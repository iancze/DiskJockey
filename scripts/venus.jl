#!/usr/bin/env julia

# Uses the EnsembleSampler to explore the posterior.
# In contrast to codes like mach_three.jl, this architecture means that within each likelihood call, the global lnprob for all channels is evaluated in *serial*, while the walkers themselves are parallelized.

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--p", "-p"
    help = "number of processes to add"
    arg_type = Int
    default = 0
    "--run_index", "-r"
    help = "Output run index"
    arg_type = Int
    "--plotly"
    help = "Send the chain samples to plot.ly?"
    action = :store_true
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
    "--cpus"
    help = "Which CPUS to add"
    arg_type = Array{Int, 1}
    eval_arg = true
end

parsed_args = parse_args(ARGS, s)

cpus = parsed_args["cpus"]

if cpus != nothing
    using ClusterManagers
    addprocs(LocalAffinityManager(;affinities=cpus))
end

# Since we've made this a #!/usr/bin/env julia script, we can no longer specify the extra
# processors via the julia executable, so we need to ghost this behavior back in.
if parsed_args["p"] > 0
    addprocs(parsed_args["p"])
end

import YAML
config = YAML.load(open(parsed_args["config"]))

outfmt(run_index::Int) = config["out_base"] * @sprintf("run%02d/", run_index)

# This code is necessary for multiple simultaneous runs on odyssey
# so that different runs do not write into the same output directory
if parsed_args["run_index"] == nothing
    run_index = 0
    outdir = outfmt(run_index)
    while ispath(outdir)
        println(outdir, " exists")
        run_index += 1
        outdir = outfmt(run_index)
    end
else
    run_index = parsed_args["run_index"]
    outdir = outfmt(run_index)
    println("Deleting old $outdir")
    run(`rm -rf $outdir`)
end

# make the output directory
println("Creating ", outdir)
mkdir(outdir)

@everywhere using JudithExcalibur.constants
@everywhere using JudithExcalibur.visibilities
@everywhere using JudithExcalibur.image
@everywhere using JudithExcalibur.gridding
@everywhere using JudithExcalibur.model
@everywhere using Base.Test

# load data and figure out how many channels
dvarr = DataVis(config["data_file"])
nchan = length(dvarr)

if haskey(config, "exclude")
    exclude = config["exclude"]
    # which channels of the dset to fit
    keylist = filter(x->(!in(x, exclude)), Int[i for i=1:nchan])
else
    keylist = Int[i for i=1:nchan]
end

# This program is meant to be started with the -p option.
nchild = nworkers()
println("Workers allocated ", nchild)

# make the values of run_index and config available on all processes
for process in procs()
    @spawnat process global run_id=run_index
    @spawnat process global cfg=config
    @spawnat process global kl=keylist
end
println("Mapped variables to all processes")

# Add radmc3d everywhere, since SLURM seems to dislike inheriting it from the system PATH
@everywhere ENV["PATH"] = ENV["PATH"] * cfg["RADMC_PATH"]

# Set the correct model
@everywhere function convert_p(p::Vector{Float64})
    gamma = cfg["parameters"]["gamma"][1]
    if cfg["fix_d"]
        return convert_vector(p, cfg["model"], cfg["parameters"]["dpc"], gamma)
    else
        return convert_vector(p, cfg["model"], -1.0, gamma)
    end
end

# Now, redo this to only load the dvarr for the keys that we need, and conjugate
@everywhere dvarr = DataVis(cfg["data_file"], kl)
@everywhere visibilities.conj!(dvarr)
@everywhere nchan = length(dvarr)

@everywhere const global species = cfg["species"]

@everywhere const global npix = cfg["npix"] # number of pixels, can alternatively specify x and y separately

# Keep track of the current home directory
@everywhere const global homedir = pwd() * "/"

# Create the model grid
@everywhere grd = cfg["grid"]
@everywhere grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)

@everywhere dpc_mu = cfg["dpc_prior"]["mu"]
@everywhere dpc_sig = cfg["dpc_prior"]["sig"]

# Only calculate the interpolation closures if we are fixing distance.
if cfg["fix_d"]

    angular_width = (1.1 * 2 * grd["r_out"])/cfg["parameters"]["dpc"][1] # [radians]

    # Simply calculate pix_AU as 1.1 * (2 * r_out) / npix
    # This is assuming that RADMC always calculates the image as 110% the full extent of the grid
    @everywhere pix_AU = (1.1 * 2 * grd["r_out"]) / npix # [AU/pixel]

    # Ignore the sin, since we use small angle approximation
    @everywhere dl = pix_AU/cfg["parameters"]["dpc"][1] * arcsec

    @everywhere uu = fftshift(fftfreq(npix, dl)) * 1e-3 # [kλ]
    @everywhere vv = fftshift(fftfreq(npix, dl)) * 1e-3 # [kλ]

    # For each channel, also calculate the interpolation closures
    @everywhere int_arr = Array(Function, nchan)
    @everywhere for (i, dset) in enumerate(dvarr)
        int_arr[i] = plan_interpolate(dset, uu, vv)
    end
end

# This function is fed to the EnsembleSampler
# That means, using the currently available global processes, like the data visibilities,
# and it must create it's own temporary directory to write the necessary files for
# RADMC to run.
@everywhere function fprob(p::Vector{Float64})

    # Convert the vector to a specific type of pars, e.g. ParametersStandard, ParametersTruncated, etc, which will be used for multiple dispatch from here on out.
    pars = convert_p(p)

    lnpr = lnprior(pars, dpc_mu, dpc_sig, grid)
    if lnpr == -Inf
        return -Inf
    end

    # Each walker needs to create it's own temporary directory
    # where all RADMC-3D files will reside and be driven from
    # It only needs to last for the duration of this function, so we use a tempdir
    keydir = mktempdir() * "/"

    # Copy all relevant configuration scripts to this subdirectory so that RADMC-3D can run.
    # these are mainly setup files that will be static throughout the run
    # they were written by JudithInitialize.jl and write_grid()
    for fname in ["radmc3d.inp", "wavelength_micron.inp", "lines.inp", "molecule_" * molnames[species] * ".inp"]
        run(`cp $(homedir)$fname $keydir`)
    end

    # change the subprocess to reside in this directory for the remainder of the run
    # where it will drive its own independent RADMC3D process for a subset of channels
    cd(keydir)

    write_grid(keydir, grid)

    # Compute disk structure files using model.jl, write to disk in current directory
    # Based upon typeof(pars), the subroutines within this function will dispatch to the correct model.
    write_model(pars, keydir, grid, species)

    # Doppler shift the dataset wavelengths to rest-frame wavelength
    beta = pars.vel/c_kms # relativistic Doppler formula
    lams = Array(Float64, nchan)
    for i=1:nchan
        lams[i] =  dvarr[i].lam * sqrt((1. - beta) / (1. + beta)) # [microns]
    end

    write_lambda(lams, keydir) # write into current directory

    # Run RADMC-3D, redirect output to /dev/null
    run(pipeline(`radmc3d image incl $(pars.incl) posang $(pars.PA) npix $npix loadlambda`, DevNull))

    # Read the RADMC-3D images from disk
    im = try
        imread() # we should already be in the sub-directory, no path required
    # If the synthesized image is screwed up, just say there is zero probability.
    catch SystemError
        println("Failed to read synthesized image for parameters ", p)
        -Inf
    end

    # The image synthesis faild for some reason or another. Clean up by deleting the director and returning -Inf.
    if im == -Inf
        run(`rm -rf $keydir`)
        return -Inf
    end

    if cfg["fix_d"]
        # After the fact, we should be able to check that the pixel size of the image is the
        # same as the one we originally calculated from the outer disk radius.
        @assert abs((im.pixsize_x/AU  - pix_AU)/pix_AU) < 1e-5

        # If we aren't fixing distance, then the pixel size is read directly from the calculated image.
    end

    # Convert raw images to the appropriate distance
    skim = imToSky(im, pars.dpc)

    # Apply the gridding correction function before doing the FFT
    # No shift needed, since we will shift the resampled visibilities
    corrfun!(skim)

    lnprobs = Array(Float64, nchan)
    # Do the Fourier domain stuff per channel
    for i=1:nchan
        dv = dvarr[i]
        # FFT the appropriate image channel
        vis_fft = transform(skim, i)

        if cfg["fix_d"]
            # Interpolate the `vis_fft` to the same locations as the DataSet
            mvis = int_arr[i](dv, vis_fft)
        else
            mvis = ModelVis(dv, vis_fft)
        end

        # Apply the phase shift here
        phase_shift!(mvis, pars.mu_RA, pars.mu_DEC)

        # Calculate the likelihood between these two using our chi^2 function.
        lnprobs[i] = lnprob(dv, mvis)
    end

    # remove the temporary directory in which we currently reside
    run(`rm -rf $keydir`)

    # Sum them all together and feed back to the master process
    lnp = sum(lnprobs) + lnpr

    return lnp

end

using NPZ

pos0 = npzread(config["pos0"])
ndim, nwalkers = size(pos0)

# Use the EnsembleSampler to do the optimization
using JudithExcalibur.EnsembleSampler

sampler = Sampler(nwalkers, ndim, fprob)

# make sure that we've loaded a pos0 with the right dimensions.
size1, size2 = size(pos0)
@assert size1==ndim "pos0 array does not match number of input dimensions."
@assert size2==nwalkers "pos0 array does not match number of walkers."

if parsed_args["plotly"]
    function f(sampler::Sampler, outdir::AbstractString)
        # Run the chain to plotly
        println("Called plotly.")
        chain_file = "$(outdir)chain.npy"
        config_file = "config.yaml"
        name = config["name"]
        try
            spawn(`plotly_walkers.py --name $name --chain $chain_file --config $config_file`)
        catch
            println("Couldn't reach plotly server.")
        end
    end
    run_schedule(sampler, pos0, config["samples"], config["loops"], outdir, f)
else
    run_schedule(sampler, pos0, config["samples"], config["loops"], outdir)
end

write_samples(sampler, outdir)
