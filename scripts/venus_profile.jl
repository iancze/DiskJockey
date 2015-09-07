#!/usr/bin/env julia

# Use the EnsembleSampler
# This means that each likelihood evaluation is *serial*, while the walkers are parallelized.

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    # "--opt1"
    # help = "an option with an argument"
    "--p", "-p"
    help = "number of processes to add"
    arg_type = Int
    default = 0
    "--run_index", "-r"
    help = "Output run index"
    arg_type = Int
    "--chain"
    help = "Write the chain to ~/web/ directory?"
    action = :store_true
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

parsed_args = parse_args(ARGS, s)

# Since we've made this a #!/usr/bin/env julia script, we can no longer specify the extra
# processors via the julia prompt, so we need to ghost this behavior back in.
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

using JudithExcalibur.constants
using JudithExcalibur.visibilities
using JudithExcalibur.image
using JudithExcalibur.gridding
using JudithExcalibur.model


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

# go through any previously created directories and remove them
function cleardirs!(keylist::Vector{Int})
    println("Removing old directories")
    for key in keylist
        keydir = basedir * "jud$key"
        run(`rm -rf $keydir`)
    end
    println("Removed directories")
end

# make the values of run_index and config available on all processes
global run_id=run_index
global cfg=config
global kl=keylist


# Now, redo this to only load the dvarr for the keys that we need, and conjugate
dvarr = DataVis(cfg["data_file"], kl)
visibilities.conj!(dvarr)
nchan = length(dvarr)

const global species = cfg["species"]
basefmt(id::Int) = cfg["base_dir"] * @sprintf("run%02d/", id)
const global basedir = basefmt(run_id)

const global npix = cfg["npix"] # number of pixels, can alternatively specify x and y separately

# Keep track of the current home directory
const global homedir = pwd() * "/"

# make the internal Judith directory, if it doesn't exist
if !ispath(basedir)
    println("Creating ", basedir)
    mkdir(basedir)
end

# Clear all directories
cleardirs!(keylist)

# Delete the old log file (if it exists)
const logfile = outdir * "log.log"
if isfile(logfile)
    rm(logfile)
end


using Logging
# change the default logger
Logging.configure(filename=logfile, level=DEBUG)

debug("Created logfile.")

# Create the model grid
grd = cfg["grid"]
global const grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)

# Regenerate all of the static files (e.g., amr_grid.inp)
# so that they may be later copied
debug("Writing grid")
write_grid(basedir, grid)
debug("Wrote grid")

# Calculate the lnprior based upon the current parameter values
# function lnprior(pars::Parameters)
#     mu_d = 142. # [pc]
#     sig_d = 6. # [pc]
#     return -0.5 * (pars.dpc - mu_d)^2 / sig_d^2
# end


# Here, if we actually are going to be fixing distance, evaluate one image to get the interpolation
# closures
# Interpolation closures need to be global objects
# if config["fix_d"]
#     dpc = cfg["parameters"]["dpc"][1] # [pc] distance
#     M_star, r_c, T_10, q, logM_gas, ksi, incl, PA, vel, mu_RA, mu_DEC = p
# else
#     M_star, r_c, T_10, q, logM_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC = p
# end

# For each channel, also calculate the interpolation closures
# npix = cfg["npix"]
# pixsize = cfg["pixsize"] # [cm]
# pix_AU = pixsize/AU # [AU]
#
# dl = sin(pix_AU/cfg["parameters"]["dpc"][1] * arcsec)
#
# uu = fftshift(fftfreq(npix, dl)) * 1e-3 # [kλ]
# vv = fftshift(fftfreq(npix, dl)) * 1e-3 # [kλ]
#
# int_arr = Array(Function, length(keys))
# for (i, dset) in enumerate(dvarr)
#     int_arr[i] = plan_interpolate(dset, uu, vv)
# end


# This function is fed to the EnsembleSampler
# That means, using the currently available global processes, like the data visibilities,
# and it must create it's own temporary directory to write the necessary files for
# RADMC to run.
function fprob(p::Vector{Float64})

    # println("Starting in ", pwd())

    # Each walker needs to create it's own temporary directory
    # where all RADMC files will reside and be driven from
    # It only needs to last for the duration of this function, so let's use a tempdir
    keydir = mktempdir() * "/"

    # Copy all relevant configuration scripts to this subdirectory
    # these are mainly setup files that will be static throughout the run
    # they were written by JudithInitialize.jl and write_grid()
    for fname in ["radmc3d.inp", "wavelength_micron.inp", "lines.inp", "molecule_" * molnames[species] * ".inp"]
        ff = homedir * fname
        run(`cp $ff $keydir`)
    end
    # run(`cp radmc3d.inp $keydir`)
    # run(`cp wavelength_micron.inp $keydir`)
    # run(`cp lines.inp $keydir`)
    #
    ag = basedir * "amr_grid.inp"
    run(`cp $ag $keydir`)

    # mf = "molecule_" * molnames[species] * ".inp"
    # run(`cp $mf $keydir`)

    # change the subprocess to reside in this directory for the remainder of the run
    # where it will drive its own independent RADMC3D process for a subset of channels
    cd(keydir)
    # println("Changed to ", pwd())

    # Fix the following arguments: gamma, dpc
    gamma = 1.0 # surface temperature gradient exponent

    if cfg["fix_d"]
        dpc = cfg["parameters"]["dpc"][1] # [pc] distance
        M_star, r_c, T_10, q, logM_gas, ksi, incl, PA, vel, mu_RA, mu_DEC = p
    else
        M_star, r_c, T_10, q, logM_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC = p
    end

    # Enforce hard priors on physical parameters
    # Short circuit evaluation if we know the RADMC won't be valid.
    if ksi <= 0. || T_10 <= 0. || r_c <= 0.0 || M_star <= 0.0 || T_10 > 1500. || q < 0. || q > 1.0
        return -Inf
    end

    if incl < 0. || incl > 180.
        return -Inf
    end

    M_gas = 10^logM_gas

    # If we are going to fit with some parameters dropped out, here's the place to do it
    pars = Parameters(M_star, r_c, T_10, q, gamma, M_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

    # Compute parameter file using model.jl, write to disk in current directory
    write_model(pars, keydir, grid, species)

    # copied from f
    # dvarr, int_arr = data

    # We are using the following conventions: inclination ranges from
    # 0 to 180 degrees. 0 means face on, angular momentum vector pointing
    # at observer; 90 means edge on; and 180 means face on, angular momentum
    # vector pointing away from observer.
    # These are the same as the RADMC-3D conventions.
    incl = pars.incl # [deg]

    # We also adopt the RADMC-3D convention for position angle, which defines position angle
    # by the angular momentum vector.
    # A positive PA angle means the disk angular momentum vector will be
    # rotated counter clockwise (from North towards East).
    PA = pars.PA # [deg]

    vel = pars.vel # [km/s]

    # Doppler shift the dataset wavelengths to rest-frame wavelength
    beta = vel/c_kms # relativistic Doppler formula
    lams = Array(Float64, nchan)
    for i=1:nchan
        lams[i] =  dvarr[i].lam * sqrt((1. - beta) / (1. + beta)) # [microns]
    end

    write_lambda(lams, keydir) # write into current directory


    # println("Now in $keydir")
    # run(`ls`)

    # println("p ", p)
    # Run RADMC-3D, redirect output to /dev/null
    tic()
    run(`radmc3d image incl $incl posang $PA npix $npix loadlambda` |> DevNull)
    println("RADMC3D call")
    toc()

    # Read the RADMC-3D images from disk (we should already be in sub-directory)
    im = imread()

    # Convert raw images to the appropriate distance
    skim = imToSky(im, pars.dpc)

    # Apply the gridding correction function before doing the FFT
    # No shift needed, since we will shift the resampled visibilities
    corrfun!(skim)

    tic()
    lnprobs = Array(Float64, nchan)
    # Do the Fourier domain stuff per channel
    for i=1:nchan
        dv = dvarr[i]
        # FFT the appropriate image channel
        vis_fft = transform(skim, i)

        # Interpolate the `vis_fft` to the same locations as the DataSet
        # mvis = int_arr[i](dv, vis_fft)

        mvis = ModelVis(dv, vis_fft)

        # Apply the phase shift here
        phase_shift!(mvis, pars.mu_RA, pars.mu_DEC)

        # Calculate chi^2 between these two
        lnprobs[i] = lnprob(dv, mvis)
    end

    println("Interpolation time")
    toc()
    # remove the temporary directory in which we currently reside
    # println("removing $keydir")
    run(`rm -rf $keydir`)


    # Sum them all together and feed back to the master process
    lnp = sum(lnprobs)

    return lnp

end


debug("Initializing MCMC")
using Distributions
using PDMats

pp = config["parameters"]
# The parameters we'll be using
if config["fix_d"]
    params = ["M_star", "r_c", "T_10", "q", "logM_gas", "ksi", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
else
    params = ["M_star", "r_c", "T_10", "q", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
end


nparam = length(params)
starting_param = Array(Float64, nparam)
jumps = Array(Float64, nparam)

for i=1:nparam
    starting_param[i], jumps[i] = pp[params[i]]
end

# given these made up parameters, call fprob
println(fprob(starting_param))

# println("Now for the profiling.")
#
# clear_malloc_data()
#
# fprob(starting_param)
#
# @profile fprob(starting_param)
# Profile.print()

tic()
fprob(starting_param)
println("Full f Timing")
toc()
