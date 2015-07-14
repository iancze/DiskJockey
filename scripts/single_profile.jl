# lnprob evaluation
# for a single process only, designed for profiling purposes
# that means the starting point is receiving the pars::Parameters object
# and the end point is returning a lnprop

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--run_index", "-r"
    help = "Output run index"
    arg_type = Int
    "config"
    help = "a YAML configuration file"
    required = true
end

parsed_args = parse_args(ARGS, s)

import YAML
cfg = YAML.load(open(parsed_args["config"]))
outfmt(run_index::Int) = cfg["out_base"] * @sprintf("run%02d/", run_index)
basefmt(id::Int) = cfg["base_dir"] * @sprintf("run%02d/", id)

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

# make the directories
println("Creating ", outdir)
mkdir(outdir)

@everywhere basefmt(id::Int) = cfg["base_dir"] * @sprintf("run%02d/", id)

@everywhere const global basedir = basefmt(run_index)

# make the internal Judith directory, if it doesn't exist
if !ispath(basedir)
    println("Creating ", basedir)
    mkdir(basedir)
end

using constants
using visibilities
using image
using gridding
using model

# This is the likelihood function called by each individual process
function f(data, keys::Vector{Int}, p::Parameters)
    # split into datasets and interpolation functions
    dvarr, int_arr = data

    nkeys = length(keys)

    # We are using the following convention: inclination ranges from
    # 0 to 180 degrees. 0 means face on, angular momentum vector pointing
    # at observer; 90 means edge on; and 180 means face on, angular momentum
    # vector pointing away from observer.
    # These are also the RADMC conventions.
    incl = p.incl # [deg]

    # We also adopt the RADMC convention for position angle, which defines position angle
    # by the angular momentum vector.
    # A positive PA angle means the disk angular momentum vector will be
    # rotated counter clockwise (from North towards East).
    PA = p.PA # [deg]

    vel = p.vel # [km/s]
    npix = cfg["npix"] # number of pixels, can alternatively specify x and y separately

    # Doppler shift the dataset wavelengths to rest-frame wavelength
    beta = vel/c_kms # relativistic Doppler formula
    lams = Array(Float64, length(keys))
    for i=1:nkeys
        lams[i] =  dvarr[i].lam * sqrt((1. - beta) / (1. + beta)) # [microns]
    end

    # key = keys[1]
    # keydir = basedir * "jud$key/"
    # println("Writing cameralambda into $keydir")
    write_lambda(lams, "") # write into current directory

    # Run RADMC3D, redirect output to /dev/null
    run(`radmc3d image incl $incl posang $PA npix $npix loadlambda` |> DevNull)

    # Read the RADMC3D images from disk (we should already be in sub-directory)
    im = imread()

    # Convert raw images to the appropriate distance
    skim = imToSky(im, p.dpc)

    # Apply the gridding correction function before doing the FFT
    # shifts necessary as if the image were already offset
    corrfun!(skim, p.mu_RA, p.mu_DEC) # alpha = 1.0

    lnprobs = Array(Float64, nkeys)
    # Do the Fourier domain stuff per channel
    for i=1:nkeys
        dv = dvarr[i]
        # FFT the appropriate image channel
        vis_fft = transform(skim, i)

        # Apply the phase correction here, before we do the interpolation
        phase_shift!(vis_fft, p.mu_RA, p.mu_DEC)

        # Interpolate the `vis_fft` to the same locations as the DataSet
        # mvis = ModelVis(dv, vis_fft)

        mvis = int_arr[i](dv, vis_fft)

        # Calculate chi^2 between these two
        lnprobs[i] = lnprob(dv, mvis)
    end

    # Sum them all together and feed back to the master process
    return sum(lnprobs)
end

# Create the model grid
grd = cfg["grid"]
global const grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)

# Regenerate all of the static files (e.g., amr_grid.inp)
# so that they may be later copied
write_grid(basedir, grid)

# Essentially, `mach_three.initfunc` is replicated here, since we'll just be profiling a single channel

keys = [11, 12, 13, 14, 15, 16]

@everywhere function initfunc(keys)
    # In contrast to `burma_shave.jl`, now each dataset receives several
    # channels, specified by keys

    # Load the relevant chunk of the dataset, return array
    dvarr = DataVis(cfg["data_file"], keys)
    for dset in dvarr
        # Conjugation is necessary for the SMA, methinks
        visibilities.conj!(dset) # Swap UV convention
    end

    # Create a directory where all RADMC files will reside and be driven from
    # using the first key as the directory index
    key = keys[1]
    keydir = basedir * "jud$key"
    println("Creating $keydir")
    mkdir(keydir)

    # Copy all relevant configuration scripts to this subdirectory
    # these are mainly setup files which will not change throughout the run
    run(`cp radmc3d.inp $keydir`)
    ag = basedir * "amr_grid.inp"
    run(`cp $ag $keydir`)
    run(`cp lines.inp $keydir`)
    run(`cp molecule_co.inp $keydir`)
    run(`cp wavelength_micron.inp $keydir`)

    # change the subprocess to reside in this directory for the remainder of the run
    # where it will drive its own independent RADMC3D process
    cd(keydir)


    # For each channel, also calculate the interpolation closures
    npix = cfg["npix"]
    pixsize = cfg["pixsize"] # [cm]
    pix_AU = pixsize/AU # [AU]

    dl = sin(pix_AU/cfg["parameters"]["dpc"][1] * arcsec)

    uu = fftshift(fftfreq(npix, dl)) * 1e-3 # [kλ]
    vv = fftshift(fftfreq(npix, dl)) * 1e-3 # [kλ]

    int_arr = Array(Function, length(keys))
    for (i, dset) in enumerate(dvarr)
        int_arr[i] = plan_interpolate(dset, uu, vv)
    end


    # return the array of datasets
    return (dvarr, int_arr)
end


# Load the relevant chunk of the dataset, return array
data = initfunc(keys)


# READ THE PARAMETERS FROM DISK
pp = cfg["parameters"]
# The parameters we'll be using
params = ["M_star", "r_c", "T_10", "q", "logM_CO", "ksi", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
p = Array(Float64, nparam)
jumps = Array(Float64, nparam)

for i=1:nparam
    p[i], jumps[i] = pp[params[i]]
end


# FPROB START
gamma = 1.0 # surface temperature gradient exponent
dpc = cfg["parameters"]["dpc"][1] # [pc] distance

M_star, r_c, T_10, q, logM_CO, ksi, incl, PA, vel, mu_RA, mu_DEC = p

# Enforce hard priors on physical parameters
# Short circuit evaluation if we know the RADMC won't be valid.
if ksi <= 0. || T_10 <= 0. || r_c <= 0.0 || M_star <= 0.0
    return -Inf
end

if incl < 0. || incl > 180.
    return -Inf
end

M_CO = 10^logM_CO

# If we are going to fit with some parameters dropped out, here's the place to do it
pars = Parameters(M_star, r_c, T_10, q, gamma, M_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

# Compute parameter file using model.jl, write to disk
write_model(pars, basedir, grid)

# Time it again
tic()
write_model(pars, basedir, grid)
toc()

nd = basedir * "numberdens_co.inp"
gv = basedir * "gas_velocity.inp"
gt = basedir * "gas_temperature.inp"
mt = basedir * "microturbulence.inp"

# Copy new parameter files to all subdirectories

key = keys[1]
keydir = basedir * "jud$key"
run(`cp $nd $keydir`)
run(`cp $gv $keydir`)
run(`cp $gt $keydir`)
run(`cp $mt $keydir`)


# FPROB DONE


# given these made up parameters, call fprob
println(f(data, keys, pars))

println("Now for the profiling.")

clear_malloc_data()
f(data, keys, pars)


@profile f(data, keys, pars)
Profile.print()


println("f Timing")
tic()
f(data, keys, pars)
toc()
