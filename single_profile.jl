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

using constants
using visibilities
using image
using gridding
using model

# This is the likelihood function called by each individual process
function f(dv::DataVis, key::Int, p::Parameters)

    # Unpack these variables from p
    incl = p.incl # [deg]
    vel = p.vel # [km/s]
    PA = p.PA # [deg] Position angle runs counter clockwise, due to looking at sky.
    npix = 256 # number of pixels, can alternatively specify x and y separately

    # Doppler shift the dataset wavelength to rest-frame wavelength
    beta = vel/c_kms # relativistic Doppler formula
    lam0 =  dv.lam * sqrt((1. - beta) / (1. + beta)) # [microns]

    tic()
    # Run RADMC3D, redirect output to /dev/null
    run(`radmc3d image incl $incl posang $PA npix $npix lambda $lam0`) # |> DevNull)
    println("RADMC3D time")
    toc()

    tic()
    # Read the RADMC3D image from disk
    im = imread()

    # Convert raw image to the appropriate distance
    skim = imToSky(im, p.dpc)

    # Apply the gridding correction function before doing the FFT
    corrfun!(skim, 1.0, p.mu_RA, p.mu_DEC) # alpha = 1.0
    println("Image read time")
    toc()

    tic()
    # FFT the appropriate image channel
    vis_fft = transform(skim)
    println("FFT time")
    toc()

    # Interpolate the `vis_fft` to the same locations as the DataSet
    tic()
    mvis = ModelVis(dv, vis_fft)
    println("Interpolate time")
    toc()

    # Apply the phase correction here, since there are fewer data points
    phase_shift!(mvis, p.mu_RA, p.mu_DEC)

    # Calculate chi^2 between these two
    tic()
    println("lnprob time")
    lnp = lnprob(dv, mvis)
    toc()

    return lnp
end

# Regenerate all of the static files (e.g., amr_grid.inp)
# so that they may be later copied
write_grid("")

key = 12

# call the initfunc with a chosen key, returning the data
dset = DataVis("data/V4046Sgr.hdf5", key)
visibilities.conj!(dset)

#From Rosenfeld et al. 2012, Table 1
M_star = 1.75 # [M_sun] stellar mass
r_c =  45. # [AU] characteristic radius
T_10 =  115. # [K] temperature at 10 AU
q = 0.63 # temperature gradient exponent
gamma = 1.0 # surface temperature gradient exponent
M_CO = 0.933 # [M_earth] disk mass of CO
ksi = 0.14 # [km/s] microturbulence
dpc = 73.0
incl = 147. # [degrees] inclination
vel = -31.18 # [km/s]
PA = -14.
mu_RA = 0.22 # [arcsec]
mu_DEC = -0.57 # [arcsec]

# Turn the parameters in the YAML file into the parameters object
# The code will only fit the parameters listed in the file

pars = Parameters(M_star, r_c, T_10, q, gamma, M_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

tic()
write_model(pars, "")

println("Model writing time")
toc()

# don't bother copying all of the files written by write_model into a subdirectory

# given these made up parameters, call fprob
println(f(dset, key, pars))
println("Now for the profiling.")

clear_malloc_data()
f(dset, key, pars)


# @profile f(dset, key, pars)
# Profile.print()
