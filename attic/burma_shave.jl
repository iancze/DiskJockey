# The bearded devil
# Is forced
# To dwell
# In the only place
# Where they don't sell
# BURMA SHAVE

# http://burma-shave.org/jingles/

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    # "--opt1"
    # help = "an option with an argument"
    "--run_index", "-r"
    help = "Output run index"
    arg_type = Int
    # default = 0
    "--chain"
    help = "Write the chain to ~/web/ directory?"
    action = :store_true
    "config"
    help = "a YAML configuration file"
    required = true
end

parsed_args = parse_args(ARGS, s)

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

using visibilities
# load data and figure out how many channels
dvarr = DataVis(config["data_file"])
nchan = length(dvarr)

if haskey(config, "exclude")
    exclude = config["exclude"]
    const global keylist = filter(x->(!in(x, exclude)), Int[i for i=1:nchan]) # which channels of the dset to fit
else
    const global keylist = Int[i for i=1:nchan]
end

println("keylist is ", keylist)

# go through any previously created directories and remove them
function cleardirs!(keylist::Vector{Int})
    println("Removing old directories")
    for key in keylist
        keydir = basedir * "jud$key"
        run(`rm -rf $keydir`)
    end
    println("Removed directories")
end

# Since we are starting with the -p option, it is no longer necessary to
# programatically add processes.
nchild = length(keylist)
# addprocs(nchild)
println("Workers allocated ", workers())

# make the values of run_index and config available on all processes
for process in procs()
    @spawnat process global run_id=run_index
    @spawnat process global cfg=config
end

@everywhere basefmt(id::Int) = cfg["base_dir"] * @sprintf("run%02d/", id)

@everywhere const global basedir = basefmt(run_id)

# make the internal Judith directory, if it doesn't exist
if !ispath(basedir)
    println("Creating ", basedir)
    mkdir(basedir)
end

# Clear all directories
cleardirs!(keylist)

@everywhere using constants
@everywhere using parallel
@everywhere using visibilities
@everywhere using image
@everywhere using gridding
@everywhere using model

# Delete the old log file (if it exists)
const logfile = outdir * "log.log"
if isfile(logfile)
    rm(logfile)
end


using Logging
# change the default logger
Logging.configure(filename=logfile, level=DEBUG)

debug("Created logfile.")

@everywhere function initfunc(key)

    # Load the relevant chunk of the dataset
    dset = DataVis(cfg["data_file"], key)
    # Conjugation is necessary for the SMA, methinks
    visibilities.conj!(dset) # Swap UV convention

    # Create a directory where all RADMC files will reside and be driven from
    keydir = basedir * "jud$key"
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

    # return the dataset
    return dset
end

# This is the likelihood function called by each individual process
@everywhere function f(dv::DataVis, key::Int, p::Parameters)

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

    # Doppler shift the dataset wavelength to rest-frame wavelength
    beta = vel/c_kms # relativistic Doppler formula
    lam0 =  dv.lam * sqrt((1. - beta) / (1. + beta)) # [microns]

    # Run RADMC3D, redirect output to /dev/null
    run(`radmc3d image incl $incl posang $PA npix $npix lambda $lam0` |> DevNull)

    # Read the RADMC3D image from disk (we should already be in sub-directory)
    im = imread()

    # Convert raw image to the appropriate distance
    skim = imToSky(im, p.dpc)

    # Apply the gridding correction function before doing the FFT
    # shifts necessary as if the image were already offset
    corrfun!(skim, 1.0, p.mu_RA, p.mu_DEC) # alpha = 1.0 (relevant for spherical gridding function)
    # corrfun!(skim, 1.0, 0.0, 0.0) # alpha = 1.0 (relevant for spherical gridding function)

    # FFT the appropriate image channel
    vis_fft = transform(skim)

    # Apply the phase correction here, before we do the interpolation
    phase_shift!(vis_fft, p.mu_RA, p.mu_DEC)

    # Interpolate the `vis_fft` to the same locations as the DataSet
    mvis = ModelVis(dv, vis_fft)

    # Apply the phase correction here, since there are fewer data points
    # phase_shift!(mvis, p.mu_RA, p.mu_DEC)

    # Calculate chi^2 between these two
    return lnprob(dv, mvis)
end

# Create the model grid
grd = config["grid"]
global const grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)

# Regenerate all of the static files (e.g., amr_grid.inp)
# so that they may be later copied
debug("Writing grid")
write_grid(basedir, grid)
debug("Wrote grid")

debug("Initializing processes")
pipes = initialize(nchild, keylist, initfunc, f)
gather!(pipes)
debug("Initialized processes")

# Calculate the lnprior based upon the current parameter values
function lnprior(pars::Parameters)
    mu_d = 141. # [pc]
    sig_d = 7. # [pc]
    return -0.5 * (pars.dpc - mu_d)^2 / sig_d^2
end

# this function is called only on the main process, which proposes MCMC jumps
# to this function, and farms out the likelihood evaluation to all of the child
# processes
function fprob(p::Vector{Float64})

    # Here is where we make the distinction between a proposed vector of floats
    # (i.e., the parameters), and the object which defines all of the disk parameters
    # which every single subprocess will use

    # Parameters has the following definition (in model.jl)
    # M_star::Float64 # [g] stellar mass
    # r_c::Float64 # [cm] characteristic radius
    # T_10::Float64 # [K] temperature at 10 AU
    # q::Float64 # temperature gradient exponent
    # gamma::Float64 # surface temperature gradient exponent
    # M_CO::Float64 # [g] disk mass of CO
    # ksi::Float64 # [cm s^{-1}] microturbulence
    # dpc::Float64 # [pc] distance to system
    # incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    # PA::Float64 # [degrees] position angle (East of North)
    # vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    # mu_RA::Float64 # [arcsec] central offset in RA
    # mu_DEC::Float64 # [arcsec] central offset in DEC

    # Fix the following arguments: gamma, dpc
    gamma = 1.0 # surface temperature gradient exponent
    dpc = cfg["parameters"]["dpc"][1] # [pc] distance
    M_star = cfg["parameters"]["M_star"][1] # M_sun

    r_c, T_10, q, logM_CO, ksi, incl, PA, vel, mu_RA, mu_DEC = p

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

    nd = basedir * "numberdens_co.inp"
    gv = basedir * "gas_velocity.inp"
    gt = basedir * "gas_temperature.inp"
    mt = basedir * "microturbulence.inp"

    # Copy new parameter files to all subdirectories
    for key in keylist
        keydir = basedir * "jud$key"
        run(`cp $nd $keydir`)
        run(`cp $gv $keydir`)
        run(`cp $gt $keydir`)
        run(`cp $mt $keydir`)
    end

    distribute!(pipes, pars)
    return gather!(pipes) # the summed lnprob
end

# wrapper for NLopt requires gradient as an argument (even if it's not used)
function fgrad(p::Vector, grad::Vector)
    val = fprob(p)
    debug(p, " : ", val)
    return val
end

function fp(p::Vector)
    val = fprob(p)
    debug(p, " : ", val)
    return val
end


debug("Initializing MCMC")
using Distributions
using PDMats

pp = config["parameters"]
# The parameters we'll be using
params = ["r_c", "T_10", "q", "logM_CO", "ksi", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
starting_param = Array{Float64}(nparam)
jumps = Array{Float64}(nparam)

for i=1:nparam
    starting_param[i], jumps[i] = pp[params[i]]
end

jump_param = PDiagMat(jumps.^2)
jump_param = full(jump_param)

# Perturb the starting parameters
proposal = MvNormal(jump_param)
starting_param = starting_param .+ 7. * rand(proposal)

# If we've provided an empirically measured covariance matrix for the MCMC
# jump proposals, use that instead of our jumps without covariance
if haskey(config, "opt_jump")
    using NPZ
    jump_param = config["jump_scale"]^2 * npzread(config["opt_jump"])
end


# Now try optimizing the function using NLopt
# using NLopt
#
# nparam = length(starting_param)
# opt = Opt(:LN_COBYLA, nparam)
#
# max_objective!(opt, fgrad)
# ftol_abs!(opt, 0.05) # the precision we want lnprob to
#
# lower_bounds!(opt, lower)
# upper_bounds!(opt, upper)
#
# (optf,optx,ret) = optimize(opt, starting_param)
# println(optf, " ", optx, " ", ret)

using LittleMC

# Code to hook into the chain.js plot generator
if parsed_args["chain"]
    csv = open("/n/home07/iczekala/web/chain.js/mc.csv", "w")
else
    csv = open(outdir * "mc.csv", "w")
end

#Write the parameter header
writecsv(csv, params')

mc = MC(fp, config["samples"], starting_param, jump_param, csv)
debug("Initialized MCMC")

start(mc)

println(mean(mc.samples, 2))
println(std(mc.samples, 2))

runstats(mc)

debug("Acceptance: ", LittleMC.acceptance(mc))

LittleMC.write(mc, outdir * "mc.hdf5")
close(csv)

quit!(pipes)
