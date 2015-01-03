# lnprob evaluation for V4046Sgr

const global keylist = Int[i for i=1:23]

# go through any previously created directories and remove them before the start
# of the run
function cleardirs!(keylist::Vector{Int})
    println("Removing old directories")
    for key in keylist
        keydir = "jud$key"
        run(`rm -rf $keydir`)
    end
    println("Removed directories")
end

# Clear all directories
cleardirs!(keylist)

nchild = length(keylist)
addprocs(nchild)

@everywhere using constants
@everywhere using parallel
@everywhere using visibilities
@everywhere using image
@everywhere using gridding
@everywhere using model

# Delete the old log file (if it exists)
const logfile = "log.log"
if isfile(logfile)
    rm(logfile)
end

# enable the logging module
using Logging
# change the default logger
Logging.configure(filename=logfile, level=DEBUG)


@everywhere function initfunc(key)

    # Load the relevant chunk of the dataset
    dset = DataVis("data/V4046Sgr.hdf5", key)

    # Create a directory where all RADMC files will reside and be driven from
    keydir = "jud$key"
    mkdir(keydir)

    # Copy all relevant configuration scripts to this subdirectory
    # these are mainly setup files which will not change throughout the run
    run(`cp radmc3d.inp $keydir`)
    run(`cp amr_grid.inp $keydir`)
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

    # Unpack these variables from p
    incl = p.incl # [deg]
    vel = p.vel # [km/s]
    PA = 90. - p.PA # [deg] Position angle runs counter clockwise, due to looking at sky.
    npix = 128 # number of pixels, can alternatively specify x and y separately

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
    corrfun!(skim, 1.0) # alpha = 1.0 (relevant for spherical gridding function)

    # FFT the appropriate image channel
    vis_fft = transform(skim)

    # Interpolate the `vis_fft` to the same locations as the DataSet
    mvis = ModelVis(dv, vis_fft)

    # Apply the phase correction here, since there are fewer data points
    phase_shift!(mvis, p.mu_x, p.mu_y)

    # Calculate chi^2 between these two
    return lnprob(dv, mvis)
end

# Regenerate all of the static files (e.g., amr_grid.inp)
# so that they may be later copied
write_grid()

pipes = initialize(nchild, keylist, initfunc, f)
gather!(pipes)

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
    # mu_x::Float64 # [arcsec] central offset in RA
    # mu_y::Float64 # [arcsec] central offset in DEC

    # Fix the following arguments: gamma, dpc
    gamma = 1.0 # surface temperature gradient exponent
    dpc = 73.0 # [pc] distance

    # so that p coming in is
    # [M_star, r_c, T_10, dpc, incl, PA, vel]
    M_star, r_c, T_10, q, M_CO, ksi, incl, PA, vel, mu_x, mu_y = p

    # If we are going to fit with some parameters dropped out, here's the place to do it
    # the p... command "unrolls" the vector into a series of arguments
    # The parameters type carries around everything in cgs (except mu_x, mu_y)
    pars = Parameters(M_star, r_c, T_10, q, gamma, M_CO, ksi, dpc, incl, PA, vel, mu_x, mu_y)

    # Compute parameter file using model.jl, write to disk
    write_model(pars)

    # Copy new parameter files to all subdirectories
    for key in keylist
        keydir = "jud$key"
        run(`cp numberdens_co.inp $keydir`)
        run(`cp gas_velocity.inp $keydir`)
        run(`cp gas_temperature.inp $keydir`)
        run(`cp microturbulence.inp $keydir`)
    end

    distribute!(pipes, pars)
    return gather!(pipes)
end

#From Rosenfeld et al. 2012, Table 1
M_star = 1.75 # [M_sun] stellar mass
r_c =  45. # [AU] characteristic radius
T_10 =  115. # [K] temperature at 10 AU
q = 0.63 # temperature gradient exponent
gamma = 1.0 # surface temperature gradient exponent
M_CO = 0.933 # [M_earth] disk mass of CO
ksi = 0.14 # [km/s] microturbulence
dpc = 73.0
incl = 33. # [degrees] inclination
#vel = 2.87 # LSR [km/s]
vel = -31.18 # [km/s]
PA = 73.
mu_x = 0.0 # [arcsec]
mu_y = 0.0 # [arcsec]


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

using Distributions
using PDMats

starting_param = [M_star, r_c, T_10, q, M_CO, ksi, incl, PA, vel, mu_x, mu_y]
jump_param = PDiagMat([0.01, 0.3, 0.2, 0.005, 0.01, 0.01, 0.1, 0.1, 0.005, 0.02, 0.02].^2)

# println("Evaluating fprob")
# println(fprob(starting_param))
# quit!(pipes)
# quit()

# Now try optimizing the function using NLopt
# using NLopt
#
# nparam = length(starting_param)
# opt = Opt(:LN_COBYLA, nparam)
#
# max_objective!(opt, fgrad)
# ftol_abs!(opt, 0.05) # the precision we want lnprob to
#
# (optf,optx,ret) = optimize(opt, starting_param)
# println(optf, " ", optx, " ", ret)


using LittleMC

mc = MC(fp, 2000, starting_param, jump_param)

start(mc)

println(mean(mc.samples, 2))
println(std(mc.samples, 2))

runstats(mc)

write(mc, "mc.hdf5")

quit!(pipes)
