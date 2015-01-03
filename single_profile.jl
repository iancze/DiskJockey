# lnprob evaluation for V4046Sgr
# for a single process only, designed for profiling purposes
# that means the starting point is receiving the pars::Parameters object
# and the end point is returning a lnprop

using constants
using parallel
using visibilities
using image
using gridding
using model

# This is the likelihood function called by each individual process
function f(dv::DataVis, key::Int, p::Parameters)

    # Unpack these variables from p
    incl = p.incl # [deg]
    vel = p.vel # [km/s]
    PA = 90. - p.PA # [deg] Position angle runs counter clockwise, due to looking at sky.
    npix = 256 # number of pixels, can alternatively specify x and y separately

    # Doppler shift the dataset wavelength to rest-frame wavelength
    beta = vel/c_kms # relativistic Doppler formula
    lam0 =  dv.lam * sqrt((1. - beta) / (1. + beta)) # [microns]

    tic()
    # Run RADMC3D, redirect output to /dev/null
    run(`radmc3d image incl $incl posang $PA npix $npix lambda $lam0`) # |> DevNull)
    println("RADMC3D time")
    toc()

    # Read the RADMC3D image from disk
    im = imread()

    # Convert raw image to the appropriate distance
    skim = imToSky(im, p.dpc)

    # Apply the gridding correction function before doing the FFT
    corrfun!(skim, 1.0) # alpha = 1.0 (relevant for spherical gridding function)

    tic()
    im = imread()
    skim = imToSky(im, p.dpc)
    corrfun!(skim, 1.0) # alpha = 1.0 (relevant for spherical gridding function)
    println("Image reading time")
    toc()

    vis_fft = transform(skim)
    tic()
    # FFT the appropriate image channel
    vis_fft = transform(skim)
    println("FFT time")
    toc()

    # Interpolate the `vis_fft` to the same locations as the DataSet
    mvis = ModelVis(dv, vis_fft)
    tic()
    mvis = ModelVis(dv, vis_fft)
    println("Interpolate time")
    toc()

    # Apply the phase correction here, since there are fewer data points
    phase_shift!(mvis, p.mu_x, p.mu_y)

    # Calculate chi^2 between these two
    lnprob(dv, mvis)
    tic()
    println("lnprob time")
    lnp = lnprob(dv, mvis)
    toc()

    return lnp
end

# Regenerate all of the static files (e.g., amr_grid.inp)
# so that they may be later copied
write_grid()

key = 12

# call the initfunc with a chosen key, returning the data
dset = DataVis("data/V4046Sgr.hdf5", key)

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

pars = Parameters(M_star, r_c, T_10, q, gamma, M_CO, ksi, dpc, incl, PA, vel, mu_x, mu_y)

write_model(pars)
tic()
write_model(pars)
println("Model writing time")
toc()

# don't bother copying all of the files written by write_model into a subdirectory

# given these made up parameters, call fprob
println(f(dset, key, pars))
