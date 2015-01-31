push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

#Make the movie!

using constants
using image
using model

# Movie of V4046Sgr, showing 7 representative channels

global const nchan = 7
global const vels = linspace(-1.5, 1.5, nchan) # [km/s]
# CO 2-1 rest frame
lam0 = cc/230.538e9 * 1e4 # [microns]
# convert velocities to wavelengths
lams = lam0 * (vels/c_kms + 1)

function make_image(pars, id::Int)

    write_model(pars, "", grid)
    vel = pars.vel # [km/s]
    # RADMC conventions for inclination and PA
    incl = pars.incl # [deg]
    PA = pars.PA # [deg] Position angle runs counter clockwise

    run(`radmc3d image incl $incl posang $PA npix $npix loadlambda` |> DevNull)

    cp("image.out", @sprintf("image%04d.out", id))

end

M_star = 1.75 # [M_sun] stellar mass
r_c = 45. # [AU] characteristic radius
T_10 = 115. # [K] temperature at 10 AU
q = 0.63 # temperature gradient exponent
gamma = 1.0 # surface density gradient
logM_CO = 0.2 # [M_earth] disk mass of CO
ksi = 0.14 # [km/s] microturbulence
dpc = 73. # [pc] distance
incl = 147. # [degrees] inclination
PA = -17. # [degrees] position angle
vel = -31.18 # [km/s]
mu_RA = 0.2 # [arcsec] centroid location
mu_DEC = -0.6 # [arcsec]

# Parameters(M_star, r_c, T_10, q, gamma, M_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)
pars = Parameters(M_star, r_c, T_10, q, gamma, 10^logM_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

const global npix = 256 # number of pixels
const global grid = Grid(100, 32, 0.5, 800, true)
write_grid("", grid)
write_lambda(lams)

nframes = 4
ids = Int[i for i=1:nframes]
incls = linspace(0, 90., nframes)

for i in ids
    pars.incl = incls[i]
    make_image(pars, i)
end