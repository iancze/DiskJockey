# All necessary variables (ie, setup vars) should be defined here, to be
# shared between any disparate processes
# push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/DiskJockey/")
# push!(LOAD_PATH, "/n/home07/iczekala/DiskJockey/")
# push!(LOAD_PATH, "/pool/scout0/DiskJockey/")

module consonance

export pars, nframes, nframes_per_proc, nchan, vels, lams, vmax, moviedir
export vmax_d, vmax_g, vmax_v


using DiskJockey.constants
using DiskJockey.model

# For odyssey use
# homedir = "/n/home07/iczekala/DiskJockey/"
# scratchdir = "/scratch/"
# outdir = "/n/home07/iczekala/DiskJockey/output/movie/" #Where the images.out and images.png are stored

# For local use
moviedir = "/home/ian/Grad/Research/Disks/DiskJockey/movie/"

# For cfa CF use
# homedir = "/pool/scout0/DiskJockey/"
# scratchdir = "/pool/cf/iczekala/scratch/"
# outdir = "/pool/scout0/DiskJockey/output/movie/"

# How many frames per process?
nframes_per_proc = 50

global const nchan = 5
global const vels = linspace(-1.5, 1.5, nchan) # [km/s]
# CO 2-1 rest frame
lam0 = cc/230.538e9 * 1e4 # [microns]
# convert velocities to wavelengths
lams = lam0 * (vels/c_kms + 1)

M_star = 1.75 # [M_sun] stellar mass
r_c = 45. # [AU] characteristic radius
r_in = 1.0 # [AU]
r_cav = 2.0 # [AU]
delta = 1.0 # no change
T_10 = 115. # [K] temperature at 10 AU
q = 0.63 # temperature gradient exponent
gamma = 1.0 # surface density gradient
logM_CO = 0.2 # [M_earth] disk mass of CO
ksi = 0.14 # [km/s] microturbulence
dpc = 73. # [pc] distance
incl = 45. # [degrees] inclination
PA = 0. # [degrees] position angle
vel = 0.0 # [km/s]
mu_RA = 0.0 # [arcsec] centroid location
mu_DEC = 0.0 # [arcsec]

# Given a starting parameter, vary it by steps dp to the low bound, then high
# bound, then back to where we started

# In the new edition, keep it at the top and bottom values for some period
function smooth_vary(start, low, high, dp, npause)
    ndown1 = round(Integer, (start - low)/dp) + 1
    nup = round(Integer, (high - low)/dp) + 1
    ndown2 = round(Integer, (high - start)/dp) + 1

    return [linspace(start, low, ndown1)' low * ones(npause)' linspace(low, high, nup)' high * ones(npause)' linspace(high, start, ndown2)' start * ones(npause)']
end

# Create a master parameter list
# First, adjust in radius
# radiuses = smooth_vary(r_c, 25., 65., 0.5)
# nradiuses = length(radiuses)

# First, adjust in inclination
incls = smooth_vary(incl, 0., 90., 1., 15)
nincls = length(incls)

# then, adjust in mass
masses = smooth_vary(M_star, 1.0, 2.5, 0.02, 15)
nmasses = length(masses)

# nframes = nincls + nmasses
#
# # Now create a giant array of Parameters objects
# pars = Array(Parameters, nframes)
#
# for i=1:nframes
#     if i <= nincls
#         pars[i] = Parameters(M_star, r_c, r_in, r_cav, delta, T_10, q, gamma, 10^logM_CO, ksi, dpc, incls[i], PA, vel, mu_RA, mu_DEC)
#     else
#         pars[i] = Parameters(masses[i - nincls], r_c, r_in, r_cav, delta, T_10, q, gamma, 10^logM_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)
#     end
# end
#
# println("There are ", nframes, " frames to be generated.")
# println("Inclinations ", incls)
# println("Masses ", masses)
# println("Radiuses ", radiuses)

# Return a normalized instance that is symmetric about 0
# function scale(data)
#     s = maximum(abs(data))
#     return norm = plt.Normalize(vmin=-s, vmax=s, clip=false)
# end

# used as arg to imshow

# Poster plot material
nframes = 5
pars = [Parameters(1.0, r_c, r_in, r_cav, delta, T_10, q, gamma, 10^logM_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC),
Parameters(1.75, r_c, r_in, r_cav, delta, T_10, q, gamma, 10^logM_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC),
Parameters(2.50, r_c, r_in, r_cav, delta, T_10, q, gamma, 10^logM_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC),
Parameters(1.75, r_c, r_in, r_cav, delta, T_10, q, gamma, 10^logM_CO, ksi, dpc, 0.0, PA, vel, mu_RA, mu_DEC),
Parameters(1.75, r_c, r_in, r_cav, delta, T_10, q, gamma, 10^logM_CO, ksi, dpc, 90.0, PA, vel, mu_RA, mu_DEC)
]

vmin_d = 0.0
vmax_d = 6.415994694015499e9

vmin_g = 0.0
vmax_g = 2.9207242096782996e11

vmin_v = 0.0
vmax_v = 45.0193854891226

end # module
