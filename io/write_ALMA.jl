push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

# Generate a set of model visibilities and write them to HDF5
# Rather than run everything in parallel, like `burma_shave.jl`, do this
# all in serial, here.

# The follow up step is to use a the python script `write_ALMA.py` to convert
# these from HDF5 and pack into numpy arrays.

using constants
using image
using model
using HDF5
using visibilities
using gridding

# read the wavelengths for all 25 channels
fid = h5open("../data/AKSco/AKSco.hdf5", "r")
lams = read(fid["lams"]) # [μm]
close(fid)
nchan = length(lams)

# Peak of best-fitting parameters
M_star = 2.42 # [M_sun] stellar mass
r_c =  11.2 # [AU] characteristic radius
T_10 =  88.16 # [K] temperature at 10 AU
q = 0.54 # temperature gradient exponent
gamma = 1.0 # surface temperature gradient exponent
logM_CO = 0.23 # [M_earth] disk mass of CO
ksi = 0.29 # [km/s] microturbulence
dpc = 141.0
incl = 108. # [degrees] inclination
PA = 141.23
vel = -26.03 # [km/s]
mu_RA = 0.055 # [arcsec]
mu_DEC = 0.0489# [arcsec]

const global grid = Grid(64, 32, 0.5, 300, true)

M_CO = 10.^logM_CO

pars = Parameters(M_star, r_c, T_10, q, gamma, M_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

vel = pars.vel # [km/s]
# RADMC conventions
incl = pars.incl # [deg]
PA = pars.PA # [deg] Position angle runs counter clockwise, due to looking at sky.
npix = 256 # number of pixels, can alternatively specify x and y separately

# Doppler shift the dataset wavelength to rest-frame wavelength
beta = vel/c_kms # relativistic Doppler formula
shift_lams =  lams .* sqrt((1. - beta) / (1. + beta)) # [microns]

write_grid("", grid)
write_model(pars, "", grid)
write_lambda(shift_lams)

files = ["lines.inp", "molecule_co.inp", "wavelength_micron.inp"]
for file in files
    cp("../" * file, file)
end

cp("../radmc3d.inp.gas", "radmc3d.inp")


run(`radmc3d image incl $incl posang $PA npix $npix loadlambda`)

im = imread()
skim = imToSky(im, pars.dpc)
corrfun!(skim, 1.0, pars.mu_RA, pars.mu_DEC) # alpha = 1.0

dvarr = DataVis("../data/AKSco/AKSco.hdf5")
mvarr = Array(DataVis, nchan)

for i=1:nchan

    # FFT the appropriate image channel
    vis_fft = transform(skim, i)

    # Apply the phase correction here, since there are fewer data points
    phase_shift!(vis_fft, pars.mu_RA, pars.mu_DEC)
    # Interpolate the `vis_fft` to the same locations as the DataSet
    mvis = ModelVis(dvarr[i], vis_fft)

    dvis = visibilities.ModelVis2DataVis(mvis)

    # Complex conjugate for SMA convention
    visibilities.conj!(dvis)
    mvarr[i] = dvis
end

visibilities.write(mvarr, "../data/AKSco/AKSco_model.hdf5")
