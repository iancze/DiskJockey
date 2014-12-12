push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

# Using the analytic form of the FT of the Gaussian, and the u,v sampling and noise
# from a real SMA dataset, make a fake dataset.

using gridding
using visibilities
using constants
using image
using gauss_model
using constants

# Read in the dataset
dv = DataVis("../data/V4046Sgr.hdf5", 12)
nvis = length(dv.VV)

# Realistic Gaussian will have scale dimensions (fatter in x direction)
const s_x = 1.2 # [arcsec]
const s_y = 1.0 # [arcsec]
const p0 = [s_x, s_y]


# Using the analytic formula for the FT of the Gaussian, compute the true
# visibilities sampled at the u,v points in the SMA dataset
VV_fake = Array(Complex128, nvis)
invsig = Array(Float64, nvis)
for i=1:nvis
    model = FTGauss(dv.uu[i], dv.vv[i], p0)
    noise = randn() + randn()*im # Signal to noise of 1
    VV_fake[i] = model + noise
    invsig[i] = 1.
end

# Make the fake dataset using these Gaussian visibilities and noise terms
dv_fake = DataVis(dv.lam, dv.uu, dv.vv, VV_fake, invsig)

write(dv_fake, "../data/V4046Sgr_fake.hdf5")
