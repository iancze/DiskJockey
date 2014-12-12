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
const scale = 0.5
for i=1:nvis
    model = FTGauss(dv.uu[i], dv.vv[i], p0)
    noise = scale * (randn() + randn()*im) # Just adding in some noise here
    VV_fake[i] = model + noise
    invsig[i] = 1./scale
end

# Make the fake dataset using these Gaussian visibilities and noise terms
dv_fake = DataVis(dv.lam, dv.uu, dv.vv, VV_fake, invsig)

# # Instead of using the u,v points in the SMA dataset, actually sample somewhat
# # more uniformly on a grid.
# nfake = 40
# u = linspace(-150, 150, nfake)
# v = linspace(-150, 150, nfake)
# uu = Array(Float64, nfake^2)
# vv = Array(Float64, nfake^2)
# VV_fake = Array(Complex128, nfake^2)
# invsig = Array(Float64, nfake^2)
# const scale = 0.5
# for i=1:nfake
#     for j=1:nfake
#         k = (i -1) * nfake + j
#         uu[k] = u[i]
#         vv[k] = v[j]
#         model = FTGauss(u[i], v[j], p0)
#
#         noise = scale * (randn() + randn()*im) # Just adding in some noise here
#         VV_fake[k] = model + noise
#         invsig[k] = 1./scale
#     end
# end
#
# # Make the fake dataset using these Gaussian visibilities and noise terms
# dv_fake = DataVis(dv.lam, uu, vv, VV_fake, invsig)



write(dv_fake, "../data/V4046Sgr_fake.hdf5")
