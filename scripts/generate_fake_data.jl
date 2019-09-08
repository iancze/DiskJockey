#!/usr/bin/env julia
using Pkg; Pkg.activate("DiskJockey")

# Using the analytic form of the FT of the Gaussian, and the u,v sampling and noise
# from a real dataset, make a fake dataset.

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

# Load the file
parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using DiskJockey.gridding
using DiskJockey.visibilities
using DiskJockey.constants
using DiskJockey.image
using DiskJockey.constants
using DiskJockey.gauss

import PyPlot
import PyPlot.plt
using LaTeXStrings


# Load the dataset.
dvarr = DataVis(config["data_file"])
nvis = length(dvarr[1].VV)
nlam = 1 #length(dvarr)
lam0 = dvarr[1].lam

# Realistic Gaussian will have scale dimensions
const mu_RA = 1.0 # [arcsec]
const mu_DEC = -0.5 # [arcsec]
# const mu_RA = 0. # [arcsec]
# const mu_DEC = 0. # [arcsec]
const s_x = 1.0 # [arcsec]
const s_y = 3.0 # [arcsec]
const rho = 0.0
const theta = 30 # deg rotation
# const theta = 0.0 # deg rotation
const p0 = [mu_RA, mu_DEC, s_x, s_y, rho] # [arcsec]
# const p0 = [0.0, 0.0, s_x, s_y, rho] # [arcsec]

# Make an image of the fake dataset
# full span of the image
ra = fftspace(10., 256) # [arcsec]
dec = fftspace(10., 256) # [arcsec]

# convert ra and dec in [arcsec] to radians, and then take the sin to convert to ll, mm
ll = sin.(ra * arcsec) # direction cosines
mm = sin.(dec * arcsec)

# The natural, shifted Gaussian image
img = imageGauss(ll, mm, p0, 1, theta=theta)
skim = SkyImage(img, ra, dec, lam0)


# Because the sky convention is different than the way the SkyImage is stored,
# we need to flip the array for plotting
fig, ax = plt[:subplots](nrows=2, figsize=(5, 8))

ext = (skim.ra[end], skim.ra[1], skim.dec[1], skim.dec[end])
ax[1][:imshow](flipdim(skim.data[:,:,1], 2), interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext)
ax[1][:contour](flipdim(skim.data[:,:,1], 2), origin="lower", extent=ext)
ax[1][:set_title]("Sky Projection")
ax[1][:set_xlabel](L"$\alpha$ [arcsec]")
ax[1][:set_ylabel](L"$\delta$ [arcsec]")
ax[1][:set_aspect]("equal", "datalim")

ext = (ll[1], ll[end], mm[1], mm[end])
ax[2][:imshow](skim.data[:,:,1], interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext)
ax[2][:contour](skim.data[:,:,1], origin="lower", extent=ext)
ax[2][:set_title]("Raw Array")
ax[2][:set_xlabel](L"$ll$ [rad]")
ax[2][:set_ylabel](L"$mm$ [rad]")
ax[2][:set_aspect]("equal", "datalim")

fig[:subplots_adjust](left=0.15, right=0.85, hspace=0.25)
plt[:savefig]("gaussian_img_array.png")


# Using the analytic formula for the FT of the Gaussian, compute the true
# visibilities sampled at the u,v points in chosen dataset
dvarr_fake = Array{DataVis}(nlam)

# Noise scaling
scale = 1.0

# The V4046 Sgr dataset
dv = dvarr[1]
# Make a new array for each DataVis, otherwise references will be kept and duplicated
VV_fake = Array{Complex128}(nvis)
invsig = Array{Float64}(nvis)
for i=1:nvis
    model = FTGauss(dv.uu[i], dv.vv[i], p0, 1, theta=theta)
    noise = 0 # scale * (randn() + randn()*im) # Just adding in some noise here
    VV_fake[i] = model + noise
    invsig[i] = 1./scale
end
dvarr_fake[1] = DataVis(dv.lam, dv.uu, dv.vv, VV_fake, invsig)

visibilities.write(dvarr_fake, "data_fake.hdf5")
