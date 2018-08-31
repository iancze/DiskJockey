#!/usr/bin/env julia

# Calculates some limits of the problem just using the dataset and parameters in the config file.
# * the largest baseline in u and v, separately.
# * the mean velocity of the dataset.
# * the physical scale of the image, and whether it fully encapsulates the model grid for
#   all of the possible distances to explore.

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

using DiskJockey.visibilities
using DiskJockey.constants
using DiskJockey.model

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
grid = Grid(config["grid"])

function wl_to_vel{T}(wl::T)
    return c_kms * (wl - lam0)/lam0
end

function vel_to_wl{T}(vel::T)
    beta = vel/c_kms # relativistic Doppler formula
    return lam0 * sqrt((1. - beta) / (1. + beta)) # [microns]
end

dvarr = DataVis(config["data_file"])

lam_start = dvarr[1].lam
lam_end = dvarr[end].lam

vmin = wl_to_vel(lam_start)
vmax = wl_to_vel(lam_end)
println("Dataset channels are velocities from $vmin to $vmax and span ", vmax - vmin, " km/s.")
println("Midpoint is ", (vmax + vmin)/2, " km/s.")

max_base = max_baseline(dvarr)

println("Max baseline ", max_base, " kilolambda")


# Given the size of the image in arcseconds and the number of pixels, figure out the nyquist
# sampling and whether or not this is sufficiently higher than our largest baseline
dRA = config["size_arcsec"] / config["npix"] # [arcsec/pix]

# Convert this to dRA or dDEC
dRA_max = 1/(nyquist_factor * max_base * 1e3) / arcsec # [arcsec / pix]

@assert dRA < dRA_max "Nyquist sampling not satisfied: Either increase the number of pixels in your image, or decrease the angular size of your image. dRA ($dRA [arcsec/pix]) is greater than dRA_max ($dRA_max [arcsec/pix]), which is determined from the maximum baseline of your dataset."

println("Nyquist sampling satisfied. dRA: $dRA [arcsec/pix] ; dRA_max: $dRA_max [arcsec/pix]")

# Given the angular size of the image, calculate what the physical size of the image would be
# at the closest source distance, and make sure that this is still larger than outer radius of
# the model grid.

dpc = config["parameters"]["dpc"]

# Calculate the physical size of the image at the closer distance, and make sure it's still larger
# than 110% of the outer grid radius
sizeau = 0.5 * config["size_arcsec"] * dpc # [AU]
outer_radius = 1.1 * config["grid"]["r_out"]
@assert sizeau > outer_radius "The angular size of your image is too small to fully encapsulate the model grid at the distance to the source. Increase the angular size of your image via `size_arcsec` or decrease the outer radius of your model grid. Size at the distance $dpc: $sizeau [AU]; outer radius of the grid + 10% $outer_radius [AU]"

println("Image size satisfied. Half-Image size at $dpc distance: $sizeau [AU]; outer radius of the grid + 10%: $outer_radius [AU]")
