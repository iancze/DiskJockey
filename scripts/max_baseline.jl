#!/usr/bin/env julia

# Read the HDF5 file and calculate the largest baseline in u and v, separately.

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--file"
    help = "The data file."
    default = "data.hdf5"
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

# Load the file
parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using JudithExcalibur.visibilities
using JudithExcalibur.constants

dvarr = DataVis(parsed_args["file"])

max_base = max_baseline(dvarr)

println("Max baseline ", max_base, " kilolambda")

# Convert this to dRA or dDEC
dRA_max = 1/(nyquist_factor * max_base * 1e3) / arcsec # [arcsec]

mu_d = config["dpc_prior"]["mu"]
sig_d = config["dpc_prior"]["sig"]

dlow = mu_d - 3. * sig_d
dhigh = mu_d + 3. * sig_d

npix = config["npix"]
println("At a resolution of  npix=$npix, the maximum outer radii for the model grid are")

# These are upper limits on the total width [AU] of the image at each distance. If the disk is large enough that it exceedes these radii, then we will need to use more pixels in the image (ie, it's very resolved).
println("dpc ", dlow, " r_out ", dlow * dRA_max * npix/(2 * 1.1))
println("dpc ", dhigh, " r_out ", dhigh * dRA_max * npix/(2 * 1.1))

# Then to an upper limit on the physical width of the image, given the current distance.
# phys_width_lim = pars.dpc * dRA_max * npix # [AU]

# Now see if the image is larger than this
# if (1.1 * 2 * grd["r_out"]) > phys_width_lim
#     println("Proposed disk r_out too large for given distance and number of pixels. Increase number of pixels in image to sample sufficiently high spatial frequencies. ", pars.dpc, " ", pars.r_c, " ", phys_width_lim)
#     return -Inf
# end
