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

mu_d, sig_d = config["parameters"]["dpc"]

dlow = mu_d - 3. * sig_d
dhigh = mu_d + 3. * sig_d

npix = config["npix"]
println("For npix $npix, the maximum outer radii allowed are")

# These are upper limits on the total width [AU] of the image at each distance. If the disk is large enough that it exceedes these radii, then we will need to use more pixels in the image (ie, it's very resolved).
println("dpc ", dlow, " r_out ", dlow * dRA_max * npix/(2 * 1.1))
println("dpc ", dhigh, " r_out ", dhigh * dRA_max * npix/(2 * 1.1))
