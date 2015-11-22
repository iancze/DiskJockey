#!/usr/bin/env julia

# Read the HDF5 file and calculate the largest baseline in u and v, separately.

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--file"
    help = "The data file."
    default = "data.hdf5"
end

# Load the file
parsed_args = parse_args(ARGS, s)

using JudithExcalibur.visibilities
using JudithExcalibur.constants

dvarr = DataVis(parsed_args["file"])

max_baseline = 0.0

for dv in dvarr
    max_uu = maximum(dv.uu)
    if max_uu > max_baseline
        max_baseline = max_uu
    end

    max_vv = maximum(dv.vv)
    if max_vv > max_baseline
        max_baseline = max_vv
    end

end

println("Max baseline ", max_baseline, " kilolambda")

# Convert this to dRA or dDEC

dRA = 1/(2 * max_baseline * 1e3) # [radians]

println("max dRA ", dRA, " radians")

println("max dRA ", dRA /arcsec, " arcsec")
