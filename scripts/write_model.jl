#!/usr/bin/env julia

# Generate a set of model visibilities and write them to HDF5
# The follow up step is to use the python scripts `write_ALMA.py` or `write_SMA.py` to convert
# these from HDF5 and pack into numpy or FITS arrays.

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
    "--out"
    help = "The file name to write out the model visibilities."
    default = "model.hdf5"
end

parsed_args = parse_args(ARGS, s)

using JudithExcalibur.constants
using JudithExcalibur.image
using JudithExcalibur.model
using JudithExcalibur.visibilities
using JudithExcalibur.gridding

using HDF5

import YAML
config = YAML.load(open(parsed_args["config"]))

# Let the model writing be taken care of by `JudithInitialize.jl`, the RADMC-3D synthesis be taken care of by `plot_model.jl` and then this simply reads in the image, does FFT, downsamples, etc.

# read the wavelengths for all channels
fid = h5open(config["data_file"], "r")
lams = read(fid["lams"]) # [μm]
close(fid)
nchan = length(lams)

# Read the parameters from the config file
pp = config["parameters"]

a = (st)->pp[st][1]

# The parameters we'll be using
pars = Parameters(a("M_star"), a("r_c"), a("T_10"), a("q"), a("gamma"), 10.^a("logM_CO"), a("ksi"), a("dpc"), a("incl"), a("PA"), a("vel"), a("mu_RA"), a("mu_DEC") )

im = imread()
skim = imToSky(im, pars.dpc)
corrfun!(skim) # alpha = 1.0

dvarr = DataVis(config["data_file"])
mvarr = Array(DataVis, nchan)
chi2s = Array(Float64, nchan)

for i=1:nchan

    dv = dvarr[i]

    # FFT the appropriate image channel
    vis_fft = transform(skim, i)

    # Interpolate the `vis_fft` to the same locations as the DataSet
    mvis = ModelVis(dv, vis_fft)

    # Apply the phase correction here, since there are fewer data points
    phase_shift!(vis_fft, pars.mu_RA, pars.mu_DEC)

    dvis = visibilities.ModelVis2DataVis(mvis)

    # Complex conjugate for SMA convention
    visibilities.conj!(dvis)
    mvarr[i] = dvis

    chi2s[i] = visibilities.chi2(dv, mvis)
end

println("Chi^2 s ", chi2s)
N = nchan * 2 * length(dvarr[1].VV)
println("Chi^2 :", sum(chi2s))
println("Reduced Chi^2 ", sum(chi2s)/N)

visibilities.write(mvarr, parsed_args["out"])