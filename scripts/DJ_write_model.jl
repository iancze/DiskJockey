#!/usr/bin/env julia

# Generate a set of model and residual visibilities and then write them to UVHDF5 file format
# More about this format, including scripts to convert to UVFITS and CASA Measurement set,
# can be found at https://github.com/Astrochem/UVHDF5

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
    "--out-model"
    help = "The file name to write out the model visibilities."
    default = "model.hdf5"
    "--out-resid"
    help = "The file name to write out the residual visibilities."
    default = "resid.hdf5"
end

parsed_args = parse_args(ARGS, s)

using DiskJockey.constants
using DiskJockey.image
# using DiskJockey.model
using DiskJockey.visibilities
using DiskJockey.gridding

using HDF5

import YAML
config = YAML.load(open(parsed_args["config"]))

# Let the model writing be taken care of by `DJInitialize.jl`, the RADMC-3D synthesis be taken care of by `plot_model.jl` and then this simply reads in the image, does FFT, downsamples, etc.

# read the wavelengths for all channels
fid = h5open(config["data_file"], "r")
nchan = length(read(fid["freqs"]))
close(fid)

# Read the parameters from the config file
pp = config["parameters"]

a = (st)->pp[st][1]

# The parameters we'll be using
# pars = Parameters(a("M_star"), a("r_c"), 100., a("q"), a("gamma"), 10.^a("logM_CO"), a("ksi"), a("dpc"), a("incl"), a("PA"), a("vel"), a("mu_RA"), a("mu_DEC") )
mu_RA = a("mu_RA")
mu_DEC = a("mu_DEC")
dpc = a("dpc")

im = imread()
skim = imToSky(im, dpc)
corrfun!(skim) # alpha = 1.0

# For *this purpose only*, read in the flagged data, so that we can export a model for these
# visibilities
dvarr = DataVis(config["data_file"], true)
# Do this as we do in `mach_three.jl`
for dset in dvarr
    # Conjugation is necessary for the SMA and ALMA
    visibilities.conj!(dset) # Swap UV convention
end

mvarr = Array(DataVis, nchan)
chi2s = Array(Float64, nchan)

for i=1:nchan

    dv = dvarr[i]

    # FFT the appropriate image channel
    vis_fft = transform(skim, i)

    # Interpolate the `vis_fft` to the same locations as the DataSet
    mvis = ModelVis(dv, vis_fft)

    # Apply the phase correction here, since there are fewer data points
    phase_shift!(mvis, mu_RA, mu_DEC)

    dvis = visibilities.ModelVis2DataVis(mvis)

    mvarr[i] = dvis

    chi2s[i] = visibilities.chi2(dv, mvis)
end

# Now generate the residual visibilities
rvarr = ResidVis(dvarr, mvarr)

# Now swap the model and residual visibilities back to ALMA/SMA convetion
for i=1:nchan
    visibilities.conj!(mvarr[i])
    visibilities.conj!(rvarr[i])
end

println("Chi^2 s ", chi2s)
N = nchan * 2 * length(dvarr[1].VV)

# Only use the unmasked channels in the chi2 calculation
if haskey(config, "exclude")
    exclude = config["exclude"]
    # which channels of the dset to fit
    keylist = filter(x->(!in(x, exclude)), Int[i for i=1:nchan])
else
    keylist = Int[i for i=1:nchan]
end

println("Note: includes flagged visibilities!")
chi2s = chi2s[keylist]
println("Unmasked Chi^2s", chi2s)

println("Chi^2 :", sum(chi2s))
println("Reduced Chi^2 ", sum(chi2s)/N)

visibilities.write(mvarr, parsed_args["out-model"])
visibilities.write(rvarr, parsed_args["out-resid"])

# copy visibility flags
visibilities.copy_flags(config["data_file"], parsed_args["out-model"])
visibilities.copy_flags(config["data_file"], parsed_args["out-resid"])
