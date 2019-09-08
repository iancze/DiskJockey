#!/usr/bin/env julia

# Generate a set of model and residual visibilities and then write them to UVHDF5 file format
# More about this format, including scripts to convert to UVFITS and CASA Measurement set,
# can be found at https://github.com/Astrochem/UVHDF5

using Pkg; Pkg.activate("DiskJockey")

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
using DiskJockey.model
using DiskJockey.visibilities
using DiskJockey.gridding

using HDF5

import YAML
config = YAML.load(open(parsed_args["config"]))

# The model writing is taken care of by `DJ_initialize.jl`, the RADMC-3D synthesis is taken care of by `DJ_plot_model.jl` and then this simply reads in the image, does FFT, downsamples, etc.

dvarr = DataVis("resid.hdf5")
nchan = length(dvarr)

for dset in dvarr
    # Conjugation is necessary for the SMA and ALMA
    visibilities.conj!(dset) # Swap UV convention
end

chi2s = Array{Float64}(undef, nchan)

for i=1:nchan
    dvis = dvarr[i]
    chi2s[i] = sum(abs2, dvis.invsig .* dvis.VV)
end

N = nchan * 2 * length(dvarr[1].VV)

# The data are stored in increasing frequency, so
# exclude: [1] means exclude the most redshifted channel
# whereas
# exclude : [nchan] excludes the most blueshifted channel
if haskey(config, "exclude")
    exclude = config["exclude"]
    # which channels of the dset to fit
    # keylist = filter(x->(!in(x, exclude)), Int[i for i=1:nchan])

    lam0 = lam0s[config["species"] * config["transition"]]
    # calculate the velocities corresponding to dvarr
    lams = Float64[dv.lam for dv in dvarr]
    vels = c_kms * (lams .- lam0)/lam0
    # get the mask
    vel_mask = generate_vel_mask(exclude, vels)
else
    # keylist = Int[i for i=1:nchan]
    vel_mask = trues(nchan)
end

chi2s = chi2s[vel_mask]
println("Reduced Chi^2 ", sum(chi2s)/N)
