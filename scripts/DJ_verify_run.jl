#!/usr/bin/env julia

# Because we've been having so many errors with launching, run a check script.
# Load config, print the parameters that we will be sampling and the mean starting position of the walkers.

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using Statistics
using DiskJockey.model
using DiskJockey.visibilities
using DiskJockey.constants
using NPZ

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]
params = config["parameters"]
fix_params = config["fix_params"]
reg_params = registered_params[model]

println("Running model $model with the $species $transition transition.")
println()

println("We are fixing the following parameters to these values:")
for param in fix_params
  println("$param : ", params[param])
end

println()

fit_params = filter(x->∉(x,fix_params), reg_params)
n_fit = length(fit_params)

pos0 = npzread(config["pos0"])
ndim, nwalkers = size(pos0)
means = mean(pos0, dims=2) # Find the mean walker position for each parameter

@assert ndim == n_fit "Number of specified parameters ($n_fit) does not match the number of dimensions in the walker file ($ndim)."


println("We are fitting with the following $n_fit parameters:")
for (par, mu) in zip(fit_params, means)
  println("$par : $mu")
end

println("Using $nwalkers walkers in the MCMC sampling.")

# load data and figure out how many channels
dvarr = DataVis(config["data_file"])
nchan = length(dvarr)

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

println("Using the following velocity channels: ", vels[vel_mask])
