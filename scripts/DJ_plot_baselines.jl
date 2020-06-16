#!/usr/bin/env julia
using Pkg; Pkg.activate("DiskJockey")

# Read the data file and plot the location of all of the baselines, in k\lambda
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using DiskJockey.constants
using DiskJockey.image
using DiskJockey.model
using DiskJockey.visibilities
using DiskJockey.gridding

import PyPlot.plt
using LaTeXStrings

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]

pars = convert_dict(config["parameters"], config["model"])

# Load the data file so we can use it to compute lnprob
dvarr = DataVis(config["data_file"])
dv = dvarr[1]

fig,ax = plt[:subplots](figsize=(6,6))
ax[:plot](dv.uu, dv.vv, ".")
ax[:set_xlabel](L"UU [k$\lambda$]")
ax[:set_ylabel](L"VV [k$\lambda$]")

plt[:savefig]("baselines.png")
