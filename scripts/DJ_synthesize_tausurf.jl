#!/usr/bin/env julia

# Run the RADMC3D tausurf command to synthesize a visualization of the tau=1 surface.

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
    "--tau"
    help = "The optical depth to visualize. Default=1.0"
    default = 1.0
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using DiskJockey.constants
using DiskJockey.model

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]

pars = convert_dict(config["parameters"], config["model"])
npix = config["npix"] # number of pixels

tau = parsed_args["tau"]

# lambda should have already been written by DJInitialize.jl
tic()
run(`radmc3d tausurf $tau incl $(pars.incl) posang $(pars.PA) npix $npix loadlambda`)
println("Tausurf runtime")
# Move the image.out into a new file, since this conflicts with the normal channelmaps output.
run(`cp image.out image_tausurf.out`)
toc()
