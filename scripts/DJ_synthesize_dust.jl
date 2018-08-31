#!/usr/bin/env julia

# Given some model parameters, synthesize the images. Need to run DJInitialize.jl first.

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
using DiskJockey.model

model = config["model"]

pars = convert_dict(config["parameters"], config["model"])

npix = config["npix"] # number of pixels

grid = Grid(config["grid"])
(sizeau_desired, sizeau_command) = size_au(config["size_arcsec"], pars.dpc, grid) # [AU]

run(`radmc3d mctherm`)

# lambda should have already been written by DJInitialize.jl
tic()
run(`radmc3d image incl $(pars.incl) posang $(pars.PA) npix $npix sizeau $sizeau_command lambda 1300`)
println("Synthesis runtime")
toc()

run(`radmc3d sed incl $(pars.incl) posang $(pars.PA) npix $npix sizeau $sizeau_command`)
