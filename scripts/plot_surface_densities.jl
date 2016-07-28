#!/usr/bin/env julia

using ArgParse

s = ArgParseSettings(description=" Given the output from the MCMC results, plot random draws of the proposed surface-density profiles.")
@add_arg_table s begin
    "--burn"
    help = "Amount of samples to delete from beginning of chain.npy"
    arg_type = Int
    "--draws"
    help = "How many draws of surface-density profiles to show."
    default = 15
    arg_type = Int
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using DiskJockey.model
using DiskJockey.constants

import PyPlot.plt
using LaTeXStrings

model = config["model"]
pars = convert_dict(config["parameters"], config["model"])
