#!/usr/bin/env julia

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

# Load the file
parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using DiskJockey.visibilities
using DiskJockey.constants
using DiskJockey.model

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]

dvarr = DataVis(config["data_file"])

# determine q and weights for one channel
qq = get_qq(dvarr[1])
wgt = dvarr[1].invsig.^2

using Plots
# plot the weights as a function of q
p = scatter(qq, wgt, yaxis=(:log))
savefig("wgt_vs_qq_log.png")

q = scatter(qq, wgt)
savefig("wgt_vs_qq.png")
