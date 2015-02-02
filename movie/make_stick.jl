# Using matplotlib, make a stick figure showing the disk inclination
using consonance

push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")
push!(LOAD_PATH, "/n/home07/iczekala/JudithExcalibur/")

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    # "--opt1"
    # help = "an option with an argument"
    "--run_index", "-r"
    help = "Output run index"
    arg_type = Int
end

parsed_args = parse_args(ARGS, s)
run_index = parsed_args["run_index"]

using constants
using model

import PyPlot.plt
using LaTeXStrings

# let's say the x axis stretches from -70 to + 70 "AU"

arrow_length = 1.0 * pars.r_c

# Projected arrow length
# al = arrow_length * sind(pars.incl)

major_axis = pars.r_c
minor_axis = pars.r_c * cosd(incl)

# If we're nifty, maybe we can make a rotation arrow along the side of the
# disk that can grow larger when the mass increases

function plot_stick(pars::Parameters, id::Int)
    fig, ax = plt.subplots(nrows=1, figsize=(3,3))
    ax[:plot](0., 0., "ko")

    plt.savefig("stick.png")
    # plt.savefig(outdir * @sprintf("stick%04d.png", id))
end

plot_stick(pars[1], 1)

# How many frames per process?
# start = run_index * nframes_per_proc + 1
#
# ids = Int[i for i=start:(start + nframes_per_proc)]
#
# for id in ids
#     if id <= nframes
#         plot_stick(pars[id] , id)
#     else
#         break
#     end
# end
