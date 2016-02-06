# Using matplotlib, make a stick figure showing the disk inclination
using consonance

push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/DiskJockey/")
push!(LOAD_PATH, "/n/home07/iczekala/DiskJockey/")
push!(LOAD_PATH, "/pool/scout0/DiskJockey/")

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

# Based upon the mass and inclination of the disk, caluclate length of the
# velocity vector.
# Say that max vel = 0.5 r_c, for .
# Rather, can we color in dots showing projected velocity? Shades of blue/red?
# Alternative would be to shade the halves of the disk based upon the velocity.

# We could also show the central star in the disk, as an orangy thing getting
# bigger and smaller. That's probably fine.


# If we're nifty, maybe we can make a rotation arrow along the side of the
# disk that can grow larger when the mass increases

function plot_stick(pars::Parameters, id::Int)

    # Projected arrow length
    al = pars.r_c * sind(pars.incl)

    major_axis = 2 * pars.r_c
    minor_axis = 2 * pars.r_c * cosd(pars.incl)


    fig, ax = plt.subplots(nrows=1, figsize=(3,3))


    ellipse = PyPlot.matplotlib[:patches][:Ellipse]((0.,0.), major_axis, minor_axis)
    ax[:add_artist](ellipse)
    ellipse[:set_facecolor]("Chocolate")

    ellipse = PyPlot.matplotlib[:patches][:Ellipse]((0.,0.), 12/cosd(pars.incl), 12)
    ax[:add_artist](ellipse)
    ellipse[:set_facecolor]("White")

    smidge = 0.3
    space = 3.

    # This should be an ellipse
    # ax[:plot](0., 0., "o", mfc="White", mec="Black", markersize=12, zorder=100)
    ax[:plot](0., 0., "o", mfc="Gold", mec="Black", markersize=10)

    ax[:annotate]("", xy=(smidge, al), xytext=(smidge, 0.), arrowprops={"width"=>0.08, "frac"=>0.2, "headwidth"=>4., "shrink"=>0.0, "facecolor"=>"black"})

    # ax[:annotate]("", xy=(-pars.r_c - space, -al), xytext=(-pars.r_c - space, 0.), arrowprops={"width"=>0.8, "frac"=>0.2, "headwidth"=>4., "shrink"=>0.0, "facecolor"=>"blue", "edgecolor"=>"blue"})

    # ax[:annotate]("", xy=(pars.r_c + space, al), xytext=(pars.r_c + space, 0.), arrowprops={"width"=>0.8, "frac"=>0.2, "headwidth"=>4., "shrink"=>0.0, "facecolor"=>"red", "edgecolor"=>"red"})

    width = 70.
    ax[:set_xlim](-width, width)
    ax[:set_ylim](-width, width)

    plt.savefig("stick.png")
    # plt.savefig(outdir * @sprintf("stick%04d.png", id))
end

p = pars[1]
p.incl = 0.
plot_stick(p, 1)

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
