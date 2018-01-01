#!/usr/bin/env julia

# Given some model parameters, plot the channel maps and
# integrated spectrum

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
    "--linear"
    help = "Plot the regular linear channel maps."
    action = :store_true
    "--log"
    help = "Plot the channel maps with log stretch."
    action = :store_true
    "--blur"
    help = "Plot the channel maps, convolved with a fake beam."
    action = :store_true
    "--spectrum"
    help = "Plot the spatially-integrated spectrum."
    action = :store_true
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using DiskJockey.constants
using DiskJockey.image
using DiskJockey.model
using HDF5

import PyPlot.plt
using LaTeXStrings
import Images

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]

# Read the incination, position angle, and (grid) radius of the disk. Plot ellipses.
pars = config["parameters"]
r_c = pars["r_c"]
incl = pars["incl"]
PA = pars["PA"]
dpc = pars["dpc"]

# from matplotlib.patches import Ellipse
#     ax[:add_artist](PyPlot.matplotlib[:patches][:Ellipse](xy=xy, width=BMIN, height=BMAJ, angle=BPA, facecolor="0.8", linewidth=0.2))
width = r_c/dpc # [arcsec]
height = width * cosd(incl) # [arcsec]

# cmap = plt[:get_cmap]("viridis")
# cmap = plt[:get_cmap]("inferno")
# cmap = plt[:get_cmap]("plasma")
cmap = plt[:get_cmap]("Blues")

# First contour is at 3 sigma, and then contours go up (or down) in multiples of spacing
function get_levels(rms::Float64, vmax::Float64, spacing=3)
    levels = Float64[]

    val = 3 * rms
    while (val < vmax)
        append!(levels, [val])
        val += rms * spacing
    end

    return levels
end

# function plot_beam(ax, BMAJ, BMIN, xy=(1,-1))
#     BMAJ = 3600. * header["BMAJ"] # [arcsec]
#     BMIN = 3600. * header["BMIN"] # [arcsec]
#     BPA =  header["BPA"] # degrees East of North
#     # from matplotlib.patches import Ellipse
#     ax[:add_artist](PyPlot.matplotlib[:patches][:Ellipse](xy=xy, width=BMIN, height=BMAJ, angle=BPA, facecolor="0.8", linewidth=0.2))
# end

# Plot the channel maps using sky convention
"""Plot the channel maps using the sky convention. If log is true, plot intensity using
a log scale."""
function plot_chmaps(img::image.SkyImage; log=false, contours=true, fname="channel_maps_sky.png")

    if log
        ldata = log10.(img.data + 1e-99)
        vvmax = maximum(ldata)
        norm = PyPlot.matplotlib[:colors][:Normalize](vvmax - 6, vvmax)
    else
        vvmax = maximum(abs, img.data)
        # println(vmin, " ", vmax, " ", vvmax)
        norm = PyPlot.matplotlib[:colors][:Normalize](0, vvmax)

        if contours
            levels = get_levels(rms, vvmax)
        end
    end

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    # Image needs to be flipped along RA dimension
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    # Figure out how many plots we'll have.
    ncols = 8
    nrows = ceil(Int, nlam/ncols)

    xx = 1.5 * 9
    dx = 1.5
    yy = (nrows + 1) * 1.5
    dy = 1.5

    fig, ax = plt[:subplots](nrows=nrows, ncols=ncols, figsize=(xx, yy))

    for row=1:nrows
        for col=1:ncols
            iframe = col + (row - 1) * ncols

            if col != 1 || row != nrows
                ax[row, col][:xaxis][:set_ticklabels]([])
                ax[row, col][:yaxis][:set_ticklabels]([])
            else
                ax[row, col][:set_xlabel](L"$\Delta \alpha$ ('')")
                ax[row, col][:set_ylabel](L"$\Delta \delta$ ('')")
            end

            if iframe > nlam
                # Plot a blank square if we run out of channels
                ax[row, col][:imshow](zeros((im_ny, im_nx)), cmap=cmap, vmin=0, vmax=20, extent=ext, origin="lower")

            else
                #Flip the frame for Sky convention
                frame = flipdim(img.data[:,:,iframe], 2)

                if log
                    frame += 1e-15 #Add a tiny bit so that we don't have log10(0)
                    lframe = log10.(frame)
                    im = ax[row, col][:imshow](lframe, extent=ext, interpolation="none", origin="lower", cmap=cmap, norm=norm)

                else
                    im = ax[row, col][:imshow](frame, extent=ext, interpolation="none", origin="lower", cmap=cmap, norm=norm)

                    if contours
                        ax[row, col][:contour](frame, origin="lower", colors="k", levels=levels, extent=ext, linestyles="solid", linewidths=0.2)
                    end

                end

                ax[row, col][:add_artist](PyPlot.matplotlib[:patches][:Ellipse]((0,0), width, height, PA, linewidth=0.15, facecolor="none", edgecolor="w"))

                if iframe==1
                    # Plot the colorbar
                    cax = fig[:add_axes]([(xx - 0.35 * dx)/xx, (yy - 1.5 * dy)/yy, (0.1 * dx)/xx, dy/yy])
                    cbar = fig[:colorbar](mappable=im, cax=cax)

                    cbar[:ax][:tick_params](labelsize=6)
                    fig[:text](0.99, (yy - 1.7 * dy)/yy, "Jy/beam", size=8, ha="right")
                end

                ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)
            end

        end
    end

    fig[:subplots_adjust](hspace=0.00, wspace=0.00, top=(yy - 0.5 * dy)/yy, bottom=(0.5 * dy)/yy, left=(0.5 * dx)/xx, right=(xx - 0.5 * dy)/xx)

    plt[:savefig](fname, dpi=600)

end

# Plot the spatially-integrated spectrum
function plot_spectrum(img::image.SkyImage; fname="spectrum.png")

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    spec = imToSpec(img)


    ax[:plot](vels, spec[:,2], ls="steps-mid")

    ax[:set_ylabel](L"$f_\nu$ [Jy]")
    ax[:set_xlabel](L"$v$ [km/s]")

    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig](fname)

    # Calculate integrated line intensity.
    tot = image.integrateSpec(spec, lam0)
    println("Total line flux ", tot, " Jy / km / s")
end

pars = convert_dict(config["parameters"], config["model"])

# Assume that image has already been generated by synthesize_model.jl
im = imread()
skim = imToSky(im, pars.dpc)

# Do the velocity conversion here for plot labels
global nlam = length(skim.lams)
# convert wavelengths to velocities
global vels = c_kms * (skim.lams .- lam0)/lam0

# Make sure we are always plotting channels from blueshift to redshift.
if vels[2] < vels[1]
  println("Flipping channel maps to plot in order of increasing velocity.")
  skim.data = skim.data[:,:,end:-1:1]
  vels = vels[end:-1:1]
end

if parsed_args["linear"]
    plot_chmaps(skim, fname="chmaps_linear.png", contours=false)
end

if parsed_args["log"]
    plot_chmaps(skim, fname="chmaps_log.png", log=true, contours=false)
end

if parsed_args["blur"]
    beam = config["beam"]
    rms = beam["rms"] # Jy/beam
    BMAJ = beam["BMAJ"]/2 # semi-major axis [arcsec]
    BMIN = beam["BMIN"]/2 # semi-minor axis [arcsec]
    BAVG = (BMAJ + BMIN)/2
    BPA = beam["BPA"] # position angle East of North [degrees]

    println("Beam sigma ", BAVG, " [arcsec]")

    arcsec_ster = (4.25e10)
    # Convert beam from arcsec^2 to Steradians
    global rms = rms/(pi * BMAJ * BMIN) * arcsec_ster

    println("bluring maps")
    sk_blur = blur(skim, [BAVG, BAVG])

    println("Plotting blured maps")
    plot_chmaps(sk_blur, fname="chmaps_blur.png", log=false, contours=true)
end

if parsed_args["spectrum"]
    plot_spectrum(skim)
end
