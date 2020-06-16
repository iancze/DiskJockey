#!/usr/bin/env julia

# Given some model parameters, plot the channel maps and
# integrated spectrum

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
using HDF5

import PyPlot.plt
using LaTeXStrings
import Images

# used to plot masked values
using PyCall
@pyimport numpy.ma as ma
import NaNMath

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species * transition]
model = config["model"]
params = config["parameters"]

# Choose a diverging colorscheme for in front of/ behind the plane
cmap = plt[:get_cmap]("coolwarm")
cmap_img = plt[:get_cmap]("Blues")

# Use the same configuration as plot_chmaps, but we'll use different scaling.
function plot_tausurf(img::image.TausurfImg; fname = "tausurf.png")

    plotcbar = false

    data = img.data ./ AU # Convert from cm to AU

    vmin = abs(NaNMath.minimum(data)) # [AU]
    vmax = NaNMath.maximum(data) # [AU]
    vvmax = maximum([vmin, vmax])

    norm = PyPlot.matplotlib[:colors][:Normalize](-vvmax, vvmax)
    # norm = PyPlot.matplotlib[:colors][:SymLogNorm](linthresh=1,vmin=-10, vmax=10)

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    xxx = ((Float64[i for i = 0:im_nx - 1] + 0.5) - im_nx / 2.) * img.pixsize_x ./ AU
    yyy = ((Float64[i for i = 0:im_ny - 1] + 0.5) - im_ny / 2.) * img.pixsize_y ./ AU

    ext = (xxx[1], xxx[end], yyy[1], yyy[end])

    # Figure out how many plots we'll have.
    ncols = 8
    nrows = ceil(Int, nlam / ncols)

    xx = 1.5 * 9
    dx = 1.5
    yy = (nrows + 1) * 1.5
    dy = 1.5

    fig, ax = plt[:subplots](nrows = nrows, ncols = ncols, figsize = (xx, yy))

    for row = 1:nrows
        for col = 1:ncols
            iframe = col + (row - 1) * ncols

            if col != 1 || row != nrows
                ax[row, col][:xaxis][:set_ticklabels]([])
                ax[row, col][:yaxis][:set_ticklabels]([])
            else
                ax[row, col][:set_xlabel](L"$x$ (AU)")
                ax[row, col][:set_ylabel](L"$y$ (AU)")
            end

            if iframe > nlam
                continue # just pass through for empty channels

            else
                # No flip is needed for the tausurf because it was never put into a SkyImage. It can be plotted as is.
                frame = data[:,:,iframe]

                # Stuff the actual frame into a numpy MaskedArray object, where the masked values are the NaNs.
                mframe = pycall(ma.array, Any, frame, mask = isnan.(frame))

                # If the entire image is NaN, then just skip plotting it.
                if all(isnan, frame)
                    continue
                else
                    im = ax[row, col][:imshow](mframe, extent = ext, interpolation = "none", origin = "lower", cmap = cmap, norm = norm)

                    if !plotcbar
                        # Plot the colorbar
                        cax = fig[:add_axes]([(xx - 0.35 * dx) / xx, (yy - 1.5 * dy) / yy, (0.1 * dx) / xx, dy / yy])
                        cbar = fig[:colorbar](mappable = im, cax = cax)

                        cbar[:ax][:tick_params](labelsize = 6)
                        fig[:text](0.99, (yy - 1.7 * dy) / yy, "AU", size = 8, ha = "right")
                        plotcbar = true
                    end

                end

                ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords = "axes fraction", size = 8)
            end

        end
    end

    fig[:subplots_adjust](hspace = 0.00, wspace = 0.00, top = (yy - 0.5 * dy) / yy, bottom = (0.5 * dy) / yy, left = (0.5 * dx) / xx, right = (xx - 0.5 * dy) / xx)

    plt[:savefig](fname, dpi = 600)
end


# Use the same configuration as plot_chmaps, but we'll use different scaling.
function plot_taupos(pos::image.TausurfPos; fname = "taupos.png")

    data_x = pos.data_x ./ AU # Convert from cm to AU
    data_y = pos.data_y ./ AU # Convert from cm to AU
    data_z = pos.data_z ./ AU # Convert from cm to AU

    # Figure out how many plots we'll have.
    ncols = 8
    nrows = ceil(Int, nlam / ncols)

    xx = 1.5 * 9
    dx = 1.5
    yy = (nrows + 1) * 1.5
    dy = 1.5

    fig, ax = plt[:subplots](nrows = nrows, ncols = ncols, figsize = (xx, yy), sharex = true, sharey = true)

    for row = 1:nrows
        for col = 1:ncols
            iframe = col + (row - 1) * ncols

            if col != 1 || row != nrows
                ax[row, col][:xaxis][:set_ticklabels]([])
                ax[row, col][:yaxis][:set_ticklabels]([])
            else
                ax[row, col][:set_xlabel](L"$xx$ (AU)")
                ax[row, col][:set_ylabel](L"$yy$ (AU)")
            end

            if iframe > nlam
                continue # skip the plot

            else

                ax[row, col][:plot](data_x[:,:,iframe], data_y[:,:,iframe], "k.", ms = 2)
                ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords = "axes fraction", size = 8)
            end

        end
    end

    fig[:subplots_adjust](hspace = 0.00, wspace = 0.00, top = (yy - 0.5 * dy) / yy, bottom = (0.5 * dy) / yy, left = (0.5 * dx) / xx, right = (xx - 0.5 * dy) / xx)

    plt[:savefig](fname, dpi = 600)
end

"Plot the (x,y), (y,z), and (x,z) views through the disk, along with the (projected) direction of the observer. Assumes that P.A. of the synthesized disk is 0."
function plot_taupos_tryptych(pos::image.TausurfPos, img::image.SkyImage, incl::Real; fname = "taupos_tryptych.png")
    # Figure out how many plots we'll have.

    data_x = pos.data_x ./ AU # Convert from cm to AU
    data_y = pos.data_y ./ AU # Convert from cm to AU
    data_z = pos.data_z ./ AU # Convert from cm to AU

    ncols = 4
    nrows = nlam

    vvmax = maximum(abs, img.data)
    norm = PyPlot.matplotlib[:colors][:Normalize](0, vvmax)

    # Image needs to be flipped along RA dimension
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    dx = 2.5
    xx = dx * (ncols + 1)
    dy = 2.5
    yy = (nrows + 1) * dy

    fig, ax = plt[:subplots](nrows = nrows, ncols = ncols, figsize = (xx, yy))

    for row = 1:nrows

        ax[row,1][:set_xlabel](L"$\Delta \alpha$ ('')")
        ax[row,1][:set_ylabel](L"$\Delta \delta$ ('')")
        ax[row,2][:set_xlabel](L"$x$ (AU)")
        ax[row,2][:set_ylabel](L"$y$ (AU)")
        ax[row,3][:set_xlabel](L"$y$ (AU)")
        ax[row,3][:set_ylabel](L"$z$ (AU)")
        ax[row,4][:set_xlabel](L"$x$ (AU)")
        ax[row,4][:set_ylabel](L"$z$ (AU)")

        # Plot the SkyImage
        frame = flipdim(img.data[:,:,row], 2)
        im = ax[row, 1][:imshow](frame, extent = ext, interpolation = "none", origin = "lower", cmap = cmap_img, norm = norm)

        # Plot and set the bounds to have equal aspect area.
        ax[row, 2][:plot](data_x[:,:,row], data_y[:,:,row], "k.", ms = 2)
        ax[row, 2][:set_aspect]("equal", "datalim")

        ax[row, 3][:plot](data_y[:,:,row], data_z[:,:,row], "k.", ms = 2)
        ax[row, 3][:set_aspect]("equal", "datalim")

        ax[row, 4][:plot](data_x[:,:,row], data_z[:,:,row], "k.", ms = 2)
        ax[row, 4][:set_aspect]("equal", "datalim")

        ax[row, 1][:annotate](@sprintf("%.1f", vels[row]), (0.1, 0.8), xycoords = "axes fraction", size = 8)
    end

    fig[:subplots_adjust](hspace = 0.35, wspace = 0.35, top = (yy - 0.5 * dy) / yy, bottom = (0.5 * dy) / yy, left = (0.5 * dx) / xx, right = (xx - 0.5 * dy) / xx)

    plt[:savefig](fname, dpi = 120)
end

# Read the image
im = taureadImg()
# Do the velocity conversion here for plot labels
global nlam = length(im.lams)
# convert wavelengths to velocities
global vels = c_kms * (im.lams .- lam0) / lam0

pos = taureadPos()

# Make sure we are always plotting channels from blueshift to redshift.
if vels[2] < vels[1]
  println("Flipping channel maps to plot in order of increasing velocity.")
  im.data = im.data[:,:,end:-1:1]

  pos.data_x = pos.data_x[:,:,end:-1:1]
  pos.data_y = pos.data_y[:,:,end:-1:1]
  pos.data_z = pos.data_z[:,:,end:-1:1]

  vels = vels[end:-1:1]
end

pars = convert_dict(config["parameters"], config["model"])

# Assume that image has already been generated by synthesize_model.jl
img = imread()
skimg = imToSky(img, pars.dpc)

# plot_tausurf(im)
plot_taupos_tryptych(pos, skimg, params["incl"])
