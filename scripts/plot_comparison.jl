#!/usr/bin/env julia


using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "image1"
    help = "The first comparison image"
    "image2"
    help = "The second comparison image"
    "--per"
    help = "Percentage to scale the residual."
    default = 0.1 # 10% is max
    "--log"
    help = "Should the intensity be plotted in log scale?"
    action = :store_true
end

parsed_args = parse_args(ARGS, s)

# Load two separate high resolution images and compare them.

using JudithExcalibur.constants
using JudithExcalibur.image
using JudithExcalibur.model

# using HDF5

import PyPlot.plt
using LaTeXStrings


img1 = imread(parsed_args["image1"])

img2 = imread(parsed_args["image2"])

nlam = length(img1.lams)

# Make sure we have the same number of channels to compare images.
@assert nlam == length(img2.lams) "Images have a different number of channels!"

# Plot the raw channel maps directly from RADMC
function compare_maps(img1::image.RawImage, img2::image.RawImage)

    (im_ny, im_nx) = size(img1.data)[1:2] # y and x dimensions of the image

    resid = img1.data .- img2.data

    ext = (1, im_nx, 1, im_ny) # Python array convention for Matplotlib

    fig, ax = plt.subplots(nrows=3, ncols=nlam, figsize=(1.5 * nlam, 5.))

    cm = plt.get_cmap("bwr_r")

    for iframe=1:nlam

        f1 = img1.data[:,:,iframe]
        f2 = img2.data[:,:,iframe]
        r = resid[:,:,iframe]

        vvmax = maxabs(vcat(f1, f2))

        # frame = img.data[:,:,iframe]
        # frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
        # max = maximum(log10(frame))
        # ax[row, col][:imshow](log10(frame), vmin=max - 6, vmax=max, interpolation="none", origin="lower", cmap=plt.get_cmap("PuBu"), extent=ext)
        # levels = linspace(max - 0.8, max, 5)

        # ax[row, col][:contour](log10(frame), origin="lower", colors="k", levels=levels, linestyles="solid", linewidths=0.2, extent=ext)

        ax[1, iframe][:imshow](f1, interpolation="none", origin="lower", cmap=cm, extent=ext, vmin=-vvmax, vmax=vvmax)

        ax[2, iframe][:imshow](f2, interpolation="none", origin="lower", cmap=cm, extent=ext, vmin=-vvmax, vmax=vvmax)

        # Scale the residuals separately

        # Set the max scale as a percentage of the original. Say 10% is the max.
        # 10% of total flux
        vvmax_new = maxabs(r)

        ax[3, iframe][:imshow](r, interpolation="none", origin="lower", cmap=cm, extent=ext, vmin=-vvmax_new, vmax=vvmax_new)

        # Now, write down what percentage this was of the original
        scale = vvmax_new/ vvmax
        ax[3, iframe][:annotate](@sprintf("%.5f", scale), (0.1, 0.8), xycoords="axes fraction", size=8)

        ax[1, iframe][:xaxis][:set_ticklabels]([])
        ax[1, iframe][:yaxis][:set_ticklabels]([])
        ax[2, iframe][:xaxis][:set_ticklabels]([])
        ax[2, iframe][:yaxis][:set_ticklabels]([])

        if iframe != 1
            ax[3, iframe][:xaxis][:set_ticklabels]([])
            ax[3, iframe][:yaxis][:set_ticklabels]([])
        else
            ax[3, iframe][:set_xlabel](L"$X$")
            ax[3, iframe][:set_ylabel](L"$Y$")
        end


    end

    fig[:subplots_adjust](hspace=0.001, wspace=0.001, top=0.95, bottom=0.05, left=0.01, right=0.99)

    plt.savefig("channel_maps_compare.png")

end


compare_maps(img1, img2)
