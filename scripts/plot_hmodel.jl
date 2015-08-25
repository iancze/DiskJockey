#!/usr/bin/env julia

# Given some model parameters, synthesize and plot the channel maps and
# integrated spectrum

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--norad"
    help = "Use the image already here."
    action = :store_true
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using JudithExcalibur.constants
using JudithExcalibur.image
using JudithExcalibur.hmodel
# using constants
# using image
# using model
using HDF5

import PyPlot.plt
using LaTeXStrings
import Images

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]

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


# Plot the raw channel maps directly from RADMC
function plot_chmaps(img::image.RawImage)

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    # CO 2-1 rest frame
    # lam0 = cc/230.538e9 * 1e4 # [microns]
    # nlam = length(img.lams)

    # convert wavelengths to velocities
    # vels = c_kms * (img.lams .- lam0)/lam0

    fig, ax = plt.subplots(nrows=2, ncols=12, figsize=(12, 2.8))

    ext = (1, im_nx, 1, im_ny) # Python array convention for Matplotlib

    for row=1:2
        for col=1:12
            iframe = col + (row - 1) * 12

            if iframe > nlam
                # Stop if we run out of channels
                break
            end

            frame = img.data[:,:,iframe]
            frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
            max = maximum(log10(frame))
            ax[row, col][:imshow](log10(frame), vmin=max - 6, vmax=max, interpolation="none", origin="lower", cmap=plt.get_cmap("PuBu"), extent=ext)
            levels = linspace(max - 0.8, max, 5)
            ax[row, col][:contour](log10(frame), origin="lower", colors="k", levels=levels, linestyles="solid", linewidths=0.2, extent=ext)

            if col != 1 || row != 2
                ax[row, col][:xaxis][:set_ticklabels]([])
                ax[row, col][:yaxis][:set_ticklabels]([])
            else
                ax[row, col][:set_xlabel](L"$X$")
                ax[row, col][:set_ylabel](L"$Y$")
            end

            ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)

        end
    end

    fig[:subplots_adjust](hspace=0.01, wspace=0.05, top=0.9, bottom=0.1, left=0.05, right=0.95)

    plt.savefig("channel_maps_raw_image.png")

end

# Plot the channel maps using sky convention
function plot_chmaps(img::image.SkyImage, fname="channel_maps_sky.png")

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    # Image needs to be flipped along RA dimension
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    # Figure out how many plots we'll have.
    ncols = 8
    nrows = iceil(nlam/ncols)

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 1.5 * nrows))

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
                # Stop if we run out of channels
                # Plot a blank square
                ax[row, col][:imshow](zeros((im_ny, im_nx)), cmap=plt.get_cmap("PuBu"), vmin=0, vmax=20, extent=ext, origin="lower")

            else
                #Flip the frame for Sky convention
                frame = fliplr(img.data[:,:,iframe])
                frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
                lframe = log10(frame)
                max = maximum(lframe)
                ix,iy = ind2sub(size(lframe), indmax(lframe))

                # ax[row, col][:imshow](lframe, extent=ext, interpolation="none", vmin=max - 6, vmax=max, origin="lower", cmap=plt.get_cmap("PuBu"))
                ax[row, col][:imshow](frame, extent=ext, interpolation="none", origin="lower", cmap=plt.get_cmap("PuBu"), norm=norm)

                # println("frame ", iframe)
                # println("levels ", levels)
                # println("extrema ", extrema(frame))
                if length(levels) > 0
                    ax[row, col][:contour](frame, origin="lower", colors="k", levels=levels, extent=ext, linestyles="solid", linewidths=0.2)
                end

                ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)
            end

        end
    end

    fig[:subplots_adjust](hspace=0.06, wspace=0.01, top=0.9, bottom=0.1, left=0.05, right=0.95)

    plt.savefig(fname)

end

# Plot the raw array for the Sky image
function plot_chmaps_data(img::image.SkyImage)

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    ll = sin(img.ra .* arcsec)
    mm = sin(img.dec .* arcsec)
    ext = (ll[1], ll[end], mm[1], mm[end])

    nlam = length(img.lams)

    fig, ax = plt.subplots(nrows=2, ncols=12, figsize=(12, 2.8))

    for row=1:2
        for col=1:12
            iframe = col + (row - 1) * 12

            if iframe > nlam
                # Stop if we run out of channels
                break
            end

            frame = img.data[:,:,iframe]
            frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
            max = maximum(log10(frame))
            ax[row, col][:imshow](log10(frame), extent=ext, vmin=max - 6, vmax=max, interpolation="none", origin="lower", cmap=plt.get_cmap("PuBu"))
            levels = linspace(max - 0.8, max, 5)
            ax[row, col][:contour](log10(frame), origin="lower", colors="k", levels=levels, extent=ext, linestyles="solid", linewidths=0.2)

            if col != 1 || row != 2
                ax[row, col][:xaxis][:set_ticklabels]([])
                ax[row, col][:yaxis][:set_ticklabels]([])
            else
                ax[row, col][:set_xlabel](L"$ll$")
                ax[row, col][:set_ylabel](L"$mm$")
            end

            ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)

        end
    end

    fig[:subplots_adjust](hspace=0.01, wspace=0.05, top=0.9, bottom=0.1, left=0.05, right=0.95)

    plt.savefig("channel_maps_sky_raw.png")

end


# Plot the spatially-integrated spectrum
function plot_spectrum(img::image.SkyImage)

    fig = plt.figure()
    ax = fig[:add_subplot](111)

    spec = imToSpec(img)

    ax[:plot](vels, spec[:,2], ls="steps-mid")
    # ax[:plot](vels, reverse(spec[:,2]), ls="steps-mid")
    ax[:set_ylabel](L"$f_\nu$ [Jy]")
    ax[:set_xlabel](L"$v$ [km/s]")

    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("spectrum.png")
end

pp = config["parameters"]
params = ["M_star", "r_c", "T_10m", "q_m", "T_10a", "q_a", "gamma", "h", "delta", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
starting_param = Array(Float64, nparam)

for i=1:nparam
    starting_param[i] = pp[params[i]][1]
end

# Convert logM_gas to M_gas
starting_param[10] = 10^starting_param[10]

pars = Parameters(starting_param...)
#
# vel = pars.vel # [km/s]
# # RADMC conventions for inclination and PA
incl = pars.incl # [deg]
PA = pars.PA # [deg] Position angle runs counter clockwise
npix = config["npix"] # number of pixels
#
# # Doppler shift the dataset wavelength to rest-frame wavelength
# beta = vel/c_kms # relativistic Doppler formula
# shift_lams =  lams .* sqrt((1. - beta) / (1. + beta)) # [microns]
#
# grd = config["grid"]
# grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)
#
# # Comment out afterward
# # nchan = 251
# # global const vels = linspace(-10., 10, nchan) # [km/s]
# # # CO 2-1 rest frame
# # lam0 = cc/230.538e9 * 1e4 # [microns]
# # # convert velocities to wavelengths
# # shift_lams = lam0 * (vels/c_kms + 1)


if !parsed_args["norad"]
    run(`radmc3d image incl $incl posang $PA npix $npix loadlambda`)
end

im = imread()

# plot_chmaps(im)

skim = imToSky(im, pars.dpc)


# Do the velocity conversion here
global nlam = length(skim.lams)
# convert wavelengths to velocities
global vels = c_kms * (skim.lams .- lam0)/lam0

# get the colorscale corresponding to the velocities
vmin, vmax = extrema(skim.data)
vvmax = maxabs(skim.data)

ldata = log10(skim.data + 1e-99)
vlmax = maximum(ldata)

# if in(config, "beam")
beam = config["beam"]
rms = beam["rms"] # Jy/beam
BMAJ = beam["BMAJ"]/2 # semi-major axis [arcsec]
BMIN = beam["BMIN"]/2 # semi-minor axis [arcsec]
BAVG = (BMAJ + BMIN)/2
BPA = beam["BPA"] # position angle East of North [degrees]

println("Beam sigma ", BAVG, " [arcsec]")

arcsec_ster = (4.25e10)
# Convert beam from arcsec^2 to Steradians
rms = rms/(pi * BMAJ * BMIN) * arcsec_ster
# println("RMS: ", rms, " vmax: ", vmax)
global levels = get_levels(rms, vmax)

norm = PyPlot.matplotlib[:colors][:Normalize](vlmax - 8, vlmax)

# println("Plotting hires maps")
# plot_chmaps(skim)

println("bluring maps")
sk_blur = blur(skim, [BAVG, BAVG])

# Now redo all this for the blurred
vmin, vmax = extrema(sk_blur.data)
global levels = get_levels(rms, vmax)
norm = PyPlot.matplotlib[:colors][:Normalize](0, vmax)

println("Plotting blured maps")
plot_chmaps(sk_blur, "channel_maps_blur.png")

plot_spectrum(skim)
