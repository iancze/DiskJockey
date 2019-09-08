#!/usr/bin/env julia
using Pkg; Pkg.activate("DiskJockey")


using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "image1"
    help = "The first comparison image"
    "config1"
    help = "a YAML configuration file for the first image"
    "image2"
    help = "The second comparison image"
    "config2"
    help = "a YAML configuration file for the second image"
    "--per"
    help = "Percentage to scale the residual."
    default = 0.1 # 10% is max
    "--log"
    help = "Should the intensity be plotted in log scale?"
    action = :store_true
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config1"]))
config2 = YAML.load(open(parsed_args["config2"]))

# Load two separate high resolution images and compare them, computing the difference and the difference in lnprob.

using DiskJockey.constants
using DiskJockey.image
using DiskJockey.model
using DiskJockey.visibilities
using DiskJockey.gridding

import PyPlot.plt
using LaTeXStrings

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]

pars = convert_dict(config["parameters"], config["model"])
pars2 = convert_dict(config2["parameters"], config2["model"])

img1 = imread(parsed_args["image1"])
img2 = imread(parsed_args["image2"])

nchan = length(img1.lams)
# Make sure we have the same number of channels to compare images.
@assert nchan == length(img2.lams) "Images have a different number of channels!"

# Convert raw images to the appropriate distance
skimg1 = imToSky(img1, pars.dpc)
skimg2 = imToSky(img2, pars2.dpc)
skresid = skimg1 - skimg2

# Store maximums for plotting
max_img = maximum(abs, skimg1.data)
if maximum(abs, skimg2.data) > max_img
    max_img = maximum(abs, skimg2.data)
end
max_resid = maximum(abs, skresid.data)

# Apply the gridding correction function before doing the FFT
# No shift needed, since we will shift the resampled visibilities
# The skimg1, skimg2 are for plotting.
skimg01 = corrfun(skimg1)
skimg02 = corrfun(skimg2)

# To store the FFT for each channel
vis_fft1 = Array{FullModelVis}(nchan)
vis_fft2 = Array{FullModelVis}(nchan)
vis_resid = Array{FullModelVis}(nchan)

# To store the computed lnprob for each channel
lnprobs1 = Array{Float64}(nchan)
lnprobs2 = Array{Float64}(nchan)
lnprobs_diff = Array{Float64}(nchan)

# Load the data file so we can use it to compute lnprob
dvarr = DataVis(config["data_file"])
for dset in dvarr
    # Conjugation is necessary for the SMA and ALMA
    visibilities.conj!(dset) # Swap UV convention
end

half_pix = config["size_arcsec"] / (2 * config["npix"])

# Fourier transform each channel, while also noting the max scale values
for i=1:nchan
    dv = dvarr[i]
    # FFT the appropriate image channel
    vis_fft1[i] = transform(skimg01, i)
    vis_fft2[i] = transform(skimg02, i)
    vis_resid[i] = vis_fft1[i] - vis_fft2[i]

    mvis1 = ModelVis(dv, vis_fft1[i])
    mvis2 = ModelVis(dv, vis_fft2[i])

    # Apply the phase shift here
    phase_shift!(mvis1, pars.mu_RA + half_pix, pars.mu_DEC - half_pix)
    phase_shift!(mvis2, pars2.mu_RA + half_pix, pars2.mu_DEC - half_pix)

    # Calculate the likelihood between these two
    lnprobs1[i] = lnprob(dv, mvis1)
    lnprobs2[i] = lnprob(dv, mvis2)
    lnprobs_diff[i] = lnprobs1[i] - lnprobs2[i]
end

# Store the max and min values for each type of plot
max_real = maximum(abs, real(vis_fft1[1].VV))
max_real_resid = maximum(abs, real(vis_resid[1].VV))
max_imag = maximum(abs, imag(vis_fft1[1].VV))
max_imag_resid = maximum(abs, imag(vis_resid[1].VV))

for i=1:nchan

    max_real_temp = maximum(abs, real(vis_fft1[i].VV))
    if max_real_temp > max_real
        max_real = max_real_temp
    end

    max_real_temp = maximum(abs, real(vis_fft2[i].VV))
    if max_real_temp > max_real
        max_real = max_real_temp
    end

    max_real_temp = maximum(abs, real(vis_resid[i].VV))
    if max_real_temp > max_real_resid
        max_real_resid = max_real_temp
    end

    max_imag_temp = maximum(abs, imag(vis_fft1[i].VV))
    if max_imag_temp > max_imag
        max_imag = max_imag_temp
    end

    max_imag_temp = maximum(abs, imag(vis_fft2[i].VV))
    if max_imag_temp > max_imag
        max_imag = max_imag_temp
    end

    max_imag_temp = maximum(abs, imag(vis_resid[i].VV))
    if max_imag_temp > max_imag_resid
        max_imag_resid = max_imag_temp
    end

end

# Create the normalize instances for each
norm_img = PyPlot.matplotlib[:colors][:Normalize](0, max_img)
norm_resid = PyPlot.matplotlib[:colors][:Normalize](-max_resid, max_resid)

norm_real = PyPlot.matplotlib[:colors][:Normalize](-max_real, max_real)
norm_real_resid = PyPlot.matplotlib[:colors][:Normalize](-max_real_resid, max_real_resid)

norm_imag = PyPlot.matplotlib[:colors][:Normalize](-max_imag, max_imag)
norm_imag_resid = PyPlot.matplotlib[:colors][:Normalize](-max_imag_resid, max_imag_resid)


# Now go through and plot everything
fig, ax = plt[:subplots](nrows=10, ncols=nchan, figsize=(2.0 * nchan, 24.))

# Image needs to be flipped along RA dimension
ext_img = (skimg1.ra[end], skimg1.ra[1], skimg1.dec[1], skimg1.dec[end])

# extent in Fourier domain
ext_vis = (vis_fft1[1].uu[1], vis_fft1[1].uu[end], vis_fft1[1].vv[1], vis_fft1[1].vv[end])

cm = plt[:get_cmap]("PuBu")
# Reversed color map
cm_resid = plt[:get_cmap]("bwr_r")

for i=1:nchan
    # Images
    ax[1, i][:imshow](flipdim(skimg1.data[:,:,i], 2), interpolation="none", origin="lower", cmap=cm, extent=ext_img, norm=norm_img)
    ax[2, i][:imshow](flipdim(skimg2.data[:,:,i], 2), interpolation="none", origin="lower", cmap=cm, extent=ext_img, norm=norm_img)
    im = ax[3, i][:imshow](flipdim(skresid.data[:,:,i], 2), interpolation="none", origin="lower", cmap=cm_resid, extent=ext_img, norm=norm_resid)

    if i==nchan
        cax = fig[:add_axes]([0.95, 0.7, 0.01, 0.07])
        cbar = fig[:colorbar](im, cax=cax)
        # cbar.ax.tick_params(labelsize=6)
    end

    ## Fourier Domain

    # Real
    ax[4, i][:imshow](real(vis_fft1[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_real)
    ax[5, i][:imshow](real(vis_fft2[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_real)
    im = ax[6, i][:imshow](real(vis_resid[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_real_resid)

    if i==nchan
        cax = fig[:add_axes]([0.95, 0.43, 0.01, 0.07])
        cbar = fig[:colorbar](im, cax=cax)
        # cbar.ax.tick_params(labelsize=6)
    end

    # Imag
    ax[7, i][:imshow](imag(vis_fft1[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_imag)
    ax[8, i][:imshow](imag(vis_fft2[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_imag)
    im = ax[9, i][:imshow](imag(vis_resid[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_imag_resid)

    if i==nchan
        cax = fig[:add_axes]([0.95, 0.15, 0.01, 0.07])
        cbar = fig[:colorbar](im, cax=cax)
        # cbar.ax.tick_params(labelsize=6)
    end

    ax[9, i][:annotate](@sprintf("%.2f", lnprobs_diff[i]), (0.1, 0.2), xycoords="axes fraction", size=8)

    dv = dvarr[i]
    ax[10,i][:plot](dv.uu, dv.vv, "k.", ms=0.2)
    ax[10,i][:set_xlim](vis_fft1[1].uu[1], vis_fft1[1].uu[end])
    ax[10,i][:set_ylim](vis_fft1[1].vv[1], vis_fft1[1].vv[end])


    # Remove tick labels from all subsequent columns.
    if i==1
        # Do nothing
    else
        for j=1:10
            ax[j, i][:xaxis][:set_ticklabels]([])
            ax[j, i][:yaxis][:set_ticklabels]([])
        end
    end

    # for j=1:10
    #     if (j==1 && i==1) || (j==4 && i==1) || (j==10 && i==1)
    #         # Do nothing
    #     else
    #         ax[j, i][:xaxis][:set_ticklabels]([])
    #         ax[j, i][:yaxis][:set_ticklabels]([])
    #     end
    # end
end

for j=1:3
    ax[j, 1][:set_xlabel](L"$\Delta \alpha$ ['']")
    ax[j, 1][:set_ylabel](L"$\Delta \delta$ ['']")
end

for j=4:10
    ax[j, 1][:set_ylabel](L"UU $[k\lambda]$")
    ax[j, 1][:set_xlabel](L"VV $[k\lambda]$")
end

name1 = config["name"]
name2 = config2["name"]

ax[1, 1][:annotate]("image $name1", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[2, 1][:annotate]("image $name2", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[3, 1][:annotate]("image $name1 - $name2", (0.1, 0.8), xycoords="axes fraction", size=8)

ax[4, 1][:annotate]("real $name1", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[5, 1][:annotate]("real $name2", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[6, 1][:annotate]("real $name1 - $name2", (0.1, 0.8), xycoords="axes fraction", size=8)

ax[7, 1][:annotate]("imag $name1", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[8, 1][:annotate]("imag $name2", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[9, 1][:annotate]("imag $name1 - $name2", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[9, 1][:annotate]("lnprob diff", (0.1, 0.4), xycoords="axes fraction", size=8)
ax[10, 1][:annotate]("baselines", (0.1, 0.8), xycoords="axes fraction", size=8)

fig[:subplots_adjust](hspace=0.3, wspace=0.001, top=0.95, bottom=0.05, left=0.1, right=0.9)
fig[:savefig]("chmaps_compare.png", dpi=300)
