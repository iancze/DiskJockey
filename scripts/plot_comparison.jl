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
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

# Load two separate high resolution images and compare them, computing the difference and the difference in lnprob.

using JudithExcalibur.constants
using JudithExcalibur.image
using JudithExcalibur.model
using JudithExcalibur.visibilities
using JudithExcalibur.gridding

import PyPlot.plt
using LaTeXStrings

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]

pars = convert_dict(config["parameters"], config["model"])

img1 = imread(parsed_args["image1"])
img2 = imread(parsed_args["image2"])

nchan = length(img1.lams)
# Make sure we have the same number of channels to compare images.
@assert nchan == length(img2.lams) "Images have a different number of channels!"

# Convert raw images to the appropriate distance
skimg1 = imToSky(img1, pars.dpc)
skimg2 = imToSky(img2, pars.dpc)
skresid = skimg1 - skimg2

# Store maximums for plotting
max_img = maxabs(skimg1.data)
if maxabs(skimg2.data) > max_img
    max_img = maxabs(skimg2.data)
end
max_resid = maxabs(skresid.data)

# Apply the gridding correction function before doing the FFT
# No shift needed, since we will shift the resampled visibilities
# The skimg1, skimg2 are for plotting.
skimg01 = corrfun(skimg1)
skimg02 = corrfun(skimg2)

# To store the FFT for each channel
vis_fft1 = Array(FullModelVis, nchan)
vis_fft2 = Array(FullModelVis, nchan)
vis_resid = Array(FullModelVis, nchan)

# To store the computed lnprob for each channel
lnprobs1 = Array(Float64, nchan)
lnprobs2 = Array(Float64, nchan)
lnprobs_diff = Array(Float64, nchan)

# Load the data file so we can use it to compute lnprob
dvarr = DataVis(config["data_file"])

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
    phase_shift!(mvis1, pars.mu_RA, pars.mu_DEC)
    phase_shift!(mvis2, pars.mu_RA, pars.mu_DEC)

    # Calculate the likelihood between these two
    lnprobs1[i] = lnprob(dv, mvis1)
    lnprobs2[i] = lnprob(dv, mvis2)
    lnprobs_diff[i] = lnprobs1[i] - lnprobs2[i]
end

# Store the max and min values for each type of plot
max_real = maxabs(real(vis_fft1[1].VV))
max_real_resid = maxabs(real(vis_resid[1].VV))
max_imag = maxabs(imag(vis_fft1[1].VV))
max_imag_resid = maxabs(imag(vis_resid[1].VV))

for i=1:nchan

    max_real_temp = maxabs(real(vis_fft1[i].VV))
    if max_real_temp > max_real
        max_real = max_real_temp
    end

    max_real_temp = maxabs(real(vis_fft2[i].VV))
    if max_real_temp > max_real
        max_real = max_real_temp
    end

    max_real_temp = maxabs(real(vis_resid[i].VV))
    if max_real_temp > max_real_resid
        max_real_resid = max_real_temp
    end

    max_imag_temp = maxabs(imag(vis_fft1[i].VV))
    if max_imag_temp > max_imag
        max_imag = max_imag_temp
    end

    max_imag_temp = maxabs(imag(vis_fft2[i].VV))
    if max_imag_temp > max_imag
        max_imag = max_imag_temp
    end

    max_imag_temp = maxabs(imag(vis_resid[i].VV))
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
fig, ax = plt[:subplots](nrows=10, ncols=nchan, figsize=(1.2 * nchan, 16.))

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
    ax[3, i][:imshow](flipdim(skresid.data[:,:,i], 2), interpolation="none", origin="lower", cmap=cm_resid, extent=ext_img, norm=norm_resid)

    ## Fourier Domain

    # Real
    ax[4, i][:imshow](real(vis_fft1[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_real)
    ax[5, i][:imshow](real(vis_fft2[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_real)
    ax[6, i][:imshow](real(vis_resid[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_real_resid)

    # Imag
    ax[7, i][:imshow](imag(vis_fft1[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_imag)
    ax[8, i][:imshow](imag(vis_fft2[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_imag)
    ax[9, i][:imshow](imag(vis_resid[i].VV), extent=ext_vis, interpolation="none", origin="lower", cmap=cm_resid, norm=norm_imag_resid)

    ax[9, i][:annotate](@sprintf("%.2f", lnprobs_diff[i]), (0.1, 0.2), xycoords="axes fraction", size=8)

    dv = dvarr[i]
    ax[10,i][:plot](dv.uu, dv.vv, "k.", ms=0.2)
    ax[10,i][:set_xlim](vis_fft1[1].uu[1], vis_fft1[1].uu[end])
    ax[10,i][:set_ylim](vis_fft1[1].vv[1], vis_fft1[1].vv[end])


    # For all rows
    for j=1:10
        if (j==1 && i==1) || (j==4 && i==1) || (j==10 && i==1)
            # Do nothing
        else
            ax[j, i][:xaxis][:set_ticklabels]([])
            ax[j, i][:yaxis][:set_ticklabels]([])
        end
    end
end

ax[1, 1][:set_ylabel](L"$\Delta \delta$ ['']")
ax[4, 1][:set_ylabel](L"UU $[k\lambda]$")

ax[1, 1][:annotate]("image 1", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[2, 1][:annotate]("image 2", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[3, 1][:annotate]("resid", (0.1, 0.8), xycoords="axes fraction", size=8)

ax[4, 1][:annotate]("real 1", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[5, 1][:annotate]("real 2", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[6, 1][:annotate]("real resid", (0.1, 0.8), xycoords="axes fraction", size=8)

ax[7, 1][:annotate]("imag 1", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[8, 1][:annotate]("imag 2", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[9, 1][:annotate]("imag resid", (0.1, 0.8), xycoords="axes fraction", size=8)
ax[9, 1][:annotate]("lnprob diff", (0.1, 0.4), xycoords="axes fraction", size=8)
ax[10, 1][:annotate]("baselines", (0.1, 0.8), xycoords="axes fraction", size=8)

fig[:subplots_adjust](hspace=0.01, wspace=0.001, top=0.95, bottom=0.05, left=0.05, right=0.99)
fig[:savefig]("chmaps_compare.png", dpi=300)
