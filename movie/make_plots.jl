push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

#Make the movie!

using constants
using image

import PyPlot.plt
using LaTeXStrings

# Movie of V4046Sgr, showing 7 representative channels

global const nchan = 7
global const vels = linspace(-1.5, 1.5, nchan) # [km/s]
# CO 2-1 rest frame
lam0 = cc/230.538e9 * 1e4 # [microns]
# convert velocities to wavelengths
lams = lam0 * (vels/c_kms + 1)

# Plot the channel maps using sky convention
function plot_chmaps(img::image.SkyImage, id::Int)

    fig, ax = plt.subplots(ncols=nchan, figsize=(8, 1.6))
    # Image needs to be flipped along RA dimension
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    for col=1:nchan

        # Set all labels but the leftmost to blank
        if col != 1
            ax[col][:xaxis][:set_ticklabels]([])
            ax[col][:yaxis][:set_ticklabels]([])
        else
            ax[col][:set_xlabel](L"$\Delta \alpha$ ('')")
            ax[col][:set_ylabel](L"$\Delta \delta$ ('')")
            ax[col][:tick_params](axis="both", which="major", labelsize=8)
        end

        #Flip the frame for Sky convention
        frame = fliplr(img.data[:,:,col])
        frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
        lframe = log10(frame)
        max = maximum(lframe)
        ax[col][:imshow](lframe, extent=ext, vmin=max - 6, vmax=max, interpolation="none", origin="lower", cmap=plt.get_cmap("PuBu"))
        ax[col][:annotate](@sprintf("%.1f", vels[col]), (0.1, 0.8), xycoords="axes fraction", size=8)

    end

    fig[:subplots_adjust](wspace=0.08, top=0.95, bottom=0.26, left=0.1, right=0.9)

    plt.savefig(@sprintf("%04d.png", id))
    plt.close(fig)

end

# Determine which image*.out files are in this directory.
imgfunc = x -> contains(x, "image") && contains(x, ".out")
imglist = filter(imgfunc, readdir())



# Go through each one, read them into an array.

# Determine the maximum and minimum of this array


# plot frames
im = imread(file=fname)
skim = imToSky(im, pars.dpc)

plot_chmaps(skim, id)
