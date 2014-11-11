# Plot the image and spectrum

module plot_image

using image
using constants

import PyPlot.plt
using LaTeXStrings


function plot_chmaps(img::image.SkyImage)

    (im_nx, im_ny) = size(img.data)[1:2] #x and y dimensions of the image

    ext = (img.ra[1], img.ra[end], img.dec[1], img.dec[end])

    lam0 = cc/230.538e9 * 1e6 # [microns]
    nlam = length(img.lams)
    vels = linspace(-4.4, 4.4, nlam) # [km/s]

    fig, ax = plt.subplots(nrows=2, ncols=12, figsize=(12, 2.8))

    #Plot a blank img in the last frame, moment map will go here eventually
    ax[2, 12][:imshow](zeros((im_nx, im_ny)), cmap=plt.get_cmap("Greys"), vmin=0, vmax=20)
    ax[2, 12][:xaxis][:set_ticklabels]([])
    ax[2, 12][:yaxis][:set_ticklabels]([])

    # Loop through all of the different channels and plot them
    #for iframe=1:nlam

    for row=1:2
        for col=1:12
            iframe = col + (row - 1) * 12

            if iframe > nlam
                break
            end

            if col != 1 || row != 2
                ax[row, col][:xaxis][:set_ticklabels]([])
                ax[row, col][:yaxis][:set_ticklabels]([])
            else
                ax[row, col][:set_xlabel](L"$\Delta \alpha$ ('')")
                ax[row, col][:set_ylabel](L"$\Delta \delta$ ('')")
            end

            ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)

            frame = img.data[:,:,iframe]
            frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
            max = maximum(log10(frame))
            ax[row, col][:imshow](log10(frame), extent=ext, vmin=max - 7, vmax=max, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"))
        end
    end

    fig[:subplots_adjust](hspace=0.01, wspace=0.05, top=0.9, bottom=0.1, left=0.05, right=0.95)

    plt.savefig("plots/channel_maps.png")

end

# Plot the spatially-integrated spectrum
function plot_spectrum(spec::Array{Float64, 2})

    fig = plt.figure()
    ax = fig[:add_subplot](111)

    lam0 = cc/230.538e9 * 1e6 # [microns]
    nvels = length(spec[:,1])
    vels = linspace(-4.4, 4.4, nvels) # [km/s]

    ax[:plot](vels, spec[:,2], ls="steps-mid")
    ax[:set_ylabel](L"$f_\nu$ [Jy]")
    ax[:set_xlabel](L"$v$ [$\textrm{km/s}$]")


    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("plots/spectrum.png")
end

end
