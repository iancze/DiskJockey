module debug

# Home to all of the routines that are useful for debugging purposes

using ..constants
using ..model
using ..image

import PyPlot
import PyPlot.plt
using LaTeXStrings

export format_grid, synthesize_image, plot_chmaps, plot_slices

"""Format grid values to a standard string"""
function format_grid(grid::Grid, npix::Int)
    #nr, ntheta, r_in, r_out, npix
    r_in = @sprintf("%.1f", grid.rs[1]/AU)
    r_out = @sprintf("%.1f", grid.rs[end]/AU)
    return "$(grid.nr)_$(grid.ntheta)_$(r_in)_$(r_out)_$npix"
end

"""Synthesize images for these parameters, grid, and npix values. dlams provides the frequencies at which to synthesize images. Assumes that you have run DJInitialize.jl inside of homedir. When done, copies the image to outfile."""
function synthesize_image(dlams::Array{Float64, 1}, pars::Parameters, grid::Grid, npix::Int, outfile::AbstractString, species::AbstractString)

    # Assumes that all of the necessary files are already in the current working directory.
    write_grid("", grid)

    # Compute parameter file using model.jl, write to disk in current directory
    write_model(pars, "", grid, species)

    # Doppler shift the dataset wavelengths to rest-frame wavelength
    nchan = length(dlams)
    beta = pars.vel/c_kms # relativistic Doppler formula
    lams = Array(Float64, nchan)
    for i=1:nchan
        lams[i] =  dlams[i] * sqrt((1. - beta) / (1. + beta)) # [microns]
    end

    write_lambda(lams, "") # write into current directory

    # Run RADMC-3D, redirect output to /dev/null
    run(pipeline(`radmc3d image incl $(pars.incl) posang $(pars.PA) npix $npix loadlambda`, DevNull))

    # Copy the output image to the new filename
    run(`mv image.out image_$(outfile).out`)
end


# Plot the channel maps using sky convention
function plot_chmaps(img::image.SkyImage, lam0, fname="channel_maps_sky.png")

    vels = c_kms * (img.lams .- lam0)/lam0

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    # Image needs to be flipped along RA dimension
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    nlam = length(img.lams)

    # Figure out how many plots we'll have.
    ncols = 8
    nrows = ceil(Int, nlam/ncols)

    fig, ax = plt[:subplots](nrows=nrows, ncols=ncols, figsize=(12, 1.5 * nrows))

    # max = maximum(log10(img.data))

    vvmax = maxabs(img.data)
    norm = PyPlot.matplotlib[:colors][:Normalize](0, vvmax)

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
                ax[row, col][:imshow](zeros((im_ny, im_nx)), cmap=plt[:get_cmap]("PuBu"), vmin=0, vmax=20, extent=ext, origin="lower")

            else
                #Flip the frame for Sky convention
                frame = flipdim(img.data[:,:,iframe], 2)

                # frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
                # lframe = log10(frame)
                #
                # ax[row, col][:imshow](lframe, extent=ext, interpolation="none", vmin=max - 6, vmax=max, origin="lower", cmap=plt[:get_cmap]("PuBu"))

                ax[row, col][:imshow](frame, extent=ext, interpolation="none", origin="lower", cmap=plt[:get_cmap]("PuBu"), norm=norm)


                ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)
            end
        end
    end
    fig[:subplots_adjust](hspace=0.06, wspace=0.01, top=0.9, bottom=0.1, left=0.05, right=0.95)
    plt[:savefig](fname)
end

function plot_slices(img::image.SkyImage, index, fname="")
end

end #module
