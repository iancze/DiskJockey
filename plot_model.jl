# Given some model parameters, synthesize and plot the channel maps and
# integrated spectrum

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    # "--opt1"
    # help = "an option with an argument"
    # default = 0
    "--norad"
    help = "Use the image alread here."
    action = :store_true
    "config"
    help = "a YAML configuration file"
    required = true
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using constants
using image
using model
using HDF5

import PyPlot.plt
using LaTeXStrings

# Plot the raw channel maps directly from RADMC
function plot_chmaps(img::image.RawImage)

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    # CO 2-1 rest frame
    lam0 = cc/230.538e9 * 1e4 # [microns]
    nlam = length(img.lams)

    # convert wavelengths to velocities
    vels = c_kms * (img.lams .- lam0)/lam0

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

    plt.savefig("plots/channel_maps_raw_image.png")

end

# Plot the channel maps using sky convention
function plot_chmaps(img::image.SkyImage)

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    # Image needs to be flipped along RA dimension
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    # CO 2-1 rest frame
    lam0 = cc/230.538e9 * 1e4 # [microns]
    nlam = length(img.lams)

    # convert wavelengths to velocities
    vels = c_kms * (img.lams .- lam0)/lam0

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
                max = maximum(log10(frame))
                ax[row, col][:imshow](log10(frame), extent=ext, vmin=max - 6, vmax=max, interpolation="none", origin="lower", cmap=plt.get_cmap("PuBu"))
                levels = linspace(max - 0.8, max, 5)
                ax[row, col][:contour](log10(frame), origin="lower", colors="k", levels=levels, extent=ext, linestyles="solid", linewidths=0.2)



                ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)
            end

        end
    end

    fig[:subplots_adjust](hspace=0.06, wspace=0.01, top=0.9, bottom=0.1, left=0.05, right=0.95)

    plt.savefig("plots/channel_maps_sky.png")

end

# Plot the raw array for the Sky image
function plot_chmaps_data(img::image.SkyImage)

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    ll = sin(img.ra .* arcsec)
    mm = sin(img.dec .* arcsec)
    ext = (ll[1], ll[end], mm[1], mm[end])

    # CO 2-1 rest frame
    lam0 = cc/230.538e9 * 1e4 # [microns]
    nlam = length(img.lams)

    # convert wavelengths to velocities
    vels = c_kms * (img.lams .- lam0)/lam0

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

    plt.savefig("plots/channel_maps_sky_raw.png")

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

# read the wavelengths for all 23 channels
fid = h5open(config["data_file"], "r")
lams = read(fid["lams"]) # [μm]
close(fid)

pp = config["parameters"]
params = ["M_star", "r_c", "T_10", "q", "gamma", "logM_CO", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
starting_param = Array(Float64, nparam)

for i=1:nparam
    starting_param[i] = pp[params[i]][1]
end

# Convert logM_CO to M_CO
starting_param[6] = 10^starting_param[6]

pars = Parameters(starting_param...)

vel = pars.vel # [km/s]
# RADMC conventions for inclination and PA
incl = pars.incl # [deg]
PA = pars.PA # [deg] Position angle runs counter clockwise
npix = config["npix"] # number of pixels

# Doppler shift the dataset wavelength to rest-frame wavelength
beta = vel/c_kms # relativistic Doppler formula
shift_lams =  lams .* sqrt((1. - beta) / (1. + beta)) # [microns]

write_grid("")
write_model(pars, "")
write_lambda(shift_lams)

if !parsed_args["norad"]
    run(`radmc3d image incl $incl posang $PA npix $npix loadlambda`)
end

im = imread()

# plot_chmaps(im)

skim = imToSky(im, pars.dpc)

plot_chmaps(skim)
# plot_chmaps_data(skim)

# spec = imToSpec(skim)
# plot_spectrum(spec)
