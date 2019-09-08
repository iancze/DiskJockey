# Given some model parameters, synthesize and plot the channel maps and
# integrated spectrum
using Pkg; Pkg.activate("DiskJockey")

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


# Plot the channel maps using sky convention
function plot_chmaps(img::image.SkyImage)

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    # Image needs to be flipped along RA dimension
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    # CO 2-1 rest frame
    lam0 = cc/230.538e9 * 1e4 # [microns]
    nlam = length(img.lams)

    # Specify velocities

    # convert wavelengths to velocities
    vels = c_kms * (img.lams .- lam0)/lam0

    # Figure out how many plots we'll have.
    ncols = 5
    nrows = 1

    fig, ax = plt.subplots(nrows=1, ncols=ncols, figsize=(5, 1.5))


    for col=1:ncols
        iframe = col

        if iframe > nlam
            # Stop if we run out of channels
            # Plot a blank square
            ax[col][:imshow](zeros((im_ny, im_nx)), cmap=plt.get_cmap("PuBu"), vmin=0, vmax=20, extent=ext, origin="lower")

        else
            #Flip the frame for Sky convention
            frame = fliplr(img.data[:,:,iframe])
            frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
            lframe = log10(frame)
            max = maximum(lframe)
            ix,iy = ind2sub(size(lframe), indmax(lframe))
            ax[col][:imshow](lframe, extent=ext, vmin=max - 6, vmax=max, interpolation="none", origin="lower", cmap=plt.get_cmap("PuBu"))
            levels = linspace(max - 0.8, max, 8)
        end

        if col != 1
            ax[col][:xaxis][:set_ticklabels]([])
            ax[col][:yaxis][:set_ticklabels]([])
            ax[col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)
        else
            ax[col][:set_xlabel](L"$\Delta \alpha$ ('')")
            ax[col][:set_ylabel](L"$\Delta \delta$ ('')")
            ax[col][:annotate](@sprintf("%.1f km/s", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)
        end

    end


    fig[:subplots_adjust](hspace=0.06, wspace=0.01, top=0.98, bottom=0.13, left=0.05, right=0.95)

    plt.savefig("plots/channel_maps_sky.svg")

end

# Plot the spatially-integrated spectrum
function plot_spectrum(img::image.SkyImage)

    fig = plt.figure(figsize=(3.4,2.5))
    ax = fig[:add_subplot](111)

    spec = imToSpec(img)

    # CO 2-1 rest frame
    lam0 = cc/230.538e9 * 1e4 # [microns]

    # convert wavelengths to velocities
    vels = c_kms * (img.lams .- lam0)/lam0 # [km/s]

    ax[:plot](vels, spec[:,2], ls="steps-mid")
    # ax[:plot](vels, reverse(spec[:,2]), ls="steps-mid")
    ax[:set_ylabel](L"$f_\nu$ [Jy]")
    ax[:set_xlabel](L"$v$ [km/s]")
    ax[:set_xlim](minimum(vels), maximum(vels))

    fig[:subplots_adjust](left=0.2, bottom=0.2, right=0.8)

    plt.savefig("plots/spectrum_talk.svg")
end

# read the wavelengths for all 23 channels
fid = h5open(config["data_file"], "r")
lams = read(fid["lams"]) # [μm]
close(fid)


pp = config["parameters"]
params = ["M_star", "r_c", "T_10", "q", "gamma", "logM_CO", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
starting_param = Array{Float64}(nparam)

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
# beta = vel/c_kms # relativistic Doppler formula
# shift_lams =  lams .* sqrt((1. - beta) / (1. + beta)) # [microns]

grd = config["grid"]
grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)

# Comment out afterward
# For channel maps
# nchan = 5
# global const vels = linspace(-4., 4, nchan) # [km/s]

# For spectrum
nchan = 60
global const vels = linspace(-16., 16, nchan) # [km/s]
# # CO 2-1 rest frame
lam0 = cc/230.538e9 * 1e4 # [microns]
# # convert velocities to wavelengths
shift_lams = lam0 * (vels/c_kms + 1)
#
# write_grid("", grid)
# write_model(pars, "", grid)
# write_lambda(shift_lams, "")
#
# if !parsed_args["norad"]
#     run(`radmc3d image incl $incl posang $PA npix $npix loadlambda`)
# end

im = imread()

# plot_chmaps(im)

skim = imToSky(im, pars.dpc)

# plot_chmaps(skim)
plot_spectrum(skim)
