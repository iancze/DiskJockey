# Given some model parameters, synthesize and plot the channel maps and
# integrated spectrum

using constants
using image
using model
using HDF5

import PyPlot.plt
using LaTeXStrings

function plot_chmaps(img::image.RawImage)

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    # CO 2-1 rest frame
    lam0 = cc/230.538e9 * 1e4 # [microns]
    nlam = length(img.lams)
    println("Plotting ", nlam, " channels")

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

    plt.savefig("plots/channel_maps_raw.png")

end

function plot_chmaps(img::image.SkyImage)

    (im_ny, im_nx) = size(img.data)[1:2] # y and x dimensions of the image

    # Image will be flipped
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    # CO 2-1 rest frame
    lam0 = cc/230.538e9 * 1e4 # [microns]
    nlam = length(img.lams)
    println("Plotting ", nlam, " channels")

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

            #Flip the frame for Sky convention
            frame = fliplr(img.data[:,:,iframe])
            frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
            max = maximum(log10(frame))
            ax[row, col][:imshow](log10(frame), extent=ext, vmin=max - 6, vmax=max, interpolation="none", origin="lower", cmap=plt.get_cmap("PuBu"))
            levels = linspace(max - 0.8, max, 5)
            ax[row, col][:contour](log10(frame), origin="lower", colors="k", levels=levels, extent=ext, linestyles="solid", linewidths=0.2)

            if col != 1 || row != 2
                ax[row, col][:xaxis][:set_ticklabels]([])
                ax[row, col][:yaxis][:set_ticklabels]([])
            else
                ax[row, col][:set_xlabel](L"$\Delta \alpha$ ('')")
                ax[row, col][:set_ylabel](L"$\Delta \delta$ ('')")
            end

            ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)

        end
    end

    fig[:subplots_adjust](hspace=0.01, wspace=0.05, top=0.9, bottom=0.1, left=0.05, right=0.95)

    plt.savefig("plots/channel_maps_sky.png")

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
fid = h5open("data/V4046Sgr.hdf5", "r")
lams = read(fid["lams"]) # [μm]
close(fid)

#From Rosenfeld et al. 2012, Table 1
M_star = 1.75 # [M_sun] stellar mass
r_c =  45. # [AU] characteristic radius
T_10 =  115. # [K] temperature at 10 AU
q = 0.63 # temperature gradient exponent
gamma = 1.0 # surface temperature gradient exponent
M_CO = 0.933 # [M_earth] disk mass of CO
ksi = 0.14 # [km/s] microturbulence
dpc = 73.0
incl = 135. # [degrees] inclination
vel = -31.18 # [km/s]
PA = -17.
mu_RA = 0.0 # [arcsec]
mu_DEC = 0.0 # [arcsec]

pars = Parameters(M_star, r_c, T_10, q, gamma, M_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

vel = pars.vel # [km/s]

# We are using the Pietu convention, where inclination ranges from +90 to -90 degrees.
# +90 means face on, angular momentum vector pointing at observer.
# 0 means edge on
# -90 means face on, angular momentum vector pointing away from observer.
# RADMC conventions define
# 0 as face on, angular momentum towards observer.
# 90 as edge on
# 180 as face on, angular momentum away from observer.
# Therefore, we convert from Pietu convention (pars.incl) to RADMC convetion (incl)
incl = pars.incl # [deg]

# We also adopt the Pietu convention for position angle, which defines position angle
# by the angular momentum vector. No conversion for RADMC is necessary.
PA = pars.PA # [deg] Position angle runs counter clockwise, due to looking at sky.
npix = 128 # number of pixels, can alternatively specify x and y separately

# Doppler shift the dataset wavelength to rest-frame wavelength
beta = vel/c_kms # relativistic Doppler formula
shift_lams =  lams .* sqrt((1. - beta) / (1. + beta)) # [microns]

write_grid("")
write_model(pars, "")
write_lambda(shift_lams)

run(`radmc3d image incl $incl posang $PA npix $npix loadlambda`)

im = imread()

plot_chmaps(im)

skim = imToSky(im, pars.dpc)

plot_chmaps(skim)

spec = imToSpec(skim)
plot_spectrum(spec)
