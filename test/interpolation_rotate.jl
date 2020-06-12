# This file is designed to test the visibility interpolation algorithm.
import PyPlot
import PyPlot.plt
using LaTeXStrings

using DiskJockey.constants
using DiskJockey.gridding
using DiskJockey.visibilities
using DiskJockey.image
using DiskJockey.gauss

# Because of the flipped nature of the sky (but not flipped nature of the UV plane)
# there are some tricky conventions about how to pack the array.



# Test to see if we get the convolutional interpolation correct by using a 2D
# Elliptical Gaussian.
# Realistic Gaussian will have scale dimensions (fatter in x direction)
const mu_RA = 2.0 # [arcsec]
const mu_DEC = -0.5 # [arcsec]
const s_x = 1.0 # [arcsec]
const s_y = 2.0 # [arcsec]
const rho = 0.0
const theta = 20.0 # deg
const p0 = [mu_RA, mu_DEC, s_x, s_y, rho] # [arcsec]
const p_center = [0.0, 0.0, s_x, s_y, rho] # [arcsec]

# If the synthesized beam is 0.3", to oversample this by a factor of 10 would require something like 512x512 pixels
nx = 128
ny = 128

# full span of the image
ra = fftspace(10., nx) # [arcsec]
dec = fftspace(10., ny) # [arcsec]

# convert ra and dec in [arcsec] to radians, and then take the sin to convert to ll, mm
ll = sin.(ra * arcsec) # direction cosines
mm = sin.(dec * arcsec)

# The natural, rotated and shifted Gaussian image
img = imageGauss(ll, mm, p0, 1, theta=theta)

# A centered Gaussian image which we will then rotate and shift via the phase theorem
img_plain1 = imageGauss(ll, mm, p_center, 1) #theta = 0

lam0 = lam0s["12CO2-1"] # [microns] 12CO 2-1

skim = SkyImage(img, ra, dec, lam0)

skim_plain1 = SkyImage(img_plain1, ra, dec, lam0)

# Apply a centered correction function on the rotated and offset image and Fourier transform the image
corrfun!(skim, 0., 0.)
shift_fft = transform(skim)

# Now, create images that we will use as test cases for doing the interpolation before or after
# sampling

# Apply a centered correction function on the centered image (to rotate and phase shift downsampled
# visibilities later)
corrfun!(skim_plain1, 0, 0)
vis_fft_center = transform(skim_plain1)
# Apply Rotation
# Apply Phase shift


# Return a normalized instance that is symmetric about 0
function scale(data)
    s = maximum(abs.(data))
    return norm = plt[:Normalize](vmin=-s, vmax=s, clip=false)
end


# Basic plot of the shifted image

# Because the sky convention is different than the way the SkyImage is stored,
# we need to flip the array for plotting
fig, ax = plt[:subplots](nrows=2, figsize=(5, 8))

ext = (skim.ra[end], skim.ra[1], skim.dec[1], skim.dec[end])
ax[1][:imshow](flipdim(skim.data[:,:,1], 2), interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext)
ax[1][:contour](flipdim(skim.data[:,:,1], 2), origin="lower", extent=ext)
ax[1][:set_title]("Sky Projection")
ax[1][:set_xlabel](L"$\alpha$ [arcsec]")
ax[1][:set_ylabel](L"$\delta$ [arcsec]")

ext = (ll[1], ll[end], mm[1], mm[end])
ax[2][:imshow](skim.data[:,:,1], interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext)
ax[2][:contour](skim.data[:,:,1], origin="lower", extent=ext)
ax[2][:set_title]("Raw Array")
ax[2][:set_xlabel](L"$ll$")
ax[2][:set_ylabel](L"$mm$")

fig[:subplots_adjust](left=0.15, right=0.85, hspace=0.25)
plt[:savefig]("gaussian_img_array_rot.png")


# Basic plot of the *centered* Gaussian image.
fig, ax = plt[:subplots](nrows=1, figsize=(5, 5))
# Real, analytic Gaussian
ext = (skim.ra[end], skim.ra[1], skim.dec[1], skim.dec[end])
aximg = ax[:imshow](flipdim(skim_plain1.data[:,:,1], 2), interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext) #, norm = scale(img))
ax[:set_title]("image")
ax[:set_xlabel](L"$\alpha$ [arcsec]")
ax[:set_ylabel](L"$\delta$ [arcsec]")
#[left, bottom, width, height]
cax = fig[:add_axes]([0.84, 0.25, 0.03, 0.45])
cb = fig[:colorbar](aximg, cax=cax)

fig[:subplots_adjust](left=0.15, right=0.85, hspace=0.25)
plt[:savefig]("gaussian_img_center_rot.png")


# It's nice to see these plots and what they look like, but it doesn't really make sense
# to compare them, unless we are comparing *interpolated* visibilities to the analytic formula.
# Just comparing the image is WRONG.

# THIS (following) IS THE ONLY THING WE SHOULD CARE ABOUT.

# Next, choose some (u,v) points within the bounds and see how the
# interpolated values compare to what the FTGauss would return.

nu = length(shift_fft.uu)
zer = zeros(nu)

function plot_1d(analytic, approx, fname)
    fig, ax = plt[:subplots](nrows=2, figsize=(5, 8))
    ax[1][:plot](shift_fft.uu, zer, ".k", label="Grid spacing")
    ax[1][:plot](uu, real(approx), "ob", label="Interp")
    ax[1][:plot](uu, real(analytic), ".r", label="Analytic")
    # ax[1][:plot](vis_fft.uu, real(analytic_u), "or", label="Analytic")
    ax[1][:set_xlim](-100, 100)
    ax[1][:set_title]("Real")
    ax[1][:set_xlabel](L"u [k $\lambda$]")
    ax[1][:legend]()

    ax[2][:plot](shift_fft.uu, zer, ".k", label="Grid spacing")
    ax[2][:plot](uu, imag(approx), "ob", label="Interp")
    ax[2][:plot](uu, imag(analytic), ".r", label="Analytic")
    # ax[2][:plot](vis_fft.uu, imag(analytic_u), "or", label="Analytic")
    ax[2][:set_xlim](-100, 100)
    ax[2][:set_title]("Imag")
    ax[2][:set_xlabel](L"u [k $\lambda$]")

    plt[:savefig](fname)
end

n = 100
uu = linspace(-100, 100, n) # [kλ]
approx = Array{Complex128}(n)
analytic = Array{Complex128}(n)

# First, let's see how the interpolated points, with *no shift*, correspond to the analytic form
# with no shift.

for i=1:n
    u = uu[i]
    v = 0.0
    #println("Interpolating at $u, $v")
    approx[i] = interpolate_uv(u, v, vis_fft_center)
    analytic[i] = FTGauss(u, v, [0.0, 0.0, s_x, s_y, rho], 1)
end

plot_1d(analytic, approx, "interpolation_norot_noshift.png")

# Next, let's see how the interpolated points, then rotated, correspond to the analytic form with rotation.

for i=1:n
    u = uu[i]
    v = 0.0

    up = cosd(theta) * u - sind(theta) * v
    vp = sind(theta) * u - cosd(theta) * v

    approx[i] = interpolate_uv(up, vp, vis_fft_center)
    analytic[i] = FTGauss(u, v, [0.0, 0.0, s_x, s_y, rho], 1, theta=theta)
end

plot_1d(analytic, approx, "interpolation_rot_noshift.png")

# Finally, let's see how the interpolated points, rotated, then phase-shifted, compare to the analytic form with rotation and phase shift.
mu = Float64[mu_RA, mu_DEC] * arcsec # [radians]

for i=1:n
    u = uu[i]
    v = 0.0

    up = cosd(theta) * u - sind(theta) * v
    vp = sind(theta) * u - cosd(theta) * v

    # Interpolate points at rotated u', v' points, but then shift according to
    # direction and amplitude of mu offset using u,v coordinates.

    # Now shift the interpolated points according to the phase theorem
    R = Float64[u, v] * 1e3 #[λ]
    # Not actually in polar phase form
    shift = exp(-2pi * 1.0im * (R' * mu)[1])

    #println("Interpolating at $u, $v")
    approx[i] = shift * interpolate_uv(up, vp, vis_fft_center)
    analytic[i] = FTGauss(u, v, p0, 1, theta=theta)

end

plot_1d(analytic, approx, "interpolation_rot_shift.png")
