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
const s_x = 1.2 # [arcsec]
const s_y = 1.0 # [arcsec]
const rho = 0.5
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

# The natural, shifted Gaussian image
img = imageGauss(ll, mm, p0, 1)

# A centered Gaussian image which we will then shift via the phase theorem
img_plain1 = imageGauss(ll, mm, p_center, 1)
img_plain2 = imageGauss(ll, mm, p_center, 1)
img_plain3 = imageGauss(ll, mm, p_center, 1)

lam0 = lam0s["12CO2-1"] # [microns] 12CO 2-1

skim = SkyImage(img, ra, dec, lam0)


skim_plain1 = SkyImage(img_plain1, ra, dec, lam0)
skim_plain2 = SkyImage(img_plain2, ra, dec, lam0)
skim_plain3 = SkyImage(img_plain3, ra, dec, lam0)

# Apply a centered correction function on the offset image and Fourier transform the image
corrfun!(skim, 0., 0.)
shift_fft = transform(skim)

# Now, create images that we will use as test cases for doing the interpolation before or after
# sampling

# Apply a centered correction function on the centered image (to phase shift downsampled
# visibilities later)
# This is the recommended approach.
corrfun!(skim_plain1, 0, 0)
vis_fft_center = transform(skim_plain1)

# Apply an offset correction function on the centered image (then  phase shift the full visibilities)
corrfun!(skim_plain2, mu_RA, mu_DEC)
vis_fft_shift = transform(skim_plain2)
phase_shift!(vis_fft_shift, mu_RA, mu_DEC)

# Apply an offset correction function on the centered image, but later we will phase shift the visibilities *after* they have been downsampled.
# NOTE, this DOES NOT WORK, and is simply included here for full testing.
corrfun!(skim_plain3, mu_RA, mu_DEC)
vis_fft_offset = transform(skim_plain3)


# This right here is not actually a fair comparison, because we haven't yet *interpolated* the FFTed values using the gridding convolution functions.

# Complex subtraction
# println("Maximum FFT discrepancy: ", maximum(abs(plain_fft.VV - vis_analytic)))
# println("Maximum FFT discrepancy: ", maximum(abs(shift_fft.VV - vis_analytic)))



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
plt[:savefig]("gaussian_img_array.png")


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
plt[:savefig]("gaussian_img_center.png")


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
    analytic[i] = FTGauss(u, v, p_center, 1)
end

plot_1d(analytic, approx, "interpolation.png")

# Next, let's see how the image shifted, then interpolated points compare with the
# analytic form

for i=1:n
    u = uu[i]
    v = 0.0
    #println("Interpolating at $u, $v")
    approx[i] = interpolate_uv(u, v, shift_fft)
    analytic[i] = FTGauss(u, v, p0, 1)
end

plot_1d(analytic, approx, "interpolation_raw_image_offset.png")

# Next,  let's see how the FFT'ed, then interpolated, then shifted points compare with the
# analytic form
mu = Float64[mu_RA, mu_DEC] * arcsec # [radians]

for i=1:n
    u = uu[i]
    v = 0.0

    # Now shift the interpolated points according to the phase theorem
    R = Float64[u, v] * 1e3 #[λ]
    # Not actually in polar phase form
    shift = exp(-2pi * 1.0im * (R' * mu)[1])

    approx[i] = shift * interpolate_uv(u, v, vis_fft_center)
    analytic[i] = FTGauss(u, v, p0, 1)
end

plot_1d(analytic, approx, "interpolation_shift_visibilities.png")


# Finally, let's check the way we've been doing things, which is to phase shift the full visibilities first, then downsample.

for i=1:n
    u = uu[i]
    v = 0.0
    #println("Interpolating at $u, $v")
    approx[i] = interpolate_uv(u, v, vis_fft_shift)
    analytic[i] = FTGauss(u, v, p0, 1)
end

plot_1d(analytic, approx, "interpolation_shift_image.png")

# As a final, final check that the resampling is in fact linear, let's use the centered image, offset corrfun, downsample this, then apply the phase shift at the end. This is where all my confusion started. If we were just able to do the previous two tests OK, it seems like this *shouldn't work*, because if in fact everything is linear, then the only difference is that this has the corrfun applied w/ a shift.

for i=1:n
    u = uu[i]
    v = 0.0

    # Now shift the interpolated points according to the phase theorem
    R = Float64[u, v] * 1e3 #[λ]
    # Not actually in polar phase form
    shift = exp(-2pi * 1.0im * (R' * mu)[1])

    approx[i] = shift * interpolate_uv(u, v, vis_fft_offset)
    analytic[i] = FTGauss(u, v, p0, 1)
end

plot_1d(analytic, approx, "interpolation_linear_operations.png")

# Ok, this *does not work*. What is the reason? Is the resampling operation not linear apparently?
# I don't think it's that, rather, I think that in order for the gridding correction function/ gridding convolution function pair to work properly, and properly cancel each other out, then they need to correspond to the same phase center. So, either you do both fixed to 0,0, or you do both offset to some value. Once the interpolation is done, you can phase shift as much as you want, but don't mix order of operations.

# Now, we can check things on an actual 2D image if we so please.

# Create analytic function on a smaller grid
n = 256
uu = linspace(-150, 150, n)
vv = linspace(-150, 150, n)

function plot_2d(analytic::Matrix{Complex128}, approx::Matrix{Complex128}, fname::AbstractString)

    println("Analyzing $fname")

    # Plot the reals first

    fig, ax = plt[:subplots](nrows=3, figsize=(5, 11))

    ext = (uu[1], uu[end], vv[1], vv[end])

    axan = ax[1][:imshow](real(analytic), interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext)
    ax[1][:set_title]("Analytic FT")
    ax[1][:set_xlabel](L"uu [k$\lambda$]")
    ax[1][:set_ylabel](L"vv [k$\lambda$]")

    cax = fig[:add_axes]([0.84, 0.70, 0.03, 0.25])
    cb = fig[:colorbar](axan, cax=cax)

    axfft = ax[2][:imshow](real(approx), interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext)
    ax[2][:set_title]("Interpolated Visibilites from FFT")
    ax[2][:set_xlabel](L"uu [k$\lambda$]")
    ax[2][:set_ylabel](L"vv [k$\lambda$]")

    cax = fig[:add_axes]([0.84, 0.40, 0.03, 0.25])
    cb = fig[:colorbar](axfft, cax=cax)

    diff = real(analytic - approx)

    # Measure the max value in the analytic array, just as an opportunity to set scale.
    println("Maximum analytic value ", maximum(abs, analytic))

    # Then, measure the total Mod squared difference across the whole image
    println("Total sqrt(modsquared) difference across image ", sqrt(sum(abs2, diff)))


    axdif = ax[3][:imshow](diff, interpolation="none", origin="lower", cmap=plt[:get_cmap]("bwr"), extent=ext, norm=scale(diff))
    ax[3][:set_title]("Difference")
    ax[3][:set_xlabel](L"uu [k$\lambda$]")
    ax[3][:set_ylabel](L"vv [k$\lambda$]")

    cax = fig[:add_axes]([0.84, 0.10, 0.03, 0.25])
    cb = fig[:colorbar](axdif, cax=cax)

    fig[:subplots_adjust](hspace=0.25, top=0.97, bottom=0.06)

    plt[:savefig](fname * "_real.png")


    fig, ax = plt[:subplots](nrows=3, figsize=(5, 11))

    axan = ax[1][:imshow](imag(analytic), interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext)
    ax[1][:set_title]("Analytic FT")
    ax[1][:set_xlabel](L"uu [k$\lambda$]")
    ax[1][:set_ylabel](L"vv [k$\lambda$]")

    cax = fig[:add_axes]([0.84, 0.70, 0.03, 0.25])
    cb = fig[:colorbar](axan, cax=cax)

    axfft = ax[2][:imshow](imag(approx), interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext)
    ax[2][:set_title]("Interpolated Visibilites from FFT")
    ax[2][:set_xlabel](L"uu [k$\lambda$]")
    ax[2][:set_ylabel](L"vv [k$\lambda$]")

    cax = fig[:add_axes]([0.84, 0.40, 0.03, 0.25])
    cb = fig[:colorbar](axfft, cax=cax)

    diff = imag(analytic - approx)

    axdif = ax[3][:imshow](diff, interpolation="none", origin="lower", cmap=plt[:get_cmap]("bwr"), extent=ext, norm=scale(diff))
    ax[3][:set_title]("Difference")
    ax[3][:set_xlabel](L"uu [k$\lambda$]")
    ax[3][:set_ylabel](L"vv [k$\lambda$]")

    cax = fig[:add_axes]([0.84, 0.10, 0.03, 0.25])
    cb = fig[:colorbar](axdif, cax=cax)

    fig[:subplots_adjust](hspace=0.25, top=0.97, bottom=0.06)

    plt[:savefig](fname * "_imag.png")

    println()

end

# First, let's see how the interpolated points, with *no shift*, correspond to the analytic form
# with no shift.
approx = Array{Complex128}(n, n)
analytic = FTGauss(uu, vv, p_center, 1) # can take in full arrays

for i=1:n
    for j=1:n
        approx[j,i] = interpolate_uv(uu[i], vv[j], vis_fft_center)
    end
end

plot_2d(analytic, approx, "2D_interpolation")

# Next, let's see how the image shifted, then interpolated points compare with the
# analytic form

analytic = FTGauss(uu, vv, p0, 1) # can take in full arrays

for i=1:n
    for j=1:n
        approx[j,i] = interpolate_uv(uu[i], vv[j], shift_fft)
    end
end

plot_2d(analytic, approx, "2D_interpolation_raw_image_offset")

# Next,  let's see how the FFT'ed, then interpolated, then shifted points compare with the
# analytic form

for i=1:n
    for j=1:n
        R = Float64[uu[i], vv[j]] * 1e3 #[λ]
        # Not actually in polar phase form
        shift = exp(-2pi * 1.0im * (R' * mu)[1])

        approx[j,i] = shift * interpolate_uv(uu[i], vv[j], vis_fft_center)
    end
end

plot_2d(analytic, approx, "2D_interpolation_shift_visibilities")


# Finally, let's check the way we've been doing things, which is to phase shift the full visibilities first, then downsample.


for i=1:n
    for j=1:n
        approx[j,i] = interpolate_uv(uu[i], vv[j], vis_fft_shift)
    end
end

plot_2d(analytic, approx, "2D_interpolation_shift_image")



fig, ax = plt[:subplots](nrows=2, figsize=(5, 8))

axan = ax[1][:imshow](abs.(analytic), interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext)
ax[1][:set_title]("Amplitude [Analytic]")
ax[1][:set_xlabel](L"uu [k$\lambda$]")
ax[1][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.70, 0.03, 0.25])
cb = fig[:colorbar](axan, cax=cax)

axfft = ax[2][:imshow](angle.(analytic), interpolation="none", origin="lower", cmap=plt[:get_cmap]("Greys"), extent=ext)
ax[2][:set_title]("Phase [Analytic]")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.20, 0.03, 0.25])
cb = fig[:colorbar](axfft, cax=cax)

fig[:subplots_adjust](hspace=0.25, top=0.95, bottom=0.1)

plt[:savefig]("interpolation_phase.png")
