push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

using gridding
using visibilities
using constants
using image

# Test to see if we get the convolutional interpolation correct by using a 2D Elliptical Gaussian.

# Realistic Gaussian will have scale dimensions (fatter in x direction)
const s_x = 3. * arcsec # [radians]
const s_y = 1.5 * arcsec # [radians]


# Given two arrays of l and m coordinates, fill an array of the Gaussian image following the MATLAB convention.
function imageGauss(ll::Vector{Float64}, mm::Vector{Float64})
    nx = length(ll)
    ny = length(mm)
    img = Array(Float64, ny, nx)
    for i=1:nx
        for j=1:ny
            img[j, i] = exp(-pi * ((ll[i]/s_x)^2 +  (mm[j]/s_y)^2))
        end
    end
    return img
end

# Given two arrays of u and v coordinates, fill an array with the FT of aforementioned Gaussian.
# u, v should be in units of klambda.
function FTGauss(uu::Vector{Float64}, vv::Vector{Float64})
    nu = length(uu)
    nv = length(vv)
    img = Array(Float64, nv, nu)
    for i=1:nu
        for j=1:nv
            img[j, i] = abs(s_x * s_y) * exp(-pi * ((s_x * uu[i])^2 + (s_y * vv[j])^2))
        end
    end
    return img
end

# Say the beam is 0.3", then oversampling by 10 would be 512x512, or so

nx = 512
ny = 512

ra = linspace(-6, 6, nx) # [arcsec]
dec = linspace(-6, 6, ny) # [arcsec]

# convert ra and dec in [arcsec] to radians, and then take the sin to convert to ll, mm
ll = sin(ra * arcsec)
mm = sin(dec * arcsec)

# find the spacing between the elements 
dl = ll[2] - ll[1]
dm = mm[2] - mm[1]

# These are the proper coordinates in the u-v plane.
uu = fftshift(fftfreq(nx, dl)) .* 1e-3
vv = fftshift(fftfreq(ny, dm)) .* 1e-3

img = imageGauss(ll, mm)
lam0 = cc/230.538e9 * 1e4 # [microns]
skim = SkyImage(img, ra, dec, lam0)

vis = FTGauss(uu, vv)

import PyPlot
import PyPlot.plt
using LaTeXStrings

# These plots may be similar, but we would need to show logscale with some set degree of stretch. Right now I can't tell. Also, the DFT is not normallized properly. How do we normalize a 2D FFT?

fig, ax = plt.subplots(nrows=2, figsize=(5, 8))

ext = (minimum(ra), maximum(ra), minimum(dec), maximum(dec))

aximg = ax[1][:imshow](skim.data[:, :, 1], interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
ax[1][:set_title]("image")
ax[1][:set_xlabel](L"$\alpha$ [arcsec]")
ax[1][:set_ylabel](L"$\delta$ [arcsec]")

#[left, bottom, width, height] 
cax = fig[:add_axes]([0.84, 0.65, 0.03, 0.25])
cb = fig[:colorbar](aximg, cax=cax)

ext = (minimum(uu), maximum(uu), minimum(vv), maximum(vv))
axvis = ax[2][:imshow](vis, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
ax[2][:set_title]("Fourier")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.15, 0.03, 0.25])
cb = fig[:colorbar](axvis, cax=cax)

fig[:subplots_adjust](left=0.15, right=0.85, hspace=0.25)
plt.savefig("../plots/gaussian.png")


# FFT the original image and see how well it matches the visibility space
fftvis = fillModelVis(transform(skim))

fig, ax = plt.subplots(nrows=2, figsize=(5, 8))

ext = (minimum(uu), maximum(uu), minimum(vv), maximum(vv))
axan = ax[1][:imshow](vis, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
ax[1][:set_title]("Fourier")
ax[1][:set_xlabel](L"uu [k$\lambda$]")
ax[1][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.65, 0.03, 0.25])
cb = fig[:colorbar](axan, cax=cax)

ext = (minimum(fftvis.uu), maximum(fftvis.uu), minimum(fftvis.vv), maximum(fftvis.vv))
axfft = ax[2][:imshow](real(fftvis.VV), interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)

ax[2][:set_title]("FFT")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")


cax = fig[:add_axes]([0.84, 0.15, 0.03, 0.25])
cb = fig[:colorbar](axfft, cax=cax)

plt.savefig("../plots/gaussian_fft.png")

# Next, choose randomly distributed (u,v) points within some bounds and see how the interpolated values compare to what the FourierGauss would return.
