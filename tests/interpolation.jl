push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

using gridding
using visibilities
using constants
using image
using gauss_model

# Test to see if we get the convolutional interpolation correct by using a 2D Elliptical Gaussian.

# Apply Hanning tapering to the image.
function hanning!(img::SkyImage)

    r0 = maximum(img.ra) # assumes these axes are symmetric
    d0 = maximum(img.dec)
    nr = length(img.ra)
    rd = length(img.dec)
    nlam = length(img.lams)

    for l=1:nlam
        for i=1:nr
            for j=1:rd
                img.data[j, i, l] = img.data[j, i, l] *
                (0.5 + 0.5 * cos(pi * img.ra[i]/r0)) * (0.5 + 0.5 * cos(pi * img.dec[j]/d0))
            end
        end
    end

end

# If the synthesized beam is 0.3", to oversample this by a factor of 10 would require something like 512x512 pixels
nx = 512
ny = 512

ra = fftspace(8., nx) # [arcsec]
dec = fftspace(8., ny) # [arcsec]

# convert ra and dec in [arcsec] to radians, and then take the sin to convert to ll, mm
ll = sin(ra * arcsec) # direction cosines
mm = sin(dec * arcsec)

img = imageGauss(ll, mm)

lam0 = cc/230.538e9 * 1e4 # [microns]
skim = SkyImage(img, ra, dec, lam0)

# Do one FFT without the correction function
plain_fft = transform(skim)

# Apply the gridding correction function before doing the FFT
corrfun!(skim, 1.0)

# FFT the image and see how well it matches the visibility space
vis_fft = transform(skim)

println("Imaginary FFT: Minimum ", minimum(imag(plain_fft.VV)), " Maximum: ", maximum(imag(plain_fft.VV)))

# Take the vis from vis_fft
uu = vis_fft.uu # [kλ]
vv = vis_fft.vv # [kλ]

# Analytic visibilites
vis_analytic = FTGauss(uu, vv)

println("Real FFT discrepancy: Minimum ", minimum(real(plain_fft.VV) - vis_analytic),
" Maximum: ", maximum(real(plain_fft.VV) - vis_analytic))

import PyPlot
import PyPlot.plt
using LaTeXStrings

fig, ax = plt.subplots(nrows=2, figsize=(5, 8))

ext = (minimum(ra), maximum(ra), minimum(dec), maximum(dec))

# Real, analytic Gaussian
aximg = ax[1][:imshow](img, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
ax[1][:set_title]("image")
ax[1][:set_xlabel](L"$\alpha$ [arcsec]")
ax[1][:set_ylabel](L"$\delta$ [arcsec]")

#[left, bottom, width, height]
cax = fig[:add_axes]([0.84, 0.65, 0.03, 0.25])
cb = fig[:colorbar](aximg, cax=cax)

ext = (minimum(uu), maximum(uu), minimum(vv), maximum(vv))
axvis = ax[2][:imshow](vis_analytic, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
ax[2][:set_title]("Fourier")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.15, 0.03, 0.25])
cb = fig[:colorbar](axvis, cax=cax)

fig[:subplots_adjust](left=0.15, right=0.85, hspace=0.25)
plt.savefig("../plots/gaussian.png")


fig, ax = plt.subplots(nrows=2, figsize=(5, 8))

ext = (minimum(vis_fft.uu), maximum(vis_fft.uu), minimum(vis_fft.vv), maximum(vis_fft.vv))
axan = ax[1][:imshow](real(plain_fft.VV), interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
ax[1][:set_title]("Real FFT")
ax[1][:set_xlabel](L"uu [k$\lambda$]")
ax[1][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.65, 0.03, 0.25])
cb = fig[:colorbar](axan, cax=cax)

axfft = ax[2][:imshow](imag(plain_fft.VV), interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)

ax[2][:set_title]("Imag FFT")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")


cax = fig[:add_axes]([0.84, 0.15, 0.03, 0.25])
cb = fig[:colorbar](axfft, cax=cax)

plt.savefig("../plots/gaussian_fft.png")

# Plot the relative error comparing the analytic FT to the (real) FFT
fig = plt.figure()
ax = fig[:add_subplot](111)
cbimg = ax[:imshow]((real(plain_fft.VV) - vis_analytic), interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)

cax = fig[:add_axes]([0.84, 0.15, 0.03, 0.45])
cb = fig[:colorbar](cbimg, cax=cax)
ax[:set_title]("FT Difference")

fig[:subplots_adjust](left=0.15, right=0.85, hspace=0.25)
plt.savefig("../plots/gaussian_difference.png")


# Next, choose randomly distributed (u,v) points within some bounds and see how the interpolated values compare to what the FourierGauss would return.

n = 100
uu = linspace(-100, 100, n) # [kλ]
approx = Array(Complex128, n)
analytic = Array(Float64, n)

for i=1:n
    u = uu[i]
    v = 0.0
    approx[i] = interpolate_uv(u, v, vis_fft)
    analytic[i] = FTGauss(u, v)
end

#@printf("(%.2f, %.2f): ", u, v)
#print("Approximate ", approx)
#println(" Analytic ", analytic)


import PyPlot
import PyPlot.plt
using LaTeXStrings

fig = plt.figure()
ax = fig[:add_subplot](111)

nu = length(vis_fft.uu)

zer = zeros(length(vis_fft.uu))

analytic_u = Array(Float64, nu)
for i=1:nu
    analytic_u[i] = FTGauss(vis_fft.uu[i], 0.0)
end

ax[:plot](vis_fft.uu, zer, ".k", label="Grid spacing")
ax[:plot](uu, real(approx), "ob", label="Approx")
ax[:plot](uu, analytic, ".r", label="Analytic")
ax[:plot](vis_fft.uu, analytic_u, "or", label="Analytic")
ax[:set_xlim](-100, 100)
ax[:set_xlabel](L"u [k $\lambda$]")

ax[:legend]()

plt.savefig("../plots/interpolation.png")

# Let's try this over a full grid of images.

# Create analytic function on a smaller grid
n = 100
uu = linspace(-150, 150, n)
vv = linspace(-150, 150, n)

vis_analytic_small = FTGauss(uu, vv)

vis_intp = Array(Complex128, n, n)
for i=1:n
    for j=1:n
        vis_intp[j, i] = interpolate_uv(uu[i], vv[j], vis_fft)
    end
end

fig, ax = plt.subplots(nrows=3, figsize=(5, 11))

ext = (minimum(uu), maximum(uu), minimum(vv), maximum(vv))

axan = ax[1][:imshow](vis_analytic_small, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
ax[1][:set_title]("Analytic FT")
ax[1][:set_xlabel](L"uu [k$\lambda$]")
ax[1][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.70, 0.03, 0.25])
cb = fig[:colorbar](axan, cax=cax)

axfft = ax[2][:imshow](real(vis_intp), interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
ax[2][:set_title]("Interpolated Visibilites from FFT")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.40, 0.03, 0.25])
cb = fig[:colorbar](axfft, cax=cax)

axdif = ax[3][:imshow](vis_analytic_small - real(vis_intp), interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
ax[3][:set_title]("Difference")
ax[3][:set_xlabel](L"uu [k$\lambda$]")
ax[3][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.10, 0.03, 0.25])
cb = fig[:colorbar](axdif, cax=cax)

fig[:subplots_adjust](hspace=0.25, top=0.97, bottom=0.06)

plt.savefig("../plots/interpolation_difference.png")
