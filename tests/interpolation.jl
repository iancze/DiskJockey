push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

using gridding
using visibilities
using constants
using image
using gauss_model

# Test to see if we get the convolutional interpolation correct by using a 2D Elliptical Gaussian.
# Realistic Gaussian will have scale dimensions (fatter in x direction)
const mu_x = 1.0 # [arcsec]
const mu_y = 1.0 # [arcsec]
const s_x = 1.2 # [arcsec]
const s_y = 1.0 # [arcsec]
const p0 = [mu_x, mu_y, s_x, s_y] # [arcsec]

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
nx = 128
ny = 128

ra = fftspace(10., nx) # [arcsec]
dec = fftspace(10., ny) # [arcsec]

# convert ra and dec in [arcsec] to radians, and then take the sin to convert to ll, mm
ll = sin(ra * arcsec) # direction cosines
mm = sin(dec * arcsec)

img = imageGauss(ll, mm, p0)

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
vis_analytic = FTGauss(uu, vv, p0)

# Complex subtraction
println("Maximum FFT discrepancy:  ", maximum(abs(plain_fft.VV - vis_analytic)))

import PyPlot
import PyPlot.plt
using LaTeXStrings



# Return a normalized instance that is symmetric about 0
function scale(data)
    s = maximum(abs(data))
    return norm = plt.Normalize(vmin=-s, vmax=s, clip=false)
end

#=

fig, ax = plt.subplots(nrows=1, figsize=(5, 5))

ext = (minimum(ra), maximum(ra), minimum(dec), maximum(dec))
# Real, analytic Gaussian
aximg = ax[:imshow](img, interpolation="none", origin="lower", cmap=plt.get_cmap("Greys"), extent=ext) #, norm = scale(img))
ax[:set_title]("image")
ax[:set_xlabel](L"$\alpha$ [arcsec]")
ax[:set_ylabel](L"$\delta$ [arcsec]")
#[left, bottom, width, height]
cax = fig[:add_axes]([0.84, 0.25, 0.03, 0.45])
cb = fig[:colorbar](aximg, cax=cax)

fig[:subplots_adjust](left=0.15, right=0.85, hspace=0.25)
plt.savefig("../plots/gaussian_img.png")


# from here on out, since we are only showing visibilities, this stays same
ext = (minimum(vis_fft.uu), maximum(vis_fft.uu), minimum(vis_fft.vv), maximum(vis_fft.vv))

# Real and imaginary components of analytic Gaussian
fig, ax = plt.subplots(nrows=2, figsize=(5, 8))

re = ax[1][:imshow](real(vis_analytic), interpolation="none", origin="lower", cmap=plt.get_cmap("bwr"), extent=ext, norm = scale(real(vis_analytic)))
ax[1][:set_title]("Real Analytic")
ax[1][:set_xlabel](L"uu [k$\lambda$]")
ax[1][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.65, 0.03, 0.25])
cb = fig[:colorbar](re, cax=cax)

im = ax[2][:imshow](imag(vis_analytic), interpolation="none", origin="lower", cmap=plt.get_cmap("bwr"), extent=ext, norm = scale(imag(vis_analytic)))

ax[2][:set_title]("Imag Analytic")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.15, 0.03, 0.25])
cb = fig[:colorbar](im, cax=cax)

plt.savefig("../plots/gaussian_analytic.png")



# Real and imaginary components of the FFT Gaussian
fig, ax = plt.subplots(nrows=2, figsize=(5, 8))

re = ax[1][:imshow](real(plain_fft.VV), interpolation="none", origin="lower", cmap=plt.get_cmap("bwr"), extent=ext, norm = scale(real(plain_fft.VV)))
ax[1][:set_title]("Real FFT")
ax[1][:set_xlabel](L"uu [k$\lambda$]")
ax[1][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.65, 0.03, 0.25])
cb = fig[:colorbar](re, cax=cax)

im = ax[2][:imshow](imag(plain_fft.VV), interpolation="none", origin="lower", cmap=plt.get_cmap("bwr"), extent=ext, norm = scale(imag(plain_fft.VV)))

ax[2][:set_title]("Imag FFT")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.15, 0.03, 0.25])
cb = fig[:colorbar](im, cax=cax)

plt.savefig("../plots/gaussian_fft.png")



# Difference between the analytic Gaussian and FFT Gaussian
fig, ax = plt.subplots(nrows=2, figsize=(5, 8))

diff = vis_analytic - plain_fft.VV

re = ax[1][:imshow](real(diff), interpolation="none", origin="lower", cmap=plt.get_cmap("bwr"), extent=ext, norm=scale(real(diff)))
ax[1][:set_title]("Real difference")
ax[1][:set_xlabel](L"uu [k$\lambda$]")
ax[1][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.65, 0.03, 0.25])
cb = fig[:colorbar](re, cax=cax)

im = ax[2][:imshow](imag(diff), interpolation="none", origin="lower", cmap=plt.get_cmap("bwr"), extent=ext, norm=scale(imag(diff)))

ax[2][:set_title]("Imag difference")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.15, 0.03, 0.25])
cb = fig[:colorbar](im, cax=cax)


plt.savefig("../plots/gaussian_difference.png")

=#

# Next, choose some (u,v) points within the bounds and see how the interpolated values compare to what the FTGauss would return.

n = 100
uu = linspace(-100, 100, n) # [kλ]
approx = Array(Complex128, n)
analytic = Array(Complex128, n)

for i=1:n
    u = uu[i]
    v = 0.0
    #println("Interpolating at $u, $v")
    approx[i] = interpolate_uv(u, v, vis_fft)
    analytic[i] = FTGauss(u, v, p0)
end

fig = plt.figure()
ax = fig[:add_subplot](111)

nu = length(vis_fft.uu)

zer = zeros(length(vis_fft.uu))

analytic_u = Array(Complex128, nu)
for i=1:nu
    analytic_u[i] = FTGauss(vis_fft.uu[i], 0.0, p0)
end

ax[:plot](vis_fft.uu, zer, ".k", label="Grid spacing")
ax[:plot](uu, real(approx), "ob", label="Approx")
ax[:plot](uu, real(analytic), ".r", label="Analytic")
ax[:plot](vis_fft.uu, real(analytic_u), "or", label="Analytic")
ax[:set_xlim](-100, 100)
ax[:set_xlabel](L"u [k $\lambda$]")

ax[:legend]()

plt.savefig("../plots/interpolation.png")

# Let's try this over a full grid of images.

# Create analytic function on a smaller grid
n = 100
uu = linspace(-150, 150, n)
vv = linspace(-150, 150, n)

vis_analytic_small = FTGauss(uu, vv, p0)

vis_intp = Array(Complex128, n, n)
for i=1:n
    for j=1:n
        vis_intp[j, i] = interpolate_uv(uu[i], vv[j], vis_fft)
    end
end

fig, ax = plt.subplots(nrows=3, figsize=(5, 11))

ext = (minimum(uu), maximum(uu), minimum(vv), maximum(vv))

axan = ax[1][:imshow](real(vis_analytic_small), interpolation="none", origin="lower", cmap=plt.get_cmap("Greys"), extent=ext)
ax[1][:set_title]("Analytic FT")
ax[1][:set_xlabel](L"uu [k$\lambda$]")
ax[1][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.70, 0.03, 0.25])
cb = fig[:colorbar](axan, cax=cax)

axfft = ax[2][:imshow](real(vis_intp), interpolation="none", origin="lower", cmap=plt.get_cmap("Greys"), extent=ext)
ax[2][:set_title]("Interpolated Visibilites from FFT")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.40, 0.03, 0.25])
cb = fig[:colorbar](axfft, cax=cax)

diff = real(vis_analytic_small - vis_intp)

axdif = ax[3][:imshow](diff, interpolation="none", origin="lower", cmap=plt.get_cmap("bwr"), extent=ext, norm=scale(diff))
ax[3][:set_title]("Difference")
ax[3][:set_xlabel](L"uu [k$\lambda$]")
ax[3][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.10, 0.03, 0.25])
cb = fig[:colorbar](axdif, cax=cax)

fig[:subplots_adjust](hspace=0.25, top=0.97, bottom=0.06)

plt.savefig("../plots/interpolation_difference_real.png")



fig, ax = plt.subplots(nrows=3, figsize=(5, 11))

axan = ax[1][:imshow](imag(vis_analytic_small), interpolation="none", origin="lower", cmap=plt.get_cmap("Greys"), extent=ext)
ax[1][:set_title]("Analytic FT")
ax[1][:set_xlabel](L"uu [k$\lambda$]")
ax[1][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.70, 0.03, 0.25])
cb = fig[:colorbar](axan, cax=cax)

axfft = ax[2][:imshow](imag(vis_intp), interpolation="none", origin="lower", cmap=plt.get_cmap("Greys"), extent=ext)
ax[2][:set_title]("Interpolated Visibilites from FFT")
ax[2][:set_xlabel](L"uu [k$\lambda$]")
ax[2][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.40, 0.03, 0.25])
cb = fig[:colorbar](axfft, cax=cax)

diff = imag(vis_analytic_small - vis_intp)

axdif = ax[3][:imshow](diff, interpolation="none", origin="lower", cmap=plt.get_cmap("bwr"), extent=ext, norm=scale(diff))
ax[3][:set_title]("Difference")
ax[3][:set_xlabel](L"uu [k$\lambda$]")
ax[3][:set_ylabel](L"vv [k$\lambda$]")

cax = fig[:add_axes]([0.84, 0.10, 0.03, 0.25])
cb = fig[:colorbar](axdif, cax=cax)

fig[:subplots_adjust](hspace=0.25, top=0.97, bottom=0.06)

plt.savefig("../plots/interpolation_difference_imag.png")
