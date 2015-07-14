# Because I haven't been able to use CFITSIO to read the visibilities, use
# astropy.io.fits to read the visibilities and convert them to HDF5

from astropy.io import fits
import h5py
import numpy as np

cc = 2.99792458e10 # [cm s^-1]

# Reading SMA dataset
fname = "../data/V4046Sgr/V4046Sgr.12CO21.final.vis.fits"
# fname = "data/V4046Sgr.12CO21.model.vis.fits"

f = fits.open(fname)

data = f[0].data
hdr = f[0].header
nfreq = hdr["NAXIS4"]

freqs = hdr["CRVAL4"] + hdr["CDELT4"] * np.arange(nfreq)  # Hz

# Convert each uu, vv cordinate from light nanoseconds to kilolambda,
# depending on which channel frequency we are using
uu = 1e-3 * (freqs * np.tile(data["UU"], (nfreq, 1)).T).T
vv = 1e-3 * (freqs * np.tile(data["VV"], (nfreq, 1)).T).T
# uu, vv are now (nfreq, nvis) shape arrays

shape = uu.shape

vis = np.squeeze(data["DATA"])  # Remove all of the "zombie" 1D columns
# On disk, stored as an (npoints, nfreqs, 3) array, where last dimension is
# (real, imag, weight)

# Read and convert all of these to (nfreq, nvis) arrays
real = vis[:, :, 0].T
imag = vis[:, :, 1].T
weight = vis[:, :, 2].T

f.close()
#Try plotting the locations of the UV points
#
import matplotlib.pyplot as plt

#plt.errorbar(np.sqrt(uu[12,:100]**2 + vv[12,:100]**2), real[12,:100], yerr=1./np.sqrt(weight[12,:100]), ls="", fmt="o")
#plt.show()
#

#plt.hist(np.abs(imag[12] * np.sqrt(weight[12])))
#plt.show()
chan = 12

fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111)
ax.plot(uu[chan], vv[chan], "k.", ms=0.5)
ax.set_xlabel(r"uu [k$\lambda$]", size=8)
ax.set_ylabel(r"vv [k$\lambda$]", size=8)
lim = 400
ax.set_xlim(lim, -lim)
ax.set_ylim(-lim, lim)

labels = ax.get_xticklabels()
for label in labels:
    label.set_rotation(40)
ax.tick_params(axis="both", which="major", labelsize=8)

fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.8)

plt.savefig("../plots/V4046Sgr/uv_spacings.svg")
plt.savefig("../plots/V4046Sgr/uv_spacings.png")
