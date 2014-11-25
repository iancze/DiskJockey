# Because I haven't been able to use CFITSIO to read the visibilities, use
# astropy.io.fits to read the visibilities and convert them to HDF5

from astropy.io import fits
# import h5py
import numpy as np

# Reading SMA dataset
fname = "data/V4046Sgr.12CO21.final.vis.fits"

f = fits.open(fname)

data = f[0].data
hdr = f[0].header
nfreq = hdr["NAXIS4"]

freqs = hdr["CRVAL4"] + hdr["CDELT4"] * np.arange(nfreq)  # Hz

# Convert each uu, vv cordinate from radians^{-1} to kilolambda
uu = 1e-3 * (freqs * np.tile(data["UU"], (nfreq, 1)).T).T
vv = 1e-3 * (freqs * np.tile(data["VV"], (nfreq, 1)).T).T
# uu, vv are now (nfreq, nvis) shape arrays

vis = np.squeeze(data["DATA"])  # Remove all of the excess 1D columns
# On disk, stored as an (npoints, nfreqs, 3) array, where last dimension is
# (real, imag, weight)

# Read and convert all of these to (nfreq, nvis) arrays
real = vis[:, :, 0].T
imag = vis[:, :, 1].T
weight = vis[:, :, 2].T

#Try plotting the locations of the UV points

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.plot(uu[10], vv[10], ".")
ax.set_xlabel(r"UU [k$\lambda$]")
ax.set_ylabel(r"YY [k$\lambda$]")
fig.subplots_adjust(left=0.2, right=0.8, bottom=0.15)

plt.savefig("plots/uv_spacings.png")

