#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Convert SMA FITS files into an HDF5 file for JudithExcalibur.")
parser.add_argument("FITS", help="The SMA FITS file.")
parser.add_argument("--out", default="data.hdf5", help="The output file.")
parser.add_argument("--plot", action="store_true", help="Make a plot of the UV coverage.")
args = parser.parse_args()

# Because I haven't been able to use Julia's CFITSIO to read the visibilities, use
# astropy.io.fits to read the visibilities and convert them to HDF5
from astropy.io import fits
import h5py
import numpy as np

cc = 2.99792458e10 # [cm s^-1]

# Reading SMA dataset
fname = args.FITS

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

# Now, stuff each of these into an HDF5 file.
fid = h5py.File(args.out, "w")

# Convert the frequencies from Hz to micron.
fid.create_dataset("lams", (nfreq,), dtype="float64")[:] = cc/freqs*1e4 #[microns]

fid.create_dataset("uu", shape, dtype="float64")[:,:] = uu
fid.create_dataset("vv", shape, dtype="float64")[:,:] = vv

fid.create_dataset("real", shape, dtype="float64")[:,:] = real
fid.create_dataset("imag", shape, dtype="float64")[:,:] = imag

fid.create_dataset("invsig", shape, dtype="float64")[:,:] = np.sqrt(weight)

fid.close()


if args.plot:

    #Try plotting the locations of the UV points
    import matplotlib.pyplot as plt

    #plt.errorbar(np.sqrt(uu[12,:100]**2 + vv[12,:100]**2), real[12,:100], yerr=1./np.sqrt(weight[12,:100]), ls="", fmt="o")
    #plt.show()
    #

    #plt.hist(np.abs(imag[12] * np.sqrt(weight[12])))
    #plt.show()
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot(uu[0], vv[0], ".")
    ax.set_xlabel(r"UU [k$\lambda$]")
    ax.set_ylabel(r"VV [k$\lambda$]")
    ax.set_xlim(max(uu[0]), min(uu[0]))
    #ax.set_ylim(-75, 75)
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.15)

    plt.savefig("uv_coverage.png")