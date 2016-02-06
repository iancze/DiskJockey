#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Convert HDF5 files into SMA FITS files.")
parser.add_argument("FITS", help="The original FITS data set, so that we can copy it to stuff in new values.")
parser.add_argument("--fname-model", default="model.hdf5", help="The name of the model visibilities HDF5 file.")
parser.add_argument("--fname-resid", default="resid.hdf5", help="The name of the model visibilities HDF5 file.")
parser.add_argument("--out-model", default="model.vis.fits", help="The output file for the model.")
parser.add_argument("--out-resid", default="resid.vis.fits", help="The output file for the residuals.")

args = parser.parse_args()

# Picking up where `write_model.jl` left off

from astropy.io import fits
import h5py
import numpy as np
import shutil

cc = 2.99792458e10 # [cm s^-1]

# Read the model from the HDF5 file
fid = h5py.File(args.fname_model, "r")

freqs = cc/fid["lams"][:]*1e4 # [Hz]
uu = fid["uu"][:,:] # [klam]
vv = fid["vv"][:,:] # [klam]
real = fid["real"][:,:] # [Jy]
imag = fid["imag"][:,:] # [Jy]
weight = fid["invsig"][:,:]**2
fid.close()

# (nfreq, nvis) arrays

# Copy the original dataset to something new, so that we can update the new data set.
shutil.copy(args.FITS, args.out_model)
shutil.copy(args.FITS, args.out_resid)


# Open the old file so we can check if we had the wavelengths increasing or decreasing
# Reading SMA dataset
f = fits.open(args.FITS)
data = f[0].data
hdr = f[0].header
nfreq = hdr["NAXIS4"]
ofreqs = hdr["CRVAL4"] + hdr["CDELT4"] * np.arange(nfreq)  # Hz
f.close()

print("Original frequencies are", ofreqs)


# Overwrite a copy of the original dataset with these values.
hdulist = fits.open(args.out_model, mode="update")

# data["DATA"] is originally (23302, 1, 1, 25, 1, 3)
# (nvis, 1, 1, nchan, 1, 3)
# 3 is real, imag, weight)

# print(hdulist[0].data["DATA"].shape)

# Concatenate the different parts of the visibility
D = np.array([real, imag, weight]).T

if ofreqs[-1] > ofreqs[0]:
    print("Originally stored with decreasing wavelength, reversing model order.")
    D = D[:, ::-1, :]

else:
    print("Originally stored with increasing wavelength, keeping model order the same.")

# Add the zombie dimensions back in to the visibilities
vis = D[:, np.newaxis, np.newaxis, :, np.newaxis, :]

hdulist[0].data["DATA"][:] = vis

# Write the changes to disk
hdulist.flush()
hdulist.close()

# Now do the same for the residuals
# Read the model from the HDF5 file
fid = h5py.File(args.fname_resid, "r")

freqs = cc/fid["lams"][:]*1e4 # [Hz]
uu = fid["uu"][:,:] # [klam]
vv = fid["vv"][:,:] # [klam]
real = fid["real"][:,:] # [Jy]
imag = fid["imag"][:,:] # [Jy]
weight = fid["invsig"][:,:]**2
fid.close()

hdulist = fits.open(args.out_resid, mode="update")
D = np.array([real, imag, weight]).T

if ofreqs[-1] > ofreqs[0]:
    print("Originally stored with decreasing wavelength, reversing residual order.")
    D = D[:, ::-1, :]

else:
    print("Originally stored with increasing wavelength, keeping residual order the same.")

vis = D[:, np.newaxis, np.newaxis, :, np.newaxis, :]
hdulist[0].data["DATA"][:] = vis

hdulist.flush()
hdulist.close()
