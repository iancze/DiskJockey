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

# Overwrite a copy of the original dataset with these values.
hdulist = fits.open(args.out_model, mode="update")

# data["DATA"] is originally (23302, 1, 1, 25, 1, 3)
# (nvis, 1, 1, nchan, 1, 3)
# 3 is real, imag, weight)

# print(hdulist[0].data["DATA"].shape)

# Concatenate the different parts of the visibility
D = np.array([real, imag, weight]).T
# (nvis, nfreq, 3)

# Add the zombie dimensions back in
vis = D[:, np.newaxis, np.newaxis, :, np.newaxis, :]

# print(vis.shape)

# Stuff the new visibilities into the existing dataset
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
vis = D[:, np.newaxis, np.newaxis, :, np.newaxis, :]
hdulist[0].data["DATA"][:] = vis
hdulist.flush()
hdulist.close()
