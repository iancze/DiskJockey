#!/usr/bin/env python

# Picking up where `write_model.jl` left off

import argparse

parser = argparse.ArgumentParser(description="Convert model visibilities in HDF5 format to ALMA NPZ save files.")
parser.add_argument("--fname", default="model.hdf5", help="The name of the model visibilities HDF5 file.")
parser.add_argument("--flip", action="store_true", help="Should the frequencies be packed in a descinding order (e.g., 13CO)?")
parser.add_argument("--out", default="model.vis.npz", help="The output file.")
args = parser.parse_args()


from astropy.io import fits
import h5py
import numpy as np
import shutil

cc = 2.99792458e10 # [cm s^-1]

# Read all of the data from the HDF5 file
fid = h5py.File(args.fname, "r")

lams = fid["lams"][:] * 1e-6 # [m]
uu = fid["uu"][:,:] # [klam]
vv = fid["vv"][:,:] # [klam]
real = fid["real"][:,:] # [Jy]
imag = fid["imag"][:,:] # [Jy]
weight = fid["invsig"][:,:]**2
fid.close()

# Convert u and v from kilo-lambda back to meters
u = uu[0,:] * 1e3 * lams[0] # [m]
v = vv[0,:] * 1e3 * lams[0] # [m]


# This means we will have to reverse the order of the real, imaginary, and weights

# This file has categories
# ['Re', 'Wt', 'u', 'Im', 'v']

# len(data["u"]) => 24555
# len(data["v"]) => 24555
# data['Re'].shape => (50, 24555)
# data['Im'].shape => (50, 24555)
# data['Wt'].shape => (50, 24555)

# Sean delivered the data set as an NPZ file.
# I kept everything in increasing *wavelength* order
# He kept everything in increasing *frequency* order

if args.flip:
    # Therefore, if he gave me a dataset that is with frequency decreasing, don't need to do anything (e.g. 13CO).
    np.savez(args.out, u=u, v=v, Re=real[:, :], Im=imag[:, :], Wt=weight[:, :] )
else:
    # But if he gave me a dataset with frequency increasing, then I need to flip the order here (e.g., 12CO).
    np.savez(args.out, u=u, v=v, Re=real[::-1, :], Im=imag[::-1, :], Wt=weight[::-1, :] )
