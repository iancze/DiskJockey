# Picking up where `write_ALMA.jl` left off

from astropy.io import fits
import h5py
import numpy as np
import shutil

cc = 2.99792458e10 # [cm s^-1]

# Read all of the data from the HDF5 file
fid = h5py.File("../data/AKSco/AKSco_model.hdf5", "r")

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

# Sean delivered the data set as an NPZ file.
# I kept everything in increasing *wavelength* order
# He kept everything in increasing *frequency* order

# This means we will have to reverse the order of the real, imaginary, and weights

# This file has categories
# ['Re', 'Wt', 'u', 'Im', 'v']

# len(data["u"]) => 24555
# len(data["v"]) => 24555
# data['Re'].shape => (50, 24555)
# data['Im'].shape => (50, 24555)
# data['Wt'].shape => (50, 24555)

np.savez("../data/AKSco/AKSco_model.vis.npz", u=u, v=v, Re=real[::-1, :], Im=imag[::-1, :], Wt=weight[::-1, :] )
