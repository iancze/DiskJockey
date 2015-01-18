# Picking up where `write_SMA.jl` left off

from astropy.io import fits
import h5py
import numpy as np
import shutil

cc = 2.99792458e10 # [cm s^-1]

# Read all of the data from the HDF5 file
fid = h5py.File("data/V4046Sgr_model.hdf5", "r")

freqs = cc/fid["lams"][:]*1e4 # [Hz]
uu = fid["uu"][:,:] # [klam]
vv = fid["vv"][:,:] # [klam]
real = fid["real"][:,:] # [Jy]
imag = fid["imag"][:,:] # [Jy]
weight = fid["invsig"][:,:]**2
fid.close()

# (nfreq, nvis) arrays

# New SMA dataset
fname = "data/V4046Sgr.12CO21.model.vis.fits"

# Copy the original dataset to something new
shutil.copy("data/V4046Sgr.12CO21.final.vis.fits", fname)

# Overwrite a copy of the original dataset with these values.
hdulist = fits.open(fname, mode="update")

# data["DATA"] is originally (23302, 1, 1, 25, 1, 3)
# (nvis, 1, 1, nchan, 1, 3)
# 3 is real, imag, weight)

# Concatenate the different parts of the visibility
D = np.array([real, imag, weight]).T
# (nvis, nfreq, 3)

# Add the zombie dimensions back in
vis = D[:, np.newaxis, np.newaxis, :, np.newaxis, :]

# Stuff the new visibilities into the existing dataset
hdulist[0].data["DATA"][:] = vis

# Write the changes to disk
hdulist.flush()
hdulist.close()
