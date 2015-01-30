import h5py
import numpy as np

cc = 2.99792458e10 # [cm s^-1]

data = np.load("../data/AKSco/AKSco.vis.npz")

# This file has categories
# ['Re', 'Wt', 'u', 'Im', 'v']

# len(data["u"]) => 24555
# len(data["v"]) => 24555

# data['Re'].shape => (50, 24555)
# data['Im'].shape => (50, 24555)
# data['Wt'].shape => (50, 24555)

# There are 24555 visibilities in each channel
# There are 50 channels
nchan, nvis = data["Re"].shape

# The index 0 channel has a frequency of 230.550055 GHz.
# Each channel has a width of 305.176 kHz
nu0 = 230.550055e9 # [Hz]
dnu = 305.176e3 # [Hz]
freqs = nu0 + np.arange(nchan) * dnu # [Hz]
lams = cc/freqs * 1e4 # [microns]

# The data set is time-averaged into 30s intervals and binned spectrally by a
# factor of 5 (original data has emission in 250 channels!).
# Two polarizations are averaged already.

# u and v are measured in meters, convert to microns and then
# convert these to kilo-lambda
uu = 1e-3 * (np.tile(data["u"] * 1e6, (nchan, 1)).T / lams).T
vv = 1e-3 * (np.tile(data["v"] * 1e6, (nchan, 1)).T / lams).T

# Plot the UV samples for a given channel
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.plot(uu[25], vv[25], ".")
ax.set_xlabel(r"UU [k$\lambda$]")
ax.set_ylabel(r"VV [k$\lambda$]")
ax.set_xlim(max(uu[25]), min(uu[25]))
# ax.set_ylim(-75, 75)
fig.subplots_adjust(left=0.2, right=0.8, bottom=0.15)

plt.savefig("../plots/uv_spacings_ALMA.png")

# uu, vv are now (nchan, nvis) shape arrays
shape = uu.shape

# Convert these to (nchan, nvis) arrays
real = data["Re"]
imag = data["Im"]
weight = data["Wt"]

# Now, stuff each of these into an HDF5 file.
fid = h5py.File("../data/AKSco/AKSco.hdf5", "w")

#Currently, everything is stored in decreasing wavelength order, lets flip this.
fid.create_dataset("lams", (nchan,), dtype="float64")[:] = lams[::-1]

fid.create_dataset("uu", shape, dtype="float64")[:,:] = uu[::-1, :]
fid.create_dataset("vv", shape, dtype="float64")[:,:] = vv[::-1, :]

fid.create_dataset("real", shape, dtype="float64")[:,:] = real[::-1, :]
fid.create_dataset("imag", shape, dtype="float64")[:,:] = imag[::-1, :]

fid.create_dataset("invsig", shape, dtype="float64")[:,:] = np.sqrt(weight)[::-1, :]

fid.close()
