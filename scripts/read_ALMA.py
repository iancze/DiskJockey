#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Convert ALMA NPZ save files into an HDF5 file for JudithExcalibur.", prefix_chars="@")
parser.add_argument("NPZ", help="The ALMA NPZ file produced by Sean.")
parser.add_argument("nu0", type=float, help="The starting frequency for index 0, in Hz")
parser.add_argument("dnu", type=float, help="The change in frequency with array index. Can be negative. In Hz.")
parser.add_argument("@out", default="data.hdf5", help="The output file.")
parser.add_argument("@plot", action="store_true", help="Make a plot of the UV coverage.")
args = parser.parse_args()

import h5py
import numpy as np

cc = 2.99792458e10 # [cm s^-1]

data = np.load(args.NPZ)


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

#####################################
# # For AKSco.vis.npz
# # The index 0 channel has a frequency of 230.550055 GHz.
# # Each channel has a width of 305.176 kHz
# nu0 = 230.550055e9 # [Hz]
# dnu = 305.176e3 # [Hz]
#####################################

#####################################
# For AKSco.12CO.305kHz.vis.npz
# The index 0 channel has a frequency of 230.546088 GHz.
# Each channel has a width of 305.176 kHz
# nu0 = 230.546088e9 # [Hz]
# dnu = 305.176e3 # [Hz]
#####################################

#####################################
# # For AKSco.12CO.122kHz.vis.npz
# # The index 0 channel has a frequnecy of 230.545996 GHz
# nu0 = 230.545996e9 # [Hz]
# dnu = 122.070e3 # [Hz]
#####################################

#####################################
# For GWOri.13CO.data.vis.npz
# nu0 = 220.371112e9 # [Hz]
# dnu = -122.070e3 # [Hz]
#####################################

#####################################
# For GWOri C18O
# nu0 = 219.532855e9 # [Hz]
# dnu = -122.070e3 # [Hz]
####################################

#####################################
# For DQ Tau 12CO
# nu0 = 230.514180e9  # [Hz]
# dnu = 61.035e3 # [Hz]
####################################


nu0 = args.nu0
dnu = args.dnu

print("nu0: ", nu0)
print("dnu: ", dnu)

freqs = nu0 + np.arange(nchan) * dnu # [Hz]
lams = cc/freqs * 1e4 # [microns]

# The data set is time-averaged into 30s intervals and binned spectrally by a
# factor of 5 (original data has emission in 250 channels!).
# Two polarizations are averaged already.

# u and v are measured in meters, convert to microns and then
# convert these to kilo-lambda
uu = 1e-3 * (np.tile(data["u"] * 1e6, (nchan, 1)).T / lams).T
vv = 1e-3 * (np.tile(data["v"] * 1e6, (nchan, 1)).T / lams).T

# uu, vv are now (nchan, nvis) shape arrays
shape = uu.shape

# Convert these to (nchan, nvis) arrays
real = data["Re"]
imag = data["Im"]
weight = data["Wt"]

# Test to see if there are any autocorrelation baselines present
ant1 = data["ant1"]
ant2 = data["ant2"]
same = (ant1 == ant2)
print("There are {} autocorrelation visibilities.".format(np.sum(same)))

# Now, stuff this into an HDF5 file.
fid = h5py.File(args.out, "w")

# The HDF5 file must be packed with wavelength increasing in order.

if dnu > 0:
    # e.g, like AK Sco 12CO 2-1
    print("Swapping order to be increasing wavelength.")
    # Reverse the order
    #Currently, everything is stored in decreasing wavelength order, lets flip this.
    fid.create_dataset("lams", (nchan,), dtype="float64")[:] = lams[::-1]

    fid.create_dataset("uu", shape, dtype="float64")[:,:] = uu[::-1, :]
    fid.create_dataset("vv", shape, dtype="float64")[:,:] = vv[::-1, :]

    fid.create_dataset("real", shape, dtype="float64")[:,:] = real[::-1, :]
    fid.create_dataset("imag", shape, dtype="float64")[:,:] = imag[::-1, :]

    fid.create_dataset("invsig", shape, dtype="float64")[:,:] = np.sqrt(weight)[::-1, :]

else:
    # e.g., like GW Ori 13CO 2-1
    print("Order is already in increasing wavelength.")
    fid.create_dataset("lams", (nchan,), dtype="float64")[:] = lams

    fid.create_dataset("uu", shape, dtype="float64")[:,:] = uu
    fid.create_dataset("vv", shape, dtype="float64")[:,:] = vv

    fid.create_dataset("real", shape, dtype="float64")[:,:] = real
    fid.create_dataset("imag", shape, dtype="float64")[:,:] = imag

    fid.create_dataset("invsig", shape, dtype="float64")[:,:] = np.sqrt(weight)


fid.close()

if args.plot:

    # Plot the UV samples for the first channel
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot(uu[0], vv[0], ".")
    ax.set_xlabel(r"UU [k$\lambda$]")
    ax.set_ylabel(r"VV [k$\lambda$]")
    ax.set_xlim(max(uu[0]), min(uu[0]))
    # ax.set_ylim(-75, 75)
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.15)

    plt.savefig("uv_spacings.png")
