#!/usr/bin/env python

# Original script written by Jane Huang, CfA

# Reads in an image.out file from RADMC-3D and creates a new FITS file.
# Ideal for conversion from RADMC output to CASA simobserve, for ALMA proposals

import argparse

parser = argparse.ArgumentParser(description="Convert RADMC-3D image.out into a FITS file. Optionally provide information that will be added to the header of the FITS file.")
parser.add_argument("--image", default="image.out", help="The name of the file created by RADMC-3D.")
parser.add_argument("--fits", default="image.fits", help="The name of the FITS file to which you want to export the image.")
parser.add_argument("--dpc", default=140., type=float, help="At what distance [pc] is this source? Assumes Taurus distance by default.")
parser.add_argument("--RA", default=0, type=float, help="Assign this as the RA to the object.")
parser.add_argument("--DEC", default=0, type=float, help="Assign this as the DEC to the object.")
args = parser.parse_args()

import numpy as np
from astropy.io import fits, ascii
import os

AU = 1.49597870700e13 # [cm]
pc = 3.08567758e18 # [cm]
cc = 2.99792458e10 # [cm s^-1]

dpc = args.dpc # [pc]
RA = args.RA
DEC = args.DEC


# Read in the file from the RADMC-3D format
imagefile = open(args.image)
iformat = imagefile.readline()

im_nx, im_ny = imagefile.readline().split() #number of pixels along x and y axes
im_nx = int(im_nx)
im_ny = int(im_ny)

nlam = int(imagefile.readline())

pixsize_x, pixsize_y = imagefile.readline().split() #pixel sizes in cm for observer at infinity
pixsize_x = float(pixsize_x)
pixsize_y = float(pixsize_y)

# Differential RA and DEC
# ra = ((np.arange(im_nx) + 0.5) - im_nx/2.) * pixsize_x/(AU*dpc)
# dec = ((np.arange(im_ny) + 0.5) - im_ny/2.) * pixsize_y/(AU*dpc)

imvals = ascii.read(args.image, format = 'fast_csv', guess = False, data_start = 4, fast_reader = {'use_fast_converter':True})['1']

lams = imvals[:nlam]

# Convert lams (in microns) into frequencies in Hz
freqs =  cc / (lams * 1e-4) # [Hz]

CRVAL3 = freqs[0]

single_chan = False

if len(lams) > 1:
    dnu = freqs[1] - freqs[0]
    assert np.all((np.diff(freqs) - dnu)/dnu < 0.01), "This conversion script does not support unequal channel widths for multi-channel imagecubes."
else:
    "Single channel image, assuming bandwidth (CDELT3) is 2GHz"
    dnu = 2e9 #[GHz]
    single_chan = True
CDELT3 = dnu

pixsize = pixsize_x*pixsize_y/(dpc*pc)**2 #pixel size in steradians

#RADMC gives intensities in erg cm^(-2) s^(-1) Hz^(-1) ster^(-1); convert to Jy/pixel
intensities = np.reshape(imvals[nlam:],[nlam, im_ny, im_nx])* pixsize*10**23

if single_chan:
    total_flux = np.sum(intensities)
    freq = freqs[0]
    print("Total (spatially-integrated) flux density at {:.3e} Hz flux in image: {:.3e} Jy".format(freq, total_flux))
else:
    channel_flux = np.sum(intensities, axis=(1,2))
    # Make a spectrum plot
    # Hz on one side
    # km/s on the other
    # Use central channel as 0 velocity.
    nu0 = np.average(freqs)
    vels = cc * (freqs - nu0)/nu0 * 1e-5

    integrated_flux = np.trapz(channel_flux, -vels)

    print("Total (spatially and velocity) integrated flux is {:.3e} Jy - km/s, integrated over {:.1f} km/s".format(integrated_flux, vels[0] - vels[-1]))
    print("See spectrum.png for visual representaion.")

    import matplotlib.pyplot as plt
    # plt.plot(freqs, channel_flux)
    plt.plot(vels, channel_flux)
    plt.xlabel(r"$\nu$ [Hz]")
    plt.xlabel("v [km/s]")
    plt.ylabel(r"$f_\nu$ [Jy]")
    plt.savefig("spectrum.png")


# Estimate the total flux in the image. Sum all the pixels in each channel to get the flux density (e.g., Jy, measured at each channel.)

# If there is only one channel, print out a message saying that this is the flux density measured at the frequency of the image.

# Might also want to make a plot of the spectrum itself, in units of Jy (flux density).

# Then, if there is more than one channel, integrate along the frequency dimension to get a measure of the integrated line Jy-km/s, and say what it is and over what velocity range it was integrated.

# Estimate the total flux in each channel, to make a spectrum.

# Convert to float32 to store in FITS?
intensities = intensities.astype('float32')


# Now, export the image to a FITS file

#check back later to make this more general (i.e., deal with the bug in cvel)
hdu = fits.PrimaryHDU(intensities)
header = hdu.header
header['EPOCH'] = 2000.
header['EQUINOX'] = 2000.
# Latitude and Longitude of the pole of the coordinate system.
header['LATPOLE'] = -1.436915713634E+01
header['LONPOLE'] = 180.

# Define the RA coordinate
header['CTYPE1'] = 'RA---SIN'
header['CUNIT1'] = 'DEG'

cdelt1 = -pixsize_x/(pc*dpc)*180/np.pi
header['CDELT1'] = cdelt1
# Pixel coordinates of the reference point. For example, if the image is 256 pixels wide, then this
# would refer to the center of the 127th pixel.
if im_nx % 2 == 0:
    header['CRPIX1'] = int(0.5*im_nx + 1)
else:
    header['CRPIX1'] = int(0.5*im_nx+0.5)
header['CRVAL1'] = RA

# Define the DEC coordinate
header['CTYPE2'] = 'DEC--SIN'
header['CUNIT2'] = 'DEG'
header['CDELT2'] = -1*cdelt1 #assumes square image
if im_ny % 2 == 0:
    header['CRPIX2'] = int(0.5*im_ny + 1)
else:
    header['CRPIX2'] = int(0.5*im_ny+0.5)
header['CRVAL2'] = DEC

# Define the frequency coordiante
header['CTYPE3'] = 'FREQ'
header['CUNIT3'] = 'Hz'
header['CRPIX3'] = 1.
header['CDELT3'] = CDELT3
header['CRVAL3'] = CRVAL3

header['SPECSYS'] = 'LSRK'
header['VELREF'] = 257
header['BSCALE'] = 1.
header['BZERO'] = 0.
header['BUNIT'] = 'JY/PIXEL'
header['BTYPE']='Intensity'

hdu.writeto(args.fits, overwrite = True)
