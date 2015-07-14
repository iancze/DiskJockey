#!/usr/bin/env python2

import radmc3dPy as r

imag = r.image.readImage()

# print("There are {} wavelengths in this data cube.".format(imag.nwav))
# print("They are {}".format(imag.wav))

#r.image.plotImage(imag, arcsec=True, dpc=73., log=True, maxlog=5, ifreq=12)

# r.image.plotImage(imag, au=True, log=True, maxlog=5, ifreq=12)

# for viewing dust
r.image.plotImage(imag, au=True, log=True, maxlog=5)
