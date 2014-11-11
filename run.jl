#!/usr/local/bin julia

# Run the RADMC3D FORTRAN program to make spectra and images

# CO J = 2-1 transition at 230.538 GHz

using constants

# iline = 2
incl = 33. #33. # deg. 0 deg = face on, 90 = edge on.
vel = 0.0 # km/s
PA = 90 - 73. # 73 deg. Position angle, runs counter clockwise, due to looking at sky.
npix = 96 # number of pixels, can alternatively specify x and y separately

lam0 = cc/230.538e9 * 1e6 # [microns]


#loads the camera_wavelength_micron.inp file
#run(`radmc3d image incl $incl posang $PA vkms $vel npix $npix loadlambda`)
#run(`radmc3d image incl $incl posang $PA vkms $vel npix $npix lambda $lam0`)
run(`radmc3d spectrum incl $incl posang $PA vkms $vel npix $npix loadlambda`)

#optionally add `doppcatch` to enable Doppler Catching. Seems to give weird results, though.

#run(`radmc3d image iline $iline vkms $vel incl $incl phi $phi`)
