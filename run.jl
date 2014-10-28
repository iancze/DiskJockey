#!/usr/local/bin julia

# Run the RADMC3D FORTRAN program to make spectra and images

# CO J = 2-1 transition at 230.538 GHz

#Compute the various wavelengths at which we want to synthesize images (23 images spaced -4.40 -0- 4.40 km/s).

iline = 2
incl = 33. # deg. 0 deg = face on, 90 = edge on.
vel = 0.0 # km/s
phi = 0.0 # deg. Position angle.

#loads the camera_wavelength_micron.inp file
run(`radmc3d image incl $incl phi $phi loadlambda doppcatch`)

#optionally add `doppcatch` to enable Doppler Catching. Seems to give weird results, though.

#run(`radmc3d image iline $iline vkms $vel incl $incl phi $phi`)
