#!/usr/bin/bash

for f in amr_grid.inp camera_wavelength_micron.inp chmaps_hires_linear.png chmaps_hires_log.png chmaps_blur_linear.png config.yaml gas_temperature.inp gas_velocity.inp image.out lines.inp microturbulence.inp molecule_co.inp numberdens_co.inp radmc3d.inp radmc3d.out wavelength_micron.inp chain.npy lnprob.npy triangle.png walkers.png config.yaml InitializeWalkers.ipynb model.hdf5 resid.hdf5 pos0.npy
do
  rm -v $f
done

rm -v *png
rm -rf output
