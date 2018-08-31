# This notebook is designed to allow you to tweak how you might like your walkers initialized. Edit the cells as you see fit and then proceed to evaluate each cell and save the final `pos0.npy` file.

import numpy as np


# Generally, you want at least a few walkers for each dimension you may be exploring. To start with we are only using `4`, but you may want to eventually increase this to `8` or `10` if you have more cores available for computation. For the **standard** model, there are 12 parameters.

nparam = 14
nwalkers = 4 * nparam
print(nwalkers)


# If you are fixing dparmeters, then you should change the previous line to reflect the number of parameters you will be sampling.

# Below, we create an array of starting walker positions, similar to how `emcee` is initialized. You should tweak the `low` and `high` ranges to correspond to a small guess around your starting position.

p0 = np.array([np.random.uniform(1.03, 1.05, nwalkers), # mass [M_sun]
              np.random.uniform(20., 21.0, nwalkers), #r_c [AU]
              np.random.uniform(110., 115, nwalkers), #T_10 [K]
              np.random.uniform(0.50, 0.55, nwalkers), # q
              np.random.uniform(-0.5, 0.0, nwalkers), # gamma
              np.random.uniform(2.0, 3.0, nwalkers), # alpha
              np.random.uniform(6.0, 7.0, nwalkers), # beta
              np.random.uniform(-3.4, -3.1, nwalkers), #log10 Sigma_c [log10 g/cm^2]
              np.random.uniform(0.17, 0.18, nwalkers), #xi [km/s]
              np.random.uniform(145.0, 150.0, nwalkers), #inc [degrees]
              np.random.uniform(40.0, 41.0, nwalkers), #PA [degrees]
              np.random.uniform(-0.1, 0.1, nwalkers), #vz [km/s]
              np.random.uniform(-0.1, 0.1, nwalkers), #mu_a [arcsec]
              np.random.uniform(-0.1, 0.1, nwalkers)]) #mu_d [arcsec]


# Save the new position file to disk
np.save("pos0.npy", p0)
