import numpy as np


# Generally, you want at least a few walkers for each dimension you may be exploring. For the **truncated** model, there are 13 parameters.

nparam = 13
nwalkers = 4 * nparam
print(nwalkers)

# If you are fixing dparmeters, then you should change the previous line to reflect the number of parameters you will be sampling.

# Below, we create an array of starting walker positions, similar to how `emcee` is initialized. You should tweak the `low` and `high` ranges to correspond to a small guess around your starting position.

p0 = np.array([np.random.uniform(1.03, 1.05, nwalkers), # mass [M_sun]
              np.random.uniform(30., 50.0, nwalkers), #r_c [AU]
              np.random.uniform(110., 115, nwalkers), #T_10 [K]
              np.random.uniform(0.70, 0.71, nwalkers), # q
               np.random.uniform(0.0, 1.5, nwalkers), # gamma_e
              np.random.uniform(-3.4, -3.5, nwalkers), #log10 Sigma_c [log10 g/cm^2]
              np.random.uniform(0.17, 0.18, nwalkers), #xi [km/s]
               np.random.uniform(144.0, 145.0, nwalkers), #dpc [pc]
              np.random.uniform(159.0, 160.0, nwalkers), #inc [degrees]
              np.random.uniform(40.0, 41.0, nwalkers), #PA [degrees]
              np.random.uniform(-0.1, 0.1, nwalkers), #vz [km/s]
              np.random.uniform(-0.1, 0.1, nwalkers), #mu_a [arcsec]
              np.random.uniform(-0.1, 0.1, nwalkers)]) #mu_d [arcsec]


# Save the new position file to disk
np.save("pos0.npy", p0)
