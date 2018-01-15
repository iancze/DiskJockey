import numpy as np

nparam = 14
nwalkers = 4 * nparam

p0 = np.array([np.random.uniform(1.74, 1.78, nwalkers), # mass [M_sun]
    np.random.uniform(44., 46.0, nwalkers), #r_c [AU]
    np.random.uniform(110., 115, nwalkers), #T_10 [K]
    np.random.uniform(0.60, 0.65, nwalkers), # q
    np.random.uniform(-1.0, 0.0, nwalkers), # gamma 
    np.random.uniform(0.11, 1.0, nwalkers), # log10alpha
    np.random.uniform(6.0, 7.0, nwalkers), # beta
    np.random.uniform(0.0, 0.2, nwalkers), #log10 Sigma_c [log10 g/cm^2]
    np.random.uniform(0.17, 0.18, nwalkers), #xi [km/s]
    np.random.uniform(145.0, 150.0, nwalkers), #inc [degrees]
    np.random.uniform(350.0, 359.0, nwalkers), #PA [degrees]
    np.random.uniform(-31.2, -31.1, nwalkers), #vz [km/s]
    np.random.uniform(0.1, 0.2, nwalkers), #mu_a [arcsec]
    np.random.uniform(-0.6, -0.4, nwalkers)]) #mu_d [arcsec]

np.save("pos0.npy", p0)
