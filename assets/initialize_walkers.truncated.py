
# coding: utf-8

# This notebook is designed to allow you to tweak how you might like your walkers initialized. Edit the cells as you see fit and then proceed to evaluate each cell and save the final `pos0.npy` file.

# In[1]:


import numpy as np


# Generally, you want at least a few walkers for each dimension you may be exploring. For the **truncated** model, there are 13 parameters.

# In[4]:


nparam = 13
nwalkers = 4 * nparam
print(nwalkers)


# If you are fixing distance, then you should change the previous line to 
#     
#     nparam = 12
# 
# and then comment out the `dpc` row below. Below, we create an array of starting walker positions, similar to how `emcee` is initialized. You should tweak the `low` and `high` ranges to correspond to a small guess around your starting position.

# In[5]:


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


# In[6]:


# Just to check we have the right shape
p0.shape


# In[7]:


# Save the new position file to disk
np.save("pos0.npy", p0)


# In[8]:


# Just to check that we have written the file, you can read it back in and check that it has the proper shape.
np.load("pos0.npy").shape

