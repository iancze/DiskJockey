
# coding: utf-8

# This notebook is designed to allow you to tweak how you might like your walkers initialized. Edit the cells as you see fit and then proceed to evaluate each cell and save the final `pos0.npy` file.

# In[3]:


import numpy as np


# Generally, you want at least a few walkers for each dimension you may be exploring. To start with we are only using `4`, but you may want to eventually increase this to `8` or `10` if you have more cores available for computation. For the **vertical** model, there are 13 parameters.

# In[4]:


nparam = 13
nwalkers = 4 * nparam
print(nwalkers)


# If you are allowing distance to float, then you should change the previous line to
# 
#     nparam = 14
# 
# and then uncomment out the `dpc` row below. Below, we create an array of starting walker positions, similar to how `emcee` is initialized. You should tweak the `low` and `high` ranges to correspond to a small guess around your starting position.

# In[5]:


p0 = np.array([np.random.uniform(1.03, 1.05, nwalkers), # mass [M_sun]
              np.random.uniform(20., 21.0, nwalkers), #r_c [AU]
              np.random.uniform(30., 40., nwalkers), #T_10m [K]
              np.random.uniform(0.50, 0.55, nwalkers), # q_m 
              np.random.uniform(110., 115, nwalkers), #T_10a [K]
              np.random.uniform(0.50, 0.55, nwalkers), # q_a 
              np.random.uniform(-3.4, -3.1, nwalkers), #log10 Sigma_c [log10 g/cm^2]
              np.random.uniform(0.17, 0.18, nwalkers), #xi [km/s]
#                np.random.uniform(144.0, 146.0, nwalkers), #dpc [pc]
              np.random.uniform(44.0, 46.0, nwalkers), #inc [degrees]
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

