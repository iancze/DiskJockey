import numpy as np
from emcee import EnsembleSampler

def lnprob(p):
    x, y = p
    lnp = -((1.0 - x)**2 + 100 * (y - x**2)**2)
    return lnp


ndim, nwalkers = 2, 40

p0 = np.array([np.random.rand(ndim) for i in range(nwalkers)])

sampler = EnsembleSampler(nwalkers, ndim, lnprob)
p0, prob, state = sampler.run_mcmc(p0, 10000)
#
sampler.reset()
#
p0, prob, state = sampler.run_mcmc(p0, 10000)
#
np.save("chain.npy", sampler.chain)
