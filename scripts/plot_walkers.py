#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Measure statistics across multiple chains.")
parser.add_argument("--burn", type=int, default=0, help="How many samples to discard from the beginning of the chain for burn in.")
parser.add_argument("--draw", type=int, help="If specified, print out a random sample of N draws from the posterior, after burn in.")

args = parser.parse_args()

import numpy as np

chain = np.load("chain.npy")

# Truncate burn in from chain
chain = chain[:, args.burn:, :]

# Convention within the Julia EnsembleSampler is
# ndim, niter, nwalkers = chain.shape

# However, when written to disk, we should have been following the emcee convention
nwalkers, niter, ndim = chain.shape

nsamples = nwalkers * niter
# Flatchain is made after the walkers have been burned
flatchain = np.reshape(chain, (nsamples, ndim))


if args.draw is not None:
    # draw samples from the posterior
    for i in range(args.draw):
        ind = np.random.randint(nsamples)
        print(flatchain[ind])

    import sys
    sys.exit()


import matplotlib.pyplot as plt
import triangle

fig, ax = plt.subplots(nrows=ndim, ncols=1, figsize=(10, 1.5 * ndim))

iterations = np.arange(niter)

for i in range(ndim):
    for j in range(nwalkers):
        ax[i].plot(iterations, chain[j, :, i], lw=0.2, color="k")


ax[-1].set_xlabel("Iteration")

fig.savefig("walkers.png")


# labels = ["a", "b"]
#     # Fixed distance
labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]",
r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]",
# r"$d$ [pc]",
r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]",
r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
figure = triangle.corner(flatchain, labels=labels, quantiles=[0.16, 0.5, 0.84], plot_contours=True, plot_datapoints=False, show_titles=True)

figure.savefig("triangle.png")
