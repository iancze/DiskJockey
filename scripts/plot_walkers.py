#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Measure statistics across multiple chains.")
parser.add_argument("--burn", type=int, default=0, help="How many samples to discard from the beginning of the chain for burn in.")
parser.add_argument("--draw", type=int, help="If specified, print out a random sample of N draws from the posterior, after burn in.")
parser.add_argument("--new_pos", help="If specified, create a new pos0 array with this filename using the number of walkers contained in draw.")


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

# Save it again
print("Overwriting flatchain.npy")
np.save("flatchain.npy", flatchain)

if args.draw is not None:
    # draw samples from the posterior

    inds = np.random.randint(nsamples, size=args.draw)
    pos0 = flatchain[inds]

    for i in range(args.draw):
        print(pos0[i])

    if args.new_pos:
        np.save(args.new_pos, pos0.T)

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

flatchain = np.load("flatchain.npy")

# labels = ["a", "b"]
#     # Fixed distance
labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$r_\textrm{in}$ [AU]",
r"$r_\textrm{cav}$ [AU]", r"$\delta$", r"$T_{10}$ [K]", r"$q$",
r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$d$ [pc]",
r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']",
r"$\mu_\delta$ ['']"]
figure = triangle.corner(flatchain, labels=labels, quantiles=[0.16, 0.5, 0.84], plot_contours=True, plot_datapoints=False, show_titles=True)

figure.savefig("triangle.png")
