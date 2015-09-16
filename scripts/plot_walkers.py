#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Measure statistics across multiple chains.")
parser.add_argument("--burn", type=int, default=0, help="How many samples to discard from the beginning of the chain for burn in.")

args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt
import triangle


chain = np.load("chain.npy")

# Truncate burn in from chain

chain = chain[:, args.burn:, :]

ndim, niter, nwalkers = chain.shape

fig, ax = plt.subplots(nrows=ndim, ncols=1, figsize=(10, 1.5 * ndim))

iterations = np.arange(niter)

for i in range(ndim):
    for j in range(nwalkers):
        ax[i].plot(iterations, chain[i, :, j], lw=0.2, color="k")


ax[-1].set_xlabel("Iteration")

fig.savefig("walkers.png")

flatchain = np.reshape(chain, (ndim, niter * nwalkers))

# flatchain = np.load("flatchain.npy")

# labels = ["a", "b"]
figure = triangle.corner(flatchain.T, quantiles=[0.16, 0.5, 0.84], plot_contours=True, plot_datapoints=False, show_titles=True)

figure.savefig("triangle.png")
