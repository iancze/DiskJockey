#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Measure statistics across multiple chains.")
parser.add_argument("--burn", type=int, default=0, help="How many samples to discard from the beginning of the chain for burn in.")
parser.add_argument("--draw", type=int, help="If specified, print out a random sample of N draws from the posterior, after burn in.")
parser.add_argument("--new_pos", help="If specified, create a new pos0 array with this filename using the number of walkers contained in draw.")
parser.add_argument("--config", help="name of the config file used for the run.", default="../../config.yaml")

args = parser.parse_args()
import numpy as np


# This first bit of code is run for every invocation of the script
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

# Save it after cutting out burn-in
print("Overwriting flatchain.npy")
np.save("flatchain.npy", flatchain)

# If we can tell what type of model we were sampling, we can give everything appropriate labels.
# Otherwise, we'll just use default indexes.
import yaml
try:
    f = open(args.config)
    config = yaml.load(f)
    f.close()

    model = config["model"]
    fix_d = config["fix_d"]

    if model == "standard":
        if fix_d:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
        else:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
    elif model == "truncated":
        if fix_d:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_in$ [AU]", r"$r_out$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
        else:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_in$ [AU]", r"$r_out$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
    elif model == "cavity":
        if fix_d:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$r_\textrm{cav}$ [AU]",       r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
        else:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$r_\textrm{cav}$ [AU]",       r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
    else:
        print("Model type not found, using default labels.")
        labels = ["par {}".format(i) for i in range(ndim)]

except:
    print("No appropriate config file found, using default labels.")
    labels = ["par {}".format(i) for i in range(ndim)]

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
        ax[i].plot(iterations, chain[j, :, i], lw=0.15, color="k")

    ax[i].set_xlabel(labels[i])

ax[-1].set_xlabel("Iteration")

fig.savefig("walkers.png")

flatchain = np.load("flatchain.npy")

def hdi(samples, bins=80):

    hist, bin_edges = np.histogram(samples, bins=bins, density=True)
    # convert bin_edges into bin centroids
    bin_centers = bin_edges[:-1] + np.diff(bin_edges)

    dbin = bin_edges[1] - bin_edges[0]
    nbins = len(bin_centers)

    # Now, sort all of the bin heights in decreasing order of probability
    indsort = np.argsort(hist)[::-1]
    histsort = hist[indsort]
    binsort = bin_centers[indsort]

    binmax = binsort[0]

    prob = histsort[0] * dbin
    i = 0
    while prob < 0.683:
        i += 1
        prob = np.sum(histsort[:i] * dbin)

    level = histsort[i]

    indHDI = hist > level
    binHDI = bin_centers[indHDI]

    # print("Ranges: low: {}, max: {}, high: {}".format(binHDI[0], binmax, binHDI[-1]))
    # print("Diffs: max:{}, low:{}, high:{}, dbin:{}".format(binmax, binmax - binHDI[0], binHDI[-1]-binmax, dbin))

    # Now, return everything necessary to make a plot
    # "lower": lower confidence interval
    # "upper": upper confidence interval
    #
    # "plus": + errorbar
    # "minus": - errorbar

    plus = binHDI[-1]-binmax
    minus = binmax - binHDI[0]

    return {"bin_centers":bin_centers, "hist":hist, "max":binmax, "lower":binHDI[0], "upper":binHDI[-1], "plus":plus, "minus":minus, "dbin":dbin, "level":level}


def plot_hdis(flatchain, fname="hdi.png"):

    # Plot the bins with highlighted ranges
    fig,ax = plt.subplots(ncols=1, nrows=ndim, figsize=(6, 1.5 * ndim))

    for i in range(flatchain.shape[1]):
        vals = hdi(flatchain[:,i])

        ax[i].plot(vals["bin_centers"], vals["hist"], ls="steps-mid")
        ax[i].axhline(vals["level"], ls=":", color="k")
        ax[i].set_xlabel(labels[i])
        ax[i].set_ylabel("probability")

        print(labels[i])
        print("Ranges: low: {}, max: {}, high: {}".format(vals["lower"], vals["max"], vals["upper"]))
        print("Diffs: max:{}, low:{}, high:{}, dbin:{}".format(vals["max"], vals["minus"], vals["plus"], vals["dbin"]))
        print()

    fig.subplots_adjust(hspace=0.1)
    fig.savefig(fname)

plot_hdis(flatchain)

# Make the triangle plot
figure = triangle.corner(flatchain, labels=labels, quantiles=[0.16, 0.5, 0.84], plot_contours=True, plot_datapoints=False, show_titles=True)
figure.savefig("triangle.png")
