#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Measure statistics across multiple chains.")
parser.add_argument("--burn", type=int, default=0, help="How many samples to discard from the beginning of the chain for burn in.")
parser.add_argument("--draw", type=int, help="If specified, print out a random sample of N draws from the posterior, after burn in.")
parser.add_argument("--filter", type=float, default=0.0, help="Used in conjunction with --draw. Only take samples with lnprob this percentile[0 - 100) and above. Default is to take all.")
parser.add_argument("--new_pos", help="If specified, create a new pos0 array with this filename using the number of walkers contained in draw.")
parser.add_argument("--config", help="name of the config file used for the run.", default="config.yaml")
parser.add_argument("--hdis", action="store_true", help="Plot the highest density intervals.")
parser.add_argument("--tri", help="Plot the triangle too.", action="store_true")
parser.add_argument("--drop", action="store_true", help="Drop the samples which have lnp==-np.inf")
parser.add_argument("--interactive", action="store_true", help="Pop up the walker window so that you can zoom around.")

args = parser.parse_args()
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# This first bit of code is run for every invocation of the script
chain = np.load("chain.npy")

# Load the lnprobabilities
lnprobs = np.load("lnprob.npy")
# Truncate for burn in, shape (nwalkers, niter)
lnprobs = lnprobs[:, args.burn:]

flat_lnprobs = lnprobs.flatten()

# Set a colorscale for the lnprobs
cmap = matplotlib.cm.get_cmap("brg")

final_lnprobs = lnprobs[:, -1]
norm = matplotlib.colors.Normalize(vmin=np.min(final_lnprobs), vmax=np.max(final_lnprobs))

# Determine colors based on the ending lnprob of each walker
colors = [cmap(norm(val)) for val in final_lnprobs]

# Truncate burn in from chain
chain = chain[:, args.burn:, :]

# Convention within the Julia EnsembleSampler is
# ndim, niter, nwalkers = chain.shape
# However, when written to disk, we should have been following the emcee convention
nwalkers, niter, ndim = chain.shape

print("Previous run used {} walkers.".format(nwalkers))

nsamples = nwalkers * niter
# Flatchain is made after the walkers have been burned
flatchain = np.reshape(chain, (nsamples, ndim))

# Calculate the AIC and BIC using the maximum value of the likelihood
AIC = 2 * ndim - 2 * np.max(flat_lnprobs)
print("AIC : ", AIC)
print("a lower AIC means the model is closer to the 'truth'... as far as you trust diagnostics like these.")
# BIC = np.log(n_data) * ndim - 2 * np.max(flat_lnprobs)

# Keep only the samples which haven't evaluated to -np.inf (the prior disallows them). This usually originates from using a starting position which is already outside the prior.
if args.drop:
    ind = flat_lnprobs > -np.inf
    flat_lnprobs = flat_lnprobs[ind]
    flatchain = flatchain[ind]

nsamples = flatchain.shape[0]

# Save it after cutting out burn-in and -np.inf samples
print("Overwriting flatchain.npy")
np.save("flatchain.npy", flatchain)


# A dictionary of parameter lists for conversion.
registered_params = {"standard": ["M_star", "r_c", "T_10", "q", "gamma", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"],
"truncated": ["M_star", "r_c", "T_10", "q", "gamma", "gamma_e", "logSigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"],
"cavity": ["M_star", "r_c", "r_cav", "T_10", "q", "gamma", "gamma_cav", "logSigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"],
"vertical": ["M_star", "r_c", "T_10m", "q_m", "T_10a", "q_a", "T_freeze", "X_freeze", "sigma_s", "gamma", "h", "delta", "logSigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"],
"verticalEta": ["M_star", "r_c", "T_10m", "q_m", "T_freeze", "X_freeze", "sigma_s", "gamma", "h", "eta", "delta", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"],
"nuker": ["M_star", "r_c", "T_10", "q", "gamma", "log_alpha", "beta", "log_M_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]}

registered_labels = {"standard": [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\gamma$", r"$\log M_\mathrm{gas} \quad \log [M_\odot]$", r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"],
"truncated": [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\gamma$", r"$\gamma_e$", r"$\log \Sigma_c \quad \log [\mathrm{g/cm}^2]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"],
"cavity": [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$r_\mathrm{cav}$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\gamma$", r"$\gamma_\mathrm{cav}$", r"$\log \Sigma_c \quad \log [\mathrm{g/cm}^2]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"],
"vertical": [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10m}$ [K]", r"$q_m$", r"$T_{10a}$ [K]", r"$q_a$", r"T_\mathrm{freeze}", r"$X_\mathrm{freeze}$", r"$\sigma_s$", r"$\gamma$", r"$h$", r"$\delta$", r"$\log \Sigma_c \quad \log [\mathrm{g/cm}^2]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"],
"verticalEta": [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10m}$ [K]", r"$q_m$", r"T_\mathrm{freeze}", r"$X_\mathrm{freeze}$", r"$\sigma_s$", r"$\gamma$", r"$h$", r"$\eta$", r"$\delta$", r"$\log M_\mathrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"],
"nuker": [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\gamma$", r"$\log \alpha$", r"$\beta$", r"$\log M_\mathrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]}

def get_labels(model, fix_params):
    '''
    Collect the labels for each model, so that we can plot.
    '''
    reg_params = registered_params[model]
    reg_labels = registered_labels[model]
    fit_index = [i for (i,param) in enumerate(reg_params) if param not in fix_params]

    labels = [reg_labels[i] for i in fit_index]

    return labels

# If we can tell what type of model we were sampling, we can give everything appropriate labels.
# Otherwise, we'll just use default indexes.
import yaml
f = open(args.config)
config = yaml.load(f, Loader=yaml.FullLoader)
f.close()

model = config["model"]
labels = get_labels(model, config["fix_params"])

# except:
#     print("No appropriate config file found, using default labels.")
#     labels = ["par {}".format(i) for i in range(ndim)]

if args.draw is not None:
    # draw samples from the posterior

    # Take only those samples above some lnprob floor.
    # This functionality allows us to discard problematic walkers carried over between runs.
    floor = np.percentile(flat_lnprobs, args.filter)
    flatchain_filtered = flatchain[(flat_lnprobs > floor)]

    inds = np.random.randint(len(flatchain_filtered), size=args.draw)
    pos0 = flatchain_filtered[inds]

    for i in range(args.draw):
        print(pos0[i])

    if args.new_pos:
        np.save(args.new_pos, pos0.T)

    import sys
    sys.exit()


# if we have a filter argument, remove these walkers from the plot, so that we can get
# better look at the progress of the chains
if args.filter > 0:
    # calculate the median lnprob for each walker
    walker_lnprobs = np.median(lnprobs, axis=1)
    floor = np.percentile(walker_lnprobs, args.filter)
    walker_filtered_ind = (walker_lnprobs > floor)

    lnprobs = lnprobs[walker_filtered_ind]
    chain = chain[walker_filtered_ind]
    nwalkers, niter, ndim = chain.shape
    print("After filtering, there are {} walkers.".format(nwalkers))



import matplotlib.pyplot as plt

iterations = np.arange(niter)
step = 100

# Plot a stream plot of the position of each walker as a function of ensemble iteration.
fig, ax = plt.subplots(nrows=(ndim + 1), ncols=1, figsize=(10, 1.5 * ndim))

#Plot the lnprob in the first panel
for j in range(nwalkers):
    ax[0].plot(iterations, lnprobs[j], lw=0.15, color=colors[j])

avg = np.average(lnprobs, axis=0)
ax[0].plot(iterations, avg, lw=1.1, color="w")
ax[0].plot(iterations, avg, lw=0.9, color="b")
ax[0].set_ylabel("lnprob")

# Plot the walkers for each dimension in the subsequent panels.
for i in range(ndim):
    for j in range(nwalkers):
        ax[i +1].plot(iterations, chain[j, :, i], lw=0.15, color=colors[j])

    # also plot the instanteous walker average to watch for drift
    avg = np.average(chain[:, :, i], axis=0)
    ax[i+1].plot(iterations, avg, lw=1.1, color="w")
    ax[i+1].plot(iterations, avg, lw=0.9, color="b")
    ax[i+1].set_ylabel(labels[i])

ax[-1].set_xlabel("Iteration")

if args.interactive:
    fig.subplots_adjust(bottom=0.0, top=1.0, right=1.0, hspace=0.0)
    plt.show()
else:
    fig.savefig("walkers.png", dpi=300)


def hdi(samples, bins=40):
    '''
    Compute the highest density interval for each 1D marginal distribution.
    '''
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

    plus = binHDI[-1]-binmax
    minus = binmax - binHDI[0]

    return {"bin_centers":bin_centers, "hist":hist, "max":binmax, "lower":binHDI[0], "upper":binHDI[-1], "plus":plus, "minus":minus, "dbin":dbin, "level":level}


def plot_hdis(flatchain, fname="hdi.png"):
    '''
    Plot the bins with highlighted highest density interval ranges, and print their ranges to the terminal.
    '''
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

    fig.subplots_adjust(hspace=0.5, bottom=0.05, top=0.99)
    fig.savefig(fname)

if args.hdis:
    plot_hdis(flatchain)


# Make the triangle plot
if args.tri:
    import corner
    figure = corner.corner(flatchain, bins=30, labels=labels, quantiles=[0.16, 0.5, 0.84], plot_contours=True, plot_datapoints=True, show_titles=True)
    figure.savefig("triangle.png")
else:
    print("Not plotting triangle, no --tri flag.")

# Compute the autocorrelation time, following emcee
# Notes here: http://dfm.io/posts/autocorr/
print("Integrated autocorrelation time")
from emcee import autocorr
print(autocorr.integrated_time(chain))
