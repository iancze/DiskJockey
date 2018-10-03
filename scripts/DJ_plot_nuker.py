#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Load the config file
import argparse

parser = argparse.ArgumentParser(description="Measure statistics across multiple chains.")
parser.add_argument("--burn", type=int, default=0, help="How many samples to discard from the beginning of the chain for burn in.")
parser.add_argument("--draw", type=int, default=10, help="If specified, print out a random sample of N draws from the posterior, after burn in.")
parser.add_argument("--config", help="name of the config file used for the run.", default="config.yaml")
args = parser.parse_args()


# A dictionary of parameter lists for conversion.
registered_params = {"standard": ["M_star", "r_c", "T_10", "q", "gamma", "logSigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"],
"truncated": ["M_star", "r_c", "T_10", "q", "gamma", "gamma_e", "logSigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"],
"cavity": ["M_star", "r_c", "r_cav", "T_10", "q", "gamma", "gamma_cav", "logSigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"],
"vertical": ["M_star", "r_c", "T_10m", "q_m", "T_10a", "q_a", "T_freeze", "X_freeze", "sigma_s", "gamma", "h", "delta", "logSigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"],
"verticalEta": ["M_star", "r_c", "T_10m", "q_m", "T_freeze", "X_freeze", "sigma_s", "gamma", "h", "eta", "delta", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"],
"nuker": ["M_star", "r_c", "T_10", "q", "gamma", "alpha", "beta", "Sigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]}

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


# If we can tell what type of model we were sampling, we can give everything appropriate labels.
# Otherwise, we'll just use default indexes.
import yaml
f = open(args.config)
config = yaml.load(f)
f.close()

# Get the radii cells from the config file.
rs = np.logspace(np.log10(config["grid"]["r_in"]), np.log10(config["grid"]["r_out"]), config["grid"]["nr"])

inds = np.random.randint(len(flatchain), size=args.draw)

r_c = flatchain[inds,1]
gamma = flatchain[inds,4]
alpha = 10**flatchain[inds,5]
beta = flatchain[inds,6]
Sigma_c = 10**flatchain[inds,7]

# print(r_c, gamma, alpha, beta, Sigma_c)

def nuker(r_c, gamma, alpha, beta, Sigma_c):
    return Sigma_c * (rs/r_c)**(-gamma) * (1 + (rs/r_c)**alpha)**((gamma - beta)/alpha)


Ss = []
Mtot = []

AU = 1.4959787066e13 # [cm]
M_sun = 1.99e33 # [g]

for (r, g, a, b, S) in zip(r_c, gamma, alpha, beta, Sigma_c):
    S = nuker(r, g, a, b, S)
    M = 2 * np.pi * np.trapz(S * rs * AU, rs * AU) / M_sun

    # Convert from g / AU^2 to Msun

    Ss.append(S)
    Mtot.append(M)

fig, ax = plt.subplots(nrows=1, figsize=(5,5))
for S in Ss:
    ax.loglog(rs, S, "k", lw=0.2, alpha=0.5)

ax.set_ylabel(r"$\Sigma(r)$")
ax.set_xlabel(r"$r$ [au]")
fig.savefig("nuker.png", dpi=120)

print("M_disk values [Msun]", Mtot)
