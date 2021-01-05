#!/usr/bin/env python

# Convert old chain.npy using Sigma_c to logM_gas
import argparse
import yaml
import numpy as np
from gas_mass_conversions import logsigma_to_logM

parser = argparse.ArgumentParser(description="Convert the value of logM_gas to Sigma_c in MCMC sample files (chain.npy).")
parser.add_argument("--config", help="name of the config file used for the run.", default="config.yaml")
parser.add_argument("--old-chain", default="chain.npy", help="Name of the file storing the previous walker positions.")
parser.add_argument("--new-chain", default="chain_logM.npy", help="Name of the file storing the logM")

args = parser.parse_args()

f = open(args.config)
config = yaml.load(f)
f.close()

model = config["model"]

p = config["parameters"]

gamma = p["gamma"]

chain = np.load(args.old_chain)
nwalkers, niter, ndim = chain.shape

# Based upon the parameter order, convert the logMgas samples to logSigma_c and resave the new file. Because all of these parameters appear before distance, we don't need to worry whether distance was fixed or not.

if (model == "standard"):
    r_c = chain[:,:,1]
    logSigma_c = chain[:,:,3]

    chain[:,:,3] = logsigma_to_logM["standard"](r_c, gamma, logSigma_c)

else:
    print("Model type not found.")

# Save the new positions file
np.save(args.new_chain, chain)
