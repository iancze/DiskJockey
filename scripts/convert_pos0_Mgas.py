#!/usr/bin/env python

# Convert old pos0.npy using logM_gas to a new value of Sigma_c
import argparse
import yaml
import numpy as np
from gas_mass_conversions import logM_to_logsigma

parser = argparse.ArgumentParser(description="Convert the value of logM_gas to Sigma_c in walker position (pos0.npy) files.")
parser.add_argument("--config", help="name of the config file used for the run.", default="config.yaml")
parser.add_argument("--old-pos0", default="pos0.npy", help="Name of the file storing the previous walker positions.")
parser.add_argument("--new-pos0", default="pos0_Sigma_c.npy", help="Name of the file storing the Sigma_c ")

args = parser.parse_args()

f = open(args.config)
config = yaml.load(f)
f.close()

model = config["model"]
fix_d = config["fix_d"]
p = config["parameters"]

gamma = p["gamma"]

pos = np.load(args.old_pos0)
ndim, nwalkers = pos.shape

# Based upon the parameter order, convert the logMgas samples to logSigma_c and resave the new file. Because all of these parameters appear before distance, we don't need to worry whether distance was fixed or not.

if (model == "standard") or (model == "vertical"):
    r_c = pos[1]
    logM_gas = pos[4]

    pos[4] = logM_to_logsigma["standard"](r_c, gamma, logM_gas)

elif model == "truncated":
    print("truncated not implemented yet")

elif model == "cavity":
    r_c = pos[1]
    r_cav = pos[2]
    gamma_cav = pos[5]
    logM_gas = pos[6]

    pos[6] = logM_to_logsigma["cavity"](r_c, r_cav, gamma, gamma_cav, logM_gas)

else:
    print("Model type not found.")

# Save the new positions file
np.save(args.new_pos0, pos)
