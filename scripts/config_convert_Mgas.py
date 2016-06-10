#!/usr/bin/env python

# Convert old config files using logM_gas to a new value of Sigma_c
import argparse
import yaml
from gas_mass_conversions import logM_to_logsigma

parser = argparse.ArgumentParser(description="Convert the value of logM_gas to Sigma_c in config.yaml files.")
parser.add_argument("--config", help="name of the config file used for the run.", default="config.yaml")

args = parser.parse_args()



f = open(args.config)
config = yaml.load(f)
f.close()

model = config["model"]


p = config["parameters"]

if model == "standard":
    logSigma_c = logM_to_logsigma[model](p["r_c"], p["gamma"], p["logM_gas"])

elif model == "vertical":
    logSigma_c = logM_to_logsigma[model](p["r_c"], p["gamma"], p["logM_gas"])

elif model == "truncated":
    pass

elif model == "cavity":
    logSigma_c = logM_to_logsigma[model](p["r_c"], p["r_cav"], p["gamma"], p["gamma_cav"], p["logM_gas"])

else:
    print("Model type not found.")


print("logSigma_c", logSigma_c, "[g/cm^2]")
