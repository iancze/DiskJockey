#!/usr/bin/env python

# Convert old config files using logM_gas to a new value of Sigma_c
import argparse
import yaml
from gas_mass_conversions import logsigma_to_logM

parser = argparse.ArgumentParser(description="Convert the value of logSigma_c to logM_gas in config.yaml files.")
parser.add_argument("--config", help="name of the config file used for the run.", default="config.yaml")

args = parser.parse_args()



f = open(args.config)
config = yaml.load(f)
f.close()

model = config["model"]


p = config["parameters"]

if model == "standard":
    log_M = logsigma_to_logM[model](p["r_c"], p["gamma"], p["logSigma_c"])

elif model == "vertical":
    log_M = logsigma_to_logM[model](p["r_c"], p["gamma"], p["logSigma_c"])

elif model == "truncated":
    pass

elif model == "cavity":
    log_M = logsigma_to_logM[model](p["r_c"], p["r_cav"], p["gamma"], p["gamma_cav"], p["logSigma_c"])

else:
    print("Model type not found.")


print("logM", log_M, "[log10 M_sun]")
