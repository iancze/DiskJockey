#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(prog="parallel.py", description="Run Starfish fitting model in parallel.")
parser.add_argument("--fix-distance", help="Flip fixing the distance.", action="store_true")
parser.add_argument("--model", choices=["standard", "truncated", "cavity"], help="Set the model to this type.")
args = parser.parse_args()


# Use PyYAML to edit the config file

import yaml

f = open("config.yaml")
config = yaml.load(f)
f.close()


if args.fix_distance
    # Just change the distance
    config["fix_d"] = True
else:
    config["fix_d"] = False

if args.model == "standard":
    config["model"] = "standard"
    # Also set parameters to a reasonable dictionary
    config["parameters"] = {"M_star": 1.03, "PA": 152.0, "T_10": 91.85, "dpc": 145.0, "gamma": 1.0, "incl": 45.0, "ksi": 0.2, "logM_gas": -3.8, "mu_DEC": 0.0, "mu_RA": 0.0, "q": 0.5, "r_c": 500.0, "vel": 0.0}

elif args.model == "truncated":
    config["model"] = "truncated"
    pass

elif args.model == "cavity":
    config["model"] = "cavity"
    config["parameters"] = {"M_star": 1.03, "PA": 152.0, "T_10": 91.85, "dpc": 145.0, "gamma": 1.0, "incl": 45.0, "ksi": 0.2, "logM_gas": -3.8, "mu_DEC": 0.0, "mu_RA": 0.0, "q": 0.5, "r_c": 500.0, "r_cav":30, "vel": 0.0}


g = open("config.yaml", mode="w")
yaml.dump(config, g)
g.close()
