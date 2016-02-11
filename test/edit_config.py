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

import numpy as np
nwalkers = 4

if args.fix_distance:
    # Just change the distance
    config["fix_d"] = True
else:
    config["fix_d"] = False

if args.model == "standard":
    config["model"] = "standard"
    # Also set parameters to a reasonable dictionary
    config["parameters"] = {"M_star": 1.03, "PA": 152.0, "T_10": 91.85, "dpc": 145.0, "gamma": 1.0, "incl": 45.0, "ksi": 0.2, "logM_gas": -3.8, "mu_DEC": 0.0, "mu_RA": 0.0, "q": 0.5, "r_c": 500.0, "vel": -1.6}

    if args.fix_distance:
        p0 = np.array([np.random.uniform(1.03, 1.05, nwalkers), # mass [M_sun]
              np.random.uniform(300., 420.0, nwalkers), #r_c [AU]
              np.random.uniform(80., 90, nwalkers), #T_10 [K]
              np.random.uniform(0.55, 0.65, nwalkers), # q
              np.random.uniform(-4.0, -3.8, nwalkers), #log10 M_gas [log10 M_sun]
              np.random.uniform(0.4, 0.5, nwalkers), #xi [km/s]
              np.random.uniform(43.0, 47.0, nwalkers), #inc [degrees]
              np.random.uniform(152.0, 156.0, nwalkers), #PA [degrees]
              np.random.uniform(-2.0, -1.9, nwalkers), #vz [km/s]
              np.random.uniform(0.2, 0.3, nwalkers), #mu_a [arcsec]
              np.random.uniform(-0.3, -0.1, nwalkers)]) #mu_d [arcsec]
    else:
        p0 = np.array([np.random.uniform(1.03, 1.05, nwalkers), # mass [M_sun]
              np.random.uniform(300., 420.0, nwalkers), #r_c [AU]
              np.random.uniform(80., 90, nwalkers), #T_10 [K]
              np.random.uniform(0.55, 0.65, nwalkers), # q
              np.random.uniform(-4.0, -3.8, nwalkers), #log10 M_gas [log10 M_sun]
              np.random.uniform(0.4, 0.5, nwalkers), #xi [km/s]
              np.random.uniform(135.0, 145.0, nwalkers), #dpc [pc]
              np.random.uniform(43.0, 47.0, nwalkers), #inc [degrees]
              np.random.uniform(152.0, 156.0, nwalkers), #PA [degrees]
              np.random.uniform(-2.0, -1.9, nwalkers), #vz [km/s]
              np.random.uniform(0.2, 0.3, nwalkers), #mu_a [arcsec]
              np.random.uniform(-0.3, -0.1, nwalkers)]) #mu_d [arcsec]

elif args.model == "truncated":
    config["model"] = "truncated"
    pass

elif args.model == "cavity":
    config["model"] = "cavity"
    config["parameters"] = {"M_star": 1.03, "PA": 152.0, "T_10": 91.85, "dpc": 145.0, "gamma": 1.0, "incl": 45.0, "ksi": 0.2, "logM_gas": -3.8, "mu_DEC": 0.0, "mu_RA": 0.0, "q": 0.5, "r_c": 500.0, "r_cav":30, "vel": 0.0}

    if args.fix_distance:
        p0 = np.array([np.random.uniform(1.03, 1.05, nwalkers), # mass [M_sun]
              np.random.uniform(300., 420.0, nwalkers), #r_c [AU]
              np.random.uniform(10., 20.0, nwalkers), #r_cav [AU]
              np.random.uniform(80., 90, nwalkers), #T_10 [K]
              np.random.uniform(0.55, 0.65, nwalkers), # q
              np.random.uniform(-4.0, -3.8, nwalkers), #log10 M_gas [log10 M_sun]
              np.random.uniform(0.4, 0.5, nwalkers), #xi [km/s]
              np.random.uniform(43.0, 47.0, nwalkers), #inc [degrees]
              np.random.uniform(152.0, 156.0, nwalkers), #PA [degrees]
              np.random.uniform(-2.0, -1.9, nwalkers), #vz [km/s]
              np.random.uniform(0.2, 0.3, nwalkers), #mu_a [arcsec]
              np.random.uniform(-0.3, -0.1, nwalkers)]) #mu_d [arcsec]
    else:
        p0 = np.array([np.random.uniform(1.03, 1.05, nwalkers), # mass [M_sun]
              np.random.uniform(300., 420.0, nwalkers), #r_c [AU]
              np.random.uniform(10., 20.0, nwalkers), #r_cav [AU]
              np.random.uniform(80., 90, nwalkers), #T_10 [K]
              np.random.uniform(0.55, 0.65, nwalkers), # q
              np.random.uniform(-4.0, -3.8, nwalkers), #log10 M_gas [log10 M_sun]
              np.random.uniform(0.4, 0.5, nwalkers), #xi [km/s]
              np.random.uniform(135.0, 145.0, nwalkers), #dpc [pc]
              np.random.uniform(43.0, 47.0, nwalkers), #inc [degrees]
              np.random.uniform(152.0, 156.0, nwalkers), #PA [degrees]
              np.random.uniform(-2.0, -1.9, nwalkers), #vz [km/s]
              np.random.uniform(0.2, 0.3, nwalkers), #mu_a [arcsec]
              np.random.uniform(-0.3, -0.1, nwalkers)]) #mu_d [arcsec]


g = open("config.yaml", mode="w")
yaml.dump(config, g)
g.close()

# Save the new position file to disk
np.save("pos0.npy", p0)
