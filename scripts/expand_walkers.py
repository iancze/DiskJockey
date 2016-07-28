#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Measure statistics across multiple chains.")
parser.add_argument("--draw", type=int, help="If specified, print out a random sample of N draws from the posterior, after burn in. Otherwise use the same number of walkers/dim as the previous run.")
parser.add_argument("--frac", help="The uncertainty by which to expand the distance posterior [%]", default=1.0)
parser.add_argument("--burn", type=int, default=0, help="How many samples to discard from the beginning of the chain for burn in.")
parser.add_argument("--new_pos", help="Create a new pos0 array with this filename using the number of walkers contained in draw.", default="pos0_d.npy")
parser.add_argument("--config", help="name of the config file used for the run.", default="config.yaml")

args = parser.parse_args()

frac = args.frac/100.

import numpy as np

# This first bit of code is run for every invocation of the script
chain = np.load("chain.npy")

# Truncate burn in from chain
chain = chain[:, args.burn:, :]

# Convention within the Julia EnsembleSampler is
# ndim, niter, nwalkers = chain.shape
# However, when written to disk, we should have been following the emcee convention
nwalkers, niter, ndim = chain.shape

draw = args.draw if args.draw is not None else int((ndim + 1) * (nwalkers/ndim))
print("Using {} walkers.".format(draw))

nsamples = nwalkers * niter
# Flatchain is made after the walkers have been burned
flatchain = np.reshape(chain, (nsamples, ndim))

# Save it after cutting out burn-in
print("Overwriting flatchain.npy")
np.save("flatchain.npy", flatchain)

# If we can tell what type of model we were sampling, we can give everything appropriate labels.
# Otherwise, we'll just use default indexes.
import yaml

f = open(args.config)
config = yaml.load(f)
f.close()

model = config["model"]
assert config["fix_d"] == True, "Apparantly the config file says that this previous run was made with distance floating?"

# index at which to insert d before (python style) 6
if model == "standard":
    d_index = 6

elif model == "truncated":
    d_index = 7

elif model == "cavity":
    d_index = 8

elif model == "vertical":
    d_index = 8

else:
    print("Model type not found.")


# draw samples from the posterior

inds = np.random.randint(nsamples, size=draw)
pos0 = flatchain[inds]

# Create a sample of distance positions that we can insert into the flatchain
dpc_s = config["parameters"]["dpc"] * np.random.uniform(low=1 - frac, high=1 + frac, size=(draw, 1))


first = pos0[:,:d_index]
middle = dpc_s
last = pos0[:,d_index:]

# pos0 has shape nsamples, ndim
pos0 = np.concatenate((first, middle, last), axis=1)

print("New draws are")
for i in range(draw):
    print(pos0[i])


np.save(args.new_pos, pos0.T)
print("New positions saved to", args.new_pos)
