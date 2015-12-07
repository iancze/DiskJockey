#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Use plot.ly to visualize walkers.")
parser.add_argument("--burn", type=int, default=0, help="How many samples to discard from the beginning of the chain for burn in.")
parser.add_argument("--chain", default="chain.npy", help="The name of the file storing the walker positions.")
parser.add_argument("--name", default="hist", help="The name of the object that we are fitting. The plot.ly plots will show up under this label.")

args = parser.parse_args()

import numpy as np

chain = np.load(args.chain)

# Convention within the Julia EnsembleSampler is
# ndim, niter, nwalkers = chain.shape
# However, when written to disk, we should have been following the emcee convention
nwalkers, niter, ndim = chain.shape

# nsamples = nwalkers * niter
# Flatchain is made after the walkers have been burned
# flatchain = np.reshape(chain, (nsamples, ndim))


from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go

# First, let's try to make a single 2D contour and a single histogram



fig = tools.make_subplots(rows=ndim, cols=1, shared_xaxes=True, vertical_spacing=0.005)

x = np.arange(niter)

# Label the axes appropriately based upon how many parameters we have
if ndim == 10:
    labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
elif ndim == 11:
    labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
else:
    labels = None
