#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Measure statistics across multiple chains.")
parser.add_argument("--opt_jump", default="opt_jump.npy", help="The numpy save file holding the estimated covariance.")
parser.add_argument("--N", default=10000, help="Number of samples to draw from the covariance matrix.")

args = parser.parse_args()
N = args.N

# Use the opt_jump.npy, sample a bunch, and then plot it later.

import numpy as np
from numpy.random import multivariate_normal as mvn

Sigma = np.load(args.opt_jump)

ndim = Sigma.shape[0]
samples = mvn(np.zeros(ndim), Sigma, N)

np.save("samples.npy", samples)


import triangle

# Fixed distance
labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]",
r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]",
# r"$d$ [pc]",
r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]",
r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]

# Vertical temperature gradient
# labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10,m}$ [K]",
# r"$q_m$", r"$T_{10,a}$ [K]", r"$q_a$", r"$\gamma$", r"$h$", r"$\delta$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]",
# # r"$d$ [pc]",
# r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]",
# r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]

figure = triangle.corner(samples, quantiles=[0.16, 0.5, 0.84],
    plot_contours=True, plot_datapoints=False, labels=labels, show_titles=True)
figure.savefig("fake_triangle.png")
