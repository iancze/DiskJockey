#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Use plot.ly to visualize walkers.")
parser.add_argument("--burn", type=int, default=0, help="How many samples to discard from the beginning of the chain for burn in.")
parser.add_argument("--chain", default="chain.npy", help="The name of the file storing the walker positions.")
parser.add_argument("--name", default="walkers", help="The name of the object that we are fitting. The plot.ly plots will show up under this label.")
parser.add_argument("--config", help="name of the config file used for the run.", default="../../config.yaml")
parser.add_argument("--open", help="Pop up the graph when finished in your browser?", action="store_true")

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

fig = tools.make_subplots(rows=ndim, cols=1, shared_xaxes=True, vertical_spacing=0.005)

# If we can tell what type of model we were sampling, we can give everything appropriate labels.
# Otherwise, we'll just use default indexes.
import yaml
try:
    f = open(args.config)
    config = yaml.load(f)
    f.close()

    model = config["model"]
    fix_d = config["fix_d"]

    if model == "standard":
        if fix_d:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
        else:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
    elif model == "truncated":
        if fix_d:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_in$ [AU]", r"$r_out$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
        else:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_in$ [AU]", r"$r_out$ [AU]", r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
    elif model == "cavity":
        if fix_d:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$r_\textrm{cav}$ [AU]",       r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
        else:
            labels = [r"$M_\ast\quad [M_\odot]$", r"$r_c$ [AU]", r"$r_\textrm{cav}$ [AU]",       r"$T_{10}$ [K]", r"$q$", r"$\log M_\textrm{gas} \quad \log [M_\odot]$",  r"$\xi$ [km/s]", r"$d$ [pc]", r"$i_d \quad [{}^\circ]$", r"PA $[{}^\circ]$", r"$v_r$ [km/s]", r"$\mu_\alpha$ ['']", r"$\mu_\delta$ ['']"]
    else:
        print("Model type not found, using default labels.")
        labels = ["par {}".format(i) for i in range(ndim)]

except:
    print("No appropriate config file found, using default labels.")
    labels = ["par {}".format(i) for i in range(ndim)]

# Only plot 12 walkers per dimension
nlines = 12 if nwalkers > 12 else nwalkers

# Only plot every 10th sample.
stride = 10
if niter > 100:
    chain = chain[:,::stride,:]
    # Iteration #
    x = stride * np.arange(chain.shape[1])
else:
    x = np.arange(chain.shape[1])

for i in range(ndim):
    # To save memory, at most plot 12 walkers per dimension.

    for j in range(nlines):
        y = chain[j,:,i]

        trace = go.Scatter(x=x, y=y, line={"color":"black", "width":0.4}, hoverinfo="none", mode="lines")

        fig.append_trace(trace, i + 1, 1)

        fig['layout']['yaxis{}'.format(i + 1)].update(title=labels[i])

# print(fig['layout'])
# fig['layout']['yaxis1'].update(title='yaxis 1 title')


fig['layout'].update(height=(200 * ndim), width=900, title=args.name, showlegend=False)

plot_url = py.plot(fig, filename=args.name, auto_open=args.open)
