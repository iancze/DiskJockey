#!/usr/bin/env python

# Use plot.ly to visualize walkers.

import numpy as np

chain = np.load("chain.npy")

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

x = np.arange(niter)

# Perhaps only send the last 200 samples to see if this speeds things up?


for i in range(ndim):
    for j in range(nwalkers):
        if niter > 500:
            trace = go.Scatter(x=x, y=chain[j,niter-500:,i], line={"color":"black", "width":0.4}, hoverinfo="none")
        else:
            trace = go.Scatter(x=x, y=chain[j,:,i], line={"color":"black", "width":0.4}, hoverinfo="none")
        fig.append_trace(trace, i + 1, 1)


fig['layout'].update(height=(200 * ndim), width=900, title='Walkers', showlegend=False)
plot_url = py.plot(fig, filename='Walkers')
