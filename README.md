JudithExcalibur
===============

Copyright Ian Czekala 2014-15

email: iancze@gmail.com

This package derives dynamical masses for T Tauri stars using the Keplerian motion of their circumstellar disks. Applied to SMA and ALMA CO line data.

See an explanation of [how dynamical mass measurements](http://iancze.github.io/dynamical/) work.

# Organization

JudithExcalibur is designed to forward model inteferometric observations of protoplanetary disks, for the purpose of deriving a precise measurement of the central (sub-)stellar mass.

This package relies upon the excellent radiative synthesis package RADMC-3D to perform the radiative transfer of the disk model. Comprehensive documentation for RADMC-3D can be found [here](http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/). Before using JudithExcalibur, you (obviously) need to [install](http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/) RADMC-3D first.

There are several stages to evaluating the best-fitting parameters of a disk structure model:

1. model specification
2. radiative synthesis
3. visibility interpolation
4. likelihood evaluation

Because spectral line datasets are large, and synthesizing models is computationally expensive, I have designed this package to work in a parallel environment, where the calculations for each frequency channel can be distributed to independent processors. Therefore, please keep this architecture in mind when navigating the source code.

The following is a description of the most important files in the package

**model.jl**: Contains the actual specification of the parametric disk model, as well as the tools to write to disk the synthesis files RADMC-3D requires.

**image.jl**: Contains type definitions to read images produced by RADMC-3D, as well as convert from physical coordinates to sky coordinates.

**visibilities.jl**: Contains type definitions to hold the dataset and the model visibilities. Additionally contains functions to apply phase shifts to the visibilities corresponding to shifts in the image plane. Also contains functions to FFT images to the visibility plane.

**gridding.jl**: Contains the prolate-spheroidal wave function definitions from Schwab 1984, used when doing the visibility interpolations.

**parallel.jl**: A simple pipe-like parallel implementation designed to farm parameters to individual RADMC-3D processes in order to synthesize chunks of channels, as opposed to the full spectrum all at once.

**LittleMC.jl**: A simple Metropolis-Hastings Markov-Chain Monte Carlo implementation designed to sample the posterior distribution of parameters.

**mach_three.jl**: The driver script for performing MCMC with a disk model. The main disk-related functions are `fprob` and `f`. The organization of this file is designed so that the likelihood call can be parallelized using the framework written in `parallel.jl`.

Parameter files are specified using YAML. Some example parameter files are located in the `scripts/` directory.

To run this code on a new dataset using 4 cores, you would run

    mach_three.jl -p 3 scripts/my_data.yaml
