JudithExcalibur
===============

Copyright Ian Czekala 2014-15

email: iancze@gmail.com

This package derives dynamical masses for T Tauri stars using the Keplerian motion of their circumstellar disks, applied to radio interferometric data from the Atacama Large Millimeter Array (ALMA) and the Submillimeter Array (SMA). **If you use this code or a derivative of it in your research, you must cite [Czekala et al. 2015 ApJ, 806 154C](http://adsabs.harvard.edu/abs/2015ApJ...806..154C).**

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

Recently, I have organized this code into a proper Julia package.

## Installation

First, make sure you have installed RADMC-3D and you can successfully run one of the example scripts.

Because `JudithExcalibur` isn't yet an official Julia package (nor will it likely be), for now, installation involves simply cloning the repository. First, open up a Julia prompt in the REPL, then type

    julia> Pkg.clone("https://github.com/iancze/JudithExcalibur.git")

Lastly, I have written several "driver" command line scripts that are used to perform the actual mass fitting. To complete the installation, you should add these files to your PATH. To figure out where the package is installed

    julia> Pkg.dir("JudithExcalibur")
    "/home/ian/.julia/JudithExcalibur"

Your PATH will vary. The scripts are located inside of a `scripts` directory, so if you are using bash or Z-shell, you will want to add the PATH that looks something like

    export PATH="/home/ian/.julia/JudithExcalibur/scripts:$PATH"

inside of your `.bashrc` or `.zshrc` file. Finally,

    $ source ~/.zshrc

To check that you have properly added the scripts, you can try

    $ JudithInitialize.jl --test
    Your JudithExcalibur scripts are successfully linked.
    Exiting

Note, if you would like to use the plotting scripts (`scripts/plot_model.jl`) or the IO routines for SMA and ALMA data (`scripts/read_SMA.py`), you will need a Python installation with the following packages installed: `numpy`, `scipy`, `matplotlib`, `h5py`, and `astropy`. The anaconda Python distribution would be a great way to take care of this.

## Use

Now, the `JudithExcalibur` package should be installed globally on your system. Because it is likely that you will want to fit more than just one disk, or use different model specifications for a particular disk, the code structure is organized so that you will have a separate directory for each run. For example,

    $ mkdir ExcitingDisk
    $ cd ExcitingDisk

Once you've created a working directory for your project, then you'll want to initialize this directory with a config file

    $ JudithInitialize.jl --new-config
    Copied default config.yaml file to current working directory.
    Exiting

Now, open up `config.yaml` with your favorite text editor and change the fields as you see fit. Next, we'll need to initialize the directory with the appropriate RADMC-3D files specific to your dataset.

    $ JudithInitialize.jl

To run this code on a new dataset using 4 cores, you would run

    mach_three.jl -p 3

where the `-p` flag specifies how many additional processes to spawn, therefore the total number of processes used is `p + 1`. For maximal performance, set `p = ncores - 1` where `ncores` is the number of cores on your machine.
