JudithExcalibur
===============

Copyright Ian Czekala 2014-16

email: iancze@gmail.com

This package derives dynamical masses for T Tauri stars using the Keplerian motion of their circumstellar disks, applied to radio interferometric data from the Atacama Large Millimeter Array (ALMA) and the Submillimeter Array (SMA). **If you use this code or a derivative of it in your research, you must cite [Czekala et al. 2015 ApJ, 806 154C](http://adsabs.harvard.edu/abs/2015ApJ...806..154C).**

See an explanation of [how dynamical mass measurements](http://iancze.github.io/dynamical/) work.

Papers published using JudithExcalibur:

* *A Disk-based Dynamical Constraint on the Mass of the Young Binary AK Sco*: [Czekala et al. 2015 ApJ, 806 154C](http://adsabs.harvard.edu/abs/2015ApJ...806..154C)
* *A Disk-based Dynamical Constraint on the Mass of the Young Binary DQ Tau* : [Czekala et al. 2016 2016arXiv160103806C](http://adsabs.harvard.edu/abs/2016arXiv160103806C)

# Organization

JudithExcalibur is designed to forward model inteferometric observations of protoplanetary disks, for the purpose of deriving a precise measurement of the central (sub-)stellar mass.

This package relies upon the excellent radiative synthesis package RADMC-3D to perform the radiative transfer of the disk model. Comprehensive documentation for RADMC-3D can be found [here](http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/). Before using JudithExcalibur, you (obviously) need to [install](http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/) RADMC-3D first.

There are several stages to evaluating the best-fitting parameters of a disk structure model:

1. model specification
2. radiative synthesis
3. visibility interpolation
4. likelihood evaluation

I have organized the code in this package into source (in the `src` directory) and scripts (in the `scripts`) directory. After installation (see below), if you do

    using JudithExcalibur

You will be able to access any of the modules within `src`. These provide the base functionality for the package in common tasks, such as reading in images produced by RADMC-3D, Fourier transforms, and visibility interpolation. My hope is to make these components as general as possible so if you would like to extend this package to fit a novel type of disk, it will be easy to reuse many of the core functionality.

The scripts are not technically part of the Julia package (you cannot import them like the modules) but instead provide "driver scripts" that utilize the core modules to address a certain research question. For example, `sythesize_model.jl` synthesizes and a model for your current parameters. These are run from your system shell after adding them to your `PATH`

    $ sythesize_model.jl

## Computational demand

Because spectral line datasets are large, and synthesizing models is computationally expensive, I have designed this package to work in a parallel environment, where the calculations for each frequency channel can be distributed to independent processors. Therefore, please keep this architecture in mind when navigating the source code. Due to the computationally expensive nature of the radiative synthesis, fitting sizable datasets (e.g., SMA and ALMA) will *require a moderately sized number of CPU cores for computation in a reasonable timeframe*. For example, to fully explore the posterior for the AK Sco example dataset will require a few days of computing time on ~32 cores or more.

## Installation

First, make sure you have installed RADMC-3D and you can successfully run one of the example scripts contained within this package. My current RADMC-3D installation is v0.38, although I expect that everything should work well on other recent versions of the package. Second, I have now migrated this package to run on Julia v0.4+, and I plan to keep it up to date with current Julia releases.

Because `JudithExcalibur` is not an official Julia package (nor will it likely be), for now, installation involves simply cloning the repository. First, open up a Julia prompt in the REPL, then type

    julia> Pkg.clone("https://github.com/iancze/JudithExcalibur.git")

As mentioned previously, there are several "driver" command line scripts that are used to perform the actual mass fitting. To complete the installation, you should add these files to your PATH. To figure out where the package is installed

    julia> Pkg.dir("JudithExcalibur")
    "/home/ian/.julia/JudithExcalibur"

Your PATH will vary. The scripts are located inside of the `scripts` directory, so if you are using bash or Z-shell, you will want to add the PATH that looks something like

    export PATH="/home/ian/.julia/JudithExcalibur/scripts:$PATH"

inside of your `.bashrc` or `.zshrc` file. Finally,

    $ source ~/.zshrc

To check that you have properly added the scripts, you can try in your system shell

    $ JudithInitialize.jl --test
    Your JudithExcalibur scripts are successfully linked.
    Exiting

Second, most of the scripts and IO routines require Python and several Python packages. Please install the following packages via your own package manager or a distribution like anaconda:

* numpy
* scipy
* matplotlib
* Jupyter/IPython
* h5py
* astropy
* Jupyter/IJulia

With the package successfully installed, see the documentation in the `docs/` folder on how to get started fitting a specific source, in particular the [cookbook](docs/cookbook.md).

For those interested in the source code, the most important files to start browsing are

**model.jl**: Contains the actual specification of the parametric disk model, as well as the tools to write to disk the synthesis files RADMC-3D requires.

**image.jl**: Contains type definitions to read images produced by RADMC-3D, as well as convert from physical coordinates to sky coordinates.

**visibilities.jl**: Contains type definitions to hold the dataset and the model visibilities. Additionally contains functions to apply phase shifts to the visibilities corresponding to shifts in the image plane. Also contains functions to FFT images to the visibility plane.

**gridding.jl**: Contains the prolate-spheroidal wave function definitions from Schwab 1984, used when doing the visibility interpolations.

**venus.jl**: This implementation uses the Ensemble Sampler (a Julia port from Dan Foreman-Mackey's `emcee` python package) to sample the posterior distribution using parallelized walkers.
