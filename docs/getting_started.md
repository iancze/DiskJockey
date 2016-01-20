# AK Sco Cookbook

After you have installed both `RADMC-3D` and `JudithExcalibur`, it's time to get up and running!

## Initialization

Now, the `JudithExcalibur` package should be successfully installed user-wide. Because it is likely that you will want to fit more than just one protoplanetary disk, or perhaps try different model specifications for a particular disk, the code structure is organized so that you will have a separate directory for each run. As an example to get you started, we have included the AK Sco dataset from ALMA. You can download the dataset in HDF5 format [here](https://figshare.com/articles/ALMA_AK_Sco_12CO_J_2_1_Visibilities/2066022).

    $ mkdir AKSco
    $ cd AKSco

Once you've created a working directory for your project, then you'll want to initialize this directory with a config file

    $ JudithInitialize.jl --new-project=standard
    Copied default config.yaml file to current working directory.
    Exiting

Now, open up `config.yaml` with your favorite text editor and change the fields as you see fit, including which transition of CO you would like to fit. For now, I have only included functionality for 12CO, 13CO, and C18O in LTE. Please create an issue on github if you would like a new species added. In particular, I recommend checking to see if the transition and species are correct.

## Synthesizing a model image

Next, we'll need to initialize the directory with the appropriate RADMC-3D files specific to your dataset.

    $ JudithInitialize.jl

This command copies the appropriate molecular data to your current directory and writes out the input files for RADMC-3D, using the parameters defined in `config.yaml`.

Before starting running MCMC chains, it is a good idea to take a guess at some parameters and try synthesizing a set of model channel maps. Then, you can play around with values in `config.yaml` to see what works well.

Now, synthesize a model using RADMC-3D:

    $ synthesize_model.jl

During this process, the output from RADMC-3D is piped to `STDOUT`. This may be a good place to debug if anything looks fishy. Note that if you make changes to `config.yaml`, you'll need to rerun `JudithInitialize.jl` to recreate the input files for RADMC-3D before running `synthesize_model.jl`, otherwise you'll be synthesizing stale input files.

Next, you may like to plot up channel maps of the synthesized image and a spatially-integrated spectrum. This is done via

    $ plot_chmaps.jl

And, if you'd like to plot up a zeroth-moment map

    $ plot_moments.jl

For me, a typical workflow for playing around with channel maps is

1. Edit `config.yaml` to parameters that might make sense
2. at the command line, run

    $ JudithInitialize.jl && synthesize_model.jl && plot_chmaps.jl && plot_moments.jl

the `&&` ensures that the previous command completes before moving on to the next command.
3. Inspect the resulting plots of the data, and if I am not satisfied go back to 1.

It is a very good idea to inspect your channel maps to make sure that there isn't any weird structure, that you have enough pixels to resolve the disk structure, and that your model grid appears to be at high enough resolution. A few extra minutes or hours spent debugging your images during this step can save you days (of supercomputer time) in the steps ahead.

## Setting up a parallelized MCMC exploration of the parameters

Now, open up `InitilializeWalkers.ipynb` with your Jupyter notebook. To save you some computational time burning in the walkers, I found that for AK Sco, the following walker starting positions worked well.

    p0 = np.array([np.random.uniform(low=1.03, high=1.05, nwalkers), # mass
              np.random.uniform(20., 21.0, nwalkers), #r_c
              np.random.uniform(110., 115, nwalkers), #T_10
              np.random.uniform(0.70, 0.71, nwalkers), # q
              np.random.uniform(-3.4, -3.5, nwalkers), #log10 M_gas
              np.random.uniform(0.17, 0.18, nwalkers), #xi
               np.random.uniform(144.0, 145.0, nwalkers), #dpc
              np.random.uniform(159.0, 160.0, nwalkers), #inc
              np.random.uniform(40.0, 41.0, nwalkers), #PA
              np.random.uniform(-0.1, 0.1, nwalkers), #vz
              np.random.uniform(-0.1, 0.1, nwalkers), #mu_a
              np.random.uniform(-0.1, 0.1, nwalkers)]) #mu_d

The following is a description of the most important files in the package

**model.jl**: Contains the actual specification of the parametric disk model, as well as the tools to write to disk the synthesis files RADMC-3D requires.

**image.jl**: Contains type definitions to read images produced by RADMC-3D, as well as convert from physical coordinates to sky coordinates.

**visibilities.jl**: Contains type definitions to hold the dataset and the model visibilities. Additionally contains functions to apply phase shifts to the visibilities corresponding to shifts in the image plane. Also contains functions to FFT images to the visibility plane.

**gridding.jl**: Contains the prolate-spheroidal wave function definitions from Schwab 1984, used when doing the visibility interpolations.

**venus.jl**: This implementation uses the Ensemble Sampler (a Julia port from Dan Foreman-Mackey's `emcee` python package) to sample the posterior distribution using parallelized walkers.
