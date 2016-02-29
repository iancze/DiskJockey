# Changelog

The following are the changes that have been implemented since the previous version. All versions are tagged as releases and available as a clone off of the master branch.

# Version 0.1.3

## Visualizing the tau=1 surface

The functionality originally stored in RADMC is now exposed via DiskJockey to visualize the tau = 1 (or other value) surface.

  $ DJinit.jl && tausurf_model.jl && plot_tausurf.jl

This also introduces the `TausurfImage` type in `src/image.jl`. This is useful to see whether the tau=1 surface is in front of or behind the midplane of the disk, projected on the sky.

The more useful plot is in 3D, however.

## Vertical temperature gradient

Now includes support for fitting a vertical temperature gradient, following the parameterization in Williams and Best 14.

## model.convert_vector

Has been simplified to ingest *kwargs*, enabling easier selection of model types and parameters.

# Version 0.1.2

The package is now named DiskJockey (previously JudithExcalibur). Check out the new logo!

## Installation

The README now provides installation instructions for tagged releases. Rather than requiring the user to install it separately, RADMC-3D is now installed automatically as part of the installation process.

## UVHDF5

Conversion from UVFITS and CASA measurement set was previously handled by scripts in this repository. However, this import and export capability has be standardized via the [UVHDF5](https://github.com/Astrochem/UVHDF5) package. Please see that package for any issues with import and export.

## Travis integration

Now the builds are tested on travis-ci, which should hopefully increase stability for the development cycle.

## EnsembleSampler

There is now an `expand_walkers.py` script, designed to add the `dpc` dimension to the sampling routine from an ensemble of walkers used on a posterior with distance fixed.

Now the ensemble sampler takes an optional function, to be called at the end of each loop as `func(sampler, outdir)`

    function run_schedule(sampler::Sampler, pos0, N::Int, loops::Int, outdir, func::Function=nothing)

See the `src/EnsembleSampler.jl` file for more details. This is primarily in support of the next feature...

## Plotly script

If you run the sampling script as

    $ venus.jl --plotly

Then after the completion of each sampling loop, this will create a plotly walkers plot corresponding to the `name` entry in `config.yaml`. This allows easy monitoring of many different chains that might be running on a cluster.

## DJInitialize.jl

(Previously JudithInitialize.jl). Now this initialization script can easily spit out an `exclude` array so that fewer channels can be fit during initial testing.

## Cavity model

By initializing a directory with

    $ DJInitialize.jl --new-project=cavity

you can start exploring the `cavity` model, which has an exponential taper inside of some radius, `r_cav`.

## gridding.jl

`gridding.jl` now exports a `corrfun(img::SkyImage)` routine, in addition to `corrfun!(img::SkyImage)`. This new function returns a corrected image as a copy, leaving the original image in the arguments unchanged. This is useful for plotting and debugging scripts so that you don't need to copy the image manually.

## visibilities.jl

`FullModelVis` now has a `-` method, which can be used to subtract two sets of dense visibilities.

## image.jl

`Image` now has a `-` method, which can be used to subtract two images.

## config.yaml copied to output directory

Now `venus.jl` will move a copy of your `config.yaml` file to the output directory. This creates an automatic record of what parameters you ran with, which will undoubtedly be useful when reviewing previous runs sometime in the near future. Default values in the initial `config.yaml` files have also been updated.

## plot_baselines.jl

Will generate a plot of where your dataset has sampled the UV plane. Also will print out an average velocity of the dataset, to allow a good starting guess for the `vel` parameter.

# Version 0.1.1

## Model specification

New implementation of models through parameter types in `model.jl`. Instead of separate model files for each new parameterization, which duplicated a lot of code and made things difficult to maintain, we are now collating all model types in `model.jl` by their parameters. We have created the `AbstractParameters` umbrella type, with subtypes `ParametersStandard`, `ParametersTruncated`, and `ParametersCavity`. Previous model specification routines like `Sigma(r::Float64, pars::Parameters)` (surface density) are now have overloaded methods for different models based upon these parameter types

    Sigma(r::Float64, pars::ParametersStandard)
    Sigma(r::Float64, pars::ParametersTruncated)
    Sigma(r::Float64, pars::ParametersCavity)

And routines that are general to all models (e.g. velocity specification) are denoted by

    velocity{T}(r::T, pars::AbstractParameters)

Note that model types `cavity` and `truncated` are still experimental and likely to change.

To go along with this change, I have updated the automatic generation of the `config.yaml` file to include new fields like `model: standard`. Within the `config.yaml` file, items in the `parameters` dictionary are now simply a single Float64 number, not an array of `[starting, jump]` like it was previously.

`model.jl` now handles the implementation of priors, that dispatch off of the parameter types. This greatly streamlines the function `fprob` in `venus.jl`.


## One-off model synthesis and plotting

I've simplified synthesizing and plotting of models. What was previously `plot_model.jl` is now split into `synthesize_model.jl` and `plot_chmaps.jl`.

Added new `plot_moments.jl` script which plots the zeroth-moment image. Soon it will plot first moment as well.

## MCMC Sampling

The specification of parameter types allowed me to greatly simplify the MCMC sampling code into a single script, `venus.jl`. I have moved previous sampling scripts to the `attic/` directory.

Created `InitializeWalkers.ipynb` that is copied to new directory to help specify walker starting positions. The user edits this with a Jupyter/IPython notebook.

## Cookbook

We now have a [cookbook example](docs/cookbook.md) for AK Sco, check it out to get started!

`plot_walkers.py` now includes ability to determine highest density interval for quoting credible intervals and should automatically label the parameters after reading from `config.yaml`.

# Version 0.1.0

Initial commit on the new versioning roadmap.
