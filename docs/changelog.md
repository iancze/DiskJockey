# Changelog

The following are the changes that have been implemented since the previous version. All versions are tagged as releases.

# Version 0.1.2

## Installation

Now provides installation instructions for tagged releases. RADMC-3D is installed automatically as a local version.

## Travis integration

Now the builds are tested on travis-ci, which should hopefully increase stability for future versions.

## fake_dataset.jl

Creation of fake data for plotting and testing purposes.

## EnsembleSampler

Now takes an optional function, to be called at the end of each loop as `func(sampler, outdir)`

    function run_schedule(sampler::Sampler, pos0, N::Int, loops::Int, outdir, func::Function=nothing)

See the `src/EnsembleSampler.jl` file for more details. This is primarily in support of the next feature...

## Plotly script

If you run the sampling script as

    $ venus.jl --plotly

After each loop, creates/updates a plotly walkers plot corresponding to the `name` entry in `config.yaml`. This allows easy monitoring of many different chains that might be running on a cluster.

## JudithInitialize.jl

Now can easily spit out an `exclude` array so that fewer channels can be fit during initial testing.

## Cavity model

By initializing a directory with

    $ JudithInitialize.jl --new-project=cavity

you can start exploring the `cavity` model, which has an exponential taper inside of some radius, `r_cav`.

## Overwrites

Now `venus.jl` should guard against overwriting previous MCMC samples in the case of a job restart on HPC. This requires changing

    run_schedule() to include a write to chain.npy file.

Anyway, not actually implemented yet.

## Channel maps

plots a circle over the image at the radius of the grid. If emission gets close, you're in trouble.

## gridding.jl

`gridding.jl` now exports a `corrfun(img::SkyImage)` routine, in contrast to `corrfun!(img::SkyImage)`. This new function returns a corrected image as a copy, leaving the original image in the arguments unchanged. This is useful for plotting and debugging scripts so that you don't need to copy the image manually.

## visibilities.jl

`FullModelVis` now has a `-` method, which can be used to subtract two sets of dense visibilities.

## image.jl

`Image` now has a `-` method, which can be used to subtract two images.

## config.yaml copied to output directory

Now `venus.jl` will move a copy of your `config.yaml` file to the output directory. This creates an automatic record of what parameters you ran with, which will undoubtedly be useful when reviewing previous runs sometime in the near future.

Also, default values in the initial `config.yaml` files have been updated.

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
