# Changelog

The following are the changes that have been implemented since the previous version. All versions are tagged as releases.

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
