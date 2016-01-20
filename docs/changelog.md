# Changelog

The following are the changes that have been implemented since the previous version. All versions are tagged as releases.

# Version 0.1.1

Now includes a cookbook example for AK Sco.

New implemenation of models. Created the `AbstractParameters` type, with subtypes `ParametersStandard` and `ParametersTruncated`. Added new parameters type

Note that model types `cavity` and `truncated` are still experimental and likely to change.

Updated the `config.yaml` format to include new fields like `model: standard`. Changed the items in the `parameters` dictionary to be a single Float64 number, not an array of `[starting, jump]` like it was previously. An example is `assets/config.yaml`.

Simplified synthesizing and plotting of models. Split `plot_model.jl` into `synthesize_model.jl` and separately `plot_chmaps.jl`.

Added new `plot_moments.jl` script which plots zeroth-moment.

Created `InitializeWalkers.ipynb` that is copied to new directory to help specify walker starting positions.

Streamlined various model implementations into the `venus.jl` code.

`model.jl` now handles the implementation of priors. This greatly streamlines the function `fprob` in `venus.jl`.



# Version 0.1.0

Initial commit.
