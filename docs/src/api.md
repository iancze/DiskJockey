# API

For those interested in the source code, the most important files to start browsing are

**model.jl**: Contains the specification of the parametric disk model, as well as the tools to write to disk the synthesis files RADMC-3D requires.

**image.jl**: Contains type definitions to read images produced by RADMC-3D, as well as convert from physical coordinates to sky coordinates.

**visibilities.jl**: Contains type definitions to hold the dataset and the model visibilities. Additionally contains functions to apply phase shifts to the visibilities corresponding to shifts in the image plane. Also contains functions to FFT images to the visibility plane.

**gridding.jl**: Contains the prolate-spheroidal wave function definitions from Schwab 1984, used when doing the visibility interpolations.

**venus.jl**: This implementation uses the Ensemble Sampler (a Julia port of Dan Foreman-Mackey's `emcee` python package) to sample the posterior distribution using parallelized walkers.


# Constants

```@docs
constants.fftspace(width::Real, N::Int)
```

# model.jl

```@autodocs
Modules = [model]
```

# image.jl

```@autodocs
Modules = [image]
```

# visibilities.jl

```@autodocs
Modules = [visibilities]
```

# gridding.jl

```@autodocs
Modules = [gridding]
```

# gauss.jl

```@autodocs
Modules = [gauss]
```

# EnsembleSampler.jl

This is a Julia port of the Python package `emcee`, written by Dan Foreman-Mackey and collaborators.

```@autodocs
Modules = [EnsembleSampler]
```


## Index
```@index
```
