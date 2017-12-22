# DiskJockey.jl Documentation

`DiskJockey.jl` is a Julia package designed for forward-modeling radio observations of the dust and gas in protoplanetary disks using interferometers like the Submillimeter Array (SMA) and the Atacama Large Millimeter/Submillimeter Array (ALMA). Its main features consist of utilities to generate common protoplanetary disk structures, generate channel maps (using radiative transfer via `RADMC-3D`), Fourier transform and sample these at the baselines of the observations, and evaluate a visibility-based likelihood function. On top of this, there is a Julia version of the Affine Invariant Ensemble sampler to perform Markov Chain Monte Carlo Sampling of the posterior distribution.

As a first step, I recommend reading [Understanding disk-based dynamical mass measurements](@ref).
