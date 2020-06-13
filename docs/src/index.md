# DiskJockey.jl Documentation

`DiskJockey.jl` is a Julia package designed for forward-modeling radio observations of the dust and gas in protoplanetary disks using interferometers like the Submillimeter Array (SMA) and the Atacama Large Millimeter/Submillimeter Array (ALMA). Its main features consist of utilities to generate common protoplanetary disk structures, generate channel maps (using radiative transfer via `RADMC-3D`), Fourier transform and sample these at the baselines of the observations, and evaluate a visibility-based likelihood function. On top of this, there is a Julia version of the Affine Invariant Ensemble sampler to perform Markov Chain Monte Carlo Sampling of the posterior distribution.

# Citation 

**If you use this code or a derivative of it in your research, we would really appreciate it if you cited the first paper in our series, [Czekala et al. 2015 ApJ, 806 154C](http://adsabs.harvard.edu/abs/2015ApJ...806..154C).**

# Papers using DiskJockey

* *A Disk-based Dynamical Constraint on the Mass of the Young Binary AK Sco* : [Czekala et al. 2015 ApJ, 806 154C](http://adsabs.harvard.edu/abs/2015ApJ...806..154C)
* *A Disk-based Dynamical Constraint on the Mass of the Young Binary DQ Tau* : [Czekala et al. 2016 ApJ, 818 156C](http://adsabs.harvard.edu/abs/2016ApJ...818..156C)
* *ALMA Measurements of Circumstellar Material in the GQ Lup System* : [Macgregor et al. 2017, ApJ, 835, 17M](http://adsabs.harvard.edu/abs/2017ApJ...835...17M)
* *ALMA Observations of the Young Substellar Binary System 2M1207* : [Ricci et al. 2017 AJ, 154, 24R](http://adsabs.harvard.edu/abs/2017AJ....154...24R)
* *The Architecture of the GW Ori Young Triple Star System and Its Disk: Dynamical Masses, Mutual Inclinations, and Recurrent Eclipses* : [Czekala et al. 2017, ApJ, 851, 132](https://ui.adsabs.harvard.edu/abs/2017ApJ...851..132C/abstract)
* *The Degree of Alignment between Circumbinary Disks and Their Binary Hosts*: [Czekala et al. 2019, ApJ 883, 1, 22](https://ui.adsabs.harvard.edu/abs/2019ApJ...883...22C/abstract)