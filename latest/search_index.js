var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#DiskJockey.jl-Documentation-1",
    "page": "Home",
    "title": "DiskJockey.jl Documentation",
    "category": "section",
    "text": "DiskJockey.jl is a Julia package designed for forward-modeling radio observations of the dust and gas in protoplanetary disks using interferometers like the Submillimeter Array (SMA) and the Atacama Large Millimeter/Submillimeter Array (ALMA). Its main features consist of utilities to generate common protoplanetary disk structures, generate channel maps (using radiative transfer via RADMC-3D), Fourier transform and sample these at the baselines of the observations, and evaluate a visibility-based likelihood function. On top of this, there is a Julia version of the Affine Invariant Ensemble sampler to perform Markov Chain Monte Carlo Sampling of the posterior distribution.As a first step, I recommend reading Understanding disk-based dynamical mass measurements."
},

{
    "location": "dynamical_mass_intro.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "dynamical_mass_intro.html#Understanding-disk-based-dynamical-mass-measurements-1",
    "page": "Introduction",
    "title": "Understanding disk-based dynamical mass measurements",
    "category": "section",
    "text": "Mass is the fundamental property that determines a star's fate. Compared to the masses of their older cousins on the main sequence, however, the masses of young pre-main sequence stars are relatively uncertain. Fortunately, young stars have a unique advantage that we can exploit to precisely infer their mass: many often host a protoplanetary disk made of gas and dust, the site of future and ongoing planet formation. By modeling the rotation of this disk, we can dynamically \"weigh\" the host star(s). Circumstellar disks are typically mm-bright, which means that interferometers like the Atacama Large Millimeter Array (ALMA) and the Submillimeter Array (SMA) are ideal instruments for deriving stellar masses.In this short intro, I will explain how we can use the kinematic fingerprint of the disk–imprinted by disk rotation–to derive precise estimates of the stellar mass."
},

{
    "location": "dynamical_mass_intro.html#Disk-geometry-and-rotation-1",
    "page": "Introduction",
    "title": "Disk geometry and rotation",
    "category": "section",
    "text": "Astronomers typically model the structure of these protoplanetary disks as an axisymmetric disk with radial surface density defined by a power law and vertical structure set by hydrostatic equilibrium. As we will see, the orientation of the disk relative to us has a profound effect on what we observe. We define the inclination of the disk i as the angle between the disk angular momentum axis and the observer: 0^circ inclination means the disk is viewed face-on (appearing on the sky as a circle), while 90^circ inclination means the disk is viewed edge-on (appearing on the sky as a thin line).(Image: disk)The disk orbits the central stellar mass under the force of gravity, yielding velocities dictated by Kepler's law. Although there is evidence that due to pressure support gas disks rotate at slightly sub-Keplerian speeds, creating a headwind for any larger dust aggregates and solids (the meter-sized barrier), this effect is negligible when determining stellar mass. For more information about disk structure, check out the lecture notes by Phil Armitage."
},

{
    "location": "dynamical_mass_intro.html#Disk-emission-mechanisms-1",
    "page": "Introduction",
    "title": "Disk emission mechanisms",
    "category": "section",
    "text": "Active young stars and their circumstellar disks can emit radiation at a host of wavelengths from the X-ray to the radio. When observing at sub-mm and radio wavelengths, there are two main sources of emission: Thermal continuum emission from the dust and spectral line emission from molecular species like carbon monoxide (CO).Owing to the high spatial resolution of an interferometer, it is possible to spatially resolve many circumstellar disks. If we set the frequency sampling finely enough, then we will have a dataset that is both spatially resolved (an image showing emission on many scales) and spectrally resolved (an image for every frequency we observe). This set of images is called a data cube or a series of channel maps."
},

{
    "location": "dynamical_mass_intro.html#Dust-emission-1",
    "page": "Introduction",
    "title": "Dust emission",
    "category": "section",
    "text": "The intensity of continuum emission from the dust is related to the density and temperature of the emitting material. Because dust emission spans a continuum of frequencies, channel maps with fine frequency spacing yield little additional information–-typically these are summed together to yield greater sensitivity. Below is a movie showing what the dust emission looks like as we vary some of the key disk parameters: the radius of the disk, the inclination of the disk, and the stellar mass (stellar flux is kept fixed). In the movie, each parameter starts at an intermediate value, is then decreased to a minimum value, then increased to a maximum value, and then decreased back to the starting value. In this movie, the changes you will notice are due to variations in the radius and inclination of the disk.(Image: dust video)"
},

{
    "location": "dynamical_mass_intro.html#Spectral-line-emission-1",
    "page": "Introduction",
    "title": "Spectral line emission",
    "category": "section",
    "text": "When we observe the disk in narrow frequency channels at the location of spectral line emission, like CO J=2-1, we see the kinematic fingerprint of the disk. Each image shows the disk emission at a specific frequency, corresponding to the blue-shift (-velocity) or red-shift (+velocity) of the CO line, due to the projected velocity along the line of sight.Why the morphology of the line emission appears as it does isn't immediately obvious. Primarily, the shape of the emission is a function of the Keplerian rotation: the inner part of the disk rotates faster than the outer part of the disk and thus imparts a more substantial velocity shift to the line emission from these regions of the disk.As the parameters of the disk change, the projected velocity can change dramatically. For example, compare the difference between i approx 0^circ and i approx 90^circ! If the disk is more inclined towards face on (0^circ), the projection of the velocity along the line of sight will reduced and we will see more emission concentrated around 0 km/s.(Image: Gas Video)"
},

{
    "location": "dynamical_mass_intro.html#Interferometers-measure-complex-visibilities-1",
    "page": "Introduction",
    "title": "Interferometers measure complex visibilities",
    "category": "section",
    "text": "Although we may be used to looking at emission in the image plane, interferometers measure the Fourier transform of the sky brightness in the visibility plane. This next movie shows what the Fourier transform of each channel map looks like. As you can see, the Fourier transform looks quite different from the image.(Image: Visibilities Video)This Fourier plane is sampled according to the baselines formed between antennae separations. Below is a typical sampling pattern for an observation with the SMA.(Image: UV Spacings)Every dot represents a sampling of the Fourier transform in the UV plane. Some of the samplings are spaced so close together in space that they smear together into a line."
},

{
    "location": "models.html#",
    "page": "Models",
    "title": "Models",
    "category": "page",
    "text": ""
},

{
    "location": "models.html#Models-1",
    "page": "Models",
    "title": "Models",
    "category": "section",
    "text": "DiskJockey supports many different disk models, and is designed in such a way so as to be easily extensible to any new model you might like (pull requests welcome)!The current models include"
},

{
    "location": "models.html#Standard-1",
    "page": "Models",
    "title": "Standard",
    "category": "section",
    "text": "Keyword: standardThis is a standard model that is used in Czekala et al. 2015 and 2016, and Rosenfeld et al. 2012."
},

{
    "location": "models.html#Cavity-1",
    "page": "Models",
    "title": "Cavity",
    "category": "section",
    "text": "Keyword: cavityThis model includes an adjustable exponential taper for the inner surface density profile. This is meant to mimic the appearance of large inner gas gaps in disks. Relevant notebook sketches of this idea can be found in notebooks/Cavity Surface Density.ipynb. The main difference from the standard model is the introduction of two new parameters: a cavity radius and decay profile exponent."
},

{
    "location": "models.html#Vertical-1",
    "page": "Models",
    "title": "Vertical",
    "category": "section",
    "text": "Keyword: verticalThis model includes a more realistic temperature profile as outlined in Dartois et al. 2003, Rosenfeld et al. 2013, and Williams and Best 2014, among many others. We follow the specific parameterization in Williams and Best 2014.The primary additional step for this model is the necessity to numerically solve the hydrostatic equilibrium equation (WB14 Eqn 1).for a given radius, r, come up with a vertical grid of z points that stretches from the midplane to an appropriate height where the gas density is approximately zero.\ncalculate ln(rho(r,z)) as a function of (r, z) by integrating WB14 Eqn 1 using quadgk. This is the unnormalized density of the disk, which means the relative values at a fixed radius should be correct, but no guarantees about anything else. We will call the unnormalized density un_rho(r,z).\nnow, we know that Sigma(r) specifies the surface density of the disk, i.e., or the integral of rho(r, z) from z = -inf to +inf. Therefore, we can find rho(r,z) by vertically integrating un_rho(r,z) and finding the normalization constant, norm(r). Because we are dealing with integrals of very big and very small numbers, we need to do some tricks to avoid overflow and underflow errors.\nfinally, we know that CO is photodissociated when it is not shielded by enough gas, i.e., there is not enough gas column density above this height to block harmful radiation. Therefore, we need to find the photodissociation height, z_phot(r), above which the molecule CO can no longer exist. To do this requires integrating from z = +inf towards the midplane (z = 0) to find the height at which we have accumulated enough gas column density to be shielded.This sounds like a lot of steps just to evaluate a single rho(r,z) point. Because RADMC-3D solves the radiative transfer on a spherical grid and the disk model is defined on a cylindrical grid, there is a careful order of operations necessary to achieve the appropriate accuracy in the shortest amount of computational time. We address this by first solving everything on a cylindrical grid to find norm(r) and z_phot(r) as a function of disk radius. Then, for a given (r_spherical, theta) point, we convert to cylindrical coordinates and solve WB14 Eqn 1 to find un_rho(r,z)."
},

{
    "location": "models.html#VerticalEta-1",
    "page": "Models",
    "title": "VerticalEta",
    "category": "section",
    "text": "Characterized with ParametersVerticalEta. In addition to the parameters described in the Vertical model, this extension has an additional parameter eta, designed to vary the height of the atmosphere with radius."
},

{
    "location": "models.html#Conventions-1",
    "page": "Models",
    "title": "Conventions",
    "category": "section",
    "text": ""
},

{
    "location": "models.html#Inclination-1",
    "page": "Models",
    "title": "Inclination",
    "category": "section",
    "text": "Disk inclination ranges from 0 to 180 degrees. 0 degrees means face on, angular momentum vector pointing at observer; 90 means edge on; and 180 means face on, angular momentum vector pointing away from observer. These are the same as the RADMC-3D conventions."
},

{
    "location": "models.html#Position-Angle-1",
    "page": "Models",
    "title": "Position Angle",
    "category": "section",
    "text": "We also adopt the RADMC-3D convention for position angle, which defines position angle by the angular momentum vector. A positive PA angle means the disk angular momentum vector will be rotated counter clockwise (from North towards East)."
},

{
    "location": "formats.html#",
    "page": "Formats",
    "title": "Formats",
    "category": "page",
    "text": ""
},

{
    "location": "formats.html#File-Formats-1",
    "page": "Formats",
    "title": "File Formats",
    "category": "section",
    "text": "Interferometric data from reduction software like CASA or MIRIAD is generally stored in measurement sets (*.ms) or FITS files (*.fits), and usually contains a lot of ancillary data unnecessary for a dynamical mass measurement (or fitting in the UV plane in general). To reduce dependency on these outside programs and overall memory footprint, we have made it so that all of the Julia code only interfaces with the minimal necessary visibility data stored in an HDF5 file.The key Julia code used to read this file is provided in /src/visibilities.jl which reads the HDF5 file and stores each channel in an array of DataVis instances. For example, this data.hdf5 file contains the following datasets in arrays of (nrows, ncols) form. nchan is the number of channels in the dataset, nvis is the number of complex visibilities in the measurement set.The format of this datafile is as follows:data.hdf5\n  lams # [μm] (nchan) Wavelength (in microns) corresponding to channel\n  uu # [kλ] (nchan, nvis) Vectors of the u locations in kilolambda\n  vv # [kλ] (nchan, nvis) Vectors of the v locations in kilolambda\n  real # [Jy] (nchan, nvis) real component of the complex visibilities\n  imag # [Jy] (nchan, nvis) imaginary component of the complex visibilities\n  invsig # [1/Jy] (nchan, nvis) the inverse of sigma for each visibility (1/sigma)and an example datafile is provided for AK Sco, available for download here. The invsig field may seem like an non-traditional way to store visibility weights (typically measured in [1/Jy^2]), but enables quick computation of the chi^2 statistic via Julia's sum(abs2,...) routine.Example scripts for reading in data from the SMA FITS format are available in the scripts/ directory. We are working on a more centralized way to easily extract ALMA visibilities from a measurement set. Right now it tends to be a bit object-specific."
},

{
    "location": "RADMC3D_setup.html#",
    "page": "RADMC Setup",
    "title": "RADMC Setup",
    "category": "page",
    "text": ""
},

{
    "location": "RADMC3D_setup.html#RADMC3D-setup-1",
    "page": "RADMC Setup",
    "title": "RADMC3D setup",
    "category": "section",
    "text": "This is a brief overview of some of the options that are set for RADMC-3D."
},

{
    "location": "RADMC3D_setup.html#Gridding-1",
    "page": "RADMC Setup",
    "title": "Gridding",
    "category": "section",
    "text": "Most models are specified (mathematically) as 2D axisymmetric in cylindrical coordinates. RADMC-3D takes files in spherical coordinates, and so this conversion is done in model.jl as the model is written out.Generally, the grid cells are logarithmically spaced in radius. The elevation (theta) cells are logarithmically spaced in theta, with cells clustered more denesly near the midplane than at the poles. Because of the axisymmetry, there is only one phi (azimuth) cell.Although RADMC-3D has the ability to determine dust temperatures through radiative Monte Carlo, we set the gas temperatures directly though a simple radial temperature profile.Then, the user also chooses the number of pixels npix for the resulting image to have. A good rule of thumb is that pixels need to be roughly 5-10x smaller than the resolution of the interferometric observations."
},

{
    "location": "RADMC3D_setup.html#Files-required-for-line-transfer-and-what-they-are-1",
    "page": "RADMC Setup",
    "title": "Files required for line transfer and what they are",
    "category": "section",
    "text": "lines.inp: control file for line transfer\nmolecule_co.inp contains properties of atom or molecule\nnumberdens_co.inp number density of molecule in units of cm^{-3}\ngas_temperature.inp gas temperature at each grid cell"
},

{
    "location": "RADMC3D_setup.html#radmc.inp-1",
    "page": "RADMC Setup",
    "title": "radmc.inp",
    "category": "section",
    "text": "tgas_eq_tdust=0 #this means we will specify the gas temperature at each grid cell, rather than setting it equal to the dust temperature."
},

{
    "location": "RADMC3D_setup.html#lines.inp-1",
    "page": "RADMC Setup",
    "title": "lines.inp",
    "category": "section",
    "text": "Format styles1 #format style 1 (may be 2)\n1\nco    leiden    0    0This means that we will only be using the CO molecule. All of the information about this molecule is stored in the file molecule_co.inp. This is read from the Leiden database. Basically it provides all of the quantum information about the energy levels."
},

{
    "location": "RADMC3D_setup.html#numberdens_XXX.inp-1",
    "page": "RADMC Setup",
    "title": "numberdens_XXX.inp",
    "category": "section",
    "text": "The number density of each species per cubic centimeter."
},

{
    "location": "RADMC3D_setup.html#gas_temperature.inp-1",
    "page": "RADMC Setup",
    "title": "gas_temperature.inp",
    "category": "section",
    "text": "Specifies the gas temperature at each cell."
},

{
    "location": "RADMC3D_setup.html#gas_velocity.inp-1",
    "page": "RADMC Setup",
    "title": "gas_velocity.inp",
    "category": "section",
    "text": "Three numbers at each grid point, which are the components of velocity along each direction. All components have units of cm/s."
},

{
    "location": "RADMC3D_setup.html#microturbulence.inp-1",
    "page": "RADMC Setup",
    "title": "microturbulence.inp",
    "category": "section",
    "text": "We can also specify a microturbulent broadening for each cell."
},

{
    "location": "RADMC3D_setup.html#partitionfunction_XXX.inp-1",
    "page": "RADMC Setup",
    "title": "partitionfunction_XXX.inp",
    "category": "section",
    "text": "This can be specified for CO or other species. Otherwise RADMC3D will calculate the partition function itself using the molecular data file."
},

{
    "location": "RADMC3D_setup.html#amr_grid.inp-1",
    "page": "RADMC Setup",
    "title": "amr_grid.inp",
    "category": "section",
    "text": "We can use different coordinate systems here.**coordsystem**:\n< 100: cartesian\n100<= - <200: spherical\n200<= - <300: cylindricalIf we want to make a spherical grid, using 2D symmetry, we can specify that we wish one dimension to be non-active.The files that need to be generated by our package areamr_grid.inp: location of grid cells\nmicroturbulence.inp: microturbulence at each cell\ngas_velocity.inp: three component velocity vector at each cell\nnumberdens_co.inp: the number density of the CO molecule at each cellAll of the other files can be pre-generated or lifted from some database."
},

{
    "location": "cookbook.html#",
    "page": "Cookbook",
    "title": "Cookbook",
    "category": "page",
    "text": ""
},

{
    "location": "cookbook.html#Cookbook-1",
    "page": "Cookbook",
    "title": "Cookbook",
    "category": "section",
    "text": "After you have installed DiskJockey, it's time to get up and running!"
},

{
    "location": "cookbook.html#Initialization-1",
    "page": "Cookbook",
    "title": "Initialization",
    "category": "section",
    "text": "Now the DiskJockey package should be successfully installed user-wide. Because it is likely that you will want to fit more than just one protoplanetary disk, or perhaps try different model specifications for a particular disk, the code structure is organized so that you will have a separate directory for each disk model. The following is an example to get you started fitting AK Sco.$ mkdir AKSco\n$ cd AKScoMake sure to download the dataset in HDF5 format here.Now, you'll want to initialize this directory with a config file. This config file will store all of the options that are specific to fitting this disk and is frequently used by many of the scripts in this package. To initialize,$ DJ_initialize.jl --new-project=standard\nCopied default config.yaml, InitializeWalkers.ipynb, and Makefile for the standard model to current working directory.\nExiting--new-project also has other than standard, such as cavity and vertical in order to fit more exotic models. Now, open up config.yaml with your favorite text editor and change the fields as you see fit, including which transition of CO you would like to fit. Currently (v0.1.3), this package only includes functionality for 12CO, 13CO, and C18O in LTE. Please create an issue on the github repository if you would like a new species added.To help get you started, here are some reasonable fields for the config.yaml file for AK Sco:General synthesis parameters:name: AKSco\ngas: true\nspecies : 12CO # Possible choices: 12CO, 13CO, C18O. Future may include HCN, etc...\ntransition: 2-1 # J =The model grid setup:grid:\n  nr: 128\n  ntheta: 40 # if mirrored about the equator, total of 80\n  nphi: 1\n  r_in: 0.1 # [AU] # Inner edge of model grid\n  r_out: 300. # [AU] # Outer edge of model gridThe distance parameters. For now, we will keep distance fixed:fix_params : [\"dpc\"]Even though we are keeping the distance fixed, we need to specify these prior parameters:dpc_prior:\n	mu: 142.\n	sig: 20.Choose what type of model will we be fitting:model : standard # choices of {standard, truncated, vertical, cavity, etc..}Now come parameters that can be used to synthesize and plot models. Due to a quirk of how YAML files are read, make sure that each of these parameter values is a float and not an int (i.e., 1.0 vs. 1).parameters:\n  M_star: 2.49 # [M_sun] stellar mass\n  r_c: 14.00 # [AU] characteristic radius\n  T_10: 91.85 # [K] temperature at 10 AU\n  q: 0.51 # temperature gradient exponent\n  gamma: 1.0 # surface density gradient\n  logSigma_c: 3.0 # log surface density at char. radius\n  ksi: 0.31 # [km/s] microturbulence\n  dpc: 142. # [pc] distance\n  incl: 109.4 # [degrees] inclination\n  PA: 141.1 # [degrees] position angle\n  vel: -26.1 # [km/s]\n  mu_RA: 0.053 # [arcsec] centroid location\n  mu_DEC: 0.045 # [arcsec]Of all of these parameters, it might be hardest to guess correctly at the systemic velocity of the source. Generally, this is best done in the data reduction stages, for example taking a quick look at the central frequency of the spectral line. Note that the AK Sco dataset is provided in the raw topocentric frame, so the velocity quoted here is not the same as the LSRK quoted in the paper.Now, we need to specify how big we want our image to be and how many pixels it should have.# Image setup\nsize_arcsec : 12.0 # [arcsec] full width/height of image\nnpix: 512The final section is parameters that roughly describe the RMS in the observation and the approximate beam size. These parameters are only used in the channel map plotting script to help make the channel maps look more comparable to the observation, and you don't need to worry about them now.beam :\n  rms : 0.01 # Jy/beam\n  BMAJ: 0.9 # arcsec # Major axis of ellipse\n  BMIN: 0.9 # arcsec # Minor axis of ellipse\n  BPA: 1.0 # degrees east of North of the semi-major axis.Now, a good thing to check is that our setup parameters actually satisfy the Nyquist theorem. There is a helper script for this$ max_baseline.jl\nDataset channels are velocities from -10.517295788189767 to -42.656605289812504 and span -32.13930950162273 km/s.\nMidpoint is -26.586950539001137 km/s.\nMax baseline 338.5121039138373 kilolambda\nNyquist sampling satisfied. dRA: 0.0234375 [arcsec/pix] ; dRA_max: 0.2769671424693894 [arcsec/pix]\nImage size satisfied. Image size at the closest distances: 510.0 [AU]; outer radius of the grid + 10%: 330.0 [AU]It looks like everything is OK to start!  (If you see very different velocities here for the channels, check that you have correctly specified species and transition in config.yaml to match the spectral line actually observed in your dataset.) "
},

{
    "location": "cookbook.html#Makefile-1",
    "page": "Cookbook",
    "title": "Makefile",
    "category": "section",
    "text": "New in v0.1.3, I've written a Makefile which should simplify a lot of the necessary tasks within the directory for a single object. You can generally do everything you need to via make <target>, where the various targets will now be described."
},

{
    "location": "cookbook.html#Plotting-up-the-model-structure-1",
    "page": "Cookbook",
    "title": "Plotting up the model structure",
    "category": "section",
    "text": "To get a first pass glimpse at what the model will look like, you can make plots of the key quantities as a function of disk position.$ make structurewill create plots in your current working directory of velocity, temperature, and surface density. If want to play around with this, change a parameter in config.yaml and then rerun make structure(Image: Temperature)(Image: Surface Density)(Image: Velocity)"
},

{
    "location": "cookbook.html#Synthesizing-a-model-image-1",
    "page": "Cookbook",
    "title": "Synthesizing a model image",
    "category": "section",
    "text": "Before jumping into running any MCMC chains, it is a good idea to take a guess at some parameters and try synthesizing a set of model channel maps to see if your model looks remotely close to the dataset. Then, you can play around with values in config.yaml to see what works well.To jump right in, just try$ make chmapsAnd the code will start synthesizing channel maps. Because the AK Sco dataset contains a lot of channels, this may take about 5 minutes to get everything done. During this process, the output from RADMC-3D is piped to STDOUT. This may be a good place to debug if anything looks fishy.When complete, this should leave you with several chmaps_*.png files in your current directory. Take a look and see if these appear reasonable.For me, a typical workflow for playing around with channel maps isedit config.yaml to parameters that might make sense\nrun make structure to see that the disk properties look reasonable\nrun make chmaps to actually synthesize images\ninspect the resulting plots of the data (chmaps_linear.png), and if I am not satisfied go back to 1It is a very good idea to inspect your channel maps to make sure that there isn't any weird structure, that you have enough pixels to resolve the disk structure, and that your model grid appears to be at high enough resolution. A few extra minutes or hours spent debugging your images during this step can save you days (of supercomputer time) in the steps ahead.(Image: High Resolution Image of AK Sco)If you'd like to make a spatially-integrated spectrum, you can also do$ make spectrum.png(Image: Spectrum of AK Sco)"
},

{
    "location": "cookbook.html#Setting-up-a-parallelized-MCMC-exploration-of-the-parameters-1",
    "page": "Cookbook",
    "title": "Setting up a parallelized MCMC exploration of the parameters",
    "category": "section",
    "text": "As you just experienced, model synthesis can take a very long time, generally 1 - 5 minutes per model in the case of AK Sco. In order to explore the posterior in a reasonable amount of time, we need to parallelize the synthesis and evaluation of the likelihood function across multiple compute cores. This is done using a Julia port of the Ensemble Sampler by Goodman and Weare 2010, implemented in Python by Foreman-Mackey et al. as emcee. For more information about this great sampler, see here.Much like emcee, starting out requires deciding upon the positions of the walkers. To aid in placing these, the DJ_initialize.jl script copied over a Jupyter/Python notebook to your current directory. Now, open up InitializeWalkers.ipynb with a Jupyter notebook. We will change these following values to correspond to your disk of choice.To save you some computational time otherwise spent on burn-in, I found that the following walker starting positions worked well for AK Scop0 = np.array([np.random.uniform(2.4, 2.5, nwalkers), # mass [M_sun]\n          np.random.uniform(14., 15.0, nwalkers), #r_c [AU]\n          np.random.uniform(92., 93., nwalkers), #T_10 [K]\n          np.random.uniform(0.51, 0.55, nwalkers), # q\n          np.random.uniform(-3.5, -3.4, nwalkers), #log10 M_gas [log10 M_sun]\n          np.random.uniform(0.3, 0.32, nwalkers), #xi [km/s]\n          np.random.uniform(140.0, 144.0, nwalkers), #dpc [pc]\n          np.random.uniform(110.0, 112.0, nwalkers), #inc [degrees]\n          np.random.uniform(140.0, 141.0, nwalkers), #PA [degrees]\n          np.random.uniform(-26.1, -26.0, nwalkers), #vz [km/s]\n          np.random.uniform(0.0, 0.05, nwalkers), #mu_a [arcsec]\n          np.random.uniform(0.0, 0.4, nwalkers)]) #mu_d [arcsec]Then finish evaluating the rest of the cells so that you save the file pos0.npy into your current working directory.How many walkers should I use? Due the the way the Ensemble Sampler advances, you can only evaluate half of the walkers simultaneously. That means that if you are running with more cores than nwalkers/2, you will have several cores idle throughout the sampling. Of course, you could now increase the number of walkers to be 2 * ncores."
},

{
    "location": "cookbook.html#Launching-the-run-1",
    "page": "Cookbook",
    "title": "Launching the run",
    "category": "section",
    "text": "The exploration of the posterior is done via the scripts/venus.jl script. It is worth exploring this piece of code to see the various moving parts. How you invoke this script depends on your cluster environment. For a dataset the size of AK Sco, it's not worth your time to start run this script (except for debugging purposes) unless you have access to 20 or more cores."
},

{
    "location": "cookbook.html#Local-machine-with-20-cores-1",
    "page": "Cookbook",
    "title": "Local machine with 20 cores",
    "category": "section",
    "text": "If you have your own 20-core machine, you can launch the script via  $ venus.jl -p 19Much like the Julia interpreter itself, the -p argument will add an addition 19 workers to the master process for a total of 20 workers."
},

{
    "location": "cookbook.html#High-Performance-Cluster-1",
    "page": "Cookbook",
    "title": "High Performance Cluster",
    "category": "section",
    "text": "These Julia scripts can take advantage many possible cores spread across multiple nodes. This may require some custom script writing for your specific cluster situation, but the main ideas are as followswrite a submission script that specifies total number of cores, time, memory, etc\nupon submission, determine how many cores you have been allocated on which nodes\ncreate a hosts.txt file which contains this information, following the Julia spec here.\nthen start your job with\njulia –machinefile hosts.txt venus.jl"
},

{
    "location": "cookbook.html#Examining-the-output-1",
    "page": "Cookbook",
    "title": "Examining the output",
    "category": "section",
    "text": "Because the MCMC run is so expensive, the code is designed to periodically write out snapshots of the samples, to both save your progress and allow you to check up on the chains mid-run. You can set the cadence in the config.yaml file under# MCMC setup\nsamples: 10\nloops: 1Basically, in each loop, the code advances N samples. After each loop, the full chain (chain.npy) is written to the output/ directory along with a current snapshot of the walker positions, pos0.npy. This way if your venus.jl script ends (either by design or cluster failure), you can copy this pos0.npy back to your project directory and start from the last known walker positions.At this point you can take the samples in chain.npy and analyze them as you would normal MCMC samples. To save you the trouble, however, we included a script to help with these tasks$ plot_walkers.pywill plot the walker positions as a function of iteration (walkers.png). Examining walkers.png is a decent way to estimate if your chains are done with burn-in.When you are ready, you can burn off these first say 300 (or more) iterations and make a corner plot (triangle.png) by$ plot_walkers.py --burn 300 --tri"
},

{
    "location": "priors.html#",
    "page": "Priors",
    "title": "Priors",
    "category": "page",
    "text": ""
},

{
    "location": "priors.html#Priors-1",
    "page": "Priors",
    "title": "Priors",
    "category": "section",
    "text": "When fitting certain disks, it may be worthwhile to include new information from separate analyses. For example, when fitting multiple CO transitions, it may be worthwhile to make priors that constrain the temperature profile.Because flexible priors like these are difficult to \"hard-code\" into the package, there is an additional functionality that allows the user to write their own prior, which will overwrite the default prior in the package.You can copy a sample stub to your current working directory via$ DJInitialize.jl --priorThen, open this file with your favorite text editor. It is important that you mimic the exact same function call as in src/model.jl."
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#DiskJockey.constants.fftspace-Tuple{Real,Int64}",
    "page": "API",
    "title": "DiskJockey.constants.fftspace",
    "category": "Method",
    "text": "Oftentimes it is necessary to get a symmetric coordinate array that spans N  elements from -width to +width, but makes sure that the middle point lands  on 0. The indices go from 0 to N -1  linspace returns  the end points inclusive, wheras we want to leave out the  right endpoint, because we are sampling the function in a cyclic manner.\n\n\n\n"
},

{
    "location": "api.html#Functions-1",
    "page": "API",
    "title": "Functions",
    "category": "section",
    "text": "constants.fftspace(width::Real, N::Int)"
},

{
    "location": "api.html#Index-1",
    "page": "API",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "changelog.html#",
    "page": "Changelog",
    "title": "Changelog",
    "category": "page",
    "text": ""
},

{
    "location": "changelog.html#Changelog-1",
    "page": "Changelog",
    "title": "Changelog",
    "category": "section",
    "text": "The following are the changes that have been implemented since the previous version. All versions are tagged as releases and available as a clone off of the master branch."
},

{
    "location": "changelog.html#Version-0.1.4-1",
    "page": "Changelog",
    "title": "Version 0.1.4",
    "category": "section",
    "text": "Upgraded to Julia version 0.6. To upgrade to this version of DiskJockey, you will first need to upgrade your own Julia distribution to 0.6 as well. I apologize for the frequent upgrades of Julia, but since Julia is itself a fast-developing language, it makes the most sense to just bite the bullet and upgrade every time the programming language upgrades as well. Once Julia reaches v1.0 (sometime in 2018, perhaps), these changes should become less frequent.After you install the new version of Julia, make sure that it actually loads as$ julia\n               _\n   _       _ _(_)_     |  A fresh approach to technical computing\n  (_)     | (_) (_)    |  Documentation: https://docs.julialang.org\n   _ _   _| |_  __ _   |  Type \"?help\" for help.\n  | | | | | | |/ _` |  |\n  | | |_| | | | (_| |  |  Version 0.6.1-pre.92 (2017-10-07 01:18 UTC)\n _/ |\\__'_|_|_|\\__'_|  |  Commit 389b23cf6e* (6 days old release-0.6)\n|__/                   |  x86_64-pc-linux-gnu\n\njulia>And then you will need to reinstall DiskJockey in this new version, which should automatically pull down all of the other packages.julia> Pkg.clone(\"https://github.com/iancze/DiskJockey.git\")If you no longer plan on using Julia v0.4 or v0.5 or any of the packages, I would recommend deleting it from your system so as not to cause any PATH conflicts. Speaking of which, don't forget to update your PATH to point to the new v0.6 version scripts."
},

{
    "location": "changelog.html#Updated-RADMC-3D-1",
    "page": "Changelog",
    "title": "Updated RADMC-3D",
    "category": "section",
    "text": "We're now running the latest version of RADMC-3D (0.41), which is automatically downloaded and installed with DiskJockey. Although not strictly required for running DiskJockey, it might help to familiarize yourself with the RADMC-3D manual, in particular the sections on input files and LTE line transfer."
},

{
    "location": "changelog.html#Renaming-1",
    "page": "Changelog",
    "title": "Renaming",
    "category": "section",
    "text": "The script max_baseline.jl is now DJ_max_baseline.jl. plot_walkers.py is now DJ_plot_walkers.py."
},

{
    "location": "changelog.html#DJ_plot_walkers.py-1",
    "page": "Changelog",
    "title": "DJ_plot_walkers.py",
    "category": "section",
    "text": "Now has proper labeling for models that fix parameters."
},

{
    "location": "changelog.html#DJ_verify_run.jl-1",
    "page": "Changelog",
    "title": "DJ_verify_run.jl",
    "category": "section",
    "text": "A new script for checking that you've specified everything properly before launching a MCMC run. This way you can avoid syntax errors in the config scripts after queuing for a cluster job."
},

{
    "location": "changelog.html#Version-0.1.3-1",
    "page": "Changelog",
    "title": "Version 0.1.3",
    "category": "section",
    "text": ""
},

{
    "location": "changelog.html#Fix-parameters-1",
    "page": "Changelog",
    "title": "Fix parameters",
    "category": "section",
    "text": "Now there is a flexible way to fix or free parameters in a fit.To do this requires new config files (which you may need to copy over from assets/) that include a params_fixed field, where you list the parameters to fix. Then, it may also require modifying the InitializeWalkers.ipynb notebooks to incorporate more or fewer parameters than are currently available."
},

{
    "location": "changelog.html#Cavity-model-1",
    "page": "Changelog",
    "title": "Cavity model",
    "category": "section",
    "text": "Now supports a cavity model, with parameters for the size of the cavity (r_cav), as well as the steepness (gamma_cav)."
},

{
    "location": "changelog.html#Truncated-model-1",
    "page": "Changelog",
    "title": "Truncated model",
    "category": "section",
    "text": "Now includes a model that has a variable outer exponential taper, allowing it to be steeper or shallower than the traditional (2 - gamma)."
},

{
    "location": "changelog.html#Change-to-Sigma_c-1",
    "page": "Changelog",
    "title": "Change to Sigma_c",
    "category": "section",
    "text": "Instead of fitting with the parameter total gas mass, we now fit with the surface density normalization at the critical radius, Sigma_c. Note that this only truly the normalization constant for the standard, truncated, and vertical models. The reason for the switch from M_gas to Sigma_c is that for more complicated models, there is no analytic formula for the total gas mass, meaning that a numerical integral would be needed for each model evaluation. It is simpler and more accurate to sample in Sigma_c and then later convert the samples to M_gas if desired."
},

{
    "location": "changelog.html#FITS-export-1",
    "page": "Changelog",
    "title": "FITS export",
    "category": "section",
    "text": "Thanks to Jane Huang (@j6626), you can now export a RADMC image to a FITS file via the script DJ_image_to_FITS.py. You can then inspect these FITS files with ds9, or use them in CASA simobserve."
},

{
    "location": "changelog.html#spectrum.png-1",
    "page": "Changelog",
    "title": "spectrum.png",
    "category": "section",
    "text": "The plotting routine now outputs integrated line flux, to be used as a sanity check against the measured value from a dataset."
},

{
    "location": "changelog.html#Half-pixel-offset-1",
    "page": "Changelog",
    "title": "Half-pixel offset",
    "category": "section",
    "text": "Has now been added to the code. This is a result from that RADMC synthesizes an image centered on (0,0), while the FFT routine expects the image to be centered on the middle of the central pixel. The only change you should see is a small, half-pixel sized offset in mu_RA and mu_DEC."
},

{
    "location": "changelog.html#Visualizing-the-tau1-surface-1",
    "page": "Changelog",
    "title": "Visualizing the tau=1 surface",
    "category": "section",
    "text": "This functionality originally provided in RADMC is now exposed via DiskJockey to visualize the tau = 1 (or other value) surface.$ DJ_initialize.jl && tausurf_model.jl && plot_tausurf.jlThis also introduces the TausurfImage type in src/image.jl. This is useful to see whether the tau=1 surface is in front of or behind the midplane of the disk, projected on the sky.It will be more useful to plot this in 3D, however. ParaView export is hopefully coming in a future version."
},

{
    "location": "changelog.html#Vertical-temperature-gradient-1",
    "page": "Changelog",
    "title": "Vertical temperature gradient",
    "category": "section",
    "text": "We now include support for fitting a vertical temperature gradient, via the vertical model type. We follow the parameterization in Williams and Best 14."
},

{
    "location": "changelog.html#model.convert_vector-1",
    "page": "Changelog",
    "title": "model.convert_vector",
    "category": "section",
    "text": "This routine is to simplify to ingesting kwargs, enabling easier selection of model types and parameters."
},

{
    "location": "changelog.html#User-defined-priors-1",
    "page": "Changelog",
    "title": "User-defined priors",
    "category": "section",
    "text": "As an experimental feature, it is now possible to assign user defined priors. This will require you to write some Julia code, but shouldn't be that difficult. An example is in assets/prior.jl, which will be read automatically from the current working directory. The main idea is to make it easier to assign disk-specific priors (i.e., constrain the disk gas mass to a specific value, or disk position, etc...)."
},

{
    "location": "changelog.html#Makefiles-1",
    "page": "Changelog",
    "title": "Makefiles",
    "category": "section",
    "text": "Now, we include a makefile that should be copied to each source directory. In a typical analysis workflow, there isn't need to regenerate each intermediate product repeatedly. Details are in the cookbook."
},

{
    "location": "changelog.html#Script-renaming-1",
    "page": "Changelog",
    "title": "Script renaming",
    "category": "section",
    "text": "To prevent these scripts from cluttering a users namespace when adding to their PATH, we've prefixed most things with DJ_. Moreover, for most tasks, the users will not need to run these tasks directly, but will just use the Makefile."
},

{
    "location": "changelog.html#model.Grid-1",
    "page": "Changelog",
    "title": "model.Grid",
    "category": "section",
    "text": "Grid type initialization now relies upon a dictionary."
},

{
    "location": "changelog.html#Version-0.1.2-1",
    "page": "Changelog",
    "title": "Version 0.1.2",
    "category": "section",
    "text": "The package is now named DiskJockey (previously JudithExcalibur). Check out the new logo!"
},

{
    "location": "changelog.html#Installation-1",
    "page": "Changelog",
    "title": "Installation",
    "category": "section",
    "text": "The README now provides installation instructions for tagged releases. Rather than requiring the user to install it separately, RADMC-3D is now installed automatically as part of the installation process."
},

{
    "location": "changelog.html#UVHDF5-1",
    "page": "Changelog",
    "title": "UVHDF5",
    "category": "section",
    "text": "Conversion from UVFITS and CASA measurement set was previously handled by scripts in this repository. However, this import and export capability has be standardized via the UVHDF5 package. Please see that package for any issues with import and export."
},

{
    "location": "changelog.html#Travis-integration-1",
    "page": "Changelog",
    "title": "Travis integration",
    "category": "section",
    "text": "Now the builds are tested on travis-ci, which should hopefully increase stability for the development cycle."
},

{
    "location": "changelog.html#EnsembleSampler-1",
    "page": "Changelog",
    "title": "EnsembleSampler",
    "category": "section",
    "text": "There is now an expand_walkers.py script, designed to add the dpc dimension to the sampling routine from an ensemble of walkers used on a posterior with distance fixed.Now the ensemble sampler takes an optional function, to be called at the end of each loop as func(sampler, outdir)function run_schedule(sampler::Sampler, pos0, N::Int, loops::Int, outdir, func::Function=nothing)See the src/EnsembleSampler.jl file for more details. This is primarily in support of the next feature..."
},

{
    "location": "changelog.html#Plotly-script-1",
    "page": "Changelog",
    "title": "Plotly script",
    "category": "section",
    "text": "If you run the sampling script as$ venus.jl --plotlyThen after the completion of each sampling loop, this will create a plotly walkers plot corresponding to the name entry in config.yaml. This allows easy monitoring of many different chains that might be running on a cluster."
},

{
    "location": "changelog.html#DJInitialize.jl-1",
    "page": "Changelog",
    "title": "DJInitialize.jl",
    "category": "section",
    "text": "(Previously JudithInitialize.jl). Now this initialization script can easily spit out an exclude array so that fewer channels can be fit during initial testing."
},

{
    "location": "changelog.html#Cavity-model-2",
    "page": "Changelog",
    "title": "Cavity model",
    "category": "section",
    "text": "By initializing a directory with$ DJInitialize.jl --new-project=cavityyou can start exploring the cavity model, which has an exponential taper inside of some radius, r_cav."
},

{
    "location": "changelog.html#gridding.jl-1",
    "page": "Changelog",
    "title": "gridding.jl",
    "category": "section",
    "text": "gridding.jl now exports a corrfun(img::SkyImage) routine, in addition to corrfun!(img::SkyImage). This new function returns a corrected image as a copy, leaving the original image in the arguments unchanged. This is useful for plotting and debugging scripts so that you don't need to copy the image manually."
},

{
    "location": "changelog.html#visibilities.jl-1",
    "page": "Changelog",
    "title": "visibilities.jl",
    "category": "section",
    "text": "FullModelVis now has a - method, which can be used to subtract two sets of dense visibilities."
},

{
    "location": "changelog.html#image.jl-1",
    "page": "Changelog",
    "title": "image.jl",
    "category": "section",
    "text": "Image now has a - method, which can be used to subtract two images."
},

{
    "location": "changelog.html#config.yaml-copied-to-output-directory-1",
    "page": "Changelog",
    "title": "config.yaml copied to output directory",
    "category": "section",
    "text": "Now venus.jl will move a copy of your config.yaml file to the output directory. This creates an automatic record of what parameters you ran with, which will undoubtedly be useful when reviewing previous runs sometime in the near future. Default values in the initial config.yaml files have also been updated."
},

{
    "location": "changelog.html#plot_baselines.jl-1",
    "page": "Changelog",
    "title": "plot_baselines.jl",
    "category": "section",
    "text": "Will generate a plot of where your dataset has sampled the UV plane. Also will print out an average velocity of the dataset, to allow a good starting guess for the vel parameter."
},

{
    "location": "changelog.html#Version-0.1.1-1",
    "page": "Changelog",
    "title": "Version 0.1.1",
    "category": "section",
    "text": ""
},

{
    "location": "changelog.html#Model-specification-1",
    "page": "Changelog",
    "title": "Model specification",
    "category": "section",
    "text": "New implementation of models through parameter types in model.jl. Instead of separate model files for each new parameterization, which duplicated a lot of code and made things difficult to maintain, we are now collating all model types in model.jl by their parameters. We have created the AbstractParameters umbrella type, with subtypes ParametersStandard, ParametersTruncated, and ParametersCavity. Previous model specification routines like Sigma(r::Float64, pars::Parameters) (surface density) are now have overloaded methods for different models based upon these parameter typesSigma(r::Float64, pars::ParametersStandard)\nSigma(r::Float64, pars::ParametersTruncated)\nSigma(r::Float64, pars::ParametersCavity)And routines that are general to all models (e.g. velocity specification) are denoted byvelocity{T}(r::T, pars::AbstractParameters)Note that model types cavity and truncated are still experimental and likely to change.To go along with this change, I have updated the automatic generation of the config.yaml file to include new fields like model: standard. Within the config.yaml file, items in the parameters dictionary are now simply a single Float64 number, not an array of [starting, jump] like it was previously.model.jl now handles the implementation of priors, that dispatch off of the parameter types. This greatly streamlines the function fprob in venus.jl."
},

{
    "location": "changelog.html#One-off-model-synthesis-and-plotting-1",
    "page": "Changelog",
    "title": "One-off model synthesis and plotting",
    "category": "section",
    "text": "I've simplified synthesizing and plotting of models. What was previously plot_model.jl is now split into synthesize_model.jl and plot_chmaps.jl.Added new plot_moments.jl script which plots the zeroth-moment image. Soon it will plot first moment as well."
},

{
    "location": "changelog.html#MCMC-Sampling-1",
    "page": "Changelog",
    "title": "MCMC Sampling",
    "category": "section",
    "text": "The specification of parameter types allowed me to greatly simplify the MCMC sampling code into a single script, venus.jl. I have moved previous sampling scripts to the attic/ directory.Created InitializeWalkers.ipynb that is copied to new directory to help specify walker starting positions. The user edits this with a Jupyter/IPython notebook."
},

{
    "location": "changelog.html#Cookbook-1",
    "page": "Changelog",
    "title": "Cookbook",
    "category": "section",
    "text": "We now have a cookbook for AK Sco, check it out to get started!plot_walkers.py now includes ability to determine highest density interval for quoting credible intervals and should automatically label the parameters after reading from config.yaml."
},

{
    "location": "changelog.html#Version-0.1.0-1",
    "page": "Changelog",
    "title": "Version 0.1.0",
    "category": "section",
    "text": "Initial commit on the new versioning roadmap."
},

]}
