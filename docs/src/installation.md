# Installation

## Julia

First, you should to install the [Julia programming language](http://julialang.org/downloads/) on your machine. 

Depending on how you choose to install Julia, you may need to take the additional step of adding the Julia executable to your `PATH`. [Here](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started#Running_directly_from_terminal) are some suggestions for OS X if you are experiencing difficulty.

Julia uses "environments" to manage package dependencies. If you are unfamiliar, it is worth reading through the documentation on how Julia's package manager works with environments [here](https://docs.julialang.org/en/v1/stdlib/Pkg/). 

## DiskJockey

`DiskJockey` is not (yet) registered an official Julia package, so it needs to be installed by using the link to the github repository. First, open up a Julia prompt in the REPL, then type `]` to enter the package manager, then use the `add` command to clone the DiskJockey repository.

    julia> ] 
    (@v1.4) pkg> add https://github.com/iancze/DiskJockey.git
    
This process may take a few minutes as the relevant packages (including RADMC-3D) are downloaded from the web and installed. If you run into errors in this build process, please file an [issue](https://github.com/iancze/DiskJockey/issues) on the github repository so that we may try to fix this. If you already have RADMC-3D installed on your system, this process won't interfere with that installation, `DiskJockey` will use the version of RADMC-3D it downloaded during the build process. If you anticipate editing the source code, rather than installing the package with the `add` command, you can install using the `develop` command

    julia> ] 
    (@v1.4) pkg> develop https://github.com/iancze/DiskJockey.git

Included within the DiskJockey github repo are also several "driver" command line scripts that are used to perform the actual mass fitting. To utilize these files, you should add the scripts directory to your system PATH. To figure out where the package is installed

    julia> using DiskJockey
    julia> dirname(pathof(DiskJockey))
    "/home/ian/.julia/dev/DiskJockey/src"

Your PATH will vary. The scripts are located inside of the `DiskJockey/scripts` directory, so if you are using bash, you will want to add the PATH that looks something like

    export PATH="/home/ian/.julia/dev/DiskJockey/scripts:$PATH"

inside of your `.bashrc` file.

Some of the Python scripts also depend on imports from this directory, so you will also need to add

    export PYTHONPATH="/home/ian/.julia/dev/DiskJockey/scripts:$PYTHONPATH"

inside of your `.bashrc` file. 

For archival purposes, tagged release versions of this package are available [here](https://github.com/iancze/DiskJockey/releases).

## Python scripts

Some of the analysis scripts and IO routines require Python 3.x and several Python packages. Please install the following packages via your own package manager or a distribution like anaconda:

* numpy
* scipy
* matplotlib
* h5py
* PyYAML
* Jupyter/IPython
* [corner.py](https://github.com/dfm/corner.py)

With the package successfully installed, see the documentation in the `docs/` folder on how to get started fitting a specific source, in particular the [Cookbook](@ref).

If you'd like to run the test suite to make sure everything checks out, start Julia and run

    julia> ]
    (@v1.4) pkg> activate .
    (@v1.4) pkg> test

this may take about 10 minutes or so. If you catch any errors for your specific machine, please report them via an Issue on the github repository.

The code package is designed to interface with visibilities in the UVHDF5 format, described [here](https://github.com/Astrochem/UVHDF5). UVHDF5 also provides scripts to convert to and from UVFITS and CASA measurement sets.
