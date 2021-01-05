# Installation

## Julia

First, you should to install the [Julia programming language](http://julialang.org/downloads/) on your machine. Depending on how you choose to install Julia, you may need to take the additional step of adding the Julia executable to your `PATH`.

Julia uses "environments" to manage package dependencies. These are similar in concept to Python's virtual environments, and very useful in practice. If you are unfamiliar, it is highly recommended to read through the Julia environment documentation [here](https://docs.julialang.org/en/v1/stdlib/Pkg/), since they are used to install DiskJockey in the recommended configuration. 

To get started, create a new folder that you will contain your work for a particular project.

    $ mkdir myproject 

Then, start Julia, and use the interactive package manager `Pkg` (enter via `]`) to create an environment specific to this project. 

    $ cd myproject
    $ julia 
    julia> ]
    pkg> activate . 
    Activating new environment at `~/myproject/Project.toml`
    (myproject) pkg>

## DiskJockey Package

Next, we will install the DiskJockey package to this same Julia environment. `DiskJockey` is not (yet) registered an official Julia package, so it needs to be installed by using the link to the github repository. 

    (myproject) pkg> add https://github.com/iancze/DiskJockey.git
    
This process may take a few minutes as the relevant packages (including RADMC-3D) are downloaded from the web and installed. If you run into errors in this build process, please file an [issue](https://github.com/iancze/DiskJockey/issues) on the github repository so that we may try to fix this. If you already have RADMC-3D installed on your system, this process won't interfere with that installation, `DiskJockey` will use the version of RADMC-3D it downloaded during the build process. 

It is a good idea to set the environment variable `JULIA_PROJECT` to the location for this project

    export JULIA_PROJECT='/home/ian/myproject'

You can add this line to your `.bashrc` or shell startup, or local directory initialization script. That way, whenever you start Julia, it will already activate the myproject environment. This makes running the scripts in the next section much easier. Advanced users with multiple active Julia projects may wish to explore alternate package activation patterns.

## DiskJockey Scripts

### Julia Scripts

Included within the DiskJockey github repo are several "driver" command line scripts in the `DiskJockey/scripts` directory that you can use to perform the actual mass fitting and analysis plotting. To utilize these files, you should add the `scripts` directory to your system PATH. To figure out where the package is installed

    julia> using DiskJockey
    julia> dirname(pathof(DiskJockey))
    "/home/ian/.julia/dev/DiskJockey/src"

Your PATH will vary. The scripts are located inside of the `DiskJockey/scripts` directory, so if you are using bash, you will want to add the PATH that looks something like

    export PATH="/home/ian/.julia/dev/DiskJockey/scripts:$PATH"

inside of your `.bashrc` file.

The scripts depend on other Julia packages for reading configuration files (in YAML) and plotting. You can either install these packages to your `myproject` environment just like you did DiskJockey, or you can use the convenience script 

    $ DJ_initialize_environment.jl 


### Python Scripts

Some of the Python scripts also depend on imports from this directory, so you will also need to add

    export PYTHONPATH="/home/ian/.julia/dev/DiskJockey/scripts:$PYTHONPATH"

inside of your `.bashrc` file. 


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


## Development workflow

If you anticipate editing the source code, rather than installing the package with the `add` command, you can install using the `develop` command

    julia> ] 
    (@v1.4) pkg> develop https://github.com/iancze/DiskJockey.git


For archival purposes, tagged release versions of this package are available [here](https://github.com/iancze/DiskJockey/releases).
