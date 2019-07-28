# Installation

## Julia

First, you should to install the [Julia programming language](https://github.com/Astrochem/UVHDF5) on your machine. Instructions can be found [here](http://julialang.org/downloads/). I have found it exciting to install from source via the [github repo](https://github.com/JuliaLang/julia/), but it is usually quicker to use a pre-compiled binary. If you've successfully installed everything, you should be able to open up a Julia interpreter and type

    $ julia
                _
    _       _ _(_)_     |  Documentation: https://docs.julialang.org
    (_)     | (_) (_)    |
    _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
    | | | | | | |/ _` |  |
    | | |_| | | | (_| |  |  Version 1.1.1 (2019-05-16)
    _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
    |__/                   |

    julia> 

Depending on how you choose to install Julia, you may need to take the additional step of adding the Julia executable to your `PATH`. [Here](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started#Running_directly_from_terminal) are some suggestions for OS X if you are experiencing difficulty.

## Fortran

Some of the packages that DiskJockey requires (e.g. [Dierckx.jl](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started#Running_directly_from_terminal)) require a Fortran compiler. On linux, try looking to make sure you have `gcc/gfortran` or something like it installed. For OS X, you can download these packages from [here](http://hpc.sourceforge.net/).

## DiskJockey

Next, we will install the `DiskJockey` package itself. Because this is not an official Julia package (for now) installation involves simply cloning the repository. Full details on Julia's (new) package manager are available [online](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html). Here we summarize the basics. First, open up a Julia prompt in the REPL, then type `]` to enter the package manager

    julia> ] 
    (v1.1) pkg> dev https://github.com/iancze/DiskJockey.git
    
This process may take a few minutes as the relevant packages (including RADMC-3D) are downloaded from the web and installed. So far, I have only been able to extensively test this installation process on Linux machines. If you run into errors in this build process, please file an [issue](https://github.com/iancze/DiskJockey/issues) on the github repository so that we may try to fix this. If you already have RADMC-3D installed on your system, this process won't interfere with that installation, `DiskJockey` will use the version of RADMC-3D it downloaded. Note that we are using the `dev` command to clone the package, since you may eventually want to edit the source code.

If your internet interrupts or something else goes wrong during this process, you may need to restart the build using 

    (v1.1) pkg > build DiskJockey

If this process completes, congratulations you've succesfully installed the DiskJockey package, and you should be able to do 

    julia> using DiskJockey 

With the release of Julia 1.x, package dependencies are now handled in a slightly more sophisticated manner. While the installation of the DiskJockey package did install all of its dependencies, those installed dependencies are only available to the package itself. Included within the DiskJockey github repo are also several "driver" command line scripts that are used to perform the actual mass fitting. To utilize these files, you should add the scripts directory to your system PATH. To figure out where the package is installed

    julia> dirname(pathof(DiskJockey))
    "/home/ian/.julia/dev/DiskJockey/src"

Your PATH will vary. The scripts are located inside of the `DiskJockey/scripts` directory, so if you are using bash, you will want to add the PATH that looks something like

    export PATH="/home/ian/.julia/dev/DiskJockey/scripts:$PATH"

inside of your `.bashrc` file.

Some of the Python scripts also depend on imports from this directory, so you will also need to add

    export PYTHONPATH="/home/ian/.julia/v0.6/DiskJockey/scripts:$PYTHONPATH"

inside of your `.bashrc` file. Finally,

    $ source ~/.bashrc

Lastly, many of these scripts have dependencies on Julia packages other than DiskJockey. To install these packages to the general environment, you will need to do 

    julia > ]
    (v1.1) pkg > add AbstractFFTs ArgParse ClusterManagers HDF5 LaTeXStrings NPZ PyPlot YAML Test

As much as I understand it, this is mostly making these packages available to the current "general" environment from where we'll run the scripts, rather than completely reinstalling them. To check that you have properly added the scripts (and check that you installed the correct version), you can try in your system shell

    $ DJInitialize.jl --version
    Your DiskJockey scripts are successfully linked.
    You are running DiskJockey 0.1.4
    Exiting

Due to the ongoing development of this package, it is easiest to keep your version current by updating off of the git master branch. For example,

    # Navigate to this directory
    $ cd /home/ian/.julia/dev/DiskJockey
    $ git pull

For archival purposes, tagged release versions of this package are available [here](https://github.com/iancze/DiskJockey/releases).

## Python scripts

Lastly, some of the analysis scripts and IO routines require Python and several Python packages. I have only tested the scripts on Python 3.x, so they are unlikely to work on Python 2.7. Please install the following packages via your own package manager or a distribution like anaconda:

* numpy
* scipy
* matplotlib
* h5py
* PyYAML
* Jupyter/IPython
* [corner.py](https://github.com/dfm/corner.py)

With the package successfully installed, see the documentation in the `docs/` folder on how to get started fitting a specific source, in particular the [Cookbook](@ref).

If you'd like to run the test suite to make sure everything checks out, start Julia and run

    julia> Pkg.test("DiskJockey")

this may take about 10 minutes or so. If you catch any errors for your specific machine, please report them via an Issue on the github repository.

The code package is designed to interface with visibilities in the UVHDF5 format, described [here](https://github.com/Astrochem/UVHDF5). There are also scripts within this repository to convert to and from UVFITS and CASA measurement sets.
