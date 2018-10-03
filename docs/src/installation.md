# Installation

## Julia

First, you should to install the [Julia programming language](https://github.com/Astrochem/UVHDF5) on your machine. Instructions can be found [here](http://julialang.org/downloads/). I have found it exciting to install from source via the [github repo](https://github.com/JuliaLang/julia/), but it may be a bit quicker to use a pre-compiled binary. If you've successfully installed everything, you should be able to open up a Julia interpreter and type

    $ julia
                   _
       _       _ _(_)_     |  A fresh approach to technical computing
      (_)     | (_) (_)    |  Documentation: http://docs.julialang.org
       _ _   _| |_  __ _   |  Type "?help" for help.
      | | | | | | |/ _` |  |
      | | |_| | | | (_| |  |  Version 0.6.2-pre+12
     _/ |\__'_|_|_|\__'_|  |  Commit 9e598b6* (83 days old release-0.4)
    |__/                   |  x86_64-unknown-linux-gnu

    julia> println("Hello world")
    Hello world

Depending on how you choose to install Julia, you may need to take the additional step of adding the Julia executable to your `PATH`. [Here](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started#Running_directly_from_terminal) are some suggestions for OS X if you are experiencing difficulty.

## Fortran

Some of the packages that DiskJockey requires (e.g. [Dierckx.jl](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started#Running_directly_from_terminal)) require a Fortran compiler. On linux, try looking to make sure you have `gcc/gfortran` or something like it installed. For OS X, you can download these packages from [here](http://hpc.sourceforge.net/).

## DiskJockey

Next, we will install the `DiskJockey` package itself. Because this is not yet an official Julia package, for now, installation involves simply cloning the repository. First, open up a Julia prompt in the REPL, then type

    julia> Pkg.clone("https://github.com/iancze/DiskJockey.git")
    julia> Pkg.build("DiskJockey")

This process may take a few minutes as the relevant packages (including RADMC-3D) are downloaded from the web and installed. So far, I have only been able to extensively test this installation process on Linux machines. If you run into errors in this build process, please file an [issue](https://github.com/iancze/DiskJockey/issues) on the github repository so that we may try to fix this. If you already have RADMC-3D installed on your system, this process won't interfere with that installation, `DiskJockey` will use the version of RADMC-3D it downloaded.

As mentioned previously, there are several "driver" command line scripts that are used to perform the actual mass fitting. To complete the installation, you should add these files to your system PATH. To figure out where the package is installed

    julia> Pkg.dir("DiskJockey")
    "/home/ian/.julia/v0.6/DiskJockey"

Your PATH will vary. The scripts are located inside of the `scripts` directory, so if you are using bash, you will want to add the PATH that looks something like

    export PATH="/home/ian/.julia/v0.6/DiskJockey/scripts:$PATH"

inside of your `.bashrc` file.

Some of the Python scripts also depend on imports from this directory, so you will also need to add

    export PYTHONPATH="/home/ian/.julia/v0.6/DiskJockey/scripts:$PYTHONPATH"

inside of your `.bashrc` file. Finally,

    $ source ~/.bashrc

To check that you have properly added the scripts (and check that you installed the correct version), you can try in your system shell

    $ DJInitialize.jl --version
    Your DiskJockey scripts are successfully linked.
    You are running DiskJockey 0.1.4
    Exiting

Due to the ongoing development of this package, it is easiest to keep your version current by updating off of the git master branch. For example,

    julia> Pkg.dir("DiskJockey")
    "/home/ian/.julia/v0.6/DiskJockey

    # Exit Julia
    # Navigate to this directory
    $ cd /home/ian/.julia/v0.6/DiskJockey
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
