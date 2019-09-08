#!/usr/bin/env julia

# The main purpose of this file is to read the config file and copy all of the relevant static
# files into this directory, write out the proper wavelength files, and make sure we can continue
# with the run properly.

using Pkg; Pkg.activate("DiskJockey")

using ArgParse

s = ArgParseSettings(description="Initialize a new project directory with the appropriate files. Can also be used to update RADMC-3D input files after making changes to config.yaml")
@add_arg_table s begin
    "--version"
    help = "Print a the version number and exit."
    action = :store_true
    "-M"
    help = "Print out the list of files this program generates."
    action = :store_true
    "--new-project"
    help = "Copy a stock configuration file to this directory. Can be 'standard', 'truncated', 'cavity', 'nuker'"
    default = "no"
    "--fit-every"
    help = "Print out an exclude array that should be manually copied to config.yaml so that only every n-th channel is fit in venus.jl"
    arg_type = Int
    default = 0
    "--prior"
    help = "Copy a default prior.jl template to this current directory, for user modification. CAREFUL! See the docs about this functionality before invoking."
    action = :store_true
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
    "--vstart", "-s"
    help = "Starting velocity in km/s."
    arg_type = Float64
    "--vend", "-e"
    help = "Ending velocity in km/s."
    arg_type = Float64
    "--nvel", "-n"
    help = "Number of velocity channels."
    arg_type = Int
end

parsed_args = parse_args(ARGS, s)

# rather kludgey substitute for previous Pkg.dir() call
using DiskJockey
srcdir = dirname(pathof(DiskJockey))
assets_dir = joinpath(srcdir, "../assets/")

# The user is going to start modeling a new disk, so copy in the new configuration file.
if parsed_args["new-project"] != "no"
    model = parsed_args["new-project"]

    cp(assets_dir * "config.$(model).yaml", pwd() * "/config.yaml")
    cp(assets_dir * "initialize_walkers.$(model).py", pwd() * "/initialize_walkers.py")
    cp(assets_dir * "Makefile", pwd() * "/Makefile")

    println("Copied default config.yaml, initialize_walkers.py, and Makefile for the $model model to current working directory.")
    println("Exiting")
    exit()
end

if parsed_args["prior"]
    cp(assets_dir * "prior.jl", pwd() * "/prior.jl")
    println("prior.jl copied to current directory. If you didn't mean to do this, delete this file now. Otherwise, edit this file with your favorite text editor to enforce the prior you would like to see for this specific disk. Please see the docs for more information.")
    println("Exiting")
    exit()
end

using DiskJockey.constants

# This is just for new users to test that the package is successfully installed.
if parsed_args["version"]
    println("Your DiskJockey scripts are successfully linked.")
    println("You are running DiskJockey $DISKJOCKEY_VERSION")
    println("Exiting")
    exit()
end

import YAML
config = YAML.load(open(parsed_args["config"]))

if parsed_args["M"]

    str = "amr_grid.inp camera_wavelength_micron.inp gas_temperature.inp gas_velocity.inp lines.inp microturbulence.inp radmc3d.inp wavelength_micron.inp"

    species = config["species"]
    transition = config["transition"]

    if species == "12CO"
        str *= " molecule_co.inp numberdens_co.inp"
    elseif species ==  "13CO"
        str *= " molecule_13co.inp numberdens_13co.inp"
    elseif species == "C18O"
        str *= " molecule_c18o.inp numberdens_c18o.inp"
    end

    println(str)
    exit()
end

using DiskJockey.model
using DiskJockey.visibilities
using HDF5

fit_every = parsed_args["fit-every"]
if fit_every != 0
    dvarr = DataVis(config["data_file"])
    nchan = length(dvarr)
    println("Dataset contains a total of $nchan channels.")

    to_fit = Int[i for i=1:fit_every:nchan]

    # which channels of the dset to exclude
    exclude = filter(x->(!in(x, to_fit)), Int[i for i=1:nchan])

    println("To fit only every $(fit_every)-nd/rd/th channel, place the following line in your `config.yaml file.`")
    println("exclude : $exclude")
    println("Exiting")
    exit()
end

# Create an output directory that will be useful for the results.
out_base = config["out_base"]
if !ispath(out_base)
    println("Creating ", out_base)
    mkdir(out_base)
end


# vstart = parsed_args["vstart"]
# vend = parsed_args["vend"]
# nvel = parsed_args["nvel"]
#
# # If we have specified the velocities, ignore the data files and Doppler shift
# # and create an evenly spaced array of velocities
# vels = linspace(vstart, vend, nvel) # [km/s]
#
# lam0 = lam0s[species*transition]
#
# # convert velocities to wavelengths
# shift_lams = lam0 * (vels/c_kms + 1)

# Do everything common to both models
pars = convert_dict(config["parameters"], config["model"])
grd = config["grid"]
grid = Grid(grd)


# Are we using dust or gas?
# gas
if config["gas"]
    run(`cp $(assets_dir)radmc3d.inp.gas radmc3d.inp`)

    # Need this as a dummy file, apparently.
    run(`cp $(assets_dir)wavelength_micron.inp .`)

    # Copy appropriate molecule files.

    # Which species?
    species = config["species"]
    transition = config["transition"]
    if species == "12CO"
        run(`cp $(assets_dir)molecule_co.inp .`)
        run(`cp $(assets_dir)lines_co.inp lines.inp`)
    elseif species ==  "13CO"
        run(`cp $(assets_dir)molecule_13co.inp .`)
        run(`cp $(assets_dir)lines_13co.inp lines.inp`)
    elseif species == "C18O"
        run(`cp $(assets_dir)molecule_c18o.inp .`)
        run(`cp $(assets_dir)lines_c18o.inp lines.inp`)
    end

    println("Copied over RADMC-3D gas input files.")

    vel = pars.vel # [km/s]
    npix = config["npix"] # number of pixels

    vstart = parsed_args["vstart"]
    vend = parsed_args["vend"]
    nvel = parsed_args["nvel"]

    if vstart != nothing
        # If we have specified the velocities, ignore the data files and Doppler shift
        # and create an evenly spaced array of velocities
        vels = linspace(vstart, vend, nvel) # [km/s]

        lam0 = lam0s[species*transition]

        # convert velocities to wavelengths
        shift_lams = lam0 * (vels/c_kms + 1)
    else
        # read the wavelengths for all data channels
        fid = h5open(config["data_file"], "r")
        freqs = read(fid["freqs"]) # [Hz]
        # Convert from Hz to wavelengths in μm
        lams = cc ./freqs * 1e4 # [μm]
        close(fid)

        # Doppler shift the dataset wavelength according to the velocity in the parameter file
        beta = vel/c_kms # relativistic Doppler formula
        shift_lams =  lams .* sqrt((1. - beta) / (1. + beta)) # [microns]
    end

    # Doppler shift the dataset wavelength according to the velocity in the parameter file
    # beta = vel/c_kms # relativistic Doppler formula
    # shift_lams =  lams .* sqrt((1. - beta) / (1. + beta)) # [microns]

    # Write everything to the current working directory
    write_grid("", grid)
    write_model(pars, "", grid, species)
    write_lambda(shift_lams, "")
    println("Wrote gas model input files for RADMC-3D.")

# dust
elseif config["dust"]

    # Copy the appropriate files for dust synthesis
    run(`cp $(assets_dir)radmc3d.inp.dust radmc3d.inp`)

    # Write the dust distribution
    write_grid("", grid)
    write_dust(pars, "", grid)

    println("Wrote dust grid and density file for RADMC-3D.")


    # ff = assets_dir * "radmc3d.inp.dust"
    # run(`cp $ff radmc3d.inp`)
    # ff = assets_dir * "stars.inp"
    # run(`cp $ff .`)
    # ff = assets_dir * "wavelength_micron.inp"
    # run(`cp $ff .`)
    # ff = assets_dir * "dustkappa_silicate.inp"
    # run(`cp $ff .`)
    # ff = assets_dir * "dustopac.inp"
    # run(`cp $ff .`)
    #
    # println("Copied over RADMC-3D dust input files.")
end
