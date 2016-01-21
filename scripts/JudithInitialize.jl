#!/usr/bin/env julia

# The main purpose of this file is to read the config file and copy all of the relevant static
# files into this directory, write out the proper wavelength files, and make sure we can continue
# with the run properly.

using ArgParse

s = ArgParseSettings(description="Initialize a new project directory with the appropriate files. Can also be used to update RADMC-3D input files after making changes to config.yaml")
@add_arg_table s begin
    "--version"
    help = "Print a the version number and exit."
    action = :store_true
    "--new-project"
    help = "Copy a stock configuration file to this directory. Can be 'standard', 'truncated', 'cavity'"
    default = "no"
    "--config"
    help = "a YAML configuration file"
    default = "config.yaml"
    # Options for specifying only model velocities.
    "--vstart", "-s"
    help = "Starting velocity in km/s."
    arg_type = Float64
    "--vend", "-e"
    help = "Ending velocity in km/s."
    arg_type = Float64
    "--nvel", "-n"
    help = "Number of velocity channels."
    arg_type = Int
    "--vsys"
    help = "Shift the systemic velocity"
    arg_type = Float64
end

parsed_args = parse_args(ARGS, s)

assets_dir = Pkg.dir("JudithExcalibur") * "/assets/"


# The user is going to start modeling a new disk, so copy in the new configuration file.
if parsed_args["new-project"] != "no"
    model = parsed_args["new-project"]

    cp(assets_dir * "config.$(model).yaml", pwd() * "/config.yaml")
    cp(assets_dir * "InitializeWalkers.$(model).ipynb", pwd() * "/InitializeWalkers.ipynb")

    println("Copied default model specific config.yaml and InitializeWalkers.ipynb to current working directory.")
    println("Exiting")
    quit()
end

import YAML
config = YAML.load(open(parsed_args["config"]))

using HDF5

using JudithExcalibur.constants
using JudithExcalibur.model

# This is just for new users to test that the package is successfully installed.
if parsed_args["version"]
    println("Your JudithExcalibur scripts are successfully linked.")
    println("You are running JudithExcalibur $JUDITHEXCALIBUR_VERSION")
    println("Exiting")
    quit()
end

# Create an output directory that will be useful for the results.
out_base = config["out_base"]
if !ispath(out_base)
    println("Creating ", out_base)
    mkdir(out_base)
end

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

    pars = convert_dict(config["parameters"], config["model"])
    vel = pars.vel # [km/s]
    incl = pars.incl # [deg]
    PA = pars.PA # [deg]
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
        lams = read(fid["lams"]) # [μm]
        close(fid)

        # Doppler shift the dataset wavelength according to the velocity in the parameter file
        beta = vel/c_kms # relativistic Doppler formula
        shift_lams =  lams .* sqrt((1. - beta) / (1. + beta)) # [microns]
    end

    grd = config["grid"]
    grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)

    # Write everything to the current working directory
    write_grid("", grid)
    write_model(pars, "", grid, species)
    write_lambda(shift_lams, "")
    println("Wrote gas model input files for RADMC-3D.")

# dust
else
    throw(ErrorException("Dust features not yet implemented. Only gas."))
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
