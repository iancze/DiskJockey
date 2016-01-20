#!/usr/bin/env julia

# The main purpose of this file is to read the config file and copy all of the relevant static
# files into this directory, write out the proper wavelength files, and make sure we can continue
# with the run properly.

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--test"
    help = "Print a test message to see you linked things correctly"
    action = :store_true
    "--new-config"
    help = "Copy a stock configuration file to this directory."
    action = :store_true
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

if parsed_args["test"]
    println("Your JudithExcalibur scripts are successfully linked.")
    println("Exiting")
    quit()
end

if parsed_args["new-config"]
    cp(assets_dir * "config.yaml", pwd() * "/config.yaml")
    println("Copied default config.yaml file to current working directory.")
    println("Exiting")
    quit()
end

import YAML
config = YAML.load(open(parsed_args["config"]))

using HDF5

using JudithExcalibur.constants
using JudithExcalibur.tmodel

# To initialize, we basically want to do everything that would be necessary to then go in and
# run a RADMC-3D model.

# Create an output directory that will be useful for the results.
out_base = config["out_base"]
if !ispath(out_base)
    println("Creating ", out_base)
    mkdir(out_base)
end

# Are we using dust or gas?

# gas
if config["gas"]
    ff = assets_dir * "radmc3d.inp.gas"
    run(`cp $ff radmc3d.inp`)

    # Need this as a dummy file, apparently.
    ff = assets_dir * "wavelength_micron.inp"
    run(`cp $ff .`)

    # Copy appropriate molecule files.

    # Which species?
    species = config["species"]
    transition = config["transition"]
    if species == "12CO"
        ff = assets_dir * "molecule_co.inp"
        run(`cp $ff .`)
        ff = assets_dir * "lines_co.inp"
        run(`cp $ff lines.inp`)
    elseif species ==  "13CO"
        ff = assets_dir * "molecule_13co.inp"
        run(`cp $ff .`)
        ff = assets_dir * "lines_13co.inp"
        run(`cp $ff lines.inp`)
    elseif species == "C18O"
        ff = assets_dir * "molecule_c18o.inp"
        run(`cp $ff .`)
        ff = assets_dir * "lines_c18o.inp"
        run(`cp $ff lines.inp`)
    end

    println("Copied over RADMC-3D gas input files.")

    # Now run the correct setup for the grid

    # If we want just want to generate a model without fitting to data, we can specify
    # the wavelengths arbitrarily and then use `plot_model.jl`

    # Load the starting parameters
    pp = config["parameters"]
    params = ["M_star", "r_in", "r_out", "T_10", "q", "gamma", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
    nparam = length(params)
    starting_param = Array(Float64, nparam)

    for i=1:nparam
        starting_param[i] = pp[params[i]][1]
    end

    # Convert logM_gas to M_gas
    starting_param[7] = 10^starting_param[7]


    pars = Parameters(starting_param...)

    vel = pars.vel # [km/s]
    # RADMC conventions for inclination and PA
    incl = pars.incl # [deg]
    PA = pars.PA # [deg] Position angle runs counter clockwise
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
    ff = assets_dir * "radmc3d.inp.dust"
    run(`cp $ff radmc3d.inp`)
    ff = assets_dir * "stars.inp"
    run(`cp $ff .`)
    ff = assets_dir * "wavelength_micron.inp"
    run(`cp $ff .`)
    ff = assets_dir * "dustkappa_silicate.inp"
    run(`cp $ff .`)
    ff = assets_dir * "dustopac.inp"
    run(`cp $ff .`)

    println("Copied over RADMC-3D dust input files.")
end
