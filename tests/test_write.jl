push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/DiskJockey/")
# Given some model parameters, synthesize and plot the channel maps and
# integrated spectrum

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    # "--opt1"
    # help = "an option with an argument"
    # default = 0
    "--norad"
    help = "Use the image alread here."
    action = :store_true
    "config"
    help = "a YAML configuration file"
    required = true
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using constants
using image
using emodel
using HDF5

# read the wavelengths for all 23 channels
fid = h5open("../" * config["data_file"], "r")
lams = read(fid["lams"]) # [μm]
close(fid)


pp = config["parameters"]
params = ["M_star", "a_c", "T_10", "q", "gamma", "logM_CO", "ksi", "dpc", "incl", "PA", "e", "w", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
starting_param = Array{Float64}(nparam)

for i=1:nparam
    starting_param[i] = pp[params[i]][1]
end

# Convert logM_CO to M_CO
starting_param[6] = 10^starting_param[6]

pars = Parameters(starting_param...)

vel = pars.vel # [km/s]
# RADMC conventions for inclination and PA
incl = pars.incl # [deg]
PA = pars.PA # [deg] Position angle runs counter clockwise
npix = config["npix"] # number of pixels

# Doppler shift the dataset wavelength to rest-frame wavelength
beta = vel/c_kms # relativistic Doppler formula
shift_lams =  lams .* sqrt((1. - beta) / (1. + beta)) # [microns]

grd = config["grid"]
grid = Grid(grd["nr"], grd["ntheta"], grd["nphi"], grd["r_in"], grd["r_out"], grd["na"], true)

write_grid("", grid)
write_model(pars, "", grid)
write_lambda(shift_lams, "")

println("Write model")
tic()
write_model(pars, "", grid)
toc()

println("Write lams")
tic()
write_lambda(shift_lams, "")
toc()
