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
using model

import PyPlot.plt
using LaTeXStrings

pp = config["parameters"]
params = ["M_star", "r_c", "T_10", "q", "gamma", "logM_CO", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
starting_param = Array(Float64, nparam)

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

write_grid("")
write_model(pars, "")
write_lambda(shift_lams)

# Temporarily overwrite to dust mode
cp("radmc3d.inp.dust", "radmc3d.inp")
# NEEDS to be edited

# Write the dust density structure
# wavelength_micron.inp # this seems to be relevant for GAS too? It shouldn't.
# dust_density.inp # can we use a gas-to-dust ratio of 100?
# dustopac.inp
# dust_kappa_XXX.inp

if !parsed_args["norad"]
    run(`radmc3d therm`)
    run(`radmc3d image incl $incl posang $PA npix $npix`)
end

# When done, revert to gas only mode
cp("radmc3d.inp.gas", "radmc3d.inp")
