# Test framework to see if config files and fixed/free parameters work

using DiskJockey.model

import YAML
config = YAML.load(open("config.yaml"))

# Read parameters from a config file
pars = config["parameters"]
model = config["model"]

# Test the dictionary conversion routine
pars_type = convert_dict(pars, model)
println(pars_type)
println(pars_type.Sigma_c)

# Test the vector conversion routine

# Convert pars dictionary to dictionary of symbols
xargs = convert(Dict{Symbol}{Float64}, pars)


# According to model.jl, these are the total parameters. Let's fix some
# "standard", ["M_star", "r_c", "T_10", "q", "gamma", "Sigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]

# Separate out the fixed values
fixed = ["gamma"]

# Separate out the free values
free = ["M_star", "r_c", "T_10", "q", "Sigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]

p = [1.0, 50.0, 40.0, 0.6, -1.0, 0.4, 145.0, 130.0, 130.0, 0.0, 0.0, 0.0]

# Convert a proposal (vector) of free parameters with to a parameter type using fixed parameters to supply missing fields.
# Convert these to a parameter type
pars_type = convert_vector(p, model, fixed; xargs...)

println(pars_type)
println(pars_type.Sigma_c)
