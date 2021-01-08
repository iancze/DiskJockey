module model

export generate_vel_mask, write_grid, write_model, write_lambda, write_dust, Grid, size_au
export AbstractParameters, ParametersStandard, ParametersVertical, ParametersVerticalEta, ParametersNuker, convert_vector, convert_dict, registered_params
export lnprior

# The double dot is because we are now inside the model module, and we want to import the
# constants module, which is part of the enclosing DiskJockey package.
using ..constants
using Printf
using QuadGK
# using Dierckx
# using Optim
# using ODE

"Write the wavelength sampling file. Only run on setup"
function write_lambda(lams::AbstractArray, basedir::AbstractString)
    fcam = open(basedir * "camera_wavelength_micron.inp", "w")
    nlam = length(lams)
    @printf(fcam, "%d\n", nlam)
    for lam in lams
        @printf(fcam, "%.9e\n", lam) # [microns]
    end

    close(fcam)
end

"Generate a boolean mask corresponding to the velocities which fall inside the user-specified mask ranges.
arr is an array of [start, end] pairs, like arr = [[1.0, 3.0], [4.0, 5.2]]
vels are the velocities corresponding to the dataset, like vels = [1.2, 2.4, 3.6], etc."
function generate_vel_mask(arr, vels)
    # initialize a boolean array of all trues
    mask = trues(length(vels))

    # iterate over the array of exclude values
    for (left, right) in arr
        # set those velocities which fall inside the range of channels to be false
        mask[(left .< vels) .& (vels .< right)] .= false
    end

    return mask
end

const eqmirror = true # mirror the grid about the z=0 midplane ?
# if we decide to mirror, then ncells = 1/2 of the true value

"Define a grid object which stores all of these variables
This will not change for the duration of the run."
struct Grid
    nr::Int
    ntheta::Int
    nphi::Int
    ncells::Int
    # cell edges
    Rs::Vector{Float64} # [cm]
    Thetas::Vector{Float64}
    Phis::Vector{Float64}
    # cell centers
    rs::Vector{Float64} # [cm]
    thetas::Vector{Float64}
    phis::Vector{Float64}
end

"Hold the ray-tracing grid."
function Grid(nr::Int, ntheta::Int, r_in::Real, r_out::Real) #, eqmirror::Bool=true)
    # Specify a 2D axisymmetric *separable* grid in spherical coordinates:
    # {r, theta, phi}, where theta is angle from zenith, phi is azimuth

    # Number of cells in each dimension
    nphi = 1 # axisymmetric disk
    ncells = nr * ntheta * nphi
    r_in = convert(Float64, r_in) * AU # [cm] Inner extent of disk
    r_out = convert(Float64, r_out) * AU # [cm] Outer extent of disk

    #Define the cell *walls*
    # Rs = logspace(log10(r_in), log10(r_out), nr+1) # [cm] logarithmically spaced
    Rs = 10 .^ range(log10(r_in), stop=log10(r_out), length=nr+1) # [cm] logarithmically spaced

    eqmirror = true

    if eqmirror
        ped = 0.1
        #Thetas = LinRange(0, pi/2., ntheta+1)
        # [rad] Angles are internally defined in radians, not degrees
        # Thetas = pi/2. - (logspace(log10(ped), log10(pi/2. + ped), ntheta+1) - ped)[end:-1:1]
        Thetas = pi/2.0 .- (10 .^ range(log10(ped), stop=log10(pi/2.0 + ped), length=ntheta+1) .- ped)[end:-1:1]
        #Logarithmically spaced closer near the z=0
    else
        Thetas = LinRange(0, pi, ntheta+1)
        # [rad] Angles are internally defined in radians, not degrees
        # Equally spaced in theta.
    end

    Phis = Float64[0.0, 0.0] # [rad] cell walls for inactive coordinate

    #Define the cell centers as the average between walls
    rs = 0.5 * (Rs[1:end-1] + Rs[2:end]) # [cm]
    thetas = 0.5 * (Thetas[1:end-1] + Thetas[2:end])
    phis = Float64[0.0]

    return Grid(nr, ntheta, nphi, ncells, Rs, Thetas, Phis, rs, thetas, phis)

end

"Hold the ray-tracing grid."
function Grid(nr::Int, ntheta::Int, nphi::Int, r_in::Real, r_out::Real)
    # Specify a 2D axisymmetric *separable* grid in spherical coordinates:
    # {r, theta, phi}, where theta is angle from zenith, phi is azimuth

    # Number of cells in each dimension
    ncells = nr * ntheta * nphi
    r_in = convert(Float64, r_in) * AU # [cm] Inner extent of disk
    r_out = convert(Float64, r_out) * AU # [cm] Outer extent of disk

    #Define the cell *walls*
    # Rs = logspace(log10(r_in), log10(r_out), nr+1) # [cm] logarithmically spaced
    Rs = 10 .^ range(log10(r_in), stop=log10(r_out), length=nr+1) # [cm] logarithmically spaced

    eqmirror = true

    if eqmirror
        ped = 0.1
        #Thetas = LinRange(0, pi/2., ntheta+1)
        # [rad] Angles are internally defined in radians, not degrees
        # Thetas = pi/2. - (logspace(log10(ped), log10(pi/2. + ped), ntheta+1) - ped)[end:-1:1]
        Thetas = pi/2.0 .- (10 .^ range(log10(ped), stop=log10(pi/2. + ped), length=ntheta+1) .- ped)[end:-1:1]
        #Logarithmically spaced closer near the z=0
    else
        Thetas = LinRange(0, pi, ntheta+1)
        # [rad] Angles are internally defined in radians, not degrees
        # Equally spaced in theta.
    end

    Phis = LinRange(0, 2pi, nphi+1) # [rad] cell walls for inactive coordinate

    #Define the cell centers as the average between walls
    rs = 0.5 * (Rs[1:end-1] + Rs[2:end]) # [cm]
    thetas = 0.5 * (Thetas[1:end-1] + Thetas[2:end])
    phis = 0.5 * (Phis[1:end-1] + Phis[2:end])

    return Grid(nr, ntheta, nphi, ncells, Rs, Thetas, Phis, rs, thetas, phis)

end


"Create a grid object using a logarithmic then linear then logarithmic radial spacing"
function Grid(r_in::Real, r_linstart::Real, r_linend::Real, r_out::Real, n_in::Int, n_mid::Int, n_out::Int, ntheta::Int, eqmirror::Bool=true)
    # Number of cells in each dimension
    nphi = 1 # axisymmetric disk
    nr = n_in + n_mid + n_out
    ncells = nr * ntheta * nphi
    r_in = convert(Float64, r_in) * AU # [cm] Inner extent of disk
    r_linstart = convert(Float64, r_linstart) * AU # [cm] Start of linear section
    r_linend = convert(Float64, r_linend) * AU # [cm] End of linear section
    r_out = convert(Float64, r_out) * AU # [cm] Outer extent of disk

    #Define the cell *walls*
    # logarithmically spaced inner grid
    # Rs_in = logspace(log10(r_in), log10(r_linstart), n_in + 1) # [cm]
    Rs_in = 10 .^ range(log10(r_in), stop=log10(r_linstart), length=n_in + 1) # [cm]

    # linearly spaced middle grid
    Rs_mid = LinRange(r_linstart, r_linend, n_mid + 1) # [cm]

    # logarithmically spaced outer grid
    # Rs_out = logspace(log10(r_linend), log10(r_out), n_out + 1) # [cm]
    Rs_out = 10 .^ range(log10(r_linend), stop=log10(r_out), length=n_out + 1) # [cm]

    Rs = cat(1, Rs_in[1:end-1], Rs_mid, Rs_out[2:end])

    if eqmirror
        ped = 0.1
        #Thetas = LinRange(0, pi/2., ntheta+1)
        # [rad] Angles are internally defined in radians, not degrees
        # Thetas = pi/2. - (logspace(log10(ped), log10(pi/2. + ped), ntheta+1) - ped)[end:-1:1]
        Thetas = pi/2.0 .- (10 .^ range(log10(ped), stop=log10(pi/2. + ped), length=ntheta+1) - ped)[end:-1:1]
        #Spaced closer near the z=0
    else
        Thetas = LinRange(0, pi, ntheta+1)
        # [rad] Angles are internally defined in radians, not degrees
    end

    Phis = Float64[0.0, 0.0] # [rad] cell walls for inactive coordinate

    #Define the cell centers as the average between walls
    rs = 0.5 * (Rs[1:end-1] + Rs[2:end]) # [cm]
    thetas = 0.5 * (Thetas[1:end-1] + Thetas[2:end])
    phis = Float64[0.0]

    return Grid(nr, ntheta, nphi, ncells, Rs, Thetas, Phis, rs, thetas, phis)
end

"Read from a dictionary, then choose how to make the grid based upon the arguments."
function Grid(d::Dict)
    if "r_linstart" in keys(d)
        # We're going for a log-linear-log grid
        names = ["r_in", "r_linstart", "r_linend", "r_out", "n_in", "n_mid", "n_out", "ntheta"]
    elseif d["nphi"] > 1
        names = ["nr", "ntheta", "nphi", "r_in", "r_out"]
    else
        # We're just going for a log spaced grid
        names = ["nr", "ntheta", "r_in", "r_out"]
    end
    vec = [d[name] for name in names]
    # Unroll these into an actual parameter
    return Grid(vec...)
end

"This function only needs to be run once, upon setup."
function write_grid(basedir::AbstractString, grid::Grid)
    #amr_grid.inp
    f = open(basedir * "amr_grid.inp", "w")

    #Write the header
    @printf(f, "%d\n", 1) #iformat
    @printf(f, "%d\n", 0) #regular grid (no AMR or Oct-tree)
    @printf(f, "%d\n", 100) #spherical coordiantes
    @printf(f, "%d\n", 0) #gridinfo (none needed for now)
    #incl_r incl_phi incl_z #use this axis?
    if length(grid.Phis) > 2
        @printf(f, "%d %d %d \n", 1, 1, 1) # 2D axisymmetric
    else
        @printf(f, "%d %d %d \n", 1, 1, 0) # 2D axisymmetric
    end

    #n_r    n_phi   n_z #number of cells in this dimension
    @printf(f, "%d %d %d \n", grid.nr, grid.ntheta, grid.nphi)

    for R in grid.Rs
        @printf(f, "%.9e\n", R)
    end

    for Theta in grid.Thetas
        @printf(f, "%.9e\n", Theta)
    end

    for Phi in grid.Phis
        @printf(f, "%.9e\n", Phi)
    end

    close(f)
end


# Define an abstract Parameters type, then subset of parameters that can be used to dispatch specifics

abstract type AbstractParameters end

"Parameters for the standard model."
mutable struct ParametersStandard <: AbstractParameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] characteristic radius
    T_10::Float64 # [K] temperature at 10 AU
    q::Float64 # temperature gradient exponent
    gamma::Float64 # surface density gradient exponent
    log_M_gas::Float64 # [log10 M_sun] total disk mass
    ksi::Float64 # [cm s^{-1}] microturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
    Sigma_c::Float64 # [g/cm^2] surface density at r_c
end

# constructor for ParametersStandard takes all but Sigma_c
function ParametersStandard(M_star::Float64, r_c::Float64, T_10::Float64, q::Float64, gamma::Float64, log_M_gas::Float64, ksi::Float64, dpc::Float64,incl::Float64, PA::Float64, vel::Float64, mu_RA::Float64, mu_DEC::Float64)
  
  # calculate the normalization constant 
  Sigma_c::Float64 = 10^log_M_gas * M_sun * (2 - gamma) / (2 * pi * (r_c * AU)^2)

  # initialize the actual type
  return ParametersStandard(M_star, r_c, T_10, q, gamma, log_M_gas, ksi, dpc,incl, PA, vel, mu_RA, mu_DEC, Sigma_c)
end


"Parameters for the vertical model."
mutable struct ParametersVertical <: AbstractParameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] characteristic radius
    T_10m::Float64 # [K] temperature at 10 AU, midplane
    q_m::Float64 # midplane temperature gradient exponent
    T_10a::Float64 # [K] temperature at 10 AU, atmosphere
    q_a::Float64 # atmosphere temperature gradient exponent
    T_freeze::Float64 # [K] temperature below which to reduce CO abundance
    X_freeze::Float64 # [ratio] amount to reduce CO abundance
    sigma_s::Float64 # Photodissociation boundary in units of A_V.
    gamma::Float64 # surface temperature gradient exponent
    h::Float64 # Number of scale heights that z_q is at, typically fixed to 4
    delta::Float64 # Shape exponent, currently fixed to 2
    log_M_gas::Float64 # [log10 M_sun] total disk mass
    ksi::Float64 # [cm s^{-1}] micsroturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
    Sigma_c::Float64 # [g/cm^2] surface density at r_c
end


"Parameters for the vertical model with variable slope."
mutable struct ParametersVerticalEta <: AbstractParameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] characteristic radius
    T_10m::Float64 # [K] temperature at 10 AU, midplane
    q_m::Float64 # midplane temperature gradient exponent
    T_freeze::Float64 # [K] temperature below which to reduce CO abundance
    X_freeze::Float64 # [ratio] amount to reduce CO abundance
    sigma_s::Float64 # Photodissociation boundary in units of A_V.
    gamma::Float64 # surface density gradient exponent
    h::Float64 # Number of scale heights that z_q is at, typically fixed to 4
    eta::Float64 # Exponent for zq profile
    delta::Float64 # Shape exponent, currently fixed to 2
    log_M_gas::Float64 # [log10 M_sun] total disk mass
    ksi::Float64 # [cm s^{-1}] micsroturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
    Sigma_c::Float64 # [g/cm^2] surface density at r_c
end


"Parameters for the NUKER model."
mutable struct ParametersNuker <: AbstractParameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] characteristic radius
    T_10::Float64 # [K] temperature at 10 AU
    q::Float64 # temperature gradient exponent
    gamma::Float64 # surface density gradient within r_c (negative values yield holes)
    log_alpha::Float64 #  log10 sharpness of transition (2 = smooth, 16 = sharp)
    beta::Float64 # gradient power law outside r_c (~7)
    log_M_gas::Float64 # [log10 M_sun] total disk mass
    ksi::Float64 # [cm s^{-1}] microturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
    Sigma_c::Float64 # [g/cm^2] surface density at r_c
end


# constructor for ParametersStandard takes all but Sigma_c
function ParametersNuker(M_star::Real, r_c::Real, T_10::Real, q::Real, gamma::Real, log_alpha::Real, beta::Real, log_M_gas::Real, ksi::Real, dpc::Real,incl::Real, PA::Real, vel::Real, mu_RA::Real, mu_DEC::Real)

    if r_c <= 0.0
        throw(ModelException("r_c cannot be negative, $r_c"))
    end

    M_gas = 10^log_M_gas * M_sun
    alpha = 10^log_alpha

    function integrand(r)
        return (r/r_c)^(-gamma) * (1 + (r/r_c)^alpha)^((gamma - beta)/alpha) * r 
    end

    # calculate the normalization constant via an integral
    # break it into two pieces bounded at r_c, since this might help accuracy
    res, err = quadgk(integrand, 1e-3, r_c, 1e4) # [AU^2]
    # convert area to cm^2 
    res_cm = res * AU^2 # [cm^2]
    Sigma_c = M_gas / (2 * pi * res_cm) # [g/cm^2]

     # calculate the normalization constant
    # we need to do an integral 
    # M_gas = 2 * pi * Sigma_c * integrate((r/r_c)^(-gamma) * (1 + (r/r_c)^alpha)^((gamma - beta)/alpha) * r * dr)

    # initialize the actual type
    return ParametersNuker(M_star, r_c, T_10, q, gamma, log_alpha, beta, log_M_gas, ksi, dpc,incl, PA, vel, mu_RA, mu_DEC, Sigma_c)
end

"A dictionary of parameter lists for conversion."
registered_params = Dict([("standard", ["M_star", "r_c", "T_10", "q", "gamma", "log_M_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]),
("vertical", ["M_star", "r_c", "T_10m", "q_m", "T_10a", "q_a", "T_freeze", "X_freeze", "sigma_s", "gamma", "h", "delta", "log_M_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]),
("verticalEta", ["M_star", "r_c", "T_10m", "q_m", "T_freeze", "X_freeze", "sigma_s", "gamma", "h", "eta", "delta", "log_M_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]),
("nuker", ["M_star", "r_c", "T_10", "q", "gamma", "log_alpha", "beta", "log_M_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"])])

registered_types = Dict([("standard", ParametersStandard), ("vertical", ParametersVertical), ("verticalEta", ParametersVerticalEta), ("nuker", ParametersNuker)])

"Unroll a vector of parameter values into a parameter type."
function convert_vector(p::Vector{Float64}, model::AbstractString, fix_params::Vector; args...)
    args = Dict{Symbol}{Float64}(args)

    # The goal is to assemble a length-nparam vector that can be unrolled into a parameter type
    # e.g., ParametersStandard(M_star, r_c, T_10, q, gamma, log_M_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

    # Select the registerd parameters corresponding to this model
    # These are the names listed in this file (model.jl)
    reg_params = registered_params[model]

    # fix_params is a list of strings from a config file that list the names of parameters to be fixed to the config
    # value. These names correspond to the registered names.

    # fit_params are the ones in reg_params that are not in (∉) fix_params
    fit_params = filter(x->∉(x,fix_params), reg_params)

    nparams = length(reg_params)

    # Make an empty vector of this same length
    par_vec = Array{Float64}(undef, nparams)

    # This requires assigning p to fit_params
    # Find the indexes that correspond to fit_params
    # par_indexes = findin(reg_params, fit_params)
    par_indexes = findall((in)(fit_params), reg_params)
    # Stuff p directly into these
    par_vec[par_indexes] = p

    # Then reading fix_params from args
    # First, create an array of fixed values analogous to p
    p_fixed = Float64[args[Symbol(par)] for par in fix_params]
    par_indexes = findall((in)(fix_params), reg_params)
    par_vec[par_indexes] = p_fixed

    # Then assembling these in the same orignial order as registered_params, into the parameter
    # type corresponding to the model.
    parameters = registered_types[model](par_vec...)

    return parameters

end

"Used to turn a dictionary of parameter values (from config.yaml) directly into a parameter type. Generally used for synthesis and plotting command line scripts."
function convert_dict(p::Dict, model::AbstractString)
    # Select the registerd parameters corresponding to this model
    reg_params = registered_params[model]
    nparams = length(reg_params)

    # Using this order of parameters, unpack the dictionary p into a vector
    # reg_params reads Sigma_c, not logSigma_c; alpha, not logalpha
    par_vec = Float64[p[par_name] for par_name in reg_params]

    return registered_types[model](par_vec...)

end


"The common sense priors that apply to all parameter values"
function lnprior_base(pars::AbstractParameters)
    # Create a giant short-circuit or loop to test for sensical parameter values.
    if pars.M_star <= 0.0 || pars.ksi <= 0.0 || pars.T_10 <= 0.0 || pars.r_c <= 0.0  || pars.T_10 > 1500. || pars.q < 0. || pars.q > 1.0 || pars.incl < 0. || pars.incl > 180. || pars.PA < -180. || pars.PA > 520. 
        # println("M_star ", pars.M_star)
        # println("r_c ", pars.r_c)
        # println("T_10 ", pars.T_10)
        # println("q ", pars.q)
        # println("incl ", pars.incl)
        # println("PA ", pars.PA)
        # println("Exiting in lnprior base")
        throw(ModelException("Parameters outside of prior range $pars"))
    end

    # If we've passed all the hard-cut offs by this point, return the geometrical inclination prior.
    # natural log, not log10
    return log(0.5 * sind(pars.incl))

end

function lnprior(pars::ParametersStandard, grid::Grid)

    lnp = lnprior_base(pars)

    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.
    if (3 * pars.r_c) > r_out
        throw(ModelException("Model radius ($(pars.r_c)) too large for grid size ($r_out)"))
    else
        return lnp
    end

end

function lnprior(pars::ParametersVertical, grid::Grid)
    # Create a giant short-circuit or loop to test for sensical parameter values.
    if pars.M_star <= 0.0 || pars.ksi <= 0.0 || pars.T_10a <= 0.0 || pars.T_10m <= 0.0 || pars.r_c <= 0.0  || pars.T_10a > 1500.0 || pars.q_m < 0.0 || pars.q_a < 0.0 || pars.q_m > 1.0 || pars.q_a > 1.0 || pars.incl < 0. || pars.incl > 180.0 || pars.PA < -180.0 || pars.PA > 520.0 || pars.X_freeze > 1.0 || pars.sigma_s < 0.0
        throw(ModelException("Parameters outside of prior range $pars"))
    end

    # Check to see that the temperatures make sense
    if pars.T_10m > pars.T_10a
        throw(ModelException("Atmosphere ($(pars.T_10a)) is cooler than midplane ($(pars.T_10m)."))
    end


    # If we've passed all the hard-cut offs by this point, return the geometrical inclination prior.
    # natural log, not log10
    lnp = log(0.5 * sind(pars.incl))

    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.
    if (3 * pars.r_c) > r_out
        throw(ModelException("Model radius ($(pars.r_c)) too large for grid size ($r_out)."))
    else
        return lnp
    end

end

function lnprior(pars::ParametersVerticalEta, grid::Grid)
    # Create a giant short-circuit or loop to test for sensical parameter values.
    if pars.M_star <= 0.0 || pars.ksi <= 0.0 || pars.T_10m <= 60.0 || pars.r_c <= 0.0  || pars.q_m < 0.0 || pars.q_m > 1.0 || pars.incl < 0.0 || pars.incl > 180.0 || pars.PA < -180.0 || pars.PA > 520.0 || pars.X_freeze > 1.0 || pars.sigma_s < 0.0 || pars.eta < 0.2 || pars.eta > 6.0  || pars.h < 0.5 || pars.h > 6.0 || pars.delta < 0.5 || pars.delta > 6.0 || pars.gamma < 0.5 || pars.gamma > 3.0
        throw(ModelException("Parameters outside of prior range $pars"))
    end


    # If we've passed all the hard-cut offs by this point, return the geometrical inclination prior.
    # natural log, not log10
    lnp = log(0.5 * sind(pars.incl))

    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.
    if (3 * pars.r_c) > r_out
        throw(ModelException("Model radius ($(pars.r_c)) too large for grid size ($r_out)."))
    else
        return lnp
    end

end


function lnprior(pars::ParametersNuker, grid::Grid)

    # println("In lnprior(Nuker)")
    lnp = lnprior_base(pars)

    if (pars.log_alpha < 0.0 ) || (pars.log_alpha > 2.0) || (pars.beta < 2) || (pars.beta > 10)
        # println("alpha: ", pars.alpha)
        # println("beta: ", pars.beta)
        # println("Alpha or Beta outside of prior range.")
        throw(ModelException("log_alpha ($(pars.log_alpha)) or beta ($(pars.beta)) outside of prior range."))
    end

    # Add on the prior on gamma (log, not log10)
    lnp_gamma = log(1/(1 + exp(-5 * (pars.gamma + 3))) - 1/(1 + exp(-15 * (pars.gamma - 2))))
    lnp += lnp_gamma

    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.
    if (3 * pars.r_c) > r_out
        # println("Model radius too large for grid size.")
        throw(ModelException("Model radius ($(pars.r_c)) too large for grid size ($r_out)."))
    else
        # println("Returning lnprior as ", lnp)
        return lnp
    end

end


# Determine the physical size of the image. Size_arcsec is the full width of the image,
# and so sizeau is the full size of the image as well (RADMC3D conventions).
# Because there seems to be a small but constant shift in the size of the image
function size_au(size_arcsec::Real, dpc::Real, grid::Grid)
    outer_model = (1.1 * grid.rs[end] / AU) # [AU]
    sizeau_desired = size_arcsec * dpc # [AU]

    if sizeau_desired < outer_model
        # We don't want this error to be caught by the lnprob routine, since this should
        # always halt execution of the script.
        throw(ErrorException("Image size ($size_au AU) must be larger than 110% the model grid size ($outer_model AU)."))
    end

    # because there seems to be a slight offset between what sizeau is specified and the actual
    # size of the image, we will also specifiy sizeau_command, which is the value to give to RADMC
    sizeau_command = sizeau_desired/(1.0 + RADMC_SIZEAU_SHIFT)

    return (sizeau_desired, sizeau_command) # [AU]

end

# Assume all inputs to these functions are in CGS units and in *cylindrical* coordinates.
# Parametric type T allows passing individual Float64 or Vectors.
# # Alternate functions accept pars passed around, where pars is in M_star, AU, etc...
function velocity(r::T, M_star::Float64) where {T}
    sqrt.(G * M_star ./ r)
end
velocity(r::T, pars::AbstractParameters) where {T} = velocity(r, pars.M_star * M_sun)

# For the vertical temperature gradient
function velocity(r::Float64, z::Float64, M_star::Float64)
    sqrt.(G * M_star / (r^2 + z^2)^(3.0/2)) * r
end
velocity(r::Float64, z::Float64, pars::ParametersVertical) = velocity(r, z, pars.M_star * M_sun)
velocity(r::Float64, z::Float64, pars::ParametersVerticalEta) = velocity(r, z, pars.M_star * M_sun)

function temperature(r::T, T_10::Float64, q::Float64) where {T}
    T_10 * (r ./ (10.0 * AU)).^(-q)
end
temperature(r::T, pars::AbstractParameters) where {T} = temperature(r, pars.T_10, pars.q)

function Hp(r::T, M_star::Float64, T_10::Float64, q::Float64) where {T}
    temp = temperature(r, T_10, q)
    sqrt.(kB * temp .* r.^3.0/(mu_gas * m_H * G * M_star))
end
Hp(r::T,  pars::AbstractParameters) where {T} = Hp(r, pars.M_star * M_sun, pars.T_10, pars.q)
# Scale height computed from midplane temperature for vertical temperature gradient model
Hp(r::T,  pars::Union{ParametersVertical,ParametersVerticalEta}) where {T} = Hp(r, pars.M_star * M_sun, pars.T_10m, pars.q_m)

# For the vertical temperature gradient model
T_mid(r::T, pars::ParametersVertical) where {T} = temperature(r, pars.T_10m, pars.q_m)
T_atm(r::T, pars::ParametersVertical) where {T} = temperature(r, pars.T_10a, pars.q_a)

T_mid(r::T, pars::ParametersVerticalEta) where {T} = temperature(r, pars.T_10m, pars.q_m)

function T_atm(r::T, pars::ParametersVerticalEta) where {T}
    Tm = T_mid(r, pars)

    # the atmosphere temperature is the maximum of 500 or twice the midplane
    # in order to get the structure close to the star correct
    Ta = maximum([500, 2 * Tm])
    return Ta
end


# Atmosphere height computed as multiple of scale height (computed at midplane)
function z_q(r::T, pars::ParametersVertical) where {T}
    return pars.h * Hp(r, pars)
end

function z_q(r::T, pars::ParametersVerticalEta) where {T}

    # Calculate the scale height at this radius
    H = Hp(r, pars) # [cm]
    r0 = 10 * AU # [cm]

    if r < r0
      return pars.h * H
    else
      return pars.h * H * (r / r0)^pars.eta
    end
end

function temperature(r::Float64, z::Float64, pars::ParametersVertical)
    zq = z_q(r, pars)
    Ta = T_atm(r, pars)
    Tm = T_mid(r, pars)
    if Tm > Ta
      return Ta
    end

    if z >= zq
        return Ta
    else
        # return Ta + (Tm - Ta) * cos(pi * z/ (2 * zq))^pars.delta # Dartois
        return Tm + (Ta - Tm) * sin(pi * z/ (2 * zq))^(2 * pars.delta) # Williams and Best 14
    end
end

function temperature(r::Float64, z::Float64, pars::ParametersVerticalEta)
    zq = z_q(r, pars)
    Tm = T_mid(r, pars)
    Ta = T_atm(r, pars)

    H = Hp(r, pars) # scale height

    # If we are less than one scale height from the midplane, just return the midplane temperature
    if z < H
      return Tm
    end

    if Tm > Ta
      return Ta
    end

    if z >= zq
        return Ta
    else
        return Tm + (Ta - Tm) * sin(pi * (z - H)/ (2 * zq))^(2 * pars.delta)
    end
end


"Calculate the gas surface density"
function Sigma(r::Float64, pars::Union{ParametersStandard, ParametersVertical, ParametersVerticalEta})
    r_c = pars.r_c * AU

    gamma = pars.gamma
    Sigma_c = pars.Sigma_c

    S = Sigma_c * (r/r_c)^(-gamma) * exp(-(r/r_c)^(2 - gamma))

    return S
end


"
    Sigma(r::Float64, pars::ParametersNuker)

Calculate the gas surface density using the Nuker profile."
function Sigma(r::Float64, pars::ParametersNuker)
    r_c = pars.r_c * AU

    gamma = pars.gamma
    Sigma_c = pars.Sigma_c

    alpha = 10^pars.log_alpha
    beta = pars.beta

    S = Sigma_c * (r/r_c)^(-gamma) * (1 + (r/r_c)^alpha)^((gamma - beta)/alpha)

    return S
end


# Delivers a gas density in g/cm^3
function rho_gas(r::Float64, z::Float64, pars::AbstractParameters)
    H = Hp(r, pars)
    S = Sigma(r, pars)

    # Calculate the density
    rho = S/(sqrt(2.0 * pi) * H) * exp(-0.5 * (z/H)^2)

    return rho
end


# Determines the midplane density assuming vertically isothermal
function rho_gas_mid(r::Float64, pars::AbstractParameters)
    H = Hp(r, pars)
    S = Sigma(r, pars)

    # Calculate the density at the midplane
    rho = S/(sqrt(2.0 * pi) * H)

    return rho

end

# The temperature change, dT/dz
function dT(r::Float64, z::Float64, pars::ParametersVertical)
    zq = z_q(r, pars)
    Ta = T_atm(r, pars)

    if z >= zq
        return 0.0
    else
        Tm = T_mid(r, pars)
        # return - pi * pars.delta / (2 * zq) * (Tm - Ta) * (cos(pi * z/(2 * zq)))^(pars.delta - 1) * sin(pi * z / (2 * zq))
        return pars.delta * pi/zq * (Ta - Tm) * cos(pi * z/(2 * zq)) * (sin(pi * z/(2 * zq)))^(2 * pars.delta - 1)
    end
end

# The temperature change, dT/dz
function dT(r::Float64, z::Float64, pars::ParametersVerticalEta)

    zq = z_q(r, pars)
    Ta = T_atm(r, pars)

    H = Hp(r, pars)

    if z >= zq
        return 0.0
    elseif z < H
      return 0.0
    else
        Tm = T_mid(r, pars)

        return pars.delta * pi/zq * (Ta - Tm) * cos(pi * (z - H)/(2 * zq)) * (sin(pi * (z - H)/(2 * zq)))^(2 * pars.delta - 1)
    end
end

# The integrand
function dlnrho(r::Float64, z::Float64, pars::Union{ParametersVertical,ParametersVerticalEta})
    T = temperature(r, z, pars)
    dTemp = dT(r, z, pars)
    return -(dTemp/T + (mu_gas * m_H)/(kB * T) * (G * pars.M_star * M_sun * z)/(r^2 + z^2)^1.5)
end


# The integrand for a vertically isothermal testcase
function dlnrho(r::Float64, z::Float64, pars::ParametersStandard)
    return -(mu_gas * m_H)/(kB * T) * (G * pars.M_star * M_sun * z)/r^3
end

function z_top(r::Float64, pars::Union{ParametersVertical, ParametersVerticalEta})
    zq = z_q(r, pars) # The height of the disk "atmosphere"

    return zq * 5
end


# function rho_gas(r::Float64, z::Float64, pars::Union{ParametersVertical, ParametersVerticalEta})
#
#     # Calculate the "top" of the atmosphere,
#     # a height where we are sure we can enforce that density = 0
#     ztop = z_top(r, pars)
#
#     # If we are querying a height above this, just return 0 now
#     if z > ztop
#         return constants.rho_gas_zero
#     end
#
#     # Calculate what the maximum surface density would be if we integrated all the way to  the midplane
#     sigma = Sigma(r, pars) / 2
#
#     # Calculate the photodissociation height
#     # Threshold column density
#     thresh_H2 = pars.sigma_s * Av_sigmaH / 2 # [n_H2/cm^2]
#     thresh = thresh_H2 * (mu_gas * amu / X_H2) # [g/cm^2] of gas
#
#     # If there is not even enough column density at the midplane in order to exceed the threshold
#     # to protect against photodissociation of CO, just return 0 now
#     # if sigma < thresh
#     #     return constants.rho_gas_zero
#     # end
#
#     # Create an array of heights that goes from the midplane (0, 0.001, ..., ztop)
#     zs = cat(1, [0], logspace(log10(0.001 * AU), log10(ztop), 64 - 1))
#
#     # Define the ODE at this radius r
#     function f(z, y)
#         # calculate y'
#         dlnrho(r, z, pars) * y
#     end
#
#     # Choose a starting guess corresponding to the midplane density for an isothermal model.
#     # This helps us get in the ballpark and avoid a lot of the numerical errors that result
#     start = rho_gas_mid(r, pars)
#
#     # Now solve the ODE on the grid of postulated values
#     # For the inner disk, we require a larger value of abstol
#     # We need to experiment with what is best
#     abstol=1e-22
#     z_out, y_out = ode45(f, start, zs; abstol=1e-22)
#
#     # where y_out is less than the minimum accuracy we can trust, just set it to the minimum value
#     y_out[y_out .< abstol] = abstol
#
#     # Integrate the output from the ODE solver to find the total (un-normalized) surface density
#     spl = Spline1D(z_out, y_out)
#     tot = integrate(spl, z_out[1], z_out[end])
#
#     # Apply a correction factor to the un-normalized densities
#     cor = sigma / tot
#     rhos = y_out .* cor
#
#     # This is the function that will be minimized
#     function g(zvar::Float64)
#         return abs(thresh - cor * integrate(spl, zvar, z_out[end]))
#     end
#
#     z_phot = optimize(g, z_out[1], z_out[end]).minimum
#
#     # Evaluate rho at this point
#     rho = cor * evaluate(spl, z)
#
#     # If we are above this height, return the gas density reduced by a factor of 100
#     if z > z_phot
#         rho = 1e-2 * rho
#     end
#
#     if rho < 0.0
#       println("r ", r/AU, " z ", z/AU, " z_top ", ztop/AU, " z_phot ", z_phot/AU, " y_out ", y_out)
#       return (start, z_out, y_out)
#       throw(ModelException("Rho less than 0.0"))
#     end
#
#     return rho
# end

# Now, replace these functions to simply multiply rho_gas by X_12CO/m_12CO, or X_13CO/m_13CO, etc.
n_12CO(r::Float64, z::Float64, pars::AbstractParameters) = number_densities["12CO"] * rho_gas(r, z, pars)

n_13CO(r::Float64, z::Float64, pars::AbstractParameters) = number_densities["13CO"] * rho_gas(r, z, pars)

n_C18O(r::Float64, z::Float64, pars::AbstractParameters) = number_densities["C18O"] * rho_gas(r, z, pars)

# It is realistic to include freezout of CO onto dust grains.
# This is the amount by which the number density of the CO is reduced (X_freeze) relative to
# the nominal value.
function X_freeze(temp::Float64, pars::AbstractParameters)
    # If it's cooler than the freezout temperature, reduce the number density by the given factor
    if temp <= pars.T_freeze
        return pars.X_freeze
    # Otherwise, just keep it as is.
    else
        return 1.0
    end
end

function rho_dust(r::Float64, z::Float64, pars::AbstractParameters)

    # Use the rho_gas function to get the total density in [g/cm^3]
    mGas = rho_gas(r, z, pars)

    # Convert from mGas to mDust using Gas/Dust ratio of 100
    mDust = mGas * 0.01 # [g]

    return mDust
end

# Ksi is microturbulent broadining width in units of km/s. Output of this function
# is in cm/s for RADMC (RADMC manual, eqn 7.12)
function microturbulence(ksi::Float64)
    return ksi * 1.e5 # convert from km/s to cm/s
end

microturbulence(pars::AbstractParameters) = microturbulence(pars.ksi)

function write_model(pars::AbstractParameters, basedir::AbstractString, grid::Grid, species::AbstractString)

    funcs = Dict([("12CO", n_12CO), ("13CO", n_13CO), ("C18O", n_C18O)])
    n_CO = funcs[species]

    # numberdens_co.inp
    fdens = open(basedir * "numberdens_" * molnames[species] * ".inp", "w")
    @printf(fdens, "%d\n", 1) #iformat
    @printf(fdens, "%d\n", grid.ncells)

    # gas_velocity.inp
    fvel = open(basedir * "gas_velocity.inp", "w")
    @printf(fvel, "%d\n", 1) #iformat
    @printf(fvel, "%d\n", grid.ncells)

    # gas_temperature.inp
    ftemp = open(basedir * "gas_temperature.inp", "w")
    @printf(ftemp, "%d\n", 1) #iformat
    @printf(ftemp, "%d\n", grid.ncells)

    # microturbulence.inp
    fmicro = open(basedir * "microturbulence.inp", "w")
    @printf(fmicro, "%d\n", 1) #iformat
    @printf(fmicro, "%d\n", grid.ncells)

    # Now, we will need to write the three other files as a function of grid position.
    # Therefore we will do *one* loop over these indices, calculate the required value,
    # and write it to the appropriate file.

    #Looping over the cell centers
    for phi in grid.phis
        for theta in grid.thetas
            for r in grid.rs
                #Convert from spherical to cylindrical coordinates
                z = r * cos(theta)
                r_cyl = r * sin(theta)

                @printf(fdens, "%.9e\n", n_CO(r_cyl, z, pars))
                @printf(fvel, "0 0 %.9e\n", velocity(r_cyl, pars))
                @printf(ftemp, "%.9e\n", temperature(r_cyl, pars))
                @printf(fmicro, "%.9e\n", microturbulence(pars))
            end
        end
    end

    close(fdens)
    close(fvel)
    close(ftemp)
    close(fmicro)

end

function P_x(var)
    mat =  Float64[[1, 0, 0]  [0, cos(var), -sin(var)] [0, sin(var), cos(var)]]
    return mat
end

function P_z(var)
    mat =   Float64[[cos(var), -sin(var), 0] [sin(var), cos(var),   0] [0,        0,          1]]
    return mat
end

# function to transform spherical to cartesian
function P_project(theta, phi)
    mat = Float64[ [sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)] [cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta)] [-sin(phi), cos(phi), 0]]
    return mat
end


# # Get the components of the vector in RADMC coordinates
# function velocity_inner(pars, r, theta, phi)
#
#     Pproj = P_project(theta, phi)
#     # define rotation matrices
#     Px = P_x(pars.incl_inner * pi/180)
#     Pz = P_z(pars.Omega * pi/180)
#
#
#     # get primed coordinates
#     cart, cart_prime, sphere_prime = convert_position(pars, r, theta, phi)
#
#     # query the velocity field of the inner disk, in the primed frame,
#     # which will be projected along the hat{\phi^\prime} direction.
#     # for now, just call this v_phi_prime
#     r_prime, theta_prime, phi_prime = sphere_prime
#     r_cyl_prime = r_prime * sin(theta_prime)
#     z_cyl_prime = r_prime * cos(theta_prime)
#     vel_sphere_prime = Float64[0.0, 0.0, velocity(r_cyl_prime, z_cyl_prime, pars)]
#
#     # project it to the primed cartesian unit vectors
#     vel_cart_prime = Pproj * vel_sphere_prime
#
#     # now rotate this from the primed cartesian frame to the RADMC-3D cartesian frame
#     vel_cart = Pz * Px * vel_cart_prime
#
#     # now project it on to the spherical unit vectors in the RADMC-3D frame
#     vel_sphere = Pproj.' * vel_cart
#
#     # assert that the norm of all of these vectors is always zero
#     @assert isapprox(norm(vel_cart_prime), norm(vel_sphere_prime)) "Cartesian norm doesn't match "
#     @assert isapprox(norm(vel_cart), norm(vel_sphere_prime)) "Rotated cartesion norm doesn't match"
#     @assert isapprox(norm(vel_sphere), norm(vel_sphere_prime)) "Spherical vel doesn't match"
#
#     return vel_sphere
#
# end




function write_model(pars::Union{ParametersVertical, ParametersVerticalEta}, basedir::AbstractString, grid::Grid, species::AbstractString)

    function n_CO(r_cyl, z)
        return number_densities[species] * rho_gas(r_cyl, z, pars)
    end

    # numberdens_co.inp
    fdens = open(basedir * "numberdens_" * molnames[species] * ".inp", "w")
    @printf(fdens, "%d\n", 1) #iformat
    @printf(fdens, "%d\n", grid.ncells)

    # gas_velocity.inp
    fvel = open(basedir * "gas_velocity.inp", "w")
    @printf(fvel, "%d\n", 1) #iformat
    @printf(fvel, "%d\n", grid.ncells)

    # gas_temperature.inp
    ftemp = open(basedir * "gas_temperature.inp", "w")
    @printf(ftemp, "%d\n", 1) #iformat
    @printf(ftemp, "%d\n", grid.ncells)

    # microturbulence.inp
    fmicro = open(basedir * "microturbulence.inp", "w")
    @printf(fmicro, "%d\n", 1) #iformat
    @printf(fmicro, "%d\n", grid.ncells)

    # Now, we will need to write the three other files as a function of grid position.
    # Therefore we will do *one* loop over these indices, calculate the required value,
    # and write it to the appropriate file.

    #Looping over the cell centers
    for phi in grid.phis
        for theta in grid.thetas
            for r in grid.rs
                #Convert from spherical to cylindrical coordinates
                z = r * cos(theta)
                r_cyl = r * sin(theta)

                temp = temperature(r_cyl, z, pars)

                @printf(fdens, "%.9e\n", X_freeze(temp, pars) * n_CO(r_cyl, z))
                @printf(fvel, "0 0 %.9e\n", velocity(r_cyl, z, pars))
                @printf(ftemp, "%.9e\n", temp)
                @printf(fmicro, "%.9e\n", microturbulence(pars))
            end
        end
    end

    close(fdens)
    close(fvel)
    close(ftemp)
    close(fmicro)

end

function write_dust(pars::AbstractParameters, basedir::AbstractString, grid::Grid)
    fdens = open(basedir * "dust_density.inp", "w")
    @printf(fdens, "%d\n", 1) #iformat
    @printf(fdens, "%d\n", grid.ncells)
    @printf(fdens, "%d\n", 1) # number of dust species

    for phi in grid.phis
        for theta in grid.thetas
            for r in grid.rs
                #Convert from spherical to cylindrical coordinates
                z = r * cos(theta)
                r_cyl = r * sin(theta)

                @printf(fdens, "%.9e\n", rho_dust(r_cyl, z, pars))
            end
        end
    end

    close(fdens)
end


end
