module model

export write_grid, write_model, write_lambda, write_dust, Grid, size_au
export AbstractParameters, ParametersStandard, ParametersTruncated, ParametersCavity, ParametersVertical, ParametersVerticalEta, convert_vector, convert_dict
export lnprior

# The double dot is because we are now inside the model module, and we want to import the
# constants module, which is part of the enclosing DiskJockey package.
using ..constants
using Dierckx
using Optim
using ODE

# Write the wavelength sampling file. Only run on setup
function write_lambda(lams::AbstractArray, basedir::AbstractString)
    fcam = open(basedir * "camera_wavelength_micron.inp", "w")
    nlam = length(lams)
    @printf(fcam, "%d\n", nlam)
    for lam in lams
        @printf(fcam, "%.9e\n", lam) # [microns]
    end

    close(fcam)
end

const eqmirror = true # mirror the grid about the z=0 midplane ?
# if we decide to mirror, then ncells = 1/2 of the true value

# Define a grid object which stores all of these variables
# This will not change for the duration of the run
immutable Grid
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

function Grid(nr::Int, ntheta::Int, r_in::Real, r_out::Real, eqmirror::Bool=true)
    # Specify a 2D axisymmetric *separable* grid in spherical coordinates:
    # {r, theta, phi}, where theta is angle from zenith, phi is azimuth

    # Number of cells in each dimension
    nphi = 1 # axisymmetric disk
    ncells = nr * ntheta * nphi
    r_in = convert(Float64, r_in) * AU # [cm] Inner extent of disk
    r_out = convert(Float64, r_out) * AU # [cm] Outer extent of disk

    #Define the cell *walls*
    Rs = logspace(log10(r_in), log10(r_out), nr+1) # [cm] logarithmically spaced

    if eqmirror
        ped = 0.1
        #Thetas = linspace(0, pi/2., ntheta+1)
        # [rad] Angles are internally defined in radians, not degrees
        Thetas = pi/2. - (logspace(log10(ped), log10(pi/2. + ped), ntheta+1) - ped)[end:-1:1]
        #Spaced closer near the z=0
    else
        Thetas = linspace(0, pi, ntheta+1)
        # [rad] Angles are internally defined in radians, not degrees
    end

    Phis = Float64[0.0, 0.0] # [rad] cell walls for inactive coordinate

    #Define the cell centers as the average between walls
    rs = 0.5 * (Rs[1:end-1] + Rs[2:end]) # [cm]
    thetas = 0.5 * (Thetas[1:end-1] + Thetas[2:end])
    phis = Float64[0.0]

    return Grid(nr, ntheta, nphi, ncells, Rs, Thetas, Phis, rs, thetas, phis)

end

# Create a grid object using a logarithmic then linear then logarithmic radial spacing
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
    Rs_in = logspace(log10(r_in), log10(r_linstart), n_in + 1) # [cm]

    # linearly spaced middle grid
    Rs_mid = linspace(r_linstart, r_linend, n_mid + 1) # [cm]

    # logarithmically spaced outer grid
    Rs_out = logspace(log10(r_linend), log10(r_out), n_out + 1) # [cm]

    Rs = cat(1, Rs_in[1:end-1], Rs_mid, Rs_out[2:end])

    if eqmirror
        ped = 0.1
        #Thetas = linspace(0, pi/2., ntheta+1)
        # [rad] Angles are internally defined in radians, not degrees
        Thetas = pi/2. - (logspace(log10(ped), log10(pi/2. + ped), ntheta+1) - ped)[end:-1:1]
        #Spaced closer near the z=0
    else
        Thetas = linspace(0, pi, ntheta+1)
        # [rad] Angles are internally defined in radians, not degrees
    end

    Phis = Float64[0.0, 0.0] # [rad] cell walls for inactive coordinate

    #Define the cell centers as the average between walls
    rs = 0.5 * (Rs[1:end-1] + Rs[2:end]) # [cm]
    thetas = 0.5 * (Thetas[1:end-1] + Thetas[2:end])
    phis = Float64[0.0]

    return Grid(nr, ntheta, nphi, ncells, Rs, Thetas, Phis, rs, thetas, phis)
end

# Read from a dictionary, then choose how to make the grid based upon the arguments
function Grid(d::Dict)
    if "r_linstart" in keys(d)
        # We're going for a log-linear-log grid
        names = ["r_in", "r_linstart", "r_linend", "r_out", "n_in", "n_mid", "n_out", "ntheta"]
    else
        # We're just going for a log spaced grid
        names = ["nr", "ntheta", "r_in", "r_out"]
    end
    vec = [d[name] for name in names]
    # Unroll these into an actual parameter
    return Grid(vec...)
end

#This function only needs to be run once, upon setup.
function write_grid(basedir::AbstractString, grid::Grid)
    #amr_grid.inp
    f = open(basedir * "amr_grid.inp", "w")

    #Write the header
    @printf(f, "%d\n", 1) #iformat
    @printf(f, "%d\n", 0) #regular grid (no AMR or Oct-tree)
    @printf(f, "%d\n", 100) #spherical coordiantes
    @printf(f, "%d\n", 0) #gridinfo (none needed for now)
    #incl_r incl_phi incl_z #use this axis?
    @printf(f, "%d %d %d \n", 1, 1, 0) # 2D axisymmetric
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

abstract AbstractParameters

type ParametersStandard <: AbstractParameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] characteristic radius
    T_10::Float64 # [K] temperature at 10 AU
    q::Float64 # temperature gradient exponent
    gamma::Float64 # surface temperature gradient exponent
    Sigma_c::Float64 # [g/cm^2] surface density at characteristic radius
    ksi::Float64 # [cm s^{-1}] microturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
end

type ParametersTruncated <: AbstractParameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] Characteristic radius
    T_10::Float64 # [K] temperature at 10 AU
    q::Float64 # temperature gradient exponent
    gamma::Float64 # surface temperature gradient exponent
    gamma_e::Float64 # exponent for outer exponential surface density taper
    Sigma_c::Float64 # [g/cm^2] surface density at characteristic radius
    ksi::Float64 # [cm s^{-1}] microturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
end

type ParametersCavity <: AbstractParameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] characteristic radius
    r_cav::Float64 # [AU] inner radius of the disk, where an exponentially depleted cavity starts
    T_10::Float64 # [K] temperature at 10 AU
    q::Float64 # temperature gradient exponent
    gamma::Float64 # surface density gradient exponent
    gamma_cav::Float64 # surface density gradient exponent for inner cavity
    Sigma_c::Float64 # [g/cm^2] surface density at characteristic radius
    ksi::Float64 # [cm s^{-1}] microturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
end

type ParametersVertical <: AbstractParameters
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
    Sigma_c::Float64 # [g/cm^2] surface density at characteristic radius
    ksi::Float64 # [cm s^{-1}] micsroturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
end

type ParametersVerticalEta <: AbstractParameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] characteristic radius
    T_10m::Float64 # [K] temperature at 10 AU, midplane
    q_m::Float64 # midplane temperature gradient exponent
    T_freeze::Float64 # [K] temperature below which to reduce CO abundance
    X_freeze::Float64 # [ratio] amount to reduce CO abundance
    sigma_s::Float64 # Photodissociation boundary in units of A_V.
    gamma::Float64 # surface temperature gradient exponent
    h::Float64 # Number of scale heights that z_q is at, typically fixed to 4
    eta::Float64 # Exponent for zq profile
    delta::Float64 # Shape exponent, currently fixed to 2
    Sigma_c::Float64 # [g/cm^2] surface density at characteristic radius
    ksi::Float64 # [cm s^{-1}] micsroturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
end

"""A dictionary of parameter lists for conversion."""
registered_params = Dict([("standard", ["M_star", "r_c", "T_10", "q", "gamma", "Sigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]),
("truncated", ["M_star", "r_c", "T_10", "q", "gamma", "gamma_e", "Sigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]),
("cavity", ["M_star", "r_c", "r_cav", "T_10", "q", "gamma", "gamma_cav", "Sigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]),
("vertical", ["M_star", "r_c", "T_10m", "q_m", "T_10a", "q_a", "T_freeze", "X_freeze", "sigma_s", "gamma", "h", "delta", "Sigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]),
("verticalEta", ["M_star", "r_c", "T_10m", "q_m", "T_freeze", "X_freeze", "sigma_s", "gamma", "h", "eta", "delta", "Sigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"])])

registered_types = Dict([("standard", ParametersStandard), ("truncated", ParametersTruncated), ("cavity", ParametersCavity), ("vertical", ParametersVertical), ("verticalEta", ParametersVerticalEta)])

"""Unroll a vector of parameter values into a parameter type."""
function convert_vector(p::Vector{Float64}, model::AbstractString, fix_params::Vector; args...)
    args = Dict{Symbol}{Float64}(args)

    # The goal is to assemble a length-nparam vector that can be unrolled into a parameter type
    # e.g., ParametersStandard(M_star, r_c, T_10, q, gamma, Sigma_c, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

    # Select the registerd parameters corresponding to this model
    # These are the names listed in this file (model.jl)
    reg_params = registered_params[model]

    # fix_params is a list of strings from a config file that list the names of parameters to be fixed to the config
    # value. These names correspond to the registered names.

    # fit_params are the ones in reg_params that are not in (∉) fix_params
    fit_params = filter(x->∉(x,fix_params), reg_params)

    nparams = length(reg_params)

    # Make an empty vector of this same length
    par_vec = Array(Float64, nparams)

    # This requires assigning p to fit_params
    # Find the indexes that correspond to fit_params
    par_indexes = findin(reg_params, fit_params)
    # Stuff p directly into these
    par_vec[par_indexes] = p

    # Then reading fix_params from args
    # First, create an array of fixed values analogous to p
    p_fixed = Float64[args[convert(Symbol, par)] for par in fix_params]
    par_indexes = findin(reg_params, fix_params)
    par_vec[par_indexes] = p_fixed

    # Now that we are sampling for log10M_gas for the verticalEta model, this part gets tricky.

    # Find the location of logSigma_c and make it Sigma_c
    # Even if we are using verticalEta, this will still be in here because it is in reg_params
    # Only though it currently corresponds to log10M_gas instead of log10Sigma_c
    indSigma_c = findin(reg_params, ["Sigma_c"])
    @assert length(indSigma_c) == 1 "Could not find Sigma_c in order to convert from logSigma_c or logM_gas."

    if model == "verticalEta"
      # Convert from log10M_gas to Sigma_c

      M_gas = 10.^par_vec[indSigma_c] * M_sun # [g]

      # Find gamma and r_c
      r_c = par_vec[2] * AU # [cm]
      gamma = par_vec[8]

      Sigma_c = M_gas * (2 - gamma) / (2 * pi * r_c^2)
      par_vec[indSigma_c] = Sigma_c
    elseif model == "standard"
      # Convert from log10M_gas to Sigma_c
      M_gas = 10.^par_vec[indSigma_c] * M_sun # [g]

      # Find gamma and r_c
      r_c = par_vec[2] * AU # [cm]
      gamma = par_vec[5]

      Sigma_c = M_gas * (2 - gamma) / (2 * pi * r_c^2)
      par_vec[indSigma_c] = Sigma_c

    else
      par_vec[indSigma_c] = 10.^par_vec[indSigma_c]
    end
    # Then assembling these in the same orignial order as registered_params, into the parameter
    # type corresponding to the model.
    return registered_types[model](par_vec...)

end

"""Used to turn a dictionary of parameter values (from config.yaml) directly into a parameter type. Generally used for synthesis and plotting command line scripts."""
function convert_dict(p::Dict, model::AbstractString)
    # Select the registerd parameters corresponding to this model
    reg_params = registered_params[model]
    nparams = length(reg_params)

    if model == "verticalEta"

      M_gas = 10.^p["logM_gas"] * M_sun # [g]

      # Find gamma and r_c
      r_c = p["r_c"] * AU # [cm]
      gamma = p["gamma"]

      Sigma_c = M_gas * (2 - gamma) / (2 * pi * r_c^2)
      p["Sigma_c"] = Sigma_c

    else
      # add a new field, which is the conversion of logSigma_c to Sigma_c
      p["Sigma_c"] = 10^p["logSigma_c"]
    end

    # Using this order of parameters, unpack the dictionary p into a vector
    # reg_params reads Sigma_c, not logSigma_c
    par_vec = Float64[p[par_name] for par_name in reg_params]

    return registered_types[model](par_vec...)

end


"""The common sense priors that apply to all parameter values"""
function lnprior_base(pars::AbstractParameters, dpc_mu::Float64, dpc_sig::Float64)
    # Create a giant short-circuit or loop to test for sensical parameter values.
    if pars.M_star <= 0.0 || pars.ksi <= 0. || pars.T_10 <= 0. || pars.r_c <= 0.0  || pars.T_10 > 1500. || pars.q < 0. || pars.q > 1.0 || pars.incl < 0. || pars.incl > 180. || pars.PA < -180. || pars.PA > 520.
        # println("M_star ", pars.M_star)
        # println("r_c ", pars.r_c)
        # println("T_10 ", pars.T_10)
        # println("q ", pars.q)
        # println("incl ", pars.incl)
        # println("PA ", pars.PA)
        throw(ModelException("Parameters outside of prior range."))
    end

    # Impose distance prior
    dlow = dpc_mu - 3. * dpc_sig
    dhigh = dpc_mu + 3. * dpc_sig

    # hard +/- 3 sigma cutoff
    if (pars.dpc < dlow) || (pars.dpc > dhigh)
        throw(ModelException("Distance outside of 3 sigma prior range."))
    end

    # If we've passed all the hard-cut offs by this point, return the sum of the distance prior and the geometrical inclination prior.
    return -0.5 * (pars.dpc - dpc_mu)^2 / dpc_sig^2 + log(0.5 * sind(pars.incl))

end

function lnprior(pars::ParametersStandard, dpc_mu::Float64, dpc_sig::Float64, grid::Grid)

    lnp = lnprior_base(pars, dpc_mu, dpc_sig)

    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.
    if (3 * pars.r_c) > r_out
        throw(ModelException("Model radius too large for grid size."))
    else
        return lnp
    end

end

function lnprior(pars::ParametersTruncated, dpc_mu::Float64, dpc_sig::Float64, grid::Grid)
    lnp = lnprior_base(pars, dpc_mu, dpc_sig)

    r_out = grid.Rs[end]/AU # [AU]

    if (3 * pars.r_c) > r_out || pars.gamma_e > 2.0
        throw(ModelException("Disk radii outside model grid radius, or outer power law too flat."))
    else
        return lnp
    end

end

function lnprior(pars::ParametersCavity, dpc_mu::Float64, dpc_sig::Float64, grid::Grid)
    lnp = lnprior_base(pars, dpc_mu, dpc_sig)

    r_in = grid.Rs[1]/AU # [AU]
    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.

    # Also check to make sure that r_cav is less than r_c but larger than r_in.
    if (3 * pars.r_c) > r_out || pars.r_cav < r_in || pars.r_cav > pars.r_c || pars.gamma_cav < 0.0
        throw(ModelException("Model radius too large, grid size or cavity too large, or gamma_cav less than zero."))
    else
        return lnp
    end

end

function lnprior(pars::ParametersVertical, dpc_mu::Float64, dpc_sig::Float64, grid::Grid)
    # Create a giant short-circuit or loop to test for sensical parameter values.
    if pars.M_star <= 0.0 || pars.ksi <= 0. || pars.T_10a <= 0. || pars.T_10m <= 0. || pars.r_c <= 0.0  || pars.T_10a > 1500. || pars.q_m < 0. || pars.q_a < 0. || pars.q_m > 1.0 || pars.q_a > 1.0 || pars.incl < 0. || pars.incl > 180. || pars.PA < -180. || pars.PA > 520. || pars.X_freeze > 1.0 || pars.sigma_s < 0.0
        throw(ModelException("Parameters outside of prior range."))
    end

    # Check to see that the temperatures make sense
    if pars.T_10m > pars.T_10a
        throw(ModelException("Atmosphere is cooler than midplane."))
    end

    # Impose distance prior
    dlow = dpc_mu - 3. * dpc_sig
    dhigh = dpc_mu + 3. * dpc_sig

    # hard +/- 3 sigma cutoff
    if (pars.dpc < dlow) || (pars.dpc > dhigh)
        throw(ModelException("Distance outside of 3 sigma prior range."))
    end

    # If we've passed all the hard-cut offs by this point, return the sum of the distance prior and the geometrical inclination prior.
    lnp = -0.5 * (pars.dpc - dpc_mu)^2 / dpc_sig^2 + log(0.5 * sind(pars.incl))

    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.
    if (3 * pars.r_c) > r_out
        throw(ModelException("Model radius too large for grid size."))
    else
        return lnp
    end

end

function lnprior(pars::ParametersVerticalEta, dpc_mu::Float64, dpc_sig::Float64, grid::Grid)
    # Create a giant short-circuit or loop to test for sensical parameter values.
    if pars.M_star <= 0.0 || pars.ksi <= 0. || pars.T_10m <= 60. || pars.r_c <= 0.0  || pars.q_m < 0. || pars.q_m > 1.0 || pars.incl < 0. || pars.incl > 180. || pars.PA < -180. || pars.PA > 520. || pars.X_freeze > 1.0 || pars.sigma_s < 0.0 || pars.eta < 0.2 || pars.eta > 6.0  || pars.h < 0.5 || pars.h > 6.0 || pars.delta < 0.5 || pars.delta > 6.0 || pars.gamma < 0.5 || pars.gamma > 3.0
        throw(ModelException("Parameters outside of prior range."))
    end

    # Impose distance prior
    dlow = dpc_mu - 3. * dpc_sig
    dhigh = dpc_mu + 3. * dpc_sig

    # hard +/- 3 sigma cutoff
    if (pars.dpc < dlow) || (pars.dpc > dhigh)
        throw(ModelException("Distance outside of 3 sigma prior range."))
    end

    # If we've passed all the hard-cut offs by this point, return the sum of the distance prior and the geometrical inclination prior.
    lnp = -0.5 * (pars.dpc - dpc_mu)^2 / dpc_sig^2 + log(0.5 * sind(pars.incl))

    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.
    if (3 * pars.r_c) > r_out
        throw(ModelException("Model radius too large for grid size."))
    else
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
    sizeau_command = sizeau_desired/(1. + RADMC_SIZEAU_SHIFT)

    return (sizeau_desired, sizeau_command) # [AU]

end

# Assume all inputs to these functions are in CGS units and in *cylindrical* coordinates.
# Parametric type T allows passing individual Float64 or Vectors.
# # Alternate functions accept pars passed around, where pars is in M_star, AU, etc...
function velocity{T}(r::T, M_star::Float64)
    sqrt(G * M_star ./ r)
end
velocity{T}(r::T, pars::AbstractParameters) = velocity(r, pars.M_star * M_sun)

# For the vertical temperature gradient
function velocity(r::Float64, z::Float64, M_star::Float64)
    sqrt(G * M_star / (r^2 + z^2)^(3./2)) * r
end
velocity(r::Float64, z::Float64, pars::ParametersVertical) = velocity(r, z, pars.M_star * M_sun)
velocity(r::Float64, z::Float64, pars::ParametersVerticalEta) = velocity(r, z, pars.M_star * M_sun)

function temperature{T}(r::T, T_10::Float64, q::Float64)
    T_10 * (r ./ (10. * AU)).^(-q)
end
temperature{T}(r::T, pars::AbstractParameters) = temperature(r, pars.T_10, pars.q)

function Hp{T}(r::T, M_star::Float64, T_10::Float64, q::Float64)
    temp = temperature(r, T_10, q)
    sqrt(kB * temp .* r.^3./(mu_gas * m_H * G * M_star))
end
Hp{T}(r::T,  pars::AbstractParameters) = Hp(r, pars.M_star * M_sun, pars.T_10, pars.q)
# Scale height computed from midplane temperature for vertical temperature gradient model
Hp{T}(r::T,  pars::Union{ParametersVertical,ParametersVerticalEta}) = Hp(r, pars.M_star * M_sun, pars.T_10m, pars.q_m)

# For the vertical temperature gradient model
T_mid{T}(r::T, pars::ParametersVertical) = temperature(r, pars.T_10m, pars.q_m)
T_atm{T}(r::T, pars::ParametersVertical) = temperature(r, pars.T_10a, pars.q_a)

T_mid{T}(r::T, pars::ParametersVerticalEta) = temperature(r, pars.T_10m, pars.q_m)

function T_atm{T}(r::T, pars::ParametersVerticalEta)
    Tm = T_mid(r, pars)

    # the atmosphere temperature is the maximum of 500 or twice the midplane
    # in order to get the structure close to the star correct
    Ta = maximum([500, 2 * Tm])
    return Ta
end


# Atmosphere height computed as multiple of scale height (computed at midplane)
function z_q{T}(r::T, pars::ParametersVertical)
    return pars.h * Hp(r, pars)
end

function z_q{T}(r::T, pars::ParametersVerticalEta)

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


# Calculate the gas surface density
function Sigma(r::Float64, pars::Union{ParametersStandard, ParametersVertical, ParametersVerticalEta})
    r_c = pars.r_c * AU

    gamma = pars.gamma
    Sigma_c = pars.Sigma_c

    S = Sigma_c * (r/r_c)^(-gamma) * exp(-(r/r_c)^(2 - gamma))

    return S
end

function Sigma(r::Float64, pars::ParametersTruncated)
    r_c = pars.r_c * AU # [cm]
    gamma = pars.gamma
    gamma_e = pars.gamma_e
    Sigma_c = pars.Sigma_c

    S = Sigma_c * (r/r_c)^(-gamma) * exp(-(r/r_c)^(2 - gamma_e))
    return S
end

function Sigma(r::Float64, pars::ParametersCavity)
    r_c = pars.r_c * AU
    r_cav = pars.r_cav * AU

    gamma = pars.gamma
    gamma_cav = pars.gamma_cav
    Sigma_c = pars.Sigma_c

    inner_taper = exp(-(r_cav/r)^gamma_cav)
    outer_taper = exp(-(r/r_c)^(2 - gamma))
    power_law = (r/r_c)^(-gamma)

    S = Sigma_c * inner_taper * power_law * outer_taper

    return S
end

# Delivers a gas density in g/cm^3
function rho_gas(r::Float64, z::Float64, pars::AbstractParameters)
    H = Hp(r, pars)
    S = Sigma(r, pars)

    # Calculate the density
    rho = S/(sqrt(2. * pi) * H) * exp(-0.5 * (z/H)^2)

    return rho
end

# Determines the midplane density assuming vertically isothermal
function rho_gas_mid(r::Float64, pars::AbstractParameters)
    H = Hp(r, pars)
    S = Sigma(r, pars)

    # Calculate the density at the midplane
    rho = S/(sqrt(2. * pi) * H)

    return rho

end

# The temperature change, dT/dz
function dT(r::Float64, z::Float64, pars::ParametersVertical)
    zq = z_q(r, pars)
    Ta = T_atm(r, pars)

    if z >= zq
        return 0.
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
        return 0.
    elseif z < H
      return 0.
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


function rho_gas(r::Float64, z::Float64, pars::Union{ParametersVertical, ParametersVerticalEta})

    # Calculate the "top" of the atmosphere,
    # a height where we are sure we can enforce that density = 0
    ztop = z_top(r, pars)

    # If we are querying a height above this, just return 0 now
    if z > ztop
        return constants.rho_gas_zero
    end

    # Calculate what the maximum surface density would be if we integrated all the way to  the midplane
    sigma = Sigma(r, pars) / 2

    # Calculate the photodissociation height
    # Threshold column density
    thresh_H2 = pars.sigma_s * Av_sigmaH / 2 # [n_H2/cm^2]
    thresh = thresh_H2 * (mu_gas * amu / X_H2) # [g/cm^2] of gas

    # If there is not even enough column density at the midplane in order to exceed the threshold
    # to protect against photodissociation of CO, just return 0 now
    # if sigma < thresh
    #     return constants.rho_gas_zero
    # end

    # Create an array of heights that goes from the midplane (0, 0.001, ..., ztop)
    zs = cat(1, [0], logspace(log10(0.001 * AU), log10(ztop), 64 - 1))

    # Define the ODE at this radius r
    function f(z, y)
        # calculate y'
        dlnrho(r, z, pars) * y
    end

    # Choose a starting guess corresponding to the midplane density for an isothermal model.
    # This helps us get in the ballpark and avoid a lot of the numerical errors that result
    start = rho_gas_mid(r, pars)

    # Now solve the ODE on the grid of postulated values
    # For the inner disk, we require a larger value of abstol
    # We need to experiment with what is best
    abstol=1e-22
    z_out, y_out = ode45(f, start, zs; abstol=1e-22)

    # where y_out is less than the minimum accuracy we can trust, just set it to the minimum value
    y_out[y_out .< abstol] = abstol

    # Integrate the output from the ODE solver to find the total (un-normalized) surface density
    spl = Spline1D(z_out, y_out)
    tot = integrate(spl, z_out[1], z_out[end])

    # Apply a correction factor to the un-normalized densities
    cor = sigma / tot
    rhos = y_out .* cor

    # This is the function that will be minimized
    function g(zvar::Float64)
        return abs(thresh - cor * integrate(spl, zvar, z_out[end]))
    end

    z_phot = optimize(g, z_out[1], z_out[end]).minimum

    # Evaluate rho at this point
    rho = cor * evaluate(spl, z)

    # If we are above this height, return the gas density reduced by a factor of 100
    if z > z_phot
        rho = 1e-2 * rho
    end

    if rho < 0.0
      println("r ", r/AU, " z ", z/AU, " z_top ", ztop/AU, " z_phot ", z_phot/AU, " y_out ", y_out)
      return (start, z_out, y_out)
      throw(ModelException("Rho less than 0.0"))
    end

    return rho
end

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
