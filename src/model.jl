module model

export write_grid, write_model, write_lambda, write_dust, Grid
export AbstractParameters, ParametersStandard, ParametersTruncated, ParametersCavity, ParametersVertical, convert_vector, convert_dict
export lnprior

# The double dot is because we are now inside the model module, and we want to import the
# constants module, which is part of the enclosing DiskJockey package.
using ..constants
using Dierckx

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

function Grid(nr::Int, ntheta::Int, r_in::Real, r_out::Real, eqmirror::Bool)
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
    M_gas::Float64 # [M_Sun] disk mass of gas
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
    r_in::Float64 # [AU] hard inner radius of the disk
    r_out::Float64 # [AU] radius at which density is depleted by delta
    T_10::Float64 # [K] temperature at 10 AU
    q::Float64 # temperature gradient exponent
    gamma::Float64 # surface temperature gradient exponent
    M_gas::Float64 # [M_Sun] disk mass of gas
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
    gamma::Float64 # surface temperature gradient exponent
    M_gas::Float64 # [M_Sun] disk mass of gas
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
    M_gas::Float64 # [M_Sun] disk mass of gas
    ksi::Float64 # [cm s^{-1}] micsroturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
end

"""Used to unroll a vector of parameter values into a parameter type."""
function convert_vector(p::Vector{Float64}, model::AbstractString, fix_d::Bool; args...)
    args = Dict{Symbol}{Float64}(args)
    if model == "standard"
        # Take gamma from the args
        gamma = args[:gamma]
        if fix_d
            # distance is fixed and provided by the args
            dpc = args[:dpc]
            M_star, r_c, T_10, q, logM_gas, ksi, incl, PA, vel, mu_RA, mu_DEC = p
        else
            # we are fitting with distance as a parameter
            M_star, r_c, T_10, q, logM_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC = p
        end

        M_gas = 10^logM_gas
        return ParametersStandard(M_star, r_c, T_10, q, gamma, M_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

    elseif model == "truncated"
        # Take gamma from the args
        gamma = args[:gamma]
        if fix_d
            dpc = args[:dpc]
            M_star, r_in, r_out, T_10, q, logM_gas, ksi, incl, PA, vel, mu_RA, mu_DEC = p
        else
            M_star, r_in, r_out, T_10, q, logM_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC = p
        end

        M_gas = 10^logM_gas
        return ParametersTruncated(M_star, r_in, r_out, T_10, q, gamma, M_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

    elseif model == "cavity"
        gamma = args[:gamma]
        if fix_d
            dpc = args[:dpc]
            M_star, r_c, r_cav, T_10, q, logM_gas, ksi, incl, PA, vel, mu_RA, mu_DEC = p
        else
            M_star, r_c, r_cav, T_10, q, logM_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC = p
        end

        M_gas = 10^logM_gas
        return ParametersCavity(M_star, r_c, r_cav, T_10, q, gamma, M_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

    elseif model == "vertical"
        gamma = args[:gamma]
        T_freeze = args[:T_freeze]
        X_freeze = args[:X_freeze]
        sigma_s = args[:sigma_s]
        h = args[:h]
        delta = args[:delta]
        if fix_d
            dpc = args[:dpc]
            M_star, r_c, T_10m, q_m, T_10a, q_a, logM_gas, ksi, incl, PA, vel, mu_RA, mu_DEC = p
        else
            M_star, r_c, T_10m, q_m, T_10a, q_a, logM_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC = p
        end

        M_gas = 10^logM_gas
        return ParametersVertical(M_star, r_c, T_10m, q_m, T_10a, q_a, T_freeze, X_freeze, sigma_s, gamma, h, delta, M_gas, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)
    else
        # Raise an error that we don't know this model.
        throw(ErrorException("Model type $model not yet implemented in model.jl"))
    end
end

"""Used to turn a dictionary of parameter values (from config.yaml) into a parameter type."""
function convert_dict(p::Dict, model::AbstractString)
    if model == "standard"
        names = ["M_star", "r_c", "T_10", "q", "gamma", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]

        vec = Float64[p[name] for name in names]
        vec[6] = 10^vec[6]# 10^ gas

        # Unroll these into an actual parameter
        return ParametersStandard(vec...)

    elseif model == "truncated"
        names = ["M_star", "r_in", "r_out", "T_10", "q", "gamma", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]

        vec = Float64[p[name] for name in names]
        vec[7] = 10^vec[7]# 10^ gas

        # Unroll these into an actual parameter
        return ParametersTruncated(vec...)
    elseif model == "cavity"
        names = ["M_star", "r_c", "r_cav", "T_10", "q", "gamma", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]

        vec = Float64[p[name] for name in names]
        vec[7] = 10^vec[7]# 10^ gas

        # Unroll these into an actual parameter
        return ParametersCavity(vec...)
    elseif model == "vertical"
        names = ["M_star", "r_c", "T_10m", "q_m", "T_10a", "q_a", "T_freeze", "X_freeze", "sigma_s", "gamma", "h", "delta", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]

        vec = Float64[p[name] for name in names]
        vec[13] = 10^vec[13]# 10^ gas

        # Unroll these into an actual parameter
        return ParametersVertical(vec...)

    else
        # Raise an error that we don't know this model.
        throw(ErrorException("Model type $model not yet implemented in model.jl"))
    end
end

"""The common sense priors that apply to all parameter values"""
function lnprior_base(pars::AbstractParameters, dpc_mu::Float64, dpc_sig::Float64)
    # Create a giant short-circuit or loop to test for sensical parameter values.
    if pars.M_star <= 0.0 || pars.ksi <= 0. || pars.T_10 <= 0. || pars.r_c <= 0.0  || pars.T_10 > 1500. || pars.q < 0. || pars.q > 1.0 || pars.incl < 0. || pars.incl > 180. || pars.PA < -180. || pars.PA > 520.
        return -Inf
    end

    # Impose distance prior
    dlow = dpc_mu - 3. * dpc_sig
    dhigh = dpc_mu + 3. * dpc_sig

    # hard +/- 3 sigma cutoff
    if (pars.dpc < dlow) || (pars.dpc > dhigh)
        return -Inf
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
        return -Inf
    else
        return lnp
    end

end

function lnprior(pars::ParametersTruncated, dpc_mu::Float64, dpc_sig::Float64, grid::Grid)
    lnp = lnprior_base(pars, dpc_mu, dpc_sig)

    # Hard cutoff on inner and outer edges.
    r_in = grid.Rs[1]/AU # [AU]
    r_out = grid.Rs[end]/AU # [AU]

    if pars.r_in < r_in || pars.r_out > r_out
        return -Inf
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
    if (3 * pars.r_c) > r_out || pars.r_cav < r_in || pars.r_cav > pars.r_c
        return -Inf
    else
        return lnp
    end

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
Hp{T}(r::T,  pars::ParametersVertical) = Hp(r, pars.M_star * M_sun, pars.T_10m, pars.q_m)

# For the vertical temperature gradient model
T_mid{T}(r::T, pars::ParametersVertical) = temperature(r, pars.T_10m, pars.q_m)
T_atm{T}(r::T, pars::ParametersVertical) = temperature(r, pars.T_10a, pars.q_a)

# Atmosphere height computed as multiple of scale height (computed at midplane)
function z_q{T}(r::T, pars::ParametersVertical)
    return pars.h * Hp(r, pars)
end

function temperature(r::Float64, z::Float64, pars::ParametersVertical)
    zq = z_q(r, pars)
    Ta = T_atm(r, pars)
    if z >= zq
        return Ta
    else
        Tm = T_mid(r, pars)
        # return Ta + (Tm - Ta) * cos(pi * z/ (2 * zq))^pars.delta # Dartois
        return Tm + (Ta - Tm) * sin(pi * z/ (2 * zq))^(2 * pars.delta) # Williams and Best 14
    end
end

# Calculate the gas surface density
function Sigma(r::Float64, pars::Union{ParametersStandard, ParametersVertical})
    r_c = pars.r_c * AU

    gamma = pars.gamma
    M_gas = pars.M_gas * M_sun

    Sigma_c = M_gas * (2 - pars.gamma) / (2 * pi * r_c^2)

    S = Sigma_c * (r/r_c)^(-gamma) * exp(-(r/r_c)^(2 - gamma))

    return S
end

function Sigma(r::Float64, pars::ParametersTruncated)
    r_c = 1. * AU # [cm]
    r_in = pars.r_in * AU
    r_out = pars.r_out * AU

    if r > r_in && r < r_out
        Sigma_c = pars.M_gas * M_sun * (2 - pars.gamma) / (2 * pi * r_c^2 * ((r_out/r_c)^(2 - pars.gamma) - (r_in/r_c)^(2 - pars.gamma)))
        return Sigma_c * (r/r_c)^(-pars.gamma)
    else
        return 0.0
    end
end

function Sigma(r::Float64, pars::ParametersCavity)
    r_c = pars.r_c * AU
    r_cav = pars.r_cav * AU

    gamma = pars.gamma
    M_gas = pars.M_gas * M_sun

    Sigma_c = M_gas * (2 - pars.gamma) / (2 * pi * r_c^2)

    inner_taper = exp(-(r_cav/r)^(2 - gamma))
    outer_taper = exp(-(r/r_c)^(2 - gamma))
    power_law = (r/r_c)^(-gamma)

    S = Sigma_c * inner_taper * power_law * outer_taper

    return S
end

# Delivers a gas density in g/cm^3
function rho_gas(r::Float64, z::Float64, pars::AbstractParameters)
    H = Hp(r, pars)
    S = Sigma(r, pars)

    rho = S/(sqrt(2. * pi) * H) * exp(-0.5 * (z/H)^2)
    if rho < rho_gas_critical
        return 0.0
    else
        return rho
    end
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

# The integrand
function dlnrho(r::Float64, z::Float64, pars::ParametersVertical)
    T = temperature(r, z, pars)
    dTemp = dT(r, z, pars)
    return -(dTemp/T + (mu_gas * m_H)/(kB * T) * (G * pars.M_star * M_sun * z)/(r^2 + z^2)^1.5)
end

# The integrand for a vertically isothermal testcase
function dlnrho(r::Float64, z::Float64, pars::ParametersStandard)
    return -(mu_gas * m_H)/(kB * T) * (G * pars.M_star * M_sun * z)/r^3
end

# Delivers an unnormalized gas density. Needs to be multiplied by correction factor.
function un_lnrho(r::Float64, z::Float64, ztop::Float64, pars::ParametersVertical)
    f(x) = dlnrho(r, x, pars)

    val, err = quadgk(f, ztop, z)
    return val
end

# Calculate a vertical 1-D density column at a given radius
function rho_column(r::Float64, pars::ParametersVertical)

    zq = z_q(r, pars) # The height of the disk "atmosphere"

    # 10 AU is an arbitrary number, but this is to make sure that we have a high enough bound
    # in the inner disk
    # It is assumed that there is 0 gas above this upper bound, ztop.
    if (5 * zq) >= 10 * AU
        ztop = zq * 5
    else
        ztop = 10 * AU
    end

    # Create an array of heights that goes from the midplane (0, 0.01, ..., ztop)
    zs = cat(1, [0], logspace(log10(0.01 * AU), log10(ztop), n_z_interpolator-1))

    un_lnrhos = Array(Float64, n_z_interpolator)
    # For each height, numerically integrate dlnrho/dz to yeild a yet-to-be normalized un_lnrho
    for i=1:n_z_interpolator
        un_lnrhos[i] = un_lnrho(r, zs[i], ztop, pars)
    end

    # Find the largest value from this array and subtract it to avoid numberical overflow
    # in subsequent steps
    bar_lnrho = un_lnrhos[1]
    un_lnrhos -= bar_lnrho

    # Convert from ln(rho) to rho; still unnormalized.
    un_rhos = exp(un_lnrhos)

    # Now that we have the shape of the density structure but not the normalization,
    # integrate to find the normalization factor.
    # The density is a unitful quantity, so we need to exponentiate (can't do this with un_lnrhos).
    spl = Spline1D(zs, un_rhos)

    # Integral is from -inf to inf, or 2 x (0 to inf).
    # This is the total unnormalized surface density.
    tot = 2 * integrate(spl, 0., zs[end])

    # The actual value of the surface density
    S = Sigma(r, pars)

    # The normalized density columns
    rhos = S/tot * un_rhos # [g/cm^3]

    return (zs, rhos)
end

# Translate from a gas density to a CO number density
function rho_column_CO(r::Float64, zs::Vector{Float64}, rhos::Vector{Float64}, pars::ParametersVertical)

    # Threshold column density
    thresh_H2 = pars.sigma_s * Av_sigmaH # [n_H2/cm^2]
    thresh = thresh_H2 * (mu_gas * amu / X_H2) # [g/cm^2] of gas

    # Normalize the previous spline and integrate to get zboundary
    # f(x) = S/tot * integrate(spl, x, zs[end])

    # Go through each of the z spacings, and use Riemann integration in a while loop to find
    # at which z point we've finally accumulated the threshold column density.
    column = 0.0 # total acumulated column
    icolumn::Int = n_z_interpolator # index of the z point
    while (column <= thresh) && (icolumn > 1)
        dz = zs[icolumn] - zs[icolumn - 1]
        column += dz *  rhos[icolumn]
        icolumn -= 1
    end

    # All densities below this z height remain intact, while all densities above this
    # height are set to (effectively) zero (CO is dissociated)
    rhos[icolumn:end] = constants.rho_gas_zero

    # Implement CO freezout onto grains
    for i=1:n_z_interpolator
        rhos[i] *= X_freeze(temperature(r, zs[i], pars), pars)
    end

    return (zs, rhos)
end


#
# # Calculate the multiplicative factor needed to normalize `un_rho` based upon the boundary
# # condition that the integral of rho(r, z) with z from -inf to inf should be equal to Sigma(r)
# function lncorrection_factor(r::Float64, pars::Parameters)
#     S = Sigma(r, pars)
#     zq = z_q(r, pars)
#     # Most of the total density is near the midplane. This is necessary to create a logspaced array with 0 as the beginning.
#
#     nz = 50
#
#     # How best to choose the outer bound? Some multiple of zq is probably a good idea, but this
#     # doesn't really go to high enough altitude in the inner disk. Instead, let's choose a minimum threshold.
#     if (10 * zq) >= 10 * AU
#         outer_reaches = zq * 10
#     else
#         outer_reaches = 10 * AU
#     end
#
#     zints = cat(1, [0], logspace(log10(0.01 * AU), log10(outer_reaches), nz-1))
#
#     un_lnrhos = Array(Float64, nz)
#     for i=1:nz
#         un_lnrhos[i] = un_lnrho(r, zints[i], pars)
#     end
#
#     bar_lnrho = un_lnrhos[1]
#
#     # Remove the largest value so we don't get an overflow error
#     un_lnrhos -= bar_lnrho
#
#     # Convert from ln(rho) to rho, but this is still unnormalized.
#     un_rhos = exp(un_lnrhos)
#
#     # Now that we have the shape of the density structure but not the normalization,
#     # integrate to get the normalization
#     # The density should be unitful, so we need to exponentiate.
#     spl = Spline1D(zints, un_rhos)
#
#     # Integral is from -inf to inf
#     tot = 2 * integrate(spl, 0., zints[end])
#
#     # println("Norm. tot: ", tot, " lntot:", log(tot))
#     # println("bar_lnrho. exp: ", exp(bar_lnrho), " bar_lnrho: ", bar_lnrho)
#     # println("With bar_lnrho :", tot * exp(bar_lnrho), " lntot :", log(tot) + bar_lnrho)
#
#     # The correction factor is the ratio of the actual surface density to the value of the integral of our unnormalized density.
#     lncor = log(S) - (log(tot) + bar_lnrho)
#
#     # cor = S/tot
#     # apply this as
#     # rhos = cor * unnormed_rhos
#     return lncor
# end
#
# # Now, create an interpolator for the the correction function. In spherical coordinates, that way, first evaluate it at every point (r, z=0) in the grid. Then, when we want to query the density of the disk at (r, z), which will correspond to a different r_cyl than the one used to compute the normalization, we can just use the spline interpolation and not have to redo the time consuming normalization integral.
# function make_correction_interpolator(pars::Parameters, grid::Grid)
#
#     nr = 64
#     rs = logspace(log10(0.5 * AU), log10(grid.Rs[end]), nr)
#
#     cors = Array(Float64, nr)
#
#     for i=1:nr
#         cors[i] = lncorrection_factor(rs[i], pars)
#     end
#
#     # An interpolating spline for this correction factor.
#     spl = Spline1D(rs, cors)
#
#     return spl
# end

# Because the calculation for the vertical temperature gradient is a bit more complex than other
# models, we can save time by caching the evaluation in a grid interpolator.
# This function returns an interpolator function that can be used to determine rho_gas(r, z)
function rho_gas_interpolator(pars::ParametersVertical, grid::Grid)

    # Create a 2D spline interpolator for this surface using a density column at each radius
    rs = Array(Float64, (n_z_interpolator, grid.nr))
    zs = Array(Float64, (n_z_interpolator, grid.nr))
    rhos = Array(Float64, (n_z_interpolator, grid.nr))
    for i=1:grid.nr
        r = grid.rs[i]
        rs[:,i] = r

        ZS, RHOS = rho_column(r, pars)
        zs[:,i], rhos[:,i] = rho_column_CO(r, ZS, RHOS, pars)
    end

    # Reshape all af rs, zs, and rhos into 1D arrays
    rs = reshape(rs, prod(size(rs)))
    zs = reshape(zs, prod(size(zs)))
    rhos = reshape(rhos, prod(size(rhos)))

    log10rhos = log10(rhos)

    return

    println("interpolator rho_gas extrema ", extrema(rhos))
    spl = Spline2D(rs, zs, log10(rhos)) #; kx=1, ky=1, s=length(rs))

    function interpolator(r::Float64, z::Float64)
        # Santity checks to see that we're inside the grid
        if r > grid.rs[end] || z > grid.rs[end]
            return constants.rho_gas_zero
        else
            return 10^evaluate(spl, r, z)
        end
    end

    return interpolator
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
    nCO = n_CO(r, z, pars) # number of CO molecules per cm^3

    # Convert from nCO to nH2
    nH2 = nCO / 7.e-5 # number density ratio

    # Convert from nH2 (assuming nH2 ~ nGas ) to mGas
    mGas = mu_gas * amu * nH2 # [g]

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

# Because we need to do the density caching for speed purposes, create a separate write_model
function write_model(pars::ParametersVertical, basedir::AbstractString, grid::Grid, species::AbstractString)

    # Prime the density interpolator
    n_CO = number_densities[species] * rho_gas_interpolator(pars, grid)

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

                @printf(fdens, "%.9e\n", n_CO(r_cyl, z))
                @printf(fvel, "0 0 %.9e\n", velocity(r_cyl, z, pars))
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
