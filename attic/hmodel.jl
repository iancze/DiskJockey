module hmodel

# Now for hydrostatic equilibrium, vertical temperature gradient.

export write_grid, write_model, write_lambda, write_dust, Parameters, Grid

# The double dot is because we are now inside the model module, and we want to import the
# constants module, which is part of the enclosing DiskJockey package.
using Dierckx
using ..constants

# Write the wavelength sampling file. Only run on setup
function write_lambda(lams::Array{Float64, 1}, basedir::AbstractString)
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
    Rs::Vector{Float64}
    Thetas::Vector{Float64}
    Phis::Vector{Float64}
    # cell centers
    rs::Vector{Float64}
    thetas::Vector{Float64}
    phis::Vector{Float64}
end

function Grid(nr::Int, ntheta::Int, r_in::Real, r_out::Real, eqmirror::Bool)
    # Specify a 2D axisymmetric *separable* grid in spherical coordinates:
    # {r, theta, phi}, where theta is angle from zenith, phi is azimuth

    # Number of cells in each dimension
    nphi = 1 # axisymmetric disk
    ncells = nr * ntheta * nphi
    r_in = convert(Float64, r_in) * AU # Inner extent of disk
    r_out = convert(Float64, r_out) * AU # Outer extent of disk

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
    rs = 0.5 * (Rs[1:end-1] + Rs[2:end])
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

#Let's try defining a parameters type, an object of which gets passed around.

type Parameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] characteristic radius
    r_in::Float64 # [AU] inner radius of the grid
    r_out::Float64 # [AU] outer radius of the grid
    T_10m::Float64 # [K] temperature at 10 AU, midplane
    q_m::Float64 # midplane temperature gradient exponent
    T_10a::Float64 # [K] temperature at 10 AU, atmosphere
    q_a::Float64 # atmosphere temperature gradient exponent
    T_freeze::Float64 # [K] temperature below which to reduce CO abundance
    X_freeze::Float64 # [ratio] amount to reduce CO abundance
    sigma_s::Float64 # Photodissociation boundary in units of A_V.
    gamma::Float64 # surface temperature gradient exponent
    h::Float64 # Number of scale heights that z_q is at, currently fixed to 4
    delta::Float64 # Shape exponent, currently fixed to 2
    M_gas::Float64 # [M_Sun] disk mass of gas
    delta_gas::Float64 # Fraction by which gas is reduced inside the cavity
    r_cav::Float64 # [AU] the radius interior to which the gas density is reduced by delta_gas
    ksi::Float64 # [cm s^{-1}] micsroturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
end

# Assume all inputs to these functions are in CGS units and in *cylindrical* coordinates.
# Parametric type T allows passing individual Float64 or Vectors.
# Alternate functions accept pars passed around, where pars is in M_star, AU, etc...
function velocity(r::Float64, z::Float64, M_star::Float64)
    sqrt(G * M_star / (r^2 + z^2)^(3./2)) * r
end
velocity(r::Float64, z::Float64, pars::Parameters) = velocity(r, z, pars.M_star * M_sun)

# Midplane temperature
function T_mid{T}(r::T, T_10::Float64, q::Float64)
    T_10 * (r ./ (10. * AU)).^(-q)
end
T_mid{T}(r::T, pars::Parameters) = T_mid(r, pars.T_10m, pars.q_m)

# Atmosphere temperature
function T_atm{T}(r::T, T_10::Float64, q::Float64)
    T_10 * (r ./ (10. * AU)).^(-q)
end
T_atm{T}(r::T, pars::Parameters) = T_atm(r, pars.T_10a, pars.q_a)

# Scale height computed from midplane temperature
function Hp{T}(r::T, M_star::Float64, T_10::Float64, q::Float64)
    temp = T_mid(r, T_10, q)
    return sqrt((kB * temp .* r.^3.)/(mu_gas * m_H * G * M_star))
end
Hp{T}(r::T,  pars::Parameters) = Hp(r, pars.M_star * M_sun, pars.T_10m, pars.q_m)

# Calculate the gas surface density
function Sigma{T}(r::T, pars::Parameters)
    r_c = pars.r_c * AU
    Sigma_c = pars.M_gas * M_sun * (2 - pars.gamma) / (2 * pi * r_c^2)
    Sigma_c .* (r./r_c).^(-pars.gamma) .* exp(-(r./r_c).^(2 - pars.gamma))
end

# Atmosphere height computed as multiple of scale height (computed at midplane)
function z_q{T}(r::T, pars::Parameters)
    return pars.h * Hp(r, pars)
end

# No parametric type for temperature, because it is a 2D function.
# This is according to the Williams and Best 14 equations
function temperature(r::Float64, z::Float64, pars::Parameters)
    zq = z_q(r, pars)
    Ta = T_atm(r, pars)
    if z >= zq
        return Ta
    else
        Tm = T_mid(r, pars)
        # return Ta + (Tm - Ta) * cos(pi * z/ (2 * zq))^pars.delta
        return Tm + (Ta - Tm) * sin(pi * z/ (2 * zq))^(2 * pars.delta)
    end
end

# The temperature change
function dT(r::Float64, z::Float64, pars::Parameters)
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
function dlnrho(r::Float64, z::Float64, pars::Parameters)
    T = temperature(r, z, pars)
    dTemp = dT(r, z, pars)
    return -(dTemp/T + (mu_gas * m_H)/(kB * T) * (G * pars.M_star * M_sun * z)/(r^2 + z^2)^1.5)
    # Vertically isothermal testcase
    # return -(mu_gas * m_H)/(kB * T) * (G * pars.M_star * M_sun * z)/r^3

end

# Delivers an unnormalized gas density. Needs to be multiplied by correction factor.
function un_lnrho(r::Float64, z::Float64, ztop::Float64, pars::Parameters)
    f(x) = dlnrho(r, x, pars)

    val, err = quadgk(f, ztop, z)
    return val
end

# Calculate a slice of density at a given radius
function density_slice(r::Float64, pars::Parameters)
    zq = z_q(r, pars) # The height of the disk "atmosphere"

    # 10 AU is an arbitrary number, but this is to make sure that we have a high enough bound
    # in the inner disk
    if (10 * zq) >= 10 * AU
        ztop = zq * 10
    else
        ztop = 10 * AU
    end

    # It is assumed that there is 0 gas above this upper bound, ztop.

    # Create an array of heights that goes from the midplane (0, 0.01, ..., ztop)
    nz = 64
    zs = cat(1, [0], logspace(log10(0.01 * AU), log10(ztop), nz-1))

    un_lnrhos = Array(Float64, nz)
    # For each height, integrate dlnrho/dz to yeild a yet-to-be normalized un_lnrho
    for i=1:nz
        un_lnrhos[i] = un_lnrho(r, zs[i], ztop, pars)
    end

    # Subtract the largest value from this array to avoid numberical overflow in subsequent steps
    bar_lnrho = un_lnrhos[1]

    # Remove the largest value so we don't get an overflow error
    un_lnrhos -= bar_lnrho

    # Convert from ln(rho) to rho, but this is still unnormalized.
    un_rhos = exp(un_lnrhos)

    # Now that we have the shape of the density structure but not the normalization,
    # integrate to get the normalization.
    # The density should be unitful, so we need to exponentiate (can't do this with un_lnrhos).
    spl = Spline1D(zs, un_rhos)

    # Integral is from -inf to inf
    tot = 2 * integrate(spl, 0., zs[end])

    # The actual value of the surface density
    S = Sigma(r, pars)

    # The normalized density slices
    rhos = S/tot * un_rhos # [g/cm^3]

    # Threshold column
    thresh = pars.sigma_s * Av_sigmaH # [g/cm^2]

    # Normalize the previous spline and integrate to get zboundary
    f(x) = S/tot * integrate(spl, x, zs[end])

    # Go through each of the z spacings, and use Riemann integration in a while loop to find at which z point we've finally
    column = 0.0 # total acumulated column
    icolumn::Int = 0 # index of the z point
    while column <= thresh

    end

    # fsolve, or some other root finding algorithm to find where f(x) = thresh.

    # Make all rhos above this equal to zero.

    return (zs, rhos)
end

function make_density_interpolator(pars::Parameters, grd::Grid)
    # Using the cell centers (in radius), go through and get a density slice.
    # Append the slice to the current list.
    # npoints = 64 * grd.nr
    rs = Array(Float64, (64, grid.nr))
    zs = Array(Float64, (64, grid.nr))
    rhos = Array(Float64, (64, grid.nr))
    for i=1:grid.nr
        r = grid.rs[i]
        rs[:,i] = r
        zs[:,i], rhos[:,i] = density_slice(r, pars)
    end

    # Now flatten these arrays and feed them to a 2D spline with Dierckx



end


# Calculate the multiplicative factor needed to normalize `un_rho` based upon the boundary
# condition that the integral of rho(r, z) with z from -inf to inf should be equal to Sigma(r)
function lncorrection_factor(r::Float64, pars::Parameters)
    S = Sigma(r, pars)
    zq = z_q(r, pars)
    # Most of the total density is near the midplane. This is necessary to create a logspaced array with 0 as the beginning.

    nz = 50

    # How best to choose the outer bound? Some multiple of zq is probably a good idea, but this
    # doesn't really go to high enough altitude in the inner disk. Instead, let's choose a minimum threshold.
    if (10 * zq) >= 10 * AU
        outer_reaches = zq * 10
    else
        outer_reaches = 10 * AU
    end

    zints = cat(1, [0], logspace(log10(0.01 * AU), log10(outer_reaches), nz-1))

    un_lnrhos = Array(Float64, nz)
    for i=1:nz
        un_lnrhos[i] = un_lnrho(r, zints[i], pars)
    end

    bar_lnrho = un_lnrhos[1]

    # Remove the largest value so we don't get an overflow error
    un_lnrhos -= bar_lnrho

    # Convert from ln(rho) to rho, but this is still unnormalized.
    un_rhos = exp(un_lnrhos)

    # Now that we have the shape of the density structure but not the normalization,
    # integrate to get the normalization
    # The density should be unitful, so we need to exponentiate.
    spl = Spline1D(zints, un_rhos)

    # Integral is from -inf to inf
    tot = 2 * integrate(spl, 0., zints[end])

    # println("Norm. tot: ", tot, " lntot:", log(tot))
    # println("bar_lnrho. exp: ", exp(bar_lnrho), " bar_lnrho: ", bar_lnrho)
    # println("With bar_lnrho :", tot * exp(bar_lnrho), " lntot :", log(tot) + bar_lnrho)

    # The correction factor is the ratio of the actual surface density to the value of the integral of our unnormalized density.
    lncor = log(S) - (log(tot) + bar_lnrho)

    # cor = S/tot
    # apply this as
    # rhos = cor * unnormed_rhos
    return lncor
end

# Now, create an interpolator for the the correction function. In spherical coordinates, that way, first evaluate it at every point (r, z=0) in the grid. Then, when we want to query the density of the disk at (r, z), which will correspond to a different r_cyl than the one used to compute the normalization, we can just use the spline interpolation and not have to redo the time consuming normalization integral.
function make_correction_interpolator(pars::Parameters, grid::Grid)

    nr = 64
    rs = logspace(log10(0.5 * AU), log10(grid.Rs[end]), nr)

    cors = Array(Float64, nr)

    for i=1:nr
        cors[i] = lncorrection_factor(rs[i], pars)
    end

    # An interpolating spline for this correction factor.
    spl = Spline1D(rs, cors)

    return spl
end

# spl is the spline object produced by a call to make_correction_interpolator()
function rho_gas(r::Float64, z::Float64, pars::Parameters, spl)

    # If we are inner to the radius used to compute the correction factors (in this case 0.5 AU),
    # just return 0.0 gas density
    if r < (0.5 * AU)
        return 0

    else
        # Evaluate the correction factor for this radius.
        lncor = evaluate(spl, r)

        # Calculate the unnormalized ln rho
        un = un_lnrho(r, z, pars)

        # rho = cor * exp(un)
        # apply the correction factor
        lnrho = lncor + un

        # Exponentiate it
        rho = exp(lnrho)

        # Check to make sure that this didn't return NaN, due to overflow or underflow errors somewhere in the integration pipeline. RADMC-3D will run with a NaN input to synthesize a blank image, and so the MCMC chain will gladly go on synthesizing blank images, even though the model is nonsense.
        if isnan(rho)
            println("Returned NaN for density.")
            throw(OverflowError)
        else
            return rho
        end
    end
end

n_12CO(r::Float64, z::Float64, pars::Parameters, spl) = number_densities["12CO"] * rho_gas(r, z, pars, spl)

n_13CO(r::Float64, z::Float64, pars::Parameters, spl) = number_densities["13CO"] * rho_gas(r, z, pars, spl)

n_C18O(r::Float64, z::Float64, pars::Parameters, spl) = number_densities["C18O"] * rho_gas(r, z, pars, spl)

# It is realistic to include freezout of CO onto dust grains.
# This is the amount by which the number density of the CO is reduced (X_freeze) relative to
# the nominal value.
function X_freeze(temp::Float64, pars::Parameters)
    # If it's cooler than the freezout temperature, reduce the number density by the given factor
    if temp <= pars.T_freeze
        return pars.X_freeze
    # Otherwise, just keep it as is.
    else
        return 1.0
    end
end


# Ksi is microturbulent broadining width in units of km/s. Output of this function
# is in cm/s for RADMC (RADMC manual, eqn 7.12)
function microturbulence(ksi::Float64)
    return ksi * 1.e5 # convert from km/s to cm/s
end

microturbulence(pars::Parameters) = microturbulence(pars.ksi)

function write_model(pars::Parameters, basedir::AbstractString, grid::Grid, species::AbstractString)

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

    # Calculate the correction function interpolator
    spl = make_correction_interpolator(pars, grid)

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
                XF = X_freeze(temp, pars)
                @printf(fdens, "%.9e\n", XF * n_CO(r_cyl, z, pars, spl))
                @printf(fvel, "0 0 %.9e\n", velocity(r_cyl, pars))
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

function write_dust(pars::Parameters, basedir::AbstractString, grid::Grid)
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
