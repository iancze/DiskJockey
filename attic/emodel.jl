module emodel

# Uses an eccentric disk

export write_grid, write_model, write_lambda, write_dust, Parameters, Grid

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
    As::Vector{Float64}
    # cell centers
    rs::Vector{Float64}
    thetas::Vector{Float64}
    phis::Vector{Float64}
    as::Vector{Float64}
    a_widths::Vector{Float64}
end

function Grid(nr::Int, ntheta::Int, nphi::Int, r_in::Real, r_out::Real, na::Int, eqmirror::Bool)
    # Specify a 2D axisymmetric *separable* grid in spherical coordinates:
    # {r, theta, phi}, where theta is angle from zenith, phi is azimuth

    # na is the number of eccentric bins that we want to have. Generally, I
    # think this should be less than nr.

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

    Phis = linspace(0, 2pi, nphi + 1) # [rad] cell walls for inactive coordinate

    # Define the semi-major axis bin edges such that we fully encapsulate the full range of radii
    # Assume 0 <= e <= 0.7
    eps = 0.1 * AU # [cm], tiny extra bit of comfort
    e = 0.7 # Maximum e we'll tolerate
    A_1 = (Rs[1] - eps)/(1 + e)
    A_N = (Rs[end] + eps)/(1 - e)

    As = logspace(log10(A_1), log10(A_N), na+1)


    #Define the cell centers as the average between walls
    rs = 0.5 * (Rs[1:end-1] + Rs[2:end])
    thetas = 0.5 * (Thetas[1:end-1] + Thetas[2:end])
    phis = 0.5 * (Phis[1:end-1] + Phis[2:end])

    as = 0.5 * (As[1:end-1] + As[2:end])
    a_widths = As[2:end] - As[1:end-1]

    return Grid(nr, ntheta, nphi, ncells, Rs, Thetas, Phis, As, rs, thetas, phis, as, a_widths)

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
    @printf(f, "%d %d %d \n", 1, 1, 1) # Full 3D grid
    #n_r    n_phi   n_z #number of cells in this dimension
    @printf(f, "%d %d %d \n", grid.nr, grid.ntheta, grid.nphi)

    # Write out the cell walls
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
    a_c::Float64 # [AU] characteristic radius
    T_10::Float64 # [K] temperature at 10 AU
    q::Float64 # temperature gradient exponent
    gamma::Float64 # surface temperature gradient exponent
    M_CO::Float64 # [M_earth] disk mass of CO
    ksi::Float64 # [cm s^{-1}] microturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    e::Float64 # eccentricity
    w::Float64 # [degrees] argument of periastron
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
end


# Given an array of bin edges (N) in radius from the star, and an array of radii (xx), return an
# array of indices of which bins (N-1) the radii fall into.
# Assumes both bins and xx are sorted in increasing order.
function get_bins(bin_edges, xx)
    inds = Array(Int64, length(xx))
    b::Int64 = 1
    nb = length(bin_edges)
    for (i,x) in enumerate(xx)
        while true
            # Greater than the lower bin
            if x >= bin_edges[b]

                # Smaller than the upper bin
                if x < bin_edges[b+1]
                    inds[i] = b
                    break

                # larger than the upper bin, means we need to try a larger bin index
                else
                    b += 1
                    if b + 1 > nb
                        throw(error(@sprintf("x value % s out of bin range.", x)))
                    end
                end

            # Smaller than the lower bin, something went wrong?
            else
                throw(error(@sprintf("x value %s out of bin range.", x)))

            end
        end

    end
    return inds
end

# Assume all inputs to these functions are in CGS units and in *cylindrical* coordinates.

# Assume all angles are in *radians*, including argument of periapse

# phi is azimuthal angle, a is semi-major axis, e is eccentricity
function radius{T}(phi::T, a::Float64, e::Float64, w::Float64)
    a * (1 - e^2)./(1 + e .* cos(phi - w))
end

radius{T}(phi::T, a::Float64, pars::Parameters) = radius(phi, a, pars.e, pars.w * deg)

# Linear density of the ring for the mass-on-a-wire approximation
# m is the total mass of the ring
function lambda(phi::Float64, a::Float64, e::Float64, w::Float64, m::Float64)
    m * sqrt(1 - e^2) / (2pi * a * sqrt(1 + 2e * cos(phi - w) + e^2))
end

# Parametric type T allows passing individual Float64 or Vectors.
# Alternate functions accept pars passed around, where pars is in M_star, AU, etc...
function velocity(phi::Float64, a::Float64, e::Float64, w::Float64, M_star::Float64)
    v_r = sqrt(G * M_star / a) * e * sin(phi - w) / sqrt(1 - e^2)
    v_theta = 0.0
    v_phi = sqrt(G * M_star / a) * (1 + e * cos(phi - w)) / sqrt(1 - e^2)
    return Float64[v_r, v_theta, v_phi]
end
velocity(phi::Float64, a::Float64, pars::Parameters) = velocity(phi, a, pars.e, pars.w * deg, pars.M_star * M_sun)


function temperature{T}(r::T, T_10::Float64, q::Float64)
    T_10 * (r ./ (10. * AU)).^(-q)
end
temperature{T}(r::T, pars::Parameters) = temperature(r, pars.T_10, pars.q)

function Hp{T}(r::T, M_star::Float64, T_10::Float64, q::Float64)
    temp = temperature(r, T_10, q)
    sqrt(kB * temp .* r.^3./(mu_gas * m_H * G * M_star))
end
Hp{T}(r::T,  pars::Parameters) = Hp(r, pars.M_star * M_sun, pars.T_10, pars.q)

# No parametric type for number density, because it is a 2D function.
function n_CO(r::Float64, z::Float64, M_star::Float64, T_10::Float64, q::Float64)
    H = Hp(r, M_star, T_10, q)
    return 1./ (m_CO * sqrt(2pi) * H) * exp(-0.5 * (z/ H)^2)
end

n_CO(r::Float64, z::Float64, pars::Parameters) = n_CO(r, z, pars.M_star * M_sun, pars.T_10, pars.q)


function n_nobin(r::Float64, phi::Float64, z::Float64, M_star::Float64, a_c::Float64, T_10::Float64, q::Float64, gamma::Float64, M_CO::Float64, e::Float64, w::Float64)
    n = n_CO(r, z, M_star, T_10, q)
    a = r * (1 + e * cos(phi - w))/(1 - e^2)
    Sigma = (2 - gamma) * M_CO / (2pi * a * a_c) * (a/a_c)^(1 - gamma) * exp(- (a/a_c)^(2 - gamma)) * (1 + e * cos(phi - w))/(sqrt(1 + 2 * e * cos(phi - w)) * sqrt(1 - e^2))
    return Sigma * n
end

n_nobin(r::Float64, phi::Float64, z::Float64, pars::Parameters) = n_nobin(r, phi, z, pars.M_star * M_sun, pars.a_c * AU, pars.T_10, pars.q, pars.gamma, pars.M_CO * M_earth, pars.e, pars.w * deg)

# Calculate the mass in each ring using a_c and M_co
function ring_mass(as::Vector{Float64}, a_c::Float64, M_CO::Float64, gamma::Float64)
    a = as ./ a_c
    ms = M_CO .* a.^(-gamma) .* exp(-a.^(2 - gamma)) / sum(a.^(-gamma) .* exp(-a.^(2 - gamma)))
    return ms
end

ring_mass(grid, pars) = ring_mass(grid.as, pars.a_c * AU, pars.M_CO * M_earth, pars.gamma)


function rho_dust(r::Float64, z::Float64, pars::Parameters)
    nCO = n_CO(r, z, pars) # number of CO molecules per cm^3

    # Convert from nCO to nH2
    nH2 = nCO / 7.e-5 # number density ratio

    # Convert from nH2 (assuming nH2 ~ nGas ) to mGas
    mGas = constants.m0 * nH2 # [g]

    # Convert from mGas to mDust using Gas/Dust ratio of 100
    mDust = mGas * 0.01 # [g]

    return mDust
end

# Ksi is microturbulent broadining width in units of km/s. Output of this function
# is in cm/s for RADMC (RADMC manual, eqn 7.12)
function microturbulence(ksi::Float64)
    return ksi * 1.e5 # convert from km/s to cm/s
end

microturbulence(pars::Parameters) = microturbulence(pars.ksi)

function write_model(pars::Parameters, basedir::AbstractString, grid::Grid)
    # numberdens_co.inp
    fdens = open(basedir * "numberdens_co.inp", "w")
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

    # Calculate the mass in each ring using a_c and M_co

    #Looping over the cell centers
    for phi in grid.phis
        # Calculate the ellipse bin-edges
        # R_edges = grid.As .* (1 - pars.e^2)/(1 + pars.e * cos(phi - pars.w * deg))
        #
        # # println("grid.As", grid.As)
        # # println("R_edges", R_edges)
        # # println("grid.rs", grid.rs)
        #
        # # Determine the array of bin indices which correspond to each point in grid.rs
        # inds = get_bins(R_edges, grid.rs)
        #
        # # Calculate all following phi-dependent quantities in vectorized form:
        # # m_i, a_i, dr_i, Sigma_CO_i
        #
        # m_rings = ring_mass(grid, pars)
        #
        # # Calculate arrays of these values for each radius
        # ms = [m_rings[ind] for ind in inds]
        # as = [grid.as[ind] for ind in inds]
        # das = [grid.a_widths[ind] for ind in inds]
        #
        # drs = das .* (1 - pars.e^2)/(1 + pars.e * cos(phi - pars.w * deg))
        #
        # lambdas = ms ./ (2pi .* as) .* (sqrt(1 - pars.e^2) ./ sqrt(1 + 2 * pars.e .* cos(phi - pars.w * deg) + pars.e^2))
        #
        # Sigmas = lambdas ./ drs

        for theta in grid.thetas

            for (i,r) in enumerate(grid.rs)
                #Convert from spherical to cylindrical coordinates
                z = r * cos(theta)
                r_cyl = r * sin(theta)

                a = r_cyl * (1 + pars.e * cos(phi - pars.w * deg))/(1 - pars.e^2)
                #Calculate unscaled density as function of (r_cyl, z)
                # n = Sigmas[i] * n_CO(r_cyl, z, pars)

                v = velocity(phi, a, pars)

                # @printf(fdens, "%.9e\n", n)
                @printf(fdens, "%.9e\n", n_nobin(r_cyl, phi, z, pars))
                @printf(fvel, "%.9e %.9e %.9e\n", v[1], v[2], v[3])
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
