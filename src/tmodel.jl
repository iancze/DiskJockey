module tmodel

# truncated power law model

export write_grid, write_model, write_lambda, write_dust, Parameters, Grid

# The double dot is because we are now inside the model module, and we want to import the
# constants module, which is part of the enclosing JudithExcalibur package.
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


#Define the 7 parameters listed in Rosenfeld et al., then write functions to go from those parameters to dust density, velocity, microturbulence.

#Let's try defining a parameters type, an object of which gets passed around.

type Parameters
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

# Typical values, from Rosenfeld et al. 2012, Table 1 for V4046 Sgr
M_star = 1.75 # [M_sun] stellar mass
r_c =  45. # [AU] characteristic radius
T_10 =  115. # [K] temperature at 10 AU
q = 0.63 # temperature gradient exponent
gamma = 1.0 # surface temperature gradient exponent
M_CO = 0.933 # [M_earth] disk mass of CO
ksi = 0.14e5 # [cm s^{-1}] microturbulence
dpc = 73.
incl = 33.5 # [degrees] inclination
PA = 73.
vel = 0.0
mu_RA = 0.0
mu_DEC = 0.0

# Assume all inputs to these functions are in CGS units and in *cylindrical* coordinates.
# Parametric type T allows passing individual Float64 or Vectors.
# Alternate functions accept pars passed around, where pars is in M_star, AU, etc...
function velocity{T}(r::T, M_star::Float64)
    sqrt(G * M_star ./ r)
end
velocity{T}(r::T, pars::Parameters) = velocity(r, pars.M_star * M_sun)

function temperature{T}(r::T, T_10::Float64, q::Float64)
    T_10 * (r ./ (10. * AU)).^(-q)
end
temperature{T}(r::T, pars::Parameters) = temperature(r, pars.T_10, pars.q)

function Hp{T}(r::T, M_star::Float64, T_10::Float64, q::Float64)
    temp = temperature(r, T_10, q)
    sqrt(kB * temp .* r.^3./(mu_gas * m_H * G * M_star))
end
Hp{T}(r::T,  pars::Parameters) = Hp(r, pars.M_star * M_sun, pars.T_10, pars.q)

# Calculate the gas surface density
# function Sigma(r::Float64, pars::Parameters)
#     r_c = pars.r_c * AU
#     Sigma_c = pars.M_gas * M_sun * (2 - pars.gamma) / (2 * pi * r_c^2)
#     Sigma_c * (r/r_c)^(-pars.gamma) * exp(-(r/r_c)^(2 - pars.gamma))
# end

# Calculate the gas surface density
function Sigma{T}(r::T, pars::Parameters)
    r_c = pars.r_c * AU
    Sigma_c = pars.M_gas * M_sun * (2 - pars.gamma) / (2 * pi * r_c^2)
    Sigma_c .* (r./r_c).^(-pars.gamma) .* exp(-(r./r_c).^(2 - pars.gamma))
end

# Delivers a gas density in g/cm^3
function rho_gas(r::Float64, z::Float64, pars::Parameters)
    H = Hp(r, pars)
    S = Sigma(r, pars)
    S/(sqrt(2. * pi) * H) * exp(-0.5 * (z/H)^2)
end

# Now, replace these functions to simply multiply rho_gas by X_12CO/m_12CO, or X_13CO/m_13CO, etc.

n_12CO(r::Float64, z::Float64, pars::Parameters) = number_densities["12CO"] * rho_gas(r, z, pars)

n_13CO(r::Float64, z::Float64, pars::Parameters) = number_densities["13CO"] * rho_gas(r, z, pars)

n_C18O(r::Float64, z::Float64, pars::Parameters) = number_densities["C18O"] * rho_gas(r, z, pars)

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

function rho_dust(r::Float64, z::Float64, pars::Parameters)
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
