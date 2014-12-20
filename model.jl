
module model

export write_grid, write_model, Parameters, params, lambdas, rs, thetas

using constants #Import all of my constants defined in constants.jl

#These grid parameters are fixed for the start of each run, they depend on the dataset.

# Write out the camera wavelength file
# Line frequency center, CO J = 2-1
# Compute the various wavelengths at which we want to synthesize images
# 23 images spaced -4.40 -0- 4.40 km/s
lam0 = cc/230.538e9 * 1e4 # [microns]
nvels = 23
vels = linspace(-4.4, 4.4, nvels) # [km/s]
const lambdas = vels/c_kms * lam0 + lam0

# Write the wavelength sampling file. Only run on setup
function write_lambda(lams::Array{Float64, 1})
    fcam = open("camera_wavelength_micron.inp", "w")
    nlam = length(lams)
    @printf(fcam, "%d\n", nlam)
    for lam in lams
        @printf(fcam, "%.9e\n", lam) # [microns]
    end

    close(fcam)
end

const eqmirror = true # mirror the grid about the z=0 midplane ?

# Specify a 2D axisymmetric *separable* grid in spherical coordinates, {r, theta, phi}.
# theta is angle from zenith, phi is azimuth

# Number of cells in each dimension
# if we decide to mirror, then ncells = 1/2 of the true value
const nr = 64
const ntheta = 32 # if mirror about the equator
const nphi = 1

const ncells = nr * ntheta * nphi

const r_in = 5 * AU # Inner extent of disk
const r_out = 500 * AU # Outer extent of disk

#Define the cell *walls*
const Rs = logspace(log10(r_in), log10(r_out), nr+1) # [cm] logarithmically spaced

if eqmirror
ped = 0.1
#Thetas = linspace(0, pi/2., ntheta+1)  # [rad] Angles are internally defined in radians, not degrees
const Thetas = pi/2. - (logspace(log10(ped), log10(pi/2. + ped), ntheta+1) - ped)[end:-1:1] #Spaced closer near the z=0
else
const Thetas = linspace(0, pi, ntheta+1)  # [rad] Angles are internally defined in radians, not degrees
end

const Phis = Float64[0.0, 0.0] # [rad] cell walls for inactive coordinate

#Define the cell centers as the average between walls
const rs = 0.5 * (Rs[1:end-1] + Rs[2:end])
const thetas = 0.5 * (Thetas[1:end-1] + Thetas[2:end])
const phis = Float64[0.0]

#This function only needs to be run once, upon setup.
function write_grid()
    #amr_grid.inp
    f = open("amr_grid.inp", "w")

    #Write the header
    @printf(f, "%d\n", 1) #iformat
    @printf(f, "%d\n", 0) #regular grid (no AMR or Oct-tree)
    @printf(f, "%d\n", 100) #spherical coordiantes
    @printf(f, "%d\n", 0) #gridinfo (none needed for now)
    #incl_r incl_phi incl_z #use this axis?
    @printf(f, "%d %d %d \n", 1, 1, 0) # 2D axisymmetric
    #n_r    n_phi   n_z #number of cells in this dimension
    @printf(f, "%d %d %d \n", nr, ntheta, nphi)

    for R in Rs
        @printf(f, "%.9e\n", R)
    end

    for Theta in Thetas
        @printf(f, "%.9e\n", Theta)
    end

    for Phi in Phis
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
    M_CO::Float64 # [M_earth] disk mass of CO
    ksi::Float64 # [cm s^{-1}] microturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_x::Float64 # [arcsec] central offset in RA
    mu_y::Float64 # [arcsec] central offset in DEC
end

#From Rosenfeld et al. 2012, Table 1
M_CO = 0.933 # [M_earth] disk mass of CO
r_c =  45. # [AU] characteristic radius
T_10 =  115. # [K] temperature at 10 AU
q = 0.63 # temperature gradient exponent
gamma = 1.0 # surface temperature gradient exponent
ksi = 0.14e5 # [cm s^{-1}] microturbulence
incl = 33.5 # [degrees] inclination
M_star = 1.75 # [M_sun] stellar mass
#PA = 73.

# global object which is useful for reproducing V4046Sgr
# params = Parameters(M_CO, r_c, T_10, q, gamma, ksi, i_d, M_star)

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

# No parametric type for number density, because it is a 2D function.
function n_CO(r::Float64, z::Float64, r_c::Float64, M_CO::Float64, M_star::Float64, T_10::Float64, q::Float64, gamma::Float64)
    H = Hp(r, M_star, T_10, q)
    (2. - gamma) * M_CO/(m_CO * (2. * pi)^(1.5) * r_c^2 * H) * (r/r_c)^(-gamma) * exp(-0.5 * (z/H)^2 - (r/r_c)^(2. - gamma))
end
n_CO(r::Float64, z::Float64, pars::Parameters) = n_CO(r, z, pars.r_c * AU, pars.M_CO * M_earth, pars.M_star * M_sun, pars.T_10, pars.q, pars.gamma)


function write_model(pars::Parameters)
    # numberdens_co.inp
    fdens = open("numberdens_co.inp", "w")
    @printf(fdens, "%d\n", 1) #iformat
    @printf(fdens, "%d\n", ncells)

    # gas_velocity.inp
    fvel = open("gas_velocity.inp", "w")
    @printf(fvel, "%d\n", 1) #iformat
    @printf(fvel, "%d\n", ncells)

    # gas_temperature.inp
    ftemp = open("gas_temperature.inp", "w")
    @printf(ftemp, "%d\n", 1) #iformat
    @printf(ftemp, "%d\n", ncells)

    #microturbulence.inp

    # Now, we will need to write the three other files as a function of grid position.
    # Therefore we will do *one* loop over these indices, calculate the required value,
    # and write it to the appropriate file.

    #Looping over the cell centers
    for phi in phis
        for theta in thetas
            for r in rs
                #Convert from spherical to cylindrical coordinates
                z = r * cos(theta)
                r_cyl = r * sin(theta)

                @printf(fdens, "%.9e\n", n_CO(r_cyl, z, pars))
                @printf(fvel, "0 0 %.9e\n", velocity(r_cyl, pars))
                @printf(ftemp, "%.9e\n", temperature(r_cyl, pars))
            end
        end
    end

    close(fdens)
    close(fvel)
    close(ftemp)

end


end
