# Conversion from astronomical units to CGS units
M_sun = 1.99e33 #g
AU = 1.4959787066e13 #cm
pc = 3.0856776e18 #cm
G = 6.67259e-8 #cm3 g-1 s-2
kB = 1.380658e-16 # erg K^-1 Boltzmann constant
c_ang = 2.99792458e18 #A s^-1
cc = 2.99792458e8 #m s^-1
c_kms = 2.99792458e5 #km s^-1

# Write out the camera wavelength file
#Line frequency center, CO J = 2-1
lam0 = cc/230.538e9 * 1e6 # [microns]
nvels = 23
vels = linspace(-4.4, 4.4, nvels) # [km/s]
lams = vels/c_kms * lam0 + lam0
fcam = open("camera_wavelength_micron.inp", "w")
@printf(fcam, "%d\n", nvels)
for lam in lams
    @printf(fcam, "%.9e\n", lam) # [microns]
end

close(fcam)

# Specify a 2D axisymmetric *separable* grid in spherical coordinates, {r, theta, phi}.
# theta is angle from zenith, phi is azimuth

# Number of cells in each dimension
nr = 128
ntheta = 128
nphi = 1

r_in = 5 * AU # Inner extent of disk
r_out = 400 * AU # Outer extent of disk

#Define the cell *walls*
Rs = logspace(log10(r_in), log10(r_out), nr+1) # [cm] logarithmically spaced
Thetas = linspace(0, pi, ntheta+1)  # [rad] Angles are internally defined in radians, not degrees
Phis = Float64[0.0, 0.0] # [rad] cell walls for inactive coordinate

#Define the cell centers as the average between walls
rs = 0.5 * (Rs[1:end-1] + Rs[2:end])
thetas = 0.5 * (Thetas[1:end-1] + Thetas[2:end])
phis = Float64[0.0]

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


#Define the 7 parameters listed in Rosenfeld et al., then write functions to go from those parameters to dust density, velocity, microturbulence.

#From Rosenfeld et al. 2012, Table 1
M_CO = 2.8e-6 * M_sun # disk mass of CO
r_c =  45. * AU # characteristic radius
T_10 =  115. # [K] temperature at 10 AU
q = 0.63 # temperature gradient
gamma = 1.0 # surface temperature gradient
ksi = 0.14 # [km s^{-1}]microturbulence
i_d = 33.5 # [degrees] inclination
M_star = 1.75 * M_sun # [g] stellar mass
mu_gas = 2.37 # mean molecular weight of circumstellar gas
m_H = 1.6733e-24 # g mass of hydrogen atom

# numberdens_co.inp
fdens = open("numberdens_co.inp", "w")
@printf(fdens, "%d\n", 1) #iformat
@printf(fdens, "%d\n", nr * ntheta * nphi ) #number of cells

# gas_velocity.inp
fvel = open("gas_velocity.inp", "w")
@printf(fvel, "%d\n", 1) #iformat
@printf(fvel, "%d\n", nr * ntheta * nphi ) #number of cells

# gas_temperature.inp
ftemp = open("gas_temperature.inp", "w")
@printf(ftemp, "%d\n", 1) #iformat
@printf(ftemp, "%d\n", nr * ntheta * nphi ) #number of cells

#microturbulence.inp

#Assume all inputs to these functions are in CGS units and in *cylindrical* coordinates
function velocity(r::Float64, M_star::Float64)
    sqrt(G * M_star/r)
end

function temperature(r::Float64, T_10::Float64, q::Float64)
    T_10 * (r/ (10. * AU))^(-q)
end

function Hp(r::Float64, M_star::Float64, T_10::Float64, q::Float64)
    T = temperature(r, T_10, q)
    sqrt(kB * T * r^3/(mu_gas * m_H * G * M_star))
end

function n_CO(r::Float64, z::Float64, r_c::Float64, M_CO::Float64, M_star::Float64, T_10::Float64, q::Float64, gamma::Float64)
    H = Hp(r, M_star, T_10, q)
    (2. - gamma) * M_CO/((2. * pi)^(1.5) * r_c^2 * H) * (r/r_c)^(-gamma) * exp(-0.5 * (z/H)^2 - (r/r_c)^(2. - gamma))
end

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

            @printf(fdens, "%.9e\n", n_CO(r_cyl, z, r_c, M_CO, M_star, T_10, q, gamma))
            @printf(fvel, "0 0 %.9e\n", velocity(r_cyl, M_star)) 
            @printf(ftemp, "%.9e\n", temperature(r_cyl, T_10, q))
        end
    end
end

close(fdens)
close(fvel)
close(ftemp)
