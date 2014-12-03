#Physical constants

module constants

export M_sun, AU, pc, G, kB, c_ang, cc, c_kms, mu_gas, m_H, m_CO, arcsec, fftspace

# Conversion from astronomical units to CGS units
M_sun = 1.99e33 # [g]
AU = 1.4959787066e13 # [cm]
pc = 3.0856776e18 # [cm]
G = 6.67259e-8 # [cm3 g-1 s-2]
kB = 1.380658e-16 # [erg K^-1] Boltzmann constant
c_ang = 2.99792458e18 # [A s^-1]
cc = 2.99792458e10 # [cm s^-1]
c_kms = 2.99792458e5 # [km s^-1]

mu_gas = 2.37 # mean molecular weight of circumstellar gas
m_H = 1.6733e-24 # g mass of hydrogen atom

#Atomic
amu = 1.6605402e-24 # [g]

#CO 
m_CO = 28.01 * amu #molecular weight of CO in g

#Disk specific
m0 = mu_gas * amu #mean molecular weight of gas

# convert from arcseconds to radians
arcsec = pi / (180. * 3600) # [radians]  = 1/206265 radian/arcsec

# Oftentimes it is necessary to get a symmetric coordinate array that spans N elements from -width to +width, 
# but makes sure that the middle point lands on 0. The indices go from 0 to N -1.
# `linspace` returns  the end points inclusive, wheras we want to leave out the right endpoint, 
# because we are sampling the function in a cyclic manner.
function fftspace(width::Real, N::Int)
    @assert(N % 2 == 0, "N must be even.")

    dx = width * 2. / N
    xx = Array(Float64, N)
    for i=1:N
        xx[i] = -width + (i - 1) * dx
    end
    return xx
end


end
