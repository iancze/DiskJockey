#Physical constants

module constants

export M_sun, AU, pc, G, kB, c_ang, cc, c_kms, mu_gas, m_H, m_CO

# Conversion from astronomical units to CGS units
M_sun = 1.99e33 # [g]
AU = 1.4959787066e13 # [cm]
pc = 3.0856776e18 # [cm]
G = 6.67259e-8 # [cm3 g-1 s-2]
kB = 1.380658e-16 # [erg K^-1] Boltzmann constant
c_ang = 2.99792458e18 # [A s^-1]
cc = 2.99792458e8 # [m s^-1]
c_kms = 2.99792458e5 # [km s^-1]

mu_gas = 2.37 # mean molecular weight of circumstellar gas
m_H = 1.6733e-24 # g mass of hydrogen atom

#Atomic
amu = 1.6605402e-24 # [g]

#CO 
m_CO = 28.01 * amu #molecular weight of CO in g

#Disk specific
m0 = mu_gas * amu #mean molecular weight of gas

end
