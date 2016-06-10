# This file is meant to be imported in other python scripts to help ease the transition from using logM_gas as a parameter to using Sigma_c as a parameter. This change was done so that we can efficiently use more complicated models which don't have an analytic form for the integration of the surface density profile.

import numpy as np
from scipy.integrate import quad

M_sun = 1.99e33 # [g]
AU = 1.4959787066e13 # [cm]

# The conversions for the standard and vertical models are straightforward, since the integral to find total mass can be done analytically.

# Each of these functions assumes that inputs are in units of either AU or Solar Masses

def sigma_to_logM_standard(r_c, gamma, Sigma_c):
    r_c = r_c * AU # [cm]

    M_gas = 2 * np.pi * r_c**2 * Sigma_c / (2 - gamma)
    M_gas = M_gas / M_sun # [M_sun]
    logM = np.log10(M_gas)

    return logM

def logM_to_sigma_standard(r_c, gamma, logM):
    r_c = r_c * AU

    M_gas = 10**logM * M_sun # [g]
    Sigma_c = M_gas * (2 - gamma) / (2 * np.pi * r_c**2) # [g/cm^2]

    return Sigma_c


# On the other hand, conversions for the cavity model are more difficult, because a numerical integral is needed to convert between the two.

def integral(r_c, r_cav, gamma, gamma_cav):
    r_c = r_c * AU
    r_cav = r_cav * AU

    def f(r):
        return r * (r/r_c)**(-gamma) * np.exp(-(r/r_c)**(2 - gamma) - (r_cav/r)**gamma_cav)

    # For some reason, using np.inf doesn't work properly for the integral, so we use a very large outer radius instead.
    I = quad(f, 0.0, 10000 * AU)[0] # [cm^2]
    return I

def sigma_to_logM_cavity(r_c, r_cav, gamma, gamma_cav, Sigma_c):
    I = integral(r_c, r_cav, gamma, gamma_cav)

    M_gas = 2 * np.pi * Sigma_c * I # [g]
    M_gas = M_gas / M_sun # [M_sun]
    logM = np.log10(M_gas)

    return logM

def logM_to_sigma_cavity(r_c, r_cav, gamma, gamma_cav, logM):
    I = integral(r_c, r_cav, gamma, gamma_cav)

    M_gas = 10**logM * M_sun # [g]
    Sigma_c = M_gas / (2 * np.pi * I) # [g/cm^2]

    return Sigma_c


# Use these dictionaries to get functions.

sigma_to_logM = {"standard":sigma_to_logM_standard, "vertical":sigma_to_logM_standard, "cavity":sigma_to_logM_cavity}

logM_to_sigma = {"standard":logM_to_sigma_standard, "vertical":logM_to_sigma_standard, "cavity":logM_to_sigma_cavity}


def main():
    # Some unit tests
    logM = -2.0

    r_c = 20.0
    r_cav = 10.0
    gamma = 1.0
    gamma_cav = 2.0

    sigma_c_standard = logM_to_sigma_standard(r_c, gamma, logM)

    logM_standard = sigma_to_logM_standard(r_c, gamma, sigma_c_standard)

    print("standard")
    print("logM original", logM)
    print("logM to Sigma_c", sigma_c_standard)
    print("Sigma_c to logM", logM_standard)

    sigma_c_cavity = logM_to_sigma_cavity(r_c, r_cav, gamma, gamma_cav, logM)

    logM_cavity = sigma_to_logM_cavity(r_c, r_cav, gamma, gamma_cav, sigma_c_cavity)

    print()
    print("cavity")
    print("logM original", logM)
    print("logM to Sigma_c", sigma_c_cavity)
    print("Sigma_c to logM", logM_cavity)


if __name__=="__main__":
    main()
