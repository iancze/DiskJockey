# Models

DiskJockey supports many different disk models, and is designed in such a way so as to be easily extensible to new ones.

The current models include

## Standard

Keyword: `standard`

This is a standard model that is used in Czekala et al. 2015 and 2016, and Rosenfeld et al. 2012.

## Vertical

Keyword: `vertical`

This model includes a more realistic temperature profile as outlined in Dartois et al. 2003, Rosenfeld et al. 2013, and Williams and Best 2014, among many others. We follow the specific parameterization in Williams and Best 2014.

The primary additional step for this model is the necessity to numerically solve the hydrostatic equilibrium equation (WB14 Eqn 1).

* for a given radius, `r`, come up with a vertical grid of `z` points that stretches from the midplane to an appropriate height where the gas density is approximately zero.
* calculate `ln(rho(r,z))` as a function of `(r, z)` by integrating WB14 Eqn 1 using `quadgk`. This is the *unnormalized* density of the disk, which means the relative values at a fixed radius should be correct, but no guarantees about anything else. We will call the unnormalized density `un_rho(r,z)`.
* now, we know that `Sigma(r)` specifies the *surface density* of the disk, i.e., or the integral of `rho(r, z)` from `z = -inf to +inf`. Therefore, we can find `rho(r,z)` by vertically integrating `un_rho(r,z)` and finding the normalization constant, `norm(r)`. Because we are dealing with integrals of very big and very small numbers, we need to do some tricks to avoid overflow and underflow errors.
* finally, we know that CO is photodissociated when it is not shielded by enough gas, i.e., there is not enough gas column density above this height to block harmful radiation. Therefore, we need to find the photodissociation height, `z_phot(r)`, above which the molecule CO can no longer exist. To do this requires integrating from `z = +inf` towards the midplane (`z = 0`) to find the height at which we have accumulated enough gas column density to be shielded.

This sounds like a lot of steps just to evaluate a single `rho(r,z)` point. Because RADMC-3D solves the radiative transfer on a spherical grid and the disk model is defined on a cylindrical grid, there is a careful order of operations necessary to achieve the appropriate accuracy in the shortest amount of computational time. We address this by first solving everything on a cylindrical grid to find `norm(r)` and `z_phot(r)` as a function of disk radius. Then, for a given `(r_spherical, theta)` point, we convert to cylindrical coordinates and solve WB14 Eqn 1 to find `un_rho(r,z)`.

## VerticalEta

Characterized with `ParametersVerticalEta`. In addition to the parameters described in the `Vertical` model, this extension has an additional parameter `eta`, designed to vary the height of the atmosphere with radius.

## Nuker 

The Nuker profile as used in [Tripathi et al. 2017.](https://ui.adsabs.harvard.edu/abs/2017ApJ...845...44T/abstract).
