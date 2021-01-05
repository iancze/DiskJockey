# This is file provides the ability for the user to override the default prior.

# This code is `include()`'ed inside of `/scripts/venus.jl`, so it will **literarly replace** the
# function lnprior inside `src/model.jl` (all of it!). Be careful to make sure it accepts the same arguments!

using DiskJockey.model
import DiskJockey.model.lnprior_base
import DiskJockey.constants.ModelException

# for example, a replacement for the standard model might look like
function lnprior(pars::ParametersStandard, grid::Grid)
    lnp = lnprior_base(pars)

    # Modify this code here to enforce your prior!

    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.
    if (3 * pars.r_c) > r_out
        throw(ModelException("User prior violated."))
    end

    return lnp
end

# while an example for the cavity model would look like this
function lnprior(pars::ParametersCavity, grid::Grid)
    lnp = lnprior_base(pars)

    # Modify this code here to enforce your prior!

    r_in = grid.Rs[1]/AU # [AU]
    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.

    # Also check to make sure that r_cav is less than r_c but larger than r_in.
    if (3 * pars.r_c) > r_out || pars.r_cav < r_in || pars.r_cav > pars.r_c
        throw(ModelException("User prior violated."))
    end

    return lnp

end


function lnprior(pars::ParametersVertical, grid::Grid)
    # Create a giant short-circuit or loop to test for sensical parameter values.
    if pars.M_star <= 0.0 || pars.ksi <= 0. || pars.T_10a <= 0. || pars.T_10m <= 0. || pars.r_c <= 0.0  || pars.T_10a > 1500. || pars.q_m < 0. || pars.q_a < 0. || pars.q_m > 1.0 || pars.q_a > 1.0 || pars.incl < 0. || pars.incl > 180. || pars.PA < -180. || pars.PA > 520. || pars.X_freeze > 1.0 || pars.sigma_s < 0.0
        throw(ModelException("Parameters outside of prior range."))
    end

    # Check to see that the temperatures make sense
    if pars.T_10m > pars.T_10a
        throw(ModelException("Atmosphere is cooler than midplane."))
    end


    # If we've passed all the hard-cut offs by this point, return the geometrical inclination prior.
    lnp = log(0.5 * sind(pars.incl))

    r_out = grid.Rs[end]/AU # [AU]
    # A somewhat arbitrary cutoff regarding the gridsize to prevent the disk from being too large
    # to fit on the model grid.
    if (3 * pars.r_c) > r_out
        throw(ModelException("Model radius too large for grid size."))
    else
        return lnp
    end

end
