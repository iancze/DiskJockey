# This is file provides the ability for the user to override the default prior.

# This code is `include()`'ed inside of `/scripts/venus.jl`, so it will **literarly replace** the
# function lnprior inside `src/model.jl` (all of it!). Be careful to make sure it accepts the same arguments!

using DiskJockey.model
import DiskJockey.model.lnprior_base
import DiskJockey.constants.ModelException

# for example, a replacement for the standard model might look like
function lnprior(pars::ParametersStandard, dpc_mu::Float64, dpc_sig::Float64, grid::Grid)
    lnp = lnprior_base(pars, dpc_mu, dpc_sig)

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
function lnprior(pars::ParametersCavity, dpc_mu::Float64, dpc_sig::Float64, grid::Grid)
    lnp = lnprior_base(pars, dpc_mu, dpc_sig)

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
