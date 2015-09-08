module EnsembleSampler

using NPZ

# This is a direct Julia port of [emcee](http://dan.iel.fm/emcee/current/), the ensemble sampler by Dan Foreman-Mackey et al.
#
# Right now, this is only designed to work in parallel with cores on the same node in the most straightforward example. Eventually it would be nice to incorporate Julia tasks.

export Sampler, run_mcmc, run_schedule, flatchain, reset, write_samples

# There is a type, called the sampler.
type Sampler
    nwalkers::Int
    ndim::Int
    lnprobfn::Function
    a::Float64
    chain # Array of shape (nwalkers, iterations, dim) in python.
    # For us, this is (ndim, iterations, nwalkers) in Julia, since we will be filling in
    # an entire row of ndim at once for a single iteration, making it the fastest changing dimension.
    lnprob # storing the values of the lnprob. # lnprob in python is (nwalkers, niterations)
    # in julia, it's the same
    iterations::Int # How many iterations has the ensemble taken
end

function Sampler(nwalkers::Int, ndim::Int, lnprobfn::Function)
    a = 2.0
    chain = Array(Float64, (ndim, 0, nwalkers))
    lnprob = Array(Float64, (nwalkers, 0))
    iterations = 0

    sampler = Sampler(nwalkers, ndim, lnprobfn, a, chain, lnprob, iterations)
end

# Convert the nwalkers chain into a flatchain
function flatchain(sampler::Sampler)
    # chain is stored as (ndim, iterations, nwalkers)

    # flatchain should be stored as (ndim, iterations)
    fchain = reshape(sampler.chain, (sampler.ndim, sampler.iterations * sampler.nwalkers))

    return fchain
end

# Clear the samples after burn in
function reset(sampler::Sampler)
    # self.naccepted = np.zeros(self.k)
    sampler.chain = Array(Float64, (sampler.ndim, 0, sampler.nwalkers))
    sampler.lnprob = Array(Float64, (sampler.nwalkers, 0))
    sampler.iterations = 0
end

function sample(sampler::Sampler, p0, lnprob0=nothing, iterations=1)

    p = p0

    # What is the index to divide the number of walkers in half
    halfk = ifloor(sampler.nwalkers / 2)

    # If the initial log-probabilities were not provided, calculate them
    # now.
    lnprob = lnprob0
    if lnprob == nothing
        lnprob = get_lnprob(sampler, p)
    end

    # Check to make sure that the probability function didn't return
    # ``np.nan``.
    if any(isnan(lnprob))
        println("The initial lnprob was NaN.")
        throw(DomainError())
    end

    # Store the initial size of the stored chain.
    i0 = size(sampler.chain)[2]
    # resize the chain by adding new rows to the array.

    # (ndim, iterations, nwalkers)
    sampler.chain = cat(2, sampler.chain, Array(Float64, (sampler.ndim, iterations, sampler.nwalkers)))

    #(niterations, nwalkers)
    sampler.lnprob = cat(2, sampler.lnprob, Array(Float64, sampler.nwalkers, iterations))

    for i=1:iterations
        sampler.iterations += 1
        println("Iteration ", sampler.iterations)
        # Loop over the two ensembles, calculating the proposed positions.

        # Slices for the first and second halves
        first = 1:halfk
        second = (halfk + 1):sampler.nwalkers

        # Do this for both halves, alternating
        for (S0, S1) in [(first, second), (second, first)]
            q, newlnp, acc = propose_stretch(sampler, p[:, S0], p[:, S1], lnprob[S0])

            if any(acc)
                # Update the positions and log probabilities

                # Note that the indexing needed to be done differently than emcee, since double
                # indexing in Julia does not address the same memory in the array.
                lnprob[S0[acc]] = newlnp[acc]
                p[:, S0[acc]] = q[:, acc]

            end

            ind = i0 + i
            sampler.chain[:, ind, :] = p
            sampler.lnprob[:, ind] = lnprob

        end

    end
    # Return the current position
    return p
end


# Propose a new position for one sub-ensemble given the positions of
# another.
# :param p0:
#  The positions from which to jump.
# :param p1:
#  The positions of the other ensemble.
# :param lnprob0:
#  The log-probabilities at ``p0``.
# This method returns:
# * ``q`` - The new proposed positions for the walkers in ``ensemble``.
# * ``newlnprob`` - The vector of log-probabilities at the positions
# given by ``q``.
# * ``accept`` - A vector of type ``bool`` indicating whether or not
# the proposed position for each walker should be accepted..

function propose_stretch(sampler::Sampler, p0, p1, lnprob0)
    # In python, these are (nwalkers, ndim)
    # In Julia, these are (ndim, nwalkers)

    s = p0
    Ns = size(s)[2]
    c = p1
    Nc = size(c)[2]

    # Ns and Nc refer to the number of walkers in the subset and complimentary set, respectively

    # Generate the vectors of random numbers that will produce the
    # proposal.

    # self._random.rand(Ns) provides a 1D array of values in the range [0,1.) of size Ns
    zz = ((sampler.a - 1.0) .* rand(Ns) .+ 1.0) .^ 2. ./ sampler.a

    # An array of random integers the size of the subset, designed to slice into the complimentary sample.
    rint = rand(1:Nc, Ns)

    # Index into the complementary sample, and advance based upon the scale factor.
    # c[:, rint] is a 2D array (ndim, Ns)
    q = c[:, rint] - (reshape(zz, (1, Ns)) .*  (c[:, rint] - s))

    newlnprob = get_lnprob(sampler, q)

    # Decide whether or not the proposals should be accepted.
    lnpdiff = (sampler.ndim - 1.) .* log(zz) .+ newlnprob .- lnprob0
    accept = (lnpdiff .> log(rand(Ns)))

    return q, newlnprob, accept

end

# Calculate the vector of log-probability for the walkers.
# :param pos: (optional)
#     The position vector in parameter space where the probability
#     should be calculated. This defaults to the current position
#     unless a different one is provided.

#     (ndim, nwalkers)
# This method returns:
# * ``lnprob`` - A vector of log-probabilities with one entry for each
#   walker in this sub-ensemble.
function get_lnprob(sampler::Sampler, pos)

    p = pos

    # Check that the parameters are in physical ranges.
    if any(isinf(p))
        println("At least one parameter value was infinite.")
        throw(DomainError())
    end
    if any(isnan(p))
        println("At least one parameter value was NaN.")
        throw(DomainError())
    end

    # What is the shape of p here in this subset?
    nwalkers_sub = size(pos)[2]
    lst = [p[:,i] for i=1:nwalkers_sub]

    # In Python, it seems like each row corresponds to a different walker.
    # For Julia, we really want each column to the parameters corresponding to a different walker.
    lnprob = float(pmap(sampler.lnprobfn, lst))

    # Check for lnprob returning NaN.
    if any(isnan(lnprob))
        # Print some debugging stuff.
        println("NaN value of lnprob for parameters: ")
            for pars in p[isnan(lnprob)]
                println(pars)
            end
        # Finally raise exception.
        throw(DomainError())
    end

    return lnprob

end

# Advance the Sampler for N iterations, staring from pos0.
function run_mcmc(sampler::Sampler, pos0, N::Int)
    pos = sample(sampler, pos0, nothing, N)
    return pos
end

# Write out the samples
function write_samples(sampler::Sampler, outdir="")

    fchain = flatchain(sampler)
    npzwrite(outdir * "chain.npy", sampler.chain)
    npzwrite(outdir * "flatchain.npy", fchain)
    npzwrite(outdir * "lnprob.npy", sampler.lnprob)

    # Needs to be reshaped to remove singleton dimension
    npzwrite(outdir * "pos0.npy", reshape(sampler.chain[:, end, :], (ndim, nwalkers)))

end

# Run the sampler on a periodic save schedule, to prevent losses on long-running calculations
# run N iterations for each loop
function run_schedule(sampler::Sampler, pos0, N::Int, loops::Int, outdir)
    for i=1:loops
        pos0 = sample(sampler, pos0, nothing, N)
        println("Finished loop ", i, " of ", loops)
        write_samples(sampler, outdir)
    end
    return pos0
end


end # module
