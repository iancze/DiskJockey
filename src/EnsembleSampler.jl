module EnsembleSampler

using Distributed
using NPZ

# This is a direct Julia port of [emcee](http://dan.iel.fm/emcee/current/), the ensemble sampler by Dan Foreman-Mackey et al.
#
# Right now, this is only designed to work in parallel with cores on the same node in the most straightforward example. Eventually it would be nice to incorporate Julia tasks.

export Sampler, run_mcmc, run_schedule, reset_mcmc, write_samples, emcee_chain

# There is a type, called the sampler.
mutable struct Sampler
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

function Sampler(nwalkers::Int, ndim::Int, lnprobfn::Function, test::Bool=false)
    a = 2.0

    if !test
        @assert ndim < nwalkers "I don't believe you have more parameters than walkers! Try flipping them."
    end

    chain = Array{Float64}(undef, (ndim, 0, nwalkers))
    lnprob = Array{Float64}(undef, (nwalkers, 0))
    iterations = 0

    sampler = Sampler(nwalkers, ndim, lnprobfn, a, chain, lnprob, iterations)
end

# Transpose the chain order from Julia column-major to Python row-major
function emcee_chain(sampler::Sampler)
    # ndim, niter, nwalkers = size(sampler.chain)
    # return reshape(sampler.chain, (nwalkers, niter, ndim))
    return permutedims(sampler.chain, [3, 2, 1])
end

# Convert the nwalkers chain into a flatchain
# function flatchain(sampler::Sampler)
#     # chain is stored as (ndim, iterations, nwalkers) in Julia parlance.
#
#     # flatchain should be stored as (ndim, iterations)
#     fchain = reshape(sampler.chain, (sampler.ndim, sampler.iterations * sampler.nwalkers))
#
#     return fchain
# end

# Clear the samples after burn in
function reset_mcmc(sampler::Sampler)
    # self.naccepted = np.zeros(self.k)
    sampler.chain = Array{Float64}(sampler.ndim, 0, sampler.nwalkers)
    sampler.lnprob = Array{Float64}(sampler.nwalkers, 0)
    sampler.iterations = 0
end

function sample(sampler::Sampler, p0, lnprob0=nothing, iterations=1)

    p = p0

    # What is the index to divide the number of walkers in half
    halfk = floor(Int, sampler.nwalkers / 2)

    # If the initial log-probabilities were not provided, calculate them
    # now.
    # This step is calculating the lnprob for the total number of walkers.
    lnprob = lnprob0
    if lnprob == nothing
        lnprob = get_lnprob(sampler, p)
    end

    # Check to make sure that the probability function didn't return
    # ``np.nan``.
    if any(isnan.(lnprob))
        println("The initial lnprob was NaN.")
        throw(DomainError())
    end

    # Store the initial size of the stored chain.
    i0 = size(sampler.chain)[2]
    # resize the chain by adding new rows to the array.

    # (ndim, iterations, nwalkers)
    sampler.chain = cat(2, sampler.chain, Array{Float64}(undef, sampler.ndim, iterations, sampler.nwalkers))

    #(niterations, nwalkers)
    sampler.lnprob = cat(2, sampler.lnprob, Array{Float64}(undef, sampler.nwalkers, iterations))

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
    zz = ((sampler.a - 1.0) .* rand(Ns) .+ 1.0) .^ 2.0 ./ sampler.a

    # An array of random integers the size of the subset, designed to slice into the complimentary sample.
    rint = rand(1:Nc, Ns)

    # Index into the complementary sample, and advance based upon the scale factor.
    # c[:, rint] is a 2D array (ndim, Ns)
    q = c[:, rint] - (reshape(zz, (1, Ns)) .*  (c[:, rint] - s))

    # This get_lnprob is only calculating the lnprob for the proposed positions, which is Ns, and is probably half the total number of walkers.
    newlnprob = get_lnprob(sampler, q)

    # Decide whether or not the proposals should be accepted.
    lnpdiff = (sampler.ndim - 1.0) .* log.(zz) .+ newlnprob .- lnprob0
    accept = (lnpdiff .> log.(rand(Ns)))

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
    if any(isinf.(p))
        println("At least one parameter value was infinite.")
        throw(DomainError())
    end
    if any(isnan.(p))
        println("At least one parameter value was NaN.")
        throw(DomainError())
    end

    # What is the shape of p here in this subset?
    nwalkers_sub = size(pos)[2]
    lst = [p[:,i] for i=1:nwalkers_sub]

    # In Python, it seems like each row corresponds to a different walker.
    # For Julia, we really want each column to the parameters corresponding to a different walker.
    result = pmap(sampler.lnprobfn, lst)

    # The array may contain one or two RemoteExceptions, so let's check to see what caused them.
    for (res,par) in zip(result, lst)
        if typeof(res) == RemoteException
            println("RemoteException found for parameters ", par)
            println(res)
            throw(DomainError)
        end
    end

    # If we've made it to here, everything's ok.
    lnprob = convert(Array{Float64, 1}, result)

    # Check for lnprob returning NaN.
    if any(isnan.(lnprob))
        # Print some debugging stuff.
        println("NaN value of lnprob for parameters: ")
            for pars in p[isnan.(lnprob)]
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

    # fchain = flatchain(sampler)

    npzwrite(outdir * "chain.npy", emcee_chain(sampler))
    # npzwrite(outdir * "flatchain.npy", fchain)
    npzwrite(outdir * "lnprob.npy", sampler.lnprob)

    # Needs to be reshaped to remove singleton dimension
    npzwrite(outdir * "pos0.npy", reshape(sampler.chain[:, end, :], (sampler.ndim, sampler.nwalkers)))

end

# Run the sampler on a periodic save schedule, to prevent losses on long-running calculations
# run N iterations for each loop
function run_schedule(sampler::Sampler, pos0, N::Int, loops::Int, outdir, func=nothing)
    for i=1:loops
        pos0 = sample(sampler, pos0, nothing, N)
        println("Finished loop ", i, " of ", loops)
        write_samples(sampler, outdir)
        if func != nothing
            func(sampler, outdir)
        end

    end
    return pos0
end


end # module
