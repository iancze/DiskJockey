module LittleMC

# Because there is no package to do just what I want, this is
# a simple Metropolis-Hastings sampler to just get things done.
# inspired by github.com/dfm/emcee/mh.py

export MC, sample, start, runstats

using Distributions
using PDMats
using HDF5

# The main object that contains all of the necessary MCMC information
type MC{T <: AbstractMvNormal} #Could try making this immutable later
    f::Function # Function to evaluate, generally the log-posterior
    nsamples::Int
    p0::Vector{Float64} # Starting parameters
    proposal::T # Jump proposal covariance matrix
    nparams::Int
    samples::Matrix{Float64}
    lnprobs::Vector{Float64}
    naccepted::Int
    csv::IOStream
end

# Initialization function
function MC(f::Function, nsamples::Int, p0::Vector{Float64}, propcov::Matrix{Float64}, csv::IOStream)
    # Do some custom initialization
    nparams = length(p0)

    # Make sure we didn't screw up with matrices
    @assert(size(propcov)[1] == size(propcov)[2], "propcov must be square matrix")
    @assert(nparams == size(propcov)[1] , "propcov mismatched with nparams")

    samples = Array{Float64(nparams, nsamples)
    lnprobs = Array{Float64}(nsamples)

    # Turn propcov matrix into a distribution
    proposal = MvNormal(propcov)

    MC(f, nsamples, p0, proposal, nparams, samples, lnprobs, 0, csv)
end

function sample(mc::MC, p0::Vector{Float64}, lnprob0::Float64)

    # Get parameter proposal
    q = p0 .+ rand(mc.proposal)

    # Evaluate the proposal
    lnprob_proposed = mc.f(q)

    diff = lnprob_proposed - lnprob0

    # M-H acceptance ratio
    if diff < 0
        diff = exp(diff) - rand()
    end

    if diff >= 0
        # proposal accepted
        mc.naccepted += 1
        return (q, lnprob_proposed)
    else
        # proposal rejected
        return (p0, lnprob0)
    end

    return (p, lnprob)
end

# Start running the MCMC
function start(mc::MC)
    # Evaluate at the current parameters
    p0 = mc.p0
    lnprob0 = mc.f(p0)

    # always accept this as a starting point
    mc.samples[:, 1] = p0
    mc.lnprobs[1] = lnprob0

    writecsv(mc.csv, p0')

    for i=2:mc.nsamples
        (p0, lnprob0) = sample(mc, p0, lnprob0)
        mc.samples[:, i] = p0
        mc.lnprobs[i] = lnprob0
        writecsv(mc.csv, p0')
        flush(mc.csv) # Make sure this gets written.
    end
end

# Calculate the Metropolis acceptance fraction
function acceptance(mc::MC)
    return mc.naccepted/mc.nsamples
end

# Calculate the covariance of the samples
# function covariance(mc::MC)
# end

# Wrapper function for all the analysis tools
function runstats(mc::MC)
    @printf("Acceptance fraction %.2f\n", acceptance(mc))
    # more to come
end


# Optionally write the samples to an HDF5 file
function write(mc::MC, fname::AbstractString)
    fid = h5open(fname, "w")
    fid["samples"] = mc.samples

    # Add the acceptance fraction as an attribute on the "samples" dataset
    attrs(fid["samples"])["acceptance"] = acceptance(mc)

    close(fid)
end

end #module
