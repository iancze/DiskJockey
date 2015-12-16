module ImportanceSampler

# Samples p(x), to compute f(x), using the sampling distribution q(x).

using HDF5

# Writes samples to a file.

export ISampler, sampleIS, write_samples

"""Define the important functions for the sampler.
p must be a Probability density function, that takes in a single vector of parameter values.
f must be the function you want to evaluate, that takes in a single vector of parameter values.
q must be the sampling distribution, that produces a vector of parameter values. For best results, it should be tuned to match (p * f), but oversample the tails somewhat (e.g., using a student-t distribution when sampling a Gaussian).
"""
type ISampler
    nsamples::Int
    nparam::Int
    p::Function # lnprobability density
    f::Function # Function we wish to compute over this interval. For our inference problems,
    # this is typically a constant.
    q::Function # The proposal distribution, returns lnprobability
    xx::Matrix{Float64}
    p_x::Vector{Float64}
    f_x::Vector{Float64}
    q_x::Vector{Float64}
end

function ISampler(nsamples::Int, nparam::Int, p::Function, q::Function)
    # A dummy function that returns a constant
    f(x) = 1
    return ISampler(nsamples, nparam, p, f, q, Array(Float64, (nparam, nsamples)), Array(Float64, nsamples), Array(Float64, nsamples), Array(Float64, nsamples))
end

"""Load existing samples from an HDF5 file. Put dummy functions for functions since we won't be using them now."""
function ISampler(fname::AbstractString)
    fid = h5open(fname, "r")

    xx = read(fid["xx"])
    p_x = read(fid["p_x"])
    f_x = read(fid["f_x"])
    q_x = read(fid["q_x"])

    nparam, nsamples = size(xx)
    close(fid)

    f(x) = x

    return ISampler(nsamples, nparam, f, f, f, xx, p_x, f_x, q_x)

end

function sampleIS(sampler::ISampler)
    for i=1:sampler.nsamples
        # Draw a parameter combination from q.
        # Also report the likelihood at q(x)
        x, q_x = sampler.q()

        sampler.xx[:,i] = x
        sampler.q_x[i] = q_x

        # Evaluate p(x) and f(x)
        sampler.p_x[i] = sampler.p(x)
        sampler.f_x[i] = sampler.f(x)
    end

end

"""Write the values of the sampler to an HDF5 file."""
function write_samples(sampler::ISampler, fname::AbstractString)
    fid = h5open(fname, "w")
    fid["xx"] = sampler.xx
    fid["p_x"] = sampler.p_x
    fid["f_x"] = sampler.f_x
    fid["q_x"] = sampler.q_x
    close(fid)
end

# Functions to calculate moments.


end # module
