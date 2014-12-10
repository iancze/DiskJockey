module LittleMC

# Because there is no package to do just what I want, this is
# a simple Metropolis-Hastings sampler to just get things done.
# inspired by github.com/dfm/emcee/mh.py

export MC, sample, run, write, runstats

using Distributions
using PDMats

# The main object that contains all of the necessary MCMC information
type MC{T <: AbstractMvNormal} #Could try making this immutable later
    f::Function # Function to evaluate, generally the log-posterior
    nsamples::Int
    p0::Vector{Float64} #Starting parameters
    proposal::T #May need to be changed
    nparams::Int
    samples::Matrix{Float64}
    lnprobs::Vector{Float64}
    naccepted::Int
end

# Initialization function
function MC{T <: AbstractPDMat}(f::Function, nsamples::Int, p0::Vector{Float64}, propcov::T)
    # Do some custom initialization
    nparams = length(p0)

    # Make sure we didn't screw up with matrices
    @assert(size(propcov)[1] == size(propcov)[2], "propcov must be square matrix")
    @assert(nparams == size(propcov)[1] , "propcov mismatched with nparams")

    samples = Array(Float64, nparams, nsamples)
    lnprobs = Array(Float64, nsamples)

    #Turn propcov matrix into a distribution
    proposal = MvNormal(propcov)

    MC(f, nsamples, p0, proposal, nparams, samples, lnprobs, 0)
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

    for i=2:mc.nsamples
        (p0, lnprob0) = sample(mc, p0, lnprob0)
        mc.samples[:, i] = p0
        mc.lnprobs[i] = lnprob0
    end
end

# Calculate the Metropolis acceptance fraction
function acceptance(mc::MC)
    return mc.naccepted/mc.nsamples
end

# Calculate the covariance of the samples
#function covariance(mc::MC)
#end

# Wrapper function for all the analysis tools
function runstats(mc::MC)
    @printf("Acceptance fraction %.2f\n", acceptance(mc))
    # more to come
end

const mat = PDiagMat([2.0^2, 3.0^2])
const dist = DiagNormal([0., 0.], mat)
function Gauss(x::Vector{Float64})
    return logpdf(dist, x)
end

mc = MC(Gauss, 100000, [1.0, 1.0], PDiagMat([1.5^2, 1.5^2]))
start(mc)


println(mean(mc.samples, 2))
println(std(mc.samples, 2))

println(acceptance(mc))

# Plot the samples
import PyPlot.plt
fig, ax = plt.subplots(nrows=2, sharex=true)
ax[1][:plot](mc.samples'[:, 1])
ax[2][:plot](mc.samples'[:, 2])
plt.savefig("plots/littlemc.png")

# Optionally write the samples to an HDF5 file

end #module
