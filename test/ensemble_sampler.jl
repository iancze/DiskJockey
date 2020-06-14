using DiskJockey.EnsembleSampler
using Distributed
using LinearAlgebra: logdet


# The EnsembleSampler was written to *only* work in parallel, so in order to test it properly, we 
# need to also define all of our test lnprobs to work using Distributed,
# even though they are quick enough to work on their own.

# Rosenbrock function
@everywhere function lnprob(p)
    x, y = p

    return -((1.0 - x)^2 + 100 * (y - x^2)^2)
end

@testset "sample Rosenbrock" begin
    nwalkers = 40
    ndim = 2

    sampler = Sampler(nwalkers, ndim, lnprob)

    # pos0 is the starting position, it needs to be a (ndim, nwalkers array)
    pos0 = Array{Float64}(undef, ndim, nwalkers)
    for i = 1:nwalkers
        pos0[:,i] = randn(ndim)
    end

    run_schedule(sampler, pos0, 10, 1, outputdir)
end

# Now, let's try sampling a function with more parameters than walkers (just for testing purposes)
# A 4D Gaussian to check against
@everywhere function lnprob(p)

    mu = [1.0, -2.0, 0.0, 3.4]
    Sigma = [[0.5^2,  0.00, 0.0, 0.0] [0.00, 0.3^2, 0.0, 0.0] [0.0, 0.0, 0.7^2, 0.0] [0.0, 0.0, 0.0, 1.0^2] ]

    R = p - mu # Residual vector
    ld = logdet(Sigma)

    lnp = -0.5 * ((R' * (Sigma \ R))[1] + ld + 2 * log(2 * pi))
    return lnp

end

@testset "sample Gaussian" begin
    nwalkers = 4
    ndim = 4

    sampler = Sampler(nwalkers, ndim, lnprob, true)

    # pos0 is the starting position, it needs to be a (ndim, nwalkers array)
    pos0 = Array{Float64}(undef, ndim, nwalkers)
    for i = 1:nwalkers
        pos0[:,i] = randn(ndim)
    end

    run_schedule(sampler, pos0, 1, 1, outputdir)
end