using JudithExcalibur.EnsembleSampler
using Base.Test

# write your own tests here
@test 1 == 1

# Rosenbrock function
@everywhere function lnprob(p)
    x, y = p

    return -((1.0 - x)^2 + 100 * (y - x^2)^2)
end


# A 2D Gaussian to check against
# @everywhere function lnprob(p)
#
#     a, b = p
#
#     if a < 0.0
#         return -Inf
#     end
#
#     mu = [1.0, -2.0]
#     Sigma = [[0.5^2,  0.00] [0.00, 0.3^2]]
#
#     R = p - mu # Residual vector
#     ld = logdet(Sigma)
#
#     lnp = -0.5 * ((R' * (Sigma\R))[1] + ld + 2 * log(2 * pi))
#     return lnp
#
# end


nwalkers = 40
ndim = 2

sampler = Sampler(nwalkers, ndim, lnprob)

# pos0 is the starting position, it needs to be a (ndim, nwalkers array)
pos0 = Array(Float64, ndim, nwalkers)
for i=1:nwalkers
    pos0[:,i] = randn(ndim)
end

# println("Starting positions ", pos0)

pos = run_mcmc(sampler, pos0, 10000)

reset_mcmc(sampler)

run_mcmc(sampler, pos, 10000)

using NPZ

npzwrite("chain.npy", emcee_chain(sampler))
