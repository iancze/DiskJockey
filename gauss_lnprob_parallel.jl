# lnprob evaluation for Gaussian

nchild = 3
addprocs(nchild)

@everywhere using parallel
@everywhere using visibilities
@everywhere using image
@everywhere using gauss_model
@everywhere using gridding
@everywhere using constants

# Dimensions of the problem
@everywhere const nx = 128
@everywhere const ny = 128
@everywhere const ra = fftspace(6., nx) # [arcsec]
@everywhere const dec = fftspace(6., ny) # [arcsec]

# convert ra and dec in [arcsec] to radians, and then take the sin to convert to ll, mm
@everywhere ll = sin(ra * arcsec) # direction cosines
@everywhere mm = sin(dec * arcsec)

@everywhere function initfunc(key)
    println("Loading dataset $key")
    return DataVis("data/V4046Sgr_fake.hdf5", key)
end

@everywhere function f(dv::DataVis, key::Int, p::Vector{Float64})
    # p is naturally in arcsec

    # Compute the image
    img = imageGauss(ll, mm, p, key)
    skim = SkyImage(img, ra, dec, dv.lam)

    # Apply the gridding correction function before doing the FFT
    corrfun!(skim, 1.0) #alpha = 1.0

    # FFT the appropriate image channel
    vis_fft = transform(skim)

    # Interpolate the `vis_fft` at the same locations as the DataSet
    mvis = ModelVis(dv, vis_fft)

    # Calculate chi^2 between these two
    return lnprob(dv, mvis)
end

keylist = Int[1, 2, 3]
pipes = initialize(nchild, keylist, initfunc, f)

function fprob(p::Vector{Float64})
    distribute!(pipes, p)
    return gather!(pipes)
end

println(fprob([1.2, 1.0, 1.2, 1.0]))

# Necessary wrapper for NLopt requires gradient as an argument (even if it's not used)
function fgrad(p::Vector, grad::Vector)
    val = fprob(p)
    println(p, " : ", val)
    return val
end
#
# function fp(p::Vector)
#     val = f(p)
#     println(p, " : ", val)
#     return val
# end

# # Now try optimizing the function using NLopt
using NLopt

starting_param = [1.3, 1.2, 0.7, 0.7]

nparam = length(starting_param)
opt = Opt(:LN_COBYLA, nparam)

max_objective!(opt, fgrad)
xtol_rel!(opt,1e-4)

(optf,optx,ret) = optimize(opt, starting_param)
println(optf, " ", optx, " ", ret)



# using LittleMC
#
# using Distributions
# using PDMats
#
# mc = MC(fp, 500, [1.2, 1.0, 1.0, 1.0], PDiagMat([0.03^2, 0.03^2, 0.01^2, 0.01^2]))
#
# start(mc)
#
#
# println(mean(mc.samples, 2))
# println(std(mc.samples, 2))
#
# runstats(mc)
#
# write(mc, "mc.hdf5")

quit!(pipes)
