# lnprob evaluation for Gaussian

using visibilities
using image
using gauss_model
using gridding
using constants

const dvis = DataVis("data/V4046Sgr_fake.hdf5", 1) #Read in the fake dataset

# Dimensions of the problem
const nx = 128
const ny = 128

const ra = fftspace(6., nx) # [arcsec]
const dec = fftspace(6., ny) # [arcsec]

# convert ra and dec in [arcsec] to radians, and then take the sin to convert to ll, mm
ll = sin(ra * arcsec) # direction cosines
mm = sin(dec * arcsec)


function f(p::Vector{Float64})
    # p is naturally in arcsec
    img = imageGauss(ll, mm, p)
    lam0 = cc/230.538e9 * 1e4 # [microns]
    skim = SkyImage(img, ra, dec, lam0)

    # Apply the gridding correction function before doing the FFT
    corrfun!(skim, 1.0) #alpha = 1.0

    # FFT the image
    vis_fft = transform(skim)

    # Interpolate the `vis_fft` at the same locations as the DataSet
    mvis = ModelVis(dvis, vis_fft)

    # Calculate chi^2 between these two
    return lnprob(dvis, mvis)

end

# Necessary wrapper for NLopt requires gradient as an argument (even if it's not used)
function fgrad(p::Vector, grad::Vector)
    val = f(p)
    println(p, " : ", val)
    return val
end

function fp(p::Vector)
    val = f(p)
    println(p, " : ", val)
    return val
end

# Now try optimizing the function using NLopt
# using NLopt
#
# nparam = 2
# opt = Opt(:LN_COBYLA, nparam)
#
# max_objective!(opt, fgrad)
# xtol_rel!(opt,1e-4)
#
# (optf,optx,ret) = optimize(opt, [1.3, 1.0])
# println(optf, " ", optx, " ", ret)


using LittleMC

using Distributions
using PDMats

mc = MC(fp, 200, [1.2, 1.0], PDiagMat([0.06^2, 0.05^2]))

start(mc)


println(mean(mc.samples, 2))
println(std(mc.samples, 2))

runstats(mc)

write(mc, "mc.hdf5")
