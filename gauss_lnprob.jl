# lnprob evaluation for Gaussian

using visibilities
using image
using gauss_model
using gridding
using constants

const dvarr = DataVis("data/V4046Sgr_fake.hdf5") #Read in the whole fake dataset
nlam = length(dvarr)

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

    # Compute the array of images
    data = Array(Float64, (ny, nx, nlam))
    lams = Array(Float64, nlam)
    for k=1:nlam
        data[:, :, k] = imageGauss(ll, mm, p, k)
        lams[k] = dvarr[k].lam # [microns]
    end

    skim = SkyImage(data, ra, dec, lams)

    # Apply the gridding correction function before doing the FFT
    corrfun!(skim, 1.0) #alpha = 1.0

    # Now calculate the sum of lnprobs in a loop
    lnp = Array(Float64, nlam)
    for i=1:nlam
        # FFT the appropriate image chanel
        vis_fft = transform(skim, i)

        # Interpolate the `vis_fft` at the same locations as the DataSet
        mvis = ModelVis(dvarr[i], vis_fft)

        # Calculate chi^2 between these two
        lnp[i] = lnprob(dvarr[i], mvis)
    end
    return sum(lnp)
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

# # Now try optimizing the function using NLopt
using NLopt

starting_param = [1.3, 1.2, 0.7, 0.7]

nparam = length(starting_param)
opt = Opt(:LN_COBYLA, nparam)

max_objective!(opt, fgrad)
xtol_rel!(opt,1e-4)

(optf,optx,ret) = optimize(opt, starting_param)
println(optf, " ", optx, " ", ret)



using LittleMC

using Distributions
using PDMats

mc = MC(fp, 500, [1.2, 1.0, 1.0, 1.0], PDiagMat([0.03^2, 0.03^2, 0.01^2, 0.01^2]))

start(mc)


println(mean(mc.samples, 2))
println(std(mc.samples, 2))

runstats(mc)

write(mc, "mc.hdf5")
