# Provides FFT functionality, regridding,

module visibilities

using AbstractFFTs
using FFTW
using HDF5
using ..image
using ..gridding
using ..constants

import Base.conj! # extend this for DataVis
import Base.- # extend this for FullModelVis

export DataVis, ModelVis, RawModelVis, FullModelVis, fillModelVis, ResidVis, ModelVisRotate, DataVisReal, ModelVisReal
export get_qq
export plan_interpolate, interpolate_uv
export transform, rfftfreq, fftfreq, phase_shift, phase_shift!, max_baseline, get_nyquist_pixel
export lnprob, lnprob_iter, lnprob_real
export -

using Test

# Stores the data visibilities for a single channel
mutable struct DataVis
    lam::Float64 # [μm] Wavelength (in microns) corresponding to channel
    uu::Vector{Float64} # [kλ] Vectors of the u, v locations in kilolambda, shape (nchan, nvis)
    vv::Vector{Float64} # [kλ]
    VV::Vector{ComplexF64} # [Jy] Complex visibilities
    invsig::Vector{Float64} # 1/sigma [1/Jy]
    # NOTE that there is no flags field
end

"
    DataVis(fname::AbstractString, flagged::Bool=false)

Read all the visibilities from an HDF5 string and return them as an array of DataVis objects
`flagged` means whether we should be loading the flagged points or not."
function DataVis(fname::AbstractString, flagged::Bool=false)
    fid = h5open(fname, "r")

    freqs = read(fid["freqs"]) # [Hz]
    # Convert from Hz to wavelengths in μm
    lams = cc ./ freqs * 1e4 # [μm]

    uu = read(fid, "uu")
    vv = read(fid, "vv")
    real = read(fid, "real")
    imag = read(fid, "imag")
    VV = real + imag .* im # Complex visibility
    weight = read(fid, "weight") # [1/Jy^2]

    if !flagged
        flag = convert(Array{Bool}, read(fid, "flag"))
    end

    close(fid)

    nlam = length(lams)

    # Do we want to include the flagged visibility points? In nearly all cases, the answer
    # should be no. The one exception is when we want to export model visibilities in `write_model.jl`

    # Do the flagging channel by channel

    # Return an array of DataVis
    out = Array{DataVis}(undef, nlam)
    if !flagged
        for i=1:nlam
            ch_flag = .~flag[:,i]
            out[i] = DataVis(lams[i], uu[ch_flag,i], vv[ch_flag,i], VV[ch_flag, i], sqrt.(weight[ch_flag, i]))
        end
    else
        for i=1:nlam
            out[i] = DataVis(lams[i], uu[:,i], vv[:,i], VV[:, i], sqrt.(weight[:, i]))
        end
    end
    return out
end

"
    DataVis(fname::AbstractString, index::Int, flagged::Bool=false)

Read just one channel of visibilities from the HDF5 file."
function DataVis(fname::AbstractString, index::Int, flagged::Bool=false)
    fid = h5open(fname, "r")
    # the indexing and `vec` are necessary here because HDF5 doesn't naturally
    # squeeze trailing dimensions of length-1

    freq = read(fid["freqs"])[index][1] # [Hz]
    # Convert from Hz to wavelengths in μm
    lam = cc / freq * 1e4 # [μm]

    len = length(fid["uu"][:, index])

    if !flagged
        flag = convert(Array{Bool}, read(fid["flag"]))
        ch_flag = .~flag[:,index]
    else
        ch_flag = ones(Bool, len) # Keep all visibilities, regardless of what flag says
    end

    uu = vec(fid["uu"][:,index])[ch_flag]
    vv = vec(fid["vv"][:,index])[ch_flag]
    real = vec(fid["real"][:,index])[ch_flag]
    imag = vec(fid["imag"][:,index])[ch_flag]

    VV = real + imag .* im # Complex visibility

    invsig = vec(sqrt.(fid["weight"][:,index]))[ch_flag] # [1/Jy^2]

    close(fid)

    #Return a DataVis object
    return DataVis(lam, uu, vv, VV, invsig)
end

"
    DataVis(fname::AbstractString, indices::Vector{Int}, flagged::Bool=false)

Read just a subset of channels from the HDF5 file and return an array of DataVis."
function DataVis(fname::AbstractString, indices::Vector{Int}, flagged::Bool=false)
    nchan = length(indices)
    out = Array{DataVis}(undef, nchan)
    for i=1:nchan
        out[i] = DataVis(fname, indices[i], flagged)
    end
    return out
end

"
    DataVis(fname::AbstractString, mask::BitArray, flagged::Bool=false)

Read just a subset of channels from the HDF5 file and return an array of DataVis."
function DataVis(fname::AbstractString, mask::BitArray{1}, flagged::Bool=false)
    # load the array and mask it
    return DataVis(fname, flagged)[mask]
end

struct DataVisReal
    lam::Float64 # [μm] Wavelength (in microns) corresponding to channel
    uu::Vector{Float64} # [kλ] Vectors of the u, v locations in kilolambda, shape (nchan, nvis)
    vv::Vector{Float64} # [kλ]
    VVr::Vector{Float64} # [Jy] Complex visibilities
    VVi::Vector{Float64} # [Jy] Complex visibilities
    invsig::Vector{Float64} # 1/sigma [1/Jy]
    # NOTE that there is no flags field
end

"
    DataVisReal(dvis::DataVis)

Convert a `DataVis` object, with visibilities stored as `Comlpex128`, into an object with visibilities stored with real and complex values stored separately as `Float64`."
function DataVisReal(dvis::DataVis)
    return DataVisReal(dvis.lam, dvis.uu, dvis.vv, real(dvis.VV), imag(dvis.VV), dvis.invsig)
end

"
    conj!(dv::DataVis)

Import the complex conjugate function from Base, and extend it to work
on a DataVis.
I think this is necessary because the SMA and ALMA baseline conventions are swapped from
what I'm using in the NRAO Synthesis Summer School textbook."
function conj!(dv::DataVis)
    conj!(dv.VV)
end

"
    conj!(dvarr::Array{DataVis, 1})

Apply ``conj!`` to an entire DataVis array."
function conj!(dvarr::Array{DataVis, 1})
    for dset in dvarr
        conj!(dset) # Swap UV convention
    end
end

"
    get_qq(dv::DataVis)

Given a DataVis, convert the Cartesian spatial frequency coordinates ``(u,v)`` into a radial spatial coordinate ``q``."
function get_qq(dv::DataVis)
  sqrt(dv.uu.^2 .+ dv.vv.^2)
end

"
    write(dv::DataVis, fname::AbstractString)

Take in a visibility data set and then write it to the HDF5 file
The HDF5 file actually expects multi-channel data, so instead we will need to
store all of this information with arrays of shape (..., 1) [an extra trailing]
dimension of 1."
function write(dv::DataVis, fname::AbstractString)
    nvis = length(dv.uu)
    VV = reshape(dv.VV, (nvis, 1))
    fid = h5open(fname, "w")

    fid["freqs"] = [cc/(1e-4 * dv.lam)] # expects 1D array [Hz]
    fid["uu"] = reshape(dv.uu, (nvis, 1)) # expects 2D array
    fid["vv"] = reshape(dv.vv, (nvis, 1)) # expects 2D array
    fid["real"] = real(VV)
    fid["imag"] = imag(VV)
    fid["weight"] = reshape(dv.invsig, (nvis, 1)).^2 # expects 2D array

    close(fid)
end

"
    write(dvarr::Array{DataVis, 1}, fname::AbstractString)

Write an array of DataVis to a single HDF5 file."
function write(dvarr::Array{DataVis, 1}, fname::AbstractString)
    nvis = length(dvarr[1].VV)
    nlam = length(dvarr)
    fid = h5open(fname, "w")
    # hcat here stacks the individual channel data sets into a big block
    # of shape (nvis, nlam)

    fid["freqs"] = Float64[cc / (1e-4 * dv.lam) for dv in dvarr] # [Hz]
    fid["uu"] = hcat([dv.uu for dv in dvarr]...)
    fid["vv"] = hcat([dv.vv for dv in dvarr]...)
    fid["real"] = hcat([real(dv.VV) for dv in dvarr]...)
    fid["imag"] = hcat([imag(dv.VV) for dv in dvarr]...)
    fid["weight"] = hcat([dv.invsig.^2 for dv in dvarr]...)
    close(fid)

end

"
    copy_flags(source::AbstractString, dest::AbstractString)

Copy the flags from one dataset to another."
function copy_flags(source::AbstractString, dest::AbstractString)
    println("Copying flags from $source to $dest")

    fid_source = h5open(source, "r") # read only
    fid_dest = h5open(dest, "r+") # append mode

    @assert isapprox(read(fid_source, "uu"), read(fid_dest, "uu")) "UU spacings between datasets do not match, make sure you have loaded $source with flagged=true."
    @assert isapprox(read(fid_source, "vv"), read(fid_dest, "vv")) "VV spacings between dataset do not match, make sure you have loaded $source with flagged=true."
    @assert isapprox(read(fid_source, "weight"), read(fid_dest, "weight")) "Weights between dataset do not match, make sure you have loaded $source with flagged=true."

    fid_dest["flag"] = read(fid_source, "flag")

    close(fid_source)
    close(fid_dest)
end

"
    max_baseline(dvarr)

Determine the maximum uu or vv baseline contained in the dataset, so we know at what resolution we will need to synthesize the images.

returned in kilolambda.
"
function max_baseline(dvarr::Array{DataVis, 1})
    max = 0.0

    for dv in dvarr
        max_uu = maximum(dv.uu)
        if max_uu > max
            max = max_uu
        end

        max_vv = maximum(dv.vv)
        if max_vv > max
            max = max_vv
        end
    end

    return max
end

"
    get_nyquist_pixel(max_base::Float64, angular_width::Float64)

Determine how many pixels the image needs at a given distance in order to satisfy the Nyquist sampling theorem.

`max_base` is in kilolambda.
`angular_width` is in radians."
function get_nyquist_pixel(max_base::Float64, angular_width::Float64)

    nyquist_factor = 4 # normally 2, but we want to be extra sure.
    #Calculate the maximum dRA and dDEC from max_base
    dRA_max = 1/(nyquist_factor * max_base * 1e3) # [radians]

    npix = 64
    dRA = angular_width/npix

    while dRA > dRA_max
        npix *= 2
        dRA = angular_width/npix
    end

    return npix
end

"
    RawModelVis

Type to store the product of FFT'ing a single channel of a SkyImage."
mutable struct RawModelVis
    lam::Float64 # [μm] Wavelength (in microns) corresponding to channel
    uu::Vector{Float64} # [kλ] Vectors of the u, v locations
    vv::Vector{Float64} # [kλ]
    VV::Matrix{ComplexF64} # Output from rfft
end

# Taking the complex conjugate to make a full grid over all visibility space.
mutable struct FullModelVis
    lam::Float64 # [μm] Wavelength (in microns) corresponding to channel
    uu::Vector{Float64} # [kλ] Vectors of the u, v locations
    vv::Vector{Float64} # [kλ]
    VV::Matrix{ComplexF64} # Output from rfft
end

"
    -(vis1::FullModelVis, vis2::FullModelVis)

Subtract two `FullModelVis`'s.
"
function -(vis1::FullModelVis, vis2::FullModelVis)
    # @assert vis1.lam == vis2.lam "Visibilities must have the same wavelengths."
    @assert vis1.uu == vis2.uu "Visibilities must have same uu sampling."
    @assert vis1.vv == vis2.vv "Visibilities must have same vv sampling."
    VV = vis1.VV - vis2.VV
    return FullModelVis(vis1.lam, vis1.uu, vis1.vv, VV)
end

# Produced by gridding a RawModelVis to match the data
mutable struct ModelVis
    dvis::DataVis # A reference to the matching dataset
    VV::Vector{ComplexF64} # vector of complex visibilities directly
    # corresponding to the u, v locations in DataVis
end

"
    ModelVis(dvis::DataVis, fmvis::FullModelVis)

Given a DataSet and a FullModelVis, go through and interpolate at the
u,v locations of the DataVis."
function ModelVis(dvis::DataVis, fmvis::FullModelVis)
    nvis = length(dvis.VV)
    VV = Array{ComplexF64}(undef, nvis)
    for i=1:nvis
        VV[i] = interpolate_uv(dvis.uu[i], dvis.vv[i], fmvis)
    end

    return ModelVis(dvis, VV)
end

"
    ModelVisRotate(dvis::DataVis, fmvis::FullModelVis, PA::Real)

Given a DataSet, FullModelVis, and a desired position angle rotation, use the Fourier rotation theorem to sample the image at the rotated baselines and produce a rotated (sampled) model. A positive `PA` value (in degrees) means that the new model is rotated ``PA`` number of degrees counter-clockwise in the image plane (towards East, from North)."
function ModelVisRotate(dvis::DataVis, fmvis::FullModelVis, PA::Real)

    # Convert from degrees to radians
    PA = PA * deg

    nvis = length(dvis.VV)

    VV = Array{ComplexF64}(undef, nvis)
    for i=1:nvis
        # Rotate points according to PA
        uuprime = cos(PA) * dvis.uu[i] + sin(PA) * dvis.vv[i]
        vvprime = -sin(PA) * dvis.uu[i] + cos(PA) * dvis.vv[i]

        # Interpolate at rotated points
        VV[i] = interpolate_uv(uuprime, vvprime, fmvis)
    end

    return ModelVis(dvis, VV)
end

"
    ModelVis2DataVis(mvis::ModelVis)

Collapse a ModelVis into a DataVis in order to write to disk."
function ModelVis2DataVis(mvis::ModelVis)
    DV = mvis.dvis
    return DataVis(DV.lam, DV.uu, DV.vv, mvis.VV, DV.invsig)
end

function conj!(mv::ModelVis)
    conj!(mv.VV)
end

"
    ResidVis(dvarr::Array{DataVis, 1}, mvarr)

Return a model visibility file that contains the residuals."
function ResidVis(dvarr::Array{DataVis, 1}, mvarr)
    nchan = length(dvarr)
    @assert length(mvarr) == nchan # make sure data and model are the same length.
    rvarr = Array{DataVis}(undef, nchan)
    for i=1:nchan
        dvis = dvarr[i]
        mvis = mvarr[i]
        VV = dvis.VV - mvis.VV
        rvarr[i] = DataVis(dvis.lam, dvis.uu, dvis.vv, VV, dvis.invsig)
    end
    return rvarr
end


mutable struct ModelVisReal
    dvis::DataVis # A reference to the matching dataset
    VVr::Vector{Float64} # vector of complex visibilities directly
    VVi::Vector{Float64} # vector of complex visibilities directly
    # corresponding to the u, v locations in DataVis
end

"
    ModelVisReal(mvis::ModelVis)

Convert a `ModelVis` object, with visibilities stored as `ComplexF64`, into a `ModelVisReal`, with real and imaginary visibility values stored separately as `Float64`."
function ModelVisReal(mvis::ModelVis)
    return ModelVisReal(mvis.dvis, real(mvis.VV), imag(mvis.VV))
end

"
    lnprob(dvis::DataVis, mvis::DataVis)

Lnprob for computing the difference between two datasets. Does not check whether the baselines of the two datasets are the same. Use this function at your own discretion."
function lnprob(dvis::DataVis, mvis::DataVis)
    return -0.5 * sum(abs2, dvis.invsig .* (dvis.VV - mvis.VV)) # Basic chi2
end

"
    lnprob(dvis::DataVis, mvis::ModelVis)

Likelihood function used to calculate the likelihood of the data visibilities given a current model of them. Asserts that the model visibilities were sampled at the same baselines as the data."
function lnprob(dvis::DataVis, mvis::ModelVis)
    @assert dvis === mvis.dvis # Using the wrong ModelVis, otherwise!
    return -0.5 * sum(abs2, dvis.invsig .* (dvis.VV - mvis.VV)) # Basic chi2
end

"
    lnprob(dvis::DataVis, mvis::ModelVis, mu_RA::Real, mu_DEC::Real)

Calculate the lnlikelihood using a single loop to do the phase shifting and summation. No extra memory needed for re-copying the model.
"
function lnprob(dvis::DataVis, mvis::ModelVis, mu_RA::Float64, mu_DEC::Float64)
    lnp::Float64 = 0.0
    prefac::ComplexF64 = -2pi * 1.0im * 1e3 * arcsec
    for i=1:length(dvis.VV)
        lnp += abs2(dvis.invsig[i] * (dvis.VV[i] - mvis.VV[i] * exp(prefac * (dvis.uu[i] * mu_RA + dvis.vv[i] * mu_DEC))))
    end
    return -0.5 * lnp
end

"
    lnprob_real(dvis::DataVis, mvis::ModelVis, mu_RA::Float64, mu_DEC::Float64)

Calculate the lnlikelihood using a single loop to do the phase shifting and summation. No extra memory needed for re-copying the model. Uses Euler formula to do the calculation, rather than exponentiating a complex number.
"
function lnprob_real(dreal::DataVisReal, mreal::ModelVisReal, mu)
    lnp::Float64 = 0.0
    val::Float64 = 0.0
    sval::Float64 = 0.0 # sin of val
    cval::Float64 = 0.0 # cos of val

    # dreal and mreal both have the real and imaginary terms stored separately as Float64 vectors.
    prefac::Float64 = -2pi * 1e3 * arcsec

    for i=1:length(dreal.VVr)
        val = prefac * (mu[1] * dreal.uu[i] + mu[2] * dreal.vv[i])
        sval = sin(val)
        cval = cos(val)

        lnp += ((dreal.VVr[i] - cval * mreal.VVr[i] + sval * mreal.VVi[i]) * dreal.invsig[i])^2 + ((dreal.VVi[i] - sval * mreal.VVr[i] - cval * mreal.VVi[i]) * dreal.invsig[i])^2
    end
    return -0.5 * lnp
end


function lnprob_iter(dvis::DataVis, mvis::ModelVis)
    lnp = 0.0
    for i=1:length(dvis.VV)
        lnp += abs2(dvis.invsig[i] * (dvis.VV[i] - mvis.VV[i]))
    end
    return -0.5 * lnp
end


"
    chi2(dvis::DataVis, mvis::ModelVis)

Calculate the ``\\chi^2`` value of the current model. This function is not used in inference (see `lnprob`), but it is sometimes a useful number to report."
function chi2(dvis::DataVis, mvis::ModelVis)
    @assert dvis === mvis.dvis # Using the wrong ModelVis, otherwise!
    return sum(abs2, dvis.invsig .* (dvis.VV - mvis.VV)) # Basic chi2
end

"
    phase_shift!(mvis::ModelVis, mu_RA, mu_DEC)

Given a new model centroid in the image plane (in arcseconds), shift the model
visibilities by corresponding amount."
function phase_shift!(mvis::ModelVis, mu_RA, mu_DEC)

    mu = Float64[mu_RA, mu_DEC] * arcsec # [radians]

    nvis = length(mvis.VV)
    # Go through each visibility and apply the phase shift
    for i=1:nvis
        # Convert from [kλ] to [λ]
        R = Float64[mvis.dvis.uu[i], mvis.dvis.vv[i]] * 1e3 #[λ]
        # Not actually in polar phase form
        shift = exp(-2pi * 1.0im * (R' * mu)[1])
        mvis.VV[i] = mvis.VV[i] * shift
    end
end

function phase_shift(mvis::ModelVis, mu_RA, mu_DEC)
    mv = deepcopy(mvis)
    mv.dvis = mvis.dvis
    phase_shift!(mv, mu_RA, mu_DEC)
    return mv
end

# Shift the (0,0) data phase center to these new RA, DEC centroids. The mu_RA, mu_DEC should be
# the same sign as those in the modelVis phase shift routine.
# This is only necessary if the data is so off centered relative to the RADMC image size
# that you are running into problems with tapering cutting off a majority of the image.
# function phase_shift!(dvis::DataVis, mu_RA, mu_DEC)
#
#     mu = -Float64[mu_RA, mu_DEC] * arcsec # [radians]
#
#     nvis = length(dvis.VV)
#     # Go through each visibility and apply the phase shift
#     for i=1:nvis
#         # Convert from [kλ] to [λ]
#         R = Float64[dvis.dvis.uu[i], dvis.dvis.vv[i]] * 1e3 #[λ]
#         # Not actually in polar phase form
#         shift = exp(-2pi * 1.0im * (R' * mu)[1])
#         dvis.VV[i] = dvis.VV[i] * shift
#     end
# end

"
    phase_shift!(fvis::FullModelVis, mu_RA, mu_DEC)

Given a new model centroid in the image plane (in arcseconds), shift the model
visibilities by corresponding amount."
function phase_shift!(fvis::FullModelVis, mu_RA, mu_DEC)

    mu = Float64[mu_RA, mu_DEC] * arcsec # [radians]

    nu = length(fvis.uu)
    nv = length(fvis.vv)
    # Go through each visibility and apply the phase shift
    for i=1:nu
        for j=1:nv
        # Convert from [kλ] to [λ]
        R = Float64[fvis.uu[i], fvis.vv[j]] * 1e3 #[λ]
        # Not actually in polar phase form
        shift = exp(-2pi * 1.0im * (R' * mu)[1])
        fvis.VV[j,i] = fvis.VV[j,i] * shift
        end
    end
end

"
    transform(img::SkyImage, index::Int=1)

Transform the SkyImage produced by RADMC using FFT."
function transform(img::SkyImage, index::Int=1)

    # By default, select the first channel of any spectral hypercube. This
    # routine can only transform one channel at a time.
    data = img.data[:, :, index]

    lam = img.lams[index]

    # convert ra and dec in [arcsec] to radians, and then take the sin to
    # convert to ll, mm
    # ll = sin(img.ra * arcsec)
    # mm = sin(img.dec * arcsec)

    # Sin is not necessary, since we are using the small angle approximation
    ll = img.ra * arcsec
    mm = img.dec * arcsec

    # number of elements in each array
    nl = length(ll)
    nm = length(mm)

    # find the spacing between the elements
    dl = abs(ll[2] - ll[1]) # [radians]
    dm = abs(mm[2] - mm[1]) # [radians]

    # println("Transform using dl ", dl)

    # determine uv plane coordinates in kλ
    uu = fftshift(fftfreq(nl, dl)) * 1e-3 # [kλ]
    vv = fftshift(fftfreq(nm, dm)) * 1e-3 # [kλ]

    # properly pack the data for input using fftshift to move the component at
    # RA=0,DEC=0 to the first array element: data[1,1]
    # For RADMC images, we'll be moving the RA=~0.5, DEC=~0.5 element
    # (1/2 pixel away), and so there will need to be a corresponding phase
    # shift to correct.

    # We also want to normalize the result by the input array spacings, so that
    # they are directly comparable with the analytic transforms
    # (Numerical Recipes ed. 3, Press, Eqn 12.1.6)
    out = dl * dm * fftshift(fft(fftshift(data)))
    return FullModelVis(lam, uu, vv, out)
end


"
    fillModelVis(vis::RawModelVis)

This function is designed to copy the partial arrays in RawModelVis into a
full image for easy plotting. This means that the u axis can remain the same
but we'll need to make the complex conjugate of the v axis."
function fillModelVis(vis::RawModelVis)

    # The full image will stretch from -(n/2 - 1) to (n/2 -1) and
    # have one less point in the u,v directions

    #Swap the +u and -u quadrants so that we have u increasing from -n/2 to 0 to n/2 - 1
    VV = fftshift(vis.VV, 2)

    # Create the top chunk (v >= 0), by flipping the array such that
    # u=n/2 is in the first row of the array and u=0 on the last row.
    # then trim off the first row and column, so that the array stretches
    # 0 <= v <= n/2 - 1  and -(n/2 - 1) <= u <= (n/2 - 1)
    top = flipdim(VV, 1)[2:end, 2:end]

    # Create the bottom chunk (v < 0) by using the original orientation,
    # and trimming off rows v = 0 and v = -n/2, and column u = n/2.
    # we use the property that the FFT of a real image is Hermitian, and thus use
    # the complex conjugate property to fill in the missing quadrants.
    bottom = conj(flipdim(VV, 2))[2:end-1, 1:end-1]

    # make the coordinates
    uu = fftshift(vis.uu)[2:end] # clip off u=-n/2 frequency

    vv_top = flipdim(vis.vv, 1)[2:end] # clip off v=n/2 frequency
    vv_bottom = -1 .* vis.vv[2:end-1] # clip off v=0 and v=-n/2 frequency
    vv = vcat(vv_top, vv_bottom)

    return FullModelVis(vis.lam, uu, vv, vcat(top, bottom))
end

"
    plan_interpolate(dvis::DataVis, uu::Vector{Float64}, vv::Vector{Float64})

Return a function that is used to interpolate the visibilities, in the
spirit of interpolate_uv but *much* faster.
Closures save time and money!

The point is that if the distance to the source and the size of the sythesized image
are not changing, then we will always be interpolating from a dense Visibility grid that
has the exact same U and V spacings, meaning that the weighting terms used to evaluate the
interpolation for a specific visibility can be cached.

So, using a closure, those weighting terms are calculated once and then cached for further use."
function plan_interpolate(dvis::DataVis, uu::Vector{Float64}, vv::Vector{Float64})

    nvis = length(dvis.VV)
    uinds = Array{UnitRange{Int64}}(undef, nvis)
    vinds = Array{UnitRange{Int64}}(undef, nvis)
    uws = Array{Float64}(undef, 6, nvis) #stored along columns
    vws = Array{Float64}(undef, 6, nvis)

    for i=1:nvis
        u = dvis.uu[i]
        v = dvis.vv[i]
        iu0 = argmin(abs.(u .- uu))
        iv0 = argmin(abs.(v .- vv))

        # now find the relative distance to this nearest grid point (not absolute)
        u0 = u - uu[iu0]
        v0 = v - vv[iv0]

        # determine the uu and vv distance for 3 grid points (could be later taken out)
        du = abs(uu[4] - uu[1])
        dv = abs(vv[4] - vv[1])

        # 2. Calculate the appropriate u and v indexes for the 6 nearest pixels
        # (3 on either side)

        # Are u0 and v0 to the left or the right of the index?
        # we want to index three to the left, three to the right

        # First check that we are still in bounds of the array
        # Check to make sure that at least three grid points exist in all directions
        lenu = length(uu)
        lenv = length(vv)
        @assert iu0 >= 4
        @assert iv0 >= 4
        @assert lenu - iu0 >= 4
        @assert lenv - iv0 >= 4

        if u0 >= 0.0
            # To the right of the index
            uind = iu0-2:iu0+3
        else
            # To the left of the index
            uind = iu0-3:iu0+2
        end

        if v0 >= 0.0
            # To the right of the index
            vind = iv0-2:iv0+3
        else
            # To the left of the index
            vind = iv0-3:iv0+2
        end

        # Store these slices so we can later index into the visibility array
        uinds[i] = uind
        vinds[i] = vind

        # 3. Calculate the weights corresponding to these 6 nearest pixels (gcffun)
        etau = (uu[uind] .- u)/du
        etav = (vv[vind] .- v)/dv

        uw = gcffun(etau)
        vw = gcffun(etav)

        # 4. Normalization such that it has an area of 1. Divide by w later.
        uw = uw/sum(uw)
        vw = vw/sum(vw)

        uws[:,i] = uw
        vws[:,i] = vw
    end

    tol = 1e-5 * ones(Float64, sizeof(uu))

    # This function inherits all of the variables just defined in this scope (uu, vv)
    function interpolate(data::DataVis, fmvis::FullModelVis)
        # Assert that we calculated the same UU and VV spacings for the FT'ed image, otherwise we did something wrong!
        # The 1e-5 addition is to prevent an undetermined error from the uu = 0.0 point.
        @assert all(abs.((uu .- fmvis.uu) ./ (uu .+ 1e-5)) .< tol)
        @assert all(abs.((vv .- fmvis.vv) ./ (vv .+ 1e-5)) .< tol)

        # output array
        Vmodel = Array{ComplexF64}(undef, nvis)

        for i=1:nvis

            Vdata = fmvis.VV[vinds[i], uinds[i]] # Array is packed like the image

            # 5. Loop over all 36 grid indices and sum to find the interpolation.
            cum::ComplexF64 = 0.0 + 0.0im
            for k=1:6
                for l=1:6
                    cum += uws[k,i] * vws[l,i] * Vdata[l,k] # Array is packed like the image
                end
            end

            Vmodel[i] = cum

        end

        return ModelVis(data, Vmodel)
    end

    return interpolate
end

# called ModGrid in gridding.c (KR code) and in Model.for (MIRIAD)
# Uses spheroidal wave functions to interpolate a model to a (u,v) coordinate.
# u,v are in [kλ]

"
    interpolate_uv(u::Float64, v::Float64, vis::FullModelVis)

Interpolates a dense grid of visibilities (e.g., from FFT of an image) to a specfic (u,v) point using spheroidal functions in a band-limited manner designed to reduce aliasing.
"
function interpolate_uv(u::Float64, v::Float64, vis::FullModelVis)

    # Note that vis.uu goes from positive to negative (East-West)
    # and vis.vv goes from negative to positive (North-South)

    # 1. Find the nearest gridpoint in the FFT'd image.
    iu0 = argmin(abs.(u .- vis.uu))
    iv0 = argmin(abs.(v .- vis.vv))

    # now find the relative distance from (u,v) to this nearest grid point (not absolute)
    u0 = u - vis.uu[iu0]
    v0 = v - vis.vv[iv0]

    # determine the uu and vv distance for 3 grid points (could be later taken out)
    du = abs.(vis.uu[4] - vis.uu[1])
    dv = abs.(vis.vv[4] - vis.vv[1])

    # 2. Calculate the appropriate u and v indexes for the 6 nearest pixels
    # (3 on either side)

    # Are u0 and v0 to the left or the right of the index?
    # we want to index three to the left, three to the right

    # First check that our (u,v) point still exists within the appropriate margins of the
    # dense visibility array
    # This is to make sure that at least three grid points exist in all directions
    # If this fails, this means that the synthesized image is too large compared to the sampled visibilities (meaning that the dense FFT grid is too small).
    # The max u,v sampled is 2/dRA or 2/dDec. This means that if dRA or dDEC is too large, then this
    # will fail
    lenu = length(vis.uu)
    lenv = length(vis.vv)
    @assert iu0 >= 4
    @assert iv0 >= 4
    @assert lenu - iu0 >= 4
    @assert lenv - iv0 >= 4

    if u0 >= 0.0
        # To the right of the index
        uind = iu0-2:iu0+3
    else
        # To the left of the index
        uind = iu0-3:iu0+2
    end

    if v0 >= 0.0
        # To the right of the index
        vind = iv0-2:iv0+3
    else
        # To the left of the index
        vind = iv0-3:iv0+2
    end

    etau = (vis.uu[uind] .- u)/du
    etav = (vis.vv[vind] .- v)/dv
    VV = vis.VV[vind, uind] # Array is packed like the image

    # 3. Calculate the weights corresponding to these 6 nearest pixels (gcffun)
    uw = gcffun(etau)
    vw = gcffun(etav)

    # 4. Normalization such that it has an area of 1. Divide by w later.
    w = sum(uw) * sum(vw)

    # 5. Loop over all 36 grid indices and sum to find the interpolation.
    cumulative::ComplexF64 = 0.0 + 0.0im
    for i=1:6
        for j=1:6
            cumulative += uw[i] * vw[j] * VV[j,i] # Array is packed like the image
        end
    end

    cumulative = cumulative/w

    return cumulative
end

"
    rfftfreq(n::Int, d::Float64)

Return the frequencies corresponding to the output of the real FFT. After numpy.fft.rfftfreq
    f = [0, 1, ...,     n/2-1,     n/2] / (d*n)   if n is even
    f = [0, 1, ..., (n-1)/2-1, (n-1)/2] / (d*n)   if n is odd."
function rfftfreq(n::Int, d::Float64)
    # n even
    if n % 2 == 0
        #Array contains n/2 + 1 elements from [0, n/2]
        return linspace(0, n/2, int(n/2) + 1) / (d * n)

    # n odd
    else
        #Array contains n/2 elements from [0, (n-1)/2]
        return linspace(0, (n - 1)/2, int(n/2)) / (d * n)
    end
end

"
    fftfreq(n::Int, d::Float64)

After numpy.fft.fftfreq

    f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
    f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd."
function fftfreq(n::Int, d::Float64)
    val = 1.0/(n * d)
    results = Array{Float64}(undef, n)
    N = floor(Int, (n  - 1)/2) + 1

    p1 = Float64[i for i=0:(N-1)]
    results[1:N] = p1

    p2 = Float64[i for i=-(floor(Int, n/2)):-1]
    results[N+1:end] = p2

    return results * val
end


end



# Transform the SkyImage produced by RADMC into a RawModelVis object using rfft
#function transform(img::SkyImage, index::Int=1)
#
#    # By default, select the first channel of any spectral hypercube.
#    data = img.data[:, :, index]
#
#    lam = img.lams[index]
#
#    # convert ra and dec in [arcsec] to radians, and then take the sin to convert to ll, mm
#    ll = sin(img.ra * arcsec)
#    mm = sin(img.dec * arcsec)
#
#    # number of elements in each array
#    nl = length(ll)
#    nm = length(mm)
#
#    # find the spacing between the elements
#    dl = ll[2] - ll[1] # [radians]
#    dm = mm[2] - mm[1] # [radians]
#
#    # determine uv plane coordinates in kλ
#    uu = fftfreq(nl, dl) * 1e-3 # [kλ]
#    vv = rfftfreq(nm, dm) * 1e-3 # [kλ]
#
#    # properly pack the data for input using fftshift to move the 0,0 component to the corner.
#
#    # From rfft: If A has size (n_1, ..., n_d), the result has size (floor(n_1/2)+1, ..., n_d).
#    # This means that the `v` dimension is halved but the `u` dimension remains the same.
#    # Even though u has real symmetry, we still get the conjugate values. This transpose of an axes is due
#    # to the way images are stored. Even though we have u corresponding to the x axis and v corresponding to the
#    # y axis, the image is actually stored on disk as img[v, u], with v changing fastest.
#
#    # We also want to normalize the result by the input array spacings, so that they are directly comparable with
#    # the analytic transforms (Numerical Recipes ed. 3, Press, Eqn 12.1.6)
#    out = dl * dm * rfft(fftshift(data))
#    return RawModelVis(lam, uu, vv, out)
#end
#
