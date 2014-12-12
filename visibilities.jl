# Provides FFT functionality, regridding,

module visibilities

using HDF5
using image
using gridding
using constants

export DataVis, ModelVis, RawModelVis, FullModelVis, fillModelVis, write
export interpolate_uv
export transform, rfftfreq, fftfreq
export lnprob

# Stores the data visibilities for a single channel
type DataVis
    lam::Float64 # [μm] Wavelength (in microns) corresponding to channel
    uu::Vector{Float64} # [kλ] Vectors of the u, v locations in kilolambda
    vv::Vector{Float64} # [kλ]
    VV::Vector{Complex128} # [Jy] Complex visibilities
    invsig::Vector{Float64} # 1/sigma [1/Jy]
end

# Read all the visibilities from an HDF5 string and return them as an array of DataVis objects
function DataVis(fname::ASCIIString)
    fid = h5open(fname, "r")
    lams = read(fid["lams"]) # [μm]
    uu = read(fid["uu"])
    vv = read(fid["vv"])
    real = read(fid["real"])
    imag = read(fid["imag"])
    VV = real + imag .* im # Complex visibility
    invsig = read(fid["invsig"])
    close(fid)

    nlam = length(lams)

    #Return an array of DataVis
    out = Array(DataVis, nlam)
    for i=1:nlam
        out[i] = DataVis(lams[i], uu[:,i], vv[:,i], VV[:, i], invsig[:, i])
    end
    return out
end

# Read just one channel of visibilities from the HDF5 file
function DataVis(fname::ASCIIString, index::Int)
    fid = h5open(fname, "r")
    # the indexing and `vec` are necessary here because HDF5 doesn't naturally
    # squeeze trailing dimensions of length-1
    lam = fid["lams"][index][1] # [μm]
    uu = vec(fid["uu"][:, index])
    vv = vec(fid["vv"][:, index])
    real = vec(fid["real"][:, index])
    imag = vec(fid["imag"][:, index])
    VV = real + imag .* im # Complex visibility
    invsig = vec(fid["invsig"][:, index])
    close(fid)

    #Return a DataVis object
    return DataVis(lam, uu, vv, VV, invsig)
end

# Take in a visibility data set and then write it to the HDF5 file
# The HDF5 file actually expects multi-channel data, so instead we will need to
# store all of this information with arrays of shape (..., 1) [an extra trailing]
# dimension of 1
function write(dv::DataVis, fname::ASCIIString)
    nvis = length(dv.uu)
    VV = reshape(dv.VV, (nvis, 1))
    fid = h5open(fname, "w")
    fid["lams"] = [dv.lam] # expects 1D array
    fid["uu"] = reshape(dv.uu, (nvis, 1)) # expects 2D array
    fid["vv"] = reshape(dv.vv, (nvis, 1)) # expects 2D array
    fid["real"] = real(VV)
    fid["imag"] = imag(VV)
    fid["invsig"] = reshape(dv.invsig, (nvis, 1)) # expects 2D array

    close(fid)
end

#TODO: These visibilites need to be renamed, since the are confusing with ModelVis
# Produced by FFT'ing a single channel of a SkyImage
type RawModelVis
    lam::Float64 # [μm] Wavelength (in microns) corresponding to channel
    uu::Vector{Float64} # [kλ] Vectors of the u, v locations
    vv::Vector{Float64} # [kλ]
    VV::Matrix{Complex128} # Output from rfft
end

# Taking the complex conjugate to make a full grid over all visibility space.
type FullModelVis
    lam::Float64 # [μm] Wavelength (in microns) corresponding to channel
    uu::Vector{Float64} # [kλ] Vectors of the u, v locations
    vv::Vector{Float64} # [kλ]
    VV::Matrix{Complex128} # Output from rfft
end

# Produced by gridding a RawModelVis to match the data
type ModelVis
    dvis::DataVis # A reference to the matching dataset
    VV::Vector{Complex128} # vector of complex visibilities directly
    # corresponding to the u, v locations in DataVis
end

# Given a DataSet and a FullModelVis, go through and interpolate at the
# u,v locations of the DataVis
function ModelVis(dvis::DataVis, fmvis::FullModelVis)
    nvis = length(dvis.VV)
    VV = Array(Complex128, nvis)
    for i=1:nvis
        VV[i] = interpolate_uv(dvis.uu[i], dvis.vv[i], fmvis)
    end

    return ModelVis(dvis, VV)
end

function lnprob(dvis::DataVis, mvis::ModelVis)
    @assert dvis == mvis.dvis # Using the wrong ModelVis, otherwise!

    return -0.5 * sumabs2(dvis.invsig .* (dvis.VV - mvis.VV)) # Basic chi2
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

# Transform the SkyImage produced by RADMC into a RawModelVis object using fft
function transform(img::SkyImage, index::Int=1)

    # By default, select the first channel of any spectral hypercube.
    data = img.data[:, :, index]

    lam = img.lams[index]

    # convert ra and dec in [arcsec] to radians, and then take the sin to convert to ll, mm
    ll = sin(img.ra * arcsec)
    mm = sin(img.dec * arcsec)

    # number of elements in each array
    nl = length(ll)
    nm = length(mm)

    # find the spacing between the elements
    dl = ll[2] - ll[1] # [radians]
    dm = mm[2] - mm[1] # [radians]

    # determine uv plane coordinates in kλ
    uu = fftshift(fftfreq(nl, dl)) * 1e-3 # [kλ]
    vv = fftshift(fftfreq(nm, dm)) * 1e-3 # [kλ]

    # properly pack the data for input using fftshift to move the 0,0 component to the corner.

    # We also want to normalize the result by the input array spacings, so that they are directly comparable with
    # the analytic transforms (Numerical Recipes ed. 3, Press, Eqn 12.1.6)
    out = dl * dm * fftshift(fft(fftshift(data)))
    return FullModelVis(lam, uu, vv, out)
end


# This function is designed to copy the partial arrays in RawModelVis into a full image for easy plotting
# This means that the u axis can remain the same but we'll need to make the complex conjugate of the v axis.
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

# called ModGrid in gridding.c (KR code) and in Model.for (MIRIAD)
# Uses spheroidal wave functions to interpolate a model to a (u,v) coordinate.
# u,v are in [kλ]
function interpolate_uv(u::Float64, v::Float64, vis::FullModelVis)

    # 1. Find the nearest gridpoint in the FFT'd image.
    iu0 = indmin(abs(u - vis.uu))
    iv0 = indmin(abs(v - vis.vv))

    # now find the relative distance to this nearest grid point (not absolute)
    u0 = u - vis.uu[iu0]
    v0 = v - vis.vv[iv0]

    # determine the uu and vv distance for 3 grid points (could be later taken out)
    du = abs(vis.uu[4] - vis.uu[1])
    dv = abs(vis.vv[4] - vis.vv[1])

    # 2. Calculate the appropriate u and v indexes for the 6 nearest pixels (3 on either side)

    # Are u0 and v0 to the left or the right of the index?
    # we want to index three to the left, three to the right

    # First check that we are still in bounds of the array
    # Check to make sure that at least three grid points exist in all directions
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

    # println("Sampling at uu: ", vis.uu[uind])
    # println("Sampling at vv: ", vis.uu[vind])

    etau = (vis.uu[uind] .- u)/du
    etav = (vis.vv[vind] .- v)/dv
    VV = vis.VV[vind, uind] # Array is packed like the image

    # println("etau: ", etau)
    # println("etav: ", etav)

    # 3. Calculate the weights corresponding to these 6 nearest pixels (gcffun)
    # TODO: Explore using something other than alpha=1.0
    uw = gcffun(etau, 1.0)
    vw = gcffun(etav, 1.0)

    # 4. Normalization such that it has an area of 1. Divide by w later.
    w = sum(uw) * sum(vw)

    # 5. Loop over all 36 grid indices and sum to find the interpolation.
    cum::Complex128 = 0.0 + 0.0im
    for i=1:6
        for j=1:6
            cum += uw[i] * vw[j] * VV[j,i] # Array is packed like the image
        end
    end

    cum = cum/w

    return cum
end

# Return the frequencies corresponding to the output of the real FFT. After numpy.fft.rfftfreq
# f = [0, 1, ...,     n/2-1,     n/2] / (d*n)   if n is even
# f = [0, 1, ..., (n-1)/2-1, (n-1)/2] / (d*n)   if n is odd
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

# After numpy.fft.fftfreq
# f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
# f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
function fftfreq(n::Int, d::Float64)
    val = 1./(n * d)
    results = Array(Float64, (n,))
    N = floor((n  - 1)/2) + 1

    p1 = Float64[i for i=0:(N-1)]
    results[1:N] = p1

    p2 = Float64[i for i=-(floor(n/2)):-1]
    results[N+1:end] = p2

    return results * val
end


end
