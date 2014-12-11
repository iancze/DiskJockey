# Provides FFT functionality, regridding,

module visibilities

using HDF5
using image
using constants

export transform, DataVis, RawModelVis, FullModelVis, fillModelVis, rfftfreq, fftfreq

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
    lams = cc ./ read(fid["freqs"]) * 1e4 # Convert from Hz to μm
    uu = read(fid["uu"])
    vv = read(fid["vv"])
    real = read(fid["real"])
    imag = read(fid["imag"])
    VV = real + imag .* im # Complex visibility
    invsig = sqrt(read(fid["weight"]))
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
    lam = cc / fid["freqs"][index][1] * 1e4 # Convert from Hz to μm
    uu = vec(fid["uu"][:, index])
    vv = vec(fid["vv"][:, index])
    real = vec(fid["real"][:, index])
    imag = vec(fid["imag"][:, index])
    VV = real + imag .* im # Complex visibility
    invsig = sqrt(vec(fid["weight"][:, index]))
    close(fid)

    #Return a DataVis object
    return DataVis(lam, uu, vv, VV, invsig)
end

# Take in a visibility data set and then write it to the HDF5 file
function write(dv:DataVis, fname::ASCIIString)
    fid = h5open(fname, "w")
    fid["freqs"] = reshape(dv.lam * 1e-4, tuple(size(data)..., 1))# expects 1D array
    fid["uu"] = # expects 2D array
    fid["vv"] = # expects 2D array
    fid["real"] = # expects 2D array
    fid["imag"] = # expects 2D array
    fid["weight"] = # expects 2D array

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
    dvis::DataVis #A reference to the matching dataset

end

function lnprob(dvis::DataVis, mvis::ModelVis)
    @assert dvis == mvis.dvis, "Using the wrong ModelVis!"

    return -0.5 * sumabs2(invsig *(dvis.VV - mvis.VV)) # Basic chi2
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
