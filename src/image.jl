# Read the image format written by RADMC3D.

# Read the ascii text in `image.out` and parse things into a 3 dimensional matrix (x, y, lambda)

# The first four lines are format information
# iformat # = 1 (2 is local observer)
# im_nx   im_ny #number of pixels in x and y directions
# nlam # number of images at different wavelengths
# pixsize_x  pixsize_y # size of the pixels in cm
# lambda[1] ... lambda[nlam + 1] # wavelengths (um) correspending to images
# pixels, ordered from left to right (increasing x) in the inner loop, and from bottom to top (increasing y) in the outer loop. And wavelength is the outermost loop.

"The image module contains various data types for reading and holding images produced by the radiative transfer programs (via `RADMC-3D`), as well as routines for processing these images."
module image

export imread, imToSky, imToSpec, SkyImage, -
export taureadImg, taureadPos

using ..constants
# import Images # The Images.jl package, not affiliated w/ DiskJockey
# using Dierckx

import Base.- # extend this for Image

# Define an image type, which can store the data as well as pixel spacing

abstract type Image end

"
    RawImage(data, pixsize_x, pixsize_y, lams)

Hold the raw output from `RADMC-3D` in a 3D array (npix_y, npix_x, nlam).
RawImage reflects the RADMC convention that both x and y are increasing
with array index. This means that to display the image as RADMC intends it,
you must set the first array element to the lower left corner."
mutable struct RawImage <: Image
    data::Array{Float64,3} # [ergs/s/cm^2/Hz/ster]
    pixsize_x::Float64 # [cm]
    pixsize_y::Float64 # [cm]
    lams::Vector{Float64} # [μm]
end

"SkyImage is a holder that has both RA and DEC increasing with array index
This convention is necessary for the FFT step
However, to display this image in the traditional sky convention (North up,
East to the left), you must set the first array element to the lower left
corner *and* flip the array along the RA axis: `fliplr(data)` or flipdim(data, 2)"
mutable struct SkyImage <: Image
    data::Array{Float64,3} # [Jy/pixel]
    ra::Vector{Float64} # [arcsec]
    dec::Vector{Float64} # [arcsec]
    lams::Vector{Float64} # [μm]
end

"TausurfImage is designed to hold the results of the RADMC-3D `tausurf` operation. This is
The distance in cm above or below the plane tangent to the observer, which intersects the origin of the model.

From the RADMC3D manual:
The image output file image.out will now contain, for each pixel, the position along the ray in centimeters where τ = τs. The zero point is the surface perpendicular to the direction of observation, going through the pointing position (which is, by default the origin (0, 0, 0)). Positive values mean that the surface is closer to the observer than the plane, while negative values mean that the surface is behind the plane.
So for this datastructure, it's the same thing as RawImage, just instead of intensity, we have distance above/below plane."
mutable struct TausurfImg <: Image
    data::Array{Float64,3} # cm above/behind central projected plane of disk
    pixsize_x::Float64 # [cm]
    pixsize_y::Float64 # [cm]
    lams::Vector{Float64} # [μm]
end

"Encapsulates the 3D position of the pixels representing the tau=1 surface, in the same datashape as the image.
For each pixel, this is the x, y, or z position."
mutable struct TausurfPos
    data_x::Array{Float64,3} # [cm]
    data_y::Array{Float64,3} # [cm]
    data_z::Array{Float64,3} # [cm]
    lams::Vector{Float64} # [μm]
end

"Subtraction for `SkyImage`s"
function -(img1::SkyImage, img2::SkyImage)
    # @assert img1.lams == img2.lams "Images must have the same wavelengths."
    @assert isapprox(img1.ra, img2.ra) "Images must have same RA coordinates."
    @assert isapprox(img1.dec, img2.dec) "Images must have same DEC coordinates."
    data = img1.data - img2.data
    return SkyImage(data, img1.ra, img1.dec, img1.lams)
end

"SkyImage constructor for just a single frame"
SkyImage(data::Matrix{Float64}, ra::Vector{Float64}, dec::Vector{Float64}, lam::Float64) =
SkyImage(reshape(data, tuple(size(data)..., 1)), ra, dec, [lam])

"
    imread(file=\"image.out\")

Read the image file (default=image.out) and return it as an Image object, which contains the fluxes in Jy/pixel,
the sizes and locations of the pixels in arcseconds, and the wavelengths (in microns) corresponding to the images"
function imread(file = "image.out")

    fim = open(file, "r")
    iformat = parse(Int, readline(fim))
    im_nx, im_ny = split(readline(fim))
    im_nx = parse(Int, im_nx)
    im_ny = parse(Int, im_ny)
    nlam = parse(Int, readline(fim))
    pixsize_x, pixsize_y = split(readline(fim))
    pixsize_x = parse(Float64, pixsize_x)
    pixsize_y = parse(Float64, pixsize_y)

    # Read the wavelength array
    lams = Array{Float64}(undef, nlam)
    for i = 1:nlam
        lams[i] = parse(Float64, readline(fim))
    end

    # Create an array with the proper size, and then read the file into it
    data = Array{Float64}(undef, im_ny, im_nx, nlam)

    # According to the RADMC manual, section A.15, the pixels are ordered
    # left to right (increasing x) in the inner loop, and from bottom to top
    # (increasing y) in the outer loop.
    # Basically, pack the array in order but display with origin=lower.

    # Because of the way an image is stored as a matrix, we actually pack the
    # array indices as data[y, x, lam]
    # radmc3dPy achieves something similar by keeping indices the x,y but
    # swaping loop order (radmcPy/image.py:line 675)
    for k = 1:nlam
        readline(fim) # Junk space
        for j = 1:im_ny
            for i = 1:im_nx
                data[j,i,k] = parse(Float64, readline(fim))
            end
        end
    end

    close(fim)

    # According to the RADMC3D manual, the units are *intensity* [erg cm−2 s−1 Hz−1 ster−1]

    return RawImage(data, pixsize_x, pixsize_y, lams)
end

"
    taureadImg(file=\"image_tausurf.out\")

Like imread, but for tausurf. Pixels that have no ``\\tau`` surface are set to `NaN`."
function taureadImg(file = "image_tausurf.out")
    fim = open(file, "r")
    iformat = parse(Int, readline(fim))
    im_nx, im_ny = split(readline(fim))
    im_nx = parse(Int, im_nx)
    im_ny = parse(Int, im_ny)
    nlam = parse(Int, readline(fim))
    pixsize_x, pixsize_y = split(readline(fim))
    pixsize_x = parse(Float64, pixsize_x)
    pixsize_y = parse(Float64, pixsize_y)
    
    # Read the wavelength array
    lams = Array{Float64}(nlam)
    for i = 1:nlam
        lams[i] = parse(Float64, readline(fim))
    end

    # Create an array with the proper size, and then read the file into it
    data = Array{Float64}(im_ny, im_nx, nlam)

    for k = 1:nlam
        readline(fim) # Junk space
        for j = 1:im_ny
            for i = 1:im_nx
                val = parse(Float64, readline(fim))

                data[j,i,k] = val

                # If RADMC3D signaled that there is no tau=1 surface here, set height to NaN
                if isapprox(val, -1e91)
                    data[j,i,k] = NaN #
                else
                    data[j,i,k] = val
                end
            end
        end
    end

    close(fim)

    # According to the RADMC3D manual, the units are *intensity* [erg cm−2 s−1 Hz−1 ster−1]

    return TausurfImg(data, pixsize_x, pixsize_y, lams)

end

"Read the (x,y,z) positions of the ``\\tau=1`` pixels."
function taureadPos(file = "tausurface_3d.out")
    fim = open(file, "r")
    iformat = parse(Int, readline(fim))
    im_nx, im_ny = split(readline(fim))
    im_nx = parse(Int, im_nx)
    im_ny = parse(Int, im_ny)
    nlam = parse(Int, readline(fim))
    # pixsize_x, pixsize_y = split(readline(fim))
    # pixsize_x = parse(Float64, pixsize_x)
    # pixsize_y = parse(Float64, pixsize_y)

    # Read the wavelength array
    lams = Array{Float64}(nlam)
    for i = 1:nlam
        lams[i] = parse(Float64, readline(fim))
    end

    # Create an array with the proper size, and then read the file into it
    data_x = Array{Float64}(im_ny, im_nx, nlam)
    data_y = Array{Float64}(im_ny, im_nx, nlam)
    data_z = Array{Float64}(im_ny, im_nx, nlam)

    # In contrast to the other image formats, apparently there is only a space before the lams, not inbetween lams.
    readline(fim) # Junk space

    for k = 1:nlam
        for j = 1:im_ny
            for i = 1:im_nx
                val_x, val_y, val_z = split(readline(fim))
                val_x = parse(Float64, val_x)
                val_y = parse(Float64, val_y)
                val_z = parse(Float64, val_z)

                # If RADMC3D signaled that there is no tau=1 surface here (any of the xyz values are -1e91)
                # then set the height to NaN
                if isapprox(val_x, -1e91)
                    data_x[j,i,k] = NaN
                    data_y[j,i,k] = NaN
                    data_z[j,i,k] = NaN
                else
                    data_x[j,i,k] = val_x
                    data_y[j,i,k] = val_y
                    data_z[j,i,k] = val_z
                end
            end
        end
    end

    close(fim)

    return TausurfPos(data_x, data_y, data_z, lams)

end

"
    imToSky(img::RawImage, dpc::Float64)

Convert a `RawImage` to a `SkyImage`. Assumes dpc is parsecs."
function imToSky(img::RawImage, dpc::Float64)

    # The RawImage is oriented with North up and East increasing to the left.
    # this means that for the RawImage, the delta RA array goes from + to -

    # However, the SkyImage actually requires RA (ll) in increasing form.
    # Therefore we flip along the RA axis, fliplr(data) or flipdim(data, 2)

    # println("Min and max intensity ", minimum(img.data), " ", maximum(img.data))
    # println("Pixel size ", img.pixsize_x)
    # println("Steradians subtended by each pixel ",  img.pixsize_x * img.pixsize_y / (dpc * pc)^2)

    # convert from ergs/s/cm^2/Hz/ster to to Jy/ster
    conv = 1e23 # [Jy/ster]

    # Conversion from erg/s/cm^2/Hz/ster to Jy/pixel at 1 pc distance.
    # conv = 1e23 * img.pixsize_x * img.pixsize_y / (dpc * pc)^2

    # Flip across RA dimension
    # dataJy = flipdim(img.data, 2) .* conv
    dataJy = reverse(img.data, dims = 2) .* conv

    (im_ny, im_nx) = size(dataJy)[1:2] # y and x dimensions of the image

    # The locations of pixel centers in cm
    # if n_x = 16, goes [-7.5, -6.5, ..., -0.5, 0.5, ..., 6.5, 7.5] * pixsize
    xx = ((Float64[i for i = 0:im_nx - 1] .+ 0.5) .- im_nx / 2.) * img.pixsize_x
    yy = ((Float64[i for i = 0:im_ny - 1] .+ 0.5) .- im_ny / 2.) * img.pixsize_y

    # The locations of the pixel centers in relative arcseconds
    # Note both RA and DEC increase with array index.
    ra = xx ./ (AU * dpc)
    dec = yy ./ (AU * dpc)

    return SkyImage(dataJy, ra, dec, img.lams)

end


"Take an image and integrate all the frames to create a spatially-integrated spectrum"
function imToSpec(img::SkyImage)

    # pixels in SkyImage are Jy/ster

    # convert from Jy/str to Jy/pixel using str/pixel
    dRA = abs(img.ra[2] - img.ra[1]) * arcsec
    dDEC = abs(img.dec[2] - img.dec[1]) * arcsec

    # Add up all the flux in the pixels to create the spectrum
    flux = dropdims(sum(img.data .* dRA .* dDEC, dims = (1, 2)), dims = (1, 2))
    spec = hcat(img.lams, flux) # First column is wl, second is flux

    return spec
end

"Calculate the integrated line flux"
function integrateSpec(spec::Matrix{Float64}, lam0::Float64)
    # First column is wl, second is flux

    wl = spec[:,1]
    fl = spec[:,2]

    # Convert wl to kms
    vs = c_kms * (spec[:,1] .- lam0) / lam0

    if vs[2] - vs[1] < 0
        reverse!(vs)
        reverse!(fl)
    end

    # println(vs)
    # println(fl)

    spl = Spline1D(vs, fl)
    tot = integrate(spl, vs[1], vs[end])

    return tot
end
    
"Storage for zeroth moment map"
mutable struct ZerothMoment <: Image
   data::Array{Float64,2} # [Jy · Hz / pixel]
   ra::Vector{Float64} # [arcsec]
   dec::Vector{Float64} # [arcsec]
end

"
    convert(::Type{ZerothMoment}, img::SkyImage)

Convert a `SkyImage` to a `ZerothMoment` map."
function convert(::Type{ZerothMoment}, img::SkyImage)
    # Sum along the frequency axis
    data = squeeze(sum(img.data, (3)), 3)
    return ZerothMoment(data, img.ra, img.dec)
end

end # model
