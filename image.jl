# Read the image format written by RADMC3D.

# Read the ascii text in `image.out` and parse things into a 3 dimensional matrix (x, y, lambda)

# The first four lines are format information
# iformat # = 1 (2 is local observer)
# im_nx   im_ny #number of pixels in x and y directions
# nlam # number of images at different wavelengths
# pixsize_x  pixsize_y # size of the pixels in cm
# lambda[1] ... lambda[nlam + 1] # wavelengths (um) correspending to images
# pixels, ordered from left to right (increasing x) in the inner loop, and from bottom to top (increasing y) in the outer loop. And wavelength is the outermost loop.

module image

export imread, imToSky, imToSpec, SkyImage

using constants

# Define an image type, which can store the data as well as pixel spacing

abstract Image

# RawImage reflects the RADMC convention that both x and y are increasing
# with array index. This means that to display the image as RADMC intends it,
# you must set the first array element to the lower left corner.
type RawImage <: Image
    data::Array{Float64, 3} # [ergs/s/cm^2]
    pixsize_x::Float64
    pixsize_y::Float64
    lams::Vector{Float64}
end

# SkyImage is a holder that has both RA and DEC increasing with array index
# This convention is necessary for the FFT step
# However, to display this image in the traditional sky convention (North up,
# East to the left), you must set the first array element to the lower left
# corner *and* flip the array along the RA axis: `fliplr(data)`
type SkyImage <: Image
    data::Array{Float64, 3} # [Jy/pixel]
    ra::Vector{Float64} # [arcsec]
    dec::Vector{Float64} # [arcsec]
    lams::Vector{Float64} # [μm]
end

# SkyImage constructor for just a single frame
SkyImage(data::Matrix{Float64}, ra::Vector{Float64}, dec::Vector{Float64}, lam::Float64) =
SkyImage(reshape(data, tuple(size(data)..., 1)), ra, dec, [lam])

# Read the image file (default=image.out) and return it as an Image object, which contains the fluxes in Jy/pixel,
# the sizes and locations of the pixels in arcseconds, and the wavelengths corresponding to the images
function imread(file="image.out")

    fim = open(file, "r")
    iformat = int(readline(fim))
    im_nx, im_ny = split(readline(fim))
    im_nx = int(im_nx)
    im_ny = int(im_ny)
    nlam = int(readline(fim))
    pixsize_x, pixsize_y = split(readline(fim))
    pixsize_x = float64(pixsize_x)
    pixsize_y = float64(pixsize_y)

    # Read the wavelength array
    lams = Array(Float64, nlam)
    for i=1:nlam
        lams[i] = float64(readline(fim))
    end

    # Create an array with the proper size, and then read the file into it
    data = Array(Float64, (im_ny, im_nx, nlam))

    # According to the RADMC manual, section A.15, the pixels are ordered
    # left to right (increasing x) in the inner loop, and from bottom to top
    # (increasing y) in the outer loop.
    # Basically, pack the array in order but display with origin=lower.

    # Because of the way an image is stored as a matrix, we actually pack the
    # array indices as data[y, x, lam]
    # radmc3dPy achieves something similar by keeping indices the x,y but
    # swaping loop order (radmcPy/image.py:line 675)
    for k=1:nlam
        readline(fim) # Junk space
        for j=1:im_ny
            for i=1:im_nx
                data[j,i,k] = float64(readline(fim))
            end
        end
    end

    close(fim)

    # According to the RADMC3D manual, the units are *intensity* [erg cm−2 s−1 Hz−1 ster−1]

    return RawImage(data, pixsize_x, pixsize_y, lams)
end

# Assumes dpc is parsecs
function imToSky(img::RawImage, dpc::Float64)

    # The RawImage is oriented with North up and East increasing to the left.
    # this means that for the RawImage, the delta RA array goes from + to -

    # However, the SkyImage actually requires RA (ll) in increasing form.
    # Therefore we flip along the RA axis, fliplr(data)

    #println("Min and max intensity ", minimum(img.data), " ", maximum(img.data))
    #println("Pixel size ", img.pixsize_x)
    #println("Steradians subtended by each pixel ",  img.pixsize_x * img.pixsize_y / (dpc * pc)^2)

    #convert from ergs/s/cm^2/Hz/ster to to Jy/ster
    conv = 1e23 # [Jy/ster]

    # Conversion from erg/s/cm^2/Hz/ster to Jy/pixel at 1 pc distance.
    # conv = 1e23 * img.pixsize_x * img.pixsize_y / (dpc * pc)^2

    dataJy = fliplr(img.data) .* conv

    (im_ny, im_nx) = size(dataJy)[1:2] #y and x dimensions of the image

    # The locations of pixel centers in cm
    # if n_x = 16, goes [-7.5, -6.5, ..., -0.5, 0.5, ..., 6.5, 7.5] * pixsize
    xx = ((Float64[i for i=0:im_nx-1] + 0.5) - im_nx/2.) * img.pixsize_x
    yy = ((Float64[i for i=0:im_ny-1] + 0.5) - im_ny/2.) * img.pixsize_y

    # The locations of the pixel centers in relative arcseconds
    # Note both RA and DEC increase with array index.
    ra = xx./(AU * dpc)
    dec = yy./(AU * dpc)

    return SkyImage(dataJy, ra, dec, img.lams)

end

# Given an array of imgs, concatenate them into a single image, in wavelength order
function catImages(imgs::Array{Image, 1})

    # Assert each image has the same number of pixels, and that they are all either RawImages or SkyImages

    # Determine the stacking order based upon each lam

    # Stack each img together along the spectral dimension

    # cat the lam together
    0
end

# Take an image and integrate all the frames to create a spatially-integrated spectrum
function imToSpec(img::SkyImage)

    # pixels in SkyImage are Jy/ster

    dRA = abs(img.ra[2] - img.ra[1])
    dDEC = abs(img.dec[2] - img.dec[1])

    flux = squeeze(sum(img.data .* dRA .* dDEC, (1, 2)), (1,2)) # Add up all the flux in the pixels to create the spectrum
    spec = hcat(img.lams, flux) #First column is wl, second is flux

    return spec
end

end

# Move it into Images.jl format.

# subclass AbstractImageDirect

# Or, might be able to store this directly in an Image

# Add a new type to Images.jl, `RADMC3D` to read `image.out` files

# using Images

#type RADMC3D <: Images.ImageFileType end
#add_image_file_format(".out", b"RADMC3D Image", RADMC3D)

#import Images.imread
#function imread{S<:IO}(stream::S, ::Type{RADMC3D})
#    seek(stream, 0)
#    l = strip(readline(stream))
#    l == "RADMC3D Image" || error("Not a RADMC3D file: " * l)
