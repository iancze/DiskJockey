# Provides FFT functionality, regridding, 

module visibilities

using image

export transform, RawModelVis

# Stores the visibilities for a single channel
type DataVis
    lam::Float64 # Wavelength (in microns) corresponding to channel
    uu::Vector{Float64} # Vectors of the u, v locations in UNITS
    vv::Vector{Float64}
    VV::Array{Complex{Float64}, 2} # Complex visibilities
end

# Produced by FFT'ing a single channel of a SkyImage
type RawModelVis
    lam::Float64 # Wavelength (in microns) corresponding to channel
    uu::Vector{Float64} # Vectors of the u, v locations in UNITS
    vv::Vector{Float64}
    VV::Array{Complex{Float64}, 2} # Output from rfft
end

# Produced by gridding a RawModelVis to match the data
# type GridModelVis
# end


# Transform the SkyImage produced by RADMC and 
function transform(img::SkyImage, index::Int=1)

    # By default, select the first channel of any spectral hypercube.
    data = img.data[:, :, index]

    lam = img.lams[index]

    # converted from ra and dec
    uu = Array(Float64, (1,))
    vv = Array(Float64, (1,))

    # properly pack the data for input 
    #TODO: Check that the outer fftshift is right
    out = rfft(fftshift(data))
    return RawModelVis(lam, uu, vv, out)
end


end
