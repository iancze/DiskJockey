module gauss_model

# Provide the functionality to model Gaussians as test cases
using DiskJockey.constants

export imageGauss, FTGauss

# Because of the flipped nature of the sky (but not flipped nature of the UV plane)
# there are some tricky conventions about how to pack the array.

# Given two arrays of l and m coordinates, fill an array of the Gaussian image
# following the MATLAB convention.
# p0 is a vector of [mu_RA, mu_DEC, sigma_x, sigma_y] in units of arcseconds
# mu_RA and mu_DEC are the locations of the centroid emission relative to the
# image origin (RA=0, DEC=0).
function imageGauss(ll::Vector{Float64}, mm::Vector{Float64}, p::Vector{Float64}, k::Int)

    # Both ll and mm increase with array index
    nx = length(ll)
    ny = length(mm)

    img = Array(Float64, ny, nx)
    mu = p[1:2] * arcsec #ll and mm shifts
    Sigma = Diagonal((p[3:4] * arcsec).^2) #Convert from arcsec to radians
    pre = 1. / (2pi * sqrt(det(Sigma))) * k
    for j=1:ny
        for i=1:nx
            R = Float64[ll[i] , mm[j]] - mu
            img[j, i] = pre * exp(-0.5 * (R' * (Sigma\R))[1])
        end
    end
    return img
end

# Given u and v coordinates in [kλ], evaluate the analytic FT of the
# aforementioned Gaussian
# N.B. Here Sigma refers to the (same) covariance matrix in the *image* domain
# always return a complex value
function FTGauss(uu::Float64, vv::Float64, p::Vector{Float64}, k::Int)
    uu = uu .* 1e3 #[λ]
    vv = vv .* 1e3 #[λ]
    mu_RA, mu_DEC = p[1:2]
    mu = Float64[mu_RA, mu_DEC] * arcsec #ll and mm shifts
    R = Float64[uu, vv]
    Sigma = Diagonal((p[3:4] * arcsec).^2) #Convert from arcsec to radians
    phase_shift = exp(-2pi * 1.0im * (R' * mu)[1]) # Not actually in polar phase form
    return k * exp(-2 * (pi^2) * (R' * Sigma * R)[1]) * phase_shift
    # in this case, Sigma serves as the inverse
end

# Given two arrays of u and v coordinates in [kλ], fill an array with the
# analytic FT of aforementioned Gaussian evaluated at every pairwise (u,v) pair
function FTGauss(uu::Vector{Float64}, vv::Vector{Float64}, p::Vector{Float64}, k::Int)
    nu = length(uu)
    nv = length(vv)
    # Both uu and vv increase with array index
    img = Array(Complex128, nv, nu)
    for j=1:nv
        for i=1:nu
            img[j, i] = FTGauss(uu[i], vv[j], p, k)
        end
    end
    return img
end


end # Module
