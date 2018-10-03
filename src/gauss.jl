"Provide routines for making fake Gaussian models which have analytic Fourier transforms, for the purposes of testing."
module gauss

export imageGauss, FTGauss

using ..constants


"
    imageGauss(ll::AbstractVector{Float64}, mm::AbstractVector{Float64}, p::Vector{Float64}, k::Real; theta=0)

Given two arrays of ``l`` and ``m`` coordinates, corresponding to ``x`` and ``y``, fill an array of the Gaussian image
following the MATLAB convention for images (each row corresponds to a different y value).

`p0` is a vector of `[mu_RA, mu_DEC, sigma_x, sigma_y, rho]` in units of arcseconds.
``\\rho`` is the correlation of the Gaussian, ranging from 0 to 1.

```math
\\rho = \\frac{\\sigma_{xy}}{\\sigma_x \\sigma_y}
```

``\\mu_\\alpha`` and ``\\mu_\\delta`` are the locations of the centroid emission relative to the
image origin (RA=0, DEC=0).

``k`` is a scaling pre-factor to adjust the amplitude of the Gaussian.

The image intensity as a function of ``l`` and ``m`` (sky plane coordinates) is

```math
    I(l,m) = \\frac{k}{2 \\pi \\sqrt{|\\boldsymbol{\\Sigma}|}} \\exp \\left \\{ -\\frac{1}{2} \\boldsymbol{R}^\\mathrm{T} \\boldsymbol{\\Sigma}^{-1} \\boldsymbol{R} \\right \\}
```

Where ``\\boldsymbol{\\Sigma}`` is
```math
\\boldsymbol{\\Sigma} = \\left [ \\begin{array}{cc}
\\sigma_x & \\sigma_{xy} \\\\
\\sigma_{xy} & \\sigma_y \\\\
\\end{array} \\right ]
```

``\\theta`` is a desired angle to rotate the Gaussian about the origin. To make sense of the order of operations, the rotation must occur about the origin. A positive value means a counter-clockwise rotation of the Gaussian from North towards East.
"
function imageGauss(ll::AbstractVector{Float64}, mm::AbstractVector{Float64}, p::Vector{Float64}, k::Real; theta=0)

    # Both ll and mm increase with array index
    nx = length(ll)
    ny = length(mm)

    img = Array{Float64}(ny, nx)
    mu = p[1:2] * arcsec #ll and mm shifts
    sigma_x = p[3] * arcsec
    sigma_y = p[4] * arcsec
    rho = p[5]
    sigma_xy = rho * sigma_x * sigma_y
    Sigma = Float64[[sigma_x^2, sigma_xy] [sigma_xy, sigma_y^2]]
    rot = [[cosd(theta), +sind(theta)] [-sind(theta), cosd(theta)]]
    # Sigma = Diagonal((p[3:4] * arcsec).^2) #Convert from arcsec to radians
    pre = 1. / (2pi * sqrt(det(Sigma))) * k
    for j=1:ny
        for i=1:nx
            # Rotation of Gaussian occurs about origin, so subtract off mean first.
            R = rot * (Float64[ll[i] , mm[j]] - mu)
            img[j, i] = pre * exp(-0.5 * (R' * (Sigma\R))[1]) # backslash solves for Sigma^{-1}R
        end
    end
    return img
end

"
    FTGauss(uu::Float64, vv::Float64, p::Vector{Float64}, k::Real)

Given u and v coordinates in [kλ], evaluate the analytic FT of the
aforementioned Gaussian.

`p` is a length 5 vector of `[mu_RA, mu_DEC, sigma_x, sigma_y, rho]` in units of arcseconds, corresponding to the **image** plane.

N.B. Here Sigma refers to the (same) covariance matrix in the *image* domain.

This function always returns a complex value."
function FTGauss(uu::Float64, vv::Float64, p::Vector{Float64}, k::Real; theta=0)
    uu = uu .* 1e3 #[λ]
    vv = vv .* 1e3 #[λ]

    rot = Float64[[cosd(theta), -sind(theta)] [sind(theta), cosd(theta)]]

    mu = p[1:2] * arcsec #ll and mm shifts
    sigma_x = p[3] * arcsec
    sigma_y = p[4] * arcsec
    rho = p[5]
    sigma_xy = rho * sigma_x * sigma_y
    Sigma = Float64[[sigma_x^2, sigma_xy] [sigma_xy, sigma_y^2]]

    # mu_RA, mu_DEC = p[1:2]
    # mu = Float64[mu_RA, mu_DEC] * arcsec #ll and mm shifts

    R0 = Float64[uu, vv]

    # Position to evaluate the rotated Gaussian centered at the origin. Evaluating the Gaussian as if we were
    # in uu', vv'.
    R = rot * R0
    model = k * exp(-2 * (pi^2) * (R' * Sigma * R)[1])
    # in this case, Sigma serves as an inverse matrix, since the covariance matrix corresponds to the *image* plane.
    # Sigma = Diagonal((p[3:4] * arcsec).^2) #Convert from arcsec to radians

    # Phase shift should be evaluated using the original coordinates and mu, since this is the direction that we are actually moving things. And, the interferometer samples uu,vv, not uu', vv'.
    phase_shift = exp(-2pi * 1.0im * (R0' * mu)[1]) # Not actually in polar phase form

    return model * phase_shift
end

"
    FTGauss(uu::AbstractVector{Float64}, vv::AbstractVector{Float64}, p::Vector{Float64}, k::Int)

`p` is a length 5 vector of `mu_RA`, `mu_DEC`, `sigma_RA`, `sigma_DEC`, `rho`.

Given two arrays of u and v coordinates in [kλ], fill an array with the
analytic FT of `imageGauss` evaluated at every pairwise (u,v) pair."
function FTGauss(uu::AbstractVector{Float64}, vv::AbstractVector{Float64}, p::Vector{Float64}, k::Real; theta=0)
    nu = length(uu)
    nv = length(vv)
    # Both uu and vv increase with array index
    img = Array{Complex128}(nv, nu)
    for j=1:nv
        for i=1:nu
            img[j, i] = FTGauss(uu[i], vv[j], p, k, theta=theta)
        end
    end
    return img
end

end # module
