# We want to create a spheroid object that can deliver gcffun, corrfun, etc.

module gridding

using ..image
using Printf

import Base.Math.@horner

export spheroid, corrfun, corrfun!, gcffun

"
    spheroid(eta)

`spheroid` function which assumes ``\\alpha`` = 1.0, ``m=6``,
built for speed."
function spheroid(eta::Float64)

    # Since the function is symmetric, overwrite eta
    eta = abs(eta)

    if eta <= 0.75
        nn = eta^2 - 0.75^2

        return @horner(nn, 8.203343E-2, -3.644705E-1, 6.278660E-1, -5.335581E-1, 2.312756E-1)/
            @horner(nn, 1., 8.212018E-1, 2.078043E-1)

    elseif eta <= 1.0
        nn = eta^2 - 1.0

        return @horner(nn, 4.028559E-3, -3.697768E-2, 1.021332E-1, -1.201436E-1, 6.412774E-2)/
            @horner(nn, 1., 9.599102E-1, 2.918724E-1)

    elseif eta <= 1.0 + 1e-7
        # case to allow some floating point error
        return 0.0

    else
        # Now you're really outside of the bounds
        println("The spheroid is only defined on the domain -1.0 <= eta <= 1.0. (modulo machine precision.)")
        throw(DomainError())
    end
end

# Make this function available to call with a vector of etas as well
spheroid(etas::Vector{Float64}) = Float64[spheroid(eta) for eta in etas]


"
    spheroid(eta, alpha)

Prolate spheroidal wavefunction, assuming that ``m = 6``. This allows
arguments for ``\\eta < 1.0 + 10^{-7}``, but returns 0.0 (i.e., the spheroid window is truncated)."
function spheroid(eta::Float64, alpha::Float64)

    etalim::Float64 = 0.75 # Specifically for m = 6

    # Alpha can take on values of 0.0, 0.5, 1.0, 1.5, and 2.0
    # Take alpha and convert to a string for keying purposes
    aa = @sprintf("%.1f", alpha)

    # Since the function is symmetric, overwrite eta
    eta = abs(eta)

    if eta <= etalim

        nn = eta^2 - etalim^2

        if aa == "0.0"
            return @horner(nn, 5.613913E-2,-3.019847E-1, 6.256387E-1, -6.324887E-1, 3.303194E-1)/
                   @horner(nn, 1., 9.077644E-1, 2.535284E-1)

        elseif aa == "0.5"
            return @horner(nn, 6.843713E-2,-3.342119E-1, 6.302307E-1, -5.829747E-1, 2.765700E-1)/
                   @horner(nn, 1., 8.626056E-1, 2.291400E-1)

        elseif aa == "1.0"
            return @horner(nn, 8.203343E-2, -3.644705E-1, 6.278660E-1, -5.335581E-1, 2.312756E-1)/
                   @horner(nn, 1., 8.212018E-1, 2.078043E-1)

        elseif aa == "1.5"
            return @horner(nn, 9.675562E-2,-3.922489E-1, 6.197133E-1, -4.857470E-1, 1.934013E-1)/
                   @horner(nn, 1., 7.831755E-1, 1.890848E-1)

        elseif aa == "2.0"
            return @horner(nn, 1.124069E-1,-4.172349E-1, 6.069622E-1, -4.405326E-1, 1.618978E-1)/
                   @horner(nn, 1., 7.481828E-1, 1.726085E-1)
        else
            println("The spheroid is only defined for alpha = 0.0, 0.5, 1.0, 1.5, and 2.0")
            throw(DomainError())
        end

    elseif eta <= 1.0
        nn = eta^2 - 1.0

        if aa == "0.0"
            return @horner(nn, 8.531865E-4,-1.616105E-2, 6.888533E-2, -1.109391E-1, 7.747182E-2)/
                   @horner(nn, 1., 1.101270   , 3.858544E-1)

        elseif aa == "0.5"
            return @horner(nn, 2.060760E-3,-2.558954E-2, 8.595213E-2, -1.170228E-1, 7.094106E-2)/
                   @horner(nn, 1., 1.025431   , 3.337648E-1)

        elseif aa == "1.0"
            return @horner(nn, 4.028559E-3, -3.697768E-2, 1.021332E-1, -1.201436E-1, 6.412774E-2)/
                   @horner(nn, 1., 9.599102E-1, 2.918724E-1)

        elseif aa == "1.5"
            return @horner(nn, 6.887946E-3,-4.994202E-2, 1.168451E-1, -1.207733E-1, 5.744210E-2)/
                   @horner(nn, 1., 9.025276E-1, 2.575337E-1)

        elseif aa == "2.0"
            return @horner(nn, 1.071895E-2,-6.404749E-2, 1.297386E-1, -1.194208E-1, 5.112822E-2)/
                   @horner(nn, 1., 8.517470E-1, 2.289667E-1)

        else
            println("The spheroid is only defined for alpha = 0.0, 0.5, 1.0, 1.5, and 2.0")
            throw(DomainError())
        end

    elseif eta <= 1.0 + 1e-7
        # case to allow some floating point error
        return 0.0

    else
        # return 0.0
        # Now you're really outside of the bounds
        println("The spheroid is only defined on the domain -1.0 <= eta <= 1.0. (modulo machine precision.)")
        throw(DomainError())
    end
end


"
    corrfun(eta::T) where {T}

Gridding *correction* function, but able to be passed either floating point numbers or vectors of `Float64`."
function corrfun(eta::T) where {T}
    return spheroid(eta)
end

"
    corrfun(eta::T, alpha::Float64) where {T}

Gridding *correction* function, used to pre-divide the image to correct for the effect
of the `gcffun`. This function is also the Fourier transform of `gcffun`."
function corrfun(eta::T, alpha::Float64) where {T}
    return spheroid(eta, alpha)
end

"
    corrfun!(img::SkyImage)

Apply the correction function to (and mutate) a `SkyImage` in place."
function corrfun!(img::SkyImage)
    ny, nx, nlam = size(img.data)

    # The size of one half-of the image.
    # sometimes ra and dec will be symmetric about 0, othertimes they won't
    # so this is a more robust way to determine image half-size
    maxra = abs(img.ra[2] - img.ra[1]) * nx/2
    maxdec = abs(img.dec[2] - img.dec[1]) * ny/2

    for k=1:nlam
        for i=1:nx
            for j=1:ny
                etax = (img.ra[i])/maxra
                etay = (img.dec[j])/maxdec
                if abs(etax) > 1.0 || abs(etay) > 1.0
                    # We would be querying outside the shifted image
                    # bounds, so set this emission to 0.0
                    img.data[j, i, k] = 0.0
                else
                    img.data[j, i, k] = img.data[j, i, k] / (corrfun(etax) * corrfun(etay))
                end
            end
        end
    end
end

"
    corrfun(img::SkyImage)

Apply the correction function to a `SkyImage`, but return a copy of the image, leaving the original unchanged."
function corrfun(img::SkyImage)
    im = deepcopy(img)
    corrfun!(im)
    return im
end

"
    corrfun!(img::SkyImage, mu_RA, mu_DEC)

Apply the correction function to the `SkyImage`, with an offset in RA and DEC."
function corrfun!(img::SkyImage, mu_RA, mu_DEC)
    ny, nx, nlam = size(img.data)

    # The size of one half-of the image.
    # sometimes ra and dec will be symmetric about 0, othertimes they won't
    # so this is a more robust way to determine image half-size
    maxra = abs(img.ra[2] - img.ra[1]) * nx/2
    maxdec = abs(img.dec[2] - img.dec[1]) * ny/2

    # If the image will be later offset via a phase shift, then this means that
    # the corrfunction will need to be applied *as if the image were already
    # offset.*
    # Update 1/23/16: Although this is still mathematically correct, it's better
    # to operation using `corrfun!` without any original shift, and simply shift the visibilities

    for k=1:nlam
        for i=1:nx
            for j=1:ny
                etax = (img.ra[i] + mu_RA)/maxra
                etay = (img.dec[j] + mu_DEC)/maxdec
                if abs(etax) > 1.0 || abs(etay) > 1.0
                    # We would be querying outside the shifted image
                    # bounds, so set this emission to 0.0
                    img.data[j, i, k] = 0.0
                else
                    img.data[j, i, k] = img.data[j, i, k] / (corrfun(etax) * corrfun(etay))
                end
            end
        end
    end
end

"
    gcffun(eta::T) where {T}

The gridding *convolution* function, used to do the convolution and interpolation of the visibilities in
the Fourier domain. This is also the Fourier transform of `corrfun`."
function gcffun(eta::T) where {T}
    return abs.(1 .- eta.^2) .* spheroid(eta)
end

"
    gcffun(eta::T, alpha::Float64) where {T}

The gridding *convolution* function, with variable ``\\alpha``.
"
function gcffun(eta::T, alpha::Float64) where {T}
    return abs(1 - eta.^2).^alpha .* spheroid(eta, alpha)
end


end #Module
