# We want to create a spheroid object that can deliver gcffun, corrfun, etc.

module gridding

import Base.Math.@horner

export spheroid

# Assumes we are using m = 6. 
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

    else
        println("The spheroid is only defined on the domain -1.0 <= eta <= 1.0.")
        throw(DomainError())

    end
end

# Make this function available to call with a vector of etas as well
spheroid(etas::Vector{Float64}, alpha::Float64) = Float64[spheroid(eta, alpha) for eta in etas] 


end #Module
