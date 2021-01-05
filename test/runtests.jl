using Test

# Check to see if an environment variable is defined. If so, save plots there. 
outputdir = get(ENV, "JL_PLOT_TEMPDIR", nothing)
if isnothing(outputdir)
    # If there is no environment variable grab a temporary directory and save all plots there by default 
    outputdir = mktempdir()
else
    # if it doesn't exist, make the directory 
    if !ispath(outputdir)
        mkdir(outputdir)
    end
end

println("Saving plots in test/$outputdir")

include("interpolation.jl") # visibility interpolation

include("ensemble_sampler.jl") # sampler test


