using Test

# check to see what's on our path and what directory we are in
# println(ENV["PATH"])
# println("Current working directory ", pwd())
# println("Files in this directory ")
# run(`ls`)

##################################
### Visibility interpolation
##################################

include("interpolation.jl")


###################################
### Basic unit tests
###################################

# Do more simple unit tests here
include("ensemble_sampler.jl")

###################################
### Functional, integration tests
###################################

######################################
### Standard model, fixed distance

# Just make sure the directory is clean
run(`bash clean.sh`)

# Make a fake model so we can test off of it.
println("Initializing directory with a new standard project.")
run(`DJ_initialize.jl --new-project=standard`)

# Here test the fix parameters reading routines
include("parameters_fixed.jl")

# Set up the config file and pos0.npy
println("Editing config.yaml to use the standard model with fixed distance.")
run(`python edit_config.py --fix-distance --model=standard`)

# Run the initialization and plotting scripts
run(`make all`)

# Run venus.jl and see if we sample properly. The idea is to get a simple test with a few walkers without running forever.
run(`venus.jl --test`)

# Clean up before moving on to the next test
run(`make clean`)
