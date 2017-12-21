using Base.Test

# check to see what's on our path and what directory we are in
println(ENV["PATH"])
println("Current working directory ", pwd())
println("Files in this directory ")
run(`ls`)

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

# The tricky thing is that this sampling will go on for a long time. How can we get the essential test of venus.jl without needing to wait for a full iteration across 24 walkers?  I mean, is it possible to run with fewer walkers than parameters?

# Run venus.jl and see if we sample properly.
run(`venus.jl --test`)

# Clean up before moving on to the next test
run(`make clean`)

######################################
### Vertical model, fixed distance

# Just make sure the directory is clean
# run(`bash clean.sh`)

# Make a fake model so we can test off of it.
# println("Initializing directory with a new vertical project.")
# run(`DJ_initialize.jl --new-project=vertical`)

# Set up the config file and pos0.npy
# println("Editing config.yaml to use the vertical model with fixed distance.")
# run(`python edit_config.py --fix-distance --model=vertical`)

# Run the initialization and plotting scripts
# run(`make all`)

# Run venus.jl and see if we sample properly. --test flag limits the sampling to a few walkers.
# run(`venus.jl --test`)

# Clean up before moving on to the next test
# run(`make clean`)

######################################
### Cavity model, fixed distance

# Just make sure the directory is clean
# run(`bash clean.sh`)

# Make a fake model so we can test off of it.
# println("Initializing directory with a new cavity project.")
# run(`DJ_initialize.jl --new-project=cavity`)

# Set up the config file and pos0.npy
# println("Editing config.yaml to use the cavity model with fixed distance.")
# run(`python edit_config.py --fix-distance --model=cavity`)

# Run the initialization and plotting scripts
# run(`make all`)

# The tricky thing is that this sampling will go on for a long time. How can we get the essential test of venus.jl without needing to wait for a full iteration across 24 walkers?  I mean, is it possible to run with fewer walkers than parameters?

# Run venus.jl and see if we sample properly.
# run(`venus.jl --test`)

# Clean up before moving on to the next test
# run(`make clean`)
