using Base.Test

# check to see what's on our path and what directory we are in
println(ENV["PATH"])
println("Current working directory ", pwd())
println("Files in this directory ")
run(`ls`)

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
run(`./clean.sh`)

# Make a fake model so we can test off of it.
println("Initializing directory with a new standard project.")
run(`JudithInitialize.jl --new-project=standard`)

# Set up the config file and pos0.npy
println("Editing config.yaml to use the standard model with fixed distance.")
run(`python edit_config.py --fix-distance --model=standard`)

# Run the initialization and plotting scripts
include("initialize.jl")

# The tricky thing is that this sampling will go on for a long time. How can we get the essential test of venus.jl without needing to wait for a full iteration across 24 walkers?  I mean, is it possible to run with fewer walkers than parameters?

# Run venus.jl and see if we sample properly.
run(`venus.jl --test`)

######################################
### Standard model, floating distance

# Just make sure the directory is clean
run(`./clean.sh`)

# Make a fake model so we can test off of it.
println("Initializing directory with a new standard project.")
run(`JudithInitialize.jl --new-project=standard`)

# Set up the config file and pos0.npy
println("Editing config.yaml to use the standard model with floating distance.")
run(`python edit_config.py --model=standard`)

# Run the initialization and plotting scripts
include("initialize.jl")

# Run venus.jl and see if we sample properly.
run(`venus.jl --test`)


# Now, see if we can run venus.jl with d fixed

# Edit config.yaml (how?) to un-fix d
# Use pyyaml to do the edits.

# Edit config.yaml to exclude most channels.


# Now, with d floating
# run(`python edit_config.py --fix-distance=False --model=standard`)
# include("initialize.jl")
#
#
# # Now, run with a cavity model
# # Make a fake model so we can test off of it.
# run(`JudithInitialize.jl --new-project=cavity`)
#
# run(`python edit_config.py --fix-distance=True --model=cavity`)
# include("initialize.jl")
#
#
# run(`JudithInitialize.jl --new-project=cavity`)
# run(`python edit_config.py --fix-distance=False --model=cavity`)
# include("initialize.jl")



# include("radmc3d.jl")

# Try running a quick venus.jl

# Clean the directory
# Repeat everything with truncated model
# run(`JudithInitialize.jl --new-project=truncated`)

# Clean the directory
# Repeat everything with cavity model
# run(`JudithInitialize.jl --new-project=cavity`)
