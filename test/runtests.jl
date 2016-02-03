using Base.Test

# check to see what's on our path
println(ENV["PATH"])

# include("simple_math.jl")
include("ensemble_sampler.jl")

# Make a fake model so we can test off of it.
run(`JudithInitialize.jl --new-project=standard`)

run(`python edit_config.py --fix-distance=True --model=standard`)
include("initialize.jl")

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
