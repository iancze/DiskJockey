using Base.Test

# check to see what's on our path
println(ENV["PATH"])

include("simple_math.jl")
include("ensemble_sampler.jl")

# Make a fake model so we can test off of it.
run(`JudithInitialize.jl --new-project=standard`)
include("initialize.jl")

# Now, see if we can run venus.jl with d fixed


# Now, with d floating



# include("radmc3d.jl")

# Try running a quick venus.jl

# Clean the directory
# Repeat everything with truncated model
# run(`JudithInitialize.jl --new-project=truncated`)

# Clean the directory
# Repeat everything with cavity model
# run(`JudithInitialize.jl --new-project=cavity`)
