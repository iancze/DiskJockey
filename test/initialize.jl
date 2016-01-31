# Run a script that supposedly creates a new fake disk here.

# Now, synthesize a new set of channel maps
# First, create the necessary input files
run(`JudithInitialize.jl`)

# make plots of the structure
run(`plot_structure.jl`)

# run the synthesis
run(`synthesize_model.jl`)

# plot the channel maps
run(`plot_chmaps.jl`)

# plot the moments
run(`plot_moments.jl`)

# Try downsampling the model to the baseline locations
run(`write_model.jl`)
