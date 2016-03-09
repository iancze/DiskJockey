# Run a script that supposedly creates a new fake disk here.

# Now, synthesize a new set of channel maps
# First, create the necessary input files
run(`DJ_initialize.jl`)

# make plots of the structure
run(`DJ_plot_structure.jl`)

# run the synthesis
run(`DJ_synthesize_model.jl`)

# plot the channel maps
run(`DJ_plot_chmaps.jl`)

# plot the moments
# run(`plot_moments.jl`)

# Try downsampling the model to the baseline locations
run(`DJ_write_model.jl`)
