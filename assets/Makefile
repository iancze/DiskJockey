# The files output by the plot_structure command
STRUCT := $(shell DJ_plot_structure.jl -M)

# The files output by the DJ_initialize.jl command
# Read the dependencies from the program itself (it knows which molecule files to use.)
OBJS := $(shell DJ_initialize.jl -M)

# Plot the structure
.PHONY : structure
.NOTPARALLEL : structure
structure : $(STRUCT)

$(STRUCT) : config.yaml
	DJ_plot_structure.jl

.NOTPARALLEL : $(OBJS)
$(OBJS) : config.yaml
	DJ_initialize.jl

image.out : $(OBJS)
	DJ_synthesize_model.jl

image_tausurf.out : $(OBJS)
	DJ_synthesize_tausurf.jl

.PHONY : chmaps
chmaps : chmaps_linear.png chmaps_log.png

chmaps_linear.png : image.out
	DJ_plot_chmaps.jl --linear

chmaps_log.png : image.out
	DJ_plot_chmaps.jl --log

spectrum.png : image.out
	DJ_plot_chmaps.jl --spectrum

.PHONY : tausurf
tausurf : tausurf.png

tausurf.png : image_tausurf.out
	DJ_plot_tausurf.jl


# Make the visibilities for both the model and residuals.
.NOTPARALLEL : vis
.PHONY : vis
vis : model.hdf5 resid.hdf5

model.hdf5 resid.hdf5 : image.out
	DJ_write_model.jl

all: structure chmaps spectrum.png vis

.PHONY: clean
clean:
	$(RM) $(OBJS)
	$(RM) $(STRUCT)
	$(RM) image.out
	$(RM) radmc3d.out
	$(RM) chmaps_linear.png
	$(RM) chmaps_log.png
	$(RM) chmaps_blur.png
	$(RM) spectrum.png
	$(RM) model.hdf5
	$(RM) resid.hdf5
