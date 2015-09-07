module JudithExcalibur

export constants,
    model, # normal, axisymmetric disk. Czekala et al. 15, Rosenfeld et al. 12
    emodel, # eccentric disk
    hmodel, # vertical temperature gradient
    image,
    gridding,
    visibilities,
    parallel,
    LittleMC,
    EnsembleSampler

# These statements just straight up dump the source code directly here, making JudithExcalibur.jl
# act as one giant file with multiple module definitions inside the `module JudithExcalibur`
# definition

include("constants.jl")
include("model.jl")
include("emodel.jl")
include("hmodel.jl")
include("image.jl")
include("gridding.jl")
include("visibilities.jl")
include("parallel.jl")
include("LittleMC.jl")
include("EnsembleSampler.jl")

end
