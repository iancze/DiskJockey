module JudithExcalibur

export constants,
    model,
    image,
    gridding,
    visibilities,
    EnsembleSampler
    # ImportanceSampler

# These statements just straight up dump the source code directly here, making JudithExcalibur.jl
# act as one giant file.

include("constants.jl")
include("model.jl")
include("image.jl")
include("gridding.jl")
include("visibilities.jl")
include("EnsembleSampler.jl")
# include("ImportanceSampler.jl")

end
