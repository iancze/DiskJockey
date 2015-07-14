module JudithExcalibur

export constants,
    model,
    emodel,
    image,
    gridding,
    visibilities,
    parallel,
    LittleMC

# These statements just straight up dump the source code directly here, making JudithExcalibur.jl
# act as one giant file with multiple module definitions inside the `module JudithExcalibur`
# definition

include("constants.jl")
include("model.jl")
include("emodel.jl")
include("image.jl")
include("gridding.jl")
include("visibilities.jl")
include("parallel.jl")
include("LittleMC.jl")

end
