module DiskJockey

export constants,
    model,
    image,
    gridding,
    visibilities,
    EnsembleSampler,
    RADMC3D_PATH

# These statements just straight up dump the source code directly here, making DiskJockey.jl
# act as one giant file.
include("constants.jl")
include("model.jl")
include("image.jl")
include("gridding.jl")
include("visibilities.jl")
include("EnsembleSampler.jl")

# Add RADMC-3D to the executable path
unixpath = "../deps/src/radmc-3d/version_0.41/src"
# @__FILE__ gives us an absolute path to the location where the user installed DiskJockey
const RADMC3D_PATH = normpath(joinpath(dirname(@__FILE__), unixpath))
ENV["PATH"] = RADMC3D_PATH * ":" * ENV["PATH"]

end
