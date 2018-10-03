module DiskJockey

export constants,
    model,
    image,
    gridding,
    gauss,
    visibilities,
    EnsembleSampler,
    RADMC3D_PATH

# These statements just straight up dump the source code directly here, making DiskJockey.jl
# act as one giant file.
include("constants.jl")
include("model.jl")
include("gauss.jl")
include("image.jl")
include("gridding.jl")
include("visibilities.jl")
include("EnsembleSampler.jl")

# Add RADMC-3D to the executable path
# using DiskJockey
# srcdir = dirname(pathof(DiskJockey))
# assets_dir = srcdir * "/../assets/"
function __init__()
    unixpath = "../deps/src/radmc-3d/version_0.41/src"
    # @__FILE__ gives us an absolute path to the location where the user installed DiskJockey
    RADMC3D_PATH = normpath(joinpath(dirname(@__FILE__), unixpath))
    ENV["PATH"] = RADMC3D_PATH * ":" * ENV["PATH"]
end

end
