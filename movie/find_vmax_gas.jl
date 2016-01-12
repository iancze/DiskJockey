push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/movie/")
# push!(LOAD_PATH, "/n/home07/iczekala/JudithExcalibur/")
# push!(LOAD_PATH, "/pool/scout0/JudithExcalibur/")

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    # "--opt1"
    # help = "an option with an argument"
    "--run_index", "-r"
    help = "Output run index"
    default = 0
    arg_type = Int
end

parsed_args = parse_args(ARGS, s)
run_index = parsed_args["run_index"]


# Determine the proper vmin and vmax scaling for the dust images, channel maps,
# and transforms of channel maps
using consonance

using JudithExcalibur.image
using JudithExcalibur.visibilities

import YAML
config = YAML.load(open("config.yaml"))

# Store a global max value, and then go through each image, see if it beats
# this value, if so, record the new one, then close image to keep the memory
# low.

imgdir = moviedir * @sprintf("g%02d/", run_index)

fname = imgdir * @sprintf("gimage%04d.out", 1)
im = imread(fname)

dpc = config["parameters"]["dpc"][1]
skim = imToSky(im, dpc)
nchan = length(skim.lams)
println("$nchan channels")

global vmax_g = maximum(skim.data)
global vmin_g = minimum(skim.data)

vis_fft = transform(skim, 1)

global vmax_v = maximum(abs(vis_fft.VV))
global vmin_v = minimum(abs(vis_fft.VV))

# Load all of the images into memory.
# Determine which image*.out files are in this directory.
imgfunc = x -> contains(x, "gimage") && contains(x, ".out")
imglist = filter(imgfunc, readdir(imgdir))

nimg = length(imglist)

# Go through each one, determine the maximum and minimum, and if it beats the
# current value, replace it
for fname in imglist
    im = imread(imgdir * fname)
    skim = imToSky(im, dpc)
    min,max = extrema(skim.data)
    if max > vmax_g
        vmax_g = max
    end
    if min < vmin_g
        vmin_g = min
    end
    for i=1:nchan
        vis_fft = transform(skim, i)
        min,max = extrema(abs(vis_fft.VV))
        if max > vmax_v
            vmax_v = max
        end
        if min < vmin_v
            vmin_v = min
        end
    end


end

println("Gas minimum and maximum $vmin_g, $vmax_g")
println("Visibility minimum and maximum $vmin_v, $vmax_v")
