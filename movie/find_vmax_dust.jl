push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/DiskJockey/")
push!(LOAD_PATH, "/n/home07/iczekala/DiskJockey/")
push!(LOAD_PATH, "/pool/scout0/DiskJockey/")

# Determine the proper vmin and vmax scaling for the dust images, channel maps,
# and transforms of channel maps

using image
using consonance

# Store a global max value, and then go through each image, see if it beats
# this value, if so, record the new one, then close image to keep the memory
# low.

fname = outdir * @sprintf("dimage%04d.out", 1)
im = imread(fname)
skim = imToSky(im, 73.)


global vmax_d = maximum(skim.data)
global vmin_d = minimum(skim.data)

# Load all of the images into memory.
# Determine which image*.out files are in this directory.
imgfunc = x -> contains(x, "dimage") && contains(x, ".out")
imglist = filter(imgfunc, readdir(outdir))

nimg = length(imglist)

# Go through each one, determine the maximum and minimum, and if it beats the
# current value, replace it
for fname in imglist
    im = imread(outdir * fname)
    skim = imToSky(im, 73.)
    min,max = extrema(skim.data)
    if max > vmax_d
        vmax_d = max
    end
    if min < vmin_d
        vmin_d = min
    end
end

println("Dust minimum and maximum $vmin_d, $vmax_d")
