push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")
push!(LOAD_PATH, "/n/home07/iczekala/JudithExcalibur/")
push!(LOAD_PATH, "/pool/scout0/JudithExcalibur/")

# Determine the proper vmin and vmax scaling for the dust images, channel maps,
# and transforms of channel maps

using image

# Store a global max value, and then go through each image, see if it beats
# this value, if so, record the new one, then close image to keep the memory
# low.

fname = outdir * @sprintf("dimage%04d.out", id)
im = imread(fname)
skim = imToSky(im, 73.)


global vmax = maximum(skim.data)
global vmin = minimum(skim.data)

# Load all of the images into memory.
# Determine which image*.out files are in this directory.
imgfunc = x -> contains(x, strkey) && contains(x, ".out")
imglist = filter(imgfunc, readdir(outdir))

println(imglist)

nimg = length(imglist)

# Go through each one, determine the maximum and minimum, and if it beats the
# current value, replace it
for fname in imglist
    im = imread(outdir * fname)
    skim = imToSky(im, 73.)
    min,max = extrema(skim.data)
    if max > vmax
        vmax = max
    end
    if min < vmin
        vmin = min
    end
end

println("Dust minimum and maximum ", vmin, vmax)
