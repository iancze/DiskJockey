push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")
push!(LOAD_PATH, "/n/home07/iczekala/JudithExcalibur/")
push!(LOAD_PATH, "/pool/scout0/JudithExcalibur/")

# Determine the proper vmin and vmax scaling for the dust images, channel maps,
# and transforms of channel maps

using image
using consonance
using visibilities

# Store a global max value, and then go through each image, see if it beats
# this value, if so, record the new one, then close image to keep the memory
# low.

fname = outdir * @sprintf("gimage%04d.out", 1)
im = imread(fname)
skim = imToSky(im, 73.)
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
imglist = filter(imgfunc, readdir(outdir))

nimg = length(imglist)

# Go through each one, determine the maximum and minimum, and if it beats the
# current value, replace it
for fname in imglist
    im = imread(outdir * fname)
    skim = imToSky(im, 73.)
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
