push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")
push!(LOAD_PATH, "/n/home07/iczekala/JudithExcalibur/")

# Determine the proper vmin and vmax scaling

using image

# Load all of the images into memory.
# Determine which image*.out files are in this directory.
imgfunc = x -> contains(x, "image") && contains(x, ".out")
imglist = filter(imgfunc, readdir(outdir))

nimg = length(imglist)

# load one image to determine the dimensions
im1 = imread(imglist[1])
skim1 = imToSky(im1, 73.)

shape = size(skim1.data)
newshape = tuple(shape..., nimg)

images = Array(Float64, newshape)

# Go through each one, read them into an array.
# Determine the maximum and minimum of this array
for (i,fname) in enumerate(imglist)
    im = imread(outdir * fname)
    skim = imToSky(im, 73.)
    images[:,:,:,i] = skim.data
end

println(log10(maximum(images)))
