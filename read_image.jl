# Read the image format written by RADMC3D.

# Read the ascii text and parse things into a 3 dimensional matrix (x, y, lambda)

# Read `image.out`

# The first four lines are format information
# iformat # = 1 (2 is local observer)
# im_nx   im_ny #number of pixels in x and y directions
# nlam # number of images at different wavelengths
# pixsize_x  pixsize_y # size of the pixels in cm 
# lambda[1] ... lambda[nlam + 1] # wavelengths (um) correspending to images
# pixels, ordered from left to right (increasing x) in the inner loop, and from bottom to top (increasing y) in the outer loop. And wavelength is the outermost loop.

using constants
        
fim = open("image.out", "r")
iformat = int(readline(fim))
im_nx, im_ny = split(readline(fim))
im_nx = int(im_nx)
im_ny = int(im_ny)
nlam = int(readline(fim))
pixsize_x, pixsize_y = split(readline(fim))
pixsize_x = float64(pixsize_x)
pixsize_y = float64(pixsize_y)

line=4

#println(iformat)
#print(im_nx, im_ny)
#print(nlam)
#print(pixsize_x, pixsize_y)

# Read the wavelength array
lam = Array(Float64, nlam)
for i=1:nlam
    lam[i] = float64(readline(fim))
    #line += 1
end


# Create an array with the proper size, and then read the file into it
img = Array(Float64, (im_nx, im_ny, nlam))

# Although this order is contrary to the RADMC manual, which states x should be in the inner loop, 
# it is what is necessary to orient the image properly. Also, radmc3dPy has the same ordering in
# image.py:line 675
for k=1:nlam
    readline(fim) # Junk space
    for i=1:im_nx 
        for j=1:im_ny
            img[i,j,k] = float64(readline(fim))
        end
    end
end

close(fim)

println("Image size ", size(img))

conv = pixsize_x * pixsize_y / pc^2 * 1e23 #convert from ergs/s/cm^2/Hz to to Jy/pixel
imgJy = img * conv

#The Ra/Dec axes in cm
xx = ((Float64[i for i=0:im_nx] + 0.5) - im_nx/2.) * pixsize_x
yy = ((Float64[i for i=0:im_ny] + 0.5) - im_ny/2.) * pixsize_y

dpc = 73. #pc 

ra = xx/1.496e13/dpc
dec = yy/1.496e13/dpc

ext = (ra[1], ra[end], dec[1], dec[end])

lam0 = cc/230.538e9 * 1e6 # [microns]
nvels = 23
vels = linspace(-4.4, 4.4, nvels) # [km/s]

import PyPlot.plt
using LaTeXStrings

fig, ax = plt.subplots(nrows=2, ncols=12, figsize=(12, 2.8))

#Plot a blank img in the last frame, moment map will go here eventually
ax[2, 12][:imshow](zeros((im_nx, im_ny)), cmap=plt.get_cmap("Greys"), vmin=0, vmax=20)
ax[2, 12][:xaxis][:set_ticklabels]([])
ax[2, 12][:yaxis][:set_ticklabels]([])

# Loop through all of the different channels and plot them
#for iframe=1:nlam
for row=1:2
    for col=1:12
        iframe = col + (row - 1) * 12

        if iframe > nlam
            break
        end

        if col != 1 || row != 2
            ax[row, col][:xaxis][:set_ticklabels]([])
            ax[row, col][:yaxis][:set_ticklabels]([])
        else
            ax[row, col][:set_xlabel](L"$\Delta \alpha$ ('')")
            ax[row, col][:set_ylabel](L"$\Delta \delta$ ('')")
        end

        ax[row, col][:annotate](@sprintf("%.1f", vels[iframe]), (0.1, 0.8), xycoords="axes fraction", size=8)

        frame = img[:,:,iframe]
        frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
        max = maximum(log10(frame))
        ax[row, col][:imshow](log10(frame), extent=ext, vmin=max - 7, vmax=max, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"))
    end
end

fig[:subplots_adjust](hspace=0.01, wspace=0.05, top=0.9, bottom=0.1, left=0.05, right=0.95)

plt.savefig("plots/channel_maps.png")

# Move it into Images.jl format.

# subclass AbstractImageDirect

# Or, might be able to store this directly in an Image

# Add a new type to Images.jl, `RADMC3D` to read `image.out` files

# using Images

#type RADMC3D <: Images.ImageFileType end
#add_image_file_format(".out", b"RADMC3D Image", RADMC3D)

#import Images.imread
#function imread{S<:IO}(stream::S, ::Type{RADMC3D})
#    seek(stream, 0)
#    l = strip(readline(stream))
#    l == "RADMC3D Image" || error("Not a RADMC3D file: " * l)
