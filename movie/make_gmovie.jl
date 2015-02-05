# Using ffmpeg, create a movie

# following https://trac.ffmpeg.org/wiki/Create%20a%20video%20slideshow%20from%20images
using consonance

ins = "scratch/g%04d.png"
outs = "gas_out.mp4"
run(`ffmpeg -framerate 20.0 -i $ins -c:v libx264 -r 30 -pix_fmt yuv420p $outs`)
# run(`ffmpeg -framerate 3.0 -i %04d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4`)
