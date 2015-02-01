# Using ffmpeg, create a movie
using consonance

ints = outdir * "%04d.png"
outs = outdir * "out.mp4"
run(`ffmpeg -framerate 3.0 -i $ints -c:v libx264 -r 30 -pix_fmt yuv420p $outs`)
# run(`ffmpeg -framerate 3.0 -i %04d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4`)
