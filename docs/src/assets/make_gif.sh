#!/bin/bash

# Run this after you have already made the mp4


#ffmpeg -i gas_out.mp4 -vf scale=800:-1 -r 10 -f image2pipe -vcodec ppm - | convert -delay 5 -loop 0 - output.gif
ffmpeg -i gas.mp4 -vf scale=1200:-1 -r 10 -f image2pipe -vcodec ppm - | convert -delay 5 -loop 0 - gas.gif
ffmpeg -i dust.mp4 -vf scale=400:-1 -r 10 -f image2pipe -vcodec ppm - | convert -delay 5 -loop 0 - dust.gif
ffmpeg -i vis.mp4 -vf scale=1200:-1 -r 10 -f image2pipe -vcodec ppm - | convert -delay 5 -loop 0 - vis.gif
