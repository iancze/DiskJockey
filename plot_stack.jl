using visibilities
using image
using plot_vis

img = imread()
skim = imToSky(img, 73.);

outVis = transform(skim);
plot_RawModelVis(outVis)
