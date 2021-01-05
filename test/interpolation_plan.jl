using DiskJockey.constants
using DiskJockey.gridding
using DiskJockey.visibilities
using DiskJockey.image
using DiskJockey.gauss

# load data and get one channel
dv = DataVis("data.hdf5", 1, true) # this shouldn't be flagged=true, but because we changed the masks, etc

# now actually create a FullModelVis to be interpolated using this closure
# # just use image.out
img = imread()

# Convert raw images to the appropriate distance
skim = imToSky(img, 145.0)
# Apply the gridding correction function before doing the FFT
# No shift needed, since we will shift the resampled visibilities
corrfun!(skim)

vis_fft = transform(skim, 1)


# set up the interpolation closure using the dataset and the uu, vv coordinates of the FullModelVis
interpolate_func = plan_interpolate(dv, vis_fft.uu, vis_fft.vv)


# Interpolate the `vis_fft` to the same locations as the DataSet using the closure
mvis = interpolate_func(dv, vis_fft)
