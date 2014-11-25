# Plot the visibilities

module plot_vis

export plot_RawModelVis

using visibilities
using constants

import PyPlot.plt
using LaTeXStrings

function plot_RawModelVis(vis::RawModelVis)
    fig, ax = plt.subplots(nrows=2, figsize=(5, 5))

    #Mirror the top and bottom
    vv = fftshift(vis.VV, 2)

    #Flip the image over both axes and concatenate vertically
    #vv2 = copy(vv)
    #vv2 = flipdim(vv2, 2)
    #vv2 = flipdim(vv2, 1)
    
    #vvd = vcat(fftshift(vv2,2), fftshift(vv,2))

    rvis = real(vv)
    cvis = imag(vv)

    ax[1][:imshow](rvis, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"))
    ax[1][:set_title]("real")
    ax[2][:imshow](cvis, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"))
    ax[2][:set_title]("imag")

    plt.savefig("plots/visibilities.png")
    
end


end
