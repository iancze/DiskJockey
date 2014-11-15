# Plot the visibilities

module plot_vis

export plot_RawModelVis

using visibilities
using constants

import PyPlot.plt
using LaTeXStrings

function plot_RawModelVis(vis::RawModelVis)
    fig, ax = plt.subplots(nrows=2, figsize=(12, 2.8))

    rvis = real(vis.VV)
    cvis = imag(vis.VV)

    ax[1][:imshow](rvis, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"))
    ax[2][:imshow](cvis, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"))

    plt.savefig("plots/visibilities.png")
    
end


end
