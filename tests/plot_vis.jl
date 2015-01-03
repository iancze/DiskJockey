# Plot the visibilities

module plot_vis

export plot_RawModelVis

using visibilities
using constants

import PyPlot
import PyPlot.plt
using LaTeXStrings

FSF(str) = PyPlot.matplotlib[:ticker][:FormatStrFormatter](str)

function plot_RawModelVis(vis::RawModelVis)
    fig, ax = plt.subplots(nrows=2, figsize=(5, 8))

    FullVis = fillModelVis(vis)
    VV = FullVis.VV
    uu = FullVis.uu
    vv = FullVis.vv

    rvis = real(VV)
    cvis = imag(VV)

    ext = (minimum(uu), maximum(uu), minimum(vv), maximum(vv))

    rimg = ax[1][:imshow](rvis, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
    ax[1][:set_title]("real")
    ax[1][:set_xlabel](L"uu [k$\lambda$]")
    ax[1][:set_ylabel](L"vv [k$\lambda$]")
    ax[1][:xaxis][:set_major_formatter](FSF("%.0f"))
    ax[1][:yaxis][:set_major_formatter](FSF("%.0f"))

    #[left, bottom, width, height]
    cax = fig[:add_axes]([0.84, 0.65, 0.03, 0.25])
    cb = fig[:colorbar](rimg, cax=cax)

        #ticks = np.linspace(0, np.max(cov), num=6)
        #cb.set_ticks(ticks)
        #cb.set_ticks(MaxNLocator(nbins=5))

    cimg = ax[2][:imshow](cvis, interpolation="none", origin="upper", cmap=plt.get_cmap("Greys"), extent=ext)
    ax[2][:set_title]("imag")
    ax[2][:set_xlabel](L"uu [k$\lambda$]")
    ax[2][:set_ylabel](L"vv [k$\lambda$]")
    ax[2][:xaxis][:set_major_formatter](FSF("%.0f"))
    ax[2][:yaxis][:set_major_formatter](FSF("%.0f"))

    cax = fig[:add_axes]([0.84, 0.15, 0.03, 0.25])
    cb = fig[:colorbar](cimg, cax=cax)

    fig[:subplots_adjust](left=0.15, right=0.85, hspace=0.2)
    plt.savefig("plots/visibilities.png")

end


end
