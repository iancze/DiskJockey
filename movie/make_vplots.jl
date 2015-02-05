using consonance

push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")
push!(LOAD_PATH, "/n/home07/iczekala/JudithExcalibur/")
push!(LOAD_PATH, "/pool/scout0/JudithExcalibur/")

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    # "--opt1"
    # help = "an option with an argument"
    "--run_index", "-r"
    help = "Output run index"
    arg_type = Int
end

parsed_args = parse_args(ARGS, s)
run_index = parsed_args["run_index"]

using constants
using image
using visibilities

import PyPlot.plt
using LaTeXStrings

# Log scale the visibility data in intensity
# minv, maxv are passed in true value
function scale(data::Complex128, maxv::Float64, decades::Float64=0.4)
    if log10(abs(data)) > (log10(maxv) - decades)
        return (log10(abs(data)) - (log10(maxv) - decades)) / decades
    else
        return 0.0
    end
end

# Given a 2D matrix, determine the appropriate color scaling
function complex_to_RGB(frame::Matrix{Complex128})
    ny, nx = size(frame)
    HSV = Array(Float64, (ny, nx, 3))
    minv = minimum(abs(frame))
    maxv = maximum(abs(frame))
    delta = maxv - minv
    for j=1:ny
        for i=1:nx
            HSV[j,i,:] = Float64[(angle(frame[j, i]) + pi)/(2pi), 1.0, (abs(frame[j, i]) - minv)/delta]
            # HSV[j,i,:] = Float64[(angle(frame[j, i]) + pi)/(2pi), 1.0, scale(frame[j, i], minv, maxv)]
        end
    end
    RGB = PyPlot.matplotlib[:colors][:hsv_to_rgb](HSV)
    return RGB
end

norm = plt.Normalize(vmin=vmax - 6, vmax=vmax, clip=false)

function plot_vis(img::image.SkyImage, id::Int)

    fig, ax = plt.subplots(ncols=nchan, nrows=2, figsize=(8, 3.0))

    ext_im = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    for col=1:nchan

        # Plot channel map on top row
        if col != 1
            ax[1,col][:xaxis][:set_ticklabels]([])
            ax[1,col][:yaxis][:set_ticklabels]([])
            ax[1,col][:annotate](@sprintf("%.1f", vels[col]), (0.1, 0.8), xycoords="axes fraction", size=8)
        else
            ax[1,col][:set_xlabel](L"$\Delta \alpha$ ('')", size=8)
            ax[1,col][:set_ylabel](L"$\Delta \delta$ ('')", size=8)
            ax[1,col][:tick_params](axis="both", which="major", labelsize=8)
            ax[1,col][:annotate](@sprintf("%.1f km/s", vels[col]), (0.1, 0.8), xycoords="axes fraction", size=8)
        end

        #Flip the frame for Sky convention
        frame = fliplr(img.data[:,:,col])
        frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
        lframe = log10(frame)
        max = maximum(lframe)
        ax[1,col][:imshow](lframe, extent=ext_im, norm=norm, interpolation="none", origin="lower", cmap=plt.get_cmap("PuBu"))



        # Plot visibilities on bottom

        # Set all labels but the leftmost to blank
        ax[2,col][:tick_params](axis="both", colors="white")
        if col != 1
            ax[2,col][:xaxis][:set_ticklabels]([])
            ax[2,col][:yaxis][:set_ticklabels]([])
            ax[2,col][:xaxis][:set_ticklabels]([])
            ax[2,col][:yaxis][:set_ticklabels]([])
        else
            ax[2,col][:set_xlabel](L"uu [k$\lambda$]", size=8)
            ax[2,col][:set_ylabel](L"vv [k$\lambda$]", size=8)
            ax[2,col][:tick_params](axis="both", which="major", labelsize=8)
            labels = ax[2,col][:get_xticklabels]()
            for label in labels
                label[:set_rotation](40)
                label[:set_color]("black")
            end
            for label in ax[2,col][:get_yticklabels]()
                label[:set_color]("black")
            end
        end

        vis_fft = transform(img, col) # Transform this channel

        # Crop 128 on either side to show more detail
        ncrop = 64
        ext = (vis_fft.uu[end - ncrop], vis_fft.uu[1 + ncrop], vis_fft.vv[1 + ncrop], vis_fft.vv[end - ncrop])
        # Vis needs to be flipped along uu dimension
        frame = fliplr(vis_fft.VV)[1 + ncrop:end-ncrop, 1 + ncrop:end-ncrop]
        ax[2,col][:imshow](complex_to_RGB(frame), interpolation="none", origin="lower", extent=ext)

    end

    # Annotate the parameters
    p = pars[id]
    mass = p.M_star
    r_c = p.r_c
    incl = p.incl
    label = L"$M_\ast$: " * @sprintf("%.2f", mass) * L" $M_\odot$   $r_c$: " * @sprintf("%.0f", r_c) * L" AU   $i$: " * @sprintf("%.0f", incl) * L"${}^\circ$"
    fig[:text](0.35, 0.1, label)

    fig[:subplots_adjust](wspace=0.08, hspace=0.53, top=0.95, bottom=0.26, left=0.1, right=0.9)

    plt.savefig(outdir * @sprintf("vis%04d.png", id), dpi=150)
    plt.close("all")

end

# How many frames per process?
start = run_index * nframes_per_proc + 1

ids = Int[i for i=start:(start + nframes_per_proc)]

for id in ids
    if id <= nframes
        fname = outdir * @sprintf("gimage%04d.out", id)
        im = imread(fname)
        skim = imToSky(im, 73.)
        plot_vis(skim, id)
    else
        break
    end
end
