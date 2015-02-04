using consonance

push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")
push!(LOAD_PATH, "/n/home07/iczekala/JudithExcalibur/")

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

# matplotlib.colors.hsv_to_rgb(hsv)

# Get uu,vv spacings
# dv = DataVis("../data/V4046Sgr/V4046Sgr.hdf5", 12)
# # spacings in klambda
# uu = dv.uu
# vv = dv.vv



function scale(data)
    s = maximum(abs(data))
    return norm = plt.Normalize(vmin=-s, vmax=s, clip=false)
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
        end
    end
    RGB = PyPlot.matplotlib[:colors][:hsv_to_rgb](HSV)
    return RGB
end

function plot_vis(img::image.SkyImage, id::Int)

    # Show the real and imaginary components

    fig, ax = plt.subplots(ncols=nchan, figsize=(8, 1.6))
    # Image needs to be flipped along uu dimension

    for col=1:nchan

        # Set all labels but the leftmost to blank
        if col != 1
            ax[col][:xaxis][:set_ticklabels]([])
            ax[col][:yaxis][:set_ticklabels]([])
            ax[col][:xaxis][:set_ticklabels]([])
            ax[col][:yaxis][:set_ticklabels]([])
        else
            ax[col][:set_xlabel](L"uu [k$\lambda$]", size=8)
            ax[col][:set_ylabel](L"vv [k$\lambda$]", size=8)
            ax[col][:tick_params](axis="both", which="major", labelsize=8)
            labels = ax[col][:get_xticklabels]()
            for label in labels
                label[:set_rotation](40)
            end
        end

        vis_fft = transform(img, col) # Transform this channel

        # Crop 128 on either side
        ncrop = 64
        ext = (vis_fft.uu[end - ncrop], vis_fft.uu[1 + ncrop], vis_fft.vv[1 + ncrop], vis_fft.vv[end - ncrop])
        frame = fliplr(vis_fft.VV)[1 + ncrop:end-ncrop, 1 + ncrop:end-ncrop]
        ax[col][:imshow](complex_to_RGB(frame), interpolation="none", origin="lower", extent=ext)

    end

    # Annotate the parameters
    p = pars[id]
    mass = p.M_star
    r_c = p.r_c
    incl = p.incl
    label = L"$M_\ast$: " * @sprintf("%.2f", mass) * L" $M_\odot$   $r_c$: " * @sprintf("%.0f", r_c) * L" AU   $i$: " * @sprintf("%.0f", incl) * L"${}^\circ$"
    fig[:text](0.35, 0.1, label)

    fig[:subplots_adjust](wspace=0.08, top=0.95, bottom=0.26, left=0.1, right=0.9)

    plt.savefig(outdir * @sprintf("vis%04d.png", id))
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
