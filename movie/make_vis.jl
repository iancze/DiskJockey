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

function scale(data)
    s = maximum(abs(data))
    return norm = plt.Normalize(vmin=-s, vmax=s, clip=false)
end


function plot_vis(img::image.SkyImage, id::Int)

    # Show the real and imaginary components

    fig, ax = plt.subplots(ncols=nchan, nrows=2, figsize=(8, 3.2))
    # Image needs to be flipped along uu dimension

    for col=1:nchan

        # Set all labels but the leftmost to blank
        if col != 1
            ax[1,col][:xaxis][:set_ticklabels]([])
            ax[1,col][:yaxis][:set_ticklabels]([])
            ax[2,col][:xaxis][:set_ticklabels]([])
            ax[2,col][:yaxis][:set_ticklabels]([])
        else
            ax[1,col][:set_xlabel](L"uu [k$\lambda$]")
            ax[1,col][:set_ylabel](L"vv [k$\lambda$]")
            ax[1,col][:tick_params](axis="both", which="major", labelsize=8)
            ax[2,col][:set_xlabel](L"uu [k$\lambda$]")
            ax[2,col][:set_ylabel](L"vv [k$\lambda$]")
            ax[2,col][:tick_params](axis="both", which="major", labelsize=8)
        end

        vis_fft = transform(img, col) # Transform this channel
        ext = (vis_fft.uu[end], vis_fft.uu[1], vis_fft.vv[1], vis_fft.vv[end])

        frame = fliplr(vis_fft.VV)


        ax[1,col][:imshow](real(frame), interpolation="none", origin="lower", cmap=plt.get_cmap("bwr"), extent=ext, norm = scale(real(frame)))
        ax[1][:set_title]("Real FFT")
        ax[1][:set_xlabel](L"uu [k$\lambda$]")
        ax[1][:set_ylabel](L"vv [k$\lambda$]")

        ax[2,col][:imshow](imag(frame), interpolation="none", origin="lower", cmap=plt.get_cmap("bwr"), extent=ext, norm = scale(imag(frame)))

        ax[2][:set_title]("Imag FFT")
        ax[2][:set_xlabel](L"uu [k$\lambda$]")
        ax[2][:set_ylabel](L"vv [k$\lambda$]")

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
        fname = outdir * @sprintf("image%04d.out", id)
        im = imread(fname)
        skim = imToSky(im, 73.)
        plot_chmaps(skim, id)
    else
        break
    end
end
