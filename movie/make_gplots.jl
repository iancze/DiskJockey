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

import PyPlot.plt
using LaTeXStrings


global norm = plt.Normalize(vmin=log10(vmax_g) - 6, vmax=log10(vmax_g), clip=false)

# Plot the channel maps using sky convention
function plot_chmaps(img::image.SkyImage, id::Int)

    fig, ax = plt.subplots(ncols=nchan, figsize=(6.4, 1.6))
    # Image needs to be flipped along RA dimension
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    for col=1:nchan

        # Set all labels but the leftmost to blank
        if col != 1
            ax[col][:xaxis][:set_ticklabels]([])
            ax[col][:yaxis][:set_ticklabels]([])
        else
            ax[col][:set_xlabel](L"$\Delta \alpha$ ('')", size=12)
            ax[col][:set_ylabel](L"$\Delta \delta$ ('')", size=12)
            ax[col][:tick_params](axis="both", which="major", labelsize=10)
        end

        #Flip the frame for Sky convention
        frame = fliplr(img.data[:,:,col])
        frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
        lframe = log10(frame)
        max = maximum(lframe)
        ax[col][:imshow](lframe, extent=ext, norm=norm, interpolation="none", origin="lower", cmap=plt.get_cmap("PuBu"))
        if col==1
            ax[col][:annotate](@sprintf("%.1f km/s", vels[col]), (0.1, 0.78), xycoords="axes fraction", size=12)
        else
            ax[col][:annotate](@sprintf("%.1f", vels[col]), (0.1, 0.78), xycoords="axes fraction", size=12)
        end

    end

    # Annotate the parameters
    p = pars[id]
    mass = p.M_star
    r_c = p.r_c
    incl = p.incl
    label = L"$i$: " * @sprintf("%.0f", incl) * L"${}^\circ \quad \quad M_\ast$: " * @sprintf("%.2f", mass) * L" $M_\odot$"
    fig[:text](0.35, 0.1, label, size=16)

    fig[:subplots_adjust](wspace=0.08, top=0.97, bottom=0.26, left=0.1, right=0.9)

    plt.savefig(outdir * @sprintf("g%04d.png", id), dpi=150)
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
        plot_chmaps(skim, id)
    else
        break
    end
end
