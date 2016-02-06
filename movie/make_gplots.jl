push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/DiskJockey/movie/")
# push!(LOAD_PATH, "/n/home07/iczekala/DiskJockey/")
# push!(LOAD_PATH, "/pool/scout0/DiskJockey/")

using consonance

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    # "--opt1"
    # help = "an option with an argument"
    "--run_index", "-r"
    help = "Output run index"
    default = 0
    arg_type = Int
end

parsed_args = parse_args(ARGS, s)
run_index = parsed_args["run_index"]

imgdir = moviedir * @sprintf("g%02d/", run_index)

using DiskJockey.constants
using DiskJockey.image

import YAML
config = YAML.load(open("config.yaml"))
dpc = config["parameters"]["dpc"][1]

import PyPlot.plt
using LaTeXStrings


global norm = plt[:Normalize](vmin=log10(vmax_g) - 4, vmax=log10(vmax_g), clip=false)

# Plot the channel maps using sky convention
function plot_chmaps(img::image.SkyImage, id::Int)

    fig, ax = plt[:subplots](ncols=nchan, figsize=(6.4, 1.6))
    # Image needs to be flipped along RA dimension
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    for col=1:nchan

        ax[col][:xaxis][:set_ticklabels]([])
        ax[col][:yaxis][:set_ticklabels]([])

        # Set all labels but the leftmost to blank
        # if col != 1
        #     ax[col][:xaxis][:set_ticklabels]([])
        #     ax[col][:yaxis][:set_ticklabels]([])
        # else
        #     ax[col][:set_xlabel](L"$\Delta \alpha$ ('')", size=12)
        #     ax[col][:set_ylabel](L"$\Delta \delta$ ('')", size=12)
        #     ax[col][:tick_params](axis="both", which="major", labelsize=10)
        # end

        #Flip the frame for Sky convention
        frame = flipdim(img.data[:,:,col], 2)
        frame += 1e-10 #Add a tiny bit so that we don't have log10(0)
        lframe = log10(frame)
        max = maximum(lframe)
        ax[col][:imshow](lframe, extent=ext, norm=norm, interpolation="none", origin="lower", cmap=plt[:get_cmap]("PuBu"))
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
    label_incl = L"$i$: " * @sprintf("%.0f", incl) * L"${}^\circ$"
    label_mass = L"$M_\ast$: " * @sprintf("%.2f", mass) * L" $M_\odot$"
    fig[:text](0.01, 0.87, label_mass, size=16)
    fig[:text](0.25, 0.87, label_incl, size=16)

    fig[:subplots_adjust](wspace=0.08, top=0.9, bottom=0.01, left=0.02, right=0.98)

    #plt.savefig(outdir * @sprintf("g%04d.png", id), dpi=150)
    plt[:savefig](imgdir * @sprintf("g%04d.svg", id), transparent=true)
    plt[:close]("all")

end


# How many frames per process?
start = run_index * nframes_per_proc + 1

ids = Int[i for i=start:(start + nframes_per_proc)]

for id in ids
    if id <= nframes
        fname = imgdir * @sprintf("gimage%04d.out", id)
        im = imread(fname)
        skim = imToSky(im, dpc)
        plot_chmaps(skim, id)
    else
        break
    end
end
