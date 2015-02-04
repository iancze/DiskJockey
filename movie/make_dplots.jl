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

import PyPlot.plt
using LaTeXStrings




# Plot the channel maps using sky convention
function plot_dust(img::image.SkyImage, id::Int)

    fig, ax = plt.subplots(ncols=1, figsize=(2.0, 2.0))
    # Image needs to be flipped along RA dimension
    ext = (img.ra[end], img.ra[1], img.dec[1], img.dec[end])

    ax[:set_xlabel](L"$\Delta \alpha$ ('')", size=8)
    ax[:set_ylabel](L"$\Delta \delta$ ('')", size=8)
    ax[:tick_params](axis="both", which="major", labelsize=8)

    #Flip the frame for Sky convention
    frame = fliplr(img.data[:,:,1])
    frame += 1e-99 #Add a tiny bit so that we don't have log10(0)
    lframe = log10(frame)
    vax = maximum(lframe)
    norm = plt.Normalize(vmin=vmax - 10, vmax=vmax-2, clip=false)
    ax[:imshow](lframe, extent=ext, interpolation="none", origin="lower", norm=norm, cmap=plt.get_cmap("Greys"))

    # Annotate the parameters
    p = pars[id]
    mass = p.M_star
    r_c = p.r_c
    incl = p.incl
    label = L"$M_\ast$: " * @sprintf("%.2f", mass) * L" $M_\odot$   $r_c$: " * @sprintf("%.0f", r_c) * L" AU   $i$: " * @sprintf("%.0f", incl) * L"${}^\circ$"
    fig[:text](0.1, 0.07, label, size=8)

    fig[:subplots_adjust](top=0.95, bottom=0.26, left=0.25, right=0.75)

    plt.savefig(outdir * @sprintf("d%04d.png", id), dpi=300)
    plt.close("all")

end


# How many frames per process?
start = run_index * nframes_per_proc + 1

ids = Int[i for i=start:(start + nframes_per_proc)]

for id in ids
    if id <= nframes
        fname = outdir * @sprintf("dimage%04d.out", id)
        im = imread(fname)
        skim = imToSky(im, 73.)
        plot_dust(skim, id)
    else
        break
    end
end
