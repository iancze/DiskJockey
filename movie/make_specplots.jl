using consonance

push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/DiskJockey/")
push!(LOAD_PATH, "/n/home07/iczekala/DiskJockey/")
push!(LOAD_PATH, "/pool/scout0/DiskJockey/")

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


function spec_plot(img::image.SkyImage, id::Int)

    fig = plt.figure()
    ax = fig[:add_subplot](111)

    spec = imToSpec(img)

    # CO 2-1 rest frame
    lam0 = cc/230.538e9 * 1e4 # [microns]

    # convert wavelengths to velocities
    vels = c_kms * (img.lams .- lam0)/lam0 # [km/s]

    ax[:plot](vels, spec[:,2], ls="steps-mid")
    # ax[:plot](vels, reverse(spec[:,2]), ls="steps-mid")
    ax[:set_ylabel](L"$f_\nu$ [Jy]")
    ax[:set_xlabel](L"$v$ [km/s]")


    # Annotate the parameters
    p = pars[id]
    mass = p.M_star
    r_c = p.r_c
    incl = p.incl
    label = L"$M_\ast$: " * @sprintf("%.2f", mass) * L" $M_\odot$   $r_c$: " * @sprintf("%.0f", r_c) * L" AU   $i$: " * @sprintf("%.0f", incl) * L"${}^\circ$"


    fig[:text](0.25, 0.1, label)


    plt.savefig(outdir * @sprintf("spec%04d.png", id), dpi=150)
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
        spec_plot(skim, id)
    else
        break
    end
end
