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

imgdir = scratchdir * @sprintf("d%02d/", run_index)

# make the directory where RADMC will make all its files exist
if !ispath(imgdir)
    println("Creating ", imgdir)
    mkdir(imgdir)
end

cd(imgdir)
println("Now in ", pwd())

files = ["dustopac.inp", "dustkappa_silicate.inp", "stars.inp", "lines.inp", "molecule_co.inp", "wavelength_micron.inp"]
for file in files
    src = homedir * file
    dst = imgdir * file
    cp(src, dst)
end

cp(homedir * "radmc3d.inp.dust", imgdir * "radmc3d.inp")

using constants
using image
using model


function make_image(pars, id::Int)

    write_dust(pars, "", grid)
    write_model(pars, "", grid)

    # RADMC conventions for inclination and PA
    incl = pars.incl # [deg]
    PA = pars.PA # [deg] Position angle runs counter clockwise

    lam0 = cc/230.538e9 * 1e4 # [microns]

    run(`radmc3d mctherm` |> DevNull)
    run(`radmc3d image incl $incl posang $PA npix $npix lambda $lam0` |> DevNull)

    src = "image.out"
    dst = outdir * @sprintf("dimage%04d.out", id)

    cp(src, dst)

end

const global npix = 256 # number of pixels
const global grid = Grid(100, 32, 0.5, 800, true)
write_grid("", grid)

start = run_index * nframes_per_proc + 1
ids = Int[i for i=start:(start + nframes_per_proc)]

for i in ids
    if i <= nframes
        make_image(pars[i], i)
    else
        break
    end
end
