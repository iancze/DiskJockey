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

imgdir = scratchdir * @sprintf("img%02d/", run_index)

# make the directory where RADMC will make all its files exist
if !ispath(imgdir)
    println("Creating ", imgdir)
    mkdir(imgdir)
end

cd(imgdir)
println("Now in ", pwd())

files = ["radmc3d.inp", "lines.inp", "molecule_co.inp", "wavelength_micron.inp"]
for file in files
    src = homedir * file
    dst = scratchdir * file
    cp(src, dst)
end

println("Copied all RADMC files")

using constants
using image
using model


function make_image(pars, id::Int)

    write_model(pars, "", grid)

    # RADMC conventions for inclination and PA
    incl = pars.incl # [deg]
    PA = pars.PA # [deg] Position angle runs counter clockwise

    run(`radmc3d image incl $incl posang $PA npix $npix loadlambda` |> DevNull)
    println("RADMC finished")

    src = "image.out"
    dst = outdir * @sprintf("image%04d.out", id)
    println("Copying from $src to $dst")

    cp(src, dst)
    println("copied image to outdir")

end

const global npix = 256 # number of pixels
const global grid = Grid(100, 32, 0.5, 800, true)
write_grid("", grid)
write_lambda(lams)

start = run_index * nframes_per_proc + 1
ids = Int[i for i=start:(start + nframes_per_proc)]

for i in ids
    if i <= nframes
        make_image(pars[i], i)
    else
        break
    end
end
