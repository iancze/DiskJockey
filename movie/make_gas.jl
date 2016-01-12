# push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")
push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/movie/")
# push!(LOAD_PATH, "/n/home07/iczekala/JudithExcalibur/")
# push!(LOAD_PATH, "/pool/scout0/JudithExcalibur/")

using consonance

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--run_index", "-r"
    help = "Output run index"
    default = 0
    arg_type = Int
end

parsed_args = parse_args(ARGS, s)
run_index = parsed_args["run_index"]

imgdir = moviedir * @sprintf("g%02d/", run_index)
scratchdir = moviedir * "scratch/"

import YAML
config = YAML.load(open("config.yaml"))

# make the directory where RADMC will make all its files exist
if !ispath(imgdir)
    println("Creating ", imgdir)
    mkdir(imgdir)
end

if !ispath(scratchdir)
    mkdir(scratchdir)
end

cd(scratchdir)
println("Now in ", pwd())


files = ["lines.inp", "molecule_co.inp", "wavelength_micron.inp", "radmc3d.inp"]
for file in files
    src = moviedir * file
    dst = scratchdir * file
    cp(src, dst, remove_destination=true)
end


using JudithExcalibur.constants
using JudithExcalibur.image
using JudithExcalibur.model

species = "12CO"

function make_image(pars, id::Int)

    write_model(pars, "", grid, species)

    run(pipeline(`radmc3d image incl $(pars.incl) posang $(pars.PA) npix $npix loadlambda`, DevNull))

    src = "image.out"
    dst = imgdir * @sprintf("gimage%04d.out", id)

    cp(src, dst, remove_destination=true)

end

grd = config["grid"]

const global npix = config["npix"] # number of pixels
const global grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)
write_grid("", grid)
write_lambda(lams, "")

start = run_index * nframes_per_proc + 1
ids = Int[i for i=start:(start + nframes_per_proc)]

for i in ids
    if i <= nframes
        make_image(pars[i], i)
    else
        break
    end
end
