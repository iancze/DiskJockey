push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

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

imgdir = @sprintf("img%02d/", run_index)

# make the directory where RADMC will make all its files exist
if !ispath(imgdir)
    println("Creating ", imgdir)
    mkdir(imgdir)
end

cd(imgdir)
println("Now in ", pwd())

files = ["radmc3d.inp", "lines.inp", "molecule_co.inp", "wavelength_micron.inp"]
for file in files
    cp("../../" * file, file)
end

using constants
using image
using model


global const nchan = 7
global const vels = linspace(-1.5, 1.5, nchan) # [km/s]
# CO 2-1 rest frame
lam0 = cc/230.538e9 * 1e4 # [microns]
# convert velocities to wavelengths
lams = lam0 * (vels/c_kms + 1)

function make_image(pars, id::Int)

    write_model(pars, "", grid)
    vel = pars.vel # [km/s]
    # RADMC conventions for inclination and PA
    incl = pars.incl # [deg]
    PA = pars.PA # [deg] Position angle runs counter clockwise

    run(`radmc3d image incl $incl posang $PA npix $npix loadlambda` |> DevNull)

    cp("image.out", @sprintf("../image%04d.out", id))

end

M_star = 1.75 # [M_sun] stellar mass
r_c = 45. # [AU] characteristic radius
T_10 = 115. # [K] temperature at 10 AU
q = 0.63 # temperature gradient exponent
gamma = 1.0 # surface density gradient
logM_CO = 0.2 # [M_earth] disk mass of CO
ksi = 0.14 # [km/s] microturbulence
dpc = 73. # [pc] distance
incl = 147. # [degrees] inclination
PA = -17. # [degrees] position angle
vel = -31.18 # [km/s]
mu_RA = 0.2 # [arcsec] centroid location
mu_DEC = -0.6 # [arcsec]

pars = Parameters(M_star, r_c, T_10, q, gamma, 10^logM_CO, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

const global npix = 256 # number of pixels
const global grid = Grid(100, 32, 0.5, 800, true)
write_grid("", grid)
write_lambda(lams)

# Create a master parameter list
# First, adjust in inclination
# then, adjust in mass

nframes = 4
ids = Int[i for i=1:nframes]
incls = linspace(0, 90., nframes)

for i in ids
    pars.incl = incls[i]
    make_image(pars, i)
end
