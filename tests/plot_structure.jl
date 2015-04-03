push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    # "--opt1"
    # help = "an option with an argument"
    # default = 0
    "--norad"
    help = "Use the image already here."
    action = :store_true
    "config"
    help = "a YAML configuration file"
    required = true
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

# Make some diagnostic plots to show what the model looks like in analytic terms

using model
import PyPlot.plt
using LaTeXStrings
using constants

# velocity structure
function plot_vel(pars::Parameters, grid)
    vels = model.velocity(grid.rs, pars) .* 1e-5 # convert from cm/s to km/s

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:plot](rr, vels)
    ax[:set_ylabel](L"$v_\phi$ [km/s]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("../plots/velocity.png")
end

# temperature structure
function plot_temp(pars::Parameters, grid)
    temps = model.temperature(grid.rs, pars)

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:plot](rr, temps)
    ax[:set_ylabel](L"$T$ [K]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("../plots/temperature.png")
end

# scale height
function plot_height(pars::Parameters, grid)
    heights = model.Hp(grid.rs, pars) ./ AU

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:plot](rr, heights)
    ax[:set_ylabel](L"$H_p$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("../plots/scale_height.png")
end

# density structure
function plot_dens(pars::Parameters, grid)
    nz = 64
    zs = linspace(0, 300 * AU, nz)
    zz = zs./AU
    nr = length(grid.rs)

    xx = Array(Float64, (nz, nr))
    yy = Array(Float64, (nz, nr))
    nn = Array(Float64, (nz, nr))

    for i=1:nz
        xx[i, :] = rr
    end

    for j=1:nr
        yy[:, j] = zz
    end

    for i=1:nz
        for j=1:nr
            nn[i,j] = model.n_CO(grid.rs[j], zs[i], pars)
        end
    end

    #Add a tiny bit to nn to prevent log10(0.0)
    nlog = log10(nn + 1e-99)
    max = int(maximum(nlog))

    levels = linspace(max - 15, max, 15 + 1)

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:set_ylabel](L"$z$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.77)

    img = ax[:contourf](xx, yy, nlog, levels=levels)

    cax = fig[:add_axes]([0.82, 0.22, 0.03, 0.65])
    cb = fig[:colorbar](img, cax=cax)

    #ticks = np.linspace(0, np.max(cov), num=6)
    #cb.set_ticks(ticks)
    #cb.set_ticks(MaxNLocator(nbins=5))


    #Plot the contoured density in cylindrical coordinates, then plot the spherical grid on top of it?
    # Do this by plotting a bunch of lines
    # First, the radial lines

    # Basically, each of these originates from x = 0, z = 0, and has a slope of theta
    for theta in (pi/2 - grid.Thetas)
        slope = tan(theta)
        ax[:plot](rr, slope .* rr, "k", lw=0.1)
    end

    # Arcs
    xs = Array(Float64, (grid.ntheta + 1))
    ys = Array(Float64, (grid.ntheta + 1))

    for r in rr
        for (i, theta) in enumerate(pi/2 - grid.Thetas)
            xs[i] = cos(theta) * r
            ys[i] = sin(theta) * r
        end
        ax[:plot](xs, ys, "k", lw=0.1)
    end

    ax[:set_ylim](0, maximum(zz))

    plt.savefig("../plots/density.png")

end

pp = config["parameters"]
params = ["M_star", "r_c", "T_10", "q", "gamma", "logM_CO", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
starting_param = Array(Float64, nparam)

for i=1:nparam
    starting_param[i] = pp[params[i]][1]
end

# Convert logM_CO to M_CO
starting_param[6] = 10^starting_param[6]

pars = Parameters(starting_param...)

grd = config["grid"]
grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)
const global rr = grid.rs ./ AU # convert to AU

plot_vel(pars, grid)
plot_temp(pars, grid)
plot_height(pars, grid)
plot_dens(pars, grid)
