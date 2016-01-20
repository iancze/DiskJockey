#!/usr/bin/env julia

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

# Make some diagnostic plots to show what the model looks like in analytic terms

using JudithExcalibur.model
using JudithExcalibur.constants

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]
pars = convert_dict(config["parameters"], config["model"])

import PyPlot.plt
using LaTeXStrings

# velocity structure
function plot_vel(pars::Parameters, grid)
    vels = JudithExcalibur.model.velocity(grid.rs, pars) .* 1e-5 # convert from cm/s to km/s

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)
    ax[:semilogx](rr, vels)

    # Now, go overlay small grey lines vertically
    for cell_edge in grid.Rs/AU
        ax[:axvline](cell_edge, color="0.5", lw=0.4)
    end

    ax[:set_ylabel](L"$v_\phi$ [km/s]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig]("velocity.png")
end

# temperature structure
function plot_temp(pars::Parameters, grid)
    temps = JudithExcalibur.model.temperature(grid.rs, pars)

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)
    ax[:semilogx](rr, temps)

    # Now, go overlay small grey lines vertically
    for cell_edge in grid.Rs/AU
        ax[:axvline](cell_edge, color="0.5", lw=0.4)
    end

    ax[:set_ylabel](L"$T$ [K]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig]("temperature.png")
end

# scale height
function plot_height(pars::Parameters, grid)
    heights = JudithExcalibur.model.Hp(grid.rs, pars) ./ AU

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)
    ax[:semilogx](rr, heights)
    ax[:set_ylabel](L"$H_p$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig]("scale_height.png")
end

function plot_dlnrho(pars::Parameters)
    nz = 50
    zs = linspace(0, 200 * AU, nz)
    R = 10 * AU
    dlnrhos = Array(Float64, nz)
    for i=1:nz
        dlnrhos[i] = - G * pars.M_star * M_sun * zs[i] * mu_gas * m_H / ((R^2 + zs[i]^2)^1.5 * kB * JudithExcalibur.model.temperature(R, pars))
    end

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    ax[:plot](zs/AU, dlnrhos)

    ax[:set_ylabel](L"$d \ln \rho$")
    ax[:set_xlabel](L"$z$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    fig[:savefig]("dlnrhos.png")
end

function plot_density_1D(pars::Parameters)
    nz = 50
    zs = linspace(0, 30 * AU, nz)
    println("zs ", zs/AU)
    R = 100 * AU
    rhos = Array(Float64, nz)
    for i=1:nz
        rhos[i] = JudithExcalibur.model.rho_gas(R, zs[i], pars)
    end

    lngas = log(rhos) - log(mu_gas * m_H)

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    ax[:plot](zs/AU, lngas)

    # ax[:set_ylabel](L"$\log_{10} n_\textrm{gas} \quad [\log_{10}$ 1/cm^3]")
    ax[:set_ylabel](L"$\ln n_\textrm{gas}$")
    ax[:set_xlabel](L"$z$ [AU]")
    ax[:set_xlim](0, 30)
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    fig[:savefig]("ln_gas_slice.png")
end

function plot_surface_density(pars::Parameters, grid)

    # Sigmas = JudithExcalibur.model.Sigma(grid.rs, pars)

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    Sigmas = Array(Float64, grid.nr)

    for i=1:grid.nr
        Sigmas[i] = JudithExcalibur.model.Sigma(grid.rs[i], pars)
    end

    ax[:loglog](rr, Sigmas)

    # Now, go overlay small grey lines vertically
    for cell_edge in grid.Rs/AU
        ax[:axvline](cell_edge, color="0.5", lw=0.4)
    end

    ax[:set_ylabel](L"$\Sigma\, [\textrm{g/cm}^2]$")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig]("surface_density.png")
end

# density structure
function plot_dens(pars::Parameters, grid)

    # Instead of spherical coordinates, do this with cartesian
    nz = 64
    zs = linspace(0.0, 420 * AU, nz)
    zz = zs./AU

    nr = grid.nr
    rs = grid.rs

    xx = Array(Float64, (nz, nr))
    yy = Array(Float64, (nz, nr))
    rhos = Array(Float64, (nz, nr))

    for i=1:nz
        xx[i, :] = rr
    end

    for j=1:nr
        yy[:, j] = zz
    end

    for i=1:nz
        for j=1:nr
            rhos[i,j] = model.rho_gas(grid.rs[j], zs[i], pars)
        end
    end

    nlog = log10(rhos/(mu_gas * m_H))

    levels = Float64[0.0, 1.0, 2.0, 3.0, 4.0, 5, 6, 7, 8, 9]

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    ax[:set_ylabel](L"$z$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")



    #ticks = np.linspace(0, np.max(cov), num=6)
    #cb.set_ticks(ticks)
    #cb.set_ticks(MaxNLocator(nbins=5))


    #Plot the contoured density in cylindrical coordinates, then plot the spherical grid on top of it?
    # Do this by plotting a bunch of lines
    # First, the radial lines

    # Basically, each of these originates from x = 0, z = 0, and has a slope of theta
    # for theta in (pi/2 - grid.Thetas)
    #     slope = tan(theta)
    #     ax[:plot](rr, slope .* rr, "k", lw=0.1)
    # end
    #
    # # Arcs
    # xs = Array(Float64, (grid.ntheta + 1))
    # ys = Array(Float64, (grid.ntheta + 1))
    #
    # for r in rr
    #     for (i, theta) in enumerate(pi/2 - grid.Thetas)
    #         xs[i] = cos(theta) * r
    #         ys[i] = sin(theta) * r
    #     end
    #     ax[:plot](xs, ys, "k", lw=0.1)
    # end
    #
    # ax[:set_ylim](0, maximum(zz))


    img = ax[:contourf](xx, yy, nlog, levels=levels)

    ax[:set_xscale]("log")

    # Now, go overlay small grey lines vertically
    for cell_edge in grid.Rs/AU
        ax[:axvline](cell_edge, color="0.5", lw=0.4)
    end

    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.77)

    cax = fig[:add_axes]([0.82, 0.22, 0.03, 0.65])
    cb = fig[:colorbar](img, cax=cax)

    plt[:savefig]("density.png")

end

pp = config["parameters"]
params = ["M_star", "r_c", "r_in", "r_cav", "logdelta", "T_10", "q", "gamma", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
starting_param = Array(Float64, nparam)

for i=1:nparam
    starting_param[i] = pp[params[i]][1]
end


starting_param[5] = 10^starting_param[5]
starting_param[9] = 10^starting_param[9]

pars = Parameters(starting_param...)

grd = config["grid"]


grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)

# The cell centers
const global rr = grid.rs ./ AU # convert to AU

plot_vel(pars, grid)
plot_temp(pars, grid)
plot_height(pars, grid)
plot_dens(pars, grid)
plot_density_1D(pars)
plot_surface_density(pars, grid)
plot_dlnrho(pars)
