push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/DiskJockey/")

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

using emodel
import PyPlot.plt
using LaTeXStrings
using constants

pp = config["parameters"]
params = ["M_star", "a_c", "T_10", "q", "gamma", "logM_CO", "ksi", "dpc", "incl", "PA", "e", "w", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
starting_param = Array{Float64}(nparam)

for i=1:nparam
    starting_param[i] = pp[params[i]][1]
end

# Convert logM_CO to M_CO
starting_param[6] = 10^starting_param[6]
pars = Parameters(starting_param...)

grd = config["grid"]
grid = Grid(grd["nr"], grd["ntheta"], grd["nphi"], grd["r_in"], grd["r_out"], grd["na"], true)

# Plot the rings and ring centers for the disk.
function plot_rings(pars::Parameters, grid)
    fig = plt.figure(figsize=(8,8))
    ax = fig[:add_subplot](111, polar=true)

    phis = linspace(0, 2pi, 100)

    # Plot the grid cell edges in radius and phi
    for R in grid.Rs
        rr = R * ones(phis)
        ax[:plot](phis, rr ./AU, "b", lw=0.1)
    end

    for Phi in grid.Phis
        ax[:plot](Phi * ones(grid.Rs), grid.Rs ./AU, "b", lw=0.1)
    end

    # theta, r
    for a in grid.as
        rr = emodel.radius(phis, a, pars)
        ax[:plot](phis, rr./AU, "k-.")
    end

    for A in grid.As
        RR = emodel.radius(phis, A, pars)
        ax[:plot](phis, RR./AU, "k-")
    end

    # Plot the outermost wall radii
    rout = grid.Rs[end] .* ones(phis)
    ax[:plot](phis, rout./AU, "b-")

    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("../plots/radii_outer.png")
end

function plot_inner_rings(pars::Parameters, grid)
    fig = plt.figure(figsize=(8,8))
    ax = fig[:add_subplot](111, polar=true)

    phis = linspace(0, 2pi, 100)

    # Plot the grid cell edges in radius and phi
    for R in grid.Rs
        rr = R * ones(phis)
        ax[:plot](phis, rr ./AU, "b", lw=0.1)
    end

    for Phi in grid.Phis
        ax[:plot](Phi * ones(grid.Rs), grid.Rs ./AU, "b", lw=0.1)
    end

    # theta, r
    for a in grid.as
        rr = emodel.radius(phis, a, pars)
        ax[:plot](phis, rr./AU, "k-.")
    end

    for A in grid.As
        RR = emodel.radius(phis, A, pars)
        ax[:plot](phis, RR./AU, "k-")
    end

    # Plot the innermost radii
    rin = grid.Rs[1] .* ones(phis)
    ax[:plot](phis, rin./AU, "b-")

    ax[:set_rmax](2)

    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("../plots/radii_inner.png")
end

plot_rings(pars, grid)
plot_inner_rings(pars, grid)
