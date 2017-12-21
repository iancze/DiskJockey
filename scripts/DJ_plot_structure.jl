#!/usr/bin/env julia

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "-M"
    help = "Print out the list of files this program will create."
    action = :store_true
    "config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using DiskJockey.model
using DiskJockey.constants

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]
pars = convert_dict(config["parameters"], config["model"])

if parsed_args["M"]
    args = "velocity.png temperature.png scale_height.png surface_density.png density.png grid_topgrid.png"
    # If we have vertical gradient, add some extra to this
    println(args)
    quit()
end

import PyPlot.plt
using LaTeXStrings

# Plot looking down, in polar coordinates
function plot_topgrid(pars::AbstractParameters, grid::Grid)
    fig = plt[:figure](figsize=(8,8))
    ax = fig[:add_subplot](111, polar=true)

    # Something to span the circle
    phis = linspace(0, 2pi, 100)

    # Plot the grid cell edges in radius and phi
    for R in grid.Rs
        rr = R * ones(phis)
        ax[:plot](phis, rr ./AU, "b", lw=0.1)
    end

    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)
    plt[:savefig]("grid_topgrid.png", dpi=300)
end

function plot_sidegrid(pars::AbstractParameters, grid::Grid)
end

# velocity structure
function plot_vel(pars::AbstractParameters, grid::Grid)
    vels = DiskJockey.model.velocity(grid.rs, pars) .* 1e-5 # convert from cm/s to km/s

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

# Instead of a 1D slice, we'll need to do a 2D field
function plot_vel(pars::ParametersVertical, grid::Grid)
    # Instead of spherical coordinates, do this with cartesian
    nz = 64
    zs = linspace(0.0, 20 * AU, nz)
    zz = zs./AU

    nr = grid.nr
    rs = grid.rs

    xx = Array{Float64}(nz, nr)
    yy = Array{Float64}(nz, nr)
    vels = Array{Float64}(nz, nr)

    for i=1:nz
        xx[i, :] = rr
    end

    for j=1:nr
        yy[:, j] = zz
    end

    for i=1:nz
        for j=1:nr
            vels[i,j] = DiskJockey.model.velocity(grid.rs[j], zs[i], pars)
        end
    end


    levels = Float64[0.0, 1.0, 5.0, 10.0, 20.0, 40.]

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

    # convert from cm/s to km/s
    img = ax[:contourf](xx, yy, 1e-5 * vels, levels=levels)

    ax[:set_xscale]("log")

    # # Now, go overlay small grey lines vertically
    # for cell_edge in grid.Rs/AU
    #     ax[:axvline](cell_edge, color="0.5", lw=0.4)
    # end

    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.77)

    cax = fig[:add_axes]([0.82, 0.22, 0.03, 0.65])
    cb = fig[:colorbar](img, cax=cax)

    plt[:savefig]("velocity.png")
end

# temperature structure
function plot_temp(pars::AbstractParameters, grid::Grid)
    temps = DiskJockey.model.temperature(grid.rs, pars)

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


function plot_temp(pars::Union{ParametersVertical, ParametersVerticalEta}, grid::Grid)
    # Instead of spherical coordinates, do this with cartesian
    nz = 64
    zs = linspace(0.0, 20 * AU, nz)
    zz = zs./AU

    nr = grid.nr
    rs = grid.rs

    xx = Array{Float64}(nz, nr)
    yy = Array{Float64}(nz, nr)
    temps = Array{Float64}(nz, nr)

    for i=1:nz
        xx[i, :] = rr
    end

    for j=1:nr
        yy[:, j] = zz
    end

    for i=1:nz
        for j=1:nr
            temps[i,j] = DiskJockey.model.temperature(grid.rs[j], zs[i], pars)
        end
    end

    levels = Float64[0.0, 5.0, 10.0, 20.0, 30., 40., 50.0, 100.0]

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


    img = ax[:contourf](xx, yy, temps, levels=levels)

    ax[:set_xscale]("log")

    # # Now, go overlay small grey lines vertically
    # for cell_edge in grid.Rs/AU
    #     ax[:axvline](cell_edge, color="0.5", lw=0.4)
    # end

    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.77)

    cax = fig[:add_axes]([0.82, 0.22, 0.03, 0.65])
    cb = fig[:colorbar](img, cax=cax)

    plt[:savefig]("temperature.png")
end

# scale height
function plot_height(pars::AbstractParameters, grid::Grid)
    heights = DiskJockey.model.Hp(grid.rs, pars) ./ AU

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)
    ax[:semilogx](rr, heights)
    ax[:set_ylabel](L"$H_p$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig]("scale_height.png")
end

function plot_surface_density(pars::AbstractParameters, grid::Grid)

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    Sigmas = Array{Float64}(grid.nr)

    for i=1:grid.nr
        Sigmas[i] = DiskJockey.model.Sigma(grid.rs[i], pars)
    end

    # ax[:loglog](rr, Sigmas)
    ax[:semilogy](rr, Sigmas)

    # Now, go overlay small grey lines vertically for the radial cells
    for cell_edge in grid.Rs/AU
        ax[:axvline](cell_edge, color="0.5", lw=0.4)
    end

    ax[:set_ylabel](L"$\Sigma\, [\mathrm{g/cm}^2]$")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig]("surface_density.png")
end



# density structure
function plot_dens(pars::AbstractParameters, grid)

    # Instead of spherical coordinates, do this with cartesian
    nz = 64
    zs = linspace(0.0, 150 * AU, nz)
    zz = zs./AU

    nr = grid.nr
    rs = grid.rs

    xx = Array{Float64}(nz, nr)
    yy = Array{Float64}(nz, nr)
    rhos = Array{Float64}(nz, nr)

    for i=1:nz
        xx[i, :] = rr
    end

    for j=1:nr
        yy[:, j] = zz
    end

    for i=1:nz
        for j=1:nr
            rhos[i,j] = DiskJockey.model.rho_gas(grid.rs[j], zs[i], pars)
        end
    end

    nlog = log10.(rhos/(mu_gas * m_H))

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

# Plot a 1D slice of the vertical density structure at a particular radius
function plot_density_column(pars::ParametersVertical, grid::Grid)
    r = 10 * AU

    # Calculate the normalizing constant and height of dissocation
    norm_rho, z_phot = DiskJockey.model.rho_norm_and_phot(r, pars)

    # Calculate a range of zs
    nz = 64
    zs = linspace(0, 4 * AU, nz)

    un_lnrhos = Array{Float64}(nz)
    for i=1:nz
        un_lnrhos[i] = DiskJockey.model.un_lnrho(r, zs[i], pars)
    end

    # calculate un_rhos for these zs
    rhos = norm_rho .* exp(un_lnrhos)

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    ax[:semilogy](zs./AU, rhos/(mu_gas * m_H))

    # Now, go overlay small grey lines vertically for the radial cells
    # for cell_edge in grid.Rs/AU
    #     ax[:axvline](cell_edge, color="0.5", lw=0.4)
    # end

    ax[:set_ylabel](L"$\rho_\mathrm{gas} \, [\mathrm{n/cm}^3]$")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig]("density_column.png")

end

function plot_density_column_CO(pars::ParametersVertical, grid::Grid)
    r = 50 * AU

    # Gas density [g/cm^3]
    zs, rhos = DiskJockey.model.rho_column(r, pars)

    # Gas density of CO [g/cm^3]
    zs, rhos = DiskJockey.model.rho_column_CO(r, zs, rhos, pars)

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    ax[:semilogy](zs./AU, DiskJockey.constants.f_12CO * DiskJockey.constants.X_H2 * rhos/(mu_gas * m_H))

    # Now, go overlay small grey lines vertically for the radial cells
    # for cell_edge in grid.Rs/AU
    #     ax[:axvline](cell_edge, color="0.5", lw=0.4)
    # end

    ax[:set_ylabel](L"$\rho_\mathrm{CO} \, [\mathrm{n/cm}^3]$")
    ax[:set_xlabel](L"$z$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig]("density_column_CO.png")

end

function test_dens(pars::ParametersVertical, grid::Grid)
    tic()
    rho_gas = DiskJockey.model.rho_gas_interpolator(pars, grid)
    println("Constructing interpolator")
    toc()
    #
    # r_in = grid.rs[1]
    # println("r_in ", r_in)
    # println(rho_gas(r_in, 0.0, pars))

    r2 = 0.12 * AU
    println("r2 ", r2)
    println("un_lnrho ", DiskJockey.model.un_lnrho(r2, 0.0, pars))
    println("norm + zphot ", DiskJockey.model.rho_norm_and_phot(r2, pars))
    tic()
    println(rho_gas(r2, 0.0, pars))
    toc()
end

function plot_norm(pars::ParametersVertical, grid::Grid)

    # DiskJockey.model.rho_norm(10 * AU, pars)
    # DiskJockey.model.z_phot(10 * AU, pars)
    #
    # quit()

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    norms = Array{Float64}(grid.nr)

    for i=1:grid.nr
        norms[i] = DiskJockey.model.rho_norm(grid.rs[i], pars)
    end

    fig, ax = plt[:subplots](figsize=(6,6))

    # ax[:loglog](rr, Sigmas)
    ax[:loglog](rr, norms)

    ax[:set_ylabel](L"$norm$")

    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig]("norm.png")
end

function plot_ztop(pars::ParametersVertical, grid::Grid)

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    ztops = Array{Float64}(grid.nr)

    for i=1:grid.nr
        ztops[i] = DiskJockey.model.z_top(grid.rs[i], pars)
    end

    fig, ax = plt[:subplots](figsize=(6,6))

    # ax[:loglog](rr, Sigmas)
    ax[:semilogx](rr, ztops ./ AU)

    ax[:set_ylabel](L"$z_\mathrm{top}$ [AU]")

    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt[:savefig]("ztop.png")
end

function plot_dens(pars::ParametersVertical, grid::Grid)
    # Set up the interpolator

    DiskJockey.model.rho_gas(10 * AU, 0.1 * AU, pars)

    tic()
    DiskJockey.model.rho_gas(10 * AU, 0.1 * AU, pars)
    println("Typical evaluation time")
    toc()

    # Instead of spherical coordinates, do this with cartesian
    nz = 128
    zs = linspace(0.0, 100 * AU, nz)
    zz = zs./AU

    nr = grid.nr
    rs = grid.rs

    xx = Array{Float64}(nz, nr)
    yy = Array{Float64}(nz, nr)
    rhos = Array{Float64}(nz, nr)

    for i=1:nz
        xx[i, :] = rr
    end

    for j=1:nr
        yy[:, j] = zz
    end

    for i=1:nz
        for j=1:nr
            temp = DiskJockey.model.temperature(grid.rs[j], zs[i], pars)
            rhos[i,j] = DiskJockey.model.X_freeze(temp, pars) * DiskJockey.model.rho_gas(grid.rs[j], zs[i], pars)
        end
    end

    ztops = Array{Float64}(grid.nr)
    for i=1:grid.nr
        ztops[i] = DiskJockey.model.z_top(grid.rs[i], pars)
    end

    println("rho extrema ", extrema(rhos))

    nlog = log10(rhos/(mu_gas * m_H))

    levels = Float64[0.0, 1.0, 2.0, 3.0, 4.0, 5, 6, 7, 8, 9, 10, 11, 12]

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    ax[:set_ylabel](L"$z$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")

    img = ax[:contourf](xx, yy, nlog, levels=levels)

    ax[:plot](rr, ztops ./ AU)

    ax[:set_xscale]("log")
    ax[:set_ylim](0, 100)

    # Now, go overlay small grey lines vertically
    # for cell_edge in grid.Rs/AU
    #     ax[:axvline](cell_edge, color="0.5", lw=0.4)
    # end

    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.77)

    cax = fig[:add_axes]([0.82, 0.22, 0.03, 0.65])
    cb = fig[:colorbar](img, cax=cax)

    plt[:savefig]("density.png")
end


grid = Grid(config["grid"])

# The cell centers for plotting purposes
const global rr = grid.rs ./ AU # convert to AU


if config["model"] == "vertical"
    plot_ztop(pars, grid)
end

plot_topgrid(pars, grid)
plot_vel(pars, grid)
plot_temp(pars, grid)
plot_height(pars, grid)
plot_surface_density(pars, grid)
plot_dens(pars, grid)
