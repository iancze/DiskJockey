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

using JudithExcalibur.hmodel
using JudithExcalibur.constants
import PyPlot.plt
using LaTeXStrings
using Dierckx

# velocity structure
function plot_vel(pars::Parameters, grid)
    vels = JudithExcalibur.hmodel.velocity(grid.rs, pars) .* 1e-5 # convert from cm/s to km/s

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:plot](rr, vels)
    ax[:set_ylabel](L"$v_\phi$ [km/s]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("velocity.png")
end

# temperature structure
function plot_temp(pars::Parameters, grid)

    nr = grid.nr
    nz = 64
    rs = grid.rs
    zq = JudithExcalibur.hmodel.z_q(rs[end], pars)
    ped = 0.01 * AU
    zs = (logspace(log10(ped), log10(1.5 * zq + ped), nz) - ped)

    # println("rs ,", rs)
    # println("zs ,", zs)

    # z_qs = atm_height(rs);

    xx = Array(Float64, (nz, nr));
    yy = Array(Float64, (nz, nr));
    for i=1:nz
        xx[i, :] = rs
    end

    for j=1:nr
        yy[:, j] = zs
    end

    Ts = Array(Float64, (nz,nr));
    for j=1:nz
        for i=1:nr
            Ts[j, i] = JudithExcalibur.hmodel.temperature(rs[i], zs[j], pars)
        end
    end

    fig = plt.figure(figsize=(6,4))
    ax = fig[:add_subplot](111)

    levels = logspace(0, log10(200), 50);
    img = ax[:contourf](xx./AU, yy./AU, Ts, cmap=plt.get_cmap("Blues_r"), levels=levels)
    # ax[:plot](rs./AU, z_qs./AU, "r--")

    cax = fig[:add_axes]([0.82, 0.22, 0.03, 0.65])
    cb = fig[:colorbar](img, cax=cax)
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.77)

    ax[:set_ylabel](L"$z$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")

    fig[:savefig]("temperature.png")

end

# scale height
function plot_height(pars::Parameters, grid)


    heights = JudithExcalibur.hmodel.Hp(grid.rs, pars) ./ AU

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:plot](rr, heights)
    ax[:set_ylabel](L"$H_p$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("scale_height.png")
end

# Atmospheric height
function plot_atm_height(pars::Parameters, grid)
    heights = JudithExcalibur.hmodel.z_q(grid.rs, pars) ./ AU

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:plot](rr, heights)
    ax[:set_ylabel](L"$z_q$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("atm_height.png")
end

function plot_density_gradient(pars::Parameters)

    nz = 500
    zs = linspace(0, 5 * AU, nz)
    R = 0.5 * AU
    y = Array(Float64, nz)
    for i=1:nz
        y[i] = JudithExcalibur.hmodel.dlnrho(R, zs[i], pars)
    end

    fig = plt.figure()
    ax = fig[:add_subplot](111)

    ax[:plot](zs/AU, y * AU)

    ax[:set_ylabel](L"$\partial \ln \rho_\textrm{gas} / \partial z$ [1/AU]")
    ax[:set_xlabel](L"$z$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("density_gradient.png")

end

function plot_density_1D_unnormed(pars::Parameters)
    nz = 50
    zs = linspace(0, 14 * AU, nz)
    R = 10 * AU
    y = Array(Float64, nz)
    for i=1:nz
        y[i] = JudithExcalibur.hmodel.un_lnrho(R, zs[i], pars)
    end

    fig = plt.figure()
    ax = fig[:add_subplot](111)

    ax[:plot](zs/AU, y)

    ax[:set_ylabel](L"$\ln \rho_\textrm{gas}$ unnormalized")
    ax[:set_xlabel](L"$z$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("un_ln_gas_slice.png")

end

function plot_temperature_1D(pars::Parameters)
    nz = 50
    zs = linspace(0, 200 * AU, nz)
    R = 10 * AU
    temps = Array(Float64, nz)
    for i=1:nz
        temps[i] = JudithExcalibur.hmodel.temperature(R, zs[i], pars)
    end

    fig = plt.figure()
    ax = fig[:add_subplot](111)

    ax[:plot](zs/AU, temps)

    ax[:set_ylabel]("T [K]")
    ax[:set_xlabel](L"$z$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("temperature_1D.png")
end


function plot_dlnrho(pars::Parameters)
    nz = 50
    zs = linspace(0, 200 * AU, nz)
    R = 10 * AU
    dlnrhos = Array(Float64, nz)
    for i=1:nz
        dlnrhos[i] = JudithExcalibur.hmodel.dlnrho(R, zs[i], pars)
    end

    fig = plt.figure()
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
    R = 0.5 * AU
    y = Array(Float64, nz)
    for i=1:nz
        y[i] = JudithExcalibur.hmodel.un_lnrho(R, zs[i], pars)
    end

    # Now calculate the correction factor here
    lncor = exp(JudithExcalibur.hmodel.lncorrection_factor(R, pars))

    lngas = lncor +  y - log(mu_gas * m_H)

    fig = plt.figure()
    ax = fig[:add_subplot](111)

    ax[:plot](zs/AU, lngas)

    ax[:set_ylabel](L"$\ln n_\textrm{gas}$")
    ax[:set_xlabel](L"$z$ [AU]")
    ax[:set_xlim](0,30)
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("ln_gas_slice.png")
end

# Go through some trial points and plot the interpolator to see how it's doing.
function plot_interpolator(pars::Parameters)

    nr = grid.nr
    rs = grid.rs

    lncors = Array(Float64, nr)

    for i=1:nr
        lncors[i] = JudithExcalibur.hmodel.lncorrection_factor(rs[i], pars)
    end

    # Now, let's use the interpolator to make finer samples and see if they line up.
    spl = JudithExcalibur.hmodel.make_correction_interpolator(pars, grid)

    # Now, come up with a finer spaceing of grid points
    n_fine = 100
    rs_fine = logspace(log10(rs[1]), log10(rs[end]), n_fine)
    lncors_fine = Array(Float64, n_fine)
    for i=1:n_fine
        lncors_fine[i] = evaluate(spl, rs_fine[i])
    end


    fig = plt.figure()
    ax = fig[:add_subplot](111)

    ax[:plot](rs/AU, lncors)
    ax[:plot](rs_fine/AU, lncors_fine, "o")

    ax[:set_xlabel](L"$r$ [AU]")
    ax[:set_ylabel]("ln correction factor")

    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    # ax[:set_xlim](0, 20)

    plt.savefig("cor_factor.png")

end

# density structure
function plot_dens(pars::Parameters, grid)

    nr = grid.nr
    nz = 64
    rs = grid.rs
    zq = JudithExcalibur.hmodel.z_q(rs[end], pars)

    zs = cat(1, [0], logspace(log10(0.1 * AU), log10(1.5 * zq), nz-1))

    xx = Array(Float64, (nz, nr));
    yy = Array(Float64, (nz, nr));
    for i=1:nz
        xx[i, :] = rs
    end

    for j=1:nr
        yy[:, j] = zs
    end

    println("Making correction interpolator")
    tic()
    spl = JudithExcalibur.hmodel.make_correction_interpolator(pars, grid)
    toc()
    println("Finished correction interpolator")

    rhos = Array(Float64, (nz,nr))
    for j=1:nz
        for i=1:nr
            temp = JudithExcalibur.hmodel.temperature(rs[i], zs[j], pars)
            XF = JudithExcalibur.hmodel.X_freeze(temp, pars)
            rhos[j, i] = XF * JudithExcalibur.hmodel.rho_gas(rs[i], zs[j], pars, spl)
        end
    end


    nlog = log10(rhos/(mu_gas * m_H))
    #Add a tiny bit to nn to prevent log10(0.0)

    levels = Float64[4.0, 5, 6, 7, 8, 9]
    # levels = linspace(4, 12., 8)

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:set_ylabel](L"$z$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.77)

    img = ax[:contourf](xx./AU, yy./AU, nlog, levels=levels)

    cax = fig[:add_axes]([0.82, 0.22, 0.03, 0.65])
    cb = fig[:colorbar](img, cax=cax)

    println("Minimum ", minimum(nlog))
    println("Maximum ", maximum(nlog))

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

    # ax[:set_ylim](0, maximum(zz))

    plt.savefig("density.png")

end

pp = config["parameters"]
params = ["M_star", "r_c", "T_10m", "q_m", "T_10a", "q_a", "T_freeze", "X_freeze", "gamma", "h", "delta", "logM_gas", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]
nparam = length(params)
starting_param = Array(Float64, nparam)

for i=1:nparam
    starting_param[i] = pp[params[i]][1]
end

# Convert logM_gas to M_gas
starting_param[12] = 10^starting_param[12]

pars = Parameters(starting_param...)

grd = config["grid"]
grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)
const global rr = grid.rs ./ AU # convert to AU

# plot_vel(pars, grid)
# plot_height(pars, grid)
# plot_atm_height(pars, grid)
# plot_temp(pars, grid)
# # plot_density_gradient(pars)
# plot_density_1D_unnormed(pars)
# plot_density_1D(pars)
# plot_interpolator(pars)
plot_dens(pars, grid)
# plot_temperature_1D(pars)
# plot_dlnrho(pars)
