# Make some diagnostic plots to show what the model looks like in analytic terms

using model
import PyPlot.plt
using LaTeXStrings
using constants

rr = rs ./ AU # convert to AU

# velocity structure
function plot_vel(pars::Parameters)
    vels = model.velocity(rs, pars) .* 1e-5 # convert from cm/s to km/s

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:plot](rr, vels)
    ax[:set_ylabel](L"$v_\phi$ [km/s]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("plots/velocity.png")
end

# temperature structure
function plot_temp(pars::Parameters)
    temps = model.temperature(rs, pars)

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:plot](rr, temps)
    ax[:set_ylabel](L"$T$ [K]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("plots/temperature.png")
end

# scale height 
function plot_height(pars::Parameters)
    heights = model.Hp(rs, pars) ./ AU

    fig = plt.figure()
    ax = fig[:add_subplot](111)
    ax[:plot](rr, heights)
    ax[:set_ylabel](L"$H_p$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")
    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.85)

    plt.savefig("plots/scale_height.png")
end

# density structure 
function plot_dens(pars::Parameters)

    #Plot the contoured density in cylindrical coordinates, then plot the spherical grid on top of it?
    0
end

plot_vel(params)
plot_temp(params)
plot_height(params)
