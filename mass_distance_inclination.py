# Plot showing the degeneracies between our inferences

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# Load the M_star, d, and incl samples
mdi_samples = np.load("mdi.npy")

# Functon lifted from triangle.py: https://github.com/dfm/triangle.py/
def hist2d(ax, x, y, nxbins=50, nybins=50, *args, **kwargs):
    """
    Plot a 2-D histogram of samples.
    """

    extent = [[x.min(), x.max()], [y.min(), y.max()]]
    bins = 50
    color = "k"
    linewidths = 0.8

    cmap = cm.get_cmap("gray")
    cmap._init()
    cmap._lut[:-3, :-1] = 0.
    cmap._lut[:-3, -1] = np.linspace(1, 0, cmap.N)

    X = np.linspace(extent[0][0], extent[0][1], nxbins + 1)
    Y = np.linspace(extent[1][0], extent[1][1], nybins + 1)
    try:
        H, X, Y = np.histogram2d(x.flatten(), y.flatten(), bins=(X, Y))
    except ValueError:
        raise ValueError("It looks like at least one of your sample columns "
                         "have no dynamic range. You could try using the "
                         "`extent` argument.")

    V = 1.0 - np.exp(-0.5 * np.array([1.0, 2.0, 3.0]) ** 2)
    #V = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)
    Hflat = H.flatten()
    inds = np.argsort(Hflat)[::-1]
    Hflat = Hflat[inds]
    sm = np.cumsum(Hflat)
    sm /= sm[-1]

    for i, v0 in enumerate(V):
        try:
            V[i] = Hflat[sm <= v0][-1]
        except:
            V[i] = Hflat[0]

    X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])
    X, Y = X[:-1], Y[:-1]

    # Plot the contours
    ax.pcolor(X, Y, H.max() - H.T, cmap=cmap)
    ax.contour(X1, Y1, H.T, V, colors=color, linewidths=linewidths)

    data = np.vstack([x, y])
    mu = np.mean(data, axis=1)
    cov = np.cov(data)
    if kwargs.pop("plot_ellipse", False):
        error_ellipse(mu, cov, ax=ax, edgecolor="r", ls="dashed")

    ax.set_xlim(extent[0])
    ax.set_ylim(extent[1])


M_stars = np.linspace(2.0, 3.2, num=200)
incls = np.linspace(100, 125, num=200)

deg = np.pi/180. # [radians]. Convert from degrees to radians by mult. E.g. 110 * deg

czek = {"Mstar":[2.485, 0.1], "i":[109.35, 0.5], "d":[143.45,5.7]}
# Czekala et al measure:
# M_star = 2.485 +/- 0.1
# incl = 141.1 +/- 0.30
# d = 143.45 +/- 5.7

anth = {"Mstar":[2.8, 0.11], "i":[115., 3.], "d":[141., 7.]}
# Anthonioz et al measure:
# M_star = 2.80 +/- 0.11
# i_star = 115 +/- 3
# i_disk = 121 +/- 8
# d = 141 +/- 7



# The two actual measurements with error bars we can put up are ours and Anthonioz.

fig,ax = plt.subplots(nrows=2, figsize=(3.5, 5.5), sharex=True)

hist2d(ax[0], mdi_samples[:,0], mdi_samples[:,2], nxbins=50, nybins=30)

C_incl_disk = czek["Mstar"][0] * np.sin(czek["i"][0] * deg)**2
M_star_disks = C_incl_disk / np.sin(incls * deg)**2

C_incl_RV = 2.114 # +/- 0.01 Alencar 03
M_star_RVs = C_incl_RV / np.sin(incls * deg)**3

M_star_RVs_upper = (C_incl_RV + 0.01)/ np.sin(incls * deg)**3
M_star_RVs_lower = (C_incl_RV - 0.01)/ np.sin(incls * deg)**3

C_incl_ast = anth["Mstar"][0] * np.cos(-anth["i"][0] * deg)**3
M_star_asts = C_incl_ast / np.cos(-incls * deg)**3


ax[0].fill_betweenx(incls, M_star_RVs_lower, M_star_RVs_upper, color="0.6")

ax[0].plot(M_star_disks, incls, "k:")
ax[0].annotate("disk", (2.8, 119.1), rotation=60, size=8)

ax[0].plot(M_star_RVs, incls, "k-.")
ax[0].annotate("RV", (3.0, 118.7), rotation=42, size=8)

ax[0].plot(M_star_asts, incls, "k--")
ax[0].annotate("astrometry", (2.2, 117.5), rotation=-19, size=8)

m, m_err = czek["Mstar"]
incl, incl_err = czek["i"]
ax[0].errorbar(m, incl, xerr=m_err, yerr=incl_err, fmt="bo", ecolor="blue")

m, m_err = anth["Mstar"]
incl, incl_err = anth["i"]
ax[0].errorbar(m, incl, xerr=m_err, yerr=incl_err, color="0.2", fmt="o", ecolor="0.2")

ax[0].annotate("Anthonioz+15", (2.83, 113), size=7)

ax[0].set_ylabel("inclination [degrees]", fontsize=10)
ax[0].set_ylim(105, 120)


# For the second plot:

C_dist_disk = czek["Mstar"][0] / czek["d"][0]
d_dist_disks = M_stars / C_dist_disk

C_disk_ast = anth["Mstar"][0] / anth["d"][0]**3
d_dist_ast = (M_stars/C_disk_ast)**(1/3.)

# Our disk measurement histogram
hist2d(ax[1], mdi_samples[:,0], mdi_samples[:,1])

ax[1].plot(M_stars, d_dist_disks, "k:")
ax[1].annotate("disk", (2.755, 164.5), rotation=60, size=8)

ax[1].plot(M_stars, d_dist_ast, "k--")
ax[1].annotate("astrometry", (2.9, 148), rotation=25, size=8)

# Error bar point
m, m_err = czek["Mstar"]
d, d_err = czek["d"]
ax[1].errorbar(m, d, xerr=m_err, yerr=d_err, fmt="bo", ecolor="blue")

m, m_err = anth["Mstar"]
d, d_err = anth["d"]
ax[1].errorbar(m, d, xerr=m_err, yerr=d_err, color="0.2", fmt="o", ecolor="0.2")

ax[1].set_xlabel(r"$M_\ast$ [$M_\odot$]", fontsize=10)
ax[1].set_ylabel(r"$d$ [pc]", fontsize=10)

ax[1].set_xlim(2.0, 3.2)
ax[1].set_ylim(120, 167)


fig.subplots_adjust(left=0.15, right=0.85, hspace=0.07, top=0.97, bottom=0.07)
fig.savefig("plots/AKSco/mdi.pdf")
