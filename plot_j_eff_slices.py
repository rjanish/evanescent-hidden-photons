#!/usr/bin/env python
""" read and process output files of sample_effective_current """

import numpy as np
import scipy.special as spec

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as tick

import parse_sample_effective_current as ec


def Ez_TM010(r, z):
    out = np.zeros(r.shape)
    x01 = spec.jn_zeros(0, 1)
    return spec.j0(x01*r)


def Ephi_TE011(z, r):
    out = np.zeros(r.shape)
    x11 = spec.jn_zeros(1, 1)
    return spec.j1(x11*r)*np.sin(np.pi*z)


def make_mode_slice_plots(output_files, prefix=""):
    """
    plot slices of effective current for two modes,
    to be used for comparison in paper
    """
    out = ec.read_OLD_effective_current_output(output_files, prefix)
    fig, [left, right] = plt.subplots(1,2)

    # color plot of |j_eff| for TM010
    run = "TM010-m10"
    xs = out[run]["x"][:, 0, :]
    zs = out[run]["z"][:, 0, :]
    j = out[run]["jtotal"][:, 0, :]
    abs_max = np.abs(j[np.isfinite(j)]).max()
    cplot = left.pcolor(zs, xs, j/abs_max, cmap="Reds", vmin=0, vmax=1)
    left.set_aspect("equal")

    # j_eff stream lines for TM010
    lw = 1.3
    starts = [[0.06, 0.8],
              [0.15, 1],
              [0.35, 1],
              [0.2, 1.49]]
    stream = left.streamplot(
              out[run]["z"][0, 0, :], out[run]["x"][:, 0, 0],
              out[run]["re_jz"][:, 0, :],
              out[run]["re_jr"][:, 0, :],
              broken_streamlines=False,
              color="maroon", linewidth=lw,
              start_points=starts, arrowsize=2,
              arrowstyle="->")

    # mode lines and crosses for TM010
    Ni=200
    xi = np.linspace(-0.5, 0, Ni)
    yi = np.linspace(0.0, 1, Ni)
    r_left_inside, z_left_inside = np.meshgrid(xi, yi)
    Nlines=5
    Efield_lines = -0.02*np.ones((Nlines,2))
    Efield_lines[:, 1] = np.linspace(0,1, 2*Nlines+1)[1:-1:2]

    stream = left.streamplot(
                  xi, yi, np.ones((Ni, Ni)), np.zeros((Ni, Ni)),
                  broken_streamlines=False, linewidth=lw, arrowsize=2, cmap="Blues",
                  arrowstyle="->", start_points=Efield_lines,
                  color=Ez_TM010(z_left_inside, r_left_inside))
    left.text(-0.35, 0.17, r'$\displaystyle \vec{E}_\textrm{em} $',
            color=cm.Blues(0.8), fontsize=16)
    left.text(-0.18, 1.25, r'$\displaystyle \vec{j}_\textrm{eff} $',
            color='maroon', fontsize=16)

    # cavity boundary
    left.plot([-0.5, 0, 0], [1, 1, 0], marker='',
        linestyle='-', color='black',
        linewidth=2)
    # labels
    left.set_title(r"$ \displaystyle {\rm TM}_{010}$ {\rm Emitter}",
                   fontsize=14)
    left.set_xlabel(r'$ \displaystyle z/R$', size=15, labelpad=5)
    # left.set_ylabel(r'$ \displaystyle \frac{\rho}{R}$', size=14, rotation=0)
    left.text(-0.75, 0.725, r'$ \displaystyle \frac{\rho}{R}$', size=15, rotation=0)
    left.set_ylim([0, 1.5])
    left.set_yticks(0.25*np.arange(7))
    left.set_xlim([-0.5, 0.5])
    left.set_xticks(-0.4 + 0.2*np.arange(5))
    left.xaxis.set_minor_locator(tick.MultipleLocator(0.1))
    left.yaxis.set_minor_locator(tick.MultipleLocator(0.05))
    left.xaxis.set_tick_params(labelsize=13)
    left.yaxis.set_tick_params(labelsize=13)


    # color plot of j_eff_phi for TE011
    run = "TE011-m10"
    xs = out[run]["x"][:, 0, :]
    zs = out[run]["z"][:, 0, :]
    j = out[run]["re_jphi"][:, 0, :]
    abs_max = np.abs(j[np.isfinite(j)]).max()
    cplot = right.pcolor(zs, xs, j/abs_max, cmap="Reds", vmin=0, vmax=1)

    Nrows = 5
    Ncolumns = 2
    Efield_x_2x2, Efield_y_2x2 = np.meshgrid(np.linspace(-0.5, 0, 2*Ncolumns + 1)[1:-1:2],
                                             np.linspace(0, 1, 2*Nrows + 1)[1:-1:2])
    Efield_x = Efield_x_2x2.flatten()
    Efield_y = Efield_y_2x2.flatten()
    absE = np.absolute(Ephi_TE011(Efield_x, Efield_y))
    right.scatter(Efield_x, Efield_y,
                  marker='.', s=15, linewidth=lw,
                  c=absE, cmap="Blues")
    norm = mpl.colors.Normalize(vmin=absE.min(), vmax=absE.max())
    for xp, yp, Ep in zip(Efield_x, Efield_y, absE):
        right.add_patch(plt.Circle((xp, yp), 0.035,
                                   color=cm.Blues(norm(Ep)),
                                   fill=False, linewidth=lw))

    # j_circ_out = [[ 0.25,  1.5*1/8],
    #               [ 0.25,  1.5*3/8],
    #               [ 0.25,  1.5*5/8],
    #               [ 0.25,  1.5*7/8],
                  # [-0.25,  1.5*7/8]]
    Nrows_out = 6
    Ncolumns_out = 2
    jout_x_2x2, jout_y_2x2 = np.meshgrid(np.linspace(-0.5, 0.5, 2*Ncolumns_out + 1)[1:-1:2],
                                         np.linspace(0, 1.5, 2*Nrows_out + 1)[1:-1:2])
    for xp, yp in zip(jout_x_2x2.flatten(), jout_y_2x2.flatten()):
        if not (xp < 0 and yp < 1):
            right.scatter([xp], [yp], marker='.', color="maroon", s=15, linewidth=lw)
            right.add_patch(plt.Circle((xp, yp), 0.035, color="maroon",
                                    fill=False, linewidth=lw))
    right.text(-0.35, 0.17, r'$\displaystyle \vec{E}_\textrm{em} $',
            color=cm.Blues(0.8), size=16)
    right.text(-0.18, 1.25, r'$\displaystyle \vec{j}_\textrm{eff} $',
            color='maroon', size=16)
    right.set_aspect("equal")

    right.set_title(r"$ \displaystyle {\rm TE}_{011}$ {\rm Emitter}",
                    fontsize=14)
    right.set_ylim([0, 1.5])
    right.set_yticks([])
    right.set_xlim([-0.5, 0.5])
    right.set_xticks(-0.4 + 0.2*np.arange(5))
    right.xaxis.set_minor_locator(tick.MultipleLocator(0.1))
    right.xaxis.set_tick_params(labelsize=13)
    right.set_xlabel(r'$z/R$', size=15)
    # cavity boundary
    right.plot([-0.5, 0, 0], [1, 1, 0], marker='',
        linestyle='-', color='black',
        linewidth=2)
    # save
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.03, hspace=None)
    fig.savefig("j_eff-slice.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":

    plt.rcParams['text.usetex'] = True

    paper_plot_runs = ["output-plot-TM010-m10.dat",
                       "output-plot-TE011-m10.dat"]

    make_mode_slice_plots(paper_plot_runs, prefix="output-plot-")



