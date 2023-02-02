#!/usr/bin/env python
""" read and process output files of sample_effective_current """

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import parse_sample_effective_current as ec


plt.rcParams['text.usetex'] = True

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
    divider = make_axes_locatable(left)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # fig.colorbar(cplot, cax=cax)
    left.set_title(r"$ \displaystyle {\rm TM}_{010}$ {\rm Source}")

    # j_eff stream lines for TM010
    lw = 0.9
    starts = [[0.06, 0.8],
          # [-0.4, 1.06],
          # [0.05, 1],
          [0.15, 1],
          # [0.25, 1],
          [0.35, 1],
          # [ -0.2, 1.49],
          [0.2, 1.49],
          # [0.4, 1.49],
          ]
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
    # X, Y = np.meshgrid(xi, yi)
    # E_norm = Ez_TM010(Y)/np.max(Ez_TM010(Y))
    stream = left.streamplot(
                  xi, yi, np.ones((Ni, Ni)), np.zeros((Ni, Ni)), #np.zeros(E_norm.shape),
                  broken_streamlines=False,
                  color="navy", linewidth=lw, arrowsize=2,
                  arrowstyle="->",
                  start_points=[[-0.02, 0.05],
                                [-0.02, 0.15],
                                [-0.02, 0.3],
                                [-0.02, 0.6],
                                [-0.02, 0.9]])
    left.text(-0.4, 0.4, r'$\displaystyle \vec{E}_\textrm{source} $',
            color='navy', size=14)
    left.text(-0.2, 1.25, r'$\displaystyle \vec{j}_\textrm{eff} $',
            color='maroon', size=14)
    # circ_l, circ_r = -0.5*2/7, -0.5*5/7
    # circ_x = [circ_l]*3 + [circ_r]*3
    # circ_y = [0.15*1.5, 0.45, 0.75]*2
    # c_circ = 'darkgreen'
    # for xp, yp in zip(circ_x, circ_y):
    #     left.scatter([xp], [yp], marker='x', color=c_circ,
    #                s=6*(yp/0.15)**2, linewidth=2*yp)
    #     left.add_patch(plt.Circle((xp, yp), 0.08*yp, color=c_circ,
    #                             fill=False, linewidth=2*yp))
    # cavity boundary
    left.plot([-0.5, 0, 0], [1, 1, 0], marker='',
        linestyle='-', color='black',
        linewidth=2)
    left.set_xlabel(r'$z/R$', size=12)
    left.set_ylabel(r'$r/R$', size=12)


    # color plot of j_eff_phi for TE011
    run = "TE011-m10"
    xs = out[run]["x"][:, 0, :]
    zs = out[run]["z"][:, 0, :]
    j = out[run]["re_jphi"][:, 0, :]
    abs_max = np.abs(j[np.isfinite(j)]).max()
    cplot = right.pcolor(zs, xs, j/abs_max, cmap="Reds", vmin=0, vmax=1)
    j_cir_in = [[-0.5*1/4,  1.5*1/8],
                [-0.5*1/4,  1.5*3/8],
                [-0.5*3/4,  1.5*1/8],
                [-0.5*3/4,  1.5*3/8]]
    for xp, yp in j_cir_in:
        right.scatter([xp], [yp], marker='x', color="navy", s=80, linewidth=lw)
        right.add_patch(plt.Circle((xp, yp), 0.05, color="navy",
                                fill=False, linewidth=lw))
    j_circ_out = [[ 0.25,  1.5*1/8],
                  [ 0.25,  1.5*3/8],
                  [ 0.25,  1.5*5/8],
                  [ 0.25,  1.5*7/8],
                  [-0.25,  1.5*7/8]]
    for xp, yp in j_circ_out:
        right.scatter([xp], [yp], marker='x', color="maroon", s=80, linewidth=lw)
        right.add_patch(plt.Circle((xp, yp), 0.05, color="maroon",
                                fill=False, linewidth=lw))
    right.text(-0.4, 0.7, r'$\displaystyle \vec{E}_\textrm{source} $',
            color='navy', size=14)
    right.text(-0.175, 1.2, r'$\displaystyle \vec{j}_\textrm{eff} $',
            color='maroon', size=14)
    right.set_aspect("equal")
    # divider = make_axes_locatable(right)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # fig.colorbar(cplot)
    right.set_title(r"$ \displaystyle {\rm TE}_{011}$ {\rm Source}")
    right.set_yticks([])
    right.set_yticks([])
    right.set_xlabel(r'$z/R$', size=12)
    # cavity boundary
    right.plot([-0.5, 0, 0], [1, 1, 0], marker='',
        linestyle='-', color='black',
        linewidth=2)

    fig.savefig("j_eff-slice.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":

    paper_plot_runs = ["output-plot-TM010-m10.dat",
                       "output-plot-TE011-m10.dat"]
    make_mode_slice_plots(paper_plot_runs, prefix="output-plot-")



