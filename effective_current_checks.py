#!/usr/bin/env python
""" read and process output files of sample_effective_current """

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import parse_sample_effective_current as ec


def plot_overviews(output_files, prefix=""):
    """ plot several coarse-sampled effective currents """
    out = ec.read_effective_current_output(output_files, prefix)
    # check unfolding of coordinates:
    # plot distance in the (z,x>0) plane (ie phi=0) as a function of x and z
    for run in out:
        m = out[run]["m"][0,0,0,0]
        fig, ax = plt.subplots()
        xs = out[run]["x"][0, :, 0, :]
        zs = out[run]["z"][0, :, 0, :]
        r_zx = np.sqrt(zs**2 + xs**2)
        abs_max = np.abs(r_zx[np.isfinite(r_zx)]).max()
        cplot = ax.pcolor(zs, xs, r_zx, cmap="Blues", vmin=0, vmax=abs_max)
            # plot distance in the x-z plane
        title_text = "{}: distance in xz plane".format(run)
        ax.set_title(title_text)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(cplot, cax=cax)
        ax.set_aspect("equal")
        ax.set_title(title_text)
        fig.savefig("{}.png".format(title_text), dpi=160)
        plt.close(fig)
    # check qualitative structure of effective current
    # plot effective currents for phi = 0
    for run in out:
        masses = out[run]["m"][:,0,0,0]
        for m_index, m in enumerate(masses):
            xs = out[run]["x"][m_index, :, 0, :]
            zs = out[run]["z"][m_index, :, 0, :]
            for re_or_im in ["re", "im"]:
                for comp in ["jr", "jphi", "jz"]:
                    fig, ax = plt.subplots()
                    full_case = "{}_{}".format(re_or_im, comp)
                    j = out[run][full_case][m_index, :, 0, :]
                    abs_max = np.abs(j[np.isfinite(j)]).max()
                    cplot = ax.pcolor(zs, xs, j, cmap="seismic", vmin=-abs_max, vmax=abs_max)
                    ax.set_aspect("equal")
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    fig.colorbar(cplot, cax=cax)
                    title_text = "{}-m{}-{}-{}".format(run, m_index, comp, re_or_im)
                    ax.set_title(title_text)
                    fig.savefig("{}.png".format(title_text), dpi=160)
                    plt.close(fig)


def check_angular_dependence(output_files, prefix=""):
    """
    Check that effective current is independent of polar angle, which
    must be true for any n=0 mode (such as TE011 and TM010)
    """
    print("\n" + "-"*40 +"\ncheck for phi-independence")
    out = ec.read_effective_current_output(output_files, prefix)
    print("\ntotal variation over phi:")
    for run in out:
        angles = out[run]["phi"][0, :, 0]
        avg_dsq = 0
        for re_or_im in ["re", "im"]:
            for comp in ["jr", "jphi", "jz"]:
                full_case = "{}_{}".format(re_or_im, comp)
                j0 = out[run][full_case][:, 0, :]
                finite = np.isfinite(j0)
                for phi_index in range(len(angles)):
                    ji = out[run][full_case][:, phi_index, :]
                    d = j0[finite] - ji[finite]
                    avg_dsq += np.sum(d**2)/np.sum(finite)/6.0
        print("  [{}] average difference squared: {:0.2e}".format(run, avg_dsq))
    return


def check_longitudinal_symmetry(output_files, prefix=""):
    """
    check if the effective current is even or odd in z about the
    center of the cavity. It must be either even or odd for cylindrical
    cavities, as the cavity modes themselves are either even or odd.
    """
    print("\n" + "-"*40 +"\ncheck for parity about center of cavity")
    out = ec.read_effective_current_output(output_files, prefix)
    unique_setups = []
    for run in out:
        setup = run.split(sep="-", maxsplit=1)[1]
        if setup not in unique_setups:
            unique_setups.append(setup)
    for setup in unique_setups:
        print("\n"+setup)

        left = "left-{}".format(setup)
        center = -out[left]["length"]/2.0
        right = "right-{}".format(setup)
        sample_midpoints = 0.5*(out[left]["z"][0,0,:] +
                                out[right]["z"][0,0,:])
        print("  check that sample points are even about cavity center:\n"
              "    max|(z_left + z_right)/2 - center| = {:0.2e}".format(
               np.abs(sample_midpoints - center).max()))
        print("  check that effective current is even or odd")
        for re_or_im in ["re", "im"]:
            for comp in ["jr", "jphi", "jz"]:
                full_case = "{}_{}".format(re_or_im, comp)
                diff = out[left][full_case] - out[right][full_case]
                add = out[left][full_case] + out[right][full_case]
                print("    [{}]\n"
                      "    check even:    max(|j_left - j_right|) = {:0.2e}\n"
                      "    check odd:     max(|j_left + j_right|) = {:0.2e}"
                      "".format(full_case, np.abs(diff).max(), np.abs(add).max()))


if __name__ == "__main__":
    crude_samplings = ["crude-TM010.in.out",
                       "crude-TE011.in.out",]
    plot_overviews(crude_samplings, prefix="crude-")
    # check_angular_dependence(crude_samplings, prefix="output-crude-")

    # z_even_check = ["output-testeven-left-TM010-m10.dat",
    #                 "output-testeven-left-TE011-m10.dat",
    #                 "output-testeven-right-TM010-m10.dat",
    #                 "output-testeven-right-TE011-m10.dat"]
    # check_longitudinal_symmetry(z_even_check, prefix="output-testeven-")

