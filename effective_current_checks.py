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
        fig_filename = "{}{}-dist-xz.png".format(prefix, run)
        fig.savefig(fig_filename, dpi=160)
        print("    {}".format(fig_filename))
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
                    cplot = ax.pcolor(zs, xs, j, cmap="seismic",
                                      vmin=-abs_max, vmax=abs_max)
                    # ax.set_aspect("equal")
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    fig.colorbar(cplot, cax=cax)
                    base_text = "{}-{}-{}".format(run, comp, re_or_im)
                    ax.set_title("{}, m={:0.2e}".format(base_text, m))
                    fig_filename = ("{}{}-m{}.png"
                                    "".format(prefix, base_text, m_index))
                    fig.savefig(fig_filename, dpi=160)
                    print("    {}".format(fig_filename))
                    plt.close(fig)


def check_angular_dependence(output_files, prefix=""):
    """
    Check that effective current is independent of polar angle, which
    must be true for any n=0 mode (such as TE011 and TM010)
    """
    out = ec.read_effective_current_output(output_files, prefix)
    print("\ntotal variation over phi:")
    for run in out:
        masses = out[run]["m"][:,0,0,0]
        for m_index, m in enumerate(masses):
            angles = out[run]["phi"][m_index, 0, :, 0]
            avg_dsq = 0
            for re_or_im in ["re", "im"]:
                for comp in ["jr", "jphi", "jz"]:
                    full_case = "{}_{}".format(re_or_im, comp)
                    j0 = out[run][full_case][m_index, :, 0, :]
                    finite = np.isfinite(j0)
                    for phi_index in range(len(angles)):
                        ji = out[run][full_case][m_index, :, phi_index, :]
                        d = j0[finite] - ji[finite]
                        avg_dsq += np.sum(d**2)/np.sum(finite)/6.0
            print("  [{}, m={:0.2e}] average difference squared: "
                  "{:0.2e}".format(run, m, avg_dsq))
    return


def check_parity(output_files, prefix="", tol=1e-10):
    """
    check if the effective current is even or odd in z about the
    center of the cavity. It must be either even or odd for cylindrical
    cavities, as the cavity modes themselves are either even or odd.
    """
    out = ec.read_effective_current_output(output_files, prefix)
    unique_setups = []
    print("\nprocessing runs:")
    for run in out:
        print(run)
        setup = run.split("-")[1]
        if setup not in unique_setups:
            unique_setups.append(setup)
            print("    {}".format(setup))
    for setup in unique_setups:
        left = "left-{}".format(setup)
        right = "right-{}".format(setup)
        masses = out[left]["m"][:,0,0,0]
        print("\n{} using".format(setup))
        for m_index, m in enumerate(masses):
            print("  m = {:0.2e}".format(m))
        center = -out[left]["length"]/2.0
        sample_midpoints = 0.5*( out[left]["z"][0,0,0,:] +
                                out[right]["z"][0,0,0,:])
        print("\n  verify sample points are even about cavity center:\n"
              "    max|(z_left + z_right)/2 - center| = {:0.2e}".format(
               np.abs(sample_midpoints - center).max()))
        print("\n  check that effective current is even or odd")
        for re_or_im in ["re", "im"]:
            for comp in ["jr", "jphi", "jz"]:
                full_case = "{}_{}".format(re_or_im, comp)
                check_even = np.abs(
                    out[left][full_case] - out[right][full_case]).max()
                check_odd  = np.abs(
                    out[left][full_case] + out[right][full_case]).max()
                print("    [{}]\n"
                      "        [check even]    "
                      "max(|j_left - j_right|) = {:0.2e}\n"
                      "        [check odd]     "
                      "max(|j_left + j_right|) = {:0.2e}"
                      "".format(full_case, check_odd, check_even))
                if (check_odd < tol and check_even < tol):
                    print("        {} is zero".format(full_case))
                elif (check_odd < tol and check_even > tol):
                    print("        {} is odd".format(full_case))
                elif (check_odd > tol and check_even < tol):
                    print("        {} is even".format(full_case))


if __name__ == "__main__":

    nearfield_tests = ["check-nearfield-TM010.in.out",
                       "check-nearfield-TE011.in.out",]
    farfield_tests = ["check-farfield-TM010.in.out",
                      "check-farfield-TE011.in.out",]
    parity_tests = ["testparity-left-TM010.in.out",
                    "testparity-left-TE011.in.out",
                    "testparity-right-TM010.in.out",
                    "testparity-right-TE011.in.out"]

    # print("\n" + "-"*40 +"\nplotting j_eff...")
    # plot_overviews(nearfield_tests,
    #                prefix="check-nearfield-")
    # plot_overviews(farfield_tests,
    #                prefix="check-farfield-")

    # print("\n" + "-"*40 +"\ncheck for phi-independence")
    # check_angular_dependence(nearfield_tests, prefix="check-nearfield-")

    # print("\n" + "-"*40 +"\ncheck for parity about center of cavity")
    # check_parity(parity_tests, prefix="testparity-")


    panckake_test = ["check-nearfield-TE011pancake.in.out"]
    print("\n" + "-"*40 +"\nplotting j_eff...")
    plot_overviews(panckake_test,
                   prefix="check-nearfield-")
