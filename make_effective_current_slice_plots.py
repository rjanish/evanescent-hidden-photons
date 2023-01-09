#!/usr/bin/env python
""" read and process output files of sample_effective_current """

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def read_effective_current_samples(filenames, prefix=""):
    """ read output files of sample_effective_current """

    # param file format assumptions
    param_types = {"r_N"  :int, "r_min"  :float, "r_max"  :float,
                  "phi_N":int, "phi_min":float, "phi_max":float,
                  "z_N"  :int, "z_min"  :float, "z_max":float,
                  "atol":float, "rtol":float,
                  "length":float, "radius":float,
                  "mass":float, "mode":str}
    num_params = len(param_types)
    num_header_lines = num_params + 4
        # assume param header occupies lines 1 to num_params, then
        # data identifer header occupies lines Nparam + 1 to num_params + 3,
        # so data begins on line Nparam + 4
    column_names = ["r", "phi", "z",
                    "re_jr", "error_re_jr", "im_jr", "error_im_jr",
                    "re_jphi", "error_re_jphi", "im_jphi", "error_im_jphi",
                    "re_jz", "error_re_jz", "im_jz", "error_im_jz"]

    # read files
    if isinstance(filenames, str):
        filenames = [filenames] # handle single string filename input
    runs = {}
    for filename in filenames: # assumes a list of filenames
        label = filename.split(sep='.')[0][len(prefix):]
            # remove prefix and file ext
        runs[label] = {}
        # read params
        with open(filename) as file:
            for i in range(num_params):
                key, value = next(file).split()
                runs[label][key] = param_types[key](value)
                    # assume each line of param header contains a single
                    # param name and value separated by whitespae
        # read data
        runs[label]["rawdata"] = np.loadtxt(filename, skiprows=num_header_lines)
        shape_3x3 = (runs[label]["r_N"],
                     runs[label]["phi_N"],
                     runs[label]["z_N"])
                       # this is the order of data columns in the file
        for num, name in enumerate(column_names):
            runs[label][name] = (
                runs[label]["rawdata"][:, num].reshape(shape_3x3))
                # put data in meshgrid format

        # construct other coordinates and values
        runs[label]["x"] = runs[label]["r"]*np.cos(runs[label]["phi"])
        runs[label]["y"] = runs[label]["r"]*np.sin(runs[label]["phi"])
        runs[label]["mag_jr"]   = np.sqrt(
            runs[label]["re_jr"]**2   + runs[label]["im_jr"]**2)
        runs[label]["mag_jphi"] = np.sqrt(
            runs[label]["re_jphi"]**2 + runs[label]["im_jphi"]**2)
        runs[label]["mag_jz"]   = np.sqrt(
            runs[label]["re_jz"]**2   + runs[label]["im_jz"]**2)
        runs[label]["jtotal"]   = np.sqrt(
            runs[label]["im_jr"]**2   + runs[label]["re_jr"]**2 +
            runs[label]["im_jphi"]**2 + runs[label]["re_jphi"]**2 +
            runs[label]["im_jz"]**2   + runs[label]["re_jz"]**2)
    return runs


def plot_overviews(output_files, prefix=""):
    """ plot several coarse-sampled effective currents """
    out = read_effective_current_samples(output_files, prefix)
    # check unfolding of coordinates:
    # plot distance in the (z,x>0) plane (ie phi=0) as a function of x and z
    for run in out:
        fig, ax = plt.subplots()
        xs = out[run]["x"][:, 0, :]
        zs = out[run]["z"][:, 0, :]
        r_zx = np.sqrt(zs**2 + xs**2)
        abs_max = np.abs(r_zx[np.isfinite(r_zx)]).max()
        cplot = ax.pcolor(zs, xs, r_zx, cmap="Blues", vmin=0, vmax=abs_max)
            # plot distance in the x-z plane
        title_text = "{}-r_xz".format(run)
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
        xs = out[run]["x"][:, 0, :]
        zs = out[run]["z"][:, 0, :]
        for re_or_im in ["re", "im"]:
            for comp in ["jr", "jphi", "jz"]:
                fig, ax = plt.subplots()
                full_case = "{}_{}".format(re_or_im, comp)
                j = out[run][full_case][:, 0, :]
                abs_max = np.abs(j[np.isfinite(j)]).max()
                cplot = ax.pcolor(zs, xs, j, cmap="seismic", vmin=-abs_max, vmax=abs_max)
                ax.set_aspect("equal")
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(cplot, cax=cax)
                title_text = "{}-{}-{}".format(run, comp, re_or_im)
                ax.set_title(title_text)
                fig.savefig("{}.png".format(title_text), dpi=160)
                plt.close(fig)


def check_angular_dependence(output_files, prefix=""):
    """
    Check that effective current is independent of polar angle, which
    must be true for any n=0 mode (such as TE011 and TM010)
    """
    out = read_effective_current_samples(output_files, prefix)
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
        print("  [{}] average difference squared: {:e}".format(run, avg_dsq))
    return


def check_longitudinal_symmetry(output_files, prefix=""):
    """
    check if the effective current is even or odd in z about the
    center of the cavity. It must be either even or odd for cylindrical
    cavities, as the cavity modes themselves are either even or odd.
    """
    out = read_effective_current_samples(output_files, prefix)
    unique_setups = []
    for run in out:
        setup = run.split(sep="-", maxsplit=1)[1]
        if setup not in unique_setups:
            unique_setups.append(setup)
    for setup in unique_setups:
        print("\n"+setup)

        left = "left-{}".format(setup)
        print("\n   left z-samples\n    {}\n"
              "".format(out[left]["z"][0,0,:]))
        left_center = -0.5*out[left]["length"]
        print("    relative to cavity center {}\n    {}\n"
              "".format(left_center, out[left]["z"][0,0,:] - left_center))

        right = "right-{}".format(setup)
        print("\n   right z-samples\n    {}\n"
              "".format(out[right]["z"][0,0,:]))
        right_center = -0.5*out[right]["length"]
        print("    relative to cavity center {}\n    {}\n"
              "".format(right_center, out[right]["z"][0,0,:] - right_center))


        for re_or_im in ["re", "im"]:
            for comp in ["jr", "jphi", "jz"]:
                full_case = "{}_{}".format(re_or_im, comp)
                diff = out[left][full_case] - out[right][full_case]
                add = out[left][full_case] + out[right][full_case]
                print("  [{}]\n"
                      "check even:    max(|j_left - j_right|) = {}\n"
                      "check odd:     max(|j_left + j_right|) = {}\n"
                      "".format(full_case, np.abs(diff).max(), np.abs(add).max()))


def make_mode_plots(output_files, prefix=""):
    """
    plot slices of effective current for two modes,
    to be used for comparison in paper
    """
    out = read_effective_current_samples(output_files, prefix)

    # color plot of |j_eff| for TM010
    run = "TM010-m10"
    xs = out[run]["x"][:, 0, :]
    zs = out[run]["z"][:, 0, :]
    fig, ax = plt.subplots()
    j = out[run]["jtotal"][:, 0, :]
    abs_max = np.abs(j[np.isfinite(j)]).max()
    cplot = ax.pcolor(zs, xs, j/abs_max, cmap="Reds", vmin=0, vmax=1)
    ax.set_aspect("equal")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(cplot, cax=cax)
    ax.set_title("Effective current generated by TM010 source")
    fig.savefig("{}.png".format(run), dpi=160)
    plt.close(fig)

    # color plot of j_eff_phi for TE011
    run = "TE011-m10"
    xs = out[run]["x"][:, 0, :]
    zs = out[run]["z"][:, 0, :]
    fig, ax = plt.subplots()
    j = out[run]["re_jphi"][:, 0, :]
    abs_max = np.abs(j[np.isfinite(j)]).max()
    cplot = ax.pcolor(zs, xs, j/abs_max, cmap="Reds", vmin=0, vmax=1)
    ax.set_aspect("equal")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(cplot, cax=cax)
    ax.set_title("Effective current generated by TE011 source")
    fig.savefig("{}.png".format(run), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    # crude_samplings = ["output-crude-TM010evan.dat",
    #                    "output-crude-TM010prop.dat",
    #                    "output-crude-TE011evan.dat",
    #                    "output-crude-TE011prop.dat"]
    # plot_overviews(crude_samplings, prefix="output-crude-")
    # check_angular_dependence(crude_samplings, prefix="output-crude-")

    z_even_check = ["output-testeven-left-TM010-m10.dat",
                    "output-testeven-left-TE011-m10.dat",
                    "output-testeven-right-TM010-m10.dat",
                    "output-testeven-right-TE011-m10.dat"]
    check_longitudinal_symmetry(z_even_check, prefix="output-testeven-")


    # paper_plot_runs = ["output-plot-TM010-m10.dat",
    #                    "output-plot-TE011-m10.dat"]
    # make_mode_plots(paper_plot_runs, prefix="output-plot-")



