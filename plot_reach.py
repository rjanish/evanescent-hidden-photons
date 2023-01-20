#!/usr/bin/env python

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt


def nuGHz_TE011(Rcm, Lcm):
    xprime01 = spec.jn_zeros(1, 1)
    cm_inv_in_GHz = 29.4118 # 1 = 29.4118 cm GHz, see mathematica
    return 2*np.pi*np.sqrt((xprime01/Rcm)**2 + (np.pi/Lcm)**2)*cm_inv_in_GHz

def plot_reach_from_overlap(results, prefix=""):
    fig, ax = plt.subplots()
    colors = ["r", "b"]
    color_index = 0
    Nsub =10**3
    for filename in results: # expecting two files, TE011 and TM010
        label = filename.split(sep='.')[0][len(prefix):]
        m, overlap = np.loadtxt(filename, skiprows=16).T
        m_subsample = np.geomspace(m[0], m[-1], Nsub)
        reach = overlap**(-0.5)
        color = colors[color_index]
        ax.loglog(m, reach, marker='.', linestyle='',
                  label=label, color=color, alpha=0.7)
        if label == "TE011":
            ax.loglog(m_subsample,
                      two_slopes_valley(m_subsample, reach, -2, 0.5),
                      marker='', linestyle='--', color=color, alpha=0.4,
                      label="slopes: -2 and 0.5")
        if label == "TM010":
            ax.loglog(m_subsample,
                      two_slopes_valley(m_subsample, reach, -1, 1),
                      marker='', linestyle='--', color=color, alpha=0.4,
                      label="slopes: -1 and 1")
        color_index += 1
    ax.legend()
    ax.set_xlabel("m R")
    ax.set_ylabel("epsilon (unormalized)")
    ax.set_title("unormalized reach")
    fig.savefig("checkreach.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    plot_reach_from_overlap(["checkoverlap-TE011.in.out",
                             "checkoverlap-TM010.in.out"],
                             prefix="checkoverlap-")
