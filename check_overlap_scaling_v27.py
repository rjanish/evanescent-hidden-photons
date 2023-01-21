#!/usr/bin/env python

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt


def one_slope(x, y, power, index=0):
    return y[index]*(x/x[index])**power


def two_slopes_valley(x, y, s_l, s_r, index_l=0, index_r=-1):
    left = one_slope(x, y, s_l, index=index_l)
    right = one_slope(x, y, s_r, index=index_r)
    return np.max([left, right], axis=0)


def plot_reach_from_overlap(results, prefix=""):
    fig, ax = plt.subplots()
    colors = ["r", "b"]
    color_index = 0
    Nsub =10**3
    for filename in results: # expecting two files, TE011 and TM010
        label = filename[len(prefix):len(prefix)+5]
        m, eta = np.loadtxt(filename, skiprows=16).T
        m_subsample = np.logspace(np.log10(m[0]), np.log10(m[-1]), Nsub)
        reach = (m/eta)**(0.5) 
          # Psig = eps^4 Bemit^2 exp(-2md) \eta^2/m^2 (see paper)
          # so all else fixed, \eps \propto sqrt(m/eta)
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
    fig.savefig("check-reachscaling.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    plot_reach_from_overlap(["checkoverlap-TE011.in.out",
                             "checkoverlap-TM010.in.out"],
                             prefix="checkoverlap-")
