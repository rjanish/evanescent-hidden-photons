#!/usr/bin/env python

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.special as spec

import parse_sample_effective_current as ec


def omegaTE011(R, L):
    xprime01 = spec.jn_zeros(1, 1)
    return np.sqrt((xprime01/R)**2 + (np.pi/L)**2)

def jhat(j, m, sep):
    return j/(m*np.exp(-m*sep))

def jhat_TE011(r, R, L):
    xprime01 = spec.jn_zeros(1, 1)
    scale = np.sqrt(2.0)*np.pi/(omegaTE011(R, L)*L*spec.j0(xprime01))
    return 0.5*scale*spec.j1(xprime01*r/R)

def plot_profiles_v_mass(results, prefix=""):
    out = ec.read_effective_current_output(results, prefix)
    key = list(out.keys())[0]
    R = out[key]["radius"]
    L = out[key]["length"]
    sep = out[key]["z_start"]
    masses = out[key]["m"][:,0,0,0]
    omega = omegaTE011(R, L)[0]
    rs = out[key]["r"][0,:,0,0]
    fig, ax = plt.subplots()
    rs_fine = np.linspace(0, 1, 10**3)
    ax.plot(rs_fine, jhat_TE011(rs_fine, R, L),
            linestyle='-', marker='', alpha=0.8,
            color='k', linewidth=3, label="prediction for m>>w")
    for m_index, m in enumerate(masses):
        jhat_profile = jhat(out[key]["re_jphi"][m_index, :, 0, 0], m, sep)
        ax.plot(rs, jhat_profile,
                linestyle='-', linewidth=1, marker='', alpha=0.6,
                label="m/w = {:0.1f}".format(m/omega))
    ax.legend()
    ax.set_xlabel("r")
    ax.set_ylabel("j_hat")
    fig.savefig("profile-comparison.png", dpi=160)
    plt.close(fig)


    # fig, ax = plt.subplots()
    # for m_index, m in enumerate(masses):
    #     jhat_profile = jhat(out[key]["re_jphi"][m_index, :, 0, 0], m, sep)
    #     ax.plot(rs, jhat_profile/jhat_TE011(rs, R, L),
    #             linestyle='-', linewidth=1, marker='', alpha=0.6,
    #             label="mR = {:.1e}".format(m))
    #     finite = np.isfinite(jhat_profile/jhat_TE011(rs, R, L))
    #     print(np.mean(jhat_profile[finite]/jhat_TE011(rs[finite], R, L)))

if __name__ == "__main__":
    plot_profiles_v_mass(["profile_v_mass-TE011.in.out"],
                             prefix="profile_v_mass-")
