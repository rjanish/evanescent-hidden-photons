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
    return -0.5*scale*spec.j1(xprime01*r/R)

def plot_profiles_v_mass(results, prefix=""):
    plt.rcParams['text.usetex'] = True
    out = ec.read_effective_current_output(results, prefix)
    key0 = list(out.keys())[0]
    R = out[key0]["radius"]
    L = out[key0]["length"]
    sep = out[key0]["z_start"]
    omega = omegaTE011(R, L)[0]
    colorset = ["indianred", "peru", "forestgreen", "deepskyblue"]
    ordered_masses = []
    colorindex=0
    fig, ax = plt.subplots()
    rs_fine = np.linspace(0, 1, 10**3)
    ax.plot(rs_fine, jhat_TE011(rs_fine, R, L),
            linestyle='dotted', marker='', 
            color='0.0', linewidth=2.5, 
            label=r"analytic result ($\displaystyle m \gg \omega $)")
    for run in out:
        masses = out[run]["m"][:,0,0,0]
        rs = out[run]["r"][0,:,0,0]
        for m_index, m in enumerate(masses):
            ordered_masses.append(m)
            jhat_profile = jhat(out[run]["re_jphi"][m_index, :, 0, 0], m, sep)
            ax.plot(rs, jhat_profile,
                    linestyle='-', linewidth=1.5, marker='',
                    color=colorset[colorindex])
            colorindex += 1

    ax.text(0.2, 0.02, 
            r"$\displaystyle m = {:0.1f} \, \omega$"
            "".format(ordered_masses[0]/omega),
            color=colorset[0], size=12)

    ax.text(0.09, 0.65, 
            r"$\displaystyle m = \omega$"
            "".format(ordered_masses[1]/omega),
            rotation=60, color=colorset[1], size=12)

    ax.text(0.3, 0.715, 
            r"$\displaystyle m = {:0.0f} \, \omega$"
            "".format(ordered_masses[2]/omega),
            rotation=22.5, color=colorset[2], size=12)

    ax.text(0.8, 0.2, 
            r"$\displaystyle m = {:0.0f} \, \omega$"
            "".format(ordered_masses[3]/omega),
            rotation=-40, color=colorset[3], size=12)
    
    ax.legend(frameon=False, handletextpad=0.5)
    ax.set_xlabel(r"$ \displaystyle \rho/R $", fontsize=12)
    ax.set_ylabel(r"$ \displaystyle \hat{\, j_\phi} (\rho) $", 
                  fontsize=12, rotation=0, labelpad=15)
    fig.savefig("profile-comparison.pdf", dpi=300)
    plt.close(fig)


    # fig, ax = plt.subplots()
    # for m_index, m in enumerate(masses):
    #     jhat_profile = jhat(out[run]["re_jphi"][m_index, :, 0, 0], m, sep)
    #     ax.plot(rs, jhat_profile/jhat_TE011(rs, R, L),
    #             linestyle='-', linewidth=1, marker='', alpha=0.6,
    #             label="mR = {:.1e}".format(m))
    #     finite = np.isfinite(jhat_profile/jhat_TE011(rs, R, L))
    #     print(np.mean(jhat_profile[finite]/jhat_TE011(rs[finite], R, L)))

if __name__ == "__main__":
    plot_profiles_v_mass(["profile_v_mass-TE011-1.in.out",
                          "profile_v_mass-TE011-2.in.out"],
                        #  "profile_v_mass-TE011-3.in.out",
                          #"profile_v_mass-TE011-4.in.out",
                          #"profile_v_mass-TE011-5.in.out"],
                          prefix="profile_v_mass-")
