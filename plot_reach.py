#!/usr/bin/env python

import warnings
warnings.filterwarnings("error")

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import parse_sample_effective_current as ec


def eV_cm():
    return 5.1e4

def fit_decaying_powerlaw(x, y, k):
    """
    Given two x and y values, determine A and p such that
        y = A x^p exp(-k x).
    Returns (A, p).
    """
    sort = np.argpartition(x, 1)
    x, y = x[sort], y[sort]
    p = (np.log(y[1]/y[0]) + k*(x[1] - x[0]))/np.log(x[1]/x[0])
    A = y[0]*(x[0]**-p)*np.exp(k*x[0])
    return A, p


def darkSRF_projection(m_eV):
    omega = 5.4e-6 #eV
    eps_peak = 1e-13
    d_eV = 60*eV_cm()
    k = np.zeros(m_eV.shape)
    evanescent = m_eV > omega
    k[evanescent] = np.sqrt(m_eV[evanescent]**2 - omega**2)
    return eps_peak*(omega/m_eV)*np.exp(k*d_eV)


def alps2_projection(m_eV):
    omega = 1e-4 #eV
    eps_peak = 6e-9
    out = eps_peak*np.ones(m_eV.shape)
    slant = m_eV < omega
    out[slant] = eps_peak*(omega/m_eV[slant])**2
    return out


class epsilon_reach():
    def __init__(self, Bin_T, tint_sec, T_K, Qrec, gap_cm,
                 eta_filename, eta_prefix):
        self.overall_scale = 4.907e-9  # see paper and mathematica
        self.Bin_T = Bin_T
        self.tint_sec = tint_sec
        self.T_K = T_K
        self.Qrec = Qrec
        # read numerics output files
        self.runs = ec.read_overlap_output(eta_filename)
        key = list(self.runs.keys())[0]
        self.m_cm_numeric, self.eta_numeric = self.runs[key]["rawdata"].T
        self.d_cm = gap_cm
        self.omega = self.runs[key]["omega"]
        self.order = np.argpartition(self.m_cm_numeric, self.m_cm_numeric.size - 1)
        # set lower extrapolation (power law)
        self.left_m_cm_numeric  = self.m_cm_numeric[self.order[:2]]
        left_reaches = self.numeric(self.left_m_cm_numeric,
                                    self.eta_numeric[self.order[:2]])
        self.left_scale, self.left_power = fit_decaying_powerlaw(
            self.left_m_cm_numeric, left_reaches, 0.0)
        # set upper extrapolation (power law time exponential)
        self.right_m_cm_numeric = self.m_cm_numeric[self.order[-2:]]
        right_reaches = self.numeric(self.right_m_cm_numeric,
                                     self.eta_numeric[self.order[-2:]])
        self.right_scale, self.right_power = fit_decaying_powerlaw(
            self.right_m_cm_numeric, right_reaches, -0.5*self.d_cm)
        print("\nextrapolating for {}:".format(eta_filename))
        print(" left edge power law is: p={:0.6f}".format(self.left_power))
        print("right edge power law is: p={:0.6f}".format(self.right_power))

    def numeric(self, m_cm, eta):
        return (self.overall_scale*np.exp(0.5*m_cm*self.d_cm)*
                (self.T_K/(self.Qrec*self.tint_sec))**(0.25)*
                (m_cm/(self.Bin_T*eta))**(0.5))

    def left(self, m_cm):
        return self.left_scale*(m_cm**self.left_power)

    def right(self, m_cm, upper_arg_lim=500):
        arg = 0.5*m_cm*self.d_cm
        arg[arg > upper_arg_lim] = upper_arg_lim
        return self.right_scale*(m_cm**self.right_power)*np.exp(arg)


# def plot_reach_from_overlap(filename, setups, prefix="", name=""):

#     fig, ax = plt.subplots()
#     mins = []
#     colors = ["r", "b"]
#     for color_index, setup in enumerate(setups):
#         label = filename.split(".")[0][len(prefix):]
#         reach = epsilon_reach(setup["Bin_T"], setup["tint_sec"],
#                               setup["T_K"], setup["Qrec"], setup["gap_cm"],
#                               filename, prefix)
#         snr = setup["snr"]
#         numeric_reaches = reach.numeric(reach.m_cm_numeric, reach.eta_numeric)*(snr**0.25)
#         ax.loglog(reach.m_cm_numeric/eV_cm(), numeric_reaches,
#                   marker='', linestyle='-', linewidth=1.5,
#                   alpha=0.6, color=color])
#         m_ev_limits = np.array(
#             [10**(np.floor(np.log10(reach.omega/eV_cm())) - 1), 10.0])
#         Nextrap = 50
#         m_cm_limits = m_ev_limits*eV_cm()
#         m_left_cm = np.logspace(np.log10(m_cm_limits[0]),
#                                 np.log10(reach.left_m_cm_numeric[0]),
#                                 Nextrap*int(np.log10(reach.left_m_cm_numeric[0]/m_cm_limits[0])))
#         ax.loglog(m_left_cm/eV_cm(), reach.left(m_left_cm)*(snr**0.25),
#                   marker='', linestyle='-', alpha=0.6, color=color],
#                   linewidth=1.5)
#         m_right_cm = np.logspace(np.log10(reach.right_m_cm_numeric[1]),
#                                  np.log10(m_cm_limits[1]),
#                                  Nextrap*int(np.log10(m_cm_limits[1]/reach.left_m_cm_numeric[1])))
#         ax.loglog(m_right_cm/eV_cm(), reach.right(m_right_cm)*(snr**0.25),
#                   marker='', linestyle='-', alpha=0.6, color=color],
#                   linewidth=1.5)
#         # ax.axvline(reach.omega/eV_cm(), color=color],
#         #            linestyle='dotted', alpha=0.3)
#         mins.append(numeric_reaches.min())
#     # ax.axvline(2.0/reach.d_cm/eV_cm(), color='k', linestyle='dotted', alpha=0.3)
#     ax.set_xlabel("m [eV]")
#     ax.set_ylabel("epsilon")
#     ax.set_title("reach")
#     ax.set_xlim(m_ev_limits)
#     ax.set_ylim([0.1*min(mins), 1])
#     fig.savefig("{}{}-reach.png".format(prefix, name), dpi=160)
#     # plt.show()
#     plt.close(fig)

def ExponentExceptOne(y,pos):
    exp = int(np.rint(np.log10(y)))
    if exp == 0:
        return r'$\displaystyle 1$'.format(exp)
    else:
        return r'$\displaystyle 10^{{{}}}$'.format(exp)

if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True

    setup_current = {"Bin_T" : 0.05,
                     "tint_sec" : 600.0,
                     "T_K" : 3,
                     "Qrec" : 1e9,
                     "snr" : 5,
                     "gap_cm": 0.05,
                     "eta_file":"reachplot-TE011equal.in.out",
                     "file_prefix":"reachplot-"}

    setup_11p = {"Bin_T" : 0.2,
                "tint_sec" : 3.1e7, # year
                "T_K" : 0.1,
                "Qrec" : 1e12,
                "snr" : 5,
                "gap_cm": 0.001,
                "eta_file":"reachplot-TE011pancake50_custom.out",
                "file_prefix":"reachplot-"}

    setup_11e = {"Bin_T" : 0.2,
                "tint_sec" : 3.1e7, # year
                "T_K" : 0.1,
                "Qrec" : 1e12,
                "snr" : 5,
                "gap_cm": 0.05,
                "eta_file":"reachplot-TE011equal.in.out",
                "file_prefix":"reachplot-"}

    setups = [[setup_current, "firebrick", "solid"],
              [setup_11e, 'firebrick', "dashed"],
              [setup_11p, 'firebrick', "dotted"]]

    fig, ax = plt.subplots()

    m_ev_limits = np.array([5*1e-7, 2])
    epsilon_limits = np.array([2e-15, 2e-6])

    for setup, color, ls in setups:
        reach = epsilon_reach(setup["Bin_T"], setup["tint_sec"],
                              setup["T_K"], setup["Qrec"], setup["gap_cm"],
                              setup["eta_file"], setup["file_prefix"])
        snr = setup["snr"]
        numeric_reaches = reach.numeric(reach.m_cm_numeric, reach.eta_numeric)*(snr**0.25)
        ax.loglog(reach.m_cm_numeric/eV_cm(), numeric_reaches,
                  marker='', linestyle=ls, linewidth=1.5,
                  color=color)
        Nleft, Nright = 10**3, 10**3
        margin = 0.05
        m_cm_limits = m_ev_limits*eV_cm()
        m_left_cm = np.logspace(np.log10(m_cm_limits[0]),
                                np.log10(reach.left_m_cm_numeric[0]*(1-margin)),
                                Nleft)
        ax.loglog(m_left_cm/eV_cm(), reach.left(m_left_cm)*(snr**0.25),
                  marker='', linestyle=ls, color=color,
                  linewidth=1.5)
        m_right_cm = np.logspace(np.log10(reach.right_m_cm_numeric[1]*(1+margin)),
                                 np.log10(m_cm_limits[1]), Nright)
        ax.loglog(m_right_cm/eV_cm(), reach.right(m_right_cm)*(snr**0.25),
                  marker='', linestyle=ls, color=color,
                  linewidth=1.5)


    m_dsrf = np.logspace(np.log10(m_ev_limits[0]), np.log10(6.52e-6), 10**3)
    ax.plot(m_dsrf, darkSRF_projection(m_dsrf), marker='',
            linestyle=(0,(.8,1)), color="0.3")

    # m_alps2 = np.logspace(np.log10(1e-5), np.log10(1e-3), 10**3)
    # ax.plot(m_alps2, alps2_projection(m_alps2), marker='',
    #         linestyle=(0,(.8,1)), color="0.3")

    others = [["XENON1T_Solar_SE.txt", "0.6"],
              ["Solar.txt", "0.7"],
              ["COBEFIRAS.txt", "0.5"],
              ["CROWS.txt", "0.65"],
              ["DARKSRF.txt", "0.5"]]
    for limit_filename, color in others:
        try:
            limit = np.loadtxt(limit_filename)
        except:
            limit = np.loadtxt(limit_filename, delimiter=",")
        ax.fill_between(limit[:,0], limit[:,1], epsilon_limits[1],
                        label=limit_filename[:-4], color=color, linewidth=0)

    ax.text(1e-1, 3e-13, "XENON1T (Solar)",
              color='0.2', size=8, rotation=-28)
    ax.text(8e-1, 4e-12, "Sun",
              color='0.3', size=8, rotation=-28)
    ax.text(5.5e-5, 7e-7, r"CMB",
              color='0.8', size=8)
    ax.text(8.5e-6, 7e-7, r"CROWS",
              color='0.3', size=8)
    ax.text(8e-7, 7e-7, r"DarkSRF",
              color='0.8', size=8)
    ax.text(8e-7, 3e-7, r"Pathfinder",
              color='0.8', size=8)
    ax.text(7e-7, 9e-14, r"DarkSRF",
              color='0.3', size=8, rotation=-27)

    ax.text(3e-5, 1.5e-11, r"LSthinW I",
              color="firebrick", size=12, rotation=20)
    ax.text(4e-6, 7e-15, r"LSthinW II",
              color="firebrick", size=12, rotation=20)
    ax.text(1e-3, 2e-14, r"LSthinW III",
              color="firebrick", size=12, rotation=20)

    ax.set_xlabel(r"$\displaystyle m_{A'}\; [{\rm  eV }]$", fontsize=12)
    ax.set_ylabel(r"$\displaystyle \epsilon$", fontsize=16,
                  rotation=0, labelpad=10)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(ExponentExceptOne))
    ax.set_xlim(m_ev_limits)
    ax.set_ylim(epsilon_limits)
    ax.set_yticks([1e-14, 1e-12, 1e-10, 1e-8, 1e-6])
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    fig.savefig("reach.png", dpi=160)
    # plt.show()
    plt.close(fig)
