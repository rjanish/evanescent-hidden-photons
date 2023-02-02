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
    def __init__(self, setup):
        self.overall_scale = 4.907e-9  # see paper and mathematica
        self.Bin_T = setup["Bin_T"]
        self.tint_sec = setup["tint_sec"]
        self.T_K = setup["T_K"]
        self.Qrec = setup["Qrec"]
        self.readout = setup["readout"]
        # read numerics output files
        self.runs = ec.read_overlap_output(setup["eta_filename"])
        key = list(self.runs.keys())[0]
        self.m_cm_numeric, self.eta_numeric = self.runs[key]["rawdata"].T
        self.d_cm = setup["gap_cm"]
        self.omega_cm = self.runs[key]["omega"]
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
        print("\nextrapolating for {}:".format(setup["eta_filename"]))
        print(" left edge power law is: p={:0.6f}".format(self.left_power))
        print("right edge power law is: p={:0.6f}".format(self.right_power))

    def numeric(self, m_cm, eta):
        if self.readout == "phase-sensitive":
            scale = 4.907e-9  # see paper and mathematica
            return (scale*np.exp(0.5*m_cm*self.d_cm)*
                    (self.T_K/(self.Qrec*self.tint_sec))**(0.25)*
                    (m_cm/(self.Bin_T*eta))**(0.5))
        if self.readout == "dicke":
            scale = 7.936e-8  # see paper and mathematica
            return (scale*np.exp(0.5*m_cm*self.d_cm)*
                    self.T_K**(0.25)*
                    (m_cm/(self.Bin_T*eta))**(0.5)*
                    (self.omega_cm/(self.tint_sec*self.Qrec**3))**(0.125))

    def left(self, m_cm):
        return self.left_scale*(m_cm**self.left_power)

    def right(self, m_cm, upper_arg_lim=500):
        arg = 0.5*m_cm*self.d_cm
        arg[arg > upper_arg_lim] = upper_arg_lim
        return self.right_scale*(m_cm**self.right_power)*np.exp(arg)


def ExponentExceptOne(y, pos):
    exp = int(np.rint(np.log10(y)))
    if exp == 0:
        return r'$\displaystyle 1$'.format(exp)
    else:
        return r'$\displaystyle 10^{{{}}}$'.format(exp)


def EvenExponentOnly(y, pos):
    exp = int(np.rint(np.log10(y)))
    if exp % 2 == 1:
        return ""
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
                     "eta_filename":"reachplot-TE011equal.in.out",
                     "readout":"dicke"}

    setup_11p = {"Bin_T" : 0.2,
                "tint_sec" : 3.1e7, # year
                "T_K" : 0.1,
                "Qrec" : 1e12,
                "snr" : 5,
                "gap_cm": 0.001,
                "eta_filename":"reachplot-TE011pancake50_custom.out",
                "readout":"phase-sensitive"}

    setup_11e = {"Bin_T" : 0.2,
                "tint_sec" : 3.1e7, # year
                "T_K" : 0.1,
                "Qrec" : 1e12,
                "snr" : 5,
                "gap_cm": 0.05,
                "eta_filename":"reachplot-TE011equal.in.out",
                "readout":"phase-sensitive"}

    setups = [[setup_current, "firebrick", "solid"],
              [setup_11e, 'firebrick', "dashed"],
              [setup_11p, 'firebrick', "dotted"]]

    fig, ax = plt.subplots()

    m_ev_limits = np.array([1e-6, 1])
    epsilon_limits = np.array([1e-15, 1e-7])

    for setup, color, ls in setups:
        reach = epsilon_reach(setup)
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

    ax.text(2e-1, 1.8e-12, "XENON1T",
              color='0.3', size=8, rotation=-27)

    ax.text(1.7e-1, 3e-11, "Solar Cooling",
              color='0.3', size=8, rotation=-27)

    ax.text(2.5e-5, 5e-8, r"CMB",
              color='0.8', size=8)

    ax.text(6.25e-6, 5e-8, r"CROWS",
              color='0.3', size=7)

    ax.text(1.2e-6, 5e-8, r"DarkSRF",
              color='0.8', size=8)

    ax.text(1.2e-6, 2.8e-8, r"Pathfinder",
              color='0.8', size=8)

    ax.text(1.1e-6, 1.8e-13, r"DarkSRF",
              color='0.3', size=8, rotation=-27)

    ax.text(2.5e-5, 1e-10, r"LSthinW {\rm I}",
              color="firebrick", size=10, rotation=17)

    ax.text(5e-6, 2.5e-14, r"LSthinW {\rm II}",
              color="firebrick", size=10, rotation=17)

    ax.text(1e-3, 7e-14, r"LSthinW {\rm III}",
              color="firebrick", size=10, rotation=17)


    ax.set_xlabel(r"$\displaystyle m_{A'}\; [{\rm \tiny eV }]$", fontsize=14)

    # ax.text(3e-4, 1.3e-16, r"$\displaystyle m_{A'}$", size=18)
    # ax.text(1.15e-3, 1.65e-16, r"$\displaystyle [{ \rm eV}]$", size=10)

    ax.text(2e-7, 8e-12, r"$\displaystyle \epsilon$", fontsize=20, rotation=0)
    ax.set_xlim(m_ev_limits)
    ax.set_ylim(epsilon_limits)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(ExponentExceptOne))
    ax.set_yticks(10.0**np.arange(-15,-6))
    minor_ticks = np.array(
        [np.arange(2, 10)*(10**n) for n in range(-15, -7)]).flatten()
    ax.yaxis.set_minor_locator(ticker.FixedLocator(minor_ticks))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(EvenExponentOnly))
    ax.xaxis.set_tick_params(labelsize=12)
    ax.yaxis.set_tick_params(labelsize=12)
    ax.set_aspect(0.5)
    fig.savefig("reach.pdf", dpi=300)
    # plt.show()
    plt.close(fig)
