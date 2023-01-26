#!/usr/bin/env python

import warnings
warnings.filterwarnings("error")

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

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

class epsilon_reach():
    def __init__(self, Bin_T, tint_sec, T_K, Qrec, eta_filename, eta_prefix):
        self.overall_scale = 4.907e-9  # see paper and mathematica
        self.Bin_T = Bin_T
        self.tint_sec = tint_sec
        self.T_K = T_K
        self.Qrec = Qrec
        # read numerics output files
        self.runs = ec.read_overlap_output(eta_filename)
        key = list(self.runs.keys())[0]
        self.m_cm_numeric, self.eta_numeric = self.runs[key]["rawdata"].T
        self.d_cm = self.runs[key]["sep"]
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
        

def plot_reach_from_overlap(results, Bin_T, tint_sec, T_K, Qrec, snr, prefix="", name=""):
    fig, ax = plt.subplots()
    color = ['r', 'b']
    mins = []
    for color_index, filename in enumerate(results):
        label = filename.split(".")[0][len(prefix):] 
        reach = epsilon_reach(Bin_T, tint_sec, T_K, Qrec, filename, prefix)
        numeric_reaches = reach.numeric(reach.m_cm_numeric, reach.eta_numeric)*(snr**0.25)
        ax.loglog(reach.m_cm_numeric/eV_cm(), 
                  numeric_reaches, 
                  marker='o', linestyle='-', alpha=0.6, color=color[color_index])
        m_ev_limits = np.array([1e-11, 1e2])
        Nextrap = 3
        m_cm_limits = m_ev_limits*eV_cm()
        m_left_cm = np.logspace(np.log10(m_cm_limits[0]), 
                                np.log10(reach.left_m_cm_numeric[0]), 
                                Nextrap*np.log10(reach.left_m_cm_numeric[0]/m_cm_limits[0]))
        ax.loglog(m_left_cm/eV_cm(), reach.left(m_left_cm)*(snr**0.25), 
                  marker='', linestyle='--', alpha=0.4, color=color[color_index])
        m_right_cm = np.logspace(np.log10(reach.right_m_cm_numeric[1]), 
                                 np.log10(m_cm_limits[1]),
                                 Nextrap*np.log10(m_cm_limits[1]/reach.left_m_cm_numeric[1]))
        ax.loglog(m_right_cm/eV_cm(), reach.right(m_right_cm)*(snr**0.25), 
                  marker='', linestyle='--', alpha=0.4, color=color[color_index])
        ax.axvline(reach.omega/eV_cm(), color=color[color_index], 
                   linestyle='dotted', alpha=0.3)
        mins.append(numeric_reaches.min())
    ax.axvline(2.0/reach.d_cm/eV_cm(), color='k', linestyle='dotted', alpha=0.3)
    ax.set_xlabel("m [eV]")
    ax.set_ylabel("epsilon")
    ax.set_title("reach")
    ax.set_xlim(m_ev_limits)
    ax.set_ylim([0.1*min(mins), 1])
    fig.savefig("{}{}-reach.png".format(prefix, name), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    eta_files = ["reachplot-TE011pancake.in.out",
                 "reachplot-TE011equal.in.out"]
    eta_prefix = "reachplot-"
    Bin_T = 0.05
    tint_sec = 600.0         
    T_K = 1
    Qrec = 1e9
    snr = 5
    plot_reach_from_overlap(eta_files, Bin_T, tint_sec,
                            T_K, Qrec, snr, prefix=eta_prefix, name="weak")
    Bin_T = 0.2
    tint_sec = 3.1e7 # year       
    T_K = 0.1 
    Qrec = 1e12
    snr = 5
    plot_reach_from_overlap(eta_files, Bin_T, tint_sec,
                            T_K, Qrec, snr, prefix=eta_prefix, name="11")
