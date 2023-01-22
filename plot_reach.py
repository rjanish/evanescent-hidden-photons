#!/usr/bin/env python

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt


def m_in_cm(m_in_eV):
    return m_in_eV/(1.961e-5)

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
    def __init__(d_cm, Bin_T, tint_sec, T_K, Qrec, eta_filename, eta_prefix):
        self.overall_scale = 1.108e-6  # see paper and mathematica
        self.gap_scale = 25500.0       # see paper and mathematica
        self.d_cm = d_cm
        self.Bin_T = Bin_T
        self.tint_sec = tint_sec
        self.T_K = T_K
        self.Qrec = Qrec
        self.m_cm_numeric, self.eta_numeric = np.loadtxt(eta_filename, skiprows=16).T
        ordered_indicies = np.argpartition(self.m_cm_numeric, self.m_cm_numeric.size - 1)
        # set lower extrapolation (power law)
        left_masses  = self.m_cm_numeric(ordered_indicies[:2])
        left_reaches = self.numeric_reach(left_masses, 
                                          self.eta_numeric(ordered_indicies[:2]))
        self.left_exponent, self.left_scale = fit_decaying_powerlaw(
            left_masses, left_reaches, 0.0)
        # set upper extrapolation (power law time exponential)
        right_masses  = self.m_cm_numeric(ordered_indicies[-2:])
        right_reaches = self.numeric_reach(right_masses, 
                                           self.eta_numeric(ordered_indicies[-2:]))
        self.right_exponent, self.right_scale = fit_decaying_powerlaw(
            right_masses, right_reaches, -2*self.gap_scale*self.d_cm)

    def numeric_reach(m_ev, eta):
        return (self.overall_scale*np.exp(self.gap_scale*m_ev*self.d_cm)*
                (self.T_K/(self.Qrec*self.tint_sec))**(0.25)*
                (m_ev/(self.Bin_T*eta))**(0.5))

    def left_reach(m_ev):
        return self.left_scale*(m_eV**self.left_power)/(snr**0.25)

    def right_reach(m_ev):
        return (self.right_scale*(m_eV**self.right_power) * 
                np.exp(self.gap_scale*m_ev*self.d_cm))

    def evaulate(m_ev, snr=1.0):
        # compute reach by interpolation within numerics and using left/right outisde it
        return

def plot_reach_from_overlap(results, d_cm, Bin_T, tint_sec,
                            T_K, Qrec, snr, prefix=""):
    fig, ax = plt.subplots()
    for filename in results: # expecting two files, TE011 and TM010
        # label = filename.split(sep='.')[0][len(prefix):]
        label = filename[len(prefix):len(prefix)+5] 
        m_cm, eta = np.loadtxt(filename, skiprows=16).T
        reach = epsilon_reach_numeric(d_cm, Bin_T, tint_sec, T_K, m_in_ev(m_cm), Qrec, eta, snr)
        ax.loglog(m_in_ev(m_cm), reach, marker='.', linestyle='',
                  label=label, alpha=0.7)
    ax.legend()
    ax.set_xlabel("m [eV]")
    ax.set_ylabel("epsilon")
    ax.set_title("reach")
    fig.savefig("reach.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    eta_files = ["checkoverlap-TE011.in.out",
                 "checkoverlap-TM010.in.out"]
    eta_prefix = "checkoverlap-"
    d_cm = 0.05 
    Bin_T = 0.05
    tint_sec = 600.0         
    T_K = 3 
    Qrec = 1e10
    snr = 5
    plot_reach_from_overlap(eta_files, d_cm, Bin_T, tint_sec,
                            T_K, Qrec, snr, prefix=eta_prefix)
