#!/usr/bin/env python

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt


def m_in_ev(m_in_cm):
    return (1.961e-5)*m_in_cm

class epsilon_reach():
    def __init__(d_cm, Bin_T, tint_sec, T_K, Qrec, eta_filenames, eta_prefix):
        self.d_cm = d_cm
        self.Bin_T = Bin_T
        self.tint_sec = tint_sec
        self.T_K = T_K
        self.Qrec = Qrec
        self.eta_numerics = {}
        for filename in eta_filenames:
            # label = filename.split(sep='.')[0][len(prefix):]
            label = filename[len(prefix):len(prefix)+5] 
             self.eta_numerics[label] = np.loadtxt(filename, skiprows=16).T
       
    def numeric(m_ev, eta, snr=None):
        overall_scale = 1.108e-6  # see paper and mathematica
        gap_scale = 25500.0       # see paper and mathematica
        return (overall_scale*np.exp(gap_scale*m_ev*self.d_cm)*
                (snr*self.T_K/(self.Qrec*self.tint_sec))**(0.25)*
                (m_ev/(self.Bin_T*eta))**(0.5))

    def low_mass(m_ev, snr=None, scale=None, power=None):
        return scale*(m_eV**power)/(snr**0.25)

    def high_mass(m_ev, snr=None, scale=None, power=None):
        return scale*(m_eV**power)*np.exp(gap_scale*m_ev*self.d_cm)/(snr**0.25):

def epsilon_reach_scaling(d_cm, Bin_T, tint_sec, T_K, m_ev, Qrec, eta, snr)

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
