#!/usr/bin/env python

"""
Modulates LIS using the force field approximation.
"""

from AnisotropyCalc import jLISLangner, getJ, rundataPath
import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt
import os.path
from os import listdir

def getJFF(ek, jLIS, phiF):
    """
    Computes j(ek) at Earth given a local interstellar spectrum jLIS.  ek: GeV, jLIS: (m^2 GeV s sr)^-1,
    phiF: MV.
    """

    # Fisk potential (GV)
    phiF = phiF / 1000.0
    Z = -1
    A = 1
    # Energy at HP, according to the force field approximation.
    ekIS = ek + phiF * abs(Z) / A
    me = 5.11e-4 # Electron mass in GeV

    return (2 * me * A * ek  + A**2 * ek**2) / (2 * me * A * ekIS + A**2 * ekIS**2) * jLIS(ekIS)

def getJNewFF(ek, jLIS):
    """
    Computes j(ek) at Earth given a local interstellar spectrum jLIS.  ek: GeV, jLIS: (m^2 GeV s sr)^-1.
    Uses the potential from arXiv:1511.01507, by Cholis, Hooper and Linden.
    """

    # Model parameters
    phi0 = 0.35 # GV
    Be = 5 # Magnetic field at earth (nT)
    phi1 = 0.977 # GV
    R0 = 0.5 # Reference rigidity (GV)
    alpha = 0 # Tilt angle (rad)
    me = 5.11e-4 # Electron mass in GeV
    R = sqrt(ek**2 + 2 * ek * me) # Rigidity (GV)
    beta = R / (ek + me) # v / c
    q = -1 # Charge
    A = -1 # HMF polarity

    # Modulation potential
    Phi = phi0 * (Be / 4) + phi1 * (1 if -q*A > 0 else 0) * (Be / 4) * (1 + (R/R0)**2) / (beta * (R/R0)**3) \
            * (alpha / (pi / 2))**4

    # Find interstellar rigidity and use to compute interstellar kinetic energy
    Rism = R + Phi
    ekIS = sqrt(me**2 + Rism**2) - me

    return  R**2 / Rism**2 * jLIS(ekIS)

if __name__ == "__main__":
    print('Computing modulated spectrum (SDE)')
    meE0s, meJs, meSigmas = map(list, zip(*[getJ('me/alt0/' + os.path.splitext(f)[0]) \
            for f in listdir(rundataPath + 'me/alt0/')]))
    meE0s = [x0['e0'] for x0 in meE0s]

    # Plot spectra
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # My modulated spectrum
    ax.errorbar(meE0s, meJs, yerr = meSigmas, fmt = '.', label='SDE')

    # Unmodulated spectrum
    eksLangner = 10.0 ** np.linspace(-3.0, 2.0, 1000) # Energies in GeV
    ax.plot(eksLangner, [jLISLangner(ek) for ek in eksLangner], label = r'$j^\mathrm{LIS}$')

    print('Computing modulated spectrum (force field)')
    ax.plot(eksLangner, [getJFF(ek, jLISLangner, 200) for ek in eksLangner], '--', \
            label = '$\phi_F = 200$ MV')
    ax.plot(eksLangner, [getJFF(ek, jLISLangner, 300) for ek in eksLangner], '--', \
            label = '$\phi_F = 300$ MV')

    print('Computing modulated spectrum (force field)')
    ax.plot(eksLangner, [getJNewFF(ek, jLISLangner) for ek in eksLangner], '--', label = r'New method')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-3, 1e1)
    ax.set_ylim(1e-4, 1e4)
    ax.set_xlabel(r'T (GeV)')
    ax.set_ylabel(r'Intensity ($\mathrm{MeV}^{-1}\ \mathrm{s}^{-1}\ \mathrm{sr}^{-2}\ \mathrm{m}^{-2}$)')
    ax.set_title(r'Cosmic ray spectrum modulation ($A_c < 0$, $\alpha = 0$)')
    plt.legend(loc='upper right')

    plt.show()




