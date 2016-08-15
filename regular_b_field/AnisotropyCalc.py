#!/usr/bin/env python

'''
Uses pulsar spectrum and MC data from exit_point_analyzer.nb to compute the anisotropy at Earth.
'''

from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt

# Converts r from au to parsecs
def toPc(r):
    return r * 4.848e-6

# Converts r from parsecs to au
def toAu(r):
    return r / 4.848e-6

# Galactic cosmic ray propagation parameters from model 0 in arXiv:0905.0636.
b0 = 1.4e-16 # GeV^-1 s^-1
d0 = 3781.0e-12 # pc^2 s^-1
e0 = 1.0 # GeV
delta = 33.0e-2

# Pulsar + other  parameters from arXiv:0905.0636
ecut = 1100.0 # GeV
q0 = 1.0e13 # Arbitrary source normalization
deltaT = 6.0e4 * 365.0 * 24.0 * 3600.0 # s.  Time delay between pulsar birth and positron emission
Gamma = 1.5
c = 9716.0e-12 # pc / s.  Speed of light
emax = 1.0 / (b0 * deltaT)

# Diffusion coefficient
def d(e):
    return d0 * pow(e / e0, delta)

# Diffusion length scale
def rdiff(e):
    return 2.0 * np.sqrt(d(e) * deltaT * (1.0 - pow(1.0 - e / emax, 1.0 - delta)) / ((1.0 - delta) * e / emax))

# Returns electron or positron number density per unit energy (ie, dN/dE = d^4n/(dE dx dy dz)) from a pulsar.
# The one here is approximately monogem.  See arXiv:0905.0636.
# Arguments
#   e: measured positron energy (GeV)
#   r, th, ph: position relative to the sun at which to compute dNdE (r measured in au)
#   rP, thP, phP: pulsar position relative to the sun (r measured in parsecs)
def dNdE(e, r, th, ph, rP = 290, thP = np.pi / 2.0, phP = 0.0):
    # Get distance between pulsar and (r, th, ph)
    dist = np.sqrt(pow(toPc(r), 2) + rP*rP - 2 * toPc(r) * rP * (np.sin(th) * np.sin(thP) * np.cos(ph - phP) 
        + np.cos(th) * np.cos(thP)))

    return q0 / (pow(np.pi, 3.0/2.0) * pow(rdiff(e), 3.0)) * pow(1.0 - e / emax, Gamma - 2.0) \
            * pow(e, -Gamma) * np.exp(-e / ((1.0 - e / emax) * ecut)) * (1.0 if emax > e else 0.0) \
            * np.exp(-pow(r / rdiff(e), 2))

# Computes dipole anisotropy as a function of azimuthal angle for a given energy (GeV) and tilt angle (deg).
# Returns:
#   detector pointing angles and corresponding anisotropies in %
def getAnisos(energy, alpha):
    mcDir = "/Users/acoogan/Dropbox/heliosphere_anisotropy/regular_b_field/mc_data/exit_point_analyzer/"
    fName = str(energy) + "GeV_alpha" + str(alpha) + "deg.csv"

    # Load exit_point_analyzer MC data
    rawData = np.loadtxt(mcDir + fName, delimiter = ",")
    # Normalize pExitPh to lie in [0, 2pi] rather than [-pi, pi]
    rawData[:,5] = [r if r >= 0 else (r + 2.0 * np.pi) for r in rawData[:,5]]

    # Struct to hold mc runs.  Not needed right now.
    #Trajectory = namedtuple('Trajectory', ["pObs", "pObsTh", "pObsPh", "pExit", "pExitTh", "pExitPh"])

    '''
    # Compute polar anisotropy.  Keep track of particles observed in front that came from back and vice verse.
    top = 0.0
    bottom = 0.0
    for dat in rawData:
        if dat[1] < np.pi / 2.0:
            top = top + 1.0
        else:
            bottom = bottom + 1.0

    print "Polar anisotropy: " + str(abs(top - bottom) / (top + bottom) * 100.0) + "%"
    '''

    # Compute anisotropy as a function of forward/backward region definition.  Not looking at polar anisotropy
    # yet. TODO: figure out how to go from -pi/2 to pi/2 instead...
    phBounds = np.linspace(0.0, np.pi, 1000)
    anisos = []
    anisosNoPsr = []
    for phBound in phBounds:
        # Compute anisotropy with anisotropic LIS
        forward = 0.0
        backward = 0.0
        # Compute anisotropy with isotropic LIS
        forwardNoPsr = 0.0
        backwardNoPsr = 0.0

        for dat in rawData:
            if dat[2] <= phBound or dat[2] > phBound + np.pi: # Detected with forward momentum
                # These are obtained from taking dNdE[800, 290] * (1 + 0.5 * x * aniso[800, 290]) in
                # the pulsar_spectrum notebook
                forward = forward + 0.185571 * (4.0 + 0.00284743 * np.sin(dat[4]) * np.sin(dat[5]))
                forwardNoPsr = forwardNoPsr + 1.0
            else:
                backward = backward + 0.185571 * (4.0 + 0.00284743 * np.sin(dat[4]) * np.sin(dat[5]))
                backwardNoPsr = backwardNoPsr + 1.0

        anisos.append(abs(forward - backward) / (forward + backward))
        anisosNoPsr.append(abs(forwardNoPsr - backwardNoPsr) / (forwardNoPsr + backwardNoPsr))

    return phBounds, np.multiply(anisos, 100.0), np.multiply(anisosNoPsr, 100.0)

#############################################################################################################

if __name__ == "__main__":
    # Plot anisotropy as a function of azimuthal angle
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    energy = 800
    angle = 50 # deg

    # Get anisotropies at 800 GeV
    phBs, anisos, anisosNoPsr = getAnisos(energy, 50)
    phBs = np.multiply(phBs, 180.0 / np.pi)
    ax.plot(phBs, anisos - anisosNoPsr, "b-", linewidth = 2, label = "Helio")

    # Get x bounds
    minPhB = min(phBs)
    maxPhB = max(phBs)

    # Put in anisotropy from pulsars at heliosphere for reference
    anisoMonogem = 0.52
    ax.plot(phBs, [anisoMonogem * np.cos(phB * np.pi / 180.0) for phB in phBs], '--r', label = "No helio",
            linewidth = 2)
    #anisoGeminga = 0.07

    # Make the plot look nice
    ax.set_xlabel("Detector pointing relative to pulsar (deg)", fontsize = 16)
    ax.set_xlim(minPhB, maxPhB)
    ax.set_ylabel(r"$\Delta_{\mathrm{Pulsar}} - \Delta_{\mathrm{No\ pulsar}}$ (%)", fontsize = 16)
    #ax.set_yscale("log")

    ax.legend(loc = "upper right")

    ax.set_title(r"Heliosphere's effect on $\Delta$ ($p_{\mathrm{obs}} = " + str(energy)
            + r"$ GeV, $\alpha=" + str(angle) + r"^\circ$)", fontsize = 18)
    plt.show()



