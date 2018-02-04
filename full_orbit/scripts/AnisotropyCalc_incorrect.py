#!/usr/bin/env python

'''
Uses pulsar spectrum and MC data from exit_point_analyzer.nb to compute the anisotropy at Earth.
'''

from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    mcDir = "/Users/acoogan/Dropbox/heliosphere_anisotropy/regular_b_field/mc_data/exit_point_analyzer/"
    energy = 800 # GeV
    alpha = 50 # degrees
    fName = str(energy) + "GeV_alpha" + str(alpha) + "deg.csv"

    # Load exit_point_analyzer MC data
    rawData = np.loadtxt(mcDir + fName, delimiter = ",")
    # Normalize pExitPh to lie in [0, 2pi] rather than [-pi, pi]
    rawData[:,5] = [r if r >= 0 else (r + 2.0 * np.pi) for r in rawData[:,5]]

    print "avg(pExitTh) = " + str(np.mean(rawData[:,4]))
    print "avg(pExitPh) = " + str(np.mean(rawData[:,5])) + "\n"

    # Compute polar anisotropy.  Keep track of particles observed in front that came from back and vice verse.
    frontObsBackExit = 0
    backObsFrontExit = 0
    for dat in rawData:
        if dat[1] < np.pi / 2.0:
            if dat[4] < np.pi / 2.0:
                frontObsBackExit = frontObsBackExit + 1
        else:
            # Check if it came from the front instead
            if dat[4] >= np.pi / 2.0:
                backObsFrontExit = backObsFrontExit + 1

    print "Polar anisotropy: " + str(abs(frontObsBackExit - backObsFrontExit) / (len(rawData) + 0.0) * 100.0)\
            + "%"

    # Compute anisotropy as a function of forward/backward region definition.  Not looking at polar anisotropy
    # yet.
    phBounds = np.linspace(0, np.pi, 150)
    anisos = []
    for phBound in phBounds:
        # Compute anisotropy.  Keep track of particles observed in front that came from back and vice verse.
        frontObsBackExit = 0
        backObsFrontExit = 0
        for dat in rawData:
            if dat[2] <= phBound or dat[2] > phBound + np.pi: # Detected with forward momentum
                # Check if it came from the back instead.  Confusingly, this is the opposite of the check for
                # the momentum's direction.
                if dat[5] <= phBound or dat[5] > phBound + np.pi:
                    frontObsBackExit = frontObsBackExit + 1
            else:
                # Check if it came from the front instead
                if dat[5] > phBound and dat[5] <= phBound + np.pi:
                    backObsFrontExit = backObsFrontExit + 1

        anisos.append(abs(frontObsBackExit - backObsFrontExit) / (len(rawData) + 0.0))

    # Plot anisotropy as a function of azimuthal angle
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(phBounds + np.pi / 2.0, np.multiply(anisos, 100.0), "-")
    ax.set_xlabel("Detector pointing angle (rad)")
    ax.set_ylabel("Dipole anisotropy (%)")
    ax.set_title(r"Anisotropy given isotropic LIS ($\alpha = " + str(alpha) + "^\circ$)")
    ax.set_xlim(min(phBounds + np.pi / 2.0), max(phBounds + np.pi / 2.0))
    plt.show()


