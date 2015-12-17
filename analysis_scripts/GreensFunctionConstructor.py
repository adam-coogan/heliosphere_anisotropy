#!/usr/bin/env python

"""
Modulates energy spectrum using results from propagation code.
"""

import numpy as np
import scipy
from matplotlib import rc
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import subprocess
import os.path
from os import listdir
from os.path import isfile, join
from AnisotropyCalc import getJ, rundataPath
from itertools import tee, izip

#############################################################

def jLISDelta(ek, ekLISBin):
    """
    Delta function LIS: returns 1 if ek is in ekLISBin.
    """
    if ek >= ekLISBin[0] and ek < ekLISBin[1]:
        return 1
    else:
        return 0

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

#############################################################

if __name__ == "__main__":
    # E^LIS bins
    #ekLISBins = list(pairwise([a * 1e0 for a in range(1, 3)]))
    ekLISBins = [(0.05, 0.06), (0.1, 0.2), (0.5, 0.6), (1.0, 2.0)]
    # Store e0s, js and sigmas for each bin
    binE0s = []
    binJs = []
    binSigmas = []

    for ekb in ekLISBins:
        # Calculate spectrum at Earth for the chosen bin
        print('Computing spectrum at Earth for E^LIS \in ' + str(ekb))
        e0s, js, sigmas = map(list, zip(*[getJ('me/alt0/' + os.path.splitext(f)[0], \
                jLIS = lambda ek: jLISDelta(ek, ekb)) for f in listdir(rundataPath + 'me/alt0/')]))
        e0s = [x0['e0'] for x0 in e0s]

        binE0s.append(e0s)
        binJs.append(js)
        binSigmas.append(sigmas)

    # Plot spectrum
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # Plot each bins spectrum at Earth
    for e0s, js, sigmas, ekb in zip(binE0s, binJs, binSigmas, ekLISBins):
        ax.errorbar(e0s, js, yerr = sigmas, fmt = '.', label = '(' + str(ekb[0]) + ', ' + str(ekb[1]) + ')')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-3, 1e1)
    ax.set_ylim(1e-6, 1e0)
    ax.set_xlabel(r'T (GeV)')
    ax.set_ylabel(r'Intensity ($\mathrm{MeV}^{-1}\ \mathrm{s}^{-1}\ \mathrm{sr}^{-2}\ \mathrm{m}^{-2}$)')
    ax.set_title(r'Modulation of injection at LIS ($A_c < 0$, $\alpha = 0$)')
    plt.legend(loc='lower right')

    plt.show()


