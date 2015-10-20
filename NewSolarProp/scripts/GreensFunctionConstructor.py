#!/usr/bin/env python

"""
Use this to construct Greens functions.  G(Ek;EkLIS) gives the probability of observing a particle with
kinetic energy Ek at earth given that it had energy EkLIS at the heliopause.  Here we assume that the LIS is
isotropic.
"""

import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import basename, isfile, join, splitext
from itertools import tee, izip

runDir = 'rundata/'

def loadPairs(dataDir):
    """
    Load all the data in the given directory.  We will parse this data into (Ek, EkLIS) pairs.
        dataDir: directory containing the data files.
        Returns: (Ek, EkLIS) pairs for each run in each data file.
    """
    fullDir = runDir + '/' + dataDir + '/'
    dataFiles = listdir(fullDir)
    eks = []
    pairs = []
    minEkLIS = 1e10 # Some huge number
    maxEkLIS = 0

    # Read pairs from files
    for df in dataFiles:
        # TODO: current convention: energy is part of file name.  Fine for isotropic LIS...
        ek = float(splitext(basename(df))[0][:-3])
        eks.append(ek)
        # Load the data
        rawData = np.loadtxt(fullDir + df, delimiter = ',')
        # Make the pairs
        [pairs.append((ek, ekLIS)) for ekLIS in rawData[:, 3]]

    # Sort observed energies
    eks.sort()

    return pairs, eks

def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def listEks(pairs, eks, minBin, maxBin, nBins):
    """
    Puts Ek (the energy at Earth) into a list of Eks for the bin corresponding to EkLIS.
    Arguments:
        pairs: pairs of (Ek, EkLIS).
        eks: simulated energies at earth
        minBin, maxBin: minimum and maximum bin boundaries.
        nBins: number of bins for EkLISs.
    Returns:
        List of bin boundaries and measured energies at Earth (Ek) for particles that started at the
        heliopause with energies (EkLIS) in the bin.
    """
    # Set up the bins
    ekLISBins = pairwise([10.0 ** a for a in np.linspace(np.log10(minBin), np.log10(maxBin), nBins)])
    listedEks = [(b, np.zeros(len(eks))) for b in ekLISBins]

    # Loop over all pairs
    for pair in pairs:
        ek = pair[0]
        ekLIS = pair[1]

        # Make sure EkLIS falls into one of the bins
        if ekLIS < maxBin and ekLIS >= minBin:
            # Loop over bin lists
            for lek in listedEks:
                # Check whether EkLIS falls into the current bin
                if ekLIS >= lek[0][0] and ekLIS < lek[0][1]:
                    # Add one to the count for the bin corresponding to ek
                    lek[1][eks.index(ek)] = lek[1][eks.index(ek)] + 1
                    break

    return listedEks

def roundedStr(num, n):
    """
    Returns a string containing num rounded to n significant figures.
    """
    return '{number:.{digits}g}'.format(number = num, digits = n)

#############################################################################################################

# Bin boundaries for EkLIS
minEkLIS = 0.05
maxEkLIS = 0.5
nekLISBins = 10

# Get Ek counts for each EkLIS bin
ePairs, eks = loadPairs('hyades/me/alt0/')
ekBinLists = listEks(ePairs, eks, minEkLIS, maxEkLIS, nekLISBins)

# Plot Ek measurement probabilities for each bin
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for ekbl in ekBinLists:
    # Normalize Ek counts to get probabilities P(E | E^LIS \in bin)
    ax.plot(eks, ekbl[1] / sum(ekbl[1]), '.-', label = '(' + roundedStr(ekbl[0][0], 3) + ', ' \
            + roundedStr(ekbl[0][1], 3) + ')')

# Format the plot
ax.set_xlim(min(eks), max(eks))
ax.set_ylim(0.0, 1.0)
ax.set_xscale('log')
ax.set_xlabel(r'$E_k$ at Earth (GeV)')
ax.set_ylabel(r'$P(E_k; E_k^{\mathrm{LIS}}$')
ax.set_title(r'Green\'s function for $e^-$ heliospheric propagation')
plt.legend(loc='upper left')

plt.show()


