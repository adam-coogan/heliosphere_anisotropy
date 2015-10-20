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
    es = []
    pairs = []
    minEkLIS = 1e10 # Some huge number
    maxEkLIS = 0

    # Read pairs from files
    for df in dataFiles:
        print('Loading ' + basename(df))
        # TODO: current convention: energy is part of file name.  Fine for isotropic LIS...
        e = float(splitext(basename(df))[0][:-3])
        es.append(e)
        # Load the data
        rawData = np.loadtxt(fullDir + df, delimiter = ',')
        # Make the pairs
        [pairs.append((e, eLIS)) for eLIS in rawData[:, 3]]

    # Sort observed energies
    es.sort()

    return pairs, es

def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def listEks(pairs, es, minBin, maxBin, nBins):
    """
    Puts Ek (the energy at Earth) into a list of Eks for the bin corresponding to EkLIS.
    Arguments:
        pairs: pairs of (Ek, EkLIS).
        es: simulated energies at earth
        minBin, maxBin: minimum and maximum bin boundaries.
        nBins: number of bins for EkLISs.
    Returns:
        List of bin boundaries and measured energies at Earth (Ek) for particles that started at the
        heliopause with energies (EkLIS) in the bin.
    """
    # Set up the bins
    eLISBins = pairwise([10.0 ** a for a in np.linspace(np.log10(minBin), np.log10(maxBin), nBins)])
    listedEks = [(b, np.zeros(len(es))) for b in eLISBins]

    # Loop over all pairs
    for pair in pairs:
        e = pair[0]
        eLIS = pair[1]

        # Make sure EkLIS falls into one of the bins
        if eLIS < maxBin and eLIS >= minBin:
            # Loop over bin lists
            for le in listedEks:
                # Check whether EkLIS falls into the current bin
                if eLIS >= le[0][0] and eLIS < le[0][1]:
                    # Add one to the count for the bin corresponding to e
                    le[1][es.index(e)] = le[1][es.index(e)] + 1
                    break

    return listedEks

def roundedStr(num, n):
    """
    Returns a string containing num rounded to n significant figures.
    """
    return '{number:.{digits}g}'.format(number = num, digits = n)

#############################################################################################################

# Bin boundaries for EkLIS
minEkLIS = 0.5
maxEkLIS = 5.0
neLISBins = 10

# Get Ek counts for each EkLIS bin
ePairs, es = loadPairs('hyades/me/alt0/')
eBinLists = listEks(ePairs, es, minEkLIS, maxEkLIS, neLISBins)

# Plot Ek measurement probabilities for each bin
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for ebl in eBinLists:
    # Normalize Ek counts to get probabilities P(E | E^LIS \in bin)
    ax.plot(es, ebl[1] / sum(ebl[1]), '.-', label = '(' + roundedStr(ebl[0][0], 3) + ', ' \
            + roundedStr(ebl[0][1], 3) + ')')
    print(str(sum(ebl[1])) + ' partices in bin (' + roundedStr(ebl[0][0], 3) + ', ' \
            + roundedStr(ebl[0][1], 3) + ')')

print('Plotting P(E | E^LIS \in bin)')

# Format the plot
ax.set_xlim(min(es), max(es))
ax.set_ylim(0.0, 1.0)
ax.set_xscale('log')
ax.set_xlabel(r'$E_k$ at Earth (GeV)')
ax.set_ylabel(r'$P(E | E^{\mathrm{LIS}}$')
ax.set_title(r'Green\'s function for $e^-$ heliospheric propagation')
plt.legend(loc='upper left')

plt.show()


