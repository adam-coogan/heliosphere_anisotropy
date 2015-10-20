#!/usr/bin/env python

"""
Use this to construct Greens functions.  G(Ek;EkLIS) gives the probability of observing a particle with
kinetic energy Ek at earth given that it had energy EkLIS at the heliopause.  Here we assume that the LIS is
isotropic.
"""

from itertools import tee, izip
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter, truediv
from os import listdir
from os.path import basename, isfile, join, splitext

runDir = 'rundata/'

def loadEnergies(dataDir):
    """
    Load all the data in the given directory.
        dataDir: directory containing the data files.
        Returns: [Ek], [nRuns], [[EkLISs]], where Ek is the energy at earth, nRuns is the number of runs
        generated at this energy and [EkLISs] is a list of the exit energies for these runs.
    """

    fullDir = runDir + '/' + dataDir + '/'
    dataFiles = listdir(fullDir)
    # We use a list of tuples to make sorting by ek easy
    eTuples = []

    # Read energies from files
    for df in dataFiles:
        # Current convention: energy is part of file name.  Fine for isotropic LIS...
        ek = float(splitext(basename(df))[0][:-3])
        # Load the data
        rawData = np.loadtxt(fullDir + df, delimiter = ',')
        # Make the tuple
        eTuples.append((ek, len(rawData[:, 3]), rawData[:, 3]))

    # Sort observed energies
    eTuples.sort(key = itemgetter(0))

    return eTuples[:, 0], eTuples[:, 1], eTuples[:, 2]

def getProbabilities(dataDir, ekLISBins):
    """
    Uses Bayes' Theorem to compute P(E | E^LIS \in E^LIS_bin).
    Arguments:
        dataDir: directory containing the data files.
        ekLISBins: list of tuples specifying E^LIS_bin.
    Returns:
        [([E], [P(E | E^LIS \in E^LIS_bin)])]: probability of observing particle with energy E at Earth for each
        bin.  Ordering corresponds to bin ordering.
    """

    # Load energy data
    eks, nRuns, ekLISsList = loadEnergies('hyades/me/alt0/')
    probs = [(eks, [0] * len(eks))] * len(ekLISBins)
    # Count the number of exit energies in each bin
    binCounts = [0] * len(ekLISBins)

    # Loop over energies
    for ek, n, ekLISs in zip(eks, nRuns, ekLISsList):
        ekIndex = eks.index(ek)

        # Loop over exit energies for the run
        for ekLIS in ekLISs:
            # Put ekLIS in the correct bin
            for b in ekLISBins:
                # If the energy is in the bin...
                if ekLIS >= b[0] and ekLIS < b[1]:
                    bIndex = ekLISBins.index(b)

                    # ...update count for the bin: n(E^LIS \in bin)
                    binCounts[bIndex] = binCounts[bIndex] + 1
                    # ...update count for bin given ek: n(E^LIS \in bin | E)
                    probs[bIndex][1][ekIndex] = probs[bIndex][1][ekIndex] + 1

    # Normalize n(E^LIS \in bin | E) by n(E^LIS \in bin) to get conditional probabilities
    for prob, count in zip(probs, binCounts):
        probs[1] = map(truediv, probs[1], count)

    return probs

#############################################################################################################

# Load all runs
eks, nRuns, ekLISsList = loadEnergies('hyades/me/alt0/')


