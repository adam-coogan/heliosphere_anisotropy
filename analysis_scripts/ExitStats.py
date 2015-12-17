#!/usr/bin/env python

"""
Functions for performing statistical analyses dependent only on first and last points in MCs
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt
import itertools

if __name__ == "__main__":
    # Load exit data
    runDir = 'rundata/'
    #fName = 'hyades/alt0_0.1GeV'
    #fName = 'strauss_run5'
    #fName = 'straussunits_alt0_run5'
    fName = 'macbookpro/precision_test/strauss_alt0_runs'

    rawData = np.loadtxt(runDir + fName + '.csv', delimiter = ',')

    # Parse out the coordinates, energies and times
    exitData = {'r': rawData[:,0], \
            'th': 180.0 / np.pi * rawData[:,1], \
            'ph': 180.0 / np.pi * rawData[:,2], \
            'ek': rawData[:,3], \
            's': rawData[:,4] / (24.0 * 3600.0)}


    # Toggle log scale on the theta and energy plots' y axes
    logScale = True

    # Create figure for plotting polar angle, energy and time at exit
    fig = plt.figure(figsize = (16, 6))
    thAx = fig.add_subplot(1, 3, 1)
    ekAx = fig.add_subplot(1, 3, 2)
    sAx = fig.add_subplot(1, 3, 3)

    # Fill the histograms
    thData = []
    sData = []
    ekData = []
    for i in range(len(exitData['th'])):
        thData.append(exitData['th'][i])
        sData.append(exitData['s'][i])
        ekData.append(exitData['ek'][i])

    # Get expected propagation time
    propTime = sum(exitData['s']) / float(len(exitData['s']))
    propTimeError = sum([(s - propTime)**2 for s in exitData['s']]) / (len(exitData['s']) - 1)
    propTimeError = np.sqrt(propTimeError / float(len(exitData['s'])))
    #print('Expected propagation time: ' + str(propTime) + ' +/- ' + str(propTimeError) + ' days')
    #print('Should be ~300 days for Ac > 0, ~180 days for Ac < 0')

    # Fill the histograms
    thAx.hist(thData, bins = np.arange(0.0, 180.0 + 6.0, 6.0), log = logScale)
    ekAx.hist(ekData, bins = np.arange(0.0, 3.0 + 0.1, 0.1), log = logScale)
    sMin = 0.0
    sMax = 1000.0
    sWidth = 10.0
    sAx.hist(sData, bins = np.arange(sMin, sMax, sWidth))

    # Scale axes
    axScale = len(exitData['th']) / 10000.0

    # Set up axes and labels
    thAx.set_ylim(0.9, axScale*6000.0)
    ekAx.set_ylim(0.9, axScale*2500.0)
    sAx.set_ylim(0.9, axScale*650.0) #(700.0 if fName[1:3] == 'lt' else 500.0))
    thAx.set_ylabel('Number of occurences')
    ekAx.set_xlim(0.1, 3.0)
    ekAx.set_xlabel('Energy (GeV)')
    ekAx.set_title(r'Strauss\' code ($A_c < 0$, $\langle s \rangle = ' + str(round(propTime, 3)) + ' \pm '
            + str(round(propTimeError, 3)) + '$ days)')
    thAx.set_xlim(0.0, 180.0)
    thAx.set_xlabel('Polar angle (degrees)')
    sAx.set_xlim(0.0, 1000.0)
    sAx.set_xlabel('Propagation time (days)')

    plt.show()
    #fig.savefig('../figures/' + fName + '.pdf')


