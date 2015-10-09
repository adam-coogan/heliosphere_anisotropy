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

rundataPath = '/Users/acoogan/Dropbox/heliosphere_anisotropy/NewSolarProp/rundata/hyades/'
generatorDir = '/Users/acoogan/Dropbox/heliosphere_anisotropy/NewSolarProp/bin/'
generator = 'StraussUnits'

"""
Basic constants and functions related to propagation quantities.
"""

# Solar cycle polarity
Ac = -1
# Magnetic field
rRefB = 1.0
rSun = 0.005 # AU
B0 = 5e-9 / np.sqrt(1 + (rRefB - rSun)**2) * 4.485e10 # GV / AU
Omega = 2.0 * np.pi / (25.4 * 24.0 * 3600.0) # rad / s
Vsw = 400 * 6.685e-9 # AU / s
# Electron mass
me = 0.000511
qe = -1
# Diffusion parameters
lambda0 = 0.15 # AU
kperp_kpar = 0.01

def momentum(ek):
    return np.sqrt(ek * (ek + 2.0 * me)) # GeV

def velocity(ek):
    return 0.002 * np.sqrt(ek * (ek + 2.0 * me)) / (ek + me) # AU / s

def rigidity(ek):
    return np.sqrt(ek * (ek + 2.0 * me)) # GV

def getB(r, th, ph):    
    bFact = Ac * B0 * (rRefB / r)**2
    return {bFact, 0.0, -bFact * (r - rSun) * Omega / Vsw * np.sin(th)}

def getK(r, th, ph, ek):
    # Diffusion tensor
    k = {'rr': 0.0, 'thth': 0.0, 'phph': 0.0, 'rth': 0.0, 'rph': 0.0, 'thr': 0.0, 'thph': 0.0, 'phr': 0.0, \
            'phth': 0.0}

    # Set up convenience variables
    p = momentum(ek)
    v = velocity(ek)
    B = getB(r, th, ph, -1)

    # Asymmetric part of k in spherical coordinates
    kAFact = p * v / (3 * qe * (B['r']**2 + B['ph']**2))
    k['rth'] = kAFact * B['ph']
    k['thr'] = -kAFact * B['ph']
    k['thph'] = kAFact * B['r']
    k['phth'] = -kAFact * B['r']

    # Symmetric part of k
    kpar = v / 3 * lambda0 * (1 + r / rRefLambda)
    if p > p0:
        kpar = kpar * p / p0
        kperp = kperp_kpar * kpar

    # Convert symmetric part to spherical coordinates
    tanPsi = Omega * (r - rSun) * np.sin(th) / Vsw
    cosPsi = 1 / np.sqrt(1 + tanPsi**2)
    sinPsi = np.sqrt(1 - cosPsi**2)

    k['rr'] = kpar * cosPsi**2 + kperp * sinPsi**2
    k['phph'] = kpar * sinPsi**2 + kperp * cosPsi**2
    k['thth'] = kperp
    k['rph'] = k['rph'] + (kperp - kpar) * cosPsi * sinPsi
    k['phr'] = k['phr'] + (kperp - kpar) * cosPsi * sinPsi

    return k

"""
Local interstellar spectrum.
"""

def jLISLangner(ek):
    """
    LIS parametrization from Langner 2004, with correction from Strauss' code.  ek is in GeV.
    Returns: j_LIS in units of particles m^-2 s^-1 sr^-1 MeV^-1.
    """
    P = rigidity(ek) / 1.0 # P0 = 1.0 GV

    if P < 0.0026:
        a = 126.067948848145
        b = 0.2567973348983205
        c = 1.95129069032843
        d = 0.0171199701826333

        return 1.7 * (a + c * np.log(P)) / (1.0 + b * np.log(P) + d * np.log(P)**2)
    elif P >= 0.0026 and P < 0.1:
        return 1.7 * ((52.55 + 23.01 * P) / (1.0 + 148.62 * P))**2
    elif P >= 0.1 and P <= 10.0:
        return (1555.89 + 17.36 * P - 3.4e-3 * P**2 + 5.13e-7 * P**3)/(1.0 - 11.22 * P + 7532.93 * P**2
                + 2405.01 * P**3 + 103.87 * P**4)
    elif P > 10.0:
        return 1.7 * np.exp(-0.89 - 3.22 * np.log(P))

def getJ(runName):
    """
    Performs MC integration to find the differential intensity at earth for a given run.
    runName: csv containing trajectory exit points.
    Returns:
    intitialPoint: the point at which the cosmic rays were observed and their energy.
    j: the differential intensity at that point and energy.
    sigma: the error on j.
    """
    # Load exit points
    fullRunName = rundataPath + runName + '.csv'
    rawData = np.loadtxt(fullRunName, delimiter = ',')

    # Parse out the coordinates and energies into a list of tuples of exit points
    exitData = [{'ee': ee, 'the': the, 'phe': phe} for ee, the, phe \
            in zip(rawData[:,3], rawData[:,1], rawData[:,2])]

    # Get initial point.  --->>> TODO: store initial point in first line of csv!!!
    r0 = 1.0
    th0 = np.pi / 2.0
    ph0 = 0.0
    e0 = float(os.path.splitext(os.path.basename(fullRunName))[0][:-3]) # TODO: store in first line of csv
                                                                        # rather than using naming convention
    initialPoint = {'r0': r0, 'th0': th0, 'ph0': ph0, 'e0': e0}

    # Average LIS at exit points
    #print([jLIS(ep['ee']) for ep in exitData][:10])
    jN = sum([jLISLangner(ep['ee']) / rigidity(ep['ee'])**2 for ep in exitData]) \
            / float(len(exitData)) * rigidity(e0)**2

    # Compute error
    varjN = sum([(jLISLangner(ep['ee']) * (rigidity(e0) / rigidity(ep['ee']))**2 - jN)**2 for ep in exitData])
    varjN = varjN * 1.0 / (len(exitData) - 1.0)

    return (initialPoint, jN, np.sqrt(varjN / len(exitData)))

#############################################################

"""
# Tests
x0, j, sigma = getJ('hyades/me/alt0_0.1GeV')
print('e0 = ' + str(x0['e0']) + ' GeV, j(x0) = ' + str(round(j, 5)) + ' +/- ' + str(round(sigma, 5))
        + ' MeV^-1 s^-1 sr^-2 m^-2')

x0, j, sigma = getJ('macbookpro/precision_test/strauss_alt0_runs')
print('e0 = ' + str(x0['e0']) + ' GeV, j(x0) = ' + str(round(j, 5)) + ' +/- ' + str(round(sigma, 5))
        + ' MeV^-1 s^-1 sr^-2 m^-2')
"""

# Use runs to compute modulated spectrum
stE0s, stJs, stSigmas = map(list, zip(*[getJ('strauss/alt0/' + os.path.splitext(f)[0]) \
        for f in listdir(rundataPath + 'strauss/alt0/')]))
stE0s = [x0['e0'] for x0 in stE0s]

meE0s, meJs, meSigmas = map(list, zip(*[getJ('me/alt0/' + os.path.splitext(f)[0]) \
        for f in listdir(rundataPath + 'me/alt0/')]))
meE0s = [x0['e0'] for x0 in meE0s]

#############################################################

# Plot spectra

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Strauss' modulated spectrum
ax.errorbar(stE0s, stJs, yerr = stSigmas, fmt = '.', label='Strauss')
ax.errorbar(meE0s, meJs, yerr = meSigmas, fmt = '.', label='Me')

# Unmodulated spectrum
eksLangner = 10.0 ** np.linspace(-3.0, 2.0, 1000) # Energies in GeV
ax.plot(eksLangner, [jLISLangner(ek) for ek in eksLangner])

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-3, 1e1)
ax.set_ylim(1e-4, 1e4)
ax.set_xlabel(r'T (GeV)')
ax.set_ylabel(r'Intensity ($\mathrm{MeV}^{-1}\ \mathrm{s}^{-1}\ \mathrm{sr}^{-2}\ \mathrm{m}^{-2}$)')
ax.set_title(r'Cosmic ray spectrum modulation ($A_c < 0$)')
plt.legend(loc='upper right')

plt.show()


