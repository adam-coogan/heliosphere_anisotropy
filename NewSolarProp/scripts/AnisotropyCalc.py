#!/usr/bin/env python

"""
Computes the anisotropy using the streaming flux.  Requires a set of files taken at and near the point where
the flux is to be computed.
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
from SDEModulate import rigidity, jLISLangner, getJ, me, rundataPath

# Directory containing run data relevant to numerical derivatives

"""
Basic constants and functions related to propagation quantities.
TODO: figure these out from the files!
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
rRefLambda = 1.0 # au
P0 = 1.0 # GV

def momentum(ek):
    return np.sqrt(ek * (ek + 2.0 * me)) # GeV

def velocity(ek):
    return 0.002 * np.sqrt(ek * (ek + 2.0 * me)) / (ek + me) # AU / s

def getB(r, th, ph):    
    bFact = Ac * B0 * (rRefB / r)**2
    return {'r': bFact, 'th': 0.0, 'ph': -bFact * (r - rSun) * Omega / Vsw * np.sin(th)}

def getK(r, th, ph, ek):
    # Diffusion tensor
    k = {'rr': 0.0, 'thth': 0.0, 'phph': 0.0, 'rth': 0.0, 'rph': 0.0, 'thr': 0.0, 'thph': 0.0, 'phr': 0.0, \
            'phth': 0.0}

    # Set up convenience variables
    p = momentum(ek)
    v = velocity(ek)
    B = getB(r, th, ph)

    # Asymmetric part of k in spherical coordinates
    kAFact = p * v / (3 * qe * (B['r']**2 + B['ph']**2))
    k['rth'] = kAFact * B['ph']
    k['thr'] = -kAFact * B['ph']
    k['thph'] = kAFact * B['r']
    k['phth'] = -kAFact * B['r']

    # Symmetric part of k
    kpar = v / 3 * lambda0 * (1 + r / rRefLambda)
    if p > P0:
        kpar = kpar * p / P0
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

def getjDerivs(baseFName, drFName, dthFName, dphFName, deFName, dirName):
    """
    Arguments: file names with trajectories measured to near point at which anisotropy is to be computed.  
    Directory relative to rundataPath containing the files.
    Returns: numerical approximations of dj/dr, dj/dth, dj/dph, dj/dT at Earth.  Also gives x, the Earth
    point.
    """
    # Get j at Earth
    xEarth, jEarth, _ = getJ(dirName + '/' + baseFName)
    # Get j at each of the points shifted spatially/spectrally from the one we care about
    x_dr, j_dr, _ = getJ(dirName + '/' + drFName)
    x_dth, j_dth, _ = getJ(dirName + '/' + dthFName)
    x_dph, j_dph, _ = getJ(dirName + '/' + dphFName)
    x_de, j_de, _ = getJ(dirName + '/' + deFName)

    # Compute derivatives
    dj_dr = (j_dr - jEarth) / (x_dr['r0'] - xEarth['r0'])
    dj_dth = (j_dth - jEarth) / (x_dth['th0'] - xEarth['th0'])
    dj_dph = (j_dph - jEarth) / (x_dph['ph0'] - xEarth['ph0'])
    dj_de = (j_de - jEarth) / (x_de['e0'] - xEarth['e0'])

    return (xEarth, jEarth, {"dj_dr": dj_dr, "dj_dth": dj_dth, "dj_dph": dj_dph, "dj_de": dj_de})

def getAnisotropy(baseFName, drFName, dthFName, dphFName, deFName, dirName):
    """
    Computes anisotropy at position specified in the provided file.  Files required for computing derivatives
    must exist.  .csv extension should not be included in their names.
    (See eq 1-2b in Jokipii and Kopriva 1979.)
    """
    # Set up convenience variables
    xEarth, jEarth, dj = getjDerivs(baseFName, drFName, dthFName, dphFName, deFName, dirName)
    e0 = xEarth['e0']
    k = getK(xEarth['r0'], xEarth['th0'], xEarth['ph0'], e0)

    # Streaming flux sf
    sfFact = 4.0 * np.pi / velocity(e0)
    sf = {'r': sfFact, 'th': sfFact, 'ph': sfFact}

    # Angular components
    sf['th'] = -sf['th'] * (k['thr'] * dj['dj_dr'] + k['thth'] * dj['dj_dth'] + k['thph'] * dj['dj_dph'])    
    sf['ph'] = -sf['ph'] * (k['phr'] * dj['dj_dr'] + k['phth'] * dj['dj_dth'] + k['phph'] * dj['dj_dph'])
    # r component has solar wind contribution
    sf['r'] = sf['r'] * (-k['rr'] * dj['dj_dr'] - k['rth'] * dj['dj_dth'] - k['rph'] * dj['dj_dph'] \
            + Vsw * (jEarth - 1.0 / (3.0 * velocity(e0)) \
                * (jEarth * (1 + me**2 / (me + e0)**2) \
                    + e0 * (e0 + 2 * me) / (e0 + me) \
                    * (dj['dj_de'] - jEarth / velocity(e0) * me**2 / (velocity(e0) * (me + e0)**3)))))

    # Compute the anisotropy from the streaming flux
    delta = {'r': 3.0 / (4.0 * np.pi * jEarth) * sf['r'], \
            'th': 3.0 / (4.0 * np.pi * jEarth) * sf['th'], \
            'ph': 3.0 / (4.0 * np.pi * jEarth) * sf['ph']}

    return delta

#############################################################

if __name__ == "__main__":
    #meE0s, meJs, meSigmas = map(list, zip(*[getJ('me/alt0/' + os.path.splitext(f)[0]) \
    #        for f in listdir(rundataPath + 'me/alt0/')]))
    #meE0s = [x0['e0'] for x0 in meE0s]

    print 'Computing anisotropy...'

    delta = getAnisotropy('1.0GeV', '1.0GeV_r1.1AU', '1.0GeV_th0.55pi', '1.0GeV_ph0.2pi', '1.05GeV',
            'me/highstats/alt0/')

    print 'Anisotropy delta = ' + str(delta)


