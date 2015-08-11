#!/usr/bin/env python

"""
Functions for performing statistical analyses dependent only on first and last points in MCs
"""

import SolarPropContainers
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import xml.etree.ElementTree as ET
import os, sys

## Path to root directory for MC runs
#mcPath = sys.argv[1]

###
# Parses all XML files in the specified directory.  Recursively traverses subdirectories
# Arguments:
# Returns:
###
def parseTrajs(mcPath):
    # Store trajectories
    trajectories = []

    # Traverse the directories in mcPath to find all XML files in the tree
    for (dirPath, dirNames, fileNames) in os.walk(mcPath):
        for fileName in fileNames:
            # Parse each XML file
            if (len(fileName) > 4 and fileName[-4:] == ".xml"):
                traj = parseTraj(os.path.join(dirPath, fileName))
                
                # Make sure the trajectory didn't hit the sun
                if traj != None:
                    trajectories.append(parseTraj(os.path.join(dirPath, fileName)))
    
    return trajectories

###
# Parses a SolarProp trajectory.
# Arguments:
#   trajPath: the XML file"s path.
# Returns: a MinTrajectory containing the trajectory"s first and last points.
###
def parseTraj(trajPath):
    # Load up the XML
    root = ET.parse(trajPath).getroot()

    if root.find("boundary").text.strip() == "heliopause":
        # Put points into a list.  Get the first and last elements.
        points = root.find("points").findall("point")
        first = SolarPropContainers.PPoint(float(points[0].find("r").text), float(points[0].find("th").text),
                float(points[0].find("phi").text), float(points[0].find("e").text),
                float(points[0].find("t").text))

        last = SolarPropContainers.PPoint(float(points[-1].find("r").text), float(points[-1].find("th").text),
                float(points[-1].find("phi").text), float(points[-1].find("e").text),
                float(points[-1].find("t").text))

        return SolarPropContainers.MinTrajectory(first, last)
    else:
        return None

if __name__ == "__main__":
    # Path to root directory for MC runs
    mcPath = sys.argv[1]

    # Load all the trajectories
    trajectories = parseTrajs(mcPath)

    print("Analyzing " + str(len(trajectories)) + " MC trajectories...")

    # Count the number of galactic particles
    fromGalaxy = 0
    # Store exit energies, polar angles and propagation times
    thBelow = []
    thAbove = []
    eBelow = []
    eAbove = []
    tBelow = []
    tAbove = []

    for traj in trajectories:
        # Make sure the trajectory ended at the heliopause
        if (traj.first.r >= 140.0):
            fromGalaxy += 1

            # Check whether the trajectory is above or below the specified cut and fill the correct lists.
            if (True):
                thBelow.append(traj.first.th * 180.0 / np.pi)
                eBelow.append(traj.first.e)
                tBelow.append(traj.last.t / (24.0 * 60.0 * 60.0))
            else:
                thAbove.append(traj.first.th * 180.0 / np.pi)
                eAbove.append(traj.first.e)
                tAbove.append(traj.last.t / (24.0 * 60.0 * 60.0))

    # Create figure
    fig = plt.figure()
    polarAx = fig.add_subplot(1, 3, 1)
    energyAx = fig.add_subplot(1, 3, 2)
    propTimeAx = fig.add_subplot(1, 3, 3)

    # Create histograms
    polarAx.hist(thBelow, bins = np.arange(0.0, 180.0 + 6.0, 6.0), log = True)
    energyAx.hist(eBelow, bins = np.arange(0.0, 3.0 + 0.1, 0.1), log = True)
    # Propagation times
    tWidth = 10.0
    tMin = 0.0
    tMax = 1000.0
    propTimeAx.hist(tBelow, bins = np.arange(tMin, tMax, tWidth))
    timeCounts, timeBins = np.histogram(tBelow + tAbove, bins = np.arange(tMin, tMax, tWidth))

    # Get probabilities for each bin
    timeProbs = timeCounts / (fromGalaxy + 0.0)
    # Find centers of bins
    timeBins = timeBins[1:]
    timeBins -= tWidth / 2.0
    # Get expected propagation time
    propTime = sum(t * p for t, p in zip(timeBins, timeProbs))

    # Find standard deviation of the mean to get error estimates
    timeError = np.std([float(t) / (24.0 * 60.0 * 60.0) for t in tBelow + tAbove],
            ddof = 1) / np.sqrt(fromGalaxy)

    print("The expected propagation time is " + str(propTime) + " days")

    # Set axes and labels
    polarAx.set_ylim(0.9, 6000.0)
    energyAx.set_ylim(0.9, 9000.0)
    propTimeAx.set_ylim(0.9, 500.0)
    polarAx.set_ylabel("Number of occurences")
    energyAx.set_ylabel("Number of occurences")
    propTimeAx.set_ylabel("Number of occurences")

    energyAx.set_xlim(0.1, 3.0)
    energyAx.set_xlabel("Energy (GeV)")

    polarAx.set_xlim(0.0, 180.0)
    polarAx.set_xlabel("Polar angle (degrees)")

    propTimeAx.set_xlim(0.0, 1000.0)
    propTimeAx.set_xlabel("Propagation time (days)")

    plt.show()


