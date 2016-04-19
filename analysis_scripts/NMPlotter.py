#!/usr/bin/env python

'''
Plots simplices generated at each step by NM implementation in Wavy3D's getLHCS() function.
'''

from collections import namedtuple
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Represents a simplex and the function's value at a given timestep
Simplex = namedtuple('Simplex', ['rs', 'phs', 'ds'])

if __name__ == "__main__":
    # Load the data
    fDir = "/Users/acoogan/Dropbox/heliosphere_anisotropy/nmdata/"
    fName = "nmdata.csv"
    rawData = np.loadtxt(fDir + fName, delimiter = ",")
    rList = rawData[:,0]
    #thList = rawData[:,1]
    phList = rawData[:,2]
    dList = rawData[:,3]

    # Group simplices by timestep
    simplices = []
    for i in range(0, len(rList) / 3):
        simplices.append(Simplex(rs = rList[i+0:i+3], phs = phList[i+0:i+3], ds = dList[i+0:i+3]))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # Make a color map to shade simplices at each timestep
    norm = matplotlib.colors.Normalize(vmin=0, vmax=len(simplices))
    cmap = cm.Blues
    colorizer = cm.ScalarMappable(norm = norm, cmap = cmap)

    # Plot the simplices
    for i in range(0, len(simplices)):
        # Not sure why I need the second tolist() call, but it seems to round if I don't include it.
        ax.plot(simplices[i].rs.tolist() + [simplices[i].rs.tolist()[0]],
                simplices[i].phs.tolist() + [simplices[i].phs.tolist()[0]], color = colorizer.to_rgba(i))

    # Plot optimal point
    ax.plot(1.0 * np.sin(np.pi / 4.0), np.pi / 4.0, 'xr')

    ax.set_xlim(min([min(s.rs) for s in simplices]) - 0.1, max([max(s.rs) for s in simplices]) + 0.1)
    ax.set_ylim(min([min(s.phs) for s in simplices]) - 0.1, max([max(s.phs) for s in simplices]) + 0.1)
    ax.set_title('Nelder-Mead simplex')
    ax.set_xlabel('r (au)')
    ax.set_ylabel('ph (rad)')

    plt.show()


