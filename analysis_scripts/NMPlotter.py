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
Simplex = namedtuple('Simplex', ['xs', 'ys', 'ds'])

if __name__ == "__main__":
    # Load the data
    fDir = "/Users/acoogan/Dropbox/heliosphere_anisotropy/nmdata/"
    fName = "nmdata1.csv"
    rawData = np.loadtxt(fDir + fName, delimiter = ",")
    xList = rawData[:,0] * np.sin(rawData[:,1]) * np.cos(rawData[:,2])
    yList = rawData[:,0] * np.sin(rawData[:,1]) * np.sin(rawData[:,2])
    dList = rawData[:,3]

    # Group simplices by timestep
    simplices = []
    for i in range(0, len(dList) / 3):
        simplices.append(Simplex(xs = xList[i+0:i+3].tolist(), ys = yList[i+0:i+3].tolist(),
            ds = dList[i+0:i+3].tolist()))

    fig = plt.figure()
    # Make subplot for each timestep
    nCols = 5

    # Make a color map to shade simplices at each timestep
    norm = matplotlib.colors.Normalize(vmin=0, vmax=len(simplices))
    cmap = cm.get_cmap("viridis")
    colorizer = cm.ScalarMappable(norm = norm, cmap = cmap)

    # Plot the simplices
    for i in range(0, len(simplices)):
        ax = fig.add_subplot(len(simplices) / nCols, nCols, i + 1)
        # Not sure why I need the second tolist() call, but it seems to round if I don't include it.
        ax.plot(simplices[i].xs + [simplices[i].xs[0]],
                simplices[i].ys + [simplices[i].ys[0]], color = colorizer.to_rgba(i))

        # Figure out points' distances
        maxIdx = simplices[i].ds.index(max(simplices[i].ds))
        minIdx = simplices[i].ds.index(min(simplices[i].ds))
        midIdx = [j for j in range(0, len(simplices[i].ds)) if j != maxIdx and j != minIdx][0]

        # Show how far each point is
        print simplices[i].xs
        ax.plot(simplices[i].xs[maxIdx], simplices[i].ys[maxIdx], '.', markersize = 16,
                color = colorizer.to_rgba(i))
        ax.plot(simplices[i].xs[midIdx], simplices[i].ys[midIdx], '.', markersize = 8,
                color = colorizer.to_rgba(i))
        ax.plot(simplices[i].xs[minIdx], simplices[i].ys[minIdx], '.', markersize = 4,
                color = colorizer.to_rgba(i))

        # Plot optimal point
        ax.plot(1.0 * np.sin(np.pi / 4.0), np.pi / 4.0, 'xr')

        ax.set_xlim(min([min(s.xs) for s in simplices]) - 0.1, max([max(s.xs) for s in simplices]) + 0.1)
        ax.set_ylim(min([min(s.ys) for s in simplices]) - 0.1, max([max(s.ys) for s in simplices]) + 0.1)

        #ax[len(simplices) - 1].set_title('Nelder-Mead simplex')
        if i == len(simplices):
            ax.set_xlabel('r (au)')
            ax.set_ylabel('ph (rad)')

    plt.show()


