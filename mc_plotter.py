#!/usr/bin/env python

"""
Plots points in <sigma v>, M_chi space corresponding to gamma ray lines.
"""

import matplotlib
import numpy
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import os, collections, string, math

### Data structure for storing MC samples.
Sample = collections.namedtuple('Sample', ['p_obs', 'theta_obs', 'phi_obs', 'theta_dev'])

### Reads all samples in samplepath directory.  Stores in a dictionary with keys for each momentum/energy that index lists of Samples.
def getdata(samplepath):
    # Get MC data points
    samplefiles = [os.path.splitext(f)[0] for f in os.listdir(samplepath) if f.endswith('.xml')]
    samples = {float(string.split(s)[1].replace('_', '.')): [] for s in samplefiles}

    # Extract data from XML files created by Mathematica
    for sample in samplefiles:
        root = ET.parse(samplepath + '/' + sample + '.xml').getroot()
        # Get sample momentum
        p = float(root.attrib['p_obs'])

        for s in root.findall('sample'):
            samples[p].append(Sample(float(p), float(s.attrib['theta_obs']),
                float(s.attrib['phi_obs']), float(s.attrib['theta_dev'])))

    return samples

### Plots average deviation angle as a function of p_obs.  degs specifies whether to use radians or degrees
def plotavgdev(samples, degs = True):
    conversion = 180.0 / math.pi if degs else 1.0
    avgs = {}
    errors = {}

    for p in samples:
        # Unzip the list of Samples to compute the average of theta_div
        avgs[p] = sum(zip(*samples[p])[3]) / len(samples[p])
        errors[p] = numpy.std(zip(*samples[p])[3], ddof = 1) / len(samples[p])

    fig = plt.figure()
    avgax = fig.add_subplot(1, 1, 1)

    for p in avgs:
        avgax.scatter(p, avgs[p] * conversion)
        # Could compute errors, but they are smaller than the data points!!!
        #avgax.errorbar(p, avgs[p] * conversion, yerr = errors[p] * conversion)
        print('p = ' + str(p) + ' GeV, theta_dev = ' + str(avgs[p] * conversion))

    # Make the graph look good
    avgax.set_xlim([10, 500])
    avgax.set_ylim([0, 100])
    avgax.set_xscale("log")
    avgax.set_xlabel(r'$p_{obs} (GeV)$')
    avgax.set_ylabel(r'$\theta_{dev} ($' + ('degrees)' if degs else 'radians)'))
    avgax.set_title(r'Average deviation angle as a function of observed momentum', fontsize='medium')

    plt.show()

    ## Resize and save
    #plt.gcf().set_size_inches(6.0, 4.0)
    #plt.gcf().savefig('final/tau_limits/' + theory + '_tau.pdf', bbox_inches='tight')

### Main
if __name__ == "__main__":
    samplepath = 'mc_data/'
    samples = getdata(samplepath)
    plotavgdev(samples)



