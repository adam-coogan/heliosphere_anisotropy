#!/usr/bin/env python

"""
Concatenates files generated with the same parameters (in this case, energy).
"""

import glob
from os import listdir
from os.path import basename, isfile, join
import re
import shutil
import sys

# Specify directory in which to concatenate files using the first command line argument
rundir = sys.argv[1]

# Get list of run files
runfiles = [f for f in listdir(rundir) if isfile(join(rundir, f))]
runfiles = [re.search('(.*)GeV_run.*csv', rf) for rf in runfiles]

# Figure out the unique base file names
uniqueNames = set([rf.group(1) for rf in runfiles if rf != None])
runfiles = [rf.group(0) for rf in runfiles if rf != None]

# Concatenate the files
for un in uniqueNames:
    # Find items which have the correct base name
    unFiles = [re.search(un + 'GeV_run.*csv', rf) for rf in runfiles]
    unFiles = [rundir + uf.group(0) for uf in unFiles if uf != None]

    # TODO: write code to remove and store the first line from each file for a given energy
    """
    with open(unFiles[0], 'r') as f:
        print unFiles[0]

        firstLine = f.readline()
        print firstLine
    """

    # Make output file
    outfilename = rundir + un + 'GeV.csv'
    with open(outfilename, 'w') as outfile:
        # Write first line containing initial point information to the outfile
        outfile.write('# Initial point: r = 1.0 AU, th = 1.5707963268 rad, ph = 0.0 rad, ek = ' + un \
                + ', s = 0.0 s\n')

        # Concatenate contents from files that have the same unique basename
        for uf in unFiles:
            with open(uf, 'r') as readfile:
                shutil.copyfileobj(readfile, outfile)
                outfile.write('\n')

        print 'Wrote ' + outfilename

print 'Done concatenating runs!'


