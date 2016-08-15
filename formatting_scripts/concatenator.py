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

# Number of lines at the beginning of the file that are shared in common by runs generated with the same
# parameters
sharedLines = 3

# Specify directory in which to concatenate files using the first command line argument
rundir = sys.argv[1]

# Get list of run files
runfiles = [f for f in listdir(rundir) if isfile(join(rundir, f))]
runfiles = [re.search('(.*)_run.*csv', rf) for rf in runfiles]

# Figure out the unique base file names
uniqueNames = set([rf.group(1) for rf in runfiles if rf != None])
runfiles = [rf.group(0) for rf in runfiles if rf != None]

# Concatenate the files
for un in uniqueNames:
    # Find items which have the correct base name
    unFiles = [re.search(un + '_run.*csv', rf) for rf in runfiles]
    unFiles = [rundir + uf.group(0) for uf in unFiles if uf != None]

    # Store first lines, which contain comments and initial point
    firstLines = ''
    with open(unFiles[0], 'r') as f:
        for i in range(1, sharedLines + 1):
            firstLines = firstLines + f.readline()

    # Make output file
    outfilename = rundir + un + '.csv'
    with open(outfilename, 'w') as outfile:
        # Write first line containing initial point information to the outfile
        outfile.write(firstLines)

        # Concatenate contents from files that have the same unique basename
        for uf in unFiles:
            with open(uf, 'r') as readfile:
                # Discard first lines
                for i in range(1, sharedLines + 1):
                    readfile.readline()

                # Concatenate contents from file being read into the output file
                shutil.copyfileobj(readfile, outfile)
                outfile.write('\n')

        print 'Wrote ' + outfilename

print 'Done concatenating runs!'


