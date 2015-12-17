#!/usr/bin/env python

'''
Alters the first line of files to contain the trajectory's energy, initial spatial location and some comments.
NO LONGER NEEDED.  Output file format has been altered!
'''

from os.path import basename
from os import listdir
from os.path import isfile, join
import shutil

runfiles = [f for f in listdir('./') if isfile(join('./', f))]

filestochange = [rf for rf in runfiles]

for f in filestochange:
    with open(f, 'r') as rf:
        with open(f[:-7] + '.csv', 'w') as wf:
            rf.readline() # and discard
            wf.write('# Run exit points.  Columns are r (AU), th (rad), ph (rad), ek (GeV), s (s).\n'
                    + '# First line contains initial point of trajectory:\n'
                    + '1.0,1.5707963267948966,0,' + str(f[:-7]) + ',0.0\n')
            shutil.copyfileobj(rf, wf)



