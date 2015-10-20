#!/usr/bin/env python

"""
Generates and submits jobs.
"""

from math import floor
import subprocess

def roundedStr(num, n):                                                                                       
    """                                                                                                       
    Returns a string containing num rounded to n significant figures.                                         
    """                                                                                                       
    return '{number:.{digits}g}'.format(number = num, digits = n)

#########################################################################################

# Job will be submitted to this queue
queue = 'DEFAULT'
queueStr = '\n'
if queue != 'DEFAULT':
    queueStr = '#PBS -q ' + queue + '\n'

# PBS script contents.  Blank strings will be replaced by the run's energy.
pbsContents = ['#!/bin/bash\n\n',
        '#PBS -N solarproparray', '', 'GeV\n',
        queueStr,
        '#PBS -t 1-5\n',
        '#PBS -l walltime=0:15:00\n\n',
        'cd $PBS_O_WORKDIR\n\n',
        './StraussUnits ', '', ' 10000 hyades/me/alt0/', '', 'GeV_run$PBS_ARRAYID']

pbsName = 'solarproparray.pbs'

# Energies
#energies = [1e0 * 10 ** (a / 19.0) for a in range(19)]
#energies = [roundedStr(e, 3) for e in energies[1::2]]
energies = [1.37, 1.55, 1.75, 1.99, 2.25, 2.56, 2.88, 3.29, 4.12, 4.73, 5.41, 6.08, 6.95, 7.80, 8.93]
energies = [e * 1e-2 for e in energies] + [e * 1e-1 for e in energies] + [e * 1e0 for e in energies]
energies = [roundedStr(e, 3) for e in energies]

# Submit a small number of jobs at a time
nJobs = 5
startIndices = [nJobs * si for si in range(0, 1 + int(floor(len(energies) / nJobs)))]

# This labels sets of nJobs jobs.
nBatch = 9 # Run with this
batchIndices = [j for j in range(startIndices[nBatch], startIndices[nBatch] + nJobs) if j < len(energies)]
print('Running batch ' + str(nBatch + 1) + ' out of ' + str(len(startIndices)) + '.')
print('Consists of ' + str(len(batchIndices)) + ' job arrays.')

# Submit a job for each energy in the batch
for i in batchIndices:
    e = energies[i]
    curContents = ''.join([e if s == '' else s for s in pbsContents])

    # Write script contents to file
    with open(pbsName, 'w') as pbsScript:
        pbsScript.write(curContents)

    # Run the script
    subprocess.Popen(['qsub', pbsName], cwd = '/home/acoogan/pfs/heliosphere_anisotropy/NewSolarProp/pbs/')
    print('Submitted jobs for E = ' + e + ' GeV.')


