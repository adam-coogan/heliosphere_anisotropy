#!/usr/bin/env python

"""
Generates and submits jobs.  Must be run in the pbs directory.  TODO: change this!
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

# Configuration file.  Blank strings will be replaced by the run's energy.
configBase = ['# Minimal parameter file\n',
        'Ac=-1\n',
        'm=0.000511\n',
        'charge=-1\n',
        'ds=374\n',
        'ek0=',
        '']

configName = 'config/params.config'

# PBS script contents.  Blank strings will be replaced by the run's energy.
pbsBase = ['#!/bin/bash\n\n',
        '#PBS -N ', '', 'GeV\n',
        queueStr,
        '#PBS -t 1-5\n',
        '#PBS -l walltime=0:30:00\n\n',
        'cd $PBS_O_WORKDIR\n\n',
        '# Ensures proper dynamic linking to boost libraries\n',
        'LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/pfs/sw/serial/gcc/boost-1.57.0/lib/\n',
        'export LD_LIBRARY_PATH\n\n',
        './Helios ',
        '-c ' + configName + ' ',
        '-o rundata/alt0/${PBS_JOBNAME%-$PBS_ARRAYID}\_run$PBS_ARRAYID\.csv ',
        '-n 4000']

pbsName = 'heliosarray.pbs'

# Energies
#prefactor = 0.01 # Changes order of magnitude of energy
energies = [roundedStr(10 ** (a / 10.0), 3) for a in range(-20, 11)]

# Submit a small number of jobs at a time
nJobs = 5
startIndices = [nJobs * si for si in range(0, 1 + int(floor(len(energies) / nJobs)))]

# This labels sets of nJobs jobs.
nBatch = 3 # Run with this
batchIndices = [j for j in range(startIndices[nBatch], startIndices[nBatch] + nJobs) if j < len(energies)]
print('Running batch ' + str(nBatch + 1) + ' out of ' + str(len(startIndices)) + '.')
print('Consists of ' + str(len(batchIndices)) + ' job arrays.')

# Submit a job for each energy in the batch
for i in batchIndices:
    e = energies[i]

    # Update pbs script and configuration file to use the current energy.
    curPBS = ''.join([e if s == '' else s for s in pbsBase])
    curConfig = ''.join([e if s == '' else s for s in configBase])

    # Write script and config contents.
    with open(pbsName, 'w') as pbsScript:
        pbsScript.write(curPBS)

    with open(configName, 'w') as config:
        config.write(curConfig)

    # Run the script.  call blocks until the command has finished running to avoid errors.
    subprocess.call(['qsub', pbsName], cwd = '/home/acoogan/pfs/heliosphere_anisotropy/helios/pbs/')
    print('Submitted Helios jobs for E = ' + e + ' GeV.')


