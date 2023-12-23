"""
Date: 22 December 2023
Author: Egheosa Ogbomo

Script to automatically generate files/directories for similar conditions

- Have a source directory for files to be used in every 
- Make empty ovito_datafiles and processed folder
- Be able to toggle between generated LAMMPS file (equil or comp)
- Make QSUB files for each condition being simulated
- Make QSUB files depending on HPC being used

Will require that all experiments that we are generating files for are
at the same 'stage', will need to check this.
"""

######### Imports #########
import pandas as pd 
import numpy as np
import math
import csv
import os
import collections
from pathlib import Path
import subprocess
from copy import deepcopy

def runcmd(cmd, verbose = False, *args, **kwargs):
    #bascially allows python to run a bash command, and the code makes sure 
    #the error of the subproceess is communicated if it fails
    process = subprocess.run(
        cmd,
        text=True,
        shell=True)

######### File Generation Parameters ###########

STAGE = 'First'
STARTINGDIR = 'F:/PhD/TCPDecompositionExperiments/Completed/RoughSurfaces/Iron_Oxide/8nm/' # Home directory to launch generation from
STARTINGDIR = 'F:/PhD/TCPDecompositionExperiments/Completed/DCMP_Solvent/' # Home directory to launch generation from
SOURCEDIR = os.path.join(STARTINGDIR, 'SourceDir') # Directory where enabler files are 
System = 'DCMP_Fe_48_TCP_Mixed' # System being simulated 
EquilTime = '800000' # Equilibration time
CompTime = '4000000' # Compression and shear time
Wall_V = '0.0001' # Wall velocity
Wall_Z = '' # Wall thickness
# Atom types
HType =  '5' 
FeType = '1'
OType = '3'
PType = '2'
CType = '4'
Fix_Z = '1.2' # Fixed layer thickness
Thermo_Z = '2.4' # Thermostat layer thickness
ReaxFFTyping = 'Fe P O C H' # Order of elements for ReaxFF command
Temperatures = ['600K'] # Temperatures to be simulated
Pressures = ['1GPa', '2GPa', '3GPa', '4GPa', '5GPa'] # Pressures to be simulated
EquilPress = '10' # Equilibration Temperature, in MPa
Safezone = '800' # System memory parameter
Mincap = '1800' # System memory parameter
RestartFileFreq = '100' 

############# Calling the function #########################

os.chdir(STARTINGDIR) # Go to starting directory

for Temp in Temperatures:
    for Press in Pressures:
        os.chdir(os.path.join(STARTINGDIR, Temp, Press))
        #for dir in os.listdir():
        #    assert os.path.isdir(os.path.join(os.getcwd(), dir)), 'There are non-directory files in this directory'

        RestartList = [x for x in os.listdir() if 'Restart' in x]
        
        if len(RestartList) == 0:
            runcmd('mkdir Restart_1')
            runcmd(f'copy {SOURCEDIR} /Restart_1')

        else:
            RestartNumbers = [x.split('_')[-1] for x in RestartList]
            RestartNumber = sorted(RestartNumbers)[-1]
            print(RestartNumber)
            runcmd('mkdir Restart_1')
            runcmd(f'copy {SOURCEDIR} /Restart_1')



