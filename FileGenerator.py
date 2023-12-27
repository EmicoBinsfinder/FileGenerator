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
import HelperFunctions as HF

def runcmd(cmd, verbose = False, *args, **kwargs):
    #bascially allows python to run a bash command, and the code makes sure 
    #the error of the subproceess is communicated if it fails
    process = subprocess.run(
        cmd,
        text=True,
        shell=True)

######### File Generation Parameters ###########

STAGE = 'First'
#STARTINGDIR = 'F:/PhD/TCPDecompositionExperiments/Completed/RoughSurfaces/Iron_Oxide/8nm/' # Home directory to launch generation from
STARTINGDIR = 'F:/PhD/TCPDecompositionExperiments/Completed/RoughSurfaces/Iron/2nm/' # Home directory to launch generation from
SOURCEDIR = os.path.join(STARTINGDIR, 'SourceDir') # Directory where enabler files are 
System = '90_TCP_Fe_0_2nm' # System being simulated 
EquilTime = '800000' # Equilibration time
CompTime = '4000000' # Compression and shear time
Wall_V = '0.0001' # Wall velocity
Wall_Z = '14' # Wall thickness
# Atom types
HType =  '5' 
FeType = '1'
OType = '3'
PType = '2'
CType = '4'
Fix_Z = '1.2' # Fixed layer thickness
Thermo_Z = '2.4' # Thermostat layer thickness
ReaxFFTyping = 'Fe P O C H' # Order of elements for ReaxFF command
Temperatures = ['500K', '600K', '700K'] # Temperatures to be simulated
Pressures = ['1GPa', '2GPa', '3GPa', '4GPa', '5GPa'] # Pressures to be simulated
EquilPress = '10' # Equilibration Temperature, in MPa
EquilTemp = '300'
Safezone = '80' # System memory parameter
Mincap = '180' # System memory parameter
RestartFileFreq = '100' 

############# Calling the function #########################

FirstRun = True
copycommand = 'copy' #'cp'
os.chdir(STARTINGDIR) # Go to starting directory

for Temp in Temperatures:
    for Press in Pressures:
        os.chdir(STARTINGDIR) # Go to starting directory

        #Make directories if they don't exist
        try:
            os.mkdir(f'{Temp}')
        except FileExistsError:
            pass
        
        os.chdir(os.path.join(STARTINGDIR, Temp))
        
        try:
            os.mkdir(f'{Press}')
        except FileExistsError:
            pass

        os.chdir(os.path.join(STARTINGDIR, Temp, Press))

        # RestartList = [x for x in os.listdir() if 'Restart' in x]
        
        if FirstRun:
            runcmd('mkdir FirstRun')
            os.chdir(os.path.join(os.getcwd(), 'FirstRun')) #Enter first run directory
            CWD = os.getcwd()
            for file in os.listdir(SOURCEDIR): # Copy enabler files from source directory
                runcmd(f'{copycommand} "{os.path.join(SOURCEDIR, file)}" {os.getcwd()}')
        
            HF.MakeLAMMPSFile(CWD, Wall_V, System, EquilTime, CompTime, Wall_Z, HType,
                FeType, OType, PType, CType, ReaxFFTyping, Temp[:3], EquilTemp, Press[0],
                EquilPress, Fix_Z, Thermo_Z, Safezone, Mincap, RestartFileFreq)            
            
            HF.MakePBSFile(System, Temp, Press, CWD)

        # elif len(RestartList) == 0:
        #     runcmd('mkdir Restart_1')
        #     os.chdir(os.path.join(os.getcwd(), f'Restart_1'))
        #     for file in os.listdir(SOURCEDIR):
        #         runcmd(f'{copycommand} "{os.path.join(SOURCEDIR, file)}" {os.getcwd()}')

        #     os.chdir(os.path.join(STARTINGDIR, Temp, Press, 'FirstRun'))
        #     files = os.listdir()

        #     # Get restart file progress
        #     equilfiles = [int(x.split('.')[-1]) for x in files if 'equil.restart' in x]
        #     compfiles = [int(x.split('.')[-1]) for x in files if 'comp.restart' in x]
        #     restartfiles = equilfiles + compfiles # Concatenating restart file numbers
        #     restartfiles = sorted(restartfiles)

        #     os.chdir(os.path.join(STARTINGDIR, Temp, Press, 'Restart_1')) #Enter first run directory
        #     CWD = os.getcwd()

        #     if restartfiles[-1] <= int(EquilTime):
        #         restarttype = 'Equilibration'
        #         restartfilename = f'equil.restart.{restartfiles[-1]}'
        #         runcmd(f'{copycommand} "{os.path.join(STARTINGDIR, Temp, Press, 'FirstRun', restartfilename)}"\
        #             {os.getcwd()}')
        #         HF.MakeLAMMPSRestartFile(CWD, Wall_V, restartfilename, restarttype, System,
        #                                  EquilTime, CompTime, ReaxFFTyping, Temp, EquilTemp,
        #                                  Press, EquilPress, Fix_Z, Thermo_Z, Safezone,
        #                                  Mincap, RestartFileFreq)
        #     else:
        #         restarttype = 'CompShear'
        #         restartfilename = f'comp.restart.{restartfiles[-1]}'
        #         runcmd(f'{copycommand} "{os.path.join(STARTINGDIR, Temp, Press, 'FirstRun', restartfilename)}"\
        #             {os.getcwd()}')
        #         HF.MakeLAMMPSRestartFile(CWD, Wall_V, restartfilename, restarttype, System,
        #                                  EquilTime, CompTime, ReaxFFTyping, Temp, EquilTemp,
        #                                  Press, EquilPress, Fix_Z, Thermo_Z, Safezone,
        #                                  Mincap, RestartFileFreq)
                
        #     HF.MakePBSFile(System, Temp, Press, CWD)


        # else:
        #     RestartNumbers = [x.split('_')[-1] for x in RestartList]
        #     RestartNumber = sorted(RestartNumbers)[-1]
        #     print(RestartNumber)
        #     Restart = 'Restart_1'
        #     runcmd(f'mkdir {Restart}')
            
        #     os.chdir(os.path.join(os.getcwd(), f'{Restart}'))
        #     for file in os.listdir(SOURCEDIR):
        #         runcmd(f'{copycommand} "{os.path.join(SOURCEDIR, file)}" {os.getcwd()}')




