"""
Date: 22 December 2023
Author: Egheosa Ogbomo

Script to automatically generate files/directories for similar conditions

- Have a source directory for files to be used in every 
- Make empty ovito_datafiles and processed folder
- Be able to toggle between generated LAMMPS file (equil or comp)
- Make QSUB files for each condition being simulated
- Make QSUB files depending on HPC being used
- Need to make sure that we are creating right stage restart for eachn condition

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
STARTINGDIR = 'F:/PhD/TCPDecompositionExperiments/Completed/DCMP_Solvent/' # Home directory to launch generation from
SOURCEDIR = os.path.join(STARTINGDIR, 'SourceDir') # Directory where enabler files are 
System = 'DCMP_Fe_48_TCP_Mixed' # System being simulated 
EquilTime = '800000' # Equilibration time
CompTime = '4000000' # Compression and shear time
Wall_V = '0.0001' # Wall velocity
Wall_Z = '7' # Wall thickness
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
EquilTemp = '300'
Safezone = '800' # System memory parameter
Mincap = '1800' # System memory parameter
RestartFileFreq = '100'
HPC = "UCL"

############# Running the script #########################

FirstRun = False
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

        # Enter condition directory
        os.chdir(os.path.join(STARTINGDIR, Temp, Press))

        if FirstRun:
            runcmd('mkdir FirstRun')
            os.chdir(os.path.join(os.getcwd(), 'FirstRun')) #Enter first run directory
            CWD = os.getcwd()
            for file in os.listdir(SOURCEDIR): # Copy enabler files from source directory
                runcmd(f'{copycommand} "{os.path.join(SOURCEDIR, file)}" {os.getcwd()}')
        
            HF.MakeLAMMPSFile(CWD, Wall_V, System, EquilTime, CompTime, Wall_Z, HType,
                FeType, OType, PType, CType, ReaxFFTyping, Temp[:3], EquilTemp, Press[0],
                EquilPress, Fix_Z, Thermo_Z, Safezone, Mincap, RestartFileFreq, HPC)            
            
            HF.MakePBSFile(System, Temp, Press, CWD)

            runcmd(f'qsub {System}_{Temp}_{Press}.pbs')

        # Check how many restarts there are 
        RestartList = [x for x in os.listdir() if 'Restart' in x]
        
        # Condition for if it's the first restart
        if len(RestartList) == 0:
            runcmd('mkdir Restart_1')
            os.chdir(os.path.join(os.getcwd(), f'Restart_1'))
            for file in os.listdir(SOURCEDIR):  
                runcmd(f'{copycommand} "{os.path.join(SOURCEDIR, file)}" {os.getcwd()}')

            FirstStage = 'FirstRun' 
            NextStage = 'Restart_1'

            #Check if first run directory has all relevant files
            if len(os.listdir(os.getcwd())) > 12:
                print('Previous simulation ran, creating files for next simulation')
                
                HF.MakeFiles(STARTINGDIR, Temp, Press, FirstStage, NextStage,
                copycommand, EquilTime, Wall_V, System, CompTime,
                ReaxFFTyping, EquilTemp, EquilPress, Fix_Z, Thermo_Z,
                Safezone, Mincap, RestartFileFreq, runcmd, HPC)
            
            else:
                print('Previous similation not yet run, creating files')
                HF.MakeFiles(STARTINGDIR, Temp, Press, FirstStage, NextStage,
                copycommand, EquilTime, Wall_V, System, CompTime,
                ReaxFFTyping, EquilTemp, EquilPress, Fix_Z, Thermo_Z,
                Safezone, Mincap, RestartFileFreq, runcmd, HPC)

        # Condition for if it's after the first restart
        else:
            # Get number of restarted simulations
            RestartNumbers = [x.split('_')[-1] for x in RestartList]
            CurrentRestartNumber = int(sorted(RestartNumbers)[-1])
            NextRestartNumber = CurrentRestartNumber + 1
            PreviousRestartNumber = CurrentRestartNumber - 1

            os.chdir(os.path.join(STARTINGDIR, Temp, Press, f'Restart_{CurrentRestartNumber}'))
            #Check if previous restart has ran
            
            if len(os.listdir(os.getcwd())) > 12:
                print('Previous simulation ran') # Could add note saying how far it ran
                os.chdir(os.path.join(STARTINGDIR, Temp, Press)) 
                runcmd(f'mkdir Restart_{NextRestartNumber}')
                os.chdir(os.path.join(STARTINGDIR, Temp, Press, f'Restart_{NextRestartNumber}')) 
                # Copy enabler files into 
                for file in os.listdir(SOURCEDIR):  
                    runcmd(f'{copycommand} "{os.path.join(SOURCEDIR, file)}" {os.getcwd()}')

                FirstStage = f'Restart_{CurrentRestartNumber}' 
                NextStage = f'Restart_{NextRestartNumber}'

                HF.MakeFiles(STARTINGDIR, Temp, Press, FirstStage, NextStage,
                copycommand, EquilTime, Wall_V, System, CompTime,
                ReaxFFTyping, EquilTemp, EquilPress, Fix_Z, Thermo_Z,
                Safezone, Mincap, RestartFileFreq, runcmd, HPC)

            else:
                print('Previous similation not yet run, creating files in current directory')
                os.chdir(os.path.join(STARTINGDIR, Temp, Press, f'Restart_{CurrentRestartNumber}')) 
                for file in os.listdir(SOURCEDIR):  
                    runcmd(f'{copycommand} "{os.path.join(SOURCEDIR, file)}" {os.getcwd()}')
                
                FirstStage = f'Restart_{PreviousRestartNumber}' 
                if PreviousRestartNumber == 0:
                    FirstStage = 'FirstRun'
                
                NextStage = f'Restart_{CurrentRestartNumber}'

                HF.MakeFiles(STARTINGDIR, Temp, Press, FirstStage, NextStage,
                copycommand, EquilTime, Wall_V, System, CompTime,
                ReaxFFTyping, EquilTemp, EquilPress, Fix_Z, Thermo_Z,
                Safezone, Mincap, RestartFileFreq, runcmd, HPC)
                