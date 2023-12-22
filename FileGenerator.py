"""
Date: 22 December 2023
Author: Egheosa Ogbomo

Script to automatically generate files/directories


- Have a source directory for files to be used in every 
- Make empty ovito_datafiles and processed folder
- Be able to toggle between generated LAMMPS file (equil or comp)
- Make QSUB files for each condition being simulated

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
STARTINGDIR = 'F:/PhD/TCPDecompositionExperiments/Completed/DCMP_Solvent/' # Home directory to launch generation from
SOURCEDIR = '' # Directory where enabler files are 
System = ''
EquilTemp = []

Temperatures = ['600K']
Pressures = ['']
Speeds = ['10ms']



if STAGE

