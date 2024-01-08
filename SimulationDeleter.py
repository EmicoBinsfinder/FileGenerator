"""
Script to automatically delete simulations

"""

import subprocess

def runcmd(cmd, verbose = False, *args, **kwargs):
    #bascially allows python to run a bash command, and the code makes sure 
    #the error of the subproceess is communicated if it fails
    process = subprocess.run(
        cmd,
        text=True,
        shell=True)
    
runcmd('qstat > sims.txt')

sims = []
with open('sims.txt', 'r') as file:
    next(file)
    next(file)
    filelist = file.readlines()
    print(filelist)
    for x in filelist:
        sims.append(x.split(' ')[0])

for sim in sims:
    runcmd(f'qdel {sim}')