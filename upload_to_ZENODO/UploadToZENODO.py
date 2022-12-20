import sys
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import os
import fnmatch
import yaml
import gc
import math
import warnings
import re
import time

from optparse import OptionParser
from collections import OrderedDict
from operator import add
from linecache import getline

from datetime import date
today = date.today()



### Change the path so it points to the folder with your folders
folder_path="/media/nenciric/Ricky20201/simulations/"

### Update the dictonary so it contains all the systems you want to upload
### True - will try to upload
### False - will not try to upload
### The true/false switch is here so you can do it per partes, or keep track of what you already succesufelly uploaded
#   It was useful when ZENODO uploads were failing all the time
 
systems={        
        "SMS_POPC_CHARMM_298K_Cl_countra_40mM_25000waters":["https://zenodo.org/deposit/7014209",True],
        "SMS_POPC_CHARMM_298K_Cl_countra_41mM_20000waters_overbind":["https://zenodo.org/deposit/7014213",True],
        "dibucaine_POPC_CHARMM_298K_Cl_countra_15mM_177600waters":["https://zenodo.org/deposit/7014252",True],
        "dibucaine_POPC_CHARMM_298K_Cl_countra_123mM":["https://zenodo.org/deposit/7014254",False]}


### No need to touch the code below
for system in systems:
    folder=folder_path+system+"/"
    if  systems[system][1]:
        for file in os.listdir(folder):
            if fnmatch.fnmatch(os.fsdecode(file), "*.xtc"):
                print(os.fsdecode(file))
                with open('upload_progress.txt', 'a') as f:
                    f.write(os.fsdecode(file))
                    f.write("\n")
                repeat_upload=True
                trial=0
                while repeat_upload:
                    gc.collect()
                    start = time.time()
                    os.system("./zenodo_upload.sh " + systems[system][0] + " " + folder+file +" 2> upload.report")
                    with open("upload.report","r") as f:
                        for line in f.readlines():
                            if "%" in line and "#" in line:
                                print_line=line
                                items = line.split()
                                if items[1]=="100.0%":
                                    repeat_upload=False
                    trial+=1
                    end = time.time()
                    with open('upload_progress.txt', 'a') as f:
                        f.write("TRIAL: {}, time: {} min \n".format(trial,(end-start)/60))
                    print("TRIAL: {}, time: {} min".format(trial,(end-start)/60))
                    try:
                        print(print_line)
                        with open('upload_progress.txt', 'a') as f:
                            f.write(print_line)
                    except:
                        print("Upload failed without ever starting")
                        with open('upload_progress.txt', 'a') as f:
                            f.write("Upload failed without ever starting")
