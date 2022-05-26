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


sys.path.insert(1, '/home/ricky/Documents/from_work/git/charged_molecules_binding/simulations_list/')

import AnalysisToolbox as AT

from datetime import date
today = date.today()

gc.collect()

folder_path="/home/ricky/Documents/from_work/MD/simulations/production_run/"
systems=["hMff"]


#folder_path="/home/ricky/Documents/from_work/MD/simulations/production_run/"
#systems=["testing_density_center"]





for file in os.listdir(folder_path):
    input_corr_file = folder_path+os.fsdecode(file)
    for system in systems:
        if fnmatch.fnmatch(os.fsdecode(file), "*"+system+"*"):
            newcomer=AT.AnalysisToolbox(folder_path,os.fsdecode(file),system,["ORDER_PARAMETER","BOX_DIMENSIONS"])
            #newcomer=AT.AnalysisToolbox(folder_path,os.fsdecode(file),"etidocaine",["ORDER_PARAMETER","BOX_DIMENSIONS"])
            newcomer.add_new_folders()
            #newcomer.analysis_module()
