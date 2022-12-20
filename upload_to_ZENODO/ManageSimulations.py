#!/usr/bin/env python

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

gc.collect()

def go_through_simulation(path,name):

    folder_path = path+name+"/"
    readmeA = folder_path+ "/README.yaml"
    
    today = str(date.today())
    
    
    if not os.path.isfile(readmeA):
        readme={}
    else:
        with open(readmeA) as yaml_file:
            readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
     
    readme=check_for_latest_files(path,name,readme)
    readme=temp_comp(folder_path,readme)
    readme=reduce_trajectories(folder_path,name,readme,1000)
    
    with open(readmeA, 'w') as f:
        yaml.dump(readme,f, sort_keys=False)
        
        
    print("great success!!!")
    
    

def reduce_trajectories(folder_path,name,sim,new_freq):
    try:
        with open(folder_path+sim["FILES"]["ndx"]["NAME"]) as f:
            for line in f:
                if "[" in line:
                    print(line.split()[1])
    
    
        skip = str(int(new_freq/int(sim["FILES"]["xtc"]["SAVING_FREQUENCY"])))
        os.system('echo System | gmx trjconv -f '+ folder_path+sim["FILES"]["xtc"]["NAME"] + ' -s ' + folder_path+sim["FILES"]["tpr"]["NAME"] + ' -skip '+ skip +' -o ' +
                      folder_path+"reduced_saving_freq_"+sim["FILES"]["xtc"]["NAME"])
                      
        sim["FILES"]["xtc_reduced_freq"]={}
        sim["FILES"]["xtc_reduced_freq"]["NAME"]="reduced_saving_freq_"+sim["FILES"]["xtc"]["NAME"]
        
        
        
 

        sim["FILES"]["xtc_reduced_freq"]['SAVING_FREQUENCY'] = new_freq
        sim["FILES"]['xtc_reduced_freq']['LENGTH'] = sim["FILES"]['xtc']['LENGTH']
        sim["FILES"]['xtc_reduced_freq']['BEGIN'] = sim["FILES"]['xtc']['BEGIN']
        
        file_adress = folder_path+sim["FILES"]["xtc_reduced_freq"]["NAME"]
        timepre=os.path.getmtime(file_adress)
        file_mod = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(timepre))
        sim["FILES"]["xtc_reduced_freq"]["SIZE"]=os.path.getsize(file_adress)/1000000
        sim["FILES"]["xtc_reduced_freq"]["MODIFIED"] = file_mod
    
    
        os.system('echo non-Water | gmx trjconv -f '+ folder_path+sim["FILES"]["xtc"]["NAME"] + ' -s ' + folder_path+sim["FILES"]["tpr"]["NAME"] + '  -o ' +
                      folder_path+"no_water_"+sim["FILES"]["xtc"]["NAME"])
        
        sim["FILES"]["xtc_no_water"]={}
        sim["FILES"]["xtc_no_water"]["NAME"]="no_water_"+sim["FILES"]["xtc"]["NAME"]
        
        
        sim["FILES"]["xtc_no_water"]['SAVING_FREQUENCY'] = sim["FILES"]['xtc']['SAVING_FREQUENCY']
        sim["FILES"]['xtc_no_water']['LENGTH'] = sim["FILES"]['xtc']['LENGTH']
        sim["FILES"]['xtc_no_water']['BEGIN'] = sim["FILES"]['xtc']['BEGIN']
        
        file_adress = folder_path+sim["FILES"]["xtc_no_water"]["NAME"]
        timepre=os.path.getmtime(file_adress)
        file_mod = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(timepre))
        sim["FILES"]["xtc_no_water"]["SIZE"]=os.path.getsize(file_adress)/1000000
        sim["FILES"]["xtc_no_water"]["MODIFIED"] = file_mod
    
    except Exception as e: 
        print(e)
    
    return sim

    
    
def temp_comp(folder_path,sim):  
        
    if not 'EQILIBRATED' in sim:
        #sim['BINDINGEQ'] = input("Biding of {} eqilibrated after [ps] \n".format(recieved_self.name))
        sim['EQILIBRATED'] = "0"

    
    log_file= folder_path+ "/"+sim["FILES"]["log"]["NAME"]    
    try:
        with open(log_file, 'rb') as log_info:
            for line in log_info:
                if ':-) GROMACS - gmx mdrun,' in str(line):
                    sim["FILES"]["log"]['GMX_VERSION']=(re.sub("'",'',str(line.split()[5])[1:]))
                      
    except:
        pass
    
    topology_tpr= folder_path+ "/"+sim["FILES"]["tpr"]["NAME"]
    try:
        with open(topology_tpr, 'rb') as tpr_info:
            for line in tpr_info:
                if 'VERSION' in str(line):
                    sim["FILES"]["tpr"]['GMX_VERSION']=(str(line.split()[1])[2:8])                    
    except:
        pass
    
    if not 'TEMPERATUREf' in sim:
        #get temperature from tpr; taken from AddData.py by Anne Kiirikki
        
        file1 =  'temporary_tpr.txt'

        print("Exporting information with gmx dump")  
        
        
        
        try:
            os.system('echo System | gmx dump -s '+ topology_tpr + ' > '+file1)

            with open(file1, 'rt') as tpr_info:
                topology_line=False
                for line in tpr_info:
                    if 'ref-t' in line:
                        sim['TEMPERATURE']=float(line.split()[1])
                    if 'topology:' in line:
                        topology_line=True
                        new_entry=False
                        
                        sim["COMPOSITION"]={}
                    if topology_line:
                        if 'moltype' in line:
                            molecule_name=re.sub('"','',line.split()[3])
                        if '#molecules' in line:
                            sim["COMPOSITION"][molecule_name]=int(line.split()[2])
                                
                    
                    if 'bIntermolecularInteractions' in line:
                        topology_line=False
            #os.system('rm '+file1)
        except:
            print("Cannot read tpr and get temperature")
  
    
    if not 'COMPOSITION' in sim:
        sim["COMPOSITION"]={}
        
        try:
            top_file= sim["FILES"]["top"]["NAME"]
            with open(top_file,"r") as f:
                molecules_list = False
                for line in f.readlines():
                    if molecules_list:
                        if not line.startswith(";"):
                            items = line.split()
                            if len(items)==2:
                                sim["COMPOSITION"][items[0]]=items[1]
                    elif line.startswith("[ molecules ]"):
                        molecules_list = True
        except:
            print("Cannot read top file and assign the composition")
            #sim["COMPOSITION"]="Encountered problems"
        
    return sim


    
    
    
    
def check_for_latest_files(path,name,sim):
    files_to_consider=["xtc","edr","tpr","top","mdp","ndx","gro","cpt","log"]
    
    if not "FILES" in sim:
        sim["FILES"]={}
      
    xtc_check=True    
    if "xtc" in sim["FILES"]:
        if not sim["FILES"]["xtc"]["NAME"]== "none":
            file_adress = path+name+"/"+sim["FILES"]["xtc"]["NAME"]
            timepre=os.path.getmtime(file_adress)
            file_mod = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(timepre))
            if file_mod==sim["FILES"]["xtc"]["MODIFIED"]:
                xtc_check=False 
        
    for fileU in files_to_consider:
        if not fileU in sim["FILES"]:
            sim["FILES"][fileU]={}
            
        folder_path=  path+name+"/"  
        for file in os.listdir(folder_path):
            if fnmatch.fnmatch(os.fsdecode(file), "*."+fileU):
                if not "NAME" in sim["FILES"][fileU] or sim["FILES"][fileU]["NAME"]== "none":
                    file_adress = path+name+"/"+os.fsdecode(file)
                    sim["FILES"][fileU]["NAME"] = os.fsdecode(file)
                else:
                    file_adress = path+name+"/"+sim["FILES"][fileU]["NAME"]
                timepre=os.path.getmtime(file_adress)
                file_mod = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(timepre))
                sim["FILES"][fileU]["SIZE"]=os.path.getsize(file_adress)/1000000
                sim["FILES"][fileU]["MODIFIED"] = file_mod
    
        try:
            print(sim["FILES"][fileU]["NAME"])
        except:
            sim["FILES"][fileU]["NAME"]= "none"
            sim["FILES"][fileU]["SIZE"]= "none"
            sim["FILES"][fileU]["MODIFIED"] = "none"
    
    
    if xtc_check:
        try:
            alku = time.time()
            mol = mda.Universe(folder_path+sim["FILES"]["gro"]["NAME"],folder_path+sim["FILES"]["xtc"]["NAME"])
            loppu = time.time()
            print("load of simulation, aikaa kului", loppu-alku, "s")
            Nframes=len(mol.trajectory)
            timestep = mol.trajectory.dt
            trj_length = Nframes * timestep
            begin_time=mol.trajectory.time

            sim["FILES"]["xtc"]['SAVING_FREQUENCY'] = timestep
            sim["FILES"]['xtc']['LENGTH'] = trj_length
            sim["FILES"]['xtc']['BEGIN'] = begin_time

        except Exception as e: 
            print(e)
            print("gro or xtc do not exist in the folder, or they do not match")
    
    
    if  sim["FILES"]["gro"]["SIZE"]== "none" and not sim["FILES"]["xtc"]["SIZE"]== "none":
        try:
            os.system('echo System | gmx trjconv -f '+ folder_path+sim["FILES"]["xtc"]["NAME"] + ' -s ' + folder_path+sim["FILES"]["tpr"]["NAME"] + ' -b 0 -e 0 -o ' +
                      folder_path+name + '.gro')
            sim["FILES"]["gro"]["NAME"] = name + '.gro'
            
            file_adress = folder_path +"/"+ sim["FILES"]["gro"]["NAME"]
                
            timepre=os.path.getmtime(file_adress)
            file_mod = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(timepre))
            sim["FILES"][fileU]["SIZE"]=os.path.getsize(file_adress)/1000000
            sim["FILES"][fileU]["MODIFIED"] = file_mod
            
            
        except:
            pass
    
    
    print("Checking for new trajectories within defiened conditions is succesfully finished")
    print(path)
    
    return sim
    
    
path="/home/nenciric/Documents/MD/pokus/"    

for file in os.listdir(path):
    go_through_simulation(path,os.fsdecode(file))
    
