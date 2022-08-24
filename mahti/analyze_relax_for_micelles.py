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



from scipy import optimize
sys.path.insert(1, '/home/ricky/Documents/from_work/git/NMR_FF_tools/relaxation_times/')
sys.path.insert(1, '/home/nenciric/Documents/git/NMR_FF_tools/relaxation_times/')

import relaxation_times as rt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
#from numba import jit #

#gyromagnetic ratios for further developmant
# !!! NOTICE!!!
#values taken from matlab code and projecct work and slightly different than those in Wikipedia
#these values are also in the external file --> if change is needed; has to be changed there
#values here in jupyter just for the information and verify, why they are different
#!!! NOTICE END !!!
gammaD=41.695*10**6; #r*s^(-1)*T^(-1)
gammaH=267.513*10**6;
gammaC=67.262*10**6;
gammaN=-27.166*10**6;


mahti_path="/scratch/project_2001058/regnierm/1_Charmmff36m/5_Micelles_50SDS_310K_CHARMM_Na_Neut_OPC_Water_Model/"
final_path="/media/nenciric/7686706b-e3c8-4ac3-a074-769f122a69d6/CoryPeptides/mahti/"

mahti_folders=["02_eElaB_micelle_40SDS_CHARMM_310K_Na_Neut_OPC",
               "01_eElaB_micelle_45SDS_CHARMM_310K_Na_Neut_OPC",
               "0_eElaB_micelle_50SDS_CHARMM_310K_Na_Neut_OPC",
               "12_eYqjD_micelle_40SDS_CHARMM_310K_Na_Neut_OPC",
               "11_eYqjD_micelle_45SDS_CHARMM_310K_Na_Neut_OPC",
               "1_eYqjD_micelle_50SDS_CHARMM_310K_Na_Neut_OPC",
               "22_hMff_micelle_40SDS_CHARMM_310K_Na_Neut_OPC",
               "21_hMff_micelle_45SDS_CHARMM_310K_Na_Neut_OPC",
               "2_hMff_micelle_50SDS_CHARMM_310K_Na_Neut_OPC",
               "23_hMff_micelle_60SDS_CHARMM_310K_Na_Neut_OPC",
               "32_yFis1_micelle_40SDS_CHARMM_310K_Na_Neut_OPC",
               "31_yFis1_micelle_45SDS_CHARMM_310K_Na_Neut_OPC",
               "3_yFis1_micelle_50SDS_CHARMM_310K_Na_Neut_OPC",
               "52_Magaining2_micelle_40SDS_CHARMM_37C_Na_Neut_OPC",
               "51_Magaining2_micelle_45SDS_CHARMM_37C_Na_Neut_OPC",
               "5_Magaining2_micelle_50SDS_CHARMM_37C_Na_Neut_OPC"]



"""Parameters to be specified by the user"""
OP=0 # order parameter
smallest_corr_time=0 # enter in log scale -3 fs; 0 ps; 3 ns; 6 us;
biggest_corr_time=5 # same as above
N_exp_to_fit=100 # number of exponential functions to be fitted between the samlles and biggest corr time
analyze=1/50 # the proportin of correlation data to be used for fitting, ex. 1/2 uses first half of the data
magnetic_field=2.35 # 5.99 # 8.49 T (values used in SDS paper, J.Chem. Soc.,, Faraday Trans. 1, 1988, 84(12), 4475-4486)
#magn_field=850
#magnetic_field=magn_field*2*np.pi/gammaH*10**6
nuclei="15N" #nuclei to calculate: 2H-deutherium; 13C - carbon; 15N - nitrogen 



##############3
## CHANGE IN THE CODE 6.4.2022, not going throught the whole content of the folder anymore
###############
take_all_in_folder="number" #"yes"/"no"/"number" analyze all in folder? useful for proteins, if no, fill the following line, if yes fill the folder path
input_corr_file="alphaCF.xvg"
input_prefix="NHrotaCF_" # mostly for peptides, works with take_all_in_folder="no"
residues=40


author_name="Ricky Nencini"



for system in mahti_folders:
    folder_mahti=mahti_path+system+"/correlationFUNCTfolderNH"
    folder_local=final_path+system
    print(folder_mahti)
    
    folder_path=folder_local+"/correlationFUNCTfolderNH/"
    
    output_name=system+".out"

    residues=0
    for file in os.listdir(folder_path):
        residues+=1

    if take_all_in_folder=="yes":
        for file in os.listdir(folder_path):
            input_corr_file = folder_path+os.fsdecode(file)
            rt.GetRelaxationData(OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze,magnetic_field,input_corr_file,nuclei,output_name)
    elif take_all_in_folder=="number":
        step_exp=(biggest_corr_time-smallest_corr_time)/N_exp_to_fit
        Ctimes = 10 ** np.arange(smallest_corr_time, biggest_corr_time, step_exp)
        Ctimes = Ctimes * 0.001 * 10 ** (-9);
        Ctimes_to_save=np.zeros([len(Ctimes),residues+1])
        Ctimes_to_save[:,0]=Ctimes
        for i in range(0,residues):
            input_corr_file = folder_path+input_prefix+str(i)+".xvg"
            AA=rt.GetRelaxationData(OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze,magnetic_field,input_corr_file,nuclei,output_name)
            Ctimes_to_save[:,i+1]=AA.Coeffs
    else:
        rt.GetRelaxationData(OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze,magnetic_field,input_corr_file,nuclei,output_name)
                                                                                                                                                        
    np.savetxt(output_name+".coeff",Ctimes_to_save)


    #try:
        #os.system("mkdir " + folder_local)
        #os.system("scp -r  nencinir@mahti.csc.fi:" + folder_mahti + " " + folder_local + "/")
    #except:
        #pass
