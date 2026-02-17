#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:26:17 2024

@author: rinty
"""
import os
import subprocess
import sys

current_dir = os.getcwd()
main_dir = os.path.dirname(current_dir)
code_dir = os.path.join(main_dir, 'code')
os.chdir(code_dir)
#main_dir = f'/Users/rinty/Desktop/Research/NND/'
#code_dir = os.path.join(main_dir, 'code')
#os.chdir(code_dir)

#cat_type =  'SC_stable_rate1'#'HV_short_term_rate1'# 'original' #
data_path = 'data/'
outfolder = 'declustered_catalog'


for i in range(0,1):
    #print('running simulation %i'%i)

    j = i*1

    input_file = 'Okmok_NND_cat'
    D = 2.0 #2.0  Eruption #1.8 Okmok Historical 
 
    # List of Python scripts that require input filenames
    scripts = ['3_eta_0.py','4_NND.py', '5_dist_tau.py', '6_plot_lat_t.py',
               '7_createClust.py','8_productivity.py',] 
    
    sc1 = subprocess.run(['python3', '1_create_mat_eqCat_file.py'], input=f"{data_path}\n{input_file}\n", text=True)
    
    
    # Loop through each script and run it with the input filename
    for script in scripts:
        print(f"Running {script} with input file '{input_file}'")
        try:
            # Run each script, passing the input file using the subprocess module
            subprocess.run(['python3', script ],
                            input=f"{input_file}\n{D}\n", text=True)
    
        except Exception as e:
            print(f"An error occurred while running {script}: {e}")
    
    sc3 = subprocess.run(['python3',  '9_create_declustered_catalog.py'], input=f"{data_path}\n{input_file}\n{D}\n{outfolder}\n", text=True)
    print(f'declustering done, file saved in {outfolder}')
    
print("\a")        
    
    
    
    