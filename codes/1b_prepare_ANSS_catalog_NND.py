#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 11:15:14 2024

@author: rinty
"""


import numpy as np
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from scipy import stats
from scipy.stats import gaussian_kde, kstest, gamma
from scipy.stats import variation
from scipy.optimize import curve_fit
import random
import os

project_dir = os.getcwd()
main_dir = os.path.dirname(project_dir)
os.chdir(main_dir)

dir_in ='/Users/rinty/Desktop/Research/NND/data/'
os.chdir(f"{dir_in}")

files = 'Akutan_historical.csv'
Data = pd.DataFrame()
#df_list = []
#for i in range(0,len(files)):
file_in = files
valid_rows = []
datetime_rows = []
with open(file_in, 'r') as file:
    for i, line in enumerate(file, start=1):  # Enumerate for line numbers
        try:
            # Use np.loadtxt to process a single line
            row = np.loadtxt([line], delimiter=',', dtype=float, usecols=(1, 2, 3, 4))
            valid_rows.append(row)
            row2    = np.genfromtxt([line], delimiter=(4,1,2,1,2,1,2,1,2,1,4),
                                                 usecols=(0,2,4,6,8,10))
            datetime_rows.append(row2)
        except ValueError as e:
            print(f"Skipping line {i} due to error: {e}")

mData = np.array(valid_rows).T 
mDateTime = np.array(datetime_rows).T  
headDate = ['year', 'month', 'day', 'hour', 'minute','sec']
data =  {}
for i in range( len(headDate)):
    data[headDate[i]] = mDateTime[i]

#data['N']  = np.arange( len( data['year']))

header = ['latitude', 'longitude', 'depth', 'magnitude']#, 'Nst', 'Gap', 'Dmin', 'rms']

for i in range( len(header)):
    data[header[i]] = mData[i]
    
mData1 = pd.DataFrame.from_dict(data, orient='index').T  
mData1 = mData1.sort_values(by=['year', 'month', 'day', 'hour', 'minute', 'sec'])
mData1 = mData1.reset_index(inplace=False,drop=True)
Data = pd.concat([Data, mData1], ignore_index=True)

Data['N'] = np.arange( len( Data['year']))
new_order = ['year', 'month', 'day', 'hour', 'minute','sec','N','latitude', 'longitude', 'depth', 'magnitude']
Data = Data[new_order]
#D1 = Data[Data['year']<2018]
#D2 = Data[Data['year']>2017]
Data.to_csv(f"{dir_in}/Akutan_historical_catalog.csv")    
#D1.to_csv(f"{out_dir}/HV_catalog_1959_2017.csv")    
#D2.to_csv(f"{out_dir}/HV_catalog_2018_2024.csv")         
    
    