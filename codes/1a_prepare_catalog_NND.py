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

dir_in = '/Users/rinty/Desktop/Research/NND/data/'
os.chdir(f"{dir_in}")
#files = [f for f in os.listdir() if f.endswith('.txt')]
#files.sort()
#files = 'OK08_GrowClust.csv'
files = 'yellowstone_swarm.csv'
# #============= prepare catalog ================
  
df = pd.read_csv(files)
df.replace('', 1.0, inplace=True)       # Replace empty strings
df.fillna(1.0, inplace=True) 
# Split Date into YR, MO, DY
date_parts = pd.to_datetime(df['otime']).dt
df['year'] = date_parts.year
df['month'] = date_parts.month
df['day'] = date_parts.day

# Split Time into HR, MN, SC (seconds as integer)
time_parts = pd.to_datetime(df['otime']).dt
df['hour'] = time_parts.hour
df['minute'] = time_parts.minute
df['sec'] = time_parts.second
df['N'] = np.arange( len( df))

# # Rename and reorder columns
# df_processed = df[['YR', 'MO', 'DY', 'HR', 'MN', 'SC', 'N','Latitude', 'Longitude', 'Depth[km]', 'Magnitude(ML)']]

df = df.rename(columns={
    'lat': 'latitude',
    'lon': 'longitude',
    'mag': 'magnitude',
    'dep': 'depth'
})

# Rename and reorder columns
df_processed = df[['year', 'month', 'day', 'hour', 'minute','sec','N','latitude', 'longitude', 'depth', 'magnitude']]



# Save to new CSV file
df_processed.to_csv(f'{dir_in}/yellowstone_cat.csv' )#index=False)
    