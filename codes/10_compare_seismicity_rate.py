#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 10:57:38 2025

@author: rinty
"""


import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#=================================1==============================================
#                            dir, file, params
#================================================================================
project_dir = os.getcwd()
main_dir = os.path.dirname(project_dir)
os.chdir(main_dir)

path = '/Users/rinty/Desktop/Research/ETAS/declustering/'

catalog_name = 'SoCal_relocated_catalog_1981_2023.csv'#'HV_ANSS_full_catalog.csv'
NN_file = catalog_name.replace('.csv', '_NN_d1.6_Mc_2.5.txt')
r85_file = catalog_name.replace('.csv', '_r85.csv')
HV_parent = pd.read_csv(f'{path}/data/data_original/{catalog_name}')
HV_parent.rename(columns={'magnitude': 'mag'}, inplace=True)


HV_NN = pd.read_csv(f'{path}/data_processed/declustered_catalog_original/{NN_file}')
HV_NN = HV_NN.sort_values(by='mag')

HV_r85 = pd.read_csv(f'{path}/r85/r85_original/{r85_file}')
HV_r85.rename(columns={'magnitude': 'mag'}, inplace=True)

#======= USGS HAZARD MODEL FILES ============
usgs_path = '/Users/rinty/Desktop/Research/ETAS/declustering/usgs_processed/'
USGS_NN = pd.read_csv(f'{usgs_path}/usgs_NN.csv')
USGS_r85 = pd.read_csv(f'{usgs_path}/usgs_r85.csv')

plt.figure(1, figsize=(8,8))
ax = plt.subplot(111)
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['xtick.labelsize']=14
plt.rcParams['ytick.labelsize']=14
plt.rcParams['axes.labelsize']=14
minorTick = {'which': 'minor', 'direction': 'out', 'length': 4, 'width': 0.5}
majorTick = {'which': 'major', 'direction': 'out', 'length': 8, 'width': 1}


df_list = [HV_parent, HV_NN, HV_r85]
df_names = ['Full Catalog', 'Nearest Neighbor', 'Reasenberg']
colors = ['black', 'b', 'r']
for i, df in enumerate(df_list):
    df_mag = np.array( sorted(df['mag'])) 
    cumul = np.cumsum( np.ones( len(df_mag)))[::-1]
    time = df['year'].max() - df['year'].min()
    ax.plot(df_mag, (cumul/time), color = colors[i], lw = 2, label='%s'%df_names[i] )

ax.set_xlabel('Magnitude')
ax.set_ylabel('Seismicity Rate [year]')
ax.legend( shadow = False, numpoints=1, loc = 'upper right')#, frameon = False)
ax.set_ylim(10e-3, 10e2)
ax.set_yscale('log')
ax.set_xlim(3, 8)
ax.grid('on')
minor_x = plt.MultipleLocator(2) 
ax.xaxis.set_minor_locator(minor_x)
minor_y = plt.MultipleLocator(10)  
ax.yaxis.set_minor_locator(minor_y)
ax.minorticks_on()
ax.tick_params(axis='both', colors='black', right=False, top=False, **minorTick)
ax.tick_params(axis='both', colors='black', right=False, top=False, **majorTick)
ax.grid(visible=True, which='major', linestyle='--', color='gray', alpha=0.5)
ax.grid(visible=True, which='minor', linestyle='--', color='gray', alpha=0.2)

ax.spines['left'].set_color('black')
ax.spines['right'].set_color('black')
plt.subplots_adjust(left  = 0.15,right = 0.95,bottom = 0.1 ,top = 0.95,wspace = 0.0,hspace = 0.0)
plt.show()
plot_dir ='/Users/rinty/Desktop/Research/ETAS/Figures/others/'
plt.savefig(f'{plot_dir}/fig_compare_rates_%s.png'%catalog_name.replace('.csv', ''), dpi = 1000, transparent = True)
