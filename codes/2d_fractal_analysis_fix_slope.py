#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 14 17:16:25 2025

@author: srinty
"""

import os, sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import pickle
# add parent dir to path
p2parent = os.path.dirname( os.getcwd())
sys.path.append( p2parent)
code_dir = '/home/srinty/Desktop/Research/ETAS/declustering/fractalAnalysis'
# os.chdir("..") 
# os.chdir("..") 
os.chdir(code_dir)

with open('d_CI_SC_update.pkl', 'rb') as f:
    d_C_int = pickle.load( f)
#=================================0==================
#                set parameters
#====================================================
HV_file = "SoCal_relocated_catalog_1981_2023.csv" #"HV_ANSS_full_catalog.csv"#"AE_PM_ST_152_002.txt" 
dPar = {
          # --- C(r) and D2----------
          'x_min'     : 0.1,#0.1, #smallest distance for point density calculation and C(r)
           # sample length
           'length' : 1000, # km 
            # max. fitting range
          'xmax'  : 100,#0.5, #3,
          #-------plotting-----------------------
          'testPlot'    : False,
          'plotFormat'  : 'svg',
        }



a_X, a_Y, a_Z = np.loadtxt( f"{code_dir}/data/{HV_file}",delimiter=',',skiprows=1,
                           usecols=(9,8,10), #max_rows=10000,
                            dtype = float).T
N_AE = len( a_X)
N_max = ( N_AE*(N_AE-1))*.5
# determine max. cut-off
if 'xmax' in dPar.keys() and dPar['xmax'] is not None:
    dRmax = { 'rmax' : dPar['xmax']}
    dRmin = { 'rmin' : 2}
else:
    dRmax = { 'rmax' : d_C_int['aL_ave'][d_C_int['aN_sum'] > N_max][0]}
    dRmin = { 'rmin' : dPar['x_min']}
print( 'Nmax:', N_max, dRmax)
sel_pl   = (d_C_int['aL_ave']>dRmin['rmin']) & (d_C_int['aL_ave'] < dRmax['rmax'])
#=================================3=========================================================
#                        fractal dimension
#===========================================================================================
slope, interc, r, p, __ = scipy.stats.linregress( np.log10( d_C_int['aL_ave'][sel_pl]),
                                                  np.log10( d_C_int['aCorr'][sel_pl]))

interc = 10**interc
print( 'D2' , slope, 'rmax', dRmax['rmax'])
aC_hat = interc*d_C_int['aL_ave']**slope

color  = 'k'
plt.figure( figsize = (8,6))
ax = plt.subplot( 111)
ax.loglog( d_C_int['aL_ave'], d_C_int['aCorr'], 'o',
            ms = 7, mec = color, mfc = 'w')


#=============== Inset figure ======================

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
ax_inset = inset_axes(ax, width="2000%", height="2000%", 
                      bbox_to_anchor=(0.85, .08, 0.015, 0.015),
                      bbox_transform=ax.transAxes, loc="lower right"
                      )


rows = np.arange(0,79)
n = len(rows)
D_value = np.zeros( len(rows)-1)
Rmax = np.zeros( len(rows)-1)
colors = plt.cm.rainbow(np.linspace(0,1,n))
for i, i in enumerate(rows):
    #rmin = bins[i]
    try:
        rmin = d_C_int['aL_ave'][i]
        rmax = d_C_int['aL_ave'][i+20]#bins[i+1]
    

        #=================================3=========================================================
        #                        fractal dimension
        #===========================================================================================
        slope2, interc2, r2, p2, __ = scipy.stats.linregress( np.log10( d_C_int['aL_ave'][i:i+20]),
                                                          np.log10( d_C_int['aCorr'][i:i+20]))
        
        interc2 = 10**interc2
        print( 'D2' , slope2, 'rmax', rmax, 'rmin',rmin)
        aC_hat2 = interc2*d_C_int['aL_ave'][i:i+20]**slope2
        D_value[i] = slope2
        Rmax[i] = rmax
        ax.loglog( d_C_int['aL_ave'][i:i+20], d_C_int['aCorr'][i:i+20], 'o',
                     ms = 6, mec = colors[i], mfc = 'w')
        
        ax.loglog( d_C_int['aL_ave'][i:i+20], aC_hat2, '--', color = colors[i],alpha = 0.9,)
                  #label = f' Dx = {slope2:.1f}')
    
        #ax.axvline( rmin, color = cl, alpha = 0.5, ymin = np.log10( d_C_int['aCorr'][sel_pl].min()-10), ymax = np.log10(d_C_int['aCorr'][sel_pl].min()+10))
        #ax.axvline( rmax, color = cl, alpha = 0.5, )
        ax_inset.plot(rmax,slope2,'o',color=colors[i],ms=10)

        
    except:
        pass

ax_inset.plot(Rmax,D_value,'ko', lw = 0.8, ms = 2)
ax_inset.set_xlim([-100, 1000])
ax_inset.set_ylim([0, 2.0])
ax_inset.set_xlabel( 'Rmax (km)',fontsize = 10,labelpad=0 )
ax_inset.set_ylabel( 'D value', fontsize = 10,labelpad=0)
ax_inset.tick_params(labelsize=8)
ax_inset.grid(True, which ='major', ls='--', color = 'black')
#ax_inset.grid(True, which ='minor', ls='--', color = 'gray')
ax_inset.set_facecolor('none')
minor_x = plt.MultipleLocator(5)  
ax_inset.xaxis.set_minor_locator(minor_x)
minor_y = plt.MultipleLocator(5)  
ax_inset.yaxis.set_minor_locator(minor_y)

# Set tick parameters
ax_inset.minorticks_on()
minorTick2 = {'which': 'minor', 'direction': 'out', 'length': 2, 'width': 0.5}
majorTick2 = {'which': 'major', 'direction': 'out', 'length': 4, 'width': 1}
ax_inset.tick_params(axis='both', colors='black', right=False, top=False, **minorTick2)
ax_inset.tick_params(axis='both', colors='black', right=False, top=False, **majorTick2)

#handles, labels = ax.get_legend_handles_labels()
#unique_labels = dict(zip(labels, handles))  # Use dictionary to remove duplicates



#===================main figure=================

ax.loglog( d_C_int['aL_ave'], aC_hat, '-', color = color,lw=1.5,
          label = f'D = {slope:.2f}')
# mark upper bound
ax.axvline( dRmin['rmin'], color='black')
ax.axvline( dRmax['rmax'], color = 'black')
#------labels and legends
ax.legend( loc = 'upper left')
ax.set_xlabel( 'Distance (km)', fontsize = 14)
ax.set_ylabel( 'Correlation Integral', fontsize = 14)



ax.grid()
#ax.set_ylim(10e-3, 10e1)
plt.savefig( f"{p2parent}/plots/fixed_winset5_CorrInt3D_{HV_file.split('.')[0]}_xmin{dPar['x_min']}_xmax{dPar['xmax']}_L{dPar['length']}km.{dPar['plotFormat']}", dpi = 1000, transparent = True)
plt.show()

