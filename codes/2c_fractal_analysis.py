"""
    - load AE lab data and compute fractal dimension
    - compute fractal dimension of arbitrary point cloud in 2D and 3D,
    compute correlation integral and fit linear portion in log-log space
"""
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

# add parent dir to path
p2parent = f"{os.environ['HOME']}/Desktop/Research/NND/" #os.path.dirname( os.getcwd())
#sys.path.append( p2parent)
code_dir = f"{os.environ['HOME']}/Desktop/Research/NND/code/" 
os.chdir(code_dir)
import fractalDim
from fractalDim import C_r_Int

#=================================0==================
#                set parameters
#====================================================
SC_file =  "Okmok_NND_cat.csv"
dPar = {
          # --- C(r) and D2----------
          'x_min'     : 0.1, #smallest distance for point density calculation and C(r)
           # sample length
           'length' : 1000, # km 
            # max. fitting range
          'xmax'  : 10, #
          #-------plotting-----------------------
          'testPlot'    : False,
          'plotFormat'  : 'png',
        }

#=================================1=========================================================
#                load data
#===========================================================================================
# a_X, a_Y, a_Z = np.loadtxt( f"{p2parent}/data/{AE_file}", usecols=(1,2,3), skiprows=1,
#                             dtype = float).T

a_X, a_Y, a_Z = np.loadtxt( f"{p2parent}/data/{SC_file}",delimiter=',',skiprows=1,
                           usecols=(9,8,10), #max_rows=2000,
                            dtype = float).T


N_AE = len( a_X)
# if dPar['testPlot']:
#     fig, ax = plt.subplots( 2,1, sharex=True)
#     ax[0].set_title( f"N-E={N_AE}")
#     ax[0].plot( a_X, a_Y, 'ko', ms = 1)
#     ax[0].set_xlabel('X (km)'), ax[0].set_ylabel( 'Y (km)')
#     #ax[0].set_xlim( 10, 40)
#     #ax[0].set_ylim( 0, 15)

#     ax[1].plot( a_X, -a_Z, 'ko', ms = 1)
#     ax[1].set_xlabel('X (km)'), ax[1].set_ylabel( 'Z (km)')
#     #ax[1].set_xlim(10, 40)
#     #ax[1].set_ylim( -10, 5)
#     plt.savefig( f"{p2parent}/plots/EvDist_{SC_file.split('.')[0]}.{dPar['plotFormat']}", dpi = 1000, 
#                 transparent = True)
#     plt.show()

#=================================2=========================================================
#                    Correlation Integral
#===========================================================================================

d_C_int = fractalDim.C_r_Int( a_X, a_Y, a_Z,
                             x_min = dPar['x_min'],
                             x_max = dPar['length'],
                            logBinning = True, binFactor=1.1)
N_max = ( N_AE*(N_AE-1))*.5
# determine max. cut-off
if 'xmax' in dPar.keys() and dPar['xmax'] is not None:
    dRmax = { 'rmax' : dPar['xmax']}
else:
    dRmax = { 'rmax' : d_C_int['aL_ave'][d_C_int['aN_sum'] > N_max][0]}
print( 'Nmax:', N_max, dRmax)
sel_pl   = d_C_int['aL_ave'] < dRmax['rmax']
#=================================3=========================================================
#                        fractal dimension
#===========================================================================================
slope, interc, r, p, __ = scipy.stats.linregress( np.log10( d_C_int['aL_ave'][sel_pl]),
                                                  np.log10( d_C_int['aCorr'][sel_pl]))
interc = 10**interc
print( 'D2' , slope, 'rmax', dRmax['rmax'])
aC_hat = interc*d_C_int['aL_ave']**slope

color  = 'k'
plt.figure(2)
ax = plt.subplot( 111)
ax.loglog( d_C_int['aL_ave'], d_C_int['aCorr'], 'o',
            ms = 5, mec = color, mfc = 'w')
ax.loglog( d_C_int['aL_ave'], aC_hat, '--', color = color,
          label = f'D = {slope:.2f}')
# mark upper bound
ax.axvline( dRmax['rmax'])
#------labels and legends
ax.legend( loc = 'upper left')
ax.set_xlabel( 'Distance (km)')
ax.set_ylabel( 'Correlation Integral')
ax.grid(True)
#plt.savefig( f"{p2parent}/plots/fixedCorrInt3D_update_{SC_file.split('.')[0]}._xmin{dPar['x_min']}_xmax{dPar['xmax']}_L{dPar['length']}km.{dPar['plotFormat']}", dpi = 1000, transparent = True)
plt.show()

import pickle
with open(f"{p2parent}/data_processed/d_CI_{SC_file[:-4]}.pkl", 'wb') as f:
    pickle.dump(d_C_int, f)

with open(f"{p2parent}/data_processed/d_CI_{SC_file[:-4]}.pkl", 'rb') as f:
    d_CI_SC = pickle.load( f)

