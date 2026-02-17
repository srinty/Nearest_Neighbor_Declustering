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
#os.chdir("..") 
#os.chdir("..") 
os.chdir(code_dir)
import fractalDim
from fractalDim import C_r_Int

#=================================0==================
#                set parameters
#====================================================
SC_file = "Alaska_NND_cat.csv"#"AE_PM_ST_152_002.txt" 
dPar = {
          # --- C(r) and D2----------
          'x_min'     : 0.1, #smallest distance for point density calculation and C(r)
           # sample length
           'length' : 1000, # km 
            # max. fitting range
          'xmax'  : 10, #
          #-------plotting-----------------------
          'testPlot'    : False,
          'plotFormat'  : 'svg',
        }
plt.figure()
ax = plt.subplot( 111)
#=================================1=========================================================
#                load data
#===========================================================================================
rows = [2000, 4000,6000,10000,20000]
n = len(rows)
D_value = np.zeros( len(rows))
colors = plt.cm.rainbow(np.linspace(0,1,n))
for i, row in enumerate(rows):
    a_X, a_Y, a_Z = np.loadtxt( f"{p2parent}/data/{SC_file}",delimiter=',',skiprows=1,
                           usecols=(9,8,10), max_rows=row,
                            dtype = float).T


    N_AE = len( a_X)
    
    d_C_int = C_r_Int( a_X, a_Y, a_Z,
                                 x_min = dPar['x_min'],
                                 x_max = dPar['length'],
                                logBinning = True, binFactor=1.1)
    N_max = ( N_AE*(N_AE-1))*.5
    # determine max. cut-off
    if 'xmax' in dPar.keys() and dPar['xmax'] is not None:
        dRmin = {'rmin': 0.1}
        dRmax = { 'rmax' : 10}
    else:
        dRmax = { 'rmax' : d_C_int['aL_ave'][d_C_int['aN_sum'] > N_max][0]}
    print( 'Nmax:', N_max, dRmax)
    sel_pl = (d_C_int['aL_ave']>dRmin['rmin']) & (d_C_int['aL_ave'] < dRmax['rmax'])
    #=================================3=========================================================
    #                        fractal dimension
    #===========================================================================================
    slope, interc, r, p, __ = scipy.stats.linregress( np.log10( d_C_int['aL_ave'][sel_pl]),
                                                      np.log10( d_C_int['aCorr'][sel_pl]))
    interc = 10**interc
    print( 'D2' , slope, 'rmax', dRmax['rmax'])
    aC_hat = interc*d_C_int['aL_ave']**slope
    
    color  = 'k'

    ax.loglog( d_C_int['aL_ave'], d_C_int['aCorr'], 'o',
                ms = 5, mec = colors[i], mfc = 'w')
    ax.loglog( d_C_int['aL_ave'], aC_hat, '--', color = colors[i],
              label = f'D = {slope:.2f}, events = {row}')
    # mark upper bound
    ax.axvline( dRmin['rmin'], color = 'grey', ls = '--')
    ax.axvline( dRmax['rmax'], color = 'grey', ls = '--')
    #------labels and legends
    ax.legend( loc = 'upper left')
    ax.set_xlabel( 'Distance (km)')
    ax.set_ylabel( 'Correlation Integral')
    ax.grid(True)
    
import pickle
os.chdir(p2parent)
with open(f"{p2parent}/data_processed/d_CI_okmok.pkl", 'rb') as f:
    d_C_int = pickle.load( f)
#=================================0==================
#                set parameters
#====================================================
dPar = {
          # --- C(r) and D2----------
          'x_min'     : 0.1, #smallest distance for point density calculation and C(r)
           # sample length
           'length' : 1000, # km 
            # max. fitting range
          'xmax'  : 10, #
          #-------plotting-----------------------
          'testPlot'    : False,
          'plotFormat'  : 'svg',
        }

a_X, a_Y, a_Z = np.loadtxt( f"{p2parent}/data/{SC_file}",delimiter=',',skiprows=1,
                           usecols=(9,8,10), #max_rows=10000,
                            dtype = float).T
N_AE = len( a_X)
N_max = ( N_AE*(N_AE-1))*.5
# determine max. cut-off
if 'xmax' in dPar.keys() and dPar['xmax'] is not None:
    dRmax = { 'rmax' : 10}
    dRmin = { 'rmin' : 0.1}
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
#plt.figure(2, figsize = (8,6))
#ax = plt.subplot( 111)
ax.loglog( d_C_int['aL_ave'], d_C_int['aCorr'], 'o',
            ms = 7, mec = color, mfc = 'w')
ax.loglog( d_C_int['aL_ave'], aC_hat, '-', color = 'black',
          label = f'D = {slope:.2f}, Full catalog = {N_AE}')
# # mark upper bound
ax.axvline( dRmin['rmin'], color = 'black', ls = '--')
ax.axvline( dRmax['rmax'], color = 'black', ls = '--')
ax.legend()

plt.savefig( f"{p2parent}/plots/fixed_loop7_CorrInt3D_final_{SC_file.split('.')[0]}._xmin{dPar['x_min']}_xmax{dPar['xmax']}_L{dPar['length']}km.{dPar['plotFormat']}", dpi = 1000, transparent = True)

plt.show()

