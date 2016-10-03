#################################
# The CESM python toolbox at KUP
# --- Meridional Overturning ---
#       Circulation Toolbox
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_utils_MOC as utils_MOC
#################################
# contained functions:
#################################
# - calc_MOC_mgrd()
# - calc_Mxint_auxgrd()
# - calc_MOC_auxgrd()
#################################
# please log your changes below:
#################################
# 17-Mai-2016 - buerki@climate.unibe.ch : created this toolbox
# 24-Mai-2016 - buerki@climate.unibe.ch : seperated calc_Mxint_auxgrd() from calc_MOC_auxgrd()
# 16-Jun-2016 - buerki@climate.unibe.ch : migrated calc_MW() to utils_transports
#################################

import numpy as np
import pandas as pd
import xarray as xr
import pickle
import CESM_utils_mask as utils_mask
import UTILS_misc as utils_misc
from IPython.core.debugger import Tracer; debug_here = Tracer()

# =======================================================================================
# - MOC on model grid with np-array as input
# =======================================================================================
def calc_MOC_mgrd(transport_type, M, dump_Mxint=False):
    ''' same as in calc_MOc_mgrd()'''
    # zonal integration along model grid   
    Mxint = np.nansum(M, 2) # zonal integration
    # meridional integration along model grid
    MOC = Mxint
    for j in np.arange(1,Mxint.shape[1]): 	# meridional integration
      MOC[:,j] = MOC[:,j] + MOC[:,j-1]
    # normalisation relative to North (shift values such that zero at northern boundary)
    MOC_norm = MOC - np.tile(MOC[:,-1],(MOC.shape[1],1)).T

    if dump_Mxint == True:  return(MOC, MOC_norm, Mxint)
    else:                   return(MOC, MOC_norm)


# =======================================================================================
# - MOC on auxillary grid
# =======================================================================================
'''
    This collection contains 2 functions:
     A: integration of vertical(MW)/meridional(MV) transport in zonal direction.
     B: integration of the result of B in meridional(MW)/vertical(MV) direction.

    Comments:
    =========
    > As the np.nansum() is used for zonal integration, ridges are ignored and integrated through. 
    This might be inappropriate. --> it's just the definition of the streamfunction...
    > For speed reasons some data from xarrays is copied to np-arrays
    > For speed reasons absolute indexing is used with 
    n: latitude on aux-grid, 
    j: latitude on model-grid, 
    i: longitude on model-grid and 
    k: depth on both grids
'''

# ---------------------------------------------------------------------------------------
# - zonal integration of Volume Transport along auxillary grid
# ---------------------------------------------------------------------------------------
def calc_Mxint_auxgrd(lat_ax, zd_ax, M, fraction_mask, ATLiter):
    '''
    Input:
     > lat_ax               : meridional axis of auxgrd | nparray
     > zd_ax                : vertical or density axis of auxgrd | nparray
     > M                    : volume transport (MW, MV, dMW or dMV) | nparray of shape [nz, nlat, nlon]
     > fraction_mask        : mask1TODOC
     > ATLiter              : mask2TODOC
    Output:
     > Mxint                : zonally integrated volume transport of shape [nz, nlat] | xarray
    Steps:
     > use mask (mask_auxgrd) that is 'True' where latitudes of both grids lie in the same box
     > calculate Mxint by zonal integration along aux grid
     > n-loop: over latitude on aux grid
     >   j-loop: over latitudes on model grid
     >     check whether mask_auxgrd of current n, j is True anywhere (just for speeding up)
     >     i-loop: over those longitudes where both, mask_modgrd (regional mask) and mask_auxgrd are True.
                   the order of the iteration is rolled to begin at western boundary of Atlantic
                   but actually, this 
     >       zonal integration by summing up vertical volume transports (M) of model grid for every depth (vectorially)
    '''
    # a few variables to speed up subsequent loops
    iter_lat_ax = np.arange(len(lat_ax))
    iter_lat_M = np.arange(M.shape[1])

    print('> zonal integration')
    Mxint = np.zeros([len(zd_ax), len(lat_ax)])      # pre-allocation with zeros
    for n in iter_lat_ax[:-1]: #!!
        utils_misc.ProgBar('step', step=n, nsteps=len(iter_lat_ax))
        fract_n = fraction_mask[n]
        try: foo = np.arange(fract_n[:,0].min(), fract_n[:,0].max()+1)
        except: debug_here()
        for j in foo: # limit meridional loop to fraction_maskZ
            fract_nj = fract_n[fract_n[:,0]==j]
            for i in np.intersect1d(fract_nj[:,1], ATLiter[j]): # limit zonal loop to ATLmask and fraction_mask
                fract_nji = fract_nj[fract_nj[:,1]==i,2]
                Mxint[:,n] = np.nansum([Mxint[:,n], M[:,j,i]*fract_nji], axis=0)   # zonal integration
    utils_misc.ProgBar('done')

    return(Mxint)

# ---------------------------------------------------------------------------------------
# - normalise MOC
# ---------------------------------------------------------------------------------------
def normalise(MOC, ref):
    '''
    Input:
     > MOC              : MOC | nparray
     > ref              : reference | string ('N', 'S')
    '''
    # normalisation relative to North (shift values such that zero at northern boundary)
    if ref=='N': 
        return(MOC - np.tile(MOC[:,-1],(MOC.shape[1],1)).T)
    else:
        return()
