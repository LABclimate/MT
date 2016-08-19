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
def calc_Mxint_auxgrd(lat_ax, zd_ax, transport_type, M, ncdat, path_vars, savevar=True):
    '''
    Input:
     > lat_ax               : meridional axis of auxgrd | nparray
     > zd_ax                : vertical or density axis of auxgrd | nparray
     > transport_type       : either 'W', 'V', 'dW' or 'dV' | string
     > M                    : volume transport (MW, MV, dMW or dMV) | nparray of shape [nz, nlat, nlon]
     > ncdat                : netCDFdata to load REGION_MASK
     > path_vars            : path for saving variables | string
     > savevar              : boolean
    Output:
     > Mxint                : zonally integrated volume transport of shape [nz, nlat] | xarray
    Steps:
     > generate mask (mask_auxgrd) that is 'True' where latitudes of both grids lie in the same box
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

    # get masks and iteration-indices to speed up subsequent loops (calculate if loading from file fails)
    try:    mask_auxgrd = utils_misc.loadvar(path_vars+'mask_auxgrd')
    except: mask_auxgrd = utils_mask.gen_mask_grd_overlay_lat(lat_ax, ncdat, path_vars)
    try:    iter_maskcombo = utils_misc.loadvar(path_vars+'iter_maskcombo')
    except: iter_maskcombo = utils_mask.gen_iter_maskcombo(lat_ax, ncdat, mask_auxgrd, path_vars)

    # zonal integration along aux grid
      # ... on depth-axis
    if (transport_type == 'W') or (transport_type == 'V'):
        print('> zonal integration')
        Mxint = np.zeros([len(zd_ax), len(lat_ax)])      # pre-allocation with zeros (np-array like for speed)
        for n in iter_lat_ax:
          utils_misc.ProgBar('step', step=n, nsteps=len(iter_lat_ax))
          for j in iter_lat_M:
            for i in iter_maskcombo[n,j]:                       # limit zonal integration to Atlantic and grid-overlay
              Mxint[:,n] = np.nansum([Mxint[:,n],M[:,j,i]], axis=0)   # zonal integration
        utils_misc.ProgBar('done')

      # ... on density-axis
    elif (transport_type == 'dW') or (transport_type == 'dV'):
        print('> zonal integration')
        Mxint = np.zeros([len(zd_ax), len(lat_ax)])      # pre-allocation with zeros (np-array like for speed)
        for n in iter_lat_ax:
          utils_misc.ProgBar('step', step=n, nsteps=len(iter_lat_ax))
          for j in iter_lat_M:
            for i in iter_maskcombo[n,j]:                # limit zonal integration to Atlantic and grid-overlay
              Mxint[:,n] = np.nansum([Mxint[:,n],M[:,j,i]], axis=0)   # zonal integration
        utils_misc.ProgBar('done')
    
    # save to file and return
    if savevar == True:
        utils_misc.savevar(Mxint, path_vars+'M'+transport_type+'xint_auxgrd')  
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
