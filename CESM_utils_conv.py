#################################
# The CESM python toolbox at KUP
# ----- Converting-Toolbox -----
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_utils_conv as utils_conv
#################################
# contained functions:
#################################
# - add_cyclic()
# - project_on_auxgrd()
# - resample_1dim_lin()
#################################
# please log your changes below:
#################################
# 30-Apr-2016 - buerki@climate.unibe.ch : created this toolbox
#                                         added add_cyclic()
# ??          - buerki@climate.unibe.ch : added project_on_auxgrd()
# 15-Jun-2016 - buerki@climate.unibe.ch : added resample_1dim_lin() (migration from utils_dMOC)
#################################

import numpy as np
import xarray as xr
import CESM_utils_mask as utils_mask
import CESM_utils_plt as utils_plt
import CESM_utils_conv as utils_conv
import UTILS_misc as utils_misc

# =======================================================================================
# - add cyclic boundaries along nlon or nlat to prevent gap on pcolorplot
# =======================================================================================
def add_cyclic(varin,dim='nlon'):
    '''Add a cyclic point to CESM data. Preserve datatype: xarray'''
    if dim == 'nlon':
        return(xr.concat([varin, varin.isel(nlon=0)], dim='nlon'))
    elif dim == 'nlat':
	return(xr.concat([varin, varin.isel(nlat=0)], dim='nlat'))

# =======================================================================================
# - sinusoidal projection of data on auxiliary grid
# =======================================================================================
def project_on_auxgrd(varin, angle):
    ''' angle is measured from positive x-coord on auxgrd to positive x-coord on mgrd. '''
    return(varin*np.cos(angle*np.pi/180))

# =======================================================================================
# - resampling data on new grid using linear interpolation (along single dimension)
# =======================================================================================

def resample_dens_colwise(var_mgrd, dens, dens_bins, mask='none'):
    ''' uses resample_1dim_lin() 
    Input:
     > var_mgrd:   variable on model grid (dens)
     > dens:       in-situ or potential density
     > dens_bins:  new (equally spaced) bining of density
     > mask:       mask of shape [j,i], default: all True (no mask)
    '''
    print(' > resampling on new density grid')
    # get default for regional mask
    if mask == 'none':
      mask = np.ones(shape=var_mgrd.shape[1:], dtype=bool)
    # pre-allocation
    var_rs = np.ones(shape=[len(dens_bins), var_mgrd.shape[1], var_mgrd.shape[2]])*np.nan
    # iteration column wise
    for j in np.arange(var_mgrd.shape[1]):
      utils_misc.ProgBar('step', step=j, nsteps=var_mgrd.shape[1])
      for i in np.arange(var_mgrd.shape[2]):
        if mask[j,i]==True: # skip masked [j,i]-tuples
          var_rs[:,j,i] = resample_1dim_lin(var_mgrd[:,j,i], dens[:,j,i], dens_bins)
    utils_misc.ProgBar('done')
    
    return(var_rs)
    
def resample_1dim_lin(data_mgrd, mgrd, rsgrd):
    '''
    Assumptions:
     > mgrd is monotonically increasing.
    Input:
     > data_mgrd: data on old grid (model grid)
     > mgrd:      old grid (model grid)
     > rsgrd:     new grid (resampling grid)
    Output:
     > data_rsgrd: data on new grid (resampled data)
    Comments:
     > for 'short' loopingbehaviour the gradient of mgrd MUST be monotoneously increasing.
    '''
    
    # Pre-allocation of data_rsgrd | if rsgrd is longer than mgrd fill tail with nans
    data_rsgrd = np.ones(shape =rsgrd.shape)*np.nan

    # Smart Loopingbehaviour
    if any(np.diff(mgrd)<0):     loopingbehaviour = 'long' # complete looping for mgrd with negative gradients.
    else:                        loopingbehaviour = 'short'# efficient looping only for monotoneously increasing mgrd.

    # Resampling
    idxm = 0                                    # index on mgrd
    idxrs = 0                                   # index on rsgrd
    if loopingbehaviour == 'short':
        while (idxrs < len(rsgrd)-1) & (rsgrd[idxrs] <= np.nanmax(mgrd)): #! another restriction for non-monotonical density
          while (idxm < len(mgrd)-1) & (rsgrd[idxrs] > mgrd[idxm]): # jump to closest neighbour
            idxm += 1
          if idxm == 0:                             # border values
            data_rsgrd[idxrs] = data_mgrd[idxm]
          else:                                     # centre values
            diff_1 = mgrd[idxm] - rsgrd[idxrs-1]    # > 0
            diff_2 = mgrd[idxm] - rsgrd[idxrs]      # > 0
            diff_total = diff_1 + diff_2            # = mgrd[idxm] - mgrd[idxm-1]
            # linearly weighted interpolation
            data_rsgrd[idxrs] = data_mgrd[idxm-1]*diff_2/diff_total + data_mgrd[idxm]*diff_1/diff_total
          idxrs += 1
          
    elif loopingbehaviour == 'long':
        while (idxrs < len(rsgrd)-1):
          while (idxm < len(mgrd)-1) & (rsgrd[idxrs] > mgrd[idxm]): # jump to closest neighbour
            idxm += 1
          if idxm == 0:                             # border values
            data_rsgrd[idxrs] = data_mgrd[idxm]
          else:                                     # centre values
            diff_1 = mgrd[idxm] - rsgrd[idxrs-1]    # > 0
            diff_2 = mgrd[idxm] - rsgrd[idxrs]      # > 0
            diff_total = diff_1 + diff_2            # = mgrd[idxm] - mgrd[idxm-1]
            # linearly weighted interpolation
            data_rsgrd[idxrs] = data_mgrd[idxm-1]*diff_2/diff_total + data_mgrd[idxm]*diff_1/diff_total
          idxrs += 1
          
    return(data_rsgrd)
    