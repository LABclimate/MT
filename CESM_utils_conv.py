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
# - resample_equidist()
#################################
# please log your changes below:
#################################
# 30-Apr-2016 - buerki@climate.unibe.ch : created this toolbox
#                                         added add_cyclic()
# ??          - buerki@climate.unibe.ch : added project_on_auxgrd()
# 15-Jun-2016 - buerki@climate.unibe.ch : added resample_equidist() (migration from utils_dMOC)
#################################

import numpy as np
import xarray as xr
import CESM_utils_mask as utils_mask
import CESM_utils_plt as utils_plt
import CESM_utils_conv as utils_conv


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
# - resampling data on equidistant grid by interpolation (along single dimension)
# =======================================================================================
def resample_equidist(data_mgrd, mgrd, rsgrd):
    # Assumptions:
    #  > mgrd is monotonically increasing.
    # Input:
    #  > data_mgrd: data on old grid (model grid)
    #  > mgrd:      old grid (model grid)
    #  > rsgrd:     new grid (resampling grid)
    # Output:
    #  > data_rsgrd: data on new grid (resampled data)

    # Check monotony of mgrd
    if any(np.diff(mgrd)<0):
      print('WARNING: mgrd is NOT monotonically increasing.')
      return()

    # Pre-allocation of data_rsgrd | if rsgrd is longer than mgrd fill tail with nans
    data_rsgrd = np.ones(shape =rsgrd.shape*np.nan

    # Resampling
    idxm = 0                                    # index on mgrd
    idxrs = 0                                   # index on rsgrd
    while (idxrs < len(rsgrd)-1) & (rsgrd[idxrs] <= mgrd[-1]):
      while (idxm < len(mgrd)-1) & (rsgrd[idxrs] > mgrd[idxm]): # jump to closest neighbour
        idxm += 1
      if idxm == 0:                             # border values
        data_rsgrd[idxrs] = data_mgrd[idxm]
      else:                                     # centre values
        diff_1 = mgrd[idxm] - rsgrd[idxrs-1]      # > 0
        diff_2 = mgrd[idxm] - rsgrd[idxrs]    # > 0
        diff_total = diff_1 + diff_2            # = mgrd[idxm] - mgrd[idxm-1]
        # linearly weighted interpolation
        data_rsgrd[idxrs] = data_mgrd[idxm-1]*diff_2/diff_total + data_mgrd[idxm]*diff_1/diff_total
      idxrs += 1
    return(data_rsgrd)

