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
# - resample_1dim_weightedmean()
#################################
# please log your changes below:
#################################
# 30-Apr-2016 - buerki@climate.unibe.ch : created this toolbox
#                                         added add_cyclic()
# ??          - buerki@climate.unibe.ch : added project_on_auxgrd()
# 15-Jun-2016 - buerki@climate.unibe.ch : added resample_1dim_weightedmean() (migration from utils_dMOC)
#################################

import numpy as np
import xarray as xr
import UTILS_misc as utils_misc
from IPython.core.debugger import Tracer; debug_here = Tracer()

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
def resample_colwise(var_mgrd, mgrd, rsgrd, method, mask='none', onlyposgrad='False'):
    ''' uses resample_1dim_lin()
    Input:
     > var_mgrd:    variable on model grid
     > mgrd:        in-situ or potential density
     > rsgrd:   new (equally spaced) bining of density
     > method:      string | 'wmean' (for vars like temp etc.) or 'sum' (for transports)
     > mask:        mask of shape [j,i], default: all True (no mask)
    Comments:
     > there's an ugly workaround for densities where the shape of mgrd is assumed to be [:,j,i]
     > add warning if negative gradients occur.
    '''
    
    print(' > columnwise resampling')
    # get default for regional mask
    if mask == 'none':
      mask = np.ones(shape=var_mgrd.shape[1:], dtype=bool)

    # pre-allocation
    var_rs = np.ones(shape=[len(rsgrd), var_mgrd.shape[1], var_mgrd.shape[2]])*np.nan
    # iteration column wise
    for j in np.arange(var_mgrd.shape[1]):
      utils_misc.ProgBar('step', step=j, nsteps=var_mgrd.shape[1])
      for i in np.arange(var_mgrd.shape[2]):
        if mask[j,i]==True: # skip masked [j,i]-tuples
          
          # workaround only for densities #! needs to be changed
          if len(mgrd.shape)==3:    mgrd_ij = mgrd[:,j,i]
          elif len(mgrd.shape)==1:  mgrd_ij = mgrd

          # eliminate all neagtive and flat density gradients
          if onlyposgrad == 'True':
            if np.any(np.diff(mgrd_ij)<=0):
              for k in np.arange(1,len(mgrd_ij)):
                if (mgrd_ij[k] <= mgrd_ij[k-1]):
                  mgrd_ij[k] = mgrd_ij[k-1]+1e-10

          if method == 'sum': #! doesn't work right now!
            var_rs[:,j,i] = resample_1dim_sum(var_mgrd[:,j,i], mgrd_ij, rsgrd)
          elif method == 'MW':
            var_rs[:,j,i] = resample_1dim_MW(var_mgrd[:,j,i], mgrd_ij, rsgrd)
          elif method == 'wmean':
            var_rs[:,j,i] = resample_1dim_weightedmean(var_mgrd[:,j,i], mgrd_ij, rsgrd)
    utils_misc.ProgBar('done')
    
    return(var_rs)

def resample_1dim_MW(data_mgrd, mgrd, rsgrd):
    '''
    Input:
     > data_mgrd: data on old grid (model grid)
     > mgrd:      old grid (model grid)
     > rsgrd:     new grid (resampling grid)
    Output:
     > data_rsgrd: data on new grid (resampled data)
    Comments:
     > #! check < and <= in while loops!!
    '''
    # Pre-allocation of data_rsgrd | if rsgrd is longer than mgrd fill tail with nans
    data_rsgrd = np.ones(shape =rsgrd.shape)*np.nan
    # reset indices
    idxm = 0                                    # index on mgrd
    idxrs = 0                                   # index on rsgrd
    if any(np.diff(mgrd)<=0): 
        print('\n\n\n\nfail!!!!!!!!!!!!\n\n\n\n')

    while (idxrs < len(rsgrd)-1) & (rsgrd[idxrs] <= np.nanmax(mgrd)): #! restricted to strictly monotonically increasing density
      # jump to closest neighbour
      while (idxm < len(mgrd)-1) & (rsgrd[idxrs] > mgrd[idxm]):
        idxm += 1
      # resampling
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

    
    
    
def resample_1dim_weightedmean(data_mgrd, mgrd, rsgrd):
    '''
    Input:
     > data_mgrd: data on old grid (model grid)
     > mgrd:      old grid (model grid)
     > rsgrd:     new grid (resampling grid)
    Output:
     > data_rsgrd: data on new grid (resampled data)
    Comments:
     > for 'short' loopingbehaviour the gradient of mgrd MUST be monotoneously increasing.
     > #! check < and <= in while loops!!
    '''
    
    # Pre-allocation of data_rsgrd | if rsgrd is longer than mgrd fill tail with nans
    data_rsgrd = np.ones(shape =rsgrd.shape)*np.nan

    # Find smart Loopingbehaviour
    if any(np.diff(mgrd)<0):     loopingbehaviour = 'long' # complete looping for mgrd with negative gradients.
    else:                        loopingbehaviour = 'short'# efficient looping only for monotoneously increasing mgrd.

    # Resampling
    idxm = 0                                    # index on mgrd
    idxrs = 0                                   # index on rsgrd
    if loopingbehaviour == 'short':
        while (idxrs <= len(rsgrd)-1) & (rsgrd[idxrs] <= np.nanmax(mgrd)): #! restricted to strictly monotonically increasing density
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
        while (idxrs <= len(rsgrd)-1):
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
    
    
    
#    def resample_1dim_sum(data_mgrd, mgrd, rsgrd):
#        '''
#        Input:
#         > data_mgrd: data on old grid (model grid)
#         > mgrd:      old grid (model grid)
#         > rsgrd:     new grid (resampling grid)
#        Output:
#         > data_rsgrd: data on new grid (resampled data)
#        Comments:
#         > not sure whether sum=0 should be overwritten with nans
#        '''
#        
#        # Pre-allocation of data_rsgrd | if rsgrd is longer than mgrd fill tail with nans
#        data_rsgrd = np.zeros(shape = len(rsgrd))
#        
#        # Resampling
#        inds = np.digitize(mgrd, rsgrd)
#        for i in np.arange(1,len(rsgrd)):
#            pass #!
#            
#    def sumup(data_mgrd, inds, i):
#        data_rsgrd[i] = np.sum(data_mgrd[np.where(inds==i)])
#        vfunc = np.vectorize(sumup)
#        vfunc(data_mgrd, inds, np.arange(1,len(rsgrd)))        
#        return(data_rsgrd)
    