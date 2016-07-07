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
# 06-Jul-2016 - buerki@climate.unibe.ch : added expand_karray_to_kji()
#################################

import numpy as np
import xarray as xr
import UTILS_misc as utils_misc
import sys
from IPython.core.debugger import Tracer; debug_here = Tracer()


# =======================================================================================
# - expand 1d array (in k) to 3 dimensions keeping the input array on first dimension
# =======================================================================================
def expand_karray_to_kji(k_array, len_j, len_i):
    return(np.swapaxes(np.broadcast_to(k_array, (len_i, len_j, len(k_array))), 0,-1))


# =======================================================================================
# - roll np-ji-array along longitude such that Atlantic is in one piece
# =======================================================================================
def rollATL(varin):
    return(np.roll(varin, 54, axis=1))


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
    ''' angle is measured from positive x-coord on auxgrd to positive x-coord on ogrd. '''
    return(varin*np.cos(angle*np.pi/180))

# =======================================================================================
# - resampling data on new grid using linear interpolation (along single dimension)
# =======================================================================================
def resample_colwise(odat, ogrd, ngrd, method, fill_value=np.nan, mask='none', sort_ogrd='True'):
    ''' uses resample_1dim_lin()
    Input:
     > odat:        data on model grid
     > ogrd:        old grid
     > ngrd:        new grid 
     > method:      string | 'wmean' (weighted mean), 'dMW' or 'sum' (sum over all datapoints within bin on new grid)
     > fill_value:  float or nan | value used to fill in for requested points outside the range of ogrd.
     > mask:        mask of densityshape [j,i], default: all True (no mask)
     > sort_ogrd:   bool |  if True: ogrd (and odat) will be sorted such that ogrd is monotonically increasing (not necess. in strict sense!)
                            if False: ogrd will be manipulated such that strictly monotonically increasing (brute method, not recommended)
    Output:
     > ndat:    resampled data on new grid
    Comments:
     > there's an ugly workaround for densities where the shape of ogrd is assumed to be [:,j,i]
     > add warning if negative gradients occur.
    '''
    
    print(' > columnwise resampling')

    # shape of data-array
    if len(odat.shape)==3:
      len_j = odat.shape[-2] # assumed to be the second last entry
      len_i = odat.shape[-1] # assumed to be the last entry
    elif len(odat.shape)==2:
      len_j = 1              # set to one, such that j-loop is executed only once.
      len_i = odat.shape[-1] # assumed to be the last entry
    elif len(odat.shape)==1:
      len_j = 1              # set to one, such that j-loop is executed only once.
      len_i = 1              # set to one, such that i-loop is executed only once.

    # get default for regional mask (i.e. do not mask anything)
    if mask == 'none':
      mask = np.ones(shape=[len_j, len_i], dtype=bool)

    # expand the shape of ogrd and odat to 3 dimensions (singleton-dimensions are intended)
    ndim_ogrd = len(ogrd.shape)
    ndim_odat = len(odat.shape)
    if ndim_ogrd==1: # 1dim --> 3dim
      ogrd = np.swapaxes(np.broadcast_to(ogrd, (len_i,len_j,len(ogrd))), 0,-1)
    elif ndim_ogrd==2: # 2dim --> 3dim
      sys.exit('case of two-dimensional ogrd is not implemented yet!')
      
    if ndim_odat==1: # 1dim --> 3dim
      odat = np.swapaxes(np.broadcast_to(odat, (len_i,len_j,len(ogrd))), 0,-1)
    elif ndim_odat==2: # 2dim --> 3dim
      sys.exit('case of two-dimensional odat is not implemented yet!')
      
    # pre-allocation of ndat
    if method == 'wmean':
        ndat = fill_value * np.ones(shape=[len(ngrd), len_j, len_i])
    elif method == 'dMW':
        ndat = fill_value * np.ones(shape=[len(ngrd)-1, len_j, len_i])
        influx_highdens = np.zeros(shape=[len_j, len_i])
        
    # loop over columns
    for j in np.arange(len_j):
      utils_misc.ProgBar('step', step=j, nsteps=len_j)
      for i in np.arange(len_i):
        if mask[j,i]==True: # skip masked [j,i]-tuples
          # reduce ogrd and odat to current column
          ogrd_ji = ogrd[:,j,i]
          odat_ji = odat[:,j,i]
          # detect disjunct ogrd and ngrd and continue with next column
          if (np.nanmax(ogrd_ji) < np.nanmin(ngrd)) or (np.nanmax(ngrd) < np.nanmin(ogrd_ji)):
              print('disjunct ogrd and ngrd at (j,i)=({}, {}). (please check conservation of integrated flux!)'.format(j, i))
              continue
          # make ogrd strictly monotoneously increasing
          if any(np.diff(ogrd_ji)<=0):
            if sort_ogrd == 'True':
              # sort ogrd and odat
              idx_sort = np.argsort(ogrd_ji)
              ogrd_ji = ogrd_ji[idx_sort]
              odat_ji = odat_ji[idx_sort]
            else:
              # brute method by manipulating ogrd
              for k in np.arange(1,ogrd.shape[0]):
                if ogrd_ji[k] <= ogrd_ji[k-1]:
                  ogrd_ji[k] = ogrd_ji[k-1]+1e-10
          
          # interpolation
          if method == 'wmean': # simple weighted mean interpolation
            ndat[:,j,i] = resample_1dim_weightedmean(odat_ji, ogrd_ji, ngrd, fill_value)
            
          elif method == 'dMW': # procedure for dMW
            # a) weighted mean of closest neighbours around dens_bin border values
            MW_densbinborderval = resample_1dim_weightedmean(odat_ji, ogrd_ji, ngrd, fill_value)    
            # b) absolute influx (offset) from high-density.
            idxn_last = np.where(ngrd <= np.nanmax(ogrd_ji))[0][-1] # last idxn which is smaller than ogrd.max()
            idco_highdens = np.where(ogrd_ji > ngrd[idxn_last])[0]  # array with all idxo where ogrd is greater than last ngrd
            # b-1) variant where 1st idxo is taken
            #influx_highdens[j,i] = odat_ji[idco_highdens[0]]        # first idxo after idxn_last                 
            # b-2) variant where 2nd idxo is taken
            #      if this does not exist, influx_highdens is left on zero (preallocation)
            if len(idco_highdens) > 1:
              influx_highdens[j,i] = odat_ji[idco_highdens[1]]      # second idxo after idxn_last (the first is already used for interpolation of last ndat.)
            # c) get differences in MW_mgrd by substracting outflux from influx at each bin
            MWdiff_densbin = -1*np.diff(MW_densbinborderval)        # factor *-1 as MW is upward and diff goes downward

            # d) cumulative integration from dense to light water starting with influx_highdens as an offset.
            if np.isnan(fill_value):
              ndat[:,j,i] = influx_highdens[j,i] + np.cumsum(np.nan_to_num(MWdiff_densbin)[::-1])[::-1]
            else:
              ndat[:,j,i] = influx_highdens[j,i] + np.cumsum(MWdiff_densbin[::-1])[::-1]

          elif method == 'sum': #! doesn't work right now!
            ndat[:,j,i] = resample_1dim_sum(odat_ji, ogrd_ji, ngrd, fill_value)            
          
#          if (((np.sum(np.isnan(odat_ji))) == (np.sum(np.isnan(ogrd_ji))))==False):
#              debug_here()
     
    utils_misc.ProgBar('done')
    
    # some statistics
    if method == 'dMW':
        print('Statistics on influx_highdens:' \
              '\n mean:   {}\n median: {}\n min:    {}\n max:    {}'.format(\
              np.nanmean(influx_highdens), np.nanmedian(influx_highdens), \
              np.nanmin(influx_highdens), np.nanmax(influx_highdens)))

    return(np.squeeze(ndat)) # remove singleton dimensions (i.e. (1d --> 3d) --> back to 1d)


def resample_1dim_weightedmean(odat, ogrd, ngrd, fill_value=np.nan):
    '''
    Input:
     > odat:    data on old grid
     > ogrd:    old grid
     > ngrd:    new grid (resampling grid)
    Output:
     > ndat:    data on new grid (resampled data)
    Comments:
     > #! implementation of a sorting algorithm?
     > #! check < and <= in while loops!!
    '''
    
    # Pre-allocation of ndat | if ngrd is longer than ogrd fill tail with nans
    ndat = fill_value * np.ones(shape=ngrd.shape)

    # Resampling
    idxo = 0                                    # index on ogrd
    idxn = np.where(ngrd>=np.nanmin(ogrd))[0][0]   # index on ngrd
    
    # loop through ngrd | conditions: 
        # (1) loop until the very last enty of ngrd.
        # (2) stop as soon as remaining ngrd values are all higher than the maximum value of ogrd.
    while (idxn < len(ngrd)) and (ngrd[idxn] <= np.nanmax(ogrd)):
      # lift idxo until ogrd is one step further (deeper, i.e. higher value) than ngrd.
      while (idxo < len(ogrd)-1) and (ngrd[idxn] > ogrd[idxo]):
        idxo += 1
      # resampling
      if idxo == 0:                             # border values
        ndat[idxn] = odat[idxo]
      else:                                     # centre values
        diff_1 = -1*(ogrd[idxo-1] - ngrd[idxn]) # positive
        diff_2 = ogrd[idxo] - ngrd[idxn]        # positive
        diff_total = diff_1 + diff_2            # = ogrd[idxo] - ogrd[idxo-1]
        # linearly weighted interpolation
        ndat[idxn] = odat[idxo-1]*diff_2/diff_total + odat[idxo]*diff_1/diff_total
      idxn += 1
      
    return(ndat)
    
    
    
#    def resample_1dim_sum(odat, ogrd, ngrd, fill_value=np.nan):
#        '''
#        Input:
#         > odat: data on old grid (model grid)
#         > ogrd:      old grid (model grid)
#         > ngrd:     new grid (resampling grid)
#        Output:
#         > ndat: data on new grid (resampled data)
#        Comments:
#         > not sure whether sum=0 should be overwritten with nans
#        '''
#        
#        # Pre-allocation of ndat | if ngrd is longer than ogrd fill tail with nans
#        ndat = fill_value * np.ones(shape = len(ngrd))
#        
#        # Resampling
#        inds = np.digitize(ogrd, ngrd)
#        for i in np.arange(1,len(ngrd)):
#            pass #!
#            
#    def sumup(odat, inds, i):
#        ndat[i] = np.sum(odat[np.where(inds==i)])
#        vfunc = np.vectorize(sumup)
#        vfunc(odat, inds, np.arange(1,len(ngrd)))        
#        return(ndat)
    