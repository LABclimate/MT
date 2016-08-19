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
# - resample_1dim_lininterp()
#################################
# please log your changes below:
#################################
# 30-Apr-2016 - buerki@climate.unibe.ch : created this toolbox
#                                         added add_cyclic()
# ??          - buerki@climate.unibe.ch : added project_on_auxgrd()
# 15-Jun-2016 - buerki@climate.unibe.ch : added resample_1dim_lininterp() (migration from utils_dMOC)
# 06-Jul-2016 - buerki@climate.unibe.ch : added exp_k_to_kji()
# 07-Jul-2016 - buerki@climate.unibe.ch : added exp_ji_to_kji()
#################################

import numpy as np
import xarray as xr
import UTILS_misc as utils_misc
import CESM_utils_conv as utils_conv
import sys
from IPython.core.debugger import Tracer; debug_here = Tracer()


# =======================================================================================
# - expand (copy) 1d array (in k) to 3 dimensions keeping the input array on first dimension
# =======================================================================================
def exp_k_to_kji(k_array, len_j, len_i):
    return(np.swapaxes(np.broadcast_to(k_array, (len_i, len_j, len(k_array))), 0,-1))

# =======================================================================================
# - expand (copy) 2d array (in j and i) to 3 dimensions
# =======================================================================================
def exp_ji_to_kji(ji_array, len_k):
    return(np.broadcast_to(ji_array, (len_k, ji_array.shape[0], ji_array.shape[1])))

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
# - columnwise resampling data on new grid 
# =======================================================================================
def resample_colwise(odat, ogrd, ngrd, method, fill_value=np.nan, mask='none', mono_method='True'):
    '''
    Uses:
     > utils_conv.resample_1dim_lininterp()
    Input:
     > odat:        array of dim 1 or 3 | data on model grid
     > ogrd:        array of dim | old grid
     > ngrd:        new grid 
     > method:      string | 'lin'      (linear interpolation), 
                             'dMW_db'   (for dMW calculation with density as reference for weighting)
                             'dMW_zdb'  (for dMW calculation with isopycnal depth as reference for weighting)
                             'sum'      (sum over all datapoints within bin on new grid)
     > fill_value:  float or np.nan | value used to fill in for requested points outside the range of ogrd.
     > mask:        mask of densityshape [j,i], default: all True (no mask)
     > mono_method: string | 'filter' (will sort ogrd (and odat) such that ogrd is monotonically increasing (not necess. in strict sense!))
                             'force'  (will manipulate ogrd such that strictly monotonically increasing (brute method, not recommended))
    Output:
     > ndat:    resampled data on new grid
    Comments:
     > add warning if negative gradients occur.
     > #!! for discreapencies at the low- and high-density borders see the comment in resample_1dim_weightedmean().
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

    # expand shape ngrd, ogrd and odat to 3 dimensions 
    #   note: singleton-dimensions are intended depending on original shape of odat
    def expand_shape(varin):
        if len(varin.shape) == 1:   # 1dim --> 3dim
            return(utils_conv.exp_k_to_kji(varin, len_j, len_i))
        elif len(varin.shape) == 2: # 2dim --> 3dim
            sys.exit('case of two-dimensional ogrd is not implemented yet!')
        elif len(varin.shape) == 3: # already 3dim
            return(varin) # no expansion needed.

    ngrd = expand_shape(ngrd)
    ogrd = expand_shape(ogrd)
    odat = expand_shape(odat)
    
    # get default for regional mask (i.e. do not mask anything)
    if mask == 'none':  mask = np.ones(shape=[len_j, len_i], dtype=bool)

    # pre-allocation of ndat
    ndat = fill_value * np.ones_like(ngrd)
    
    # pre-allocation of influx_highdens
    if method == 'dMW':  influx_highdens = np.zeros(shape=[len_j, len_i])
    
    # loop over columns
    for j in np.arange(len_j):
        utils_misc.ProgBar('step', step=j, nsteps=len_j)
        for i in np.arange(len_i):
            
            # skip masked [j,i]-tuples
            if mask[j,i]==False: continue
                
            # reduce ngrd, ogrd and odat to current column
            ngrd_ji = ngrd[:,j,i]
            ogrd_ji = ogrd[:,j,i]
            odat_ji = odat[:,j,i]
            
            # detect disjunct ogrd and ngrd and continue with next column
            if (np.nanmax(ogrd_ji) < np.nanmin(ngrd_ji)) or (np.nanmax(ngrd_ji) < np.nanmin(ogrd_ji)):
                print('disjunct ogrd and ngrd at (j,i)=({}, {}). (please check conservation of integrated flux!)'.format(j, i))
                continue
            
            # function to make grd strictly monotoneously increasing
            def make_mono(grd, dat, mono_method):
                if mono_method == 'sort':       # sort grd and dat relative to grd
                    idx_sort = np.argsort(grd)
                    grd, dat = grd[idx_sort], dat[idx_sort]
                elif mono_method == 'force':    # manipulate grd
                    for k in np.arange(1,ogrd.shape[0]):
                        if ogrd_ji[k] <= ogrd_ji[k-1]:
                            ogrd_ji[k] = ogrd_ji[k-1]+1e-10
                return(grd, dat)
            
            # resampling
            if method == 'lin':     # simple linear interpolation
                # make monotoneously increasing
                if any(np.diff(ogrd_ji)<=0):
                    ogrd_ji, odat_ji = make_mono(ogrd_ji, odat_ji, mono_method)
                # linear interpolation
                try:    idxn_start = np.where(ngrd_ji > np.nanmin(ogrd_ji))[0][0]
                except: idxn_start = 0  #!
                ndat[:,j,i], gaps_border, gaps_center = utils_conv.resample_1dim_lininterp(odat_ji, ogrd_ji, ngrd_ji, idxn_start, fill_value)
            
            elif method == 'dMW_db':    # for dMW calculation with density as reference for weighting
                try:    idxn_start = np.where(ngrd_ji > np.nanmin(ogrd_ji))[0][0]
                except: idxn_start = 0  #!
                ndat[:,j,i], gaps_border, gaps_center = utils_conv.resample_1dim_lininterp(odat_ji, ogrd_ji, ngrd_ji, idxn_start, fill_value=0)
                
            elif method == 'dMW_zdb':   # for dMW calculation with isopycnal depth as reference for weighting
                try:    idxn_start = np.where(ngrd_ji > np.nanmin(ogrd_ji))[0][0]
                except: idxn_start = 0  #!
                ndat[:,j,i], gaps_border, gaps_center = utils_conv.resample_1dim_lininterp(odat_ji, ogrd_ji, ngrd_ji, idxn_start, fill_value=0)
                ndat[gaps_border,j,i] = fill_value
                
            elif method == 'sum': #! doesn't work right now!
                try:    idxn_start = np.where(ngrd_ji > np.nanmin(ogrd_ji))[0][0]
                except: idxn_start = 0  #!
                ndat[:,j,i], gaps_border, gaps_center = resample_1dim_sum(odat_ji, ogrd_ji, ngrd, idxn_start, fill_value)            
                ndat[:,j,i] = fill_gaps(ndat[:,j,i], gaps_border, gaps_center, fill_value)
            
    utils_misc.ProgBar('done')

    # some statistics
    if method == 'dMW':
        print('Statistics on influx_highdens:' \
              '\n mean:   {}\n median: {}\n min:    {}\n max:    {}'.format(\
              np.nanmean(influx_highdens), np.nanmedian(influx_highdens), \
              np.nanmin(influx_highdens), np.nanmax(influx_highdens)))

    return(np.squeeze(ndat)) # remove singleton dimensions (i.e. (1d --> 3d) --> back to 1d)

# =======================================================================================
# - 1 dimensional resampling using linear interpolation
# =======================================================================================

def resample_1dim_lininterp(odat, ogrd, ngrd, idxn_start, fill_value=np.nan):
    import CESM_utils_analysis as utils_ana
    '''
    Input:
     > odat:            data on old grid
     > ogrd:            old grid
     > ngrd:            new grid (resampling grid)
     > fill_value:      0 or np.nan | used for preallocation
    Output:
     > ndat:            data on new grid (resampled data)
     > gaps_border      bool-array, same shape as ndat | True where ngrd is outside ogrd
     > gaps_center      bool-array, same shape as ndat | True where resolution of ngrd is higher than resolution of ogrd
    Comments:
     > #!! implemented, such that both, the low- and the high-value level of ngrd are left blank (i.e. on fill_value)
       --> for dMW a shift of the cumsumvalues will correct for the high-density border discrepency 
       and a comparison with the columnwise integrated transport in depth-space for those at the low-density border.
    '''

    # Pre-allocation of ndat | if ngrd is longer than ogrd fill tail with nans
    ndat = fill_value * np.ones(shape=ngrd.shape)
    gaps_border = np.zeros(shape=ngrd.shape, dtype=bool)
    gaps_center = np.zeros(shape=ngrd.shape, dtype=bool)

    # starting values and gaps at first border
    idxo = 0                                    # index on ogrd
    idxo_old = np.nan
    idxn = idxn_start                           # index on ngrd
    gaps_border[:idxn] = True                   # mark gaps at low-value border 
    #if np.any(np.isnan(ngrd)==False) and (np.where(np.isnan(ngrd)==False)[0][0] != idxn_start):
    #    debug_here()
        
    # loop through ngrd: 
        # cond 1) loop until the very last enty of ngrd.
        # cond 2) stop as soon as remaining ngrd values are all higher than the maximum value of ogrd.
    while (idxn < len(ngrd)) and (ngrd[idxn] <= np.nanmax(ogrd)):
      # lift idxo until ogrd is one step further (deeper, i.e. higher value) than ngrd.
      while (idxo < len(ogrd)-1) and (ngrd[idxn] > ogrd[idxo]):
        idxo += 1
      # resampling
      if idxo == idxo_old:                      # ngrd-resolution is higher than ogrd-resolution
        gaps_center[idxn-1:idxn+1] = True       # --> mark idxn and idxn-1 as gaps (the very careful way! #! consider defining minimum span to also mark idxn-1)
        ndat[idxn] = ndat[idxn-1]               # --> same value as above. These value remains if center-gaps are not masked.
      else:                                     # centre values
        diff_1 = -1*(ogrd[idxo-1] - ngrd[idxn]) # --> a positive number
        diff_2 = ogrd[idxo] - ngrd[idxn]        # --> a positive number
        diff_total = diff_1 + diff_2            # --> equivalent to ogrd[idxo] - ogrd[idxo-1]
        # linearly weighted interpolation
        ndat[idxn] = odat[idxo-1]*diff_2/diff_total + odat[idxo]*diff_1/diff_total
      # increase idxn and renew idxo_old
      idxo_old = idxo                         
      idxn += 1

    gaps_border[idxn:] = True                   # mark gaps at high-value border
    return(ndat, gaps_border, gaps_center)
    