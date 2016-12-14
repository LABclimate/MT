#################################
# The CESM python toolbox at KUP
# ---- Analysis Toolbox ----
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_utils_analysis as utils_ana
#################################
# contained functions:
#################################
# - runmean()
# - normalize()
# - detrend()
# - xcorr()
#################################
# please log your changes below:
#################################
# 22-Jun-2016 - buerki@climate.unibe.ch : created this toolbox
#################################

import numpy as np
import pandas as pd
import xarray as xr
import CESM_utils_analysis as utils_ana
import CESM_utils_mask as utils_mask
import CESM_utils_conv as utils_conv
import UTILS_misc as utils_misc
#from pandas.types.missing import notnull    # for nan handling for correlations # was not available on Py2 on server, thus pd.notnull() taken instead
from scipy.special import betainc           # for pvalue of correlations (ttest)

# =======================================================================================
# - cumsum
# =======================================================================================
def nancumsum(array, axis=0):
    ''' Method: cumulative sum along axis ignoring NaNs.
                only pandas ignores NaNs in cumsum function.
        Input:
         > array : 1dimensional np-array or list
         > axis  : cumsum along this axis
    '''
    ndim = len(array.shape)
    if ndim == 1:
        return(np.array(pd.Series(array).cumsum(axis=axis)))
    elif ndim == 2:
        return(np.array(pd.DataFrame(array).cumsum(axis=axis)))
    elif ndim == 3:
        return(np.array(pd.Panel(array).cumsum(axis=axis)))
        
    
# =======================================================================================
# - canonical cumsum with span=n
# =======================================================================================
def canonical_cumsum(array, n, axis, crop=True):
    from IPython.core.debugger import Tracer; debug_here = Tracer()

    ''' Input:
         > array : np-array or list
         > n     : windowsize (n>0)
         > axis  : axis along which the cumsum will be taken
         > crop  : bool | will crop the beginning of the output array, 
                          where the sum runs over uncomplete window.
                   note: if ndim of array is >1 then the output will always be cropped.
                         --> did not find a function like np.take to /access/ an array.
    '''
    ndims = len(array.shape)
    # pandas doesn't know negative axes-allocation
    if axis<0:
        axis = ndims+axis
        
    if n<=0: 
      sys.exis('The windowsize n must be chosen greater than 0!')
    if ndims==1:
        b = utils_ana.nancumsum(array)
        b[n:] = b[n:] - b[:-n]
        if crop==True:
            b = b[n-1:]  
        return(b)
    else:
        lenax = array.shape[axis]
        b = utils_ana.nancumsum(array, axis=axis)
        b_first = np.expand_dims(np.take(b, n-1, axis=axis), axis=axis)
        c = np.take(b,np.arange(n,lenax),axis=axis) - np.take(b,np.arange(lenax-n),axis=axis)
        d = np.concatenate((b_first, c), axis=axis)
        return(d)

# =======================================================================================
# - running mean
# =======================================================================================
def runmean(datain, n, axis):
    foo = np.cumsum(np.array(datain), axis=axis, dtype=float)
    lenax = foo.shape[axis]
    return foo.take(np.arange(n,lenax), axis=axis) - foo.take(np.arange(0,lenax-n), axis=axis) / n

# =======================================================================================
# - normalisation relative to std. norm. distribution
# =======================================================================================
''' from utils on alphadata02/born'''
def normalize(datain):
    # such that series of [n x loc1 x loc2] is normalized individually for each n. (mean and std over both loc-coordinates)
    dataout = (np.subtract(datain.transpose(),np.mean(datain, axis=tuple([1,2])))/np.std(datain, axis=tuple([1,2]))).transpose()
    return dataout

# =======================================================================================
# - detrend
# =======================================================================================
''' from utils on alphadata02/born'''
def detrend(datain,keep_ave):
    datain_x = np.arange(datain.size)
    p = np.polyfit(datain_x,datain,1)
    lin_trend = np.polyval(p,datain_x)
    if(keep_ave==True):
        return datain - lin_trend + np.mean(datain)
    else:
        return datain - lin_trend

# =======================================================================================
# - cross correlation
# =======================================================================================

def nanpearsonr(a, b):
    ''' adapted from ~/anaconda3/envs/Py2/lib/python2.7/site-packages/pandas/core/nanops.py
        >> nancorr()
        calculates pearson correlation of two pd.Series a and b. 
        a and b can consist of some nans.
    '''
    
    return np.corrcoef(a, b)[0,1]

def ttest_pval(len_array, corr):  
    ''' adapted from ~/anaconda3/envs/Py2/lib/python2.7/site-packages/scipy/stats/stats.py 
        >> pearsonr()
        
        calculates the p-value of the correlation indices using a ttest.
    '''
    df = len_array - 2
    t_squared = corr**2 * (df / ((1.0 - corr) * (1.0 + corr)))
    pval = betainc(0.5*df, 0.5, df/(df+t_squared))
    return pval
    
def xpearsonr(a, b, lag=0, num_df='full', return_pval=False):
    ''' calculates timelagged correlation and pvalue of two series using
            utils_ana.nanpearsonr()* and
            utils_ana.ttest_pval()
            
            *Right now, due to speed issue pd.Sereis.corr is used.
            
        Question: check whether for pval calculation the length 
        of the series should be shortened if nans are encountered within one of the series.
        Reason: in nanpearsonr(a, b) they are also shorter.
        
        !! Important note:
            The reduction to valid entries needs
            to be performed AFTER the timelag shift.
    '''
    aa, bb = pd.Series(a), pd.Series(b)
    aash = aa.shift(-lag)
    corr = pd.Series.corr(aash, bb, 'pearson')
    #valid = notnull(aash) & notnull(bb)
    #corr = pd.Series.corr(aash[valid], bb[valid], 'pearson')
    #corr = utils_ana.nanpearsonr(aash[valid], bb[valid])
    if return_pval:
        if num_df=='full': num_df = sum(pd.notnull(aash))
        pval = utils_ana.ttest_pval(num_df, corr) # nan-values do not count for the length of aa as they are excluded inside the correlation-calculation
        return(corr, pval)
    else:
        return(corr)
    

# =======================================================================================
# - integrate 3dim data along density axis and weight with thickness of boxes
# =======================================================================================
def integrate_along_dens(dat, delta):
  # expand delta, the layer-thickness, from 1d to 3d by copying the columns
  if len(delta.shape)==1:
    delta = utils_conv.exp_k_to_kji(delta, dat.shape[-2], dat.shape[-1])
  # calculate the total thickness of each column for normalisation (only count boxes with non-nan dat value)
  delta_sum = np.nansum(delta*(np.isnan(dat)==False).astype(int), axis=0)
  delta_sum[delta_sum==0] = np.nan
  # weighted sum and normalisation with delta_sum
  dat_int = np.nansum(dat*delta, axis=0) / delta_sum
  return(dat_int)
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  