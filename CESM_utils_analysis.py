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
import xarray as xr
import CESM_utils_mask as utils_mask
import CESM_utils_conv as utils_conv
import UTILS_misc as utils_misc


# =======================================================================================
# - canonical cumsum with span=n
# =======================================================================================
def canonical_cumsum(array, n, crop=False):
    ''' Input:
         > array : 1dimensional np-array or list
         > n     : windowsize (n>0)
         > crop  : bool | will crop the beginning of the output array, 
                          where the sum runs over uncomplete window.
    '''
    if n<=0: 
      sys.exis('The windowsize n must be chosen greater than 0!')
      
    b = np.cumsum(array)
    b[n:] = b[n:] - b[:-n]
    if crop==True:
      b = b[n-1:]
      
    return(b)
    
    
# =======================================================================================
# - running mean
# =======================================================================================
''' from utils on alphadata02/born'''
def runmean(datain,length,crop):
  dataout = np.zeros(datain.size)
  for i in range(length,len(datain)-length):
    dataout[i] = np.mean(datain[i-length:i+length+1])
  if(crop==True):
    return dataout[length:i+1]
  else:
    return dataout


# =======================================================================================
# - normalization relative to std. norm. distribution
# =======================================================================================
''' from utils on alphadata02/born'''
def normalize(datain):
  dataout = (datain-np.mean(datain))/np.std(datain)
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
''' from utils on alphadata02/born'''
def xcorr(a,b,maxlag):
  #negative lag = time series "a" leads
  corr_coeff = np.ones(maxlag*2+1)*-999
  for i in range(-maxlag,maxlag+1):
    if(i<0):
      aa = a[:i]
      bb = b[-i:]
    elif(i>0):
      aa = a[i:]
      bb = b[:-i]
    else:
      aa = a
      bb = b
    #print i, aa.shape, bb.shape
    corr_coeff[i+maxlag] = np.corrcoef(normalize(aa),normalize(bb))[0,1]
  return corr_coeff

    
# =======================================================================================
# - integrate 3dim data along density axis and weight with thickness of boxes
# =======================================================================================
def integrate_along_dens(dat, delta):
  # expand delta, the layer-thickness, from 1d to 3d by copying the columns
  delta = utils_conv.expand_karray_to_kji(delta, dat.shape[-2],dat.shape[-1])
  # calculate the total thickness of each column for normalisation (only count boxes with non-nan dat value)
  delta_sum = np.nansum(delta*(np.isnan(dat)==False).astype(int), axis=0)
  delta_sum[delta_sum==0] = np.nan
  # weighted sum and normalisation with delta_sum
  dat_int = np.nansum(dat*delta, axis=0) / delta_sum
  return(dat_int)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  