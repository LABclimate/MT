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
import UTILS_misc as utils_misc

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
