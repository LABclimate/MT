# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 17:51:00 2016

@author: levyn
"""

import numpy as np
import gsw
from netCDF4 import Dataset
import xarray as xr
import pickle
import sys
sys.path.append('/home/buerki/Documents/MT/scripts/')
import CESM_utils_analysis as utils_ana
import CESM_utils_mask as utils_mask
import CESM_utils_conv as utils_conv
import UTILS_misc as utils_misc
import CESM_utils_transports as utils_transp
import CESM_utils_MOC as utils_MOC
import CESM_utils_BSF as utils_BSF
import CESM_utils_time as utils_time
import CESM_paths as paths
from IPython.core.debugger import Tracer; debug_here = Tracer()

from UTILS_misc import LGS GS LG

# suppress RuntimeWaring due to All-NaN slices.
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', r'All-NaN axis encountered')
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', r'All-NaN slice encountered')

# =======================================================================================
#  Load netcdf file
# =======================================================================================
fpath='../data/'
fname3='b30.004.pop.h.1000.ave.cdf'
fname4='b40.lm850-1850.1deg.001.pop.h.1279.ann.4.cdf'
ncdat3 = xr.open_dataset(fpath+fname3, decode_times=False)
ncdat4 = xr.open_dataset(fpath+fname4, decode_times=False)
CESMversion = 4
if CESMversion==3:      ncdat = ncdat3
elif CESMversion==4:    ncdat = ncdat4