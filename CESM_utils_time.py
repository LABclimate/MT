#################################
# The CESM python toolbox at KUP
# -------- Time Toolbox --------
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_utils_time as utils_time
#################################
# contained functions:
#################################
# - concat()
#################################
# please log your changes below:
#################################
# 22-Jun-2016 - buerki@climate.unibe.ch : created this toolbox
#################################

import xarray as xr

# =======================================================================================
# - Concatenate xarrays along time axis
# =======================================================================================

def concat(output, input):
    try:    return(xr.concat([output, input], dim='time'))
    except: return(input)
