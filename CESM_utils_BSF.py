#################################
# The CESM python toolbox at KUP
# -- Barotropic Streamfunction --
#             Toolbox
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_utils_MOC as utils_MOC
#################################
# contained functions:
#################################
# - calc_MOC_on_model_grid()
#################################
# please log your changes below:
#################################
# 17-Mai-2016 - buerki@climate.unibe.ch : created this toolboxfew
#################################

import numpy as np
import xarray as xr

# =======================================================================================
# - BSF on model grid
# =======================================================================================
def calc_BSF_mgrd(MV, dump_MVzint=False):
    '''
    Comments:
     > Think about taking np.nansum() #!
    '''
    MVzint = xr.DataArray(MV.sum(dim='z_t'), 	# vertical integration
		name='MV vertically integrated',
		attrs={'units':u'centimeter^3/s'})

    BSF = xr.DataArray(MVzint, name='BSF on model grid', attrs={'units':u'Sv'})
    for i in np.arange(1,len(MVzint.nlon)): 	# zonal integration
      BSF[dict(nlon=i)] = BSF[dict(nlon=i)]+BSF[dict(nlon=i-1)]

    if dump_MVzint == True:
      return(BSF, MVzint)
    else:
      return(BSF)

