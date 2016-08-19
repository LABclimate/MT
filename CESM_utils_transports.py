#################################
# The CESM python toolbox at KUP
# --- Meridional Overturning ---
#       Circulation Toolbox
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_utils_MOC as utils_MOC
#################################
# contained functions:
#################################
# - calc_MW()
# - calc_MOC_mgrd()
# - calc_Mxint_auxgrd()
# - calc_MOC_auxgrd()
#################################
# please log your changes below:
#################################
# 16-Jun-2016 - buerki@climate.unibe.ch : created this toolbox
#                                         migrated calc_MW() from utils_MOC
#                                         migrated calc_MV() from utils_BSF
#################################

import numpy as np
import xarray as xr
import pickle
import CESM_utils_mask as utils_mask
import UTILS_misc as utils_misc

# =======================================================================================
# - Compute vertical volume transport MW (in Sv)
# =======================================================================================
def calc_MW(ncdat):
    '''
    Comments:
     > Conversion from cgs units to Sv by multiplication with 1e-12
    '''
    wvel = utils_mask.mask_ATLANTIC(ncdat.WVEL.mean(dim='time'), ncdat.REGION_MASK)
    TAREA = ncdat.TAREA # z-area of T cells
    MW = xr.DataArray(wvel*TAREA*1e-12,
		    name='vertical volume transport',
 		    attrs={'units':u'Sv'})
    return(MW)

# =======================================================================================
# - Compute meridional volume transport MV (in Sv)
# =======================================================================================
def calc_MV(ncdat):
    '''
    Comments:
     > Conversion from cgs units to Sv by multiplication with 1e-12
     #> Think about NOT to roll #!
    '''
    vvel = utils_mask.mask_ATLANTIC(ncdat.VVEL.mean(dim='time'), ncdat.REGION_MASK)
    DXU = ncdat.DXU.values 		                    # x-spacing centered at U points
    DZU = ncdat.z_w_bot.values-ncdat.z_w_top.values  # z-spacing centered at U points
    UYAREA = np.array([DXU*ii for ii in DZU])        # y-area of U cells
    MV = xr.DataArray(vvel*UYAREA*1e-12, 
		name='meridional volume transport', 
		attrs={'units':u'Sv'})
    return(MV)
