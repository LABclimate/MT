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
# - ...
#################################
# please log your changes below:
#################################
# 30-Apr-2016 - buerki@climate.unibe.ch : created this toolbox
#                                         added add_cyclic()
#################################

import xarray as xr
import CESM_utils_mask as utils_mask
import CESM_utils_plt as utils_plt
import CESM_utils_conv as utils_conv

# add cyclic boundaries along nlon or nlat to prevent gap on pcolorplot
def add_cyclic(varin,dim='nlon'):
    '''Add a cyclic point to CESM data. Preserve datatype: xarray'''
    if dim == 'nlon':
        return(xr.concat([varin, varin.isel(nlon=0)], dim='nlon'))
    elif dim == 'nlat':
	return(xr.concat([varin, varin.isel(nlat=0)], dim='nlat'))

