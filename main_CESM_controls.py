"""
Created on Wed Apr 20 2016

@author: buerki@climate.unibe.ch
"""

import numpy as np
from netCDF4 import Dataset
import xarray as xr
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/buerki/Documents/MT/scripts/')
import CESM_utils_mask as utils_mask
import CESM_utils_plt as utils_plt
import CESM_utils_conv as utils_conv
import UTILS_specials as utils_spec 	# for Progress Bars
import CESM_utils_MOC as utils_MOC
import CESM_utils_BSF as utils_BSF

# #######################################################################################
#  GET AND PROCESS DATA
# #######################################################################################
# ---------------------------------------------------------------------------------------
# load netcdf file
fpath='/alphadata02/born/lm850-1850.1deg/no_backup/annual_data/'
fname='b40.lm850-1850.1deg.001.pop.h.1279.ann.4.cdf'
ncdat = xr.open_dataset(fpath+fname, decode_times=False)

# =======================================================================================
#  Variables contained in model output
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Temperature
T = ncdat.TEMP.mean(dim='time').isel(z_t=0)
T = utils_mask.mask_ATLANTIC(temp, ncdat.REGION_MASK)

# =======================================================================================
#  Barotropic Streamfunction
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Compute meridional volume transport MV (in Sv)
MV = utils_BSF.calc_MV(ncdat)
# ---------------------------------------------------------------------------------------
# - BSF on model grid -
BSF_mgrd, MVzint = utils_BSF.calc_BSF_on_modgrd(MV, dump_MVzint=True)

# =======================================================================================
#  Meridional Overturning Streamfunction
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Compute vertical volume transport MW (in Sv)
MW = utils_MOC.calc_MW(ncdat)
# ---------------------------------------------------------------------------------------
# - MOC on model grid
MOC_mgrd = utils_MOC.calc_MOC_on_modgrd(MW)
# ---------------------------------------------------------------------------------------
# - MOC on auxillary grid
lat_auxgrd, z_w_auxgrd = utils_MOC.get_default_auxgrd(MW) # define auxillary grid
MOC_auxgrd, MWxint_auxgrd = utils_MOC.calc_MOC_on_auxgrd(lat_auxgrd, z_w_auxgrd, MW, ncdat, dump_MWxint=True, savevar=True)

# #######################################################################################
#  PLOTTING
# #######################################################################################
plt.ion() # enable interactive mode

# BSF on model grid
fig, map = utils_plt.pcolor_basemap(BSF_mgrd, 'U', cmapstep=1)
plt.title('BSF-mgrd on U grid')
utils_plt.print2pdf(fig, 'testfigures/BSF_mgrd_U')

# BSF on geographical grid calculated by model
BSF_model = utils_mask.mask_ATLANTIC(ncdat.BSF.isel(time=0), ncdat.REGION_MASK)
fig, map = utils_plt.pcolor_basemap(BSF_model, 'T', cmapstep = 1)
plt.title('BSF-model on T grid')
utils_plt.print2pdf(fig, 'testfigures/BSF_model_T')

# MOC on model grid
xcoord = ncdat.TLAT.isel(nlon=0) # should be changed I guess
fig, ax = utils_plt.plot_slice(xcoord, ncdat.z_t, MOC_mgrd, cmapstep=1, plttype='contourf')
plt.title('MOC-mgrd')
utils_plt.print2pdf(fig, 'testfigures/BSF_mgrd')

# MOC on geographical grid calculated by model
MOC_model = ncdat.MOC.isel(time=0, transport_reg=0, moc_comp=0)
fig, ax = utils_plt.plot_slice(MOC_model.lat_aux_grid, MOC_model.moc_z, MOC_model, cmapstep = 1, plttype='contourf')
plt.title('MOC-model')
utils_plt.print2pdf(fig, 'testfigures/MOC_model')

# MOC on geographical grid (manual conversion on aux grid)
fig, ax = utils_plt.plot_slice(lat_auxgrd, z_w_auxgrd, MOC_auxgrd.T, cmapstep = 1, plttype='contourf')
plt.title('MOC-auxF')
utils_plt.print2pdf(fig, 'testfigures/MOC_auxF')

# Seafloor
fig, ax = plt.contourf(ncdat.HT.roll(nlon=54), levels = np.linspace(0,560000,100))
plt.title('Depth of Seafloor')
utils_plt.print2pfig, map = utils_plt.pcolor_basemap(MW.roll(nlon=54).mean(dim='z_w_top'), 'T', cmapstep=.1)
fig, map = utils_plt.pcolor_basemap(MW.roll(nlon=54).mean(dim='z_w_top'), 'T', cmapstep=.1)
df(fig, 'testfigures/seafloor')

# Temperature
fig, map = utils_plt.pcolor_basemap(T.roll(nlon=54), 'T', cmapstep=1)
plt.title('Temperature')
utils_plt.print2pdf(fig, 'testfigures/temp')

# MW
fig, map = utils_plt.pcolor_basemap(MW.roll(nlon=54).mean(dim='z_w_top'), 'T', cmapstep=.1)
plt.title('MW')
utils_plt.print2pdf(fig, 'testfigures/MW')
