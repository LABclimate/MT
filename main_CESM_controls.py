"""
Created on Wed Apr 20 2016

@author: buerki@climate.unibe.ch
TODO: 	 add '.values' where possible to speed up code.
"""

import numpy as np
from netCDF4 import Dataset
import xarray as xr
import pickle
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/buerki/Documents/MT/scripts/')
import CESM_utils_mask as utils_mask
import CESM_utils_plt as utils_plt
import CESM_utils_conv as utils_conv
import UTILS_misc as utils_misc
import CESM_utils_MOC as utils_MOC
import CESM_utils_BSF as utils_BSF

# #######################################################################################
#  GET AND PROCESS DATA
# #######################################################################################
# ---------------------------------------------------------------------------------------
# load netcdf file
fpath='./'
fname='b40.lm850-1850.1deg.001.pop.h.1279.ann.4.cdf'
ncdat = xr.open_dataset(fpath+fname, decode_times=False)

# =======================================================================================
#  Variables contained in model output
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Temperature
T = ncdat.TEMP.mean(dim='time').isel(z_t=0)
T = utils_mask.mask_ATLANTIC(T, ncdat.REGION_MASK)

# =======================================================================================
#  Streamfunctions
# =======================================================================================
''' BSF: Barotropic Streamfunction
    MOC: Meridional Overturning Circulation Streamfunction
    MV:  meridional volume transport
    MW:  vertical volume transport
'''
# ---------------------------------------------------------------------------------------
# - Volume transports (in Sv)
MV_mgrd = utils_BSF.calc_MV(ncdat) 					  # on model grid
MV_projauxgrd = utils_conv.project_on_auxgrd(MV_mgrd, ncdat.ANGLE.values) # on auxillary grid
MW = utils_MOC.calc_MW(ncdat)					          # valid on both grids
# ---------------------------------------------------------------------------------------
# - BSF and MOC on model grid (in Sv)
BSF_mgrd, MVzint = utils_BSF.calc_BSF_mgrd(MV_mgrd, dump_MVzint=True)
MOC_mgrd_W, MWxint_mgrd = utils_MOC.calc_MOC_mgrd('W', MW, do_normalize=True, dump_Mxint=True)
MOC_mgrd_V, MVxint_mgrd = utils_MOC.calc_MOC_mgrd('V', MV_projauxgrd, do_normalize=True, dump_Mxint=True)
# ---------------------------------------------------------------------------------------
# - Auxillary grid
lat_auxgrd, z_t_auxgrd, z_w_top_auxgrd = utils_MOC.get_default_auxgrd(ncdat)  	# define auxillary grid
lat_auxgrd = np.linspace(-80, 90, 11)
# ---------------------------------------------------------------------------------------
# - MOC on auxillary grid - WVEL (in Sv)
try:    MOC_auxgrd_W = utils_misc.loadvar('variables/MOC_auxgrd_W') 		# load from file
except: 
  try:    MWxint_auxgrd = utils_misc.loadvar('variables/MWxint_auxgrd') 	                # load from file
  except: MWxint_auxgrd = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, z_w_top_auxgrd, 'W', MW, ncdat, savevar=True)
  MOC_auxgrd_W = utils_MOC.calc_MOC_auxgrd(lat_auxgrd, z_w_top_auxgrd, 'W', MWxint_auxgrd, ncdat, savevar=True)
# ---------------------------------------------------------------------------------------
# - MOC on auxillary grid - VVEL (in Sv)
try:    MOC_auxgrd_V = utils_misc.loadvar('variables/MOC_auxgrd_V') 		# load from file
except:
  try:    MVxint_auxgrd = utils_misc.loadvar('variables/MVxint_auxgrd') 	                # load from file
  except: MVxint_auxgrd = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, z_t_auxgrd, 'V', MV_projauxgrd, ncdat, savevar=True)
  MOC_auxgrd_V = utils_MOC.calc_MOC_auxgrd(lat_auxgrd, z_t_auxgrd, 'V', MVxint_auxgrd, ncdat, savevar=True)

# ---------------------------------------------------------------------------------------
# - Zonal maximum of ocean depth
try:    HT_auxgrd_xmax = utils_misc.loadvar('variables/HT_auxgrd_xmax') 	# load from file
except: HT_auxgrd_xmax = utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'T', savevar=True)
try:    HT_mgrd_xmax = utils_misc.loadvar('variables/HT_mgrd_xmax') 	        # load from file
except: HT_mgrd_xmax = utils_mask.calc_H_mgrd_xmax(ncdat, 'T', savevar=True)
try:    HU_auxgrd_xmax = utils_misc.loadvar('variables/HU_auxgrd_xmax') 	# load from file
except: HU_auxgrd_xmax = utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'U', savevar=True)
try:    HU_mgrd_xmax = utils_misc.loadvar('variables/HU_mgrd_xmax') 	        # load from file
except: HU_mgrd_xmax = utils_mask.calc_H_mgrd_xmax(ncdat, 'U', savevar=True)

# #######################################################################################
#  PLOTTING
# #######################################################################################
plt.ion() # enable interactive mode
'''
# BSF on model grid
fig, map = utils_plt.pcolor_basemap(BSF_mgrd, 'U', nlevels=100)
plt.title('BSF mgrd on U grid')
#utils_plt.print2pdf(fig, 'testfigures/BSF_mgrd_U')

# BSF on geographical grid calculated by model
BSF_model = utils_mask.mask_ATLANTIC(ncdat.BSF.isel(time=0), ncdat.REGION_MASK)
fig, map = utils_plt.pcolor_basemap(BSF_model, 'T', cmapstep = 1)
plt.title('BSF model on T grid')
#utils_plt.print2pdf(fig, 'testfigures/BSF_model_T')
'''
'''
# MOC_mgrd_W
lat_mgrd = ncdat.TLAT.isel(nlon=0) 				# should be changed I guess
fig, ax = utils_plt.plot_slice(lat_mgrd, ncdat.z_w_top, MOC_mgrd_W, nlevels=100, plttype='contourf')
plt.plot(lat_mgrd,HT_mgrd_xmax) 				# plot seafloor #! it's the T-grid!!!
plt.title('MOC mgrd W')
plt.xlim([-36,90])
#utils_plt.print2pdf(fig, 'testfigures/MOC_mgrd_W')

# MOC_mgrd_V
lat_mgrd = ncdat.TLAT.isel(nlon=0) 				# should be changed I guess
fig, ax = utils_plt.plot_slice(lat_mgrd, ncdat.z_t, MOC_mgrd_V, nlevels=100, plttype='contourf')
plt.plot(lat_mgrd,HT_mgrd_xmax)				# plot seafloor
plt.title('MOC mgrd V')
plt.xlim([-36,90])
#utils_plt.print2pdf(fig, 'testfigures/MOC_mgrd_V')
'''
# MOC_model
MOC_model = ncdat.MOC.isel(time=0, transport_reg=1, moc_comp=0)
MOC_model = MOC_model - MOC_model[:,-1] # normalization
fig, ax = utils_plt.plot_slice(MOC_model.lat_aux_grid, MOC_model.moc_z, MOC_model, nlevels=100, plttype='contourf')
plt.plot(lat_auxgrd,HT_auxgrd_xmax) 				# plot seafloor
plt.title('MOC model')
plt.xlim([-36,90])
#utils_plt.print2pdf(fig, 'testfigures/MOC_model')
'''
# MWxint_auxgrd
fig, ax = utils_plt.plot_slice(lat_auxgrd, z_w_top_auxgrd, MWxint_auxgrd, nlevels=100, plttype='contourf')
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.plot([0]*len(MW.z_w_top), MW.z_w_top, 'x', color='red') 	# plot depth-layers
plt.xlim([-36,90])
plt.title('MWxint auxgrd')
#utils_plt.print2pdf(fig, 'testfigures/MWxint_auxgrd')

# MVxint_auxgrd
fig, ax = utils_plt.plot_slice(lat_auxgrd, z_t_auxgrd, MVxint_auxgrd, nlevels=100, plttype='contourf')
plt.plot(lat_auxgrd,HU_auxgrd_xmax)  				# plot seafloor
plt.plot([0]*len(MV_projauxgrd.z_t), MV_projauxgrd.z_t, 'x', color='red') 	# plot depth-layers
plt.xlim([-36,90])
plt.title('MVxint auxgrd')
#utils_plt.print2pdf(fig, 'testfigures/MVxint_auxgrd')
'''
# MOC_auxgrd_W
fig, ax = utils_plt.plot_slice(lat_auxgrd, z_w_top_auxgrd, MOC_auxgrd_W, nlevels=100, plttype='contourf')
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.plot([0]*len(MW.z_w_top), MW.z_w_top, 'x', color='red') 	# plot depth-layers
plt.xlim([-36,90])
plt.title('MOC auxgrd W')
#utils_plt.print2pdf(fig, 'testfigures/MOC_auxgrd_W')

# MOC_auxgrd_V
fig, ax = utils_plt.plot_slice(lat_auxgrd, z_t_auxgrd, MOC_auxgrd_V, nlevels=100, plttype='contourf')
plt.plot(lat_auxgrd,HU_auxgrd_xmax)  				# plot seafloor
plt.plot([0]*len(MV_projauxgrd.z_t), MV_projauxgrd.z_t, 'x', color='red') 	# plot depth-layers
plt.xlim([-36,90])
plt.title('MOC auxgrd V')
#utils_plt.print2pdf(fig, 'testfigures/MOC_auxgrd_V')

# ---------------------------------------------------------------------------------------
'''
# Seafloor
fig, ax = plt.contourf(ncdat.HT.roll(nlon=54), levels=np.linspace(0,560000,100))
plt.title('Depth of Seafloor')
utils_plt.print2pfig, map = utils_plt.pcolor_basemap(MW.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
fig, map = utils_plt.pcolor_basemap(MW.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
df(fig, 'testfigures/seafloor')

# Temperature
fig, map = utils_plt.pcolor_basemap(T.roll(nlon=54), 'T', nlevels=100)
plt.title('Temperature')
utils_plt.print2pdf(fig, 'testfigures/temp')

# MW
fig, map = utils_plt.pcolor_basemap(MW.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
plt.title('MW')
utils_plt.print2pdf(fig, 'testfigures/MW')

# ANGLE
plt.figure()
plt.pcolor(ncdat.ANGLE*180/np.pi)
plt.contour(ncdat.REGION_MASK)
'''



'''
fig, ax = utils_plt.plot_slice(ncdat.VVEL.TLAT[:,0], z_w_top_auxgrd, MV_projauxgrd[:,:,59])
fig, ax = utils_plt.plot_slice(ncdat.VVEL.TLAT[:,0], z_w_top_auxgrd, ncdat.VVEL[0,:,:,59])
'''
