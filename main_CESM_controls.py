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
import UTILS_specials as utils_spec
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
MW = utils_MOC.calc_MW(ncdat) 						  # valid on both grids
# ---------------------------------------------------------------------------------------
# - BSF and MOC on model grid (in Sv)
BSF_mgrd, MVzint = utils_BSF.calc_BSF_mgrd(MV_mgrd, dump_MVzint=True)
MOC_mgrd_W, MWxint = utils_MOC.calc_MOC_mgrd('W', MW, do_normalize = True, dump_Mxint=True)
MOC_mgrd_V, MVxint = utils_MOC.calc_MOC_mgrd('V', MV_projauxgrd, do_normalize = True, dump_Mxint=True)
# ---------------------------------------------------------------------------------------
# - Auxillary grid
lat_auxgrd, z_t_auxgrd, z_w_top_auxgrd = utils_MOC.get_default_auxgrd(ncdat)  	# define auxillary grid
# ---------------------------------------------------------------------------------------
# - MOC on auxillary grid - WVEL (in Sv)
try:    MOC_auxgrd_W = utils_spec.loadvar('variables/MOC_auxgrd_W') 		# load from file
except: try:    MWxint = utils_spec.loadvar('variables/MOC_auxgrd_W') 	        # load from file
        except: MWxint = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, z_w_top_auxgrd, 'W', MW, ncdat, savevar=True)
        MOC_auxgrd_W = utils_MOC.calc_MOC_auxgrd(lat_auxgrd, z_w_top_auxgrd, 'W', MWxint, ncdat, savevar=True)
# ---------------------------------------------------------------------------------------
# - MOC on auxillary grid - VVEL (in Sv)
try:    MOC_auxgrd_V = utils_spec.loadvar('variables/MOC_auxgrd_V') 		# load from file
except: try:    MVxint = utils_spec.loadvar('variables/MOC_auxgrd_V') 	        # load from file
        except: MVxint = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, z_t_auxgrd, 'V', MV_projauxgrd, ncdat, savevar=True)
        MOC_auxgrd_V = utils_MOC.calc_MOC_auxgrd(lat_auxgrd, z_t_auxgrd, 'V', MVxint, ncdat, savevar=True)

# =======================================================================================
#  Seafloor masked for ATLANTIC - these are actually z_t_top values!!
# =======================================================================================
HT_auxgrd_max = utils_spec.loadvar('variables/HT_auxgrd_max') 		# load from file
HT_mgrd_max = utils_spec.loadvar('variables/HT_mgrd_max') 		# load from file

# mask HT of model grid for ATLANTIC
HTm = utils_mask.mask_ATLANTIC(ncdat.HT, ncdat.REGION_MASK)

# a few variables to speed up subsequent loops
iter_lat_auxgrd = np.arange(len(lat_auxgrd))
iter_lat_mgrd = np.arange(len(MW.nlat))

# get i-iterators for mask for auxgrd and atlantic #! rewrite comment
try: 	iter_maskcombo = utils_spec.loadvar('variables/iter_maskcombo')     
except: iter_maskcombo = utils_mask.gen_iter_maskcombo(lat_auxgrd, MW, mask_auxgrd, ncdat.REGION_MASK)

# find maximal depth for auxgrd boxes
HT_auxgrd_max = np.zeros(len(iter_lat_auxgrd))
for n in iter_lat_auxgrd:
  utils_spec.ProgBar('step', barlen=60, step=n, nsteps=len(iter_lat_auxgrd))# initialize and update progress bar
  for j in iter_lat_mgrd:
    for i in iter_maskcombo[n,j]:
      HT_auxgrd_max[n] = np.nanmax([HT_auxgrd_max[n], ncdat.HT[j,i]])
utils_spec.ProgBar('done')
utils_spec.savevar(HT_auxgrd_max, 'variables/HT_auxgrd_max') 		# save to file

# find maximal depth for mgrd boxes
HT_mgrd_max = HTm.max(dim='nlon')
utils_spec.savevar(HT_mgrd_max, 'variables/HT_mgrd_max') 		# save to file

# #######################################################################################
#  PLOTTING
# #######################################################################################
plt.ion() # enable interactive mode

# BSF on model grid
fig, map = utils_plt.pcolor_basemap(BSF_mgrd, 'U', cmapstep=1)
plt.title('BSF mgrd on U grid')
utils_plt.print2pdf(fig, 'testfigures/BSF_mgrd_U')

# BSF on geographical grid calculated by model
BSF_model = utils_mask.mask_ATLANTIC(ncdat.BSF.isel(time=0), ncdat.REGION_MASK)
fig, map = utils_plt.pcolor_basemap(BSF_model, 'T', cmapstep = 1)
plt.title('BSF model on T grid')
utils_plt.print2pdf(fig, 'testfigures/BSF_model_T')

# MOC_mgrd_W
lat_mgrd = ncdat.TLAT.isel(nlon=0) 				# should be changed I guess
fig, ax = utils_plt.plot_slice(lat_mgrd, ncdat.z_w_top, MOC_mgrd_W, cmapstep=1, plttype='contourf')
plt.plot(lat_mgrd,HTm.max(dim='nlon')) 				# plot seafloor #! it's the T-grid!!!
plt.title('MOC mgrd W')
#utils_plt.print2pdf(fig, 'testfigures/MOC_mgrd_W')

# MOC_mgrd_V
lat_mgrd = ncdat.TLAT.isel(nlon=0) 				# should be changed I guess
fig, ax = utils_plt.plot_slice(lat_mgrd, ncdat.z_t, MOC_mgrd_V, cmapstep=1, plttype='contourf')
plt.plot(lat_mgrd,HTm.max(dim='nlon')) 				# plot seafloor
plt.title('MOC mgrd V')
#utils_plt.print2pdf(fig, 'testfigures/MOC_mgrd_V')

# MOC_model
MOC_model = ncdat.MOC.isel(time=0, transport_reg=1, moc_comp=0)
#MOC_model = MOC_model - MOC_model[:,-1] # normalization
fig, ax = utils_plt.plot_slice(MOC_model.lat_aux_grid, MOC_model.moc_z, MOC_model, cmapstep = .1, plttype='contourf')
plt.plot(lat_auxgrd,HT_auxgrd_max) 				# plot seafloor
plt.title('MOC model')
plt.xlim([-36,90])
utils_plt.print2pdf(fig, 'testfigures/MOC_model')

# MOC_auxgrd_W
fig, ax = utils_plt.plot_slice(lat_auxgrd, z_w_top_auxgrd, MOC_auxgrd_W.T, cmapstep = 1, plttype='contourf')
plt.plot(lat_auxgrd,HT_auxgrd_max)  				# plot seafloor
plt.plot([-50]*len(MW.z_w_top), MW.z_w_top, 'x', color='red') 	# plot depth-layers
plt.xlim([-36,90])
plt.title('MOC auxgrd W')
utils_plt.print2pdf(fig, 'testfigures/MOC_auxgrd_W')

# MOC_auxgrd_V
fig, ax = utils_plt.plot_slice(lat_auxgrd, z_t_auxgrd, MOC_auxgrd_V.T, cmapstep = 1, plttype='contourf')
plt.plot(lat_auxgrd,HT_auxgrd_max)  				# plot seafloor
plt.plot([-50]*len(MW.z_w_top), MW.z_w_top, 'x', color='red') 	# plot depth-layers
plt.xlim([-36,90])
plt.title('MOC auxgrd V')
utils_plt.print2pdf(fig, 'testfigures/MOC_auxgrd_V')

# ---------------------------------------------------------------------------------------

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

# ANGLE
plt.figure()
plt.pcolor(ncdat.ANGLE*180/np.pi)
plt.contour(ncdat.REGION_MASK)
