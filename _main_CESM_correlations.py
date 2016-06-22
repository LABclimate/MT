"""
Created on Wed Apr 20 2016

@author: buerki@climate.unibe.ch
TODO: 	 add '.values' where possible to speed up code.
"""

import numpy as np
import gsw
from netCDF4 import Dataset
import xarray as xr
import pickle
import matplotlib.pyplot as plt
import time
import sys
sys.path.append('/home/buerki/Documents/MT/scripts/')
import CESM_utils_mask as utils_mask
import CESM_utils_plt as utils_plt
import CESM_utils_conv as utils_conv
import UTILS_misc as utils_misc
import CESM_utils_transports as utils_transp
import CESM_utils_MOC as utils_MOC
import CESM_utils_dMOC as utils_dMOC
import CESM_utils_BSF as utils_BSF
import CESM_paths as paths

# ---------------------------------------------------------------------------------------
# paths
path_mgrd = 'vars_mgrd/'

# #######################################################################################
#  GET AND PROCESS DATA
# #######################################################################################
# ---------------------------------------------------------------------------------------
# load netcdf file
fpath='./'
fname='b40.lm850-1850.1deg.001.pop.h.1279.ann.4.cdf'
ncdat = xr.open_dataset(fpath+fname, decode_times=False)

# ---------------------------------------------------------------------------------------
# - Volume transports (in Sv)
MV_mgrd = utils_transp.calc_MV(ncdat)                                           # on model grid
MV_projauxgrd = utils_conv.project_on_auxgrd(MV_mgrd, ncdat.ANGLE.values)       # on auxiliary grid
MW = utils_transp.calc_MW(ncdat)                                                # valid on both grids
# ---------------------------------------------------------------------------------------
# - Atlantic Streamfunctions (in Sv)
BSF_mod = utils_mask.mask_ATLANTIC(ncdat.BSF, ncdat.REGION_MASK)
MOC_mod = ncdat.MOC.isel(transport_reg=1, moc_comp=0)
MOC_mod = MOC_mod - MOC_mod[:,-1] # normalisation

BSF_mgrd, MVzint = utils_BSF.calc_BSF_mgrd(MV_mgrd, dump_MVzint=True)
MOC_mgrd_W, MWxint_mgrd = utils_MOC.calc_MOC_mgrd('W', MW, do_norm=True, dump_Mxint=True)

# =======================================================================================
#  Zonal maxima of ocean depth
# =======================================================================================
try:    HT_mgrd_xmax = utils_misc.loadvar(path_mgrd+'HT_mgrd_xmax')             # load from file
except: HT_mgrd_xmax = utils_mask.calc_H_mgrd_xmax(ncdat, 'T', path_mgrd)
try:    HU_mgrd_xmax = utils_misc.loadvar(path_mgrd+'HU_mgrd_xmax')             # load from file
except: HU_mgrd_xmax = utils_mask.calc_H_mgrd_xmax(ncdat, 'U', path_mgrd)

# #######################################################################################
#  ANALYSIS
# #######################################################################################
''' check: already cut to Atlantic?'''
MOC = MOC_mod
BSF = BSF_mod
# ---------------------------------------------------------------------------------------
# time windowing

# ---------------------------------------------------------------------------------------
# select region for Streamfunction indices
MOC = MOC
BSF = BSF.where((BSF.TLAT>=45) & (BSF.TLAT<=70))
TAREA = ncdat.TAREA.where((BSF.TLAT>=45) & (BSF.TLAT<=70))
# ---------------------------------------------------------------------------------------
# calculate Streamfunction indices
MOCidx = MOC.max(dim = ['moc_z','lat_aux_grid'])
BSFidx = (BSF*TAREA).mean()/TAREA.mean()
# ---------------------------------------------------------------------------------------
# calculate correlation btw. MOCidx and BSFidx
corrIDX = utils_ana.xcorr(MOCidx, BSFidx, 0)
# calculate pointwise correlation of MOC with BSFidx
corrBSFidx = np.zeros_like(MOC)
for j in np.arange(MOC.shape[-2]):
    for i in np.arange(MOC.shape[-1]):
        corrBSFidx[j,i] = utils_ana.xcorr(MOC[:,j,i], BSFidx, 0)
# calculate pointwise correlation of BSF with MOCidx
corrMOCidx = np.zeros_like(BSF)
for j in np.arange(BSF.shape[-2]):
    for i in np.arange(BSF.shape[-1]):
        corrMOCidx = utils_ana.xcorr(BSF[:,j,i], MOCidx, 0)


# #######################################################################################
#  PLOTTING
# #######################################################################################
plt.ion() # enable interactive mode
path_fig = 'corrfigs_figures_Jun22/'
# =======================================================================================
#  CCSM4 representations
# =======================================================================================
# -----------------------------------------------------------------------------------------
# BSF on geographical grid calculated by model
fig, map = utils_plt.plot_BSF(BSF_model, 'T', nlevels = 10)
plt.title('BSF model on T grid')
 #utils_plt.print2pdf(fig, 'testfigures/BSF_model_T')
# -----------------------------------------------------------------------------------------
# MOC on geographical grid calculated by model
fig, ax = utils_plt.plot_MOC(MOC_model.lat_aux_grid, MOC_model.moc_z, MOC_model, nlevels=40, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HT_auxgrd_xmax) 				# plot seafloor
plt.title('MOC model')
plt.xlim([-36,90])
 #utils_plt.print2pdf(fig, 'testfigures/MOC_model')

# =======================================================================================
#  Calculated on auxiliary (geographical) grid
# =======================================================================================
# -----------------------------------------------------------------------------------------
# MOC_auxgrd_W
fig, ax = utils_plt.plot_MOC(lat_auxgrd, z_w_top_auxgrd, MOC_auxgrd_W_norm, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
plt.title('MOC auxgrd W')
utils_plt.print2pdf(fig, path_fig+'MOC_auxgrd_W')
# -----------------------------------------------------------------------------------------
# MOC_auxgrd_V
fig, ax = utils_plt.plot_MOC(lat_auxgrd, zT_auxgrd, MOC_auxgrd_V_norm, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HU_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
plt.title('MOC auxgrd V')
 #utils_plt.print2pdf(fig, 'testfigures/MOC_auxgrd_V')
# -----------------------------------------------------------------------------------------
# dMOC_auxgrd_W (in Sv)
fig, ax = utils_plt.plot_MOC(lat_auxgrd, PD_bins, dMOC_auxgrd_W_norm, nlevels=10, plttype='pcolor+contour')
plt.title('dMOC auxgrd W (sigma2)')
plt.suptitle('density binning from {} to {} in {} steps'.format(PD_bins.min(), PD_bins.max(), len(PD_bins)))
plt.xlim([-36,90])
utils_plt.print2pdf(fig, path_fig+'dMOC_auxgrd_W_sig2')
# -----------------------------------------------------------------------------------------
# MWxint_auxgrd
fig, ax = utils_plt.plot_MOC(lat_auxgrd, z_w_top_auxgrd, MWxint_auxgrd, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
plt.title('MWxint auxgrd')
 #utils_plt.print2pdf(fig, 'testfigures/MWxint_auxgrd')
# -----------------------------------------------------------------------------------------
# MVxint_auxgrd
fig, ax = utils_plt.plot_MOC(lat_auxgrd, zT_auxgrd, MVxint_auxgrd, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HU_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
plt.title('MVxint auxgrd')
 #utils_plt.print2pdf(fig, 'testfigures/MVxint_auxgrd')














# =======================================================================================
#  Other variables
# =======================================================================================
# -----------------------------------------------------------------------------------------
# Seafloor
fig, ax = plt.contourf(ncdat.HT.roll(nlon=54), levels=np.linspace(0,560000,100))
plt.title('Depth of Seafloor')
fig, map = utils_plt.pcolor_basemap(MW.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
fig, map = utils_plt.pcolor_basemap(MW.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
 #utils_plt.print2pdf(fig, 'testfigures/seafloor')
# -----------------------------------------------------------------------------------------
# Temperature
fig, map = utils_plt.pcolor_basemap(T.roll(nlon=54), 'T', nlevels=100)
plt.title('Temperature')
 #utils_plt.print2pdf(fig, 'testfigures/temp')
# -----------------------------------------------------------------------------------------
# Density
fig, ax = utils_plt.plot_MOC(ncdat.TLAT[:,0], ncdat.z_t, np.nanmean(sig2, 2), nlevels=10, plttype='contour')
plt.title('Density')
 #utils_plt.print2pdf(fig, 'density/temp')
# -----------------------------------------------------------------------------------------
# MW
fig, map = utils_plt.pcolor_basemap(MW.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
plt.title('MW')
 #utils_plt.print2pdf(fig, 'testfigures/MW')
# -----------------------------------------------------------------------------------------
# ANGLE
plt.figure()
plt.pcolor(ncdat.ANGLE*180/np.pi)
plt.contour(ncdat.REGION_MASK)

fig, ax = utils_plt.plot_MOC(ncdat.VVEL.TLAT[:,0], z_w_top_auxgrd, MV_projauxgrd[:,:,59])
fig, ax = utils_plt.plot_MOC(ncdat.VVEL.TLAT[:,0], z_w_top_auxgrd, ncdat.VVEL[0,:,:,59])