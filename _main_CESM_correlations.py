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
import CESM_utils_time as utils_time
import CESM_utils_analysis as utils_ana
import CESM_paths as paths

# ---------------------------------------------------------------------------------------
# paths
path_mgrd = 'vars_mgrd/'
path_corr = 'vars_corr/'
fpath = paths.get_path2data('lm_1deg', 'anndat')
fnames = ['b40.lm850-1850.1deg.001.pop.h.{:04d}.ann.4.cdf'.format(i) for i in np.arange(850, 1500)]


# =======================================================================================
#  Load time independent variables
# =======================================================================================
ncdat = xr.open_dataset(fpath+fnames[0], decode_times=False)
TAREA = ncdat.TAREA

# =======================================================================================
#  Loop over time
# =======================================================================================
# load BSF and MOC (try from pickled file, otherwise load from ncdata)
try: 
    BSF_mod = utils_misc.loadvar(path_corr+'BSF_mod') 
    MOC_mod = utils_misc.loadvar(path_corr+'MOC_mod') 
except:
    BSF_mod = []
    MOC_mod = []
    for t in np.arange(len(fnames)):
        # load netcdf file
        ncdat = xr.open_dataset(fpath+fnames[t], decode_times=False)
        # write streamfunctions to variables (in Sv)
        BSF_mod = utils_time.concat(BSF_mod, utils_mask.mask_ATLANTIC(ncdat.BSF, ncdat.REGION_MASK))
        MOC_mod = utils_time.concat(MOC_mod, ncdat.MOC.isel(transport_reg=1, moc_comp=0))
    # save
    utils_misc.savevar(BSF_mod, path_corr+'BSF_mod')
    utils_misc.savevar(MOC_mod, path_corr+'MOC_mod')
# normalisation of AMOC
for t in np.arange(len(MOC_mod)):
    MOC_mod[t,:,:] = MOC_mod[t,:,:] - MOC_mod[t,:,-1]

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
TAREA = TAREA.where((BSF.TLAT>=45) & (BSF.TLAT<=70))
# ---------------------------------------------------------------------------------------
# calculate Streamfunction indices
MOCidx = MOC.max(dim = ['moc_z','lat_aux_grid'])
BSFidx = (BSF*TAREA).mean(dim=['nlat', 'nlon'])/TAREA.mean() # weighted average
# ---------------------------------------------------------------------------------------
# calculate correlation btw. MOCidx and BSFidx
corrIDX = utils_ana.xcorr(MOCidx, BSFidx, 0)
utils_misc.savevar(corrIDX, path_corr+'corrBSFidx')
# calculate pointwise correlation of MOC with BSFidx
corrBSFidx = np.zeros(shape = MOC.shape[-2:])
for j in np.arange(MOC.shape[-1]):
    utils_misc.ProgBar('step', step = j, nsteps = MOC.shape[-1])
    for k in np.arange(MOC.shape[-2]):
        corrBSFidx[k,j] = utils_ana.xcorr(MOC[:,k,j], BSFidx, 0)
utils_misc.ProgBar('done')
utils_misc.savevar(corrBSFidx, path_corr+'corrBSFidx')
# calculate pointwise correlation of BSF with MOCidx
corrMOCidx = np.zeros(shape = BSF.shape[-2:])
for j in np.arange(BSF.shape[-2]):
    utils_misc.ProgBar('step', step = j, nsteps = MOC.shape[-1])    
    for i in np.arange(BSF.shape[-1]):
        corrMOCidx[j,i] = utils_ana.xcorr(BSF[:,j,i], MOCidx, 0)
utils_misc.ProgBar('done')
utils_misc.savevar(corrMOCidx, path_corr+'corrMOCidx')


# #######################################################################################
#  PLOTTING
# #######################################################################################
plt.ion() # enable interactive mode
path_fig = 'corrfigs_figures_Jun22/'

# =======================================================================================
#  Zonal maxima of ocean depth
# =======================================================================================
try:    HT_mgrd_xmax = utils_misc.loadvar(path_mgrd+'HT_mgrd_xmax')             # load from file
except: HT_mgrd_xmax = utils_mask.calc_H_mgrd_xmax(ncdat, 'T', path_mgrd)
try:    HU_mgrd_xmax = utils_misc.loadvar(path_mgrd+'HU_mgrd_xmax')             # load from file
except: HU_mgrd_xmax = utils_mask.calc_H_mgrd_xmax(ncdat, 'U', path_mgrd)

# =======================================================================================
#  Streamfunctions
# =======================================================================================
lat_mgrd = ncdat.TLAT.isel(nlon=0)          # mean of LAT for each j #! very inappropriate

# -----------------------------------------------------------------------------------------
# BSF on geographical grid calculated by model
fig, map = utils_plt.plot_BSF(BSF.isel(time=0), 'T', nlevels = 10)
plt.title('BSF model on T grid')
 #utils_plt.print2pdf(fig, 'testfigures/BSF_model_T')
# -----------------------------------------------------------------------------------------
# MOC on geographical grid calculated by model
MOCsel = MOC.isel(time=0)
fig, ax = utils_plt.plot_MOC(MOCsel.lat_aux_grid, MOCsel.moc_z, MOCsel, nlevels=40, plttype='pcolor+contour')
plt.plot(lat_mgrd, HT_auxgrd_xmax) 				# plot seafloor
plt.title('MOC model')
plt.xlim([-36,90])
 #utils_plt.print2pdf(fig, 'testfigures/MOC_model')

# =======================================================================================
#  Streamfunction Indices
# =======================================================================================

# -----------------------------------------------------------------------------------------
# BSF on geographical grid calculated by model
fig, map = utils_plt.plot_BSF(corrBSFidx, 'T', nlevels = 10)
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
