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
import CESM_paths as paths
from IPython.core.debugger import Tracer; debug_here = Tracer()

# #######################################################################################
#  GET AND PROCESS DATA
# #######################################################################################
# ---------------------------------------------------------------------------------------
# load netcdf file
#fpath=paths.get_path2data('lm_1deg', 'anndat')
#fname='b40.lm850-1850.1deg.001.pop.h.1499.ann.4.cdf'
fpath='../../'
fname='b40.lm850-1850.1deg.001.pop.h.1279.ann.4.cdf'
ncdat = xr.open_dataset(fpath+fname, decode_times=False)

# =======================================================================================
#  Transformation on different grids (Density, Spatial auxiliary grid)
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Spatial auxiliary grid
auxgrd_name = ['lat395model_zeq60', 'lat198model_zeq60', 'lat170eq80S90N_zeq60', 'lat340eq80S90N_zeq60'][1]       # choose aux grid
lat_auxgrd, zT_auxgrd, z_w_top_auxgrd = utils_mask.gen_auxgrd(ncdat, auxgrd_name)
lat_mgrd = ncdat.TLAT.isel(nlon=0)          # mean of LAT for each j #! very inappropriate
# ---------------------------------------------------------------------------------------
# - Density grid
SA = ncdat.SALT[0,:,:,:].values             # absolute salinity
PT = ncdat.TEMP[0,:,:,:].values             # potential temperature
CT = gsw.CT_from_pt(SA, PT)                 # conservative temperature
sig2 = gsw.sigma2(SA, CT)                   # potential density anomaly referenced to 2000dbar
RHO = ncdat.RHO[0,:,:,:].values*1000-1000   # in-situ density anomaly [SI]
#dens_bins = np.linspace(28,42,100)          # dens_bins = np.linspace(1.004,1.036,65)[np.mean(b[i-1:i+1]) for i in np.arange(1,len(b))]
dens_bins = np.concatenate((np.arange(28,35), np.linspace(35, 38, 50), np.arange(39, 43)))
dens_bins_centers = np.array([np.mean(dens_bins[i-1:i+1]) for i in np.arange(1,len(dens_bins))]) #! reasonable for non-eq-spaced dens_bins?
dens_str = 'sig2'                               # string | choice of density to use for resampling (either RHO or sig2)


# =======================================================================================
# Pathnames for temporally stored variables
# =======================================================================================
path_auxgrd = paths.get_path2vars(auxgrd_name, True)
path_grd = paths.get_path2vars('grd', mkdir=True)
path_dens = paths.get_path2vars('dens', mkdir=True)
#varname_binning = 'eqbins_{}to{}in{}steps'.format(int(dens_bins.min()), int(dens_bins.max()), int(len(dens_bins)))
varname_binning = 'spec_{}to{}in{}steps'.format(int(dens_bins.min()), int(dens_bins.max()), int(len(dens_bins)))
fname_MWdens = 'MW_'+dens_str+'_'+varname_binning

# =======================================================================================
#  Variables contained in model output
# =======================================================================================
# - temperature -------------------------------------------------------------------------
T = ncdat.TEMP.mean(dim='time').isel(z_t=0)
T = utils_mask.mask_ATLANTIC(T, ncdat.REGION_MASK)
# - in-situ density ---------------------------------------------------------------------
rho = utils_mask.mask_ATLANTIC(ncdat.RHO.mean(dim='time'), ncdat.REGION_MASK)
rho = rho.mean(dim='nlon')








# =======================================================================================
#  Streamfunctions
# =======================================================================================
''' BSF: Barotropic Streamfunction
    MOC: Meridional Overturning Circulation Streamfunction
    MW:  vertical volume transport
    MV:  meridional volume transport
'''
# ---------------------------------------------------------------------------------------
# - Volume transports (in Sv)
MV_mgrd = utils_transp.calc_MV(ncdat)                                           # on model grid
MV_projauxgrd = utils_conv.project_on_auxgrd(MV_mgrd, ncdat.ANGLE.values)       # on auxiliary grid
MW_mgrd = utils_transp.calc_MW(ncdat)                                                # valid on both grids

try:    MW_dens = utils_misc.loadvar(path_dens+fname_MWdens)                    # load from file
except:
    print(' > loading failed!')
    # resampled MW_mgrd on centre of T grid
    MW_z_t = utils_conv.resample_colwise(MW_mgrd.values, MW_mgrd.z_w_top.values, ncdat.z_t.values, method='wmean')
    # resampled MW_mgrd on density axis (still pointing in vertical direction)
    # a) weighted mean of closest neighbours around dens_bin values
    MW_dens_binborders = utils_conv.resample_colwise(MW_z_t, sig2, dens_bins, method='MW', onlyposgrad='True')
    MW_dens_binborders[np.isnan(MW_dens_binborders)] = 0 # overwrite nans with 0 --> should be done in preallocation already!
    # b) integrating the differences in MW_mgrd from the densest waters towards lighter water (assuming that there's no denser sig2 than dens_bins[-1] )
    MW_dens_diff = np.diff(MW_dens_binborders, axis=0)
    MW_dens = np.cumsum(MW_dens_diff[::-1], axis=0)[::-1]
    # saving
    utils_misc.savevar(MW_dens, path_dens+fname_MWdens)                         # save to file











# ---------------------------------------------------------------------------------------
# - Streamfunctions (in Sv)...
# ... on model grid
BSF_mgrd, MVzint = utils_BSF.calc_BSF_mgrd(MV_mgrd, dump_MVzint=True)
MOC_mgrd_W, MWxint_mgrd = utils_MOC.calc_MOC_mgrd('W', MW_mgrd, do_norm=True, dump_Mxint=True)
 #MOC_mgrd_V, MVxint_mgrd = utils_MOC.calc_MOC_mgrd('V', MV_projauxgrd, do_norm=True, dump_Mxint=True)
dMOC_mgrd_W, dMOC_mgrd_W_norm, dMWxint_mgrd = utils_MOC.calc_MOC_mgrd_nparray('W', MW_dens, dump_Mxint=True)

# ... on auxiliary grid
MWxint_auxgrd = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, z_w_top_auxgrd, 'W', MW_mgrd.values, ncdat, path_auxgrd)
MOC_auxgrd_W, MOC_auxgrd_W_norm = utils_MOC.calc_MOC_auxgrd(lat_auxgrd, z_w_top_auxgrd, 'W', MWxint_auxgrd, 'forward', path_auxgrd)
 #MVxint_auxgrd = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, zT_auxgrd, 'V', MV_projauxgrd.values, ncdat, path_auxgrd)
 #MOC_auxgrd_V, MOC_auxgrd_V_norm = utils_MOC.calc_MOC_auxgrd(lat_auxgrd, zT_auxgrd, 'V', MVxint_auxgrd, path_auxgrd)
dMWxint_auxgrd = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, dens_bins_centers, 'dW', MW_dens, ncdat, path_auxgrd)
dMOC_auxgrd_W, dMOC_auxgrd_W_norm = utils_MOC.calc_MOC_auxgrd(lat_auxgrd, dens_bins_centers, 'W', dMWxint_auxgrd, 'forward', path_auxgrd)

# ... old stuff with utils_dMOC
 #dMWxint_auxgrd = utils_dMOC.calc_dMxint_auxgrd(lat_auxgrd, zT_auxgrd, 'W', MW_mgrd.values, sig2, dens_bins, ncdat, path_auxgrd)
 #dMOC_auxgrd_W = utils_dMOC.calc_dMOC_auxgrd(lat_auxgrd, dens_bins, 'W', dMWxint_auxgrd, ncdat, path_auxgrd, do_norm=False)
 #dMOC_mgrd_W, dMxint_mgrd = utils_dMOC.calc_dMOC_mgrd('W', MW_mgrd.values, sig2, dens_bins, do_norm=False, dump_dMxint=True)
 #dMOC_mgrd_W_norm = dMOC_mgrd_W - np.tile(dMOC_mgrd_W[:,-1],(384,1)).T

# =======================================================================================
#  Zonal maxima of ocean depth
# =======================================================================================
try:    HT_auxgrd_xmax = utils_misc.loadvar(path_auxgrd+'HT_auxgrd_xmax')       # load from file
except: HT_auxgrd_xmax = utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'T', path_auxgrd)
try:    HT_mgrd_xmax = utils_misc.loadvar(path_grd+'HT_mgrd_xmax')             # load from file
except: HT_mgrd_xmax = utils_mask.calc_H_mgrd_xmax(ncdat, 'T', path_grd)
#try:    HU_auxgrd_xmax = utils_misc.loadvar(path_auxgrd+'HU_auxgrd_xmax')       # load from file
#except: HU_auxgrd_xmax = utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'U', path_auxgrd)
#try:    HU_mgrd_xmax = utils_misc.loadvar(path_grd+'HU_mgrd_xmax')             # load from file
#except: HU_mgrd_xmax = utils_mask.calc_H_mgrd_xmax(ncdat, 'U', path_grd)


# #######################################################################################
#  PLOTTING
# #######################################################################################
plt.ion() # enable interactive mode
path_fig = '../figures/discussion_160701/'
# =======================================================================================
#  CCSM4 representations
# =======================================================================================
# -----------------------------------------------------------------------------------------
# BSF on geographical grid calculated by model
BSF_model = utils_mask.mask_ATLANTIC(ncdat.BSF.isel(time=0), ncdat.REGION_MASK)
fig, map = utils_plt.plot_BSF(BSF_model, 'T', nlevels = 10)
plt.title('BSF model on T grid')
utils_plt.print2pdf(fig, path_fig+'BSF_model_T')
# -----------------------------------------------------------------------------------------
# MOC on geographical grid calculated by model
MOC_model = ncdat.MOC.isel(time=0, transport_reg=1, moc_comp=0)
#MOC_model = MOC_model - MOC_model[:,-1] # normalisation
fig, ax = utils_plt.plot_MOC(MOC_model.lat_aux_grid, MOC_model.moc_z, MOC_model, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HT_auxgrd_xmax) 				# plot seafloor
plt.title('MOC model')
plt.xlim([-36,90])
utils_plt.print2pdf(fig, path_fig+'MOC_model')

# =======================================================================================
#  Calculated on model grid
# =======================================================================================
# -----------------------------------------------------------------------------------------
# BSF on model grid
fig, map = utils_plt.plot_BSF(BSF_mgrd, 'U', nlevels=100)
plt.title('BSF mgrd on U grid')
utils_plt.print2pdf(fig, 'testfigures/BSF_mgrd_U')
# -----------------------------------------------------------------------------------------
# MOC_mgrd_W
fig, ax = utils_plt.plot_MOC(lat_mgrd, ncdat.z_w_top, MOC_mgrd_W, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_mgrd,HT_mgrd_xmax) 				# plot seafloor #! it's the T-grid!!!
plt.title('MOC mgrd W')
plt.xlim([-36,90])
utils_plt.print2pdf(fig, path_fig+'MOC_mgrd_W')
# -----------------------------------------------------------------------------------------
# MOC_mgrd_V
fig, ax = utils_plt.plot_MOC(lat_mgrd, ncdat.z_t, MOC_mgrd_V, nlevels=100, plttype='pcolor+contour')
plt.plot(lat_mgrd,HT_mgrd_xmax)				# plot seafloor
plt.title('MOC mgrd V')
plt.xlim([-36,90])
 #utils_plt.print2pdf(fig, path_fig+'MOC_mgrd_V')
# -----------------------------------------------------------------------------------------
# dMOC on model grid (in Sv)
fig, ax = utils_plt.plot_MOC(lat_mgrd, dens_bins_centers, dMOC_mgrd_W_norm, nlevels=10, plttype='pcolor+contour')
plt.title('dMOC mgrd W (sigma2)')
plt.suptitle('density binning from {} to {} in {} steps'.format(dens_bins_centers.min(), dens_bins_centers.max(), len(dens_bins_centers)))
plt.xlim([-36,73])
utils_plt.print2pdf(fig, path_fig+'dMOC_mgrd_W_sig2')

# =======================================================================================
#  Calculated on auxiliary (geographical) grid
# =======================================================================================
# -----------------------------------------------------------------------------------------
# MOC_auxgrd_W
fig, ax = utils_plt.plot_MOC(lat_auxgrd, z_w_top_auxgrd, MOC_auxgrd_W, nlevels=10, plttype='pcolor+contour')
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
fig, ax = utils_plt.plot_MOC(lat_auxgrd, dens_bins_centers, dMOC_auxgrd_W_norm, nlevels=10, plttype='pcolor+contour')
plt.title('dMOC auxgrd W (sigma2)')
plt.suptitle('density binning from {} to {} in {} steps'.format(dens_bins_centers.min(), dens_bins_centers.max(), len(dens_bins_centers)))
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
fig, map = utils_plt.pcolor_basemap(MW_mgrd.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
fig, map = utils_plt.pcolor_basemap(MW_mgrd.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
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
# MW_mgrd
fig, map = utils_plt.pcolor_basemap(MW_mgrd.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
plt.title('MW_mgrd')
 #utils_plt.print2pdf(fig, 'testfigures/MW_mgrd')
# -----------------------------------------------------------------------------------------
# ANGLE
plt.figure()
plt.pcolor(ncdat.ANGLE*180/np.pi)
plt.contour(ncdat.REGION_MASK)

fig, ax = utils_plt.plot_MOC(ncdat.VVEL.TLAT[:,0], z_w_top_auxgrd, MV_projauxgrd[:,:,59])
fig, ax = utils_plt.plot_MOC(ncdat.VVEL.TLAT[:,0], z_w_top_auxgrd, ncdat.VVEL[0,:,:,59])