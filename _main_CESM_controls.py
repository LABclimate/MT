'''
CESM controls --> computation routines

@author: buerki@climate.unibe.ch
TODO: 	 add '.values' where possible to speed up code.

Variables:
    BSF: Barotropic Streamfunction
    MOC: Meridional Overturning Circulation Streamfunction
    MW:  vertical volume transport
    MV:  meridional volume transport
'''

import numpy as np
import gsw
from netCDF4 import Dataset
import xarray as xr
import pickle
import sys
sys.path.append('/home/buerki/Documents/MT/scripts/')
import CESM_utils_analysis as utils_ana
import CESM_utils_mask as utils_mask
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
# - Mask for Atlantic
ATLboolmask = utils_mask.get_ATLbools(ncdat.REGION_MASK) # boolean mask
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
#dens_bins = np.concatenate((np.arange(28,35), np.linspace(35, 38, 50), np.arange(39, 43)))
dens_bins = np.concatenate((np.linspace(28, 33, 11), np.linspace(33.1, 37.5, 45), np.linspace(38, 43, 11)))
# variables related to density_bins
dens_bins_centers = np.array([np.mean(dens_bins[i-1:i+1]) for i in np.arange(1,len(dens_bins))]) #! reasonable for non-eq-spaced dens_bins?
ddb = utils_ana.canonical_cumsum(np.diff(dens_bins)/2, 2, crop=True)    # layer thickness of density_bins
ddb_centers = np.diff(dens_bins)                                        # layer thickness from midpoint to midpoint (#! note: it is 1 element longer than ddb)

# =======================================================================================
# Pathnames for temporally stored variables
# =======================================================================================
dens_str = 'sig2'                           # string | choice of density to use for resampling (either RHO or sig2) 

path_auxgrd = paths.get_path2vars(auxgrd_name, True)
path_grd = paths.get_path2vars('grd', mkdir=True)
path_dens = paths.get_path2vars('dens', mkdir=True)
#varname_binning = 'eqbins_{}to{}in{}steps'.format(int(dens_bins.min()), int(dens_bins.max()), int(len(dens_bins)))
varname_binning = 'spec_{}to{}in{}steps'.format(int(dens_bins.min()), int(dens_bins.max()), int(len(dens_bins)))
fname_MWdens = 'MW_'+dens_str+'_'+varname_binning
fname_MWzt = 'MW_z_t'    # MW interpolated on z_t i.e. on T-cell centres

# =======================================================================================
#  Variables contained in model output
# =======================================================================================
# - temperature -------------------------------------------------------------------------
T = ncdat.TEMP.mean(dim='time')
T = utils_mask.mask_ATLANTIC(T, ncdat.REGION_MASK)
T_dens = utils_conv.resample_colwise(T.values, sig2, dens_bins, method='wmean', fill_value=np.nan, mask=ATLboolmask, sort_ogrd='True')

# =======================================================================================
#  Streamfunctions
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Volume transports (in Sv)
MW_mgrd = utils_transp.calc_MW(ncdat)                                           # on model grid
MV_mgrd = utils_transp.calc_MV(ncdat)                                          # on model grid
#MV_projauxgrd = utils_conv.project_on_auxgrd(MV_mgrd, ncdat.ANGLE.values)      # projected on auxiliary grid (same shape as on model grid)

# - conversion on density axis
try:    MW_dens = utils_misc.loadvar(path_dens+fname_MWdens)                    # load from file
except:
    print(' > loading failed!')
    # resampled MW_mgrd on centre of T grid
    MW_z_t = utils_conv.resample_colwise(MW_mgrd.values, MW_mgrd.z_w_top.values, ncdat.z_t.values, method='wmean', mask = ATLboolmask)
    utils_misc.savevar(MW_z_t, path_dens+fname_MWzt)                           # save to file
    # resampled MW_mgrd on density axis (still pointing in vertical direction)
    MW_dens = utils_conv.resample_colwise(MW_z_t, sig2, dens_bins, method='dMW', fill_value=0, mask = ATLboolmask, sort_ogrd='True')
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

