'''
CESM controls for calculation of PSI in density coordinates

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

# suppress RuntimeWaring due to All-NaN slices.
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
    
# #######################################################################################
#  GET AND PROCESS DATA
# #######################################################################################
# ---------------------------------------------------------------------------------------
# load netcdf file
fpath='../data/'
fname='b40.lm850-1850.1deg.001.pop.h.1279.ann.4.cdf'
ncdat = xr.open_dataset(fpath+fname, decode_times=False)

# =======================================================================================
#  Transformation on different grids (Density, Spatial auxiliary grid)
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Mask for Atlantic
ATLboolmask = utils_mask.get_ATLbools(ncdat.REGION_MASK.values) # boolean mask
# ---------------------------------------------------------------------------------------
# - Spatial auxiliary grid
auxgrd_name = ['lat395model_zeq60', 'lat198model_zeq60', 'lat170eq80S90N_zeq60', 'lat340eq80S90N_zeq60'][1]       # choose aux grid
lat_auxgrd, zT_auxgrd, z_w_top_auxgrd = utils_mask.gen_auxgrd(ncdat, auxgrd_name)
lat_mgrd = ncdat.TLAT.isel(nlon=0)          # mean of LAT for each j #! very inappropriate
# ---------------------------------------------------------------------------------------
# - Density grid/bins
SA = ncdat.SALT[0,:,:,:].values             # absolute salinity
PT = ncdat.TEMP[0,:,:,:].values             # potential temperature
CT = gsw.CT_from_pt(SA, PT)                 # conservative temperature
sig2 = gsw.sigma2(SA, CT)                   # potential density anomaly referenced to 2000dbar
RHO = ncdat.RHO[0,:,:,:].values*1000-1000   # in-situ density anomaly [SI]
# density bins, center-values and thicknesses
db = np.concatenate((np.linspace(20, 33, 14), np.linspace(33.1, 36, 11), np.linspace(36.05,37.5, 11), np.linspace(38, 43, 11)))
dbc = np.array([np.mean(db[i-1:i+1]) for i in np.arange(1,len(db))]) #! center-values | reasonable for non-eq-spaced db?
ddb = utils_ana.canonical_cumsum(np.diff(db)/2, 2, crop=True)        # layer thickness of density_bins
ddbc = np.diff(db)                                                   # layer thickness from midpoint to midpoint (#! note: it is 1 element longer than ddb)
# depth of isopycnals (i.e. of density_bins at both staggered grids)
z_t_3d = utils_conv.expand_karray_to_kji(ncdat.z_t, sig2.shape[-2], sig2.shape[-1])
zdb = utils_conv.resample_colwise(z_t_3d, sig2, db, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort')
zdbc = utils_conv.resample_colwise(z_t_3d, sig2, dbc, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort')

del PT, CT, z_t_3d
# ---------------------------------------------------------------------------------------
# Volume representation for axis-scaling
#! check if and where db to be replaced by dbc
dz3d = utils_conv.expand_karray_to_kji(ncdat.dz, sig2.shape[-2], sig2.shape[-1])    # in cgs
TAREA3d = utils_conv.expand_jiarray_to_kji(ncdat.TAREA, sig2.shape[0])              # in cgs
vol3d = dz3d*TAREA3d                                                                # in cgs
inds = np.digitize(sig2, db)
vol_dbs_glob = np.zeros(shape=[len(db)])
vol_dbs_reg  = np.zeros(shape=[len(db)])
vol_dbs_col  = np.zeros(shape=[len(db), sig2.shape[-2], sig2.shape[-1]])
for b in np.arange(len(db)):
    vol_dbs_glob[b] = np.sum(vol3d[inds==b])                # global        # in cgs
    vol_dbs_reg[b]  = np.sum((vol3d*ATLboolmask)[inds==b])  # regional      # in cgs
    vol_dbs_col[b]  = np.sum(vol3d[inds==b], axis=0)        # column-wise   # in cgs

# axes and ticks
ax_vol_glob = np.cumsum(vol_dbs_glob) - vol_dbs_glob/2
ax_vol_reg = np.cumsum(vol_dbs_reg) - vol_dbs_reg/2
ticks_dens = [28, 35.13, 36, 36.485, 36.92, 37.065, 37.21, 37.355]
ticks_vol_glob = ax_vol_glob[np.in1d(db, ticks_dens)]
ticks_vol_reg = ax_vol_reg[np.in1d(db, ticks_dens)]

del dz3d, TAREA3d, vol3d, inds

# =======================================================================================
# Pathnames for temporally stored variables
# =======================================================================================
dens_str = 'sig2'   # string | choice of density to use for resampling (either RHO or sig2) 

path_auxgrd = paths.get_path2vars(auxgrd_name, mkdir=True)
path_grd = paths.get_path2vars('grd', mkdir=True)
path_dens = paths.get_path2vars('dens', mkdir=True)
varname_binning = 'spec_{}to{}in{}_{}to{}in{}_{}to{}in{}_{}to{}in{}'.format(28,33,14, 33.1,36,11, 36.05,37.5,11, 38,43,11)

fname_MWdens = 'MW_'+dens_str+'_'+varname_binning
fname_MWzt = 'MW_z_t'    # MW interpolated on z_t i.e. on T-cell centres

# =======================================================================================
#  Variables contained in model output
# =======================================================================================
# - temperature -------------------------------------------------------------------------
T = ncdat.TEMP.mean(dim='time')
T = utils_mask.mask_ATLANTIC(T, ncdat.REGION_MASK)
T_dens = utils_conv.resample_colwise(T.values, sig2, db, method='lin', fill_value=np.nan, mask=ATLboolmask, mono_method='sort')


# =======================================================================================
#  Transports and Streamfunctions - V
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Volume transports (in Sv)
MV_mgrd = utils_transp.calc_MV(ncdat)                                          # on model grid
MV_projauxgrd = utils_conv.project_on_auxgrd(MV_mgrd, ncdat.ANGLE.values)      # projected on auxiliary grid i.e. towards N (same shape as on model grid)

# ---------------------------------------------------------------------------------------
# - BSF Streamfunction (in Sv)
BSF_mgrd, MVzint = utils_BSF.calc_BSF_mgrd(MV_mgrd, dump_MVzint=True)

# ---------------------------------------------------------------------------------------
# - MOC Streamfunction (in Sv)...
# ... in depth-space
MOC_mgrd_V, MVxint_mgrd = utils_MOC.calc_MOC_mgrd('V', MV_projauxgrd, do_norm=True, dump_Mxint=True)
MVxint_auxgrd = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, zT_auxgrd, 'V', MV_projauxgrd.values, ncdat, path_auxgrd)
MOC_auxgrd_V, MOC_auxgrd_V_norm = utils_MOC.calc_MOC_auxgrd(lat_auxgrd, zT_auxgrd, 'V', MVxint_auxgrd, 'forward', path_auxgrd)
# ... in density-space
MV_cumsum = np.cumsum(MV_projauxgrd, axis=0)
dMV = utils_conv.resample_colwise(MV_cumsum.values, ncdat.z_t.values, zdb, method='dMW_zdb', fill_value=np.nan, mask = ATLboolmask, mono_method='force')

dMOC_mgrd_V, dMOC_mgrd_V_norm, dMVxint_mgrd = utils_MOC.calc_MOC_mgrd_nparray('V', dMV, dump_Mxint=True)
dMVxint_auxgrd = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, db, 'dV', dMV, ncdat, path_auxgrd)
dMOC_auxgrd_V = dMVxint_auxgrd
dMOC_auxgrd_V_norm = dMOC_auxgrd_V - np.tile(dMOC_auxgrd_V[:,-1],(dMOC_auxgrd_V.shape[1],1)).T

'''
# =======================================================================================
#  Transports and Streamfunctions - W
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Volume transports (in Sv)...
# ... on model grid
MW_mgrd = utils_transp.calc_MW(ncdat)
# ... resampled on centre of T grid
MW_z_t = utils_conv.resample_colwise(MW_mgrd.values, MW_mgrd.z_w_top.values, ncdat.z_t.values, method='lin', mask = ATLboolmask, mono_method='sort')
# ... resampled on density axis
#     --> 2 versions
dMW = utils_conv.resample_colwise(MW_z_t, ncdat.z_t.values, zdb, method='dMW_zdb', fill_value=np.nan, mask = ATLboolmask, mono_method='force')
 #dMW = utils_conv.resample_colwise(MW_z_t, sig2, db, method='dMW_db', fill_value=np.nan, mask = ATLboolmask, mono_method='force')

# ---------------------------------------------------------------------------------------
# - BSF Streamfunction (in Sv)
BSF_mgrd, MVzint = utils_BSF.calc_BSF_mgrd(MV_mgrd, dump_MVzint=True)

# ---------------------------------------------------------------------------------------
# - MOC Streamfunction (in Sv)...
# ... in depth-space
MOC_mgrd_W, MWxint_mgrd = utils_MOC.calc_MOC_mgrd('W', MW_mgrd, do_norm=True, dump_Mxint=True)
MWxint_auxgrd = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, z_w_top_auxgrd, 'W', MW_mgrd.values, ncdat, path_auxgrd)
MOC_auxgrd_W, MOC_auxgrd_W_norm = utils_MOC.calc_MOC_auxgrd(lat_auxgrd, z_w_top_auxgrd, 'W', MWxint_auxgrd, 'forward', path_auxgrd)
# ... in density-space
dMOC_mgrd_W, dMOC_mgrd_W_norm, dMWxint_mgrd = utils_MOC.calc_MOC_mgrd_nparray('W', dMW, dump_Mxint=True)
dMWxint_auxgrd = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, db, 'dW', dMW, ncdat, path_auxgrd)
dMOC_auxgrd_W, dMOC_auxgrd_W_norm = utils_MOC.calc_MOC_auxgrd(lat_auxgrd, db, 'W', dMWxint_auxgrd, 'forward', path_auxgrd)

# ---------------------------------------------------------------------------------------
# - Saving variables to files
utils_misc.savevar(MW_z_t, path_dens+fname_MWzt)    # save to file
utils_misc.savevar(dMW, path_dens+fname_MWdens)     # save to file
'''

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
