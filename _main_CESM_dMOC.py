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

from UTILS_misc import loadgetsave as LGS

# suppress RuntimeWaring due to All-NaN slices.
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', r'All-NaN axis encountered')
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', r'All-NaN slice encountered')

# =======================================================================================
#  Load netcdf file
# =======================================================================================
fpath='../data/'
fname='b40.lm850-1850.1deg.001.pop.h.1279.ann.4.cdf'
ncdat = xr.open_dataset(fpath+fname, decode_times=False)

# =======================================================================================
#  Settings and Directories
# =======================================================================================
auxgrd_name = ['lat395model', 'lat198model', 'lat170eq80S90N', 'lat340eq80S90N'][1] # choose aux grid
dbsetup = np.array([[20, 33, 14], [33.1, 36, 11], [36.05, 37.5, 11], [38, 43, 11]]) # setup for density bins [lower, upper, steps]
densChoice = 'sig2'     # string | choice of density to use for resampling (either RHO or sig2) 
boolLGSnoload = False  # boolean | False: normal loadgetsave procedure | True: noload-mode > will get and save (overwrite!)
# ---------------------------------------------------------------------------------------
# automatical generation of directory names
dir_mgrd         = paths.get_path2vars('mgrd', mkdir=True)
str_db          = 'spec_%sto%sin%s_%sto%sin%s_%sto%sin%s_%sto%sin%s' % tuple(dbsetup.flatten())
dir_dens        = paths.get_path2vars('dens', mkdir=True)
dir_auxgrd      = paths.get_path2vars(auxgrd_name, mkdir=True)
# paths (directory + filenames) to temp variables
path_MW_z_t     = dir_mgrd+'MW_z_t'
path_dMW        = dir_dens+'dMW_'+densChoice+'_'+str_db
path_dMV_proj   = dir_dens+'dMV_proj_'+densChoice+'_'+str_db
path_dMVf_proj  = dir_dens+'dMVf_proj_'+densChoice+'_'+str_db
path_zdb        = dir_dens+'zdb_'+str_db
path_zdbc       = dir_dens+'zdbc_'+str_db
path_HT_auxgrd_xmax = dir_auxgrd+'HT_auxgrd_xmax'
path_HT_mgrd_xmax   = dir_mgrd+'HT_mgrd_xmax'
path_HU_auxgrd_xmax = dir_auxgrd+'HU_auxgrd_xmax'
path_HU_mgrd_xmax   = dir_mgrd+'HU_mgrd_xmax'

# =======================================================================================
#  Transformation on different grids (Density, Spatial auxiliary grid)
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Mask for Atlantic
ATLboolmask = utils_mask.get_ATLbools(ncdat.REGION_MASK.values) # boolean mask
# ---------------------------------------------------------------------------------------
# - Spatial auxiliary grid
lat_auxgrd, zT_auxgrd, z_w_auxgrd = utils_mask.gen_auxgrd(ncdat, auxgrd_name)
lat_mgrd = ncdat.TLAT.isel(nlon=0)          # mean of LAT for each j #! very inappropriate
# ---------------------------------------------------------------------------------------
# - Density grid/bins
SA = ncdat.SALT[0,:,:,:].values             # absolute salinity
PT = ncdat.TEMP[0,:,:,:].values             # potential temperature
CT = gsw.CT_from_pt(SA, PT)                 # conservative temperature
sig2 = gsw.sigma2(SA, CT)                   # potential density anomaly referenced to 2000dbar
RHO = ncdat.RHO[0,:,:,:].values*1000-1000   # in-situ density anomaly [SI]
if densChoice == 'sig2': dens = sig2
elif densChoice == 'rho': dens = RHO
# density bins (db), center-values (dbc) and thicknesses (ddb, ddbc)
db = np.concatenate((np.linspace(dbsetup[0,0],dbsetup[0,1],dbsetup[0,2]), np.linspace(dbsetup[1,0],dbsetup[1,1],dbsetup[1,2]), np.linspace(dbsetup[2,0],dbsetup[2,1],dbsetup[2,2]), np.linspace(dbsetup[3,0],dbsetup[3,1],dbsetup[3,2])))
dbc = np.array([np.mean(db[i-1:i+1]) for i in np.arange(1,len(db))]) #! center-values | reasonable for non-eq-spaced db?
ddb = utils_ana.canonical_cumsum(np.diff(db)/2, 2, axis=0, crop=True)        # layer thickness of density_bins
ddbc = np.diff(db)                                                   # layer thickness from midpoint to midpoint (#! note: it is 1 element longer than ddb)
# depth of isopycnals (i.e. of density_bins at both staggered grids)
z_t_3d = utils_conv.exp_k_to_kji(ncdat.z_t, dens.shape[-2], dens.shape[-1])
zdb = LGS(lambda: utils_conv.resample_colwise(z_t_3d, dens, db, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort'), path_zdb, 'zdb', noload=boolLGSnoload)
zdbc = LGS(lambda: utils_conv.resample_colwise(z_t_3d, dens, dbc, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort'), path_zdbc, 'zdbc', noload=boolLGSnoload)
#dzdb = utils_ana.canonical_cumsum(np.diff(zdb, axis=0)/2, 2, axis=0, crop=True)      # layer thickness of isopycnal depths
dzdb = np.diff(zdb, axis=0)
del PT, CT, z_t_3d, dbsetup
# ---------------------------------------------------------------------------------------
# Volume representation for axis-scaling
#! check if and where db to be replaced by dbc
dz3d = utils_conv.exp_k_to_kji(ncdat.dz, dens.shape[-2], dens.shape[-1])    # in cgs
TAREA3d = utils_conv.exp_ji_to_kji(ncdat.TAREA, dens.shape[0])              # in cgs
vol3d = dz3d*TAREA3d                                                        # in cgs
inds = np.digitize(dens, db)
vol_dbs_glob = np.zeros(shape=[len(db)])
vol_dbs_reg  = np.zeros(shape=[len(db)])
vol_dbs_col  = np.zeros(shape=[len(db), dens.shape[-2], dens.shape[-1]])
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
#  Variables contained in model output
# =======================================================================================
# - temperature -------------------------------------------------------------------------
T = ncdat.TEMP.mean(dim='time')
T = utils_mask.mask_ATLANTIC(T, ncdat.REGION_MASK)
#T_dens = utils_conv.resample_colwise(T.values, dens, db, method='lin', fill_value=np.nan, mask=ATLboolmask, mono_method='sort')
# - Zonal maxima of ocean depth ---------------------------------------------------------
HT_auxgrd_xmax = LGS(lambda: utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'T', path_HT_auxgrd_xmax), path_HT_auxgrd_xmax, 'HT_auxgrd_xmax', noload=boolLGSnoload)
HT_mgrd_xmax = LGS(lambda: utils_mask.calc_H_mgrd_xmax(lat_mgrd, ncdat, 'T', dir_mgrd), path_HT_mgrd_xmax, 'HT_mgrd_xmax', noload=boolLGSnoload)
HU_auxgrd_xmax = LGS(lambda: utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'U', dir_auxgrd), path_HU_auxgrd_xmax, 'HU_auxgrd_xmax', noload=boolLGSnoload)
HU_mgrd_xmax = LGS(lambda: utils_mask.calc_H_mgrd_xmax(lat_mgrd, ncdat, 'U', path_mgrd), path_HU_mgrd_xmax, 'HU_mgrd_xmax', noload=boolLGSnoload)

# =======================================================================================
#  Volume transports (in Sv)
# =======================================================================================
# - on model grid
MV_mgrd         = utils_transp.calc_MV(ncdat)
MVf_mgrd        = utils_transp.calc_MVflat(ncdat)
MW_mgrd         = utils_transp.calc_MW(ncdat)
# - MV projected on auxiliary grid i.e. towards N (same shape as on model grid)
MV_proj         = utils_conv.project_on_auxgrd(MV_mgrd.values, ncdat.ANGLE.values)
MVf_proj        = utils_conv.project_on_auxgrd(MVf_mgrd.values, ncdat.ANGLE.values)
# - MW resampled on centre of T grid
MW_z_t          = LGS(lambda: utils_conv.resample_colwise(MW_mgrd.values, MW_mgrd.z_w_top.values, ncdat.z_t.values, method='lin', mask=ATLboolmask, mono_method='sort'), path_MW_z_t, 'MW_z_t', noload=boolLGSnoload)
# - resampled on density axis
boolLGSnoload = True  # noload-mode
 #dMV_proj        = LGS(lambda: utils_conv.resample_colwise(MV_proj, ncdat.z_t.values, zdbc, method='dMWV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMV_proj, 'dMV_proj', noload=boolLGSnoload)
dMVf_proj       = LGS(lambda: utils_conv.resample_colwise(MVf_proj, ncdat.z_t.values, zdbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMVf_proj, 'dMVf_proj', noload=boolLGSnoload)
boolLGSnoload = False  # normal LGS procedure
dMW             = LGS(lambda: utils_conv.resample_colwise(MW_z_t, ncdat.z_t.values, zdb, method='dMW_zdb', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMW, 'dMW', noload=boolLGSnoload)
 #dMW             = LGS(lambda: utils_conv.resample_colwise(MW_z_t, dens, db, method='dMW_db', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMW, noload=boolLGSnoload)

# =======================================================================================
#  Streamfunctions - V method
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - BSF Streamfunction (in Sv)
BSF_mgrd, MVzint = utils_BSF.calc_BSF_mgrd(MV_mgrd, dump_MVzint=True)

# ---------------------------------------------------------------------------------------
# - MOC Streamfunction (in depth space) (in Sv)
#   (1) zonal integration
MVxint_mgrd     = np.nansum(MV_proj, axis=2)
MVxint_auxgrd   = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, zT_auxgrd, 'V', MV_proj, ncdat, dir_auxgrd)
#   (2) vertical integration
MOC_mgrd_V      = utils_ana.nancumsum(MVxint_mgrd, axis=0)
MOC_auxgrd_V    = utils_ana.nancumsum(MVxint_auxgrd, axis=0)

# ---------------------------------------------------------------------------------------
# - A dMOC Streamfunction (in density space) (in Sv)
#   (1) vertical integration of MV (first! as weighting known in depth-space)
MV_zint         = utils_ana.nancumsum(MV_proj, axis=0)
#   (2) resampling on density axis
dMV_zint        = utils_conv.resample_colwise(MV_zint, ncdat.z_t.values, zdb, method='dMW_db', fill_value=np.nan, mask = ATLboolmask, mono_method='force')
 #dMV_zint        = utils_conv.resample_colwise(MV_zint, dens, db, method='dMV_db', fill_value=np.nan, mask = ATLboolmask, mono_method='force')
#   (3) zonal integration
dMOC_mgrd_V     = np.nansum(dMV_zint, axis=2)
dMOC_auxgrd_V   = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, db, 'dV', dMV_zint, ncdat, dir_auxgrd)
# ---------------------------------------------------------------------------------------
# - B dMOC Streamfunction (in density space) (in Sv)
#   (1) zonal integration of dMV_proj
dMVxint_mgrd    = np.nansum(dMV_proj, axis=2)
dMVxint_auxgrd  = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, db, 'dV', dMVf_proj, ncdat, dir_auxgrd)
#   (2) integration of MV along density axis weighting with ddbc
dMOC_mgrd_V     = utils_ana.nancumsum((dMVxint_mgrd.T*np.append(ddbc,1)).T, axis=0)/np.sum(ddbc) #! pfusch!!
# ---------------------------------------------------------------------------------------
# - C dMOC Streamfunction (in density space) (in Sv)
#   (1) integration of MVf_proj along density axis weighting with zdbc
dMV_proj_dint   = utils_ana.nancumsum(dMVf_proj*np.append(ddb,np.ones([1,384,320]), axis=0), axis=0)
#   (2) zonal integration of dMV_proj
dMOC_mgrd_V    = np.nansum(dMV_proj_dint, axis=2)
dMOC_auxgrd_V  = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, zdbc, 'dV', dMV_proj_dint, ncdat, dir_auxgrd)

# =======================================================================================
#  Streamfunctions - W method
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - MOC Streamfunction (in depth space) (in Sv)
#   (1) zonal integration
MWxint_mgrd     = np.nansum(MW_mgrd, axis=2)
MWxint_auxgrd   = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, z_w_auxgrd, 'W', MW_mgrd.values, ncdat, dir_auxgrd)
#   (2) meridional integration
MOC_auxgrd_W    = utils_MOC.normalise(utils_ana.nancumsum(MWxint_auxgrd, axis=1), ref='N')

# ---------------------------------------------------------------------------------------
# - dMOC Streamfunction (in density space) (in Sv)
#   (1) zonal integration
dMWxint_mgrd    = np.nansum(dMW, axis=2)
dMWxint_auxgrd  = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, db, 'dW', dMW, ncdat, dir_auxgrd)
#   (2) meridional integration
dMOC_mgrd_W     = utils_MOC.normalise(utils_ana.nancumsum(dMWxint_mgrd, axis=1), ref='N')
dMOC_auxgrd_W   = utils_MOC.normalise(utils_ana.nancumsum(dMWxint_auxgrd, axis=1), ref='N')
