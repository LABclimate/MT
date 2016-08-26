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
fname3='b30.004.pop.h.1000.ave.cdf'
fname4='b40.lm850-1850.1deg.001.pop.h.1279.ann.4.cdf'
ncdat3 = xr.open_dataset(fpath+fname3, decode_times=False)
ncdat4 = xr.open_dataset(fpath+fname4, decode_times=False)
CESMversion = 4
if CESMversion==3:      ncdat = ncdat3
elif CESMversion==4:    ncdat = ncdat4
# =======================================================================================
#  Settings and Directories
# =======================================================================================
auxgrd_name = ['lat395model', 'lat198model', 'lat170eq80S90N', 'lat340eq80S90N'][1] # choose aux grid
 #dbsetup = np.array([[20, 33, 14], [33.1, 36, 11], [36.05, 37.5, 11], [38, 43, 11]]) # setup for density bins [lower, upper, steps]
dbsetup = np.array([[12, 22, 3], [25,30, 3], [31,35,6], [35, 38, 40], [38.5,42, 3]]) # setup for density bins [lower, upper, steps]
densChoice = 'sig2'     # string | choice of density to use for resampling (either RHO or sig2) 
boolLGSnoload = True  # boolean | False: normal loadgetsave procedure | True: noload-mode > will get and save (overwrite!)
# ---------------------------------------------------------------------------------------
# automatical generation of directory names
str_db          = 'spec_%sto%sin%s_%sto%sin%s_%sto%sin%s_%sto%sin%s_%sto%sin%s' % tuple(dbsetup.flatten())

dir_mgrd        = paths.get_path2vars('mgrd',CESMversion=CESMversion, mkdir=True)
dir_dens        = paths.get_path2vars('dens', CESMversion=CESMversion, mkdir=True)
dir_auxgrd      = paths.get_path2vars(auxgrd_name, CESMversion=CESMversion, mkdir=True)
# paths (directory + filenames) to temp variables
path_dMVp       = dir_dens+'dMVp_'+densChoice+'_'+str_db
path_dMVfp_A    = dir_dens+'dMVfp_A_'+densChoice+'_'+str_db
path_dMVfp_B    = dir_dens+'dMVfp_B_'+densChoice+'_'+str_db
path_dMVfp_C    = dir_dens+'dMVfp_C_'+densChoice+'_'+str_db
path_zdbb       = dir_dens+'zdbb_'+str_db
path_zdbc       = dir_dens+'zdbc_'+str_db
path_zdbbc      = dir_dens+'zdbbc_'+str_db
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
lat_auxgrd, zT_auxgrd, z_w_top_auxgrd = utils_mask.gen_auxgrd(ncdat4, auxgrd_name)
lat_mgrd = ncdat.TLAT.isel(nlon=0)          # mean of LAT for each j #! very inappropriate
# ---------------------------------------------------------------------------------------
# - Density grid/bins
# - conversion T-->U
#! TODO
#! note that for volume representation T-grid dbc is needed!!!!
SA = ncdat.SALT[0,:,:,:].values             # absolute salinity
PT = ncdat.TEMP[0,:,:,:].values             # potential temperature
CT = gsw.CT_from_pt(SA, PT)                 # conservative temperature
sig2 = gsw.sigma2(SA, CT)                   # potential density anomaly referenced to 2000dbar
RHO = ncdat.RHO[0,:,:,:].values*1000-1000   # in-situ density anomaly [SI]
if densChoice == 'sig2': dens = sig2
elif densChoice == 'rho': dens = RHO
# density bins:  border-values (=dbb), center-values (=dbc) and thickness (=ddb)
dbb = np.concatenate((np.linspace(dbsetup[0,0],dbsetup[0,1],dbsetup[0,2]), np.linspace(dbsetup[1,0],dbsetup[1,1],dbsetup[1,2]), np.linspace(dbsetup[2,0],dbsetup[2,1],dbsetup[2,2]), np.linspace(dbsetup[3,0],dbsetup[3,1],dbsetup[3,2]), np.linspace(dbsetup[4,0],dbsetup[4,1],dbsetup[4,2])))
dbc = np.convolve(dbb, np.array([.5,.5]))[1:-1]
# depth of isopycnals (zdbbc) calculated as z(dbc) (=zdbc) and as c(zdbb) (=zdbbc)
z_t_3d = utils_conv.exp_k_to_kji(ncdat.z_t, dens.shape[-2], dens.shape[-1])
zdbc = LGS(lambda: utils_conv.resample_colwise(z_t_3d, dens, dbc, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort'), path_zdbc, 'zdbc', noload=boolLGSnoload)
zdbb = LGS(lambda: utils_conv.resample_colwise(z_t_3d, dens, dbb, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort'), path_zdbb, 'zdbb', noload=boolLGSnoload)
zdbbc = np.ones_like(zdbc) * np.nan
for j in np.arange(zdbb.shape[-2]):
    for i in np.arange(zdbb.shape[-1]):
        zdbbc[:,j,i] = np.convolve(zdbb[:,j,i], np.array([.5,.5]))[1:-1] # centre-values
# thickness of db (=ddb) and zdb (=dzdb)
ddb = np.diff(dbb)
dzdb = np.diff(zdbb, axis=0)
del PT, CT, z_t_3d, dbsetup

# ---------------------------------------------------------------------------------------
# Volume representation for axis-scaling (fully on T-grid)
dz3d = utils_conv.exp_k_to_kji(ncdat.dz, dens.shape[-2], dens.shape[-1])    # in cgs
TAREA3d = utils_conv.exp_ji_to_kji(ncdat.TAREA, dens.shape[0])              # in cgs
vol3d = dz3d*TAREA3d                                                        # in cgs
inds = np.digitize(dens, dbc)
# pre-allocation
Vdb_glob = np.zeros(shape=[len(dbc)])
Vdb_reg  = np.zeros(shape=[len(dbc)])
Vdb_col  = np.zeros(shape=[len(dbc), dens.shape[-2], dens.shape[-1]])
# summing up volumina over different basins
for b in np.arange(len(dbc)):
    Vdb_glob[b] = np.sum(vol3d[inds==b])                # global        # in cgs
    Vdb_reg[b]  = np.sum((vol3d*ATLboolmask)[inds==b])  # regional      # in cgs
    Vdb_col[b]  = np.sum(vol3d[inds==b], axis=0)        # column-wise   # in cgs
# axes and ticks for plots
ax_vol_glob = np.cumsum(Vdb_glob) - Vdb_glob/2
ax_vol_reg = np.cumsum(Vdb_reg) - Vdb_reg/2
ticks_dens = [20.5, 28.5, 35.275, 36.025, 36.2675, 37.1375, 37.75, 38.75]
ticks_dens = [33, 35, 36, dbc[35], dbc[38], dbc[40], dbc[41]]
ticks_vol_glob = ax_vol_glob[np.in1d(dbc, ticks_dens)]
ticks_vol_reg = ax_vol_reg[np.in1d(dbc, ticks_dens)]

del dz3d, TAREA3d, vol3d, inds

# =======================================================================================
#  Variables contained in model output
# =======================================================================================
# - temperature -------------------------------------------------------------------------
T = ncdat.TEMP.mean(dim='time')
T = utils_mask.mask_ATLANTIC(T, ncdat.REGION_MASK)
#T_dens = utils_conv.resample_colwise(T.values, dens, dbc, method='lin', fill_value=np.nan, mask=ATLboolmask, mono_method='sort')
# - Zonal maxima of ocean depth --------------------------------------------------------
boolLGSnoload = False  # boolean | False: normal loadgetsave procedure | True: noload-mode > will get and save (overwrite!)
HT_auxgrd_xmax  = LGS(lambda: utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'T', path_HT_auxgrd_xmax), path_HT_auxgrd_xmax, 'HT_auxgrd_xmax', noload=boolLGSnoload)
HT_mgrd_xmax    = LGS(lambda: utils_mask.calc_H_mgrd_xmax(lat_mgrd, ncdat, 'T', dir_mgrd), path_HT_mgrd_xmax, 'HT_mgrd_xmax', noload=boolLGSnoload)
HU_auxgrd_xmax  = LGS(lambda: utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'U', dir_auxgrd), path_HU_auxgrd_xmax, 'HU_auxgrd_xmax', noload=boolLGSnoload)
HU_mgrd_xmax    = LGS(lambda: utils_mask.calc_H_mgrd_xmax(lat_mgrd, ncdat, 'U', path_mgrd), path_HU_mgrd_xmax, 'HU_mgrd_xmax', noload=boolLGSnoload)
boolLGSnoload = True  # boolean | False: normal loadgetsave procedure | True: noload-mode > will get and save (overwrite!)

# =======================================================================================
#  Volume transports (in Sv)
# =======================================================================================
# - MV on model grid
MV_mgrd     = utils_transp.calc_MV(ncdat)       # = V * DX *DZ
MVf_mgrd    = utils_transp.calc_MVflat(ncdat)   # = V * DX
# - MV projected on auxiliary grid i.e. towards N (same shape as on model grid)
#MVp         = utils_conv.project_on_auxgrd(MV_mgrd.values, ncdat.ANGLE.values)
MVp         = MV_mgrd.values
MVfp        = utils_conv.project_on_auxgrd(MVf_mgrd.values, ncdat.ANGLE.values)
# - resampled on density axis
dMVp        = LGS(lambda: utils_conv.resample_colwise(MVp, ncdat.z_t.values, zdbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMVp, 'dMVp', noload=boolLGSnoload)
 #dMVfp_A     = LGS(lambda: utils_conv.resample_colwise(MVfp, dens, dbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMVfp_A, 'dMVfp_A', noload=boolLGSnoload)
dMVfp_B     = LGS(lambda: utils_conv.resample_colwise(MVfp, ncdat.z_t.values, zdbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMVfp_B, 'dMVfp_B', noload=boolLGSnoload)
 #dMVfp_C     = LGS(lambda: utils_conv.resample_colwise(MVfp, ncdat.z_t.values, zdbbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMVfp_C, 'dMVfp_C', noload=boolLGSnoload)

# =======================================================================================
#  Streamfunctions - V method
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - BSF Streamfunction (in Sv)
BSF_mgrd, MVzint = utils_BSF.calc_BSF_mgrd(MV_mgrd, dump_MVzint=True)

# ---------------------------------------------------------------------------------------
# - MOC Streamfunction (in depth space) (in Sv)
#   (1) zonal integration
MVxint_mgrd     = np.nansum(MVp, axis=2)
MVxint_auxgrd   = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, zT_auxgrd, 'V', MVp, ncdat, dir_auxgrd)
#   (2) vertical integration
MOC_mgrd_V      = utils_ana.nancumsum(MVxint_mgrd, axis=0)
MOC_auxgrd_V    = utils_ana.nancumsum(MVxint_auxgrd, axis=0)

# ---------------------------------------------------------------------------------------
# - dMOC Streamfunction (in density space) (in Sv)

# (M0) ANDREAS
# ------------
#   (1) vertical integration of MV (first! as weighting known in depth-space)
MV_zint         = utils_ana.nancumsum(MVp, axis=0)
#   (2) resampling on density axis
dMV_zint        = utils_conv.resample_colwise(MV_zint, ncdat.z_t.values, zdbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force')
 #dMV_zint        = utils_conv.resample_colwise(MV_zint, dens, db, method='dMV_db', fill_value=np.nan, mask = ATLboolmask, mono_method='force')
#   (3) zonal integration
dMOC_mgrd_V_0     = np.nansum(dMV_zint, axis=2)
dMOC_auxgrd_V_0   = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, dbc, 'dV', dMV_zint, ncdat, dir_auxgrd)


# (M1) LEVYN
# ------------
#   (1) integration of MVfp along density axis weighting with ddb or dzdb
ddb3d = utils_conv.exp_k_to_kji(ddb, dzdb.shape[-2], dzdb.shape[-1])
 #dMVfpc_AI   = utils_ana.nancumsum(dMVfp_A*ddb3d, axis=0)
 #dMVfpc_AII  = utils_ana.nancumsum(dMVfp_A*dzdb, axis=0)
 #dMVfpc_BI   = utils_ana.nancumsum(dMVfp_B*ddb3d, axis=0)
dMVfpc_BII  = utils_ana.nancumsum(dMVfp_B*dzdb, axis=0)
 #dMVfpc_CI   = utils_ana.nancumsum(dMVfp_C*ddb3d, axis=0)
 #dMVfpc_CII  = utils_ana.nancumsum(dMVfp_C*dzdb, axis=0)

#   (2) zonal integration of dMVp
dMOC_mgrd_V_BII    = np.nansum(dMVfpc_BII, axis=-1)
dMOC_auxgrd_V_BII  = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, zdbc, 'dV', dMVfpc_BII, ncdat, dir_auxgrd)


# AII, BII, CII all work quite fine
# forget ABCI --> maybe multiplication with something in meters would help with units, not with shape.


