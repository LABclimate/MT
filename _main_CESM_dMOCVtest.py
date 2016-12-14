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
import CESM_utils_BSF as utils_BSF
import CESM_utils_time as utils_time
import CESM_paths as paths
from IPython.core.debugger import Tracer; debug_here = Tracer()

from UTILS_misc import LGS, GS, LG

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
auxgrd_name = ['lat395model', 'lat198model', 'lat99model', 'lat170eq80S90N', 'lat340eq80S90N'][2] # choose aux grid
 #1 DBsetup = np.array([[20, 33, 14], [33.1, 36, 11], [36.05, 37.5, 11], [38, 43, 11]]) # setup for density bins [lower, upper, steps]
 #2 DBsetup = np.array([[12, 22, 3], [25,30, 3], [31,35,6], [35, 38, 40], [38.5,42, 3]]) # setup for density bins [lower, upper, steps]
 #3 DBsetup = np.array([[11, 30, 6], [30,35,6], [35, 36.5, 51], [36.5,37, 51], [37,43, 7]]) # setup for density bins [lower, upper, steps]
DBsetup     = np.array([[11, 30, 6], [30,35,6], [35, 36.5, 51], [36.5,38, 101], [38,43, 8]]) # setup for density bins (border values) [lower, upper, steps]
boolLGSnoload = True    # boolean | False: normal loadgetsave procedure | True: noload-mode > will get and save (overwrite!)
# ---------------------------------------------------------------------------------------
# automatical generation of directory names
dir_mgrd        = paths.get_path2vars('mgrd',CESMversion=CESMversion, mkdir=True)
dir_dens        = paths.get_path2vars('dens', CESMversion=CESMversion, mkdir=True)
dir_DB          = 'DB_%gto%gin%g_%gto%gin%g_%gto%gin%g_%gto%gin%g_%gto%gin%g/' % tuple(DBsetup.flatten())
dir_auxgrd      = paths.get_path2vars(auxgrd_name, CESMversion=CESMversion, mkdir=True)
utils_misc.mkdir(dir_dens+dir_DB)
# paths (directory + filenames) to temp variables
path_dMV       = dir_dens+dir_DB+'dMV_'
path_dMVf      = dir_dens+dir_DB+'dMVf_'
path_zDBb       = dir_dens+dir_DB+'zDBb'
path_zDBc       = dir_dens+dir_DB+'zDBc'
path_zDBbc      = dir_dens+dir_DB+'zDBbc'
path_HT_auxgrd_xmax = dir_auxgrd+'HT_auxgrd_xmax'
path_HT_mgrd_xmax   = dir_mgrd+'HT_mgrd_xmax'
path_HU_auxgrd_xmax = dir_auxgrd+'HU_auxgrd_xmax'
path_HU_mgrd_xmax   = dir_mgrd+'HU_mgrd_xmax'
path_lat_auxgrd     = '../variables/CESM_gen/lat_auxgrd_'+auxgrd_name
path_fraction_mask = dir_auxgrd+'fraction_mask'

# =======================================================================================
#  Transformation on different grids (Density, Spatial auxiliary grid)
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Mask for Atlantic
ATLboolmask = utils_mask.get_ATLbools(ncdat.REGION_MASK.values) # boolean mask
ATLiter = utils_mask.get_ATLiter(ATLboolmask)

# ---------------------------------------------------------------------------------------
# - Spatial auxiliary grid

lat_auxgrd = LGS(lambda: utils_mask.gen_auxgrd(ncdat4, auxgrd_name), path_lat_auxgrd, 'lat_auxgrd')
lat_mgrd = ncdat.TLAT.isel(nlon=0)          # mean of LAT for each j #! very inappropriate
fraction_mask = LGS(lambda: utils_mask.gen_fraction_mask(lat_auxgrd, ncdat), path_fraction_mask, 'fraction_mask')
 # fraction_mask, stats = utils_mask.gen_fraction_mask(lat_auxgrd, ncdat)
 # masks_auxgrd = dict()
 # masks_auxgrd['overlay_lat'] = LGS(lambda: utils_mask.gen_mask_auxgrd_overlay_lat(lat_auxgrd, ncdat), path_mask_auxgrd_grd_overlay_lat, 'mask_auxgrd_overlay_lat', noload=False)
 # masks_auxgrd['iter_maskcombo'] = LGS(lambda: utils_mask.gen_iter_maskcombo(lat_auxgrd, ncdat, masks_auxgrd['overlay_lat']), path_mask_auxgrd_iter_maskcombo, 'iter_maskcombo', noload=False)

# ---------------------------------------------------------------------------------------
# - Density grid/bins
#! note that for volume representation T-grid DBc is needed!!!!
SA = ncdat.SALT[0,:,:,:].values             # absolute salinity
PT = ncdat.TEMP[0,:,:,:].values             # potential temperature
CT = gsw.CT_from_pt(SA, PT)                 # conservative temperature
sig2T = gsw.sigma2(SA, CT)                   # potential density anomaly referenced to 2000dbar
RHO = ncdat.RHO[0,:,:,:].values*1000-1000   # in-situ density anomaly [SI]

# - conversion T-->U
sig2U = np.zeros_like(sig2T)
foo1 = utils_ana.canonical_cumsum(sig2T, 2, axis=-1)
sig2U[:,:-1,:-1] = .25*utils_ana.canonical_cumsum(foo1, 2, axis=-2)
sig2U[:,-1,:-1] = .5*utils_ana.canonical_cumsum(sig2T, 2, axis=-1)[:,-1,:]
sig2U[:,:-1,-1] = .5*utils_ana.canonical_cumsum(sig2T, 2, axis=-2)[:,:,-1]
sig2U[:,-1,-1] = sig2T[:,-1,-1]

# density bins:  border-values (=DBb), center-values (=DBc) and thickness (=dDB)
DBb = np.concatenate((np.linspace(DBsetup[0,0],DBsetup[0,1],DBsetup[0,2]), np.linspace(DBsetup[1,0],DBsetup[1,1],DBsetup[1,2]), np.linspace(DBsetup[2,0],DBsetup[2,1],DBsetup[2,2]), np.linspace(DBsetup[3,0],DBsetup[3,1],DBsetup[3,2]), np.linspace(DBsetup[4,0],DBsetup[4,1],DBsetup[4,2])))
DBc = np.convolve(DBb, [.5,.5])[1:-1]  # find midpoints
# depth of isopycnals (zDBbc) calculated as z(DBc) (=zDBc) and as c(zDBb) (=zDBbc)
z_t_3d = utils_conv.exp_k_to_kji(ncdat.z_t, sig2U.shape[-2], sig2U.shape[-1])

 # zDBc = LGS(lambda: utils_conv.resample_colwise(z_t_3d, sig2U, DBc, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort'), path_zDBc, 'zDBc', noload=False)
 # zDBb = LGS(lambda: utils_conv.resample_colwise(z_t_3d, sig2U, DBb, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort'), path_zDBb, 'zDBb', noload=False)
 # zDBbc = np.ones_like(zDBc) * np.nan
 # for j in np.arange(zDBb.shape[-2]):
 #     for i in np.arange(zDBb.shape[-1]):
 #         zDBbc[:,j,i] = np.convolve(zDBb[:,j,i], [.5,.5])[1:-1] # centre-values

# thickness of DB (=dDB) and zDBb (=dzDB)
dDB = np.diff(DBb)
#dzDB = np.diff(zDBb, axis=0)
del PT, CT, z_t_3d, DBsetup

# ---------------------------------------------------------------------------------------
# Volume representation for axis-scaling (fully on T-grid)
dz3d = utils_conv.exp_k_to_kji(ncdat.dz, sig2T.shape[-2], sig2T.shape[-1])    # in cgs
TAREA3d = utils_conv.exp_ji_to_kji(ncdat.TAREA, sig2T.shape[0])              # in cgs
vol3d = dz3d*TAREA3d                                                        # in cgs
inds = np.digitize(sig2T, DBc)
# pre-allocation
VDB_glob = np.zeros(shape=[len(DBc)])
VDB_reg  = np.zeros(shape=[len(DBc)])
VDB_col  = np.zeros(shape=[len(DBc), sig2T.shape[-2], sig2T.shape[-1]])
# summing up volumina over different basins
for b in np.arange(len(DBc)):
    VDB_glob[b] = np.sum(vol3d[inds==b])                # global        # in cgs
    VDB_reg[b]  = np.sum((vol3d*ATLboolmask)[inds==b])  # regional      # in cgs
    VDB_col[b]  = np.sum(vol3d[inds==b], axis=0)        # column-wise   # in cgs
# axes and ticks for plots
ax_vol_glob = np.cumsum(VDB_glob) - VDB_glob/2
ax_vol_reg = np.cumsum(VDB_reg) - VDB_reg/2
 #1 ticks_dens = [20.5, 28.5, 35.275, 36.025, 36.2675, 37.1375, 37.75, 38.75]
 #2 ticks_dens = [33, 35, 36, DBc[35], DBc[38], DBc[40], DBc[41]]
 #3 ticks_dens = [33.5, 35, 36.005, DBc[62], DBc[100], DBc[110], DBc[113]]
ticks_dens = [33.5, 35, 36.005, DBc[62], DBc[100], DBc[110]]
ticks_vol_glob = ax_vol_glob[np.in1d(DBc, ticks_dens)]
ticks_vol_reg = ax_vol_reg[np.in1d(DBc, ticks_dens)]
ticks_dens_rd = np.round(ticks_dens, 2)

del dz3d, TAREA3d, vol3d, inds

# =======================================================================================
#  Variables contained in model output
# =======================================================================================
# - temperature -------------------------------------------------------------------------
T = ncdat.TEMP.mean(dim='time')
T = utils_mask.mask_ATLANTIC(T, ncdat.REGION_MASK)
#T_dens = utils_conv.resample_colwise(T.values, sig2U, DBc, method='lin', fill_value=np.nan, mask=ATLboolmask, mono_method='sort')
# - Zonal maxima of ocean depth --------------------------------------------------------

#! #TODO REWRITE THESE FUNCTIONS
HT_auxgrd_xmax  = LGS(lambda: utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'T', fraction_mask), path_HT_auxgrd_xmax, 'HT_auxgrd_xmax')
HU_auxgrd_xmax  = LGS(lambda: utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'U', fraction_mask), path_HU_auxgrd_xmax, 'HU_auxgrd_xmax')
HT_mgrd_xmax    = LGS(lambda: utils_mask.calc_H_mgrd_xmax(ncdat, 'T'), path_HT_mgrd_xmax, 'HT_mgrd_xmax')
HU_mgrd_xmax    = LGS(lambda: utils_mask.calc_H_mgrd_xmax(ncdat, 'U'), path_HU_mgrd_xmax, 'HU_mgrd_xmax')

# =======================================================================================
#  Volume transports (in Sv)
# =======================================================================================
# - MV on model grid
MV      = utils_transp.calc_MV(ncdat).values        # = V * DX *DZ
MVf     = utils_transp.calc_MVflat(ncdat).values    # = V * DX
# - mutations
'''
MV = np.ones_like(MV)                                                  # mutation
MVf = np.ones_like(MV)                                                # mutation
MV[:,190:300,-65:-45] = np.nan*np.ones_like(MV[:,190:300,-65:-45])      # mask
MVf[:,190:300,-65:-45] = np.nan*np.ones_like(MVf[:,190:300,-65:-45])    # mask
'''

# ---------------------------------------------------------------------------------------
# - resampled on density axis
 #dMV       = LGS(lambda: utils_conv.resample_colwise(MV, ncdat.z_t.values, zDBc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMV, 'dMV', noload=boolLGSnoload)
 #dMVf      = LGS(lambda: utils_conv.resample_colwise(MVf.values, ncdat.z_t.values, zDBc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMVf, 'dMVf', noload=boolLGSnoload) #B
 dMVf      = utils_conv.resample_colwise(MVf, sig2U, DBc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force')
 #dMVf      = LGS(lambda: utils_conv.resample_colwise(MVf.values, ncdat.z_t.values, zDBbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMVfp, 'dMVfp', noload=boolLGSnoload) #C
 #dMVp      = utils_conv.project_on_auxgrd(dMV, ncdat.ANGLE.values)
 #dMVfp     = utils_conv.project_on_auxgrd(dMVf, ncdat.ANGLE.values)

# =======================================================================================
#  Streamfunctions - V method
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - BSF Streamfunction (in Sv)
'''
BSF_mgrd, MVzint = utils_BSF.calc_BSF_mgrd(MV, dump_MVzint=True)
'''

# ---------------------------------------------------------------------------------------
# - MOC Streamfunction (in depth space) (in Sv)
#   (1) projection on auxilary grid
MVp         = utils_conv.project_on_auxgrd(MV, ncdat.ANGLE.values)
#   (2) zonal integration
MVxint_mgrd     = np.nansum(MV, axis=2)
MVxint_auxgrd   = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, ncdat.z_t.values, MVp, fraction_mask, ATLiter)
#   (3) vertical integration
MOC_mgrd_V      = utils_ana.nancumsum(MVxint_mgrd, axis=0)
MOC_auxgrd_V    = utils_ana.nancumsum(MVxint_auxgrd, axis=0)

# ---------------------------------------------------------------------------------------
# - dMOC Streamfunction (in density space) (in Sv)

# (M0)
# -----
#   (1) vertical integration of MV (first! as weighting known in depth-space)
MVc         = utils_ana.nancumsum(MV, axis=0)
#   (2) resampling on density axis
dMVc        = utils_conv.resample_colwise(MVc, ncdat.z_t.values, dDB, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force') # II
 #dMVc        = utils_conv.resample_colwise(MVc, sig2U, DB, method='dMV_DB', fill_value=np.nan, mask = ATLboolmask, mono_method='force') # I
#   (3) projection on auxilary grid
dMVcp       = utils_conv.project_on_auxgrd(dMVc, ncdat.ANGLE.values)
#   (4) zonal integration
dMOC_mgrd_V_0   = np.nansum(dMVc, axis=2)
dMOC_auxgrd_V_0 = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, DBc, dMVcp, fraction_mask, ATLiter)


# (M1)dDB
# -----
#   (1) integration of MVfp along density axis weighting with dDB
dDB3d       = utils_conv.exp_k_to_kji(dDB, dMVf.shape[-2], dMVf.shape[-1])
dMVfc       = utils_ana.nancumsum(dMVf*dDB3d, axis=0)
#   (2) projection on auxilary grid
dMVfcp      = utils_conv.project_on_auxgrd(dMVfc, ncdat.ANGLE.values)
#   (3) zonal integration of dMVfc(p)
dMOC_mgrd_V_B     = np.nansum(dMVfc, axis=-1)
dMOC_auxgrd_V_B   = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, dDB, dMVfcp, fraction_mask, ATLiter)
