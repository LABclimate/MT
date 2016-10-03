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
auxgrd_name = ['lat395model', 'lat198model', 'lat99model', 'lat170eq80S90N', 'lat340eq80S90N'][2] # choose aux grid
 #1 dbsetup = np.array([[20, 33, 14], [33.1, 36, 11], [36.05, 37.5, 11], [38, 43, 11]]) # setup for density bins [lower, upper, steps]
 #2 dbsetup = np.array([[12, 22, 3], [25,30, 3], [31,35,6], [35, 38, 40], [38.5,42, 3]]) # setup for density bins [lower, upper, steps]
 #3 dbsetup = np.array([[11, 30, 6], [30,35,6], [35, 36.5, 51], [36.5,37, 51], [37,43, 7]]) # setup for density bins [lower, upper, steps]
dbsetup = np.array([[11, 30, 6], [30,35,6], [35, 36.5, 51], [36.5,38, 101], [38,43, 8]]) # setup for density bins [lower, upper, steps]
densChoice = 'sig2'     # string | choice of density to use for resampling (either RHO or sig2) 
boolLGSnoload = True  # boolean | False: normal loadgetsave procedure | True: noload-mode > will get and save (overwrite!)
# ---------------------------------------------------------------------------------------
# automatical generation of directory names
str_db          = 'db_%gto%gin%g_%gto%gin%g_%gto%gin%g_%gto%gin%g_%gto%gin%g' % tuple(dbsetup.flatten())
str2_db         = 'db %g-%g(%g) %g-%g(%g) %g-%g(%g) %g-%g(%g) %g-%g(%g)' % tuple(dbsetup.flatten())
dir_mgrd        = paths.get_path2vars('mgrd',CESMversion=CESMversion, mkdir=True)
dir_dens        = paths.get_path2vars('dens', CESMversion=CESMversion, mkdir=True)
dir_auxgrd      = paths.get_path2vars(auxgrd_name, CESMversion=CESMversion, mkdir=True)
utils_misc.mkdir(dir_dens+'/{}/'.format(str_db))
# paths (directory + filenames) to temp variables
path_dMV       = dir_dens+'/{}/'.format(str_db)+'dMV_'+densChoice
path_dMVf      = dir_dens+'/{}/'.format(str_db)+'dMVf_'+densChoice
path_zdbb       = dir_dens+'/{}/'.format(str_db)+'zdbb'
path_zdbc       = dir_dens+'/{}/'.format(str_db)+'zdbc'
path_zdbbc      = dir_dens+'/{}/'.format(str_db)+'zdbbc'
path_HT_auxgrd_xmax = dir_auxgrd+'HT_auxgrd_xmax'
path_HT_mgrd_xmax   = dir_mgrd+'HT_mgrd_xmax'
path_HU_auxgrd_xmax = dir_auxgrd+'HU_auxgrd_xmax'
path_HU_mgrd_xmax   = dir_mgrd+'HU_mgrd_xmax'
path_lat_auxgrd     = '../variables/CESM_gen/lat_auxgrd_'+auxgrd_name
path_fraction_mask = dir_auxgrd+'fraction_mask'
path_mask_auxgrd_grd_overlay_lat = dir_auxgrd+'mask_auxgrd_overlay_lat'
path_mask_auxgrd_iter_maskcombo = dir_auxgrd+'iter_maskcombo'
# =======================================================================================
#  Transformation on different grids (Density, Spatial auxiliary grid)
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - Mask for Atlantic
ATLboolmask = utils_mask.get_ATLbools(ncdat.REGION_MASK.values) # boolean mask
ATLiter = utils_mask.get_ATLiter(ATLboolmask)

# ---------------------------------------------------------------------------------------
# - Spatial auxiliary grid
lat_auxgrd = LGS(lambda: utils_mask.gen_auxgrd(ncdat4, auxgrd_name), path_lat_auxgrd, 'lat_auxgrd', noload=True)
lat_mgrd = ncdat.TLAT.isel(nlon=0)          # mean of LAT for each j #! very inappropriate
#masks_auxgrd = dict()
#masks_auxgrd['overlay_lat'] = LGS(lambda: utils_mask.gen_mask_auxgrd_overlay_lat(lat_auxgrd, ncdat), path_mask_auxgrd_grd_overlay_lat, 'mask_auxgrd_overlay_lat', noload=False)
fraction_mask = LGS(lambda: utils_mask.gen_fraction_mask(lat_auxgrd, ncdat), path_fraction_mask, 'fraction_mask', noload=True)
fraction_mask, diffAREA = utils_mask.gen_fraction_mask(lat_auxgrd, ncdat)
#masks_auxgrd['iter_maskcombo'] = LGS(lambda: utils_mask.gen_iter_maskcombo(lat_auxgrd, ncdat, masks_auxgrd['overlay_lat']), path_mask_auxgrd_iter_maskcombo, 'iter_maskcombo', noload=False)
# ---------------------------------------------------------------------------------------
# - Density grid/bins
#! note that for volume representation T-grid dbc is needed!!!!
SA = ncdat.SALT[0,:,:,:].values             # absolute salinity
PT = ncdat.TEMP[0,:,:,:].values             # potential temperature
CT = gsw.CT_from_pt(SA, PT)                 # conservative temperature
sig2 = gsw.sigma2(SA, CT)                   # potential density anomaly referenced to 2000dbar
RHO = ncdat.RHO[0,:,:,:].values*1000-1000   # in-situ density anomaly [SI]
if densChoice == 'sig2': densT = sig2
elif densChoice == 'rho': densT = RHO

# - conversion T-->U
densU = np.zeros_like(densT)
foo1 = utils_ana.canonical_cumsum(densT, 2, axis=-1)
densU[:,:-1,:-1] = .25*utils_ana.canonical_cumsum(foo1, 2, axis=-2)
densU[:,-1,:-1] = .5*utils_ana.canonical_cumsum(densT, 2, axis=-1)[:,-1,:]
densU[:,:-1,-1] = .5*utils_ana.canonical_cumsum(densT, 2, axis=-2)[:,:,-1]
densU[:,-1,-1] = densT[:,-1,-1]

# density bins:  border-values (=dbb), center-values (=dbc) and thickness (=ddb)
dbb = np.concatenate((np.linspace(dbsetup[0,0],dbsetup[0,1],dbsetup[0,2]), np.linspace(dbsetup[1,0],dbsetup[1,1],dbsetup[1,2]), np.linspace(dbsetup[2,0],dbsetup[2,1],dbsetup[2,2]), np.linspace(dbsetup[3,0],dbsetup[3,1],dbsetup[3,2]), np.linspace(dbsetup[4,0],dbsetup[4,1],dbsetup[4,2])))
dbc = np.convolve(dbb, np.array([.5,.5]))[1:-1]
# depth of isopycnals (zdbbc) calculated as z(dbc) (=zdbc) and as c(zdbb) (=zdbbc)
z_t_3d = utils_conv.exp_k_to_kji(ncdat.z_t, densU.shape[-2], densU.shape[-1])
zdbc = LGS(lambda: utils_conv.resample_colwise(z_t_3d, densU, dbc, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort'), path_zdbc, 'zdbc', noload=False)
zdbb = LGS(lambda: utils_conv.resample_colwise(z_t_3d, densU, dbb, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort'), path_zdbb, 'zdbb', noload=False)
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
dz3d = utils_conv.exp_k_to_kji(ncdat.dz, densT.shape[-2], densT.shape[-1])    # in cgs
TAREA3d = utils_conv.exp_ji_to_kji(ncdat.TAREA, densT.shape[0])              # in cgs
vol3d = dz3d*TAREA3d                                                        # in cgs
inds = np.digitize(densT, dbc)
# pre-allocation
Vdb_glob = np.zeros(shape=[len(dbc)])
Vdb_reg  = np.zeros(shape=[len(dbc)])
Vdb_col  = np.zeros(shape=[len(dbc), densT.shape[-2], densT.shape[-1]])
# summing up volumina over different basins
for b in np.arange(len(dbc)):
    Vdb_glob[b] = np.sum(vol3d[inds==b])                # global        # in cgs
    Vdb_reg[b]  = np.sum((vol3d*ATLboolmask)[inds==b])  # regional      # in cgs
    Vdb_col[b]  = np.sum(vol3d[inds==b], axis=0)        # column-wise   # in cgs
# axes and ticks for plots
ax_vol_glob = np.cumsum(Vdb_glob) - Vdb_glob/2
ax_vol_reg = np.cumsum(Vdb_reg) - Vdb_reg/2
 #1 ticks_dens = [20.5, 28.5, 35.275, 36.025, 36.2675, 37.1375, 37.75, 38.75]
 #2 ticks_dens = [33, 35, 36, dbc[35], dbc[38], dbc[40], dbc[41]]
 #3 ticks_dens = [33.5, 35, 36.005, dbc[62], dbc[100], dbc[110], dbc[113]]
ticks_dens = [33.5, 35, 36.005, dbc[62], dbc[100], dbc[110]]
ticks_vol_glob = ax_vol_glob[np.in1d(dbc, ticks_dens)]
ticks_vol_reg = ax_vol_reg[np.in1d(dbc, ticks_dens)]
ticks_dens_rd = np.round(ticks_dens, 2)

del dz3d, TAREA3d, vol3d, inds

# =======================================================================================
#  Variables contained in model output
# =======================================================================================
# - temperature -------------------------------------------------------------------------
T = ncdat.TEMP.mean(dim='time')
T = utils_mask.mask_ATLANTIC(T, ncdat.REGION_MASK)
#T_dens = utils_conv.resample_colwise(T.values, densU, dbc, method='lin', fill_value=np.nan, mask=ATLboolmask, mono_method='sort')
# - Zonal maxima of ocean depth --------------------------------------------------------
HT_auxgrd_xmax  = LGS(lambda: utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'T', masks_auxgrd['iter_maskcombo']), path_HT_auxgrd_xmax, 'HT_auxgrd_xmax', noload=False)
HT_mgrd_xmax    = LGS(lambda: utils_mask.calc_H_mgrd_xmax(ncdat, 'T'), path_HT_mgrd_xmax, 'HT_mgrd_xmax', noload=False)
HU_auxgrd_xmax  = LGS(lambda: utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'U', masks_auxgrd['iter_maskcombo']), path_HU_auxgrd_xmax, 'HU_auxgrd_xmax', noload=False)
HU_mgrd_xmax    = LGS(lambda: utils_mask.calc_H_mgrd_xmax(ncdat, 'U'), path_HU_mgrd_xmax, 'HU_mgrd_xmax', noload=False)

# =======================================================================================
#  Volume transports (in Sv)
# =======================================================================================
# - MV on model grid
MV      = utils_transp.calc_MV(ncdat).values        # = V * DX *DZ
MVf     = utils_transp.calc_MVflat(ncdat).values    # = V * DX
# - mutations
MV = np.ones_like(MV)                                                  # mutation
MVf = np.ones_like(MV)                                                # mutation
MV[:,190:300,-65:-45] = np.nan*np.ones_like(MV[:,190:300,-65:-45])      # mask
MVf[:,190:300,-65:-45] = np.nan*np.ones_like(MVf[:,190:300,-65:-45])    # mask

# ---------------------------------------------------------------------------------------
# - resampled on density axis
 #dMV         = LGS(lambda: utils_conv.resample_colwise(MV, ncdat.z_t.values, zdbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMV, 'dMV', noload=boolLGSnoload)
 #dMVf        = LGS(lambda: utils_conv.resample_colwise(MVf, ncdat.z_t.values, zdbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMVf, 'dMVf', noload=boolLGSnoload) #B
 #dMVf     = LGS(lambda: utils_conv.resample_colwise(MVf.values, densU, dbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMVfp, 'dMVfp', noload=boolLGSnoload) #A
 #dMVf     = LGS(lambda: utils_conv.resample_colwise(MVf.values, ncdat.z_t.values, zdbbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force'), path_dMVfp, 'dMVfp', noload=boolLGSnoload) #C
 #dMVp        = utils_conv.project_on_auxgrd(dMV, ncdat.ANGLE.values)
 #dMVfp       = utils_conv.project_on_auxgrd(dMVf, ncdat.ANGLE.values)

# =======================================================================================
#  Streamfunctions - V method
# =======================================================================================
# ---------------------------------------------------------------------------------------
# - BSF Streamfunction (in Sv)
BSF_mgrd, MVzint = utils_BSF.calc_BSF_mgrd(MV, dump_MVzint=True)

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
dMVc        = utils_conv.resample_colwise(MVc, ncdat.z_t.values, zdbc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force') # II
 #dMVc        = utils_conv.resample_colwise(MVc, densU, db, method='dMV_db', fill_value=np.nan, mask = ATLboolmask, mono_method='force') # I
#   (3) projection on auxilary grid
dMVcp       = utils_conv.project_on_auxgrd(dMVc, ncdat.ANGLE.values)
#   (4) zonal integration
dMOC_mgrd_V_0   = np.nansum(dMVc, axis=2)
dMOC_auxgrd_V_0 = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, dbc, dMVcp, fraction_mask, ATLiter)


# (M1)
# -----
#   (1) integration of MVfp along density axis weighting with ddb or dzdb
ddb3d       = utils_conv.exp_k_to_kji(ddb, dzdb.shape[-2], dzdb.shape[-1])
dMVfc       = utils_ana.nancumsum(dMVf*dzdb, axis=0)
#   (2) projection on auxilary grid
dMVfcp      = utils_conv.project_on_auxgrd(dMVfc, ncdat.ANGLE.values)
#   (3) zonal integration of dMVfc(p)
dMOC_mgrd_V_B     = np.nansum(dMVfc, axis=-1)
dMOC_auxgrd_V_B   = utils_MOC.calc_Mxint_auxgrd(lat_auxgrd, zdbc, 'dV', dMVfcp, masks_auxgrd['overlay_lat'], masks_auxgrd['iter_maskcombo'])
