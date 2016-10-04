# preparation of all dens-conversions


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
import CESM_utils_BSF as utils_BSF
import CESM_utils_time as utils_time
import CESM_utils_analysis as utils_ana
import CESM_paths as paths
from IPython.core.debugger import Tracer; debug_here = Tracer()


# =============================================================================
#  Get Data
# =============================================================================
# -----------------------------------------------------------------------------
#  - Paths
CESMversion = 4
dbsetup = np.array([[11, 30, 6], [30,35,6], [35, 36.5, 51], [36.5,38, 101], [38,43, 8]]) # setup for density bins [lower, upper, steps]
str_db          = 'db_%gto%gin%g_%gto%gin%g_%gto%gin%g_%gto%gin%g_%gto%gin%g' % tuple(dbsetup.flatten())
dir_dens        = paths.get_path2vars('dens', CESMversion=CESMversion, mkdir=True)

fpath = paths.get_path2data('lm_1deg', 'anndat')
fnames = ['b40.lm850-1850.1deg.001.pop.h.{:04d}.ann.4.cdf'.format(i) for i in np.arange(850, 1500)]

# =======================================================================================
#  Settings and Directories
# =======================================================================================

for ii in np.arange(850,1500):
    #  - Paths
    path_sig2       = dir_dens+'sig2_'+fnames[ii][:-4]
    path_densU      = dir_dens+'densU_'+fnames[ii][:-4]
    path_zdbc       = dir_dens+'/{}/'.format(str_db)+'zdbc_'+fnames[ii][:-4]
    path_zdbb       = dir_dens+'/{}/'.format(str_db)+'zdbb_'+fnames[ii][:-4]
    path_zdbbc      = dir_dens+'/{}/'.format(str_db)+'zdbbc_'+fnames[ii][:-4]
    
    #  - Load data
    ncdat = xr.open_dataset(fpath+fnames[0], decode_times=False)

    #  - Density Conversion
    SA = ncdat.SALT[0,:,:,:].values             # absolute salinity
    PT = ncdat.TEMP[0,:,:,:].values             # potential temperature
    CT = gsw.CT_from_pt(SA, PT)                 # conservative temperature
    sig2 = gsw.sigma2(SA, CT)                   # potential density anomaly referenced to 2000dbar
    
    # - conversion T-->U
    densU = np.zeros_like(sig2)
    foo1 = utils_ana.canonical_cumsum(sig2, 2, axis=-1)
    densU[:,:-1,:-1] = .25*utils_ana.canonical_cumsum(foo1, 2, axis=-2)
    densU[:,-1,:-1] = .5*utils_ana.canonical_cumsum(sig2, 2, axis=-1)[:,-1,:]
    densU[:,:-1,-1] = .5*utils_ana.canonical_cumsum(sig2, 2, axis=-2)[:,:,-1]
    densU[:,-1,-1] = sig2[:,-1,-1]
    
    # density bins:  border-values (=dbb), center-values (=dbc) and thickness (=ddb)
    dbb = np.concatenate((np.linspace(dbsetup[0,0],dbsetup[0,1],dbsetup[0,2]), np.linspace(dbsetup[1,0],dbsetup[1,1],dbsetup[1,2]), np.linspace(dbsetup[2,0],dbsetup[2,1],dbsetup[2,2]), np.linspace(dbsetup[3,0],dbsetup[3,1],dbsetup[3,2]), np.linspace(dbsetup[4,0],dbsetup[4,1],dbsetup[4,2])))
    dbc = np.convolve(dbb, np.array([.5,.5]))[1:-1]
    
    # depth of isopycnals (zdbbc) calculated as z(dbc) (=zdbc) and as c(zdbb) (=zdbbc)
    z_t_3d = utils_conv.exp_k_to_kji(ncdat.z_t, densU.shape[-2], densU.shape[-1])
    zdbc = utils_conv.resample_colwise(z_t_3d, densU, dbc, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort')
    zdbb = utils_conv.resample_colwise(z_t_3d, densU, dbb, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort')
    zdbbc = np.ones_like(zdbc) * np.nan
    for j in np.arange(zdbb.shape[-2]):
        for i in np.arange(zdbb.shape[-1]):
            zdbbc[:,j,i] = np.convolve(zdbb[:,j,i], np.array([.5,.5]))[1:-1] # centre-values
    
    # save variables
    utils_misc.savevar(sig2, path_sig2)
    utils_misc.savevar(densU, path_densU)
    utils_misc.savevar(zdbc, path_zdbc)
    utils_misc.savevar(zdbb, path_zdbb)
    utils_misc.savevar(zdbbc, path_zdbbc)
