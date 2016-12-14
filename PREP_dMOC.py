# preparation of all dens-conversions


import numpy as np
import gsw
from netCDF4 import Dataset
import xarray as xr

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


# -----------------------------------------------------------------------------
#  - Paths to files
CESMversion = 4
DBsetup     = np.array([[11, 30, 6], [30,35,6], [35, 36.5, 51], [36.5,38, 101], [38,43, 8]]) # setup for density bins (border values) [lower, upper, steps]
# input data
fpath       = paths.get_path2data('lm_1deg', 'anndat')
fnames      = ['b40.lm850-1850.1deg.001.pop.h.{:04d}.ann.4.cdf'.format(i) for i in np.arange(850, 1500)]
# output variables
dir_dens    = paths.get_path2vars('dens', CESMversion=CESMversion, mkdir=True)
dir_DB      = 'DB_%gto%gin%g_%gto%gin%g_%gto%gin%g_%gto%gin%g_%gto%gin%g/' % tuple(DBsetup.flatten())
dir_vars = dir_dens + dir_DB

# -----------------------------------------------------------------------------
#  - Loop over files

for ii in np.arange(255,300):
    
    print ii
    
    # Load data
    ncdat = xr.open_dataset(fpath+fnames[0], decode_times=False)
    
    # Masks for Atlantic
    ATLboolmask = utils_mask.get_ATLbools(ncdat.REGION_MASK.values) # boolean mask
    ATLiter = utils_mask.get_ATLiter(ATLboolmask)
    
# -----------------------------------------------------------------------------
# PART I // DENSITY
# -----------------------------------------------------------------------------
    #  - sigma 2 on T-grid
    SA = ncdat.SALT[0,:,:,:].values             # absolute salinity

    PT = ncdat.TEMP[0,:,:,:].values             # potential temperature
    CT = gsw.CT_from_pt(SA, PT)                 # conservative temperature
    sig2T = gsw.sigma2(SA, CT)                  # potential density anomaly referenced to 2000dbar
    
    # - conversion of sigma 2 on U-grid (using cannonical cumsum method)
    sig2U = np.zeros_like(sig2T)
    foo = utils_ana.canonical_cumsum(sig2T, 2, axis=-1)
    sig2U[:,:-1,:-1] = .25*utils_ana.canonical_cumsum(foo, 2, axis=-2)
    sig2U[:,-1,:-1] = .5*utils_ana.canonical_cumsum(sig2T, 2, axis=-1)[:,-1,:]
    sig2U[:,:-1,-1] = .5*utils_ana.canonical_cumsum(sig2T, 2, axis=-2)[:,:,-1]
    sig2U[:,-1,-1] = sig2T[:,-1,-1]
    del foo
    
    # density bins:  border-values (=DBb), center-values (=DBc) and thickness (=dDB)
    DBb = np.concatenate((np.linspace(DBsetup[0,0],DBsetup[0,1],DBsetup[0,2]), np.linspace(DBsetup[1,0],DBsetup[1,1],DBsetup[1,2]), np.linspace(DBsetup[2,0],DBsetup[2,1],DBsetup[2,2]), np.linspace(DBsetup[3,0],DBsetup[3,1],DBsetup[3,2]), np.linspace(DBsetup[4,0],DBsetup[4,1],DBsetup[4,2])))
    DBc = np.convolve(DBb, [.5,.5])[1:-1] # find midpoints
    dDB = np.diff(DBb)

    # depth of isopycnals calculated as z(DBc) = zDB
    z_t_3d = utils_conv.exp_k_to_kji(ncdat.z_t, sig2U.shape[-2], sig2U.shape[-1])
    zDBc = utils_conv.resample_colwise(z_t_3d, sig2U, DBc, 'lin', fill_value=np.nan, mask = ATLboolmask, mono_method='sort')
    
    # thickness of zDBc (=dzDB) (with first thickness from surf to first bin)
    dzDBc = np.vstack([utils_conv.exp_ji_to_kji(zDBc[0],1), np.diff(zDBc, axis=0)])
    
# -----------------------------------------------------------------------------
# PART II // TRANSPORTS
# -----------------------------------------------------------------------------
    # calculate transports
    MV      = utils_transp.calc_MV(ncdat).values        # = V * DX *DZ
    MVf     = utils_transp.calc_MVflat(ncdat).values    # = V * DX

    # -------------------------------------------------------------------------
    # - MOC Streamfunction (in depth space) (in Sv)
    MVxint_mgrd     = np.nansum(MV, axis=2)
    MOC_mgrd_V      = utils_ana.nancumsum(MVxint_mgrd, axis=0)

    # -------------------------------------------------------------------------
    # - dMOC Streamfunction (in density space) (in Sv) (M1)
    #   (1) integration of MVfp along density axis weighting with dDB or dzDB
    #   (2) projection on auxilary grid
    #   (3) zonal integration of dMVfc(p)
    '''    
    dMVf        = utils_conv.resample_colwise(MVf, ncdat.z_t.values, zDBc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force')
    dMVfc       = utils_ana.nancumsum(dMVf*dzDBc, axis=0) #! changed dDB to dzDBc
    dMVfcp      = utils_conv.project_on_auxgrd(dMVfc, ncdat.ANGLE.values)
    dMOC        = np.nansum(dMVfcp, axis=-1)
    '''

    dMVf        = utils_conv.resample_colwise(MVf, sig2U, DBc, method='dMV', fill_value=np.nan, mask = ATLboolmask, mono_method='force')
    dDB3d       = utils_conv.exp_k_to_kji(dDB, dMVf.shape[-2], dMVf.shape[-1])
    dMVfc       = utils_ana.nancumsum(dMVf*dDB3d, axis=0)
    dMVfcp      = utils_conv.project_on_auxgrd(dMVfc, ncdat.ANGLE.values)
    dMOC        = np.nansum(dMVfcp, axis=-1)

# -----------------------------------------------------------------------------
# PART III // SAVE SOME VARIABLES
# -----------------------------------------------------------------------------
    # Note: the part [:-4] removes ".cdf"

    # create folder with choice of densitybinning
    utils_misc.mkdir(dir_vars)

    utils_misc.savevar(sig2T, dir_vars+'sig2T_'+fnames[ii][:-4], 'sig2T')
    utils_misc.savevar(sig2U, dir_vars+'sig2U_'+fnames[ii][:-4])
    
    #utils_misc.savevar(zDBc, dir_vars+'zDBc_'+fnames[ii][:-4], 'zDBc')
    #utils_misc.savevar(zDBb, dir_vars+'zDBb_'+fnames[ii][:-4], 'zDBb')
    #utils_misc.savevar(zDBbc, dir_vars+'zDBbc_'+fnames[ii][:-4], 'zDBbc')
    utils_misc.savevar(dMOC, dir_vars+'dMOC_'+fnames[ii][:-4], 'dMOC')

# -----------------------------------------------------------------------------
# CLOSE NC-FILE
# -----------------------------------------------------------------------------
    ncdat.close()