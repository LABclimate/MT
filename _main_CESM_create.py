###############################################################################
#                   The CESM python toolbox at KUP
#                   ----- Creation - Toolbox -----
#                               12/12/16
###############################################################################
# ABOUT:
#   > What?:   Creation of basic variables that exist in different variants
#   > How?:    First make settings in PATHS-section to set up run and CESM-run.
#               Then execute any of the blocks below directly in the terminal.
#               These blocks are not intended to be called from any other script.
#               Rather create variables with their variants in here. 
#               Within other scripts variables can be load using their unique code.
#   > Why?:    Tracking of the settings using long filenames (combination of codes)
#               containings all information about how a variable was produced.
#               /data/ folder shall only be called from this script. All other 
#               script shall use data stored to /variables/
###############################################################################
# CONTENT:
#   > anyNCDAT
#   > Gridding parameters
#   > Stack ncdata in time
#   > MOCidx
#   > BSFidx
#   > corrMOCidxBSFidx      and p-value
#   > corrMOCBSFidx         and p-value
#   > corrBSFMOCidx         and p-value
###############################################################################

'''
cd /alphadata02/buerki/no_backup/CESM_new/scripts/
Py2
ipython
'''
import numpy as np
import xarray as xr
import UTILS_misc as utils_misc
import CESM_paths as paths
import CESM_utils_time as utils_time
import CESM_utils_mask as utils_mask
import CESM_utils_analysis as utils_ana

# ==========================================================================================================================================================
# PATHS
# ==========================================================================================================================================================
# -----------------------------------------------------------------------------
# Import paths for data
CESMversion = 4
CESMrun = 'lm_1deg'
fpath = paths.get_path2data('lm_1deg', 'anndat')
fnames = ['b40.lm850-1850.1deg.001.pop.h.{:04d}.ann.4.cdf'.format(i) for i in np.arange(850, 1851)]
# -----------------------------------------------------------------------------
# Dumping paths for variables
# basic path for variables - please choose correct CESMversion, CESMrun above
path_var = '../variables/CESM{}_{}/.'.format(CESMversion, CESMrun)


# ==========================================================================================================================================================
#  Gridding parameters
# ==========================================================================================================================================================
# -----------------------------------------------------------------------------
# alternative paths for local use:
ncdat = utils_misc.loadxvar(path_var+'/any_ncdat', None)                        # get ncdata

folder_grid = '/grid'
utils_misc.savexvar(ncdat.REGION_MASK, path_var+folder_grid+'/REGION_MASK', 'REGION_MASK')
utils_misc.savexvar(ncdat.TLAT, path_var+folder_grid+'/TLAT', 'TLAT')
utils_misc.savexvar(ncdat.TLONG, path_var+folder_grid+'/TLONG', 'TLONG')
utils_misc.savexvar(ncdat.lat_aux_grid, path_var+folder_grid+'/lat_auxgrd', 'lat')
utils_misc.savexvar(ncdat.z_t, path_var+folder_grid+'/zT', 'zT')

# ==========================================================================================================================================================
#  Stack ncdata in time (e.g. BSF, MOC, ...)
# ==========================================================================================================================================================
''' Q: good choice of transport_reg=1 #??
'''
BSF = []
MOC = []
err = []
for idx, fname in enumerate(fnames):
    utils_misc.Counter(idx+1, len(fnames), 'modelyears stacked')
    try:
        ncdat = xr.open_dataset(fpath+fname, decode_times=False)
        MOC = utils_time.concat(MOC, ncdat.MOC.isel(transport_reg=1, moc_comp=0))
        BSF = utils_time.concat(BSF, utils_mask.mask_ATLANTIC(ncdat.BSF, ncdat.REGION_MASK))
        ncdat.close()
    except:
        err = np.append(err, idx)
        print('ERROR with file ' + fpath+fname + ' (idx = {:d})'.format(idx))

#MOC[66,:,:] = np.nan*np.ones_like(MOC[66,:,:])                          # bug-fix due to crap data

utils_misc.savexvar(MOC, path_var+'/MOC/MOC__mod_850to1499', 'MOC')
utils_misc.savexvar(BSF, path_var+'/BSF/BSF__mod_850to1499', 'BSF')


# ==========================================================================================================================================================
# Streamfunction Indices
# ==========================================================================================================================================================
ncdat = utils_misc.loadxvar(path_var+'/any_ncdat',None)                         # ncdat for gridding parameters
# -----------------------------------------------------------------------------
# MOCidx
code_MOC = 'MOC__mod_850to1499'
MOC = utils_misc.loadxvar(path_var+'/MOC/'+code_MOC, 'MOC')                      # load MOC

code_MOCidx = 'MOCidx__{}__max_fullATL_20downwards'.format(code_MOC)           # CODE 1 >>>>> maximum of full atlantic excluding the 20 surface-boxes
MOCidx = MOC[:,20:,:].max(axis=tuple([1,2]))                                    # maximal value omitting surface
utils_misc.savexvar(MOCidx, path_var+'/MOCidx/'+code_MOCidx, 'MOCidx')          # save MOCidx

# -----------------------------------------------------------------------------
# BSFidx
code_BSF = 'BSF__mod_850to1499'
BSF = utils_misc.loadxvar(path_var+'/BSF/'+code_BSF, 'BSF')                     # load BSF

#               ------------------------------
code_BSFidx = 'BSFidx__{}__wmean_fullSPG_45to70_rmse8'.format(code_BSF) # CODE 1 >>>>> weighted mean over full SPG from 45N to 70N and restricted to regmask<=8 (cut at GSR)
selreg = (BSF.TLAT>=45) & (BSF.TLAT<=70) & (ncdat.REGION_MASK<=8)               # select region
BSFsel = BSF.where(selreg)                                                          # --> apply to BSF
TAREAsel = ncdat.TAREA.where(selreg)                                                # --> apply to TAREA
BSFidx = (BSFsel*TAREAsel).mean(axis=tuple([1,2]))/ TAREAsel.mean()             # weighted average over BSFsel
utils_misc.savexvar(BSFidx, path_var+'/BSFidx/'+code_BSFidx, 'BSFidx')          # save BSFidx

code_BSFidx = 'BSFidx__{}__max_fullSPG_45to70_rmse8'.format(code_BSF)    # CODE 2 >>>>> maximum of full SPG from 45N to 70N and restricted to regmask<=8 (cut at GSR)
selreg = (BSF.TLAT>=45) & (BSF.TLAT<=70) & (ncdat.REGION_MASK<=8)               # select region
BSFsel = BSF.where(selreg)                                                          # --> apply to BSF
BSFidx = BSFsel.max(axis=tuple([1,2]))                                          # weighted average over BSFsel
utils_misc.savexvar(BSFidx, path_var+'/BSFidx/'+code_BSFidx, 'BSFidx')          # save BSFidx

#               ------------------------------
code_BSFidx = 'BSFidx__{}__wmean_fullSPG_45to70_rmse8_288to314'.format(code_BSF) # CODE 1W >>>>> weighted mean over western SPG from 45N to 70N and restricted to regmask<=8 (cut at GSR)
selreg = (BSF.TLAT>=45) & (BSF.TLAT<=70) & (ncdat.REGION_MASK<=8)               # select region
selW = (BSF.TLONG>=280) & (BSF.TLONG<=314)                                # select western part of SPG
BSFsel = BSF.where(selreg & selW)                                                   # --> apply to BSF
TAREAsel = ncdat.TAREA.where(selreg & selW)                                         # --> apply to TAREA
BSFidx = (BSFsel*TAREAsel).mean(axis=tuple([1,2]))/ TAREAsel.mean()             # weighted average over BSFsel
utils_misc.savexvar(BSFidx, path_var+'/BSFidx/'+code_BSFidx, 'BSFidx')          # save BSFidx

code_BSFidx = 'BSFidx__{}__max_fullSPG_45to70_rmse8_288to314'.format(code_BSF)    # CODE 2W >>>>> maximum of western SPG from 45N to 70N and restricted to regmask<=8 (cut at GSR)
selreg = (BSF.TLAT>=45) & (BSF.TLAT<=70) & (ncdat.REGION_MASK<=8)               # select region
selW = (BSF.TLONG>=280) & (BSF.TLONG<=314)                                      # select western part of SPG
BSFsel = BSF.where(selreg & selW)                                                   # --> apply to BSF
BSFidx = BSFsel.max(axis=tuple([1,2]))                                          # weighted average over BSFsel
utils_misc.savexvar(BSFidx, path_var+'/BSFidx/'+code_BSFidx, 'BSFidx')          # save BSFidx

#               ------------------------------
code_BSFidx = 'BSFidx__{}__wmean_fullSPG_45to70_rmse8_215E'.format(code_BSF) # CODE 1E >>>>> weighted mean over eastern SPG from 45N to 70N and restricted to regmask<=8 (cut at GSR)
selreg = (BSF.TLAT>=45) & (BSF.TLAT<=70) & (ncdat.REGION_MASK<=8)               # select region
selE = BSF.TLONG>=215                                                           # select eastern part of SPG
BSFsel = BSF.where(selreg & selE)                                                   # --> apply to BSF
TAREAsel = ncdat.TAREA.where(selreg & selE)                                         # --> apply to TAREA
BSFidx = (BSFsel*TAREAsel).mean(axis=tuple([1,2]))/ TAREAsel.mean()             # weighted average over BSFsel
utils_misc.savexvar(BSFidx, path_var+'/BSFidx/'+code_BSFidx, 'BSFidx')          # save BSFidx

code_BSFidx = 'BSFidx__{}__max_fullSPG_45to70_rmse8_215E'.format(code_BSF)          # CODE 2E >>>>> maximum of eastern SPG from 45N to 70N and restricted to regmask<=8 (cut at GSR)
selreg = (BSF.TLAT>=45) & (BSF.TLAT<=70) & (ncdat.REGION_MASK<=8)               # select region
selE = BSF.TLONG>=215                                                           # select eastern part of SPG
BSFsel = BSF.where(selreg & selE)                                                   # --> apply to BSF
BSFidx = BSFsel.max(axis=tuple([1,2]))                                          # weighted average over BSFsel
utils_misc.savexvar(BSFidx, path_var+'/BSFidx/'+code_BSFidx, 'BSFidx')          # save BSFidx


# ==========================================================================================================================================================
# Correlations
# ==========================================================================================================================================================
# load streamfunctions
code_MOC = 'MOC__mod_850to1499'
code_BSF = 'BSF__mod_850to1499'
MOC = utils_misc.loadxvar(path_var+'/MOC/'+code_MOC, 'MOC')
BSF = utils_misc.loadxvar(path_var+'/BSF/'+code_BSF, 'BSF')
# load streamfunction-indices
code_MOCidx = 'MOCidx__MOC__mod_850to1499__max_fullATL_20downwards'
    #code_BSFidx = 'BSFidx__BSF__mod_850to1499__max_fullSPG_45to70_rmse8'
    #code_BSFidx = 'BSFidx__BSF__mod_850to1499__wmean_fullSPG_45to70_rmse8'
    #code_BSFidx = 'BSFidx__BSF__mod_850to1499__wmean_fullSPG_45to70_rmse8_288to314'
    #code_BSFidx = 'BSFidx__BSF__mod_850to1499__wmean_fullSPG_45to70_rmse8_215E'
    #code_BSFidx = 'BSFidx__BSF__mod_850to1499__max_fullSPG_45to70_rmse8_288to314'
    #code_BSFidx = 'BSFidx__BSF__mod_850to1499__max_fullSPG_45to70_rmse8_215E'

MOCidx = utils_misc.loadxvar(path_var+'/MOCidx/'+code_MOCidx, 'MOCidx')
BSFidx = utils_misc.loadxvar(path_var+'/BSFidx/'+code_BSFidx, 'BSFidx')

# lag in time sorted by increasing abs
lags = 0;   
for i in range(1,21): lags=np.append(lags,i); lags=np.append(lags,-i)
# number of degrees of freedom for p-value caluculation
num_df = 'full'               # set 'full' if full length of dataset shall be used

# -----------------------------------------------------------------------------
# corr of MOCidx with BSFidx
for lag in lags:
    corrMOCidxBSFidx, pcorrMOCidxBSFidx = utils_ana.xpearsonr(MOCidx, BSFidx, lag, num_df, return_pval=True)        # calculation of corr
    code_corrMOCidxBSFidx = 'corrMOCidxBSFidx__{}__{}__lag{}'.format(code_MOCidx, code_BSFidx, lag)                 # saving
    code_pcorrMOCidxBSFidx = 'corrMOCidxBSFidx__{}__{}__lag{}_df{}'.format(code_MOCidx, code_BSFidx, lag, num_df)
    utils_misc.savexvar(corrMOCidxBSFidx, path_var+'/corrMOCidxBSFidx/'+code_corrMOCidxBSFidx, 'corrMOCidxBSFidx',verbose=False)
    utils_misc.savexvar(pcorrMOCidxBSFidx, path_var+'/corrMOCidxBSFidx/'+code_pcorrMOCidxBSFidx, 'pcorrMOCidxBSFidx',verbose=False)

# -----------------------------------------------------------------------------
# pointwise corr. of MOC with BSFidx
for lag in lags:
    print '\nlag = {}'.format(lag)
    corrMOCBSFidx = xr.DataArray(np.nan*np.ones(shape = np.hstack([MOC.shape[-2:]])))                       # pre-allocation
    pcorrMOCBSFidx = xr.DataArray(np.zeros_like(corrMOCBSFidx))
    for j in np.arange(MOC.shape[-1]):                                                                      # point-wise calculation of corr
        utils_misc.Counter(j+1, MOC.shape[-1], 'Latitude')
        for k in np.arange(MOC.shape[-2]):
            corrMOCBSFidx[k,j], pcorrMOCBSFidx[k,j] = utils_ana.xpearsonr(MOC[:,k,j], BSFidx, lag, num_df, return_pval=True)
    code_corrMOCBSFidx = 'corrMOCBSFidx__{}__{}__lag{}'.format(code_MOC, code_BSFidx, lag)                   # saving
    code_pcorrMOCBSFidx = 'pcorrMOCBSFidx__{}__{}__lag{}_df{}'.format(code_MOC, code_BSFidx, lag, num_df)    
    utils_misc.savexvar(corrMOCBSFidx, path_var+'/corrMOCBSFidx/'+code_corrMOCBSFidx, 'corrMOCBSFidx',verbose=False)
    utils_misc.savexvar(pcorrMOCBSFidx, path_var+'/corrMOCBSFidx/'+code_pcorrMOCBSFidx, 'pcorrMOCBSFidx',verbose=False)

# -----------------------------------------------------------------------------
# pointwise corr. of BSF with MOCidx
for lag in lags:
    print '\nlag = {}'.format(lag)
    corrBSFMOCidx = xr.DataArray(np.zeros(shape = np.hstack([BSF.shape[-2:]])))                              # pre-allocation
    pcorrBSFMOCidx = xr.DataArray(np.zeros_like(corrBSFMOCidx))
    for j in np.arange(BSF.shape[-2]):                                                                      # point-wise calculation of corr
        utils_misc.Counter(j+1, BSF.shape[-2], 'Latitude')
        for i in np.arange(BSF.shape[-1]):
            corrBSFMOCidx[j,i], pcorrBSFMOCidx[j,i] = utils_ana.xpearsonr(BSF[:,j,i], MOCidx, lag, num_df, return_pval=True)
    code_corrBSFMOCidx = 'corrBSFMOCidx__{}__{}__lag{}'.format(code_BSF, code_MOCidx, lag)                   # saving
    code_pcorrBSFMOCidx = 'pcorrBSFMOCidx__{}__{}__lag{}_df{}'.format(code_BSF, code_MOCidx, lag, num_df)
    utils_misc.savexvar(corrBSFMOCidx, path_var+'/corrBSFMOCidx/'+code_corrBSFMOCidx, 'corrBSFMOCidx',verbose=False)
    utils_misc.savexvar(pcorrBSFMOCidx, path_var+'/corrBSFMOCidx/'+code_pcorrBSFMOCidx, 'pcorrBSFMOCidx',verbose=False)
