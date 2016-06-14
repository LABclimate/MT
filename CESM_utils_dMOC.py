#################################
# The CESM python toolbox at KUP
# ----- Density dMOC Toolbox -----
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_utils_dMOC as utils_dMOC
#################################
# contained functions:
#################################
# - calc_dMOC_mgrd()
# - calc_dMxint_auxgrd()
# - calc_dMOC_auxgrd()
#################################
# please log your changes below:
#################################
# 02-Jun-2016 - buerki@climate.unibe.ch : created this toolbox
#################################

import numpy as np
import xarray as xr
import pickle
import CESM_utils_mask as utils_mask
import UTILS_misc as utils_misc


# =======================================================================================
# - dMOC on model grid
# =======================================================================================
def calc_dMOC_mgrd(vel_comp, M, PD, PD_bins, do_norm=True, dump_dMxint=False):
    '''
    Input:
     > vel_comp 	        : either 'W' or  'V' | string
     > M 			: volume transport (MW or MV) | nparray of shape [nz, nlat, nlon]
     > PD                       : potential density | nparray of shape [nz, nlat, nlon]
     > PD_bins                  : borders of PD-bins | nparray of shape [nPDbins+1]
     > do_norm 		        : boolean
     > dump_dMxint 		: boolean
    Output:
     > dMxint                   : zonally integrated volume transport of shape [nPDbins, nlat] | nparray
     > dMOC                     : dMOC of shape [nPDbins, nlat] | nparray
    '''
    iter_dens = np.arange(len(PD_bins)-1)
    # pre-allocation
    mask_PD_bins = np.zeros(shape=(len(PD_bins)-1, PD.shape[1]), dtype=object)
    dMxint = np.zeros(shape=mask_PD_bins.shape)

    # zonal integration and conversion on density axis
    for l in iter_dens:
      utils_misc.ProgBar('step', step=l, nsteps=len(PD_bins)-1)
      for j in np.arange(PD.shape[1]):
        mask_PD_bins[l,j] = np.where( (PD[:,j,:]>PD_bins[l]) & (PD[:,j,:]<PD_bins[l+1]) )
        dMxint[l,j] = np.nansum(M[mask_PD_bins[l,j][0], j, mask_PD_bins[l,j][1]])
    utils_misc.ProgBar('done')

    # meridional integration along model grid
    dMOC = np.copy(dMxint)                          # pre-allocation with dMxint
    for j in np.arange(1,dMxint.shape[1]): 	    # meridional integration S --> N
      dMOC[:,j] = dMOC[:,j] + dMOC[:,j-1]
#    for j in np.arange(0,dMxint.shape[1]-1)[::-1]:  # meridional integration N --> S
#      dMOC[:,j] = dMOC[:,j] + dMOC[:,j+1]

    # normalisation relative to North (shift values such that zero at northern boundary)
    if do_norm == True:
      dMOC = dMOC - dMOC[:,-1]

    '''
    # write to xarray
    dMOC = xr.DataArray(dMOC, attrs={'units':u'Sv'}, 
            coords={'nsigma2':np.arange(len(PD_bins)-1), 'sigma2':PD_bins[:-1], 'nlat':np.arange(PD.shape[1]), 'TLAT':ncdat.TLAT.mean(dim='nlon')},  	#! mean is inappropriate at high latitudes!
		            dims=['nsigma2', 'nlat'])

    # naming xarrays
    if vel_comp == 'W':
      dMxint.name = 'MW zonally integrated'
      dMOC.name = 'dMOC on model grid calculated from WVEL'
    elif vel_comp == 'V':
      dMxint.name = 'MV zonally integrated'
      dMOC.name = 'dMOC on model grid calculated from VVEL'
    '''

    if dump_dMxint == True:
      return(dMOC, dMxint)
    else:
      return(dMOC)

# =======================================================================================
# - dMOC on auxillary grid
# =======================================================================================
'''
REWRITE !!!!
    This collection contains 2 functions:
     A: integration of vertical(MW)/meridional(MV) transport in zonal direction.
     B: integration of the result of B in meridional(MW)/vertical(MV) direction.

    Comments:
    =========
    > As the np.nansum() is used for zonal integration, ridges are ignored and integrated through. 
    This might be inappropriate. --> it's just the definition of the streamfunction...
    > For speed reasons some data from xarrays is copied to np-arrays
    > For speed reasons absolute indexing is used with 
    n: latitude on aux-grid, 
    j: latitude on model-grid, 
    i: longitude on model-grid and k: depth on both grids
    k: depth on both grids and
    l: density
REWRITE !!!!
'''
# ---------------------------------------------------------------------------------------
# - zonal integration of Volume Transport along auxillary grid
# ---------------------------------------------------------------------------------------
def calc_dMxint_auxgrd(lat_auxgrd, z_auxgrd, vel_comp, M, PD, PD_bins, ncdat, path_vars, savevar=True):
    '''
    Input:
     > lat_auxgrd               : meridional auxillary grid | nparray
     > z_auxgrd                 : vertical auxillary grid | nparray
     > vel_comp         	: either 'W' or 'V' | string
     > M 			: volume transport (MW or MV) | nparray of shape [nz, nlat, nlon]
     > PD                       : potential density | nparray of shape [nz, nlat, nlon]
     > PD_bins                  : borders of PD-bins | nparray of shape [nPDbins+1]
     > ncdat                    : netCDFdata to load REGION_MASK
     > path_vars                : path for saving variables | string
     > savevar 		        : boolean
    Output:
     > dMxint                   : zonally integrated volume transport of shape [nPDbins, nlat] | nparray
    Steps:
     > Todo
    '''
    # a few variables to speed up subsequent loops 
    iter_lat_auxgrd = np.arange(len(lat_auxgrd))
    iter_lat_M = np.arange(M.shape[1])
    iter_dens = np.arange(len(PD_bins)-1)
    # get masks and iteration-indices to speed up subsequent loops (calculate if loading from file fails)
    try: 	mask_auxgrd = utils_misc.loadvar(path_vars+'mask_auxgrd')
    except: 	mask_auxgrd = utils_mask.gen_mask_grd_overlay_lat(lat_auxgrd, ncdat, path_vars)
    try:	iter_maskcombo = utils_misc.loadvar(path_vars+'iter_maskcombo')
    except:     iter_maskcombo = utils_mask.gen_iter_maskcombo(lat_auxgrd, ncdat, mask_auxgrd, path_vars)
    try:	maxiter_depth = utils_misc.loadvar(path_vars+'maxiter_depth') 
    except:     maxiter_depth = utils_mask.gen_maxiter_depth(lat_auxgrd, z_auxgrd, ncdat, path_vars)
    # pre-allocation with zeros
    mask_PD_bins = np.zeros(shape=(len(PD_bins)-1, PD.shape[1]), dtype=object)
    dMxint = np.zeros(shape=mask_PD_bins.shape)
    # zonal integration along auxgrid and conversion on density axis
    print('> zonal integration')
    for n in iter_lat_auxgrd:
      utils_misc.ProgBar('step', step=n, nsteps=len(iter_lat_auxgrd))
      for l in iter_dens:
        for j in iter_lat_M:
          mask_PD_bins[l,j] = np.where( (PD[:,j,:]>PD_bins[l]) & (PD[:,j,:]<PD_bins[l+1]) ) # depth and longitude
          for i in iter_maskcombo[n,j]: # limit zonal integration to Atlantic and grid-overlay
            if i in mask_PD_bins[l,j][1]: # combine both masks: for maskcombo and for densitybinning
              try: dMxint[l,n] = np.nansum([dMxint[l,n], M[mask_PD_bins[l,j][0],j,i]])
              except: pass # just for the cases where M[...] is an empty array, which is not accepted by nansum
    utils_misc.ProgBar('done')

    if vel_comp == 'W':
      if savevar == True: utils_misc.savevar(dMxint, path_vars+'dMWxint_auxgrd')  # save to file
    elif vel_comp == 'V':
      if savevar == True: utils_misc.savevar(dMxint, path_vars+'dMVxint_auxgrd')  # save to file

    return(dMxint)



# ---------------------------------------------------------------------------------------
# - dMOC on auxillary grid
# ---------------------------------------------------------------------------------------
def calc_dMOC_auxgrd(lat_auxgrd, dens_auxgrd, vel_comp, dMxint, ncdat, path_vars, do_norm=True, savevar=True):
    '''
    Input:
     > lat_auxgrd               : meridional auxillary grid | nparray
     > dens_auxgrd              : density auxillary grid | nparray
     > vel_comp         	: either 'W' or 'V' | string
     > dMxint                   : zonally integrated volume transport
     > ncdat                    : netCDFdata to load REGION_MASK
     > path_vars                : path for saving variables | string
     > do_norm                  : do normalisation relative to northern boundary | boolean     
     > savevar                  : boolean
    Output:
     > dMOC                     : dMOC of shape [nPDbins, nlat] | nparray
    Steps:
     > calculate dMOC by meridional of dMxint integration along aux grid
     > normalisation relative to northern boundary: at every point substract northernmost value at same depth, 
       such that streamfunction closes at NP.
    '''
    # a few variables to speed up subsequent loops 
    iter_lat_auxgrd = np.arange(len(lat_auxgrd))
    iter_dens_auxgrd = np.arange(len(dens_auxgrd))

    # preallocation of dMOC as np-array
    dMOC = np.copy(dMxint) # start with dMxint, which subsequently will be summed up

    if vel_comp == 'W':
      # meridional integration along aux grid
      print('> meridional integration')
      for n in iter_lat_auxgrd[1:]:
        utils_misc.ProgBar('step', step=n, nsteps=len(iter_lat_auxgrd), forceinit=True, minbarlen=60)
        dMOC[:,n] = np.nansum([dMOC[:,n], dMOC[:,n-1]], axis=0) 	        # meridional integration
      utils_misc.ProgBar('done')

      xrname = 'dMOC on auxillary grid calculated from WVEL'             # name of xarray
      fname = 'dMOC_auxgrd_W'                                            # name for saving

    elif vel_comp == 'V':
      # vertical integration along aux grid
      print('> vertical integration')
      for k in iter_dens_auxgrd[1:]:
        utils_misc.ProgBar('step', step=k, nsteps=len(iter_dens_auxgrd), forceinit=True, minbarlen=60)
        dMOC[k,:] = np.nansum([dMOC[k,:], dMOC[k-1,:]], axis=0) 		# meridional integration
      utils_misc.ProgBar('done')

      xrname = 'dMOC on auxillary grid calculated from VVEL'             # name of xarray
      fname = 'dMOC_auxgrd_V'                                            # name for saving

    '''
    # write to xarray
    dMOC = xr.DataArray(dMOC,
		attrs={'units':u'Sv'},
                coords=[np.arange(len(dens_auxgrd)), np.arange(len(lat_auxgrd))],
		dims=['dens', 'nlat'],
                name=xrname)
    '''
    # normalisation relative to North (shift values such that zero at northern boundary)
    if do_norm == True:
      dMOC = dMOC - dMOC[:,-1]

    # save to file
    if savevar == True:
      utils_misc.savevar(dMOC, path_vars + fname)

    return(dMOC)
