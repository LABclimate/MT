#################################
# The CESM python toolbox at KUP
# --- Meridional Overturning ---
#       Circulation Toolbox
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_utils_MOC as utils_MOC
#################################
# contained functions:
#################################
# - calc_MW()
# - calc_MOC_model_grid()
# - get_default_auxgrd()
# - calc_Mxint_auxgrd()
# - calc_MOC_auxgrd()
#################################
# please log your changes below:
#################################
# 17-Mai-2016 - buerki@climate.unibe.ch : created this toolboxfew
# 24?-Mai-2016 -buerki@climate.unibe.ch : seperated calc_Mxint_auxgrd() from calc_MOC_auxgrd()
#################################

import numpy as np
import xarray as xr
import pickle
import CESM_utils_mask as utils_mask
import UTILS_misc as utils_misc

# =======================================================================================
# - Compute vertical volume transport MW (in Sv)
# =======================================================================================
def calc_MW(ncdat):
    ''' 
    Comments:
     > Conversion from cgs units to Sv by multiplication with 1e-12
    '''
    wvel = utils_mask.mask_ATLANTIC(ncdat.WVEL.mean(dim='time'), ncdat.REGION_MASK)
    TAREA = ncdat.TAREA # z-area of T cells
    MW = xr.DataArray(wvel*TAREA*1e-12,
		    name='vertical volume transport',
 		    attrs={'units':u'Sv'})
    return(MW)

# =======================================================================================
# - MOC on model grid
# =======================================================================================
def calc_MOC_mgrd(velocity_component, M, do_normalize=True, dump_Mxint=False):
    '''
    Input:
     > velocity_component 	: either 'W' or  'V' (string)
     > M 			: either MW or MV (either as xarray or nparray)
     > do_normalize 		: boolean
     > dump_Mxint 		: boolean
    Comments:
     > As latitude varies along longitude of model grid I simply set the TLAT as 
       the *mean* of M-latitudes along longitudes.
       Note, that this is very inappropriate at high latitudes!
     > Think about taking np.nansum() #!
    '''
    # zonal integration along model grid    
    Mxint = xr.DataArray(M.sum(dim='nlon'),  	# zonal integration
		    attrs={'units':u'Sv'}, 
		    coords={'TLAT':M.TLAT.mean(dim='nlon')}) 	#! mean is inappropriate at high latitudes!
    # meridional integration along model grid
    MOC = xr.DataArray(Mxint)
    for j in np.arange(1,len(Mxint.nlat)): 	# meridional integration
      MOC[dict(nlat=j)] = MOC.isel(nlat=j) + MOC.isel(nlat=j-1)
    # normalization relative to North (shift values such that zero at northern boundary)
    if do_normalize == True:
      MOC = MOC - MOC[:,-1]
    # naming xarrays
    if velocity_component == 'W':
      Mxint.name = 'MW zonally integrated'
      MOC.name = 'MOC on model grid calculated from WVEL'
    elif velocity_component == 'V':
      Mxint.name = 'MV zonally integrated'
      MOC.name = 'MOC on model grid calculated from VVEL'

    if dump_Mxint == True:
      return(MOC, Mxint)
    else:
      return(MOC)

# =======================================================================================
# - MOC on auxillary grid
# =======================================================================================
'''
    This collection contains 3 functions, which :
     A: returns a default for an auxillary grid. 
     B: integration of vertical(MW)/meridional(MV) transport in zonal direction.
     C: integration of the result of B in meridional(MW)/vertical(MV) direction.

    Comments:
    =========
    > As the np.nansum() is used for zonal integration, ridges are ignored and integrated through. 
    This might be inappropriate. --> it's just the definition of the streamfunction...
    > For speed reasons some data from xarrays is copied to np-arrays
    > For speed reasons absolute indexing is used with 
    n: latitude on aux-grid, 
    j: latitude on model-grid, 
    i: longitude on model-grid and k: depth on both grids
'''
# ---------------------------------------------------------------------------------------
# - define default for auxillary grid
# ---------------------------------------------------------------------------------------
def get_default_auxgrd(ncdat):
    lat = np.linspace(-80, 90, 170)  	# latitudes
    z_t = ncdat.z_t.values 		# depth levels
    z_w_top = ncdat.z_w_top.values 	# depth levels
    return(lat, z_t, z_w_top)

# ---------------------------------------------------------------------------------------
# - zonal integration of Volume Transport along auxillary grid
# ---------------------------------------------------------------------------------------
def calc_Mxint_auxgrd(lat_auxgrd, z_auxgrd, velocity_component, M, ncdat, do_normalize=True, savevar=True):
    '''
    Input:
     > lat_auxgrd               : vector with meridional auxillary grid
     > z_auxgrd                 : vector with vertical auxillary grid
     > velocity_component 	: either 'W' or  'V' (string)    
     > M                        : volume transport (NOT rolled!!)
     > ncdat                    : netCDFdata to load REGION_MASK

    Steps:
     > generate mask (mask_auxgrd) that is 'True' where latitudes of both grids lie in the same box
     > generate array (maxiter_depth) for seafloor detection. It contains the indices of maximal depth
       for each point on model grid (T points) in order to stop k-iteration at the seafloor.
     > calculate Mxint by zonal integration along aux grid
     > n-loop: over latitude on aux grid
     >   j-loop: over latitudes on model grid
     >     check whether mask_auxgrd of current n, j is True anywhere (just for speeding up)
     >     i-loop: over those longitudes where both, mask_modgrd (regional mask) and mask_auxgrd are True.
                   the order of the iteration is rolled to begin at western boundary of Atlantic
                   but actually, this 
     >       k-loop: over depths from surface down to depth of seafloor relative to model grid (j,i position)
     >         zonal integration by summing up vertical volume transports (M) of model grid.
    '''
    # a few variables to speed up subsequent loops 
    iter_lat_auxgrd = np.arange(len(lat_auxgrd))
    iter_lat_M = np.arange(len(M.nlat))

    # get masks and iteration-indices to speed up subsequent loops (recalculate if loading from file fails)
    try: 	mask_auxgrd = utils_misc.loadvar('variables/mask_auxgrd')
    except: 	mask_auxgrd = utils_mask.gen_mask_grd_overlay_lat(lat_auxgrd, ncdat)
    try:	iter_maskcombo = utils_misc.loadvar('variables/iter_maskcombo')
    except:     iter_maskcombo = utils_mask.gen_iter_maskcombo(lat_auxgrd, ncdat, mask_auxgrd)
    try:	maxiter_depth = utils_misc.loadvar('variables/maxiter_depth') 
    except:     maxiter_depth = utils_mask.gen_maxiter_depth(lat_auxgrd, z_auxgrd, ncdat)
    
    # zonal integration along aux grid
    print('> zonal integration')
    Mxint = np.zeros([len(z_auxgrd), len(lat_auxgrd)])  	# pre-allocation with zeros (np-array like for speed)
    for n in iter_lat_auxgrd:
      utils_misc.ProgBar('step', barlen=60, step=n, nsteps=len(iter_lat_auxgrd))# initialize and update progress bar
      for j in iter_lat_M:
#        if any(mask_auxgrd[n,j,:]):						# to speed up the code
          for i in iter_maskcombo[n,j]: 					# limit zonal integration to Atlantic and grid-overlay
            for k in np.arange(int(maxiter_depth[j,i])): 			# stop at depth of seafloor
              Mxint[k,n] = np.nansum([Mxint[k,n],M[k,j,i]]) 			# zonal integration
    utils_misc.ProgBar('done') 							# terminate progress bar
    
    # write Mxint to xarray
    Mxint= xr.DataArray(Mxint,
		attrs={'units':u'Sv'},
                coords={'z_w_top':np.arange(len(z_auxgrd)), 'nlat':np.arange(len(lat_auxgrd))},
		dims=['z_w_top', 'nlat'])

    if velocity_component == 'W':
      Mxint.name = 'MW zonally integrated'                                      # naming xarray
      if savevar == True: utils_misc.savevar(Mxint, 'variables/MWxint_auxgrd')  # save to file
    elif velocity_component == 'V':
      Mxint.name = 'MV zonally integrated'                                      # naming xarray
      if savevar == True: utils_misc.savevar(Mxint, 'variables/MVxint_auxgrd')  # save to file

    return(Mxint)

# ---------------------------------------------------------------------------------------
# - MOC on auxillary grid
# ---------------------------------------------------------------------------------------
def calc_MOC_auxgrd(lat_auxgrd, z_auxgrd, velocity_component, Mxint, ncdat, do_normalize=True, savevar=True):
    '''
    Input:
     > lat_auxgrd               : vector with meridional auxillary grid
     > z_auxgrd                 : vector with vertical auxillary grid
     > velocity_component 	: either 'W' or  'V' (string)    
     > Mxint                    : zonally integrated volume transport (NOT rolled!!)
     > ncdat                    : netCDFdata to load REGION_MASK

    Steps:
     > calculate MOC by meridional of Mxint integration along aux grid
     > normalization relative to northern boundary: at every point substract northernmost value at same depth, 
       such that streamfunction closes at NP.
    '''
    # a few variables to speed up subsequent loops 
    iter_lat_auxgrd = np.arange(len(lat_auxgrd))
    iter_z_auxgrd = np.arange(len(z_auxgrd))

    # preallocation of MOC as np-array
    MOC = np.copy(Mxint) # start with Mxint, which subsequently will be summed up

    if velocity_component == 'W':
      # meridional integration along aux grid
      print('> meridional integration')
      for n in iter_lat_auxgrd[1:]:
        utils_misc.ProgBar('step', barlen=60, step=n, nsteps=len(iter_lat_auxgrd), forceinit=True)      # initialize and update progress bar
        MOC[:,n] = np.nansum([MOC[:,n], MOC[:,n-1]], axis=0) 	        # meridional integration
      utils_misc.ProgBar('done') 					# terminate progress bar

      xrname = 'MOC on auxillary grid calculated from WVEL'             # name of xarray
      fname = 'MOC_auxgrd_W'                                            # name for saving

    elif velocity_component == 'V':
      # vertical integration along aux grid
      print('> vertical integration')
      for k in iter_z_auxgrd[1:]:
        utils_misc.ProgBar('step', barlen=60, step=k, nsteps=len(iter_z_auxgrd), forceinit=True)        # initialize and update progress bar
        MOC[k,:] = np.nansum([MOC[k,:], MOC[k-1,:]], axis=0) 		# meridional integration
      utils_misc.ProgBar('done') 					# terminate progress bar

      xrname = 'MOC on auxillary grid calculated from VVEL'             # name of xarray
      fname = 'MOC_auxgrd_V'                                            # name for saving

    # write to xarray
    MOC = xr.DataArray(MOC,
		attrs={'units':u'Sv'},
                coords=[np.arange(len(z_auxgrd)), np.arange(len(lat_auxgrd))],
		dims=['z_w_top', 'nlat'],
                name=xrname)

    # normalization relative to North (shift values such that zero at northern boundary)
    if do_normalize == True:
      MOC = MOC - MOC[:,-1]

    # save to file
    if savevar == True:
      utils_misc.savevar(MOC, 'variables/'+fname)

    return(MOC)
