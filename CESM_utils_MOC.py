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
# - calc_MOC_on_model_grid()
# - get_default_auxgrd()
# - calc_MOC_on_auxgrd()
#################################
# please log your changes below:
#################################
# 17-Mai-2016 - buerki@climate.unibe.ch : created this toolboxfew
#################################

import numpy as np
import xarray as xr
import pickle
import CESM_utils_mask as utils_mask
import UTILS_specials as utils_spec

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
def calc_MOC_mgrd(M, velocity_component, do_normalize=True, dump_Mxint=False):
    '''
    Input:
     > M 			: either MW or MV (either as xarray or nparray)
     > velocity_component 	: either 'W' or  'V' (string)
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

# ---------------------------------------------------------------------------------------
# - define default for auxillary grid
# ---------------------------------------------------------------------------------------
def get_default_auxgrd(ncdat):
    lat = np.linspace(-90, 90, 180)  	# latitudes
    z_t = ncdat.z_t.values 		# depth levels
    z_w_top = ncdat.z_w_top.values 	# depth levels
    return(lat, z_t, z_w_top)

# ---------------------------------------------------------------------------------------
# - MOC on auxillary grid based on WVEL
# ---------------------------------------------------------------------------------------
def calc_MOC_auxgrd(lat_auxgrd, z_auxgrd, MW, ncdat, do_normalize=True, dump_MWxint=False, savevar=True):
    '''
    Input:
    ======
    > MW : vertical volume transport (NOT rolled!!)

    Steps:
    ======
    > define auxillary grid
    > generate mask (mask_auxgrd) that is 'True' where latitudes of both grids lie in the same box
    > generate array (maxiter_depth) for seafloor detection. It contains the indices of maximal depth
    for each point on model grid (T points) in order to stop k-iteration at the seafloor.
    > calculate MWxint by zonal integration along aux grid
    > n-loop: over latitude on aux grid
    >   j-loop: over latitudes on model grid
    >     check whether mask_auxgrd of current n, j is True anywhere (just for speeding up)
    >     i-loop: over those longitudes where both, mask_modgrd (regional mask) and mask_auxgrd are True.
    the order of the iteration is rolled to begin at western boundary of Atlantic
    but actually, this 
    >       k-loop: over depths from surface down to depth of seafloor relative to model grid (j,i position)
    >         zonal integration by summing up vertical volume transports (MW) of model grid.
    > calculate MOC by meridional integration along aux grid
    > normalization relative to northern boundary: at every point substract northernmost value at same depth, 
    such that streamfunction closes at NP.

    Comments:
    =========
    > As the np.nansum() is used for zonal integration, ridges are ignored and integrated through. 
    This might be inappropriate.
    > For speed reasons some variables are written to new variables
    > For speed reasons absolute indexing is used with 
    n: latitude on aux-grid, 
    j: latitude on model-grid, 
    i: longitude on model-grid and k: depth on both grids
    '''
    # a few variables to speed up subsequent loops 
    lat_MW = np.array(MW.TLAT)
    iter_lat_auxgrd = np.arange(len(lat_auxgrd))
    iter_lat_MW = np.arange(len(MW.nlat))

    # get masks and iteration-indices to speed up subsequent loops (recalculate if loading from file fails)
    try: 	mask_auxgrd = utils_spec.loadvar('variables/mask_auxgrd')
    except: 	mask_auxgrd = utils_mask.gen_mask_grd_overlay_lat(lat_auxgrd, MW)
    try:	iter_maskcombo = utils_spec.loadvar('variables/iter_maskcombo')
    except:     iter_maskcombo = utils_mask.gen_iter_maskcombo(lat_auxgrd, MW, mask_auxgrd, ncdat.REGION_MASK)
    try:	maxiter_depth = utils_spec.loadvar('variables/maxiter_depth') 
    except:     maxiter_depth = utils_mask.gen_maxiter_depth(lat_auxgrd, z_auxgrd, MW, ncdat.HT.values)
    
    # zonal integration along aux grid
    print('> zonal integration')
    MWxint = np.zeros([len(lat_auxgrd), len(z_auxgrd)])  	# pre-allocation with zeros (np-array like for speed)
    for n in iter_lat_auxgrd:
      utils_spec.ProgBar('step', barlen=60, step=n, nsteps=len(iter_lat_auxgrd))# initialize and update progress bar
      for j in iter_lat_MW:
#        if any(mask_auxgrd[n,j,:]):						# to speed up the code
          for i in iter_maskcombo[n,j]: 					# limit zonal integration to Atlantic and grid-overlay
            for k in np.arange(int(maxiter_depth[j,i])): 			# stop at depth of seafloor
              MWxint[n,k] = np.nansum([MWxint[n,k],MW[k,j,i]]) 			# zonal integration
    MWxint= xr.DataArray(MWxint, 						# write MWxint to xarray
		name='MW zonally integrated along auxillary grid', 
		attrs={'units':u'Sv'},
		dims={'nlat':np.arange(len(lat_auxgrd)), 'z_w_top':np.arange(len(z_auxgrd))})
    utils_spec.ProgBar('done') 							# terminate progress bar

    # meridional integration along aux grid
    print('> meridional integration')
    MOC = xr.DataArray(MWxint, 
		name='MOC on auxillary grid', 
		attrs={'units':u'Sv'},
		dims={'nlat':np.arange(len(lat_auxgrd)), 'z_w_top':np.arange(len(z_auxgrd))})
    for n in iter_lat_auxgrd[1:]:
      utils_spec.ProgBar('step', barlen=60, step=n, nsteps=len(iter_lat_auxgrd))# initialize and update progress bar
      MOC[dict(nlat=n)] = np.nansum([MOC[n,:], MOC[n-1,:]], axis=0) 		# meridional integration
    utils_spec.ProgBar('done') 							# terminate progress bar
    
    if savevar == True:
      utils_spec.savevar(MOC, 'variables/MOC_auxgrd')				# save variable #! change name

    # normalization relative to North (shift values such that zero at northern boundary)
    if do_normalize == True:
      MOC = MOC - MOC.isel(nlat=-1)
    
    if dump_MWxint == True:
      return(MOC, MWxint)
    else:
      return(MOC)


# ---------------------------------------------------------------------------------------
# - MOC on auxillary grid based on VVEL
# ---------------------------------------------------------------------------------------
#def calc_MOC_auxgrd_VVEL(lat_auxgrd, z_auxgrd, MW, ncdat, do_normalize=True, dump_MWxint=False, savevar=True):


#################################################################################################
#################################################################################################
## OOOOLD STUFF
#################################################################################################
#################################################################################################

#MOC = xr.DataArray(np.zeros([len(lat_auxgrd),len(MW.z_w_top)]),name='MOC on aux grid', dims=['nlat_aux', 'z_w_top'], coords = {'nlat_aux':np.arange(len(lat_auxgrd)),'z_w_top':MW.z_w_top})
#for n in np.arange(1,len(lat_aux_grd)+1)
#  MOC[dict(nlat_aux=n)] = MOC.isel(nlat_aux=n-1)
#  for j in np.arange(len(MW.nlat))
#    for i in np.arange(len(MW.nlon))
#      if lat_auxgrd[n] <= MW.TLAT.isel(nlon=i, nlat=j) < lat_auxgrd[n+1]:
#        for k in np.arange(len(MW.z_w_top))
#          MOC[dict(nlat_aux=n, z_w_top=k)] = MOC[n,k] + MW.isel(nlon=i, nlat=j, z_w_top=k)
