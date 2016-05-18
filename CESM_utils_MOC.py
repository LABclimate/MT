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
import UTILS_specials as utils_spec 	# for Progress Bars

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
def calc_MOC_on_modgrd(MW, dump_MWxint=False):
    '''
    Comments:
     > As latitude varies along longitude of model grid I simply set the TLAT as 
       the *mean* of MW-latitudes along longitudes.
       Note, that this is very inappropriate at high latitudes!
     > Think about taking np.nansum() #!
    '''
    MWxint = xr.DataArray(MW.sum(dim='nlon'),  	# zonal integration
		    name='MW zonally integrated',
		    attrs={'units':u'Sv'}, 
		    coords={'TLAT':MW.TLAT.mean(dim='nlon')}) 	#! mean is inappropriate at high latitudes!

    MOC = xr.DataArray(MWxint, name='MOC on model grid', attrs={'units':u'Sv'})
    for j in np.arange(1,len(MWxint.nlat)): 	# meridional integration
      MOC[dict(nlat=j)] = MOC.isel(nlat=j) + MOC.isel(nlat=j-1)
    
    if dump_MWxint == True:
      return(MOC, MWxint)
    else:
      return(MOC)

# =======================================================================================
# - MOC on auxillary grid
# =======================================================================================

# ---------------------------------------------------------------------------------------
# - define default for auxillary grid
# ---------------------------------------------------------------------------------------
def get_default_auxgrd(MW):
    lat_auxgrd = np.linspace(-90, 90, 180)  	# latitudes
    z_w_auxgrd = MW.z_w_top.values 		# depth levels
    return(lat_auxgrd, z_w_auxgrd)

# ---------------------------------------------------------------------------------------
# - MOC on auxillary grid
# ---------------------------------------------------------------------------------------
def calc_MOC_on_auxgrd(lat_auxgrd, z_w_auxgrd, MW, ncdat, dump_MWxint=False, savevar=True):
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
    iter_lon_MW = np.array(MW.nlon) 			 			# in normal order!!	

    # get masks and iteration-indices to speed up subsequent loops
    try:
      with open('variables/mask_auxgrd') as f: mask_auxgrd = pickle.load(f)
      print('success in loading ''mask_auxgrd'' from file')
    except:
      mask_auxgrd = utils_mask.gen_mask_grd_overlay_lat(lat_auxgrd, MW)
    try:
      with open('variables/iter_maskcombo') as f: iter_maskcombo = pickle.load(f)
      print('success in loading ''iter_maskcombo'' from file')      
    except:
      iter_maskcombo = utils_mask.gen_iter_maskcombo(lat_auxgrd, MW, mask_auxgrd, ncdat.REGION_MASK)
    try:
      with open('variables/maxiter_depth') as f: maxiter_depth = pickle.load(f)
      print('success in loading ''maxiter_depth'' from file')      
    except:
      maxiter_depth = utils_mask.gen_maxiter_depth(lat_auxgrd, z_w_auxgrd, MW, ncdat.HT.values)
    
    # zonal integration along aux grid
    print('> zonal integration')
    MWxint = np.zeros([len(lat_auxgrd), len(z_w_auxgrd)])  		# pre-allocation with zeros
    for n in iter_lat_auxgrd:
      utils_spec.ProgBar('step', barlen=60, step=n, nsteps=len(iter_lat_auxgrd))# initialize and update progress bar
      for j in iter_lat_MW:
        if any(mask_auxgrd[n,j,:]):						# to speed up the code
          for i in iter_maskcombo[n,j]: 					# limit zonal integration to Atlantic and grid-overlay
            for k in np.arange(int(maxiter_depth[j,i])): 			# stop at depth of seafloor
              MWxint[n,k] = np.nansum([MWxint[n,k],MW[k,j,i]]) 			# zonal integration
              print(MWxint[n,k])
    utils_spec.ProgBar('done') 							# terminate progress bar

    # meridional integration along aux grid
    print('> meridional integration')
    MOC = xr.DataArray(MWxint, 
		name='MOC on auxillary grid', 
		attrs={'units':u'Sv'},
		dims={'nlat':np.arange(len(lat_auxgrd)), 'z_w_top':np.arange(len(z_w_auxgrd))})
    for n in iter_lat_auxgrd[1:]:
      utils_spec.ProgBar('step', barlen=60, step=n, nsteps=len(iter_lat_auxgrd))# initialize and update progress bar
      for ii in 
      MOC[dict(nlat=n)] = np.nansum([MOC[n,:], MOC[n-1,:]], axis=0) 	# meridional integration
    utils_spec.ProgBar('done') 							# terminate progress bar
    
    if savevar == True:
      with open('variables/MOC', 'w') as f: pickle.dump(MOC, f) 	# save variable

     
    if dump_MWxint == True:
      return(MOC, MWxint)
    else:
      return(MOC)



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

