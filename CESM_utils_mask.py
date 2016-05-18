#################################
# The CESM python toolbox at KUP
# ------ Masking-Toolbox -------
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_utils_mask as utils_mask
#################################
# contained functions:
#################################
# - mask_ATLANTIC()
# - vars2speedup()
# - gen_mask_grd_overlay_lat()
# - gen_iter_maskcombo()
# - get_maxiter_depth()
#################################
# please log your changes below:
#################################
# 30-Apr-2016 - buerki@climate.unibe.ch : created this toolbox
#                                         added mask_all_but_reg_6_8_9_10()
# 02-Mai-2016 - buerki@climate.unibe.ch : replaced mask_all_but_reg_6_8_9_10 by mask_ATLANTIC()
# 14-Mai-2016 - buerki@climate.unibe.ch : in mask_ATLANTIC added optional argument 'outputformat'
# 17-Mai-2016 - buerki@climate.unibe.ch : created vars2speedup()
#                                         created gen_mask_grd_overlay_lat()
#                                         created gen_iter_maskcombo()
#                                         created get_maxiter_depth()
#################################

import numpy as np
import numpy.ma as ma
import pickle
import CESM_utils_mask as utils_mask
import CESM_utils_plt as utils_plt
import CESM_utils_conv as utils_conv
import UTILS_specials as utils_spec 	# for Progress Bars

# =======================================================================================
# mask any xarray with respect to values in the regional mask of the nc-data
# =======================================================================================
''' Mask (hide) everything except the regions... 
         6 (North Atlantic),
	 7 (Mediteranian),
         8 (Labrador/Davis/Baffin Seas), 
         9 (North Sea and part of Norvegian Sea),
        10 (Norvegian Sea and Northern Seas) and
	11 (Hudson Bay)
'''
def mask_ATLANTIC(varin, mask, outputformat='xr'):
    if outputformat=='ma':
        return(ma.masked_array(varin, mask=np.array(mask>=6)))
    elif outputformat=='xr':
        return(varin.where(mask>=6))


# =======================================================================================
# generate masks that we use in MOC-functions
# =======================================================================================

# write some variables to numpy-arrays in order to speed up subsequent loops
def vars2speedup(lat_auxgrd, MW):
    lat_MW = np.array(MW.TLAT)
    iter_lat_auxgrd = np.arange(len(lat_auxgrd))
    iter_lat_MW = np.arange(len(MW.nlat))
    iter_lon_MW = np.array(MW.nlon) 			# in normal order!!
    return(lat_MW, iter_lat_auxgrd, iter_lat_MW, iter_lon_MW)

# --------------------------------------------------
# generate mask_auxgrd, a mask for grid-overlay
# --------------------------------------------------
def gen_mask_grd_overlay_lat(lat_auxgrd, MW, savevar=True):
    ''' Boolean of size (nlatAUXgrid, nlatMODELgrid, nlonMODELgrid)
        It is True where latitudes of auxillary and model grid lie in the same box. 
    '''
    print('> generating mask_auxgrd')
    lat_MW, iter_lat_auxgrd, iter_lat_MW, iter_lon_MW = vars2speedup(lat_auxgrd, MW) 	# np-arrays for speed
    mask_auxgrd = np.zeros([len(lat_auxgrd), len(MW.nlat), len(MW.nlon)],dtype=bool) 	# pre-allocation as False

    for n in iter_lat_auxgrd:
      utils_spec.ProgBar('step', barlen=30, step=n, nsteps=len(iter_lat_auxgrd)) 	# initialize and update progress bar
      for j in iter_lat_MW:
        for i in iter_lon_MW:
          if lat_auxgrd[n] <= lat_MW[j,i] < lat_auxgrd[n+1]:
            mask_auxgrd[n,j,i] = True
    utils_spec.ProgBar('done')

    if savevar == True:
      with open('variables/mask_auxgrd', 'w') as f: pickle.dump(mask_auxgrd, f) 	# save variable

    return(mask_auxgrd)

# --------------------------------------------------
# generate iter_maskcombo
# --------------------------------------------------
def gen_iter_maskcombo(lat_auxgrd, MW, mask_auxgrd, mask_modgrd, savevar=True):
    ''' Array of size (nlatAUXgrid, nlatMODELgrid)
        Each element is a np.Array containing longitude-indices (i) where 
	both, mask_auxgrd and the region mask on the model grid are True. 
    '''
    print('> generating iter_maskcombo')
    lat_MW, iter_lat_auxgrd, iter_lat_MW, iter_lon_MW = vars2speedup(lat_auxgrd, MW) 	# np-arrays for speed    
    iter_maskcombo = np.zeros([len(lat_auxgrd), len(MW.nlat)], dtype=object)  		# pre-allocation as integer

    for n in iter_lat_auxgrd:
      utils_spec.ProgBar('step', barlen=30, step=n, nsteps=len(iter_lat_auxgrd)) 	# initialize and update progress bar	
      for j in iter_lat_MW:
        iter_maskcombo[n,j] = np.where((mask_auxgrd[n,j,:]) & (mask_modgrd[j,:]>=6))[0]
    utils_spec.ProgBar('done')

    if savevar == True:
      with open('variables/iter_maskcombo', 'w') as f: pickle.dump(iter_maskcombo, f) 	# save variable

    return(iter_maskcombo)

    ''' Run following lines for visual testing:

    a = np.zeros([180,384])
    b = np.zeros([180,320])
    for nn in np.arange(180):
      for jj in np.arange(384):
        if iter_maskcombo[nn,jj].shape[0] > 0:
          a[nn,jj] = 1
	  for ii in iter_maskcombo[nn,jj]:
            b[nn,ii] = 1
    plt.figure()
    plt.pcolor(b)
    '''
# --------------------------------------------------
# generate maxiter_depth, an array for seafloor detection on model grid
# --------------------------------------------------
def gen_maxiter_depth(lat_auxgrd, z_w_auxgrd, MW, seafloor_HT, savevar=True):
    ''' Array of size (nlatMODELgrid, nlonMODELgrid)
        It contains the indices of maximal depth (model T-grid) in order to stop k-iteration at the seafloor
     
    Comments:
      > For future generalization, think about, whether function argument z_w_auxgrd shall really 
        run on aux-grid or not on model grid, as they may differ from each other.
      > Make function lat_auxgrd independent
    '''
    print('> generating maxiter_depth')
    lat_MW, iter_lat_auxgrd, iter_lat_MW, iter_lon_MW = vars2speedup(lat_auxgrd, MW) 	# np-arrays for speed    
    maxiter_depth = np.zeros([len(MW.nlat), len(MW.nlon)], dtype=int)

    for j in iter_lat_MW:
      utils_spec.ProgBar('step', barlen=32, step=j, nsteps=len(iter_lat_MW))
      for i in iter_lon_MW:
        maxiter_depth[j,i] = np.where(z_w_auxgrd <= seafloor_HT[j,i])[-1][-1] 	# index of maximal depth at j,i
    utils_spec.ProgBar('done')

    if savevar == True:
      with open('variables/maxiter_depth', 'w') as f: pickle.dump(maxiter_depth, f) 	# save variable

    return(maxiter_depth)
