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
# - get_ATLbools()
# - vars2speedup()
# - gen_auxgrd()
# - gen_mask_grd_overlay_lat()
# - gen_iter_maskcombo()
# - calc_HT_mgrd_xmax()
# - calc_HT_auxgrd_xmax()
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
# 20-Mai-2016 - buerki@climate.unibe.ch : in gen_maxiter_depth() changed '<=' to '<' 
#                                         reason: <= diggs into ground for both z_w and z_t.
# 					  added implemented new utils_misc.savevar functions
# 24-Mai-2016 - buerki@climate.unibe.ch : changed MW to ncdat in various places
#                                         created calc_HT_mgrd_xmax()
#                                         created calc_HT_auxgrd_xmax()
# 31-Mai-2016 - buerki@climate.unibe.ch : in gen_maxiter_depth() changed '<' back to '<='
# 01-Jun-2016 - buerki@climate.unibe.ch : migrated gen_auxgrd from utils_MOC to utils_mask
# 15-Jun-2016 - buerki@climate.unibe.ch : added get_ATLbools()
# 16-Jun-2016 - buerki@climate.unibe.ch : removed gen_maxiter_depth()
#################################

import numpy as np
import numpy.ma as ma
import pickle
import CESM_utils_mask as utils_mask
import UTILS_misc as utils_misc
from IPython.core.debugger import Tracer; debug_here = Tracer()

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
#(mask>=6) & (mask!=7)  # ATL without Mediteranian
def mask_ATLANTIC(varin, mask, outputformat='xr'):
    if outputformat=='ma':
        return(ma.masked_array(varin, mask=np.array((mask>=6) & (mask!=7))))
    elif outputformat=='xr':
        return(varin.where((mask>=6) & (mask!=7)))

def get_ATLbools(mask):
    return(np.array((mask>=6) & (mask!=7)))

# =======================================================================================
# generate auxillary grid and related masks and mask-like iterators 
# =======================================================================================

# write some variables to numpy-arrays in order to speed up subsequent loops
def vars2speedup(lat_auxgrd, ncdat):
    lat_mgrdT = np.array(ncdat.TLAT)
    iter_lat_auxgrd = np.arange(len(lat_auxgrd))
    iter_lat_mgrdT = np.arange(len(ncdat.nlat))
    iter_lon_mgrdT = np.array(ncdat.nlon) 			# in normal order!!
    return(lat_mgrdT, iter_lat_auxgrd, iter_lat_mgrdT, iter_lon_mgrdT)

# --------------------------------------------------
# - generate auxillary grid
# --------------------------------------------------
def gen_auxgrd(ncdat, name):
    if name == 'lat170eq80S90N':    # 170 equally spaced boxes from 80S to 90N
      lat = np.linspace(-80, 90, 170)  	        # latitudes
    elif name == 'lat340eq80S90N':  # 340 equally spaced boxes from 80S to 90N
      lat = np.linspace(-80, 90, 340)  	        # latitudes
    elif name == 'lat198model':     # as in ncdat.lat_aux_grid but only every other entry
      lat = ncdat.MOC.lat_aux_grid[::2].values  # latitudes
    elif name == 'lat395model':     # as in ncdat.lat_aux_grid
      lat = ncdat.MOC.lat_aux_grid.values       # latitudes
    return(lat)

# --------------------------------------------------
# generate mask_auxgrd_overlay_lat, a mask for grid-overlay
# --------------------------------------------------
def gen_mask_auxgrd_overlay_lat(lat_auxgrd, ncdat):
    ''' Boolean of size (nlatAUXgrid, nlatMODELgrid, nlonMODELgrid)
        It is True where latitudes of auxillary and model grid lie in the same box. 
    '''
    print('> generating mask_auxgrd_overlay_lat')
    lat_mgrdT, iter_lat_auxgrd, iter_lat_mgrdT, iter_lon_mgrdT = vars2speedup(lat_auxgrd, ncdat) 	# np-arrays for speed
    mask_auxgrd_overlay_lat = np.zeros([len(lat_auxgrd), len(ncdat.nlat), len(ncdat.nlon)],dtype=bool) 	# pre-allocation as False

    for n in iter_lat_auxgrd:
      utils_misc.ProgBar('step', step=n, nsteps=len(iter_lat_auxgrd), minbarlen=60)
      for j in iter_lat_mgrdT:
        for i in iter_lon_mgrdT:
          if lat_auxgrd[n] <= lat_mgrdT[j,i] < lat_auxgrd[n+1]:
            mask_auxgrd_overlay_lat[n,j,i] = True
    utils_misc.ProgBar('done')

    return(mask_auxgrd_overlay_lat)

# --------------------------------------------------
# generate iter_maskcombo
# --------------------------------------------------
def gen_iter_maskcombo(lat_auxgrd, ncdat, mask_auxgrd_overlay_lat):
    ''' Array of size (nlatAUXgrid, nlatMODELgrid)
        Each element is a np.Array containing longitude-indices (i) where 
	both, mask_auxgrd_overlay_lat and the region mask on the model grid are True. 
    '''
    print('> generating iter_maskcombo')
    # np-arrays for speed
    lat_mgrdT, iter_lat_auxgrd, iter_lat_mgrdT, iter_lon_mgrdT = utils_mask.vars2speedup(lat_auxgrd, ncdat)
    mask_modgrd = ncdat.REGION_MASK.values
    # pre-allocation with zeros and dtype=int
    iter_maskcombo = np.zeros([len(lat_auxgrd), len(ncdat.nlat)], dtype=object)  

    for n in iter_lat_auxgrd:
      utils_misc.ProgBar('step', step=n, nsteps=len(iter_lat_auxgrd), minbarlen=60)
      for j in iter_lat_mgrdT:
        iter_maskcombo[n,j] = np.where((mask_auxgrd_overlay_lat[n,j,:]) & (mask_modgrd[j,:]>=6))[0]
    utils_misc.ProgBar('done')

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

# =======================================================================================
#  Find maixmal depth along longitudes for different grids
# =======================================================================================

# ...for model grid
def calc_H_mgrd_xmax(ncdat, TorUgrid):
    if TorUgrid == 'T':
      Hm = utils_mask.mask_ATLANTIC(ncdat.HT, ncdat.REGION_MASK) # mask Atlantic
      fname = 'HT_mgrd_xmax'
    elif TorUgrid == 'U':
      Hm = utils_mask.mask_ATLANTIC(ncdat.HU, ncdat.REGION_MASK) # mask Atlantic
      fname = 'HU_mgrd_xmax'

    H_mgrd_xmax = Hm.max(dim='nlon')    # find max along longitudes

    return(H_mgrd_xmax)

# ...for auxillary grid
def calc_H_auxgrd_xmax(lat_auxgrd, ncdat, TorUgrid, iter_maskcombo):
    # a few variables to speed up subsequent loops
    iter_lat_auxgrd = np.arange(len(lat_auxgrd))
    iter_lat_mgrd = np.arange(len(ncdat.nlat))

    # find maximal depth along longitudes
    if TorUgrid == 'T':
      H = ncdat.HT
      fname = 'HT_auxgrd_xmax'
    elif TorUgrid == 'U':
      H = ncdat.HU
      fname = 'HU_auxgrd_xmax'

    H_auxgrd_xmax = np.zeros(len(iter_lat_auxgrd))             # pre-allocation with zeros
    for n in iter_lat_auxgrd:
      utils_misc.ProgBar('step', step=n, nsteps=len(iter_lat_auxgrd))
      for j in iter_lat_mgrd:
        for i in iter_maskcombo[n,j]:
          H_auxgrd_xmax[n] = np.nanmax([H_auxgrd_xmax[n], H[j,i]])
    utils_misc.ProgBar('done')

    return(H_auxgrd_xmax)
