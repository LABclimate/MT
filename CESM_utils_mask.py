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
    elif name == 'lat99model':     # as in ncdat.lat_aux_grid but only every other entry
      lat = ncdat.MOC.lat_aux_grid[::4].values  # latitudes
    return(lat)

# --------------------------------------------------
# generate mask_auxgrd_overlay_lat, a mask for grid-overlay
# --------------------------------------------------
def gen_mask_auxgrd_overlay_lat_old(lat_auxgrd, ncdat):
    ''' Boolean of size (nlatAUXgrid, nlatMODELgrid, nlonMODELgrid)
        It is True where latitudes of auxillary and model grid lie in the same box. 
    '''
    print('> generating mask_auxgrd_overlay_lat')
    lat_mgrdT, iter_lat_auxgrd, iter_lat_mgrdT, iter_lon_mgrdT = vars2speedup(lat_auxgrd, ncdat) 	# np-arrays for speed
    mask_auxgrd_overlay_lat = np.zeros([len(lat_auxgrd), len(ncdat.nlat), len(ncdat.nlon)],dtype=bool) 	# pre-allocation as False
    for n in iter_lat_auxgrd[:-1]:
      utils_misc.ProgBar('step', step=n, nsteps=len(iter_lat_auxgrd), minbarlen=60)
      for j in iter_lat_mgrdT:
        for i in iter_lon_mgrdT:
          if lat_auxgrd[n] <= lat_mgrdT[j,i] < lat_auxgrd[n+1]:
            mask_auxgrd_overlay_lat[n,j,i] = True
    utils_misc.ProgBar('done')

    return(mask_auxgrd_overlay_lat)

# --------------------------------------------------
# generate mask_auxgrd_overlay_lat, a mask for grid-overlay
# --------------------------------------------------
def gen_mask_auxgrd_overlay_lat(auxLAT, ncdat):
    ''' Boolean of size (nlatAUXgrid, nlatMODELgrid, nlonMODELgrid)
        It is True where latitudes of auxillary and model grid lie in the same box. 
    '''
    print('> generating mask_auxgrd_overlay_lat')
    
    TLAT = np.array(ncdat.TLAT)     # mgrd
    TLONG = np.array(ncdat.TLONG)   # mgrd
    ULAT = np.array(ncdat.ULAT)     # mgrd
    COSULAT = np.cos(np.deg2rad(ULAT))          # cos of ULAT
    SINULAT = np.sin(np.deg2rad(ULAT))          # sin of ULAT
    ULONG = np.array(ncdat.ULONG)   # mgrd
    auxLAT = np.arange(-90,90, .1)  #!
    COSauxLAT = np.cos(np.deg2rad(auxLAT))      # cos of auxLAT
    SINauxLAT = np.sin(np.deg2rad(auxLAT))      # sin of auxLAT
    iter_auxLAT = np.arange(len(auxLAT))
    iter_lat_mgrdT = np.arange(len(ncdat.nlat))
    iter_lon_mgrdT = np.array(ncdat.nlon) 			# in normal order!!   

    def dist_flat(COSLAT, SINULAT, LONG):
        ''' formula found on: https://www.math.ksu.edu/~dbski/writings/haversine.pdf'''
        R = 6378000 #! R_earth in [m]
        return(R*np.sqrt(2-2*np.prod(COSLAT)*np.cos(np.deg2rad(np.diff(LONG)))-2*np.prod(SINLAT)))
    
    def area_heron(d):
        ''' formula found on: https://de.wikipedia.org/wiki/Dreiecksfl%C3%A4che'''
        s = np.sum(d)/2
        return(np.sqrt(s*(s-d[0])*(s-d[1])*(s-d[2])))
        
        
    def get_ji(idxclose, j, i):
        out = np.zeros(len(idxclose), dtype=object) # pre-allocation of output
        for ii in np.arange(len(idxclose)):
            if idxclose[ii] == 0: out[ii] = tuple([j-1, i-1])
            if idxclose[ii] == 1: out[ii] = tuple([j-1, i])
            if idxclose[ii] == 2: out[ii] = tuple([j, i-1])
            if idxclose[ii] == 3: out[ii] = tuple([j, i])
        return(out)
        
    def argsort_WtoE(idx):
        return(np.argsort([ULONG[idx[ii]] for ii in np.arange(len(idx))]))
        
    def area_1_triangle(COSLAT, SINLAT, LONG):
        d = np.zeros(3)     # pre-allocation of distances
        d[0] = dist_flat(COSLAT[[0,1]], SINLAT[[0,1]], LONG[[0,1]]) # in [m]
        d[1] = dist_flat(COSLAT[[1,2]], SINLAT[[1,2]], LONG[[1,2]]) # in [m]
        d[2] = dist_flat(COSLAT[[2,0]], SINLAT[[2,0]], LONG[[2,0]]) # in [m]
        return(area_heron(d))
        
    def area_2_triangles(COSLAT, SINLAT, LONG):
        d = np.zeros(5)     # pre-allocation of distances
        # first triangle                    
        d[0] = dist_flat(COSLAT[[0,1]], SINLAT[[0,1]], LONG[[0,1]]) # in [m]
        d[1] = dist_flat(COSLAT[[1,2]], SINLAT[[1,2]], LONG[[1,2]]) # in [m]
        d[2] = dist_flat(COSLAT[[2,0]], SINLAT[[2,0]], LONG[[2,0]]) # in [m]        
        # second triangle
        d[3] = dist_flat(COSLAT[[2,3]], SINLAT[[2,3]], LONG[[2,3]]) # in [m]
        d[4] = dist_flat(COSLAT[[0,3]], SINLAT[[0,3]], LONG[[0,3]]) # in [m]
        # areas
        AREAa = area_heron(d[[0,1,2]])
        AREAb = area_heron(d[[2,3,4]])
        return(AREAa + AREAb)
    def area_2_triangles(COSLAT, SINLAT, LONG):
        d = np.zeros(7)     # pre-allocation of distances
        # first triangle                    
        d[0] = dist_flat(COSLAT[[0,1]], SINLAT[[0,1]], LONG[[0,1]]) # in [m]
        d[1] = dist_flat(COSLAT[[1,2]], SINLAT[[1,2]], LONG[[1,2]]) # in [m]
        d[2] = dist_flat(COSLAT[[2,0]], SINLAT[[2,0]], LONG[[2,0]]) # in [m]        
        # second triangle
        d[3] = dist_flat(COSLAT[[2,3]], SINLAT[[2,3]], LONG[[2,3]]) # in [m]
        d[4] = dist_flat(COSLAT[[0,3]], SINLAT[[0,3]], LONG[[0,3]]) # in [m]
        # third triangle
        d[5] = dist_flat(COSLAT[[0,3]], SINLAT[[0,3]], LONG[[0,3]]) # in [m]
        d[6] = dist_flat(COSLAT[[0,3]], SINLAT[[0,3]], LONG[[0,3]]) # in [m]
        # areas
        AREAa = area_heron(d[[0,1,2]])
        AREAb = area_heron(d[[2,3,4]])
        return(AREAa + AREAb)
        
        
    for j in iter_lat_mgrdT[1:]:
        for i in iter_lon_mgrdT[1:]:auxLAT<MAXLAT and auxLAT>MINLAT
            closeULAT = ULAT[j-1:j+1, i-1:i+1].flatten()
            closeULONG = ULONG[j-1:j+1, i-1:i+1].flatten()
            MAXLAT = np.max(closeULAT)
            MINLAT = np.min(closeULAT)
            idxLATbtw = np.where((auxLAT<MAXLAT) & (auxLAT>MINLAT))[-1]
            LATbtw = auxLAT[idxLATbtw]
            try: LATbtwplusone = auxLAT[idxLATbtw+1]
            except: LATbtwplusone = np.hstack(auxLAT[idxLATbtw], 100) # buffer value for overshoot at upper border of lat_auxgrd
            nLATbtw = len(idxLATbtw)
            AREA = np.zeros(nLATbtw+1) # pre-allocation of AREAfractions
            
            ## --- only in one LATbin ---
            if nLATbtw == 0:
                LATbin = np.where(auxLAT<MINLAT)[-1][-1]
                frac[j,i,LATbin] = 1
                continue
            
            ## --- multiple LATbins ---
            
            # longitueds of intersects
            nLATfloor = np.sum(closeULAT<LATbtw[0])
            nLATtop = np.sum(closeULAT>LATbtw[-1])
            
            LONGbtwW = closeULONG[0] + np.diff(closeULONG[[0,2]])/np.diff(closeULAT[[0,2]]) * (LATbtw-closeULAT[0])
            LONGbtwE = closeULONG[1] + np.diff(closeULONG[[1,3]])/np.diff(closeULAT[[1,3]]) * (LATbtw-closeULAT[1])
            
            # lower part
            idxLATbelow = get_ji(np.where(closeULAT<LATbtw[0])[-1], j, i)
            nbelow = len(idxLATbelow)
            if nbelow == 0:
                print('nbelow = 0!! What next?')
                raise ValueError('Error')
            elif nbelow == 1:
                COSLAT, SINLAT, LONG= np.zeros(3), np.zeros(3), np.zeros(3)     # pre-allocation for points
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATbelow[0]], SINULAT[idxLATbelow[0]], ULONG[idxLATbelow[0]] # Lower corner
                COSLAT[1], SINLAT[1], LONG[1] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # West intersect
                COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # East intersect
                AREA[0] = area_1_triangle(COSLAT, SINLAT, LONG)
            elif nbelow == 2:
                COSLAT, SINLAT, LONG, = np.zeros(4), np.zeros(4), np.zeros(4)   # pre-allocation for points
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATbelow[0]], SINULAT[idxLATbelow[0]], ULONG[idxLATbelow[0]] # Lower corner West
                COSLAT[1], SINLAT[1], LONG[1] = COSULAT[idxLATbelow[1]], SINULAT[idxLATbelow[1]], ULONG[idxLATbelow[1]] # Lower corner East
                COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # West intersect
                COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # East intersect
                AREA[0] = area_2_triangles(COSLAT, SINLAT, LONG)
            elif nbelow == 3:
                COSLAT, SINLAT, LONG, = np.zeros(5), np.zeros(5), np.zeros(5)   # pre-allocation for points
                idx = argsort_WtoE(idxLATabove) # indices of indices of Lower corners, sorted from W to E
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATabove[idx[0]]], SINULAT[idxLATabove[idx[0]]], ULONG[idxLATabove[idx[0]]] # Lower corner West
                COSLAT[1], SINLAT[1], LONG[1] = COSULAT[idxLATabove[idx[1]]], SINULAT[idxLATabove[idx[1]]], ULONG[idxLATabove[idx[1]]] # Lower corner centre
                COSLAT[2], SINLAT[2], LONG[2] = COSULAT[idxLATabove[idx[2]]], SINULAT[idxLATabove[idx[[2]]], ULONG[idxLATabove[idx[[2]]] # Lower corner East
                COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # West intersect
                COSLAT[4], SINLAT[4], LONG[4] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # East intersect
                AREA[nLATbtw+1] = area_3_triangles(COSLAT, SINLAT, LONG)
   
            # centre part(s)      
            if nLATbtw>=2:
                for bb in np.arange(nLATbtw[:-1]): #!
                    idxLATabove = get_ji(np.where((closeULAT>LATbtw[bb]) & (closeULAT<LATbtwplusone[bb]))[-1], j, i)
                    nabove = len(idxLATabove)
                    if nabove == 0:
                        COSLAT, SINLAT, LONG, = np.zeros(4), np.zeros(4), np.zeros(4)   # pre-allocation for points
                        COSLAT[0], SINLAT[0], LONG[0] = COSauxLAT[idxLATbtw[bb+1]], SINauxLAT[idxLATbtw[bb+1]], LONGbtwW[bb+1]  # NW intersect
                        COSLAT[1], SINLAT[1], LONG[1] = COSauxLAT[idxLATbtw[bb+1]], SINauxLAT[idxLATbtw[bb+1]], LONGbtwE[bb+1]  # NE intersect
                        COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[bb]], SINauxLAT[idxLATbtw[bb]], LONGbtwW[bb]        # SW intersect
                        COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[bb], SINauxLAT[idxLATbtw[bb]], LONGbtwE[bb]         # SE intersect
                        AREA[bb+1] = area_2_triangles(COSLAT, SINLAT, LONG)
                    elif nabove == 1:
                        COSLAT, SINLAT, LONG, = np.zeros(5), np.zeros(5), np.zeros(5)   # pre-allocation for points
                        COSLAT[0], SINLAT[0], LONG[0] = COSauxLAT[idxLATbtw[bb+1]], SINauxLAT[idxLATbtw[bb+1]], LONGbtwW[bb+1]  # NW intersect
                        COSLAT[1], SINLAT[1], LONG[1] = COSauxLAT[idxLATbtw[bb+1]], SINauxLAT[idxLATbtw[bb+1]], LONGbtwE[bb+1]  # NE intersect
                        COSLAT[2], SINLAT[2], LONG[2] = COSULAT[idxLATabove[0]], SINULAT[idxLATabove[0]], ULONG[idxLATabove[0]] # Top corner West
                        COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[bb]], SINauxLAT[idxLATbtw[bb]], LONGbtwW[bb]        # SW intersect
                        COSLAT[4], SINLAT[4], LONG[4] = COSauxLAT[idxLATbtw[bb], SINauxLAT[idxLATbtw[bb]], LONGbtwE[bb]         # SE intersect
                        AREA[bb+1] = area_3_triangles(COSLAT, SINLAT, LONG)
                    elif nabove == 2:
                        
            # upper part
            idxLATabove = get_ji(np.where((closeULAT>LATbtw[-1]) & (closeULAT<LATbtwplusone[-1]))[-1], j, i)
            nabove = len(idxLATabove)
            if nabove == 0:
                print('nbelow = 0!! What next?')
                raise ValueError('Error')
            elif nabove == 1:
                COSLAT, SINLAT, LONG= np.zeros(3), np.zeros(3), np.zeros(3)     # pre-allocation for points
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATabove[0]], SINULAT[idxLATabove[0]], ULONG[idxLATabove[0]] # Upper corner
                COSLAT[1], SINLAT[1], LONG[1] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # West intersect
                COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # East intersect
                AREA[nLATbtw+1] = area_1_triangle(COSLAT, SINLAT, LONG)
            elif nabove == 2:
                COSLAT, SINLAT, LONG, = np.zeros(4), np.zeros(4), np.zeros(4)   # pre-allocation for points
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATabove[0]], SINULAT[idxLATabove[0]], ULONG[idxLATabove[0]] # Upper corner West
                COSLAT[1], SINLAT[1], LONG[1] = COSULAT[idxLATabove[1]], SINULAT[idxLATabove[1]], ULONG[idxLATabove[1]] # Upper corner East
                COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # West intersect
                COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # East intersect
                AREA[nLATbtw+1] = area_2_triangles(COSLAT, SINLAT, LONG)
            elif nabove == 3:
                COSLAT, SINLAT, LONG, = np.zeros(5), np.zeros(5), np.zeros(5)   # pre-allocation for points
                idx = argsort_WtoE(idxLATabove) # indices of indices of Upper corners, sorted from W to E
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATabove[idx[0]]], SINULAT[idxLATabove[idx[0]]], ULONG[idxLATabove[idx[0]]] # Upper corner West
                COSLAT[1], SINLAT[1], LONG[1] = COSULAT[idxLATabove[idx[1]]], SINULAT[idxLATabove[idx[1]]], ULONG[idxLATabove[idx[1]]] # Upper corner centre
                COSLAT[2], SINLAT[2], LONG[2] = COSULAT[idxLATabove[idx[2]]], SINULAT[idxLATabove[idx[2]]], ULONG[idxLATabove[idx[2]]] # Upper corner East
                COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # West intersect
                COSLAT[4], SINLAT[4], LONG[4] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # East intersect
                AREA[nLATbtw+1] = area_3_triangles(COSLAT, SINLAT, LONG)

            LATbins = np.hstack([LATbtw[0]-1, LATbtw])
            idxMAXLAT = np.argwhere(closeULAT == MAXLAT)
            idxMINLAT = np.argwhere(closeULAT == MINLAT)
            if len(idxMAXLAT.flatten()) == 4: idxMAXLAT = idxMAXLAT[0] # at low latitudes latvalues at different i are the same
            if len(idxMINLAT.flatten()) == 4: idxMINLAT = idxMINLAT[0] # at low latitudes latvalues at different i are the same
            P1 = np.array([MAXLAT, closeULONG[tuple(idxMAXLAT)]])
            P2 = 
            while True:
                MAX                
                    
    mask_auxgrd_overlay_lat = np.zeros([len(auxLAT), len(ncdat.nlat), len(ncdat.nlon)],dtype=bool) 	# pre-allocation as False
    for n in iter_auxLAT[:-1]:
      utils_misc.ProgBar('step', step=n, nsteps=len(iter_auxLAT), minbarlen=60)
      for j in iter_lat_mgrdT:
        for i in iter_lon_mgrdT:
          if auxLAT[n] <= TLAT[j,i] < auxLAT[n+1]:
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
