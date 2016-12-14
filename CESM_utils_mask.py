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

def get_ATLiter(mask):
    return(np.array([np.where(mask[ii,:]) for ii in np.arange(mask.shape[0])]).squeeze())


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
# generate mask_auxgrd_overlay_lat, the mask for grid-overlay
# --------------------------------------------------
def gen_mask_auxgrd_overlay_lat(lat_auxgrd, ncdat):
    ''' Boolean of size (nlatAUXgrid, nlatMODELgrid, nlonMODELgrid)
        It is True where latitudes of auxillary and model grid lie in the same box. 
    '''
    print('> generating mask_auxgrd_overlay_lat')
    lat_mgrdT, iter_lat_auxgrd, iter_lat_mgrdT, iter_lon_mgrdT = vars2speedup(lat_auxgrd, ncdat) 	# np-arrays for speed
    lat_auxgrd = np.array(lat_auxgrd) # for speed
    mask_auxgrd_overlay_lat = np.zeros([len(lat_auxgrd), len(ncdat.nlat), len(ncdat.nlon)],dtype=bool) 	# pre-allocation as False
    for n in iter_lat_auxgrd[:-1]:
      utils_misc.ProgBar('step', step=n, nsteps=len(iter_lat_auxgrd), minbarlen=60)
      for j in iter_lat_mgrdT:
          good_i = np.where(np.greater_equal(lat_mgrdT[j,:], lat_auxgrd[n]) & np.less(lat_mgrdT[j,:], lat_auxgrd[n+1]))
          for i in good_i:
            mask_auxgrd_overlay_lat[n,j,i] = True
    utils_misc.ProgBar('done')

    return(mask_auxgrd_overlay_lat)

# --------------------------------------------------
# generate fraction_mask, the mask for grid-overlay
# --------------------------------------------------
def gen_fraction_mask(auxLAT, ncdat):
    ''' Boolean of size (nlatAUXgrid, nlatMODELgrid, nlonMODELgrid)
        It is True where latitudes of auxillary and model grid lie in the same box. 
        
        Important Note: check all border values! sometimes I cropped the loops 
                        in order not to get missing values, e.g. as U-grid and 
                        T-grid have the same shape
    '''
    #auxLAT = np.arange(-90,90, .1)
    
    print('> generating mask_auxgrd_overlay_lat')
    # -----------------------------------------------------------------
    # - Pre-allocation of mask
    fraction_mask = np.zeros([len(auxLAT)], dtype=object) # second dimension is going to be expanded below if needed
    for ii in np.arange(len(auxLAT)):
        fraction_mask[ii] = np.nan*np.ones(3) # dummies, which will be deleted afterwards
    stats = dict()
    stats['cases'] = np.zeros(10)
    stats['nLATbtw'] = np.nan * np.ones_like(ncdat.TLAT)
    stats['area'] = np.zeros([len(ncdat.nlat), len(ncdat.nlon), 3])
    stats['area'][:,:,0] = ncdat.TAREA
    # -----------------------------------------------------------------
    # - LATITUDES and LONGITUES
    TLAT = np.array(ncdat.TLAT)                 # mgrd
    TLONG = np.array(ncdat.TLONG)               # mgrd
    ULAT = np.array(ncdat.ULAT)                 # mgrd
    ULONG = np.array(ncdat.ULONG)               # mgrd
    COSULAT = np.cos(np.deg2rad(ULAT))          # cos of ULAT
    SINULAT = np.sin(np.deg2rad(ULAT))          # sin of ULAT
    COSauxLAT = np.cos(np.deg2rad(auxLAT))      # cos of auxLAT
    SINauxLAT = np.sin(np.deg2rad(auxLAT))      # sin of auxLAT
    
    # -----------------------------------------------------------------
    # - FUNCTIONS
    def dist_flat(COSLAT, SINULAT, LONG):
        ''' formula found on: https://www.math.ksu.edu/~dbski/writings/haversine.pdf'''
        ''' Radius_Earth is mean of Wikipedia-Values for R_pol and R_eq'''
        R = np.mean([6356752, 6371008]) # mean Earth radius in [m]
        return(R*np.sqrt(2-2*np.prod(COSLAT)*np.cos(np.deg2rad(np.diff(LONG)))-2*np.prod(SINLAT)))
    
    def area_heron(d):
        ''' formula found on: https://de.wikipedia.org/wiki/Dreiecksfl%C3%A4che'''
        s = np.sum(d)/2
        return(np.sqrt(s*(s-d[0])*(s-d[1])*(s-d[2])))
        
        
    def get_ji(idxji, j, i):
        out = np.zeros(len(idxji), dtype=object) # pre-allocation of output
        for ii in np.arange(len(idxji)):
            if idxji[ii] == 0: out[ii] = tuple([j-1, i-1])
            if idxji[ii] == 1: out[ii] = tuple([j-1, i])
            if idxji[ii] == 2: out[ii] = tuple([j, i-1])
            if idxji[ii] == 3: out[ii] = tuple([j, i])
          
        return(out)
        
    def argsort_WtoE(idx):
        return(np.argsort([ULONG[idx[ii]] for ii in np.arange(len(idx))]))
        
    def area_1_triangle(COSLAT, SINLAT, LONG):
        d = np.zeros(3)     # pre-allocation of distances
        # legs in m
        d[0] = dist_flat(COSLAT[[0,1]], SINLAT[[0,1]], LONG[[0,1]])
        d[1] = dist_flat(COSLAT[[1,2]], SINLAT[[1,2]], LONG[[1,2]])
        d[2] = dist_flat(COSLAT[[2,0]], SINLAT[[2,0]], LONG[[2,0]])
        # area
        return(area_heron(d))
        
    def area_2_triangles(COSLAT, SINLAT, LONG):
        d = np.zeros(5)     # pre-allocation of distances
        # legs in m       
        d[0] = dist_flat(COSLAT[[0,1]], SINLAT[[0,1]], LONG[[0,1]])
        d[1] = dist_flat(COSLAT[[1,2]], SINLAT[[1,2]], LONG[[1,2]])
        d[2] = dist_flat(COSLAT[[2,0]], SINLAT[[2,0]], LONG[[2,0]])        
        d[3] = dist_flat(COSLAT[[2,3]], SINLAT[[2,3]], LONG[[2,3]])
        d[4] = dist_flat(COSLAT[[0,3]], SINLAT[[0,3]], LONG[[0,3]])
        # areas
        AREAa = area_heron(d[[0,1,2]])
        AREAb = area_heron(d[[2,3,4]])
        return(AREAa + AREAb)
        
    def area_3_triangles(COSLAT, SINLAT, LONG):
        d = np.zeros(7)     # pre-allocation of distances
        # legs in m
        d[0] = dist_flat(COSLAT[[0,1]], SINLAT[[0,1]], LONG[[0,1]])
        d[1] = dist_flat(COSLAT[[1,2]], SINLAT[[1,2]], LONG[[1,2]])
        d[2] = dist_flat(COSLAT[[2,0]], SINLAT[[2,0]], LONG[[2,0]])        
        d[3] = dist_flat(COSLAT[[2,3]], SINLAT[[2,3]], LONG[[2,3]])
        d[4] = dist_flat(COSLAT[[0,3]], SINLAT[[0,3]], LONG[[0,3]])
        d[5] = dist_flat(COSLAT[[3,4]], SINLAT[[3,4]], LONG[[3,4]])
        d[6] = dist_flat(COSLAT[[0,4]], SINLAT[[0,4]], LONG[[0,4]])
        # areas
        AREAa = area_heron(d[[0,1,2]])
        AREAb = area_heron(d[[2,3,4]])
        AREAc = area_heron(d[[4,5,6]])
        return(AREAa + AREAb + AREAc)
        
    # -----------------------------------------------------------------
    # - LOOP over mgrd        
    for j in np.arange(1,len(ncdat.nlat)):
        utils_misc.ProgBar('step', step=j, nsteps=len(ncdat.nlat), minbarlen=60, forceinit = True)
        for i in np.arange(1,len(ncdat.nlon)):
            ULATji = ULAT[j-1:j+1, i-1:i+1].flatten()   # flat and thus allocatable with [0,1,2,3]
            ULONGji = ULONG[j-1:j+1, i-1:i+1].flatten() # flat and thus allocatable with [0,1,2,3]
            MAXLAT = np.max(ULATji)
            MINLAT = np.min(ULATji)
            idxLATbtw = np.where((auxLAT<MAXLAT) & (auxLAT>MINLAT))[-1] # indices of auxLAT between min, max corner of ULATji
            LATbtw = auxLAT[idxLATbtw]
            nLATbtw = len(idxLATbtw)
            try: LATbtwplusone = auxLAT[idxLATbtw+1]
            except: LATbtwplusone = np.hstack((auxLAT[idxLATbtw], 100)) #! buffer value for overshoot at upper border of lat_auxgrd
            AREA = np.zeros(nLATbtw+1) # pre-allocation of AREA
            
            stats['nLATbtw'][j,i] = nLATbtw

            ## --- only in one LATbin -----------------------------------------
            if nLATbtw == 0:
                stats['cases'][0] += 1
                idx = np.where(auxLAT<MINLAT)[-1][-1]
                fraction_mask[idx] = np.vstack((fraction_mask[idx], np.array([j,i,1])))
                # -----------------------------------------------------------------
                # - calculation of areas for check
                COSLAT, SINLAT, LONG, = np.zeros(4), np.zeros(4), np.zeros(4)   # pre-allocation for points
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[j-1,i-1], SINULAT[j-1,i-1], ULONG[j-1,i-1] # SW corner
                COSLAT[1], SINLAT[1], LONG[1] = COSULAT[j-1,i], SINULAT[j-1,i], ULONG[j-1,i] # SE corner
                COSLAT[2], SINLAT[2], LONG[2] = COSULAT[j,i], SINULAT[j,i], ULONG[j,i] # NE corner
                COSLAT[3], SINLAT[3], LONG[3] = COSULAT[j,i-1], SINULAT[j,i-1], ULONG[j,i-1] # NW corner
                A = area_2_triangles(COSLAT, SINLAT, LONG)
                # - Comparison with TAREA from ncdat
                stats['area'][j,i,1] = A
                stats['area'][j,i,2] = np.diff([A, ncdat.TAREA[j,i]]) 
                
                continue

            
            ## --- multiple LATbins -------------------------------------------
            # -----------------------------------------------------------------
            # - LONGITUDES of intersects  #!THCHK
            LONGbtwW, LONGbtwE = np.zeros(nLATbtw), np.zeros(nLATbtw)
            # lowest values
            nLATfloor = np.sum(ULATji<LATbtw[0]) # number of corners below lowest LATbtw
            if nLATfloor == 1 and ULATji[0]<LATbtw[0]: # W corner is lowest
                LONGbtwW[0] = ULONGji[0] + np.diff(ULONGji[[0,2]])/np.diff(ULATji[[0,2]]) * (LATbtw[0]-ULATji[0])
                LONGbtwE[0] = ULONGji[0] + np.diff(ULONGji[[0,1]])/np.diff(ULATji[[0,1]]) * (LATbtw[0]-ULATji[0])
            elif nLATfloor == 1 and ULATji[1]<LATbtw[0]: # E corner is lowest
                LONGbtwW[0] = ULONGji[1] + np.diff(ULONGji[[1,0]])/np.diff(ULATji[[1,0]]) * (LATbtw[0]-ULATji[1])
                LONGbtwE[0] = ULONGji[1] + np.diff(ULONGji[[1,3]])/np.diff(ULATji[[1,3]]) * (LATbtw[0]-ULATji[1])
            elif nLATfloor == 2:
                LONGbtwW[0] = ULONGji[0] + np.diff(ULONGji[[0,2]])/np.diff(ULATji[[0,2]]) * (LATbtw[0]-ULATji[0])
                LONGbtwE[0] = ULONGji[1] + np.diff(ULONGji[[1,3]])/np.diff(ULATji[[1,3]]) * (LATbtw[0]-ULATji[1])
            # highest values
            nLATtop = np.sum(ULATji>LATbtw[-1]) # number of corners above highest LATbtw
            if nLATtop == 1 and ULATji[2]>LATbtw[-1]: # W corner is highest
                LONGbtwW[-1] = ULONGji[2] + np.diff(ULONGji[[2,0]])/np.diff(ULATji[[2,0]]) * (LATbtw[-1]-ULATji[2])
                LONGbtwE[-1] = ULONGji[2] + np.diff(ULONGji[[2,3]])/np.diff(ULATji[[2,3]]) * (LATbtw[-1]-ULATji[2])
            elif nLATtop == 1 and ULATji[3]>LATbtw[-1]: # E corner is highest
                LONGbtwW[-1] = ULONGji[3] + np.diff(ULONGji[[3,2]])/np.diff(ULATji[[3,2]]) * (LATbtw[-1]-ULATji[3])
                LONGbtwE[-1] = ULONGji[3] + np.diff(ULONGji[[3,1]])/np.diff(ULATji[[3,1]]) * (LATbtw[-1]-ULATji[3])
            elif nLATtop == 2:
                LONGbtwW[-1] = ULONGji[2] + np.diff(ULONGji[[2,0]])/np.diff(ULATji[[2,0]]) * (LATbtw[-1]-ULATji[2])
                LONGbtwE[-1] = ULONGji[3] + np.diff(ULONGji[[3,1]])/np.diff(ULATji[[3,1]]) * (LATbtw[-1]-ULATji[3])
            # eventual middle values
            if nLATbtw >=3:
                LONGbtwW[1:-1] = ULONGji[0] + np.diff(ULONGji[[0,2]])/np.diff(ULATji[[0,2]]) * (LATbtw[1:-1]-ULATji[0])
                LONGbtwE[1:-1] = ULONGji[1] + np.diff(ULONGji[[1,3]])/np.diff(ULATji[[1,3]]) * (LATbtw[1:-1]-ULATji[1])           
            
            # -----------------------------------------------------------------
            # - AREA calculation...
            # ... of lower part
            idxLATbelow = get_ji(np.where(ULATji<LATbtw[0])[-1], j, i)
            nbelow = len(idxLATbelow)
            if nbelow == 0:
                print('nbelow = 0!! What next?')
                raise ValueError('Error')
            elif nbelow == 1: # Case 7
                stats['cases'][7] += 1
                COSLAT, SINLAT, LONG= np.zeros(3), np.zeros(3), np.zeros(3)     # pre-allocation for points
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATbelow[0]], SINULAT[idxLATbelow[0]], ULONG[idxLATbelow[0]] # corner
                COSLAT[1], SINLAT[1], LONG[1] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # E intersect
                COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # W intersect
                AREA[0] = area_1_triangle(COSLAT, SINLAT, LONG)
            elif nbelow == 2: # Case 8
                stats['cases'][8] += 1
                COSLAT, SINLAT, LONG, = np.zeros(4), np.zeros(4), np.zeros(4)   # pre-allocation for points
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATbelow[0]], SINULAT[idxLATbelow[0]], ULONG[idxLATbelow[0]] # W corner
                COSLAT[1], SINLAT[1], LONG[1] = COSULAT[idxLATbelow[1]], SINULAT[idxLATbelow[1]], ULONG[idxLATbelow[1]] # E corner
                COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # E intersect
                COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # W intersect
                AREA[0] = area_2_triangles(COSLAT, SINLAT, LONG)
            elif nbelow == 3: # Case 9
                stats['cases'][9] += 1
                COSLAT, SINLAT, LONG, = np.zeros(5), np.zeros(5), np.zeros(5)   # pre-allocation for points
                idx = argsort_WtoE(idxLATbelow) # indices of indices of Lower corners, sorted from W to E
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATbelow[idx[0]]], SINULAT[idxLATbelow[idx[0]]], ULONG[idxLATbelow[idx[0]]] # W corner
                COSLAT[1], SINLAT[1], LONG[1] = COSULAT[idxLATbelow[idx[1]]], SINULAT[idxLATbelow[idx[1]]], ULONG[idxLATbelow[idx[1]]] # C corner
                COSLAT[2], SINLAT[2], LONG[2] = COSULAT[idxLATbelow[idx[2]]], SINULAT[idxLATbelow[idx[2]]], ULONG[idxLATbelow[idx[2]]] # E corner
                COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # E intersect
                COSLAT[4], SINLAT[4], LONG[4] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # W intersect
                AREA[0] = area_3_triangles(COSLAT, SINLAT, LONG)
   
            # ... of eventual centre part(s)      
            if nLATbtw>=2:
                for bb in np.arange(1,nLATbtw):
                    idxLATabove = get_ji(np.where((ULATji>LATbtw[bb-1]) & (ULATji<=LATbtwplusone[bb-1]))[-1], j, i)
                    nabove = len(idxLATabove)
                    if nabove == 0: # Case 4
                        stats['cases'][4] += 1
                        COSLAT, SINLAT, LONG, = np.zeros(4), np.zeros(4), np.zeros(4)   # pre-allocation for points
                        COSLAT[0], SINLAT[0], LONG[0] = COSauxLAT[idxLATbtw[bb]], SINauxLAT[idxLATbtw[bb]], LONGbtwW[bb]  # NW intersect
                        COSLAT[1], SINLAT[1], LONG[1] = COSauxLAT[idxLATbtw[bb]], SINauxLAT[idxLATbtw[bb]], LONGbtwE[bb]  # NE intersect
                        COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[bb-1]], SINauxLAT[idxLATbtw[bb-1]], LONGbtwW[bb-1]        # SW intersect
                        COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[bb-1]], SINauxLAT[idxLATbtw[bb-1]], LONGbtwE[bb-1]         # SE intersect
                        AREA[bb] = area_2_triangles(COSLAT, SINLAT, LONG)
                    elif nabove == 1 and ULONG[idxLATabove[0]]<LONGbtwW[bb-1]: # Case 5a
                        stats['cases'][5] += 1
                        COSLAT, SINLAT, LONG, = np.zeros(5), np.zeros(5), np.zeros(5)   # pre-allocation for points
                        COSLAT[0], SINLAT[0], LONG[0] = COSauxLAT[idxLATbtw[bb]], SINauxLAT[idxLATbtw[bb]], LONGbtwE[bb]  # NE intersect
                        COSLAT[1], SINLAT[1], LONG[1] = COSauxLAT[idxLATbtw[bb]], SINauxLAT[idxLATbtw[bb]], LONGbtwW[bb]  # NW intersect
                        COSLAT[2], SINLAT[2], LONG[2] = COSULAT[idxLATabove[0]], SINULAT[idxLATabove[0]], ULONG[idxLATabove[0]] # corner
                        COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[bb-1]], SINauxLAT[idxLATbtw[bb-1]], LONGbtwW[bb-1]        # SW intersect
                        COSLAT[4], SINLAT[4], LONG[4] = COSauxLAT[idxLATbtw[bb-1]], SINauxLAT[idxLATbtw[bb-1]], LONGbtwE[bb-1]         # SE intersect
                        AREA[bb] = area_3_triangles(COSLAT, SINLAT, LONG)
                    elif nabove == 1 and ULONG[idxLATabove[0]]>LONGbtwE[bb-1]: # Case 5b
                        stats['cases'][5] += 1
                        COSLAT, SINLAT, LONG, = np.zeros(5), np.zeros(5), np.zeros(5)   # pre-allocation for points
                        COSLAT[0], SINLAT[0], LONG[0] = COSauxLAT[idxLATbtw[bb]], SINauxLAT[idxLATbtw[bb]], LONGbtwE[bb]  # NE intersect
                        COSLAT[1], SINLAT[1], LONG[1] = COSauxLAT[idxLATbtw[bb]], SINauxLAT[idxLATbtw[bb]], LONGbtwW[bb]  # NW intersect
                        COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[bb-1]], SINauxLAT[idxLATbtw[bb-1]], LONGbtwW[bb-1]        # SW intersect
                        COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[bb-1]], SINauxLAT[idxLATbtw[bb-1]], LONGbtwE[bb-1]         # SE intersect
                        COSLAT[4], SINLAT[4], LONG[4] = COSULAT[idxLATabove[0]], SINULAT[idxLATabove[0]], ULONG[idxLATabove[0]] # corner
                        AREA[bb] = area_3_triangles(COSLAT, SINLAT, LONG)                        
                    elif nabove == 2: # Case 6
                        stats['cases'][6] += 1
                        idx = argsort_WtoE(idxLATabove) # indices of indices of Upper corners, sorted from W to E
                        COSLAT, SINLAT, LONG, = np.zeros(6), np.zeros(6), np.zeros(6)   # pre-allocation for points
                        COSLAT[0], SINLAT[0], LONG[0] = COSauxLAT[idxLATbtw[bb]], SINauxLAT[idxLATbtw[bb]], LONGbtwE[bb]  # NE intersect
                        COSLAT[1], SINLAT[1], LONG[1] = COSauxLAT[idxLATbtw[bb]], SINauxLAT[idxLATbtw[bb]], LONGbtwW[bb]  # NW intersect
                        COSLAT[2], SINLAT[2], LONG[2] = COSULAT[idxLATabove[idx[0]]], SINULAT[idxLATabove[idx[0]]], ULONG[idxLATabove[idx[0]]] # W corner
                        COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[bb-1]], SINauxLAT[idxLATbtw[bb-1]], LONGbtwW[bb-1]        # SW intersect
                        COSLAT[4], SINLAT[4], LONG[4] = COSauxLAT[idxLATbtw[bb-1]]  , SINauxLAT[idxLATbtw[bb-1]], LONGbtwE[bb-1]         # SE intersect
                        COSLAT[5], SINLAT[5], LONG[5] = COSULAT[idxLATabove[idx[1]]], SINULAT[idxLATabove[idx[1]]], ULONG[idxLATabove[idx[1]]] # E corner
                        AREA[bb] = area_3_triangles(COSLAT, SINLAT, LONG)

            # ... of upper part
            idxLATabove = get_ji(np.where(ULATji>LATbtw[-1])[-1], j, i)
            nabove = len(idxLATabove)
            if nabove == 0:
                print('nbelow = 0!! What next?')
                raise ValueError('Error')
            elif nabove == 1: # Case 1
                stats['cases'][1] += 1                
                COSLAT, SINLAT, LONG= np.zeros(3), np.zeros(3), np.zeros(3)     # pre-allocation for points
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATabove[0]], SINULAT[idxLATabove[0]], ULONG[idxLATabove[0]] # corner
                COSLAT[1], SINLAT[1], LONG[1] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # W intersect
                COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # E intersect
                AREA[-1] = area_1_triangle(COSLAT, SINLAT, LONG)
            elif nabove == 2: # Case 2
                stats['cases'][2] += 1    
                COSLAT, SINLAT, LONG, = np.zeros(4), np.zeros(4), np.zeros(4)   # pre-allocation for points
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATabove[1]], SINULAT[idxLATabove[1]], ULONG[idxLATabove[1]] # E corner
                COSLAT[1], SINLAT[1], LONG[1] = COSULAT[idxLATabove[0]], SINULAT[idxLATabove[0]], ULONG[idxLATabove[0]] # W corner
                COSLAT[2], SINLAT[2], LONG[2] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # W intersect
                COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # E intersect
                AREA[-1] = area_2_triangles(COSLAT, SINLAT, LONG)
            elif nabove == 3: # Case 3
                stats['cases'][3] += 1       
                COSLAT, SINLAT, LONG, = np.zeros(5), np.zeros(5), np.zeros(5)   # pre-allocation for points
                idx = argsort_WtoE(idxLATabove) # indices of indices of Upper corners, sorted from W to E
                COSLAT[0], SINLAT[0], LONG[0] = COSULAT[idxLATabove[idx[2]]], SINULAT[idxLATabove[idx[2]]], ULONG[idxLATabove[idx[2]]] # Upper corner West
                COSLAT[1], SINLAT[1], LONG[1] = COSULAT[idxLATabove[idx[1]]], SINULAT[idxLATabove[idx[1]]], ULONG[idxLATabove[idx[1]]] # Upper corner centre
                COSLAT[2], SINLAT[2], LONG[2] = COSULAT[idxLATabove[idx[0]]], SINULAT[idxLATabove[idx[0]]], ULONG[idxLATabove[idx[0]]] # Upper corner East
                COSLAT[3], SINLAT[3], LONG[3] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwW[0]           # West intersect
                COSLAT[4], SINLAT[4], LONG[4] = COSauxLAT[idxLATbtw[0]], SINauxLAT[idxLATbtw[0]], LONGbtwE[0]           # East intersect
                AREA[-1] = area_3_triangles(COSLAT, SINLAT, LONG)
            
            # -----------------------------------------------------------------
            # - Fraction of areas
            fracAREA = AREA/np.sum(AREA)
            # - Comparison with TAREA from ncdat
            stats['area'][j,i,1] = np.sum(AREA)
            stats['area'][j,i,2] = np.diff([np.sum(AREA), ncdat.TAREA[j,i]])            
            # -----------------------------------------------------------------
            # - Write to fraction_mask
            for ii in np.arange(nLATbtw+1):
                if ii == 0: idx = idxLATbtw[0]-1
                else:       idx = idxLATbtw[ii-1]
                try: fraction_mask[idx] = np.vstack((fraction_mask[idx], np.array([j,i,fracAREA[ii]])))
                except: debug_here()
    
    # ---------------------------------------------------------------------             
    # delete dummy arrays
    for ii in np.arange(len(auxLAT)):
        fraction_mask[ii] = np.delete(fraction_mask[ii], 0,0)
    print('statistics cases: \n{} '.format(stats['cases']))
    # return(fraction_mask, stats)
    return(fraction_mask)

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
def calc_H_auxgrd_xmax(lat_auxgrd, ncdat, TorUgrid):
    # a few variables to speed up subsequent loops
    iter_lat_auxgrd = np.arange(len(lat_auxgrd))
    overlayboolmask = utils_mask.gen_mask_auxgrd_overlay_lat(lat_auxgrd, ncdat)
    iter_lat_mgrd = [np.where(np.any(overlayboolmask[n,:,:], axis=1))[0] for n in iter_lat_auxgrd]
    ATLboolmask = utils_mask.get_ATLbools(ncdat.REGION_MASK)    # boolean mask 
    mask_ATLiter = utils_mask.get_ATLiter(ATLboolmask)
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
      for j in iter_lat_mgrd[n]:
        for i in mask_ATLiter[j]:
          H_auxgrd_xmax[n] = np.nanmax([H_auxgrd_xmax[n], H[j,i]])
    utils_misc.ProgBar('done')

    return(H_auxgrd_xmax)