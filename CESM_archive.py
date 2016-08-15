#################################
# The CESM python toolbox at KUP
# -- Archive of old functions --
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_archive as utils_archive
#################################
# contained functions:
#################################
# - calc_dMOC_mgrd()
#################################
# please log your changes below:
#################################
# 15-Jun-2016 - buerki@climate.unibe.ch : created this archive
#                                         added calc_dMOC_mgrd()
#################################


# =======================================================================================
# - maxiter_depth filtering
# =======================================================================================
''' It was found out only after abandoning this rather slow looping method
    that np.arange crops last entry, therefore think about adding +1
'''

# in utils_mask:
# --------------

        # --------------------------------------------------
        # generate maxiter_depth, an array for seafloor detection on model grid
        # --------------------------------------------------
        def gen_maxiter_depth(lat_auxgrd, z_w_auxgrd, ncdat, path_vars, savevar=True):
            ''' Array of size (nlatMODELgrid, nlonMODELgrid)
                It contains the indices of maximal depth (model T-grid) in order to stop k-iteration at the seafloor
             
            Comments:
              > For future generalization, think about, whether function argument z_w_auxgrd shall really 
                run on aux-grid or not on model grid, as they may differ from each other.
              > Make function lat_auxgrd independent
            '''
            print('> generating maxiter_depth')
            # np-arrays for speed    
            lat_mgrdT, iter_lat_auxgrd, iter_lat_mgrdT, iter_lon_mgrdT = utils_mask.vars2speedup(lat_auxgrd, ncdat)
            HT = ncdat.HT.values
            # pre-allocation with zeros and dtype=object
            maxiter_depth = np.zeros([len(ncdat.nlat), len(ncdat.nlon)], dtype=object)
        
            for j in iter_lat_mgrdT:
              utils_misc.ProgBar('step', step=j, nsteps=len(iter_lat_mgrdT))
              for i in iter_lon_mgrdT:
                try:    maxiter_depth[j,i] = np.where(z_w_auxgrd <= HT[j,i])[-1][-1] 	# index of maximal depth at j,i
                except: maxiter_depth[j,i] = np.array([])
            utils_misc.ProgBar('done')
        
            if savevar == True:                                         # save to file
              utils_misc.savevar(maxiter_depth, path_vars+'maxiter_depth')
        
            return(maxiter_depth)

# in utils_MOC:
# -------------

            ''' > generate array (maxiter_depth) for seafloor detection. It contains the indices of maximal depth
                  for each point on model grid (T points) in order to stop k-iteration at the seafloor.:
            '''
            try:    maxiter_depth = utils_misc.loadvar(path_vars+'maxiter_depth') 
            except: maxiter_depth = utils_mask.gen_maxiter_depth(lat_auxgrd, z_auxgrd, ncdat, path_vars)
            # [...]
                for k in np.arange(int(maxiter_depth[j,i])):      # stop at depth of seafloor
                    Mxint[k,n] = np.nansum([Mxint[k,n],M[k,j,i]], axis=0)   # zonal integration
             
             
             
             
             
# =======================================================================================
# - dMOC on model grid
# =======================================================================================
def calc_dMOC_mgrd(vel_comp, M, PD, PD_bins, do_norm=True, dump_dMxint=False):
    '''
    Input:
     > vel_comp           : either 'W' or  'V' | string
     > M 			    : volume transport (MW or MV) | nparray of shape [nz, nlat, nlon]
     > PD                 : potential density | nparray of shape [nz, nlat, nlon]
     > PD_bins            : borders of PD-bins | nparray of shape [nPDbins+1]
     > do_norm 		    : boolean
     > dump_dMxint 	    : boolean
    Output:
     > dMxint             : zonally integrated volume transport of shape [nPDbins, nlat] | nparray
     > dMOC               : dMOC of shape [nPDbins, nlat] | nparray
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
# - resampling data on new grid using linear interpolation (along single dimension)
# =======================================================================================
def resample_colwise(odat, ogrd, ngrd, method, fill_value=np.nan, mask='none', sort_ogrd='True'):
    ''' 
    Uses:
     > resample_1dim_weightedmean()
    Input:
     > odat:        data on model grid
     > ogrd:        old grid
     > ngrd:        new grid 
     > method:      string | 'wmean' (weighted mean), 'dMW' or 'sum' (sum over all datapoints within bin on new grid)
     > fill_value:  float or nan | value used to fill in for requested points outside the range of ogrd.
     > mask:        mask of densityshape [j,i], default: all True (no mask)
     > sort_ogrd:   bool |  if True: ogrd (and odat) will be sorted such that ogrd is monotonically increasing (not necess. in strict sense!)
                            if False: ogrd will be manipulated such that strictly monotonically increasing (brute method, not recommended)
    Output:
     > ndat:    resampled data on new grid
    Comments:
     > add warning if negative gradients occur.
     > #!! for discreapencies at the low- and high-density borders see the comment in resample_1dim_weightedmean().
    '''

    def fill_gaps(ndat_ji, gaps_border, gaps_center, fill_value):
        ndat_ji[gaps_border] = fill_value
        ndat_ji[gaps_center] = fill_value
        return(ndat_ji)
        
    print(' > columnwise resampling')

    # shape of data-array
    if len(odat.shape)==3:
      len_j = odat.shape[-2] # assumed to be the second last entry
      len_i = odat.shape[-1] # assumed to be the last entry
    elif len(odat.shape)==2:
      len_j = 1              # set to one, such that j-loop is executed only once.
      len_i = odat.shape[-1] # assumed to be the last entry
    elif len(odat.shape)==1:
      len_j = 1              # set to one, such that j-loop is executed only once.
      len_i = 1              # set to one, such that i-loop is executed only once.

    # get default for regional mask (i.e. do not mask anything)
    if mask == 'none':
      mask = np.ones(shape=[len_j, len_i], dtype=bool)

    # expand the shape of ogrd and odat to 3 dimensions (singleton-dimensions are intended)
    ndim_ogrd = len(ogrd.shape)
    ndim_odat = len(odat.shape)
    if ndim_ogrd==1: # 1dim --> 3dim
      ogrd = utils_conv.expand_karray_to_kji(ogrd, len_j, len_i)
    elif ndim_ogrd==2: # 2dim --> 3dim
      sys.exit('case of two-dimensional ogrd is not implemented yet!')
      
    if ndim_odat==1: # 1dim --> 3dim
      odat = utils_conv.expand_karray_to_kji(odat, len_j, len_i)
    elif ndim_odat==2: # 2dim --> 3dim
      sys.exit('case of two-dimensional odat is not implemented yet!')
      
    # pre-allocation of ndat
    if method == 'wmean':
        ndat = fill_value * np.ones(shape=[len(ngrd), len_j, len_i])
    elif method == 'dMW':
        ndat = fill_value * np.ones(shape=[len(ngrd)-1, len_j, len_i])
        influx_highdens = np.zeros(shape=[len_j, len_i])
        
    # loop over columns
    for j in np.arange(len_j):
      utils_misc.ProgBar('step', step=j, nsteps=len_j)
      for i in np.arange(len_i):
        # skip masked [j,i]-tuples
        if mask[j,i]==False:
          continue
        # reduce ogrd and odat to current column
        ogrd_ji = ogrd[:,j,i]
        odat_ji = odat[:,j,i]
        # detect disjunct ogrd and ngrd and continue with next column
        if (np.nanmax(ogrd_ji) < np.nanmin(ngrd)) or (np.nanmax(ngrd) < np.nanmin(ogrd_ji)):
          print('disjunct ogrd and ngrd at (j,i)=({}, {}). (please check conservation of integrated flux!)'.format(j, i))
          continue
        # make ogrd strictly monotoneously increasing
        if any(np.diff(ogrd_ji)<=0):
          if sort_ogrd == 'True':
            # sort ogrd and odat (recommended)
            idx_sort = np.argsort(ogrd_ji)
            ogrd_ji = ogrd_ji[idx_sort]
            odat_ji = odat_ji[idx_sort]
          else:
            # brute method by manipulating ogrd (not recommended but used by Kwon)
            for k in np.arange(1,ogrd.shape[0]):
              if ogrd_ji[k] <= ogrd_ji[k-1]:
                ogrd_ji[k] = ogrd_ji[k-1]+1e-10
        
        # interpolation
        if method == 'wmean': # simple weighted mean interpolation
          ndat[:,j,i], gaps_border, gaps_center = resample_1dim_weightedmean(odat_ji, ogrd_ji, ngrd, fill_value)
          ndat[:,j,i] = fill_gaps(ndat[:,j,i], gaps_border, gaps_center, fill_value)
          
        elif method == 'dMW': # procedure for dMW
          # a) weighted mean of closest neighbours around dens_bin border values
          MW_densbinborderval, gaps_border, gaps_center = resample_1dim_weightedmean(odat_ji, ogrd_ji, ngrd, fill_value=0)
          # b) absolute influx (offset) from high-density.
          idxn_last = np.where(ngrd <= np.nanmax(ogrd_ji))[0][-1] # last idxn which is smaller than ogrd.max()
          idco_highdens = np.where(ogrd_ji > ngrd[idxn_last])[0]  # array with all idxo where ogrd is greater than last ngrd
          # b-1) variant where 1st idxo is taken
          #influx_highdens[j,i] = odat_ji[idco_highdens[0]]        # first idxo after idxn_last                 
          # b-2) variant where 2nd idxo is taken
          #      if this does not exist, influx_highdens is left on zero (preallocation)
          if len(idco_highdens) > 1:
            influx_highdens[j,i] = odat_ji[idco_highdens[1]]      # second idxo after idxn_last (the first is already used for interpolation of last ndat.)
          # c) get differences in MW_mgrd by substracting outflux from influx at each bin
          MWdiff_densbin = -1*np.diff(MW_densbinborderval)        # factor *-1 as MW is upward and diff goes downward
          gaps_border = gaps_border[:-1] | gaps_border[1:]        # gaps need to changed too ('|' requires no gap on both sides: the very careful way)
          gaps_center = gaps_center[:-1] | gaps_center[1:]            
          # d) cumulative integration from dense to light water starting with influx_highdens as an offset.
          # --> before and aferwards fill gaps with 0 and fill_value, respectively.
          
          MWdiff_densbin = fill_gaps(MWdiff_densbin, gaps_border, gaps_center, fill_value=0)  # mask gaps with 0
          ndat[:,j,i] = influx_highdens[j,i] + np.cumsum(MWdiff_densbin[::-1])[::-1]          # cumulative sum
          ndat[:,j,i] = fill_gaps(ndat[:,j,i], gaps_border, gaps_center, fill_value)          # mask gaps with fill_value

        elif method == 'sum': #! doesn't work right now!
          ndat[:,j,i], gaps_border, gaps_center = resample_1dim_sum(odat_ji, ogrd_ji, ngrd, fill_value)            
          ndat[:,j,i] = fill_gaps(ndat[:,j,i], gaps_border, gaps_center, fill_value)
        
    utils_misc.ProgBar('done')
    
    # some statistics
    if method == 'dMW':
        print('Statistics on influx_highdens:' \
              '\n mean:   {}\n median: {}\n min:    {}\n max:    {}'.format(\
              np.nanmean(influx_highdens), np.nanmedian(influx_highdens), \
              np.nanmin(influx_highdens), np.nanmax(influx_highdens)))

    return(np.squeeze(ndat)) # remove singleton dimensions (i.e. (1d --> 3d) --> back to 1d)


# =======================================================================================
# - resampling data on new grid using linear interpolation (along single dimension)
# =======================================================================================
def resample_colwise_on_zgrd(odat, zogrd, zngrd, method, fill_value=np.nan, mask='none', sort_zngrd='False'):
    '''
    Differences to resample_colwise():
     > zogrd is same everywhere, whereas zngrd changes btw. columns
    Uses:
     > resample_1dim_weightedmean()
    Input:
     > odat:        data on model grid
     > zogrd:        old grid
     > zngrd:        new grid 
     > method:      string | 'wmean' (weighted mean), 'dMW' or 'sum' (sum over all datapoints within bin on new grid)
     > fill_value:  float or nan | value used to fill in for requested points outside the range of zogrd.
     > mask:        mask of densityshape [j,i], default: all True (no mask)
     > sort_zngrd:   bool |  if True: zogrd (and odat) will be sorted such that zogrd is monotonically increasing (not necess. in strict sense!)
                            if False: zogrd will be manipulated such that strictly monotonically increasing (brute method, not recommended)
    Output:
     > ndat:    resampled data on new grid
    Comments:
     > add warning if negative gradients occur.
     > #!! for discreapencies at the low- and high-density borders see the comment in resample_1dim_weightedmean().
    '''

    def fill_gaps(ndat_ji, gaps_border, gaps_center, fill_value):
        ndat_ji[gaps_border] = fill_value
        #ndat_ji[gaps_center] = fill_value
        return(ndat_ji)
        
    print(' > columnwise resampling')

    # shape of data-array
    if len(odat.shape)==3:
      len_j = odat.shape[-2] # assumed to be the second last entry
      len_i = odat.shape[-1] # assumed to be the last entry
    elif len(odat.shape)==2:
      len_j = 1              # set to one, such that j-loop is executed only once.
      len_i = odat.shape[-1] # assumed to be the last entry
    elif len(odat.shape)==1:
      len_j = 1              # set to one, such that j-loop is executed only once.
      len_i = 1              # set to one, such that i-loop is executed only once.

    # get default for regional mask (i.e. do not mask anything)
    if mask == 'none':
      mask = np.ones(shape=[len_j, len_i], dtype=bool)

    # expand the shape of zngrd and odat to 3 dimensions (singleton-dimensions are intended)
    ndim_zngrd = len(zngrd.shape)
    ndim_odat = len(odat.shape)
    if ndim_zngrd==1: # 1dim --> 3dim
      zngrd = utils_conv.expand_karray_to_kji(zngrd, len_j, len_i)
    elif ndim_zngrd==2: # 2dim --> 3dim
      sys.exit('case of two-dimensional zngrd is not implemented yet!')
      
    if ndim_odat==1: # 1dim --> 3dim
      odat = utils_conv.expand_karray_to_kji(odat, len_j, len_i)
    elif ndim_odat==2: # 2dim --> 3dim
      sys.exit('case of two-dimensional odat is not implemented yet!')
      
    # pre-allocation of ndat
    if method == 'wmean':
        ndat = fill_value * np.ones(shape=[zngrd.shape[0], len_j, len_i])
    elif method == 'dMW':
        ndat = fill_value * np.ones(shape=[zngrd.shape[0]-1, len_j, len_i])
        influx_highdens = np.zeros(shape=[len_j, len_i])
    
    # loop over columns
    for j in np.arange(len_j):
      utils_misc.ProgBar('step', step=j, nsteps=len_j)
      for i in np.arange(len_i):
        if mask[j,i]==True: # skip masked [j,i]-tuples
          # reduce zngrd and odat to current column
          zngrd_ji = zngrd[:,j,i]
          odat_ji = odat[:,j,i]
          # detect disjunct zogrd and zngrd_ji and continue with next column
          if (np.nanmax(zogrd) < np.nanmin(zngrd_ji)) or (np.nanmax(zngrd_ji) < np.nanmin(zogrd)):
              print('disjunct zogrd and zngrd at (j,i)=({}, {}). (please check conservation of integrated flux!)'.format(j, i))
              continue
          # make zngrd_ji strictly monotoneously increasing
          if any(np.diff(zngrd_ji)<=0):
            if sort_zngrd_ji == 'True':
              # sort zngrd_ji and odat (recommended)
              idx_sort = np.argsort(zngrd_ji)
              zngrd_ji = zngrd_ji[idx_sort]
              odat_ji = odat_ji[idx_sort]
            else:
              # brute method by manipulating zngrd_ji (not recommended but used by Kwon)
              for k in np.arange(1,zngrd_ji.shape[0]):
                if zngrd_ji[k] <= zngrd_ji[k-1]:
                  zngrd_ji[k] = zngrd_ji[k-1]+1e-10
          
          # interpolation
          if method == 'wmean': # simple weighted mean interpolation
            ndat[:,j,i], gaps_border, gaps_center = resample_1dim_weightedmean(odat_ji, zogrd, zngrd_ji, fill_value)
            ndat[:,j,i] = fill_gaps(ndat[:,j,i], gaps_border, gaps_center, fill_value)
            
          elif method == 'dMW': # procedure for dMW
            # a) weighted mean of closest neighbours around dens_bin border values
            MW_densbinborderval, gaps_border, gaps_center = resample_1dim_weightedmean(odat_ji, zogrd, zngrd_ji, fill_value=0)
            # b) absolute influx (offset) from high-density.
            try: 
                idxn_last = np.where(zngrd_ji <= np.nanmax(zogrd))[0][-1] # last idxn which is smaller than zogrd.max()
                idco_highdens = np.where(zogrd > zngrd_ji[idxn_last])[0]  # array with all idxo where zogrd is greater than last zngrd_ji
                # b-1) variant where 1st idxo is taken
                #influx_highdens[j,i] = odat_ji[idco_highdens[0]]        # first idxo after idxn_last                 
                # b-2) variant where 2nd idxo is taken
                #      if this does not exist, influx_highdens is left on zero (preallocation)
                if len(idco_highdens) > 1:
                  influx_highdens[j,i] = odat_ji[idco_highdens[1]]      # second idxo after idxn_last (the first is already used for interpolation of last ndat.)
            except: 
                pass #! workaround as it happened that all values in zngrd_ji were nans --> idxn_last could not be found.
            # c) get differences in MW_mgrd by substracting outflux from influx at each bin
            MWdiff_densbin = -1*np.diff(MW_densbinborderval)        # factor *-1 as MW is upward and diff goes downward
            gaps_border = gaps_border[:-1] | gaps_border[1:]        # gaps need to changed too ('|' requires no gap on both sides: the very careful way)
            gaps_center = gaps_center[:-1] | gaps_center[1:]            
            # d) cumulative integration from dense to light water starting with influx_highdens as an offset.
            # --> before and aferwards fill gaps with 0 and fill_value, respectively.
            
            MWdiff_densbin = fill_gaps(MWdiff_densbin, gaps_border, gaps_center, fill_value=0)  # mask gaps with 0
            ndat[:,j,i] = influx_highdens[j,i] + np.cumsum(MWdiff_densbin[::-1])[::-1]          # cumulative sum
            ndat[:,j,i] = fill_gaps(ndat[:,j,i], gaps_border, gaps_center, fill_value)          # mask gaps with fill_value

          elif method == 'sum': #! doesn't work right now!
            ndat[:,j,i], gaps_border, gaps_center = resample_1dim_sum(odat_ji, zogrd, zngrd_ji, fill_value)            
            ndat[:,j,i] = fill_gaps(ndat[:,j,i], gaps_border, gaps_center, fill_value)
          
    utils_misc.ProgBar('done')
    
    # some statistics
    if method == 'dMW':
        print('Statistics on influx_highdens:' \
              '\n mean:   {}\n median: {}\n min:    {}\n max:    {}'.format(\
              np.nanmean(influx_highdens), np.nanmedian(influx_highdens), \
              np.nanmin(influx_highdens), np.nanmax(influx_highdens)))
        print('foo')

    return(np.squeeze(ndat)) # remove singleton dimensions (i.e. (1d --> 3d) --> back to 1d)


def resample_1dim_lininterp(odat, ogrd, ngrd, idxn_start = fill_value=np.nan):
    import CESM_utils_analysis as utils_ana
    '''
    Input:
     > odat:            data on old grid
     > ogrd:            old grid
     > ngrd:            new grid (resampling grid)
     > fill_value:      0 or np.nan | used for preallocation
    Output:
     > ndat:            data on new grid (resampled data)
     > gaps_border      bool-array, same shape as ndat | True where ngrd is outside ogrd
     > gaps_center      bool-array, same shape as ndat | True where resolution of ngrd is higher than resolution of ogrd
    Comments:
     > #!! implemented, such that both, the low- and the high-value level of ngrd are left blank (i.e. on fill_value)
       --> for dMW a shift of the cumsumvalues will correct for the high-density border discrepency 
       and a comparison with the columnwise integrated transport in depth-space for those at the low-density border.
    '''

    # Pre-allocation of ndat | if ngrd is longer than ogrd fill tail with nans
    ndat = fill_value * np.ones(shape=ngrd.shape)
    gaps_border = np.zeros(shape=ngrd.shape, dtype=bool)
    gaps_center = np.zeros(shape=ngrd.shape, dtype=bool)
    # Resampling
    idxo = 0                                        # index on ogrd
    idxo_old = np.nan
    try: 
        idxn = np.where(ngrd > np.nanmin(ogrd))[0][0]   # index on ngrd
        gaps_border[:idxn] = True                       # mark gaps at low-value border 
    except:
        gaps_border[:] = True                           # mark all points as bordergaps
        return(ndat, gaps_border, gaps_center)
    # loop through ngrd: 
        # cond 1) loop until the very last enty of ngrd.
        # cond 2) stop as soon as remaining ngrd values are all higher than the maximum value of ogrd.
    while (idxn < len(ngrd)) and (ngrd[idxn] <= np.nanmax(ogrd)):
      # lift idxo until ogrd is one step further (deeper, i.e. higher value) than ngrd.
      while (idxo < len(ogrd)-1) and (ngrd[idxn] > ogrd[idxo]):
        idxo += 1
      # resampling
      if idxo == idxo_old:                      # ngrd-resolution is higher than ogrd-resolution
        gaps_center[idxn-1:idxn+1] = True       # --> mark idxn and idxn-1 as gaps (the very careful way! #! consider defining minimum span to also mark idxn-1)
        ndat[idxn] = ndat[idxn-1]               # --> same value as above. These value remains if center-gaps are not masked.
      else:                                     # centre values
        diff_1 = -1*(ogrd[idxo-1] - ngrd[idxn]) # --> a positive number
        diff_2 = ogrd[idxo] - ngrd[idxn]        # --> a positive number
        diff_total = diff_1 + diff_2            # --> equivalent to ogrd[idxo] - ogrd[idxo-1]
        # linearly weighted interpolation
        ndat[idxn] = odat[idxo-1]*diff_2/diff_total + odat[idxo]*diff_1/diff_total
      # increase idxn and renew idxo_old
      idxo_old = idxo                         
      idxn += 1

    gaps_border[idxn:] = True                   # mark gaps at high-value border
    return(ndat, gaps_border, gaps_center)
    
    
def resample_1dim_sum(odat, ogrd, ngrd, fill_value=np.nan):
    '''
    Input:
     > odat: data on old grid (model grid)
     > ogrd:       old grid (model grid)
     > ngrd:       new grid (resampling grid)
     > fill_value:      0 or np.nan | used for preallocation
    Output:
     > ndat:       data on new grid (resampled data)
    Comments:
     > not sure whether sum=0 should be overwritten with nans
    '''
    
    # Pre-allocation of ndat | if ngrd is longer than ogrd fill tail with nans
    ndat = fill_value * np.ones(shape = len(ngrd))
    
    # Resampling
    inds = np.digitize(ogrd, ngrd)
    for i in np.arange(1,len(ngrd)):
        pass #!



# =======================================================================================
#  Color Maps
# =======================================================================================

# get_cmap interpolates the colormap on a discrete range.
def get_cmap(min, max, nlevels, scheme):
    cmlist=[];
    levels=np.linspace(min,max,nlevels)
    cmlist=np.array([int(x) for x in np.linspace(0,252,len(levels)+1)]) # interpolation on [0,252]
    cmap, norm = from_levels_and_colors(levels, scheme(cmlist) , extend='both')
    return cmap, norm

   
    
# The viridis colormap will be the default in matplotlib 2.0
def get_viridis():
    return(ml.colors.ListedColormap([[0.267004, 0.004874, 0.329415],
                 [0.268510, 0.009605, 0.335427],
                 [0.269944, 0.014625, 0.341379],
                 [0.271305, 0.019942, 0.347269],
                 [0.272594, 0.025563, 0.353093],
                 [0.273809, 0.031497, 0.358853],
                 [0.274952, 0.037752, 0.364543],
                 [0.276022, 0.044167, 0.370164],
                 [0.277018, 0.050344, 0.375715],
                 [0.277941, 0.056324, 0.381191],
                 [0.278791, 0.062145, 0.386592],
                 [0.279566, 0.067836, 0.391917],
                 [0.280267, 0.073417, 0.397163],
                 [0.280894, 0.078907, 0.402329],
                 [0.281446, 0.084320, 0.407414],
                 [0.281924, 0.089666, 0.412415],
                 [0.282327, 0.094955, 0.417331],
                 [0.282656, 0.100196, 0.422160],
                 [0.282910, 0.105393, 0.426902],
                 [0.283091, 0.110553, 0.431554],
                 [0.283197, 0.115680, 0.436115],
                 [0.283229, 0.120777, 0.440584],
                 [0.283187, 0.125848, 0.444960],
                 [0.283072, 0.130895, 0.449241],
                 [0.282884, 0.135920, 0.453427],
                 [0.282623, 0.140926, 0.457517],
                 [0.282290, 0.145912, 0.461510],
                 [0.281887, 0.150881, 0.465405],
                 [0.281412, 0.155834, 0.469201],
                 [0.280868, 0.160771, 0.472899],
                 [0.280255, 0.165693, 0.476498],
                 [0.279574, 0.170599, 0.479997],
                 [0.278826, 0.175490, 0.483397],
                 [0.278012, 0.180367, 0.486697],
                 [0.277134, 0.185228, 0.489898],
                 [0.276194, 0.190074, 0.493001],
                 [0.275191, 0.194905, 0.496005],
                 [0.274128, 0.199721, 0.498911],
                 [0.273006, 0.204520, 0.501721],
                 [0.271828, 0.209303, 0.504434],
                 [0.270595, 0.214069, 0.507052],
                 [0.269308, 0.218818, 0.509577],
                 [0.267968, 0.223549, 0.512008],
                 [0.266580, 0.228262, 0.514349],
                 [0.265145, 0.232956, 0.516599],
                 [0.263663, 0.237631, 0.518762],
                 [0.262138, 0.242286, 0.520837],
                 [0.260571, 0.246922, 0.522828],
                 [0.258965, 0.251537, 0.524736],
                 [0.257322, 0.256130, 0.526563],
                 [0.255645, 0.260703, 0.528312],
                 [0.253935, 0.265254, 0.529983],
                 [0.252194, 0.269783, 0.531579],
                 [0.250425, 0.274290, 0.533103],
                 [0.248629, 0.278775, 0.534556],
                 [0.246811, 0.283237, 0.535941],
                 [0.244972, 0.287675, 0.537260],
                 [0.243113, 0.292092, 0.538516],
                 [0.241237, 0.296485, 0.539709],
                 [0.239346, 0.300855, 0.540844],
                 [0.237441, 0.305202, 0.541921],
                 [0.235526, 0.309527, 0.542944],
                 [0.233603, 0.313828, 0.543914],
                 [0.231674, 0.318106, 0.544834],
                 [0.229739, 0.322361, 0.545706],
                 [0.227802, 0.326594, 0.546532],
                 [0.225863, 0.330805, 0.547314],
                 [0.223925, 0.334994, 0.548053],
                 [0.221989, 0.339161, 0.548752],
                 [0.220057, 0.343307, 0.549413],
                 [0.218130, 0.347432, 0.550038],
                 [0.216210, 0.351535, 0.550627],
                 [0.214298, 0.355619, 0.551184],
                 [0.212395, 0.359683, 0.551710],
                 [0.210503, 0.363727, 0.552206],
                 [0.208623, 0.367752, 0.552675],
                 [0.206756, 0.371758, 0.553117],
                 [0.204903, 0.375746, 0.553533],
                 [0.203063, 0.379716, 0.553925],
                 [0.201239, 0.383670, 0.554294],
                 [0.199430, 0.387607, 0.554642],
                 [0.197636, 0.391528, 0.554969],
                 [0.195860, 0.395433, 0.555276],
                 [0.194100, 0.399323, 0.555565],
                 [0.192357, 0.403199, 0.555836],
                 [0.190631, 0.407061, 0.556089],
                 [0.188923, 0.410910, 0.556326],
                 [0.187231, 0.414746, 0.556547],
                 [0.185556, 0.418570, 0.556753],
                 [0.183898, 0.422383, 0.556944],
                 [0.182256, 0.426184, 0.557120],
                 [0.180629, 0.429975, 0.557282],
                 [0.179019, 0.433756, 0.557430],
                 [0.177423, 0.437527, 0.557565],
                 [0.175841, 0.441290, 0.557685],
                 [0.174274, 0.445044, 0.557792],
                 [0.172719, 0.448791, 0.557885],
                 [0.171176, 0.452530, 0.557965],
                 [0.169646, 0.456262, 0.558030],
                 [0.168126, 0.459988, 0.558082],
                 [0.166617, 0.463708, 0.558119],
                 [0.165117, 0.467423, 0.558141],
                 [0.163625, 0.471133, 0.558148],
                 [0.162142, 0.474838, 0.558140],
                 [0.160665, 0.478540, 0.558115],
                 [0.159194, 0.482237, 0.558073],
                 [0.157729, 0.485932, 0.558013],
                 [0.156270, 0.489624, 0.557936],
                 [0.154815, 0.493313, 0.557840],
                 [0.153364, 0.497000, 0.557724],
                 [0.151918, 0.500685, 0.557587],
                 [0.150476, 0.504369, 0.557430],
                 [0.149039, 0.508051, 0.557250],
                 [0.147607, 0.511733, 0.557049],
                 [0.146180, 0.515413, 0.556823],
                 [0.144759, 0.519093, 0.556572],
                 [0.143343, 0.522773, 0.556295],
                 [0.141935, 0.526453, 0.555991],
                 [0.140536, 0.530132, 0.555659],
                 [0.139147, 0.533812, 0.555298],
                 [0.137770, 0.537492, 0.554906],
                 [0.136408, 0.541173, 0.554483],
                 [0.135066, 0.544853, 0.554029],
                 [0.133743, 0.548535, 0.553541],
                 [0.132444, 0.552216, 0.553018],
                 [0.131172, 0.555899, 0.552459],
                 [0.129933, 0.559582, 0.551864],
                 [0.128729, 0.563265, 0.551229],
                 [0.127568, 0.566949, 0.550556],
                 [0.126453, 0.570633, 0.549841],
                 [0.125394, 0.574318, 0.549086],
                 [0.124395, 0.578002, 0.548287],
                 [0.123463, 0.581687, 0.547445],
                 [0.122606, 0.585371, 0.546557],
                 [0.121831, 0.589055, 0.545623],
                 [0.121148, 0.592739, 0.544641],
                 [0.120565, 0.596422, 0.543611],
                 [0.120092, 0.600104, 0.542530],
                 [0.119738, 0.603785, 0.541400],
                 [0.119512, 0.607464, 0.540218],
                 [0.119423, 0.611141, 0.538982],
                 [0.119483, 0.614817, 0.537692],
                 [0.119699, 0.618490, 0.536347],
                 [0.120081, 0.622161, 0.534946],
                 [0.120638, 0.625828, 0.533488],
                 [0.121380, 0.629492, 0.531973],
                 [0.122312, 0.633153, 0.530398],
                 [0.123444, 0.636809, 0.528763],
                 [0.124780, 0.640461, 0.527068],
                 [0.126326, 0.644107, 0.525311],
                 [0.128087, 0.647749, 0.523491],
                 [0.130067, 0.651384, 0.521608],
                 [0.132268, 0.655014, 0.519661],
                 [0.134692, 0.658636, 0.517649],
                 [0.137339, 0.662252, 0.515571],
                 [0.140210, 0.665859, 0.513427],
                 [0.143303, 0.669459, 0.511215],
                 [0.146616, 0.673050, 0.508936],
                 [0.150148, 0.676631, 0.506589],
                 [0.153894, 0.680203, 0.504172],
                 [0.157851, 0.683765, 0.501686],
                 [0.162016, 0.687316, 0.499129],
                 [0.166383, 0.690856, 0.496502],
                 [0.170948, 0.694384, 0.493803],
                 [0.175707, 0.697900, 0.491033],
                 [0.180653, 0.701402, 0.488189],
                 [0.185783, 0.704891, 0.485273],
                 [0.191090, 0.708366, 0.482284],
                 [0.196571, 0.711827, 0.479221],
                 [0.202219, 0.715272, 0.476084],
                 [0.208030, 0.718701, 0.472873],
                 [0.214000, 0.722114, 0.469588],
                 [0.220124, 0.725509, 0.466226],
                 [0.226397, 0.728888, 0.462789],
                 [0.232815, 0.732247, 0.459277],
                 [0.239374, 0.735588, 0.455688],
                 [0.246070, 0.738910, 0.452024],
                 [0.252899, 0.742211, 0.448284],
                 [0.259857, 0.745492, 0.444467],
                 [0.266941, 0.748751, 0.440573],
                 [0.274149, 0.751988, 0.436601],
                 [0.281477, 0.755203, 0.432552],
                 [0.288921, 0.758394, 0.428426],
                 [0.296479, 0.761561, 0.424223],
                 [0.304148, 0.764704, 0.419943],
                 [0.311925, 0.767822, 0.415586],
                 [0.319809, 0.770914, 0.411152],
                 [0.327796, 0.773980, 0.406640],
                 [0.335885, 0.777018, 0.402049],
                 [0.344074, 0.780029, 0.397381],
                 [0.352360, 0.783011, 0.392636],
                 [0.360741, 0.785964, 0.387814],
                 [0.369214, 0.788888, 0.382914],
                 [0.377779, 0.791781, 0.377939],
                 [0.386433, 0.794644, 0.372886],
                 [0.395174, 0.797475, 0.367757],
                 [0.404001, 0.800275, 0.362552],
                 [0.412913, 0.803041, 0.357269],
                 [0.421908, 0.805774, 0.351910],
                 [0.430983, 0.808473, 0.346476],
                 [0.440137, 0.811138, 0.340967],
                 [0.449368, 0.813768, 0.335384],
                 [0.458674, 0.816363, 0.329727],
                 [0.468053, 0.818921, 0.323998],
                 [0.477504, 0.821444, 0.318195],
                 [0.487026, 0.823929, 0.312321],
                 [0.496615, 0.826376, 0.306377],
                 [0.506271, 0.828786, 0.300362],
                 [0.515992, 0.831158, 0.294279],
                 [0.525776, 0.833491, 0.288127],
                 [0.535621, 0.835785, 0.281908],
                 [0.545524, 0.838039, 0.275626],
                 [0.555484, 0.840254, 0.269281],
                 [0.565498, 0.842430, 0.262877],
                 [0.575563, 0.844566, 0.256415],
                 [0.585678, 0.846661, 0.249897],
                 [0.595839, 0.848717, 0.243329],
                 [0.606045, 0.850733, 0.236712],
                 [0.616293, 0.852709, 0.230052],
                 [0.626579, 0.854645, 0.223353],
                 [0.636902, 0.856542, 0.216620],
                 [0.647257, 0.858400, 0.209861],
                 [0.657642, 0.860219, 0.203082],
                 [0.668054, 0.861999, 0.196293],
                 [0.678489, 0.863742, 0.189503],
                 [0.688944, 0.865448, 0.182725],
                 [0.699415, 0.867117, 0.175971],
                 [0.709898, 0.868751, 0.169257],
                 [0.720391, 0.870350, 0.162603],
                 [0.730889, 0.871916, 0.156029],
                 [0.741388, 0.873449, 0.149561],
                 [0.751884, 0.874951, 0.143228],
                 [0.762373, 0.876424, 0.137064],
                 [0.772852, 0.877868, 0.131109],
                 [0.783315, 0.879285, 0.125405],
                 [0.793760, 0.880678, 0.120005],
                 [0.804182, 0.882046, 0.114965],
                 [0.814576, 0.883393, 0.110347],
                 [0.824940, 0.884720, 0.106217],
                 [0.835270, 0.886029, 0.102646],
                 [0.845561, 0.887322, 0.099702],
                 [0.855810, 0.888601, 0.097452],
                 [0.866013, 0.889868, 0.095953],
                 [0.876168, 0.891125, 0.095250],
                 [0.886271, 0.892374, 0.095374],
                 [0.896320, 0.893616, 0.096335],
                 [0.906311, 0.894855, 0.098125],
                 [0.916242, 0.896091, 0.100717],
                 [0.926106, 0.897330, 0.104071],
                 [0.935904, 0.898570, 0.108131],
                 [0.945636, 0.899815, 0.112838],
                 [0.955300, 0.901065, 0.118128],
                 [0.964894, 0.902323, 0.123941],
                 [0.974417, 0.903590, 0.130215],
                 [0.983868, 0.904867, 0.136897],
                 [0.993248, 0.906157, 0.143936]]))


