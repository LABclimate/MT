# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 15:22:04 2016

@author: buerki@climate.unibe.ch
foo
"""

import sys
import UTILS_misc as utils_misc

# ---------------------------------------------------------------------------------------
# Directories of netCF data
def get_path2data(runname, key):
    
    # pre-industrial last millenium run
    ''' ** some more comments ** '''
    
    lm_1deg_root = '/alphadata02/born/lm850-1850.1deg/'
    path_lm_1deg = dict(
      anndat      = lm_1deg_root + 'no_backup/annual_data/', 
      ppannav     = lm_1deg_root + 'no_backup/postproc_annual_data/',
      mondat      = lm_1deg_root + 'monthly_data/',
      sens        = lm_1deg_root + 'sensitivity_runs/')
    
    # control runs
    ''' ** some more comments ** '''
    # 1-deg atmospherical resolution 1300 y
    ctr_1deg_root = '/alphadata02/born/no_backup/b40.1850.track1.1deg.006/'
    path_ctr_1deg = dict(
      anndat      = ctr_1deg_root + 'annual_data/',                                
      anndattmp   = ctr_1deg_root + 'annual_data/tmp/',                              #?? what's in here?? only RHO and SALT but why tmp??
      pp          = ctr_1deg_root + 'postprocessed/',                                       
      ppannav     = ctr_1deg_root + 'postprocessed/annual_average/')                 #?? what's the difference to ctr_1deg_anndat ?
    
    # 2-deg atmospherical resolution 999 y
    ctr_2deg_root = '/alphadata02/born/no_backup/b40.1850.track1.2deg.003/'
    path_ctr_2deg = dict(
      anndat     = ctr_2deg_root + '/annual_data/',
      nobckp     = ctr_2deg_root + '/no_backup/',                                    #?? what's in here? original runs already sep. for vars?
      nobckphist = ctr_2deg_root + '/no_backup/hist',                                #?? what's in here? orig. run not even sep. for vars?
      anndat0    = ctr_2deg_root + '/restart_yr999/')                                #?? restart important?
     
    if runname == 'lm_1deg':
        try:    return(path_lm_1deg[key])
        except: sys.exit('Invalid key: ' + key + 'Valid keys are: ' + str(path_lm_1deg.keys))
    elif runname == 'ctr_1deg':
        try:    return(path_ctr_1deg[key])
        except: sys.exit('Invalid key: ' + key + 'Valid keys are: ' + str(path_ctr_1deg.keys))
    elif runname == 'ctr_2deg':
        try:    return(path_ctr_2deg[key])
        except: sys.exit('Invalid key: ' + key + 'Valid keys are: ' + str(path_ctr_2deg.keys))
    # Raise error for invalid key
    else: 
        sys.exit('Invalid name for model run')
        return()
    
# ---------------------------------------------------------------------------------------
# Directories for Variables to be stored
def get_path2vars(vartype, CESMversion, mkdir=False):
    
    vars_root = '../variables/'+'CESM_V'+str(CESMversion)+'/'

    # variables fixed to grid like maxdepths etc.
    if vartype == 'mgrd':
        dirname = vars_root+'vars_mgrd/'
    # variables resampled on densityaxis
    elif vartype == 'dens':
        dirname = vars_root+'vars_dens/'
    # correlation indices etc...
    elif vartype == 'corr':
        dirname = vars_root+'vars_corr/'
        
    # lat: 170 equally spaced boxes from 80S to 90N | z: 60 boxes
    elif vartype == 'lat170eq80S90N': 
        dirname = vars_root+'vars_auxgrd/lat170eq80S90N/'
    # lat: 340 equally spaced boxes from 80S to 90N | z: 60 boxes
    elif vartype == 'lat340eq80S90N': 
        dirname = vars_root+'vars_auxgrd/lat340eq80S90N/'
    # lat: as in ncdat.lat_aux_grid | z: 60 boxes
    elif vartype == 'lat395model':
        dirname = vars_root+'vars_auxgrd/lat395model/'
    # lat: as in ncdat.lat_aux_grid but only every other entry | z: 60 boxes
    elif vartype == 'lat198model':
        dirname = vars_root+'vars_auxgrd/lat198model/'
    # lat: as in ncdat.lat_aux_grid but only every other entry | z: 60 boxes
    elif vartype == 'lat99model':
        dirname = vars_root+'vars_auxgrd/lat99model/'        
    # Raise error for invalid key
    else:
        sys.exit('Invalid key for path2vars')
        return()

    # -----------------------------------
    # return dirname and create directory
    # -----------------------------------
    if mkdir == True:
        utils_misc.mkdir(dirname)
    return(dirname)
    
# ---------------------------------------------------------------------------------------
# Directories for Variables to be stored
def get_path2figs(vartype):
    
    figs_root = '../figures/'

    # dump miscellenian test figures here
    if vartype == 'misc':
        dirname = vars_root+'figs_misc/'
    # correlation indices etc...
    elif vartype == 'corr':
        dirname = vars_root+'figs_corr/'
    # Raise error for invalid key
    else: 
        sys.exit('Invalid key for path2vars')
        return()
    
    # -----------------------------------
    # return dirname and create directory
    # -----------------------------------
    if mkdir == True:
        utils_misc.mkdir(dirname)
    return(dirname)
    