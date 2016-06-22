# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 15:22:04 2016

@author: buerki@climate.unibe.ch
"""

import sys

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
    else: sys.exit('Invalid name for model run')
    
# ---------------------------------------------------------------------------------------
# Directories for Variables named by auxillary grid
def get_path2var(auxgrd_name):
    # lat: 170 equally spaced boxes from 80S to 90N | z: 60 boxes
    if auxgrd_name == 'lat170eq80S90N_zeq60': 
        return('vars_lat170eq80S90N_zeq60/')
    # lat: 340 equally spaced boxes from 80S to 90N | z: 60 boxes
    elif auxgrd_name == 'lat340eq80S90N_zeq60': 
        return('vars_lat340eq80S90N_zeq60/')
    # lat: as in ncdat.lat_aux_grid | z: 60 boxes
    elif auxgrd_name == 'lat395model_zeq60':
        return('vars_lat395model_zeq60/')
    # lat: as in ncdat.lat_aux_grid but only every other entry | z: 60 boxes
    elif auxgrd_name == 'lat198model_zeq60':
        return('vars_lat198model_zeq60/')
