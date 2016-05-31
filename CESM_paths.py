# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 15:22:04 2016

@author: buerki@climate.unibe.ch
"""

global path_lm_1deg, path_ctr_1deg, path_ctr_2deg


# last millenium run
''' ** some more comments ** '''

lm_1deg_root = '/alphadata02/born/lm850-1850.1deg/'
path_lm_1deg = dict(
  anndat      = lm_1deg_root + 'no_backup/annual_data/', 
  ppanndat    = lm_1deg_root + 'no_backup/postproc_annual_data/',
  mondat      = lm_1deg_root + 'montly_data/',
  sens        = lm_1deg_root + 'sensitivity_runs/')



# control runs
''' ** some more comments ** '''
# 1-deg atmospherical resolution
ctr_1deg_root = '/alphadata02/born/no_backup/b40.1850.track1.1deg.006/'
path_ctr_1deg = dict(
  anndat      = ctr_1deg_root + 'annual_data/',                                
  anndattmp   = ctr_1deg_root + 'annual_data/tmp/',                              #?? what's in here?? only RHO and SALT but why tmp??
  pp          = ctr_1deg_root + 'postprocessed/',                                       
  ppannav     = ctr_1deg_root + 'postprocessed/annual_average/')                 #?? what's the difference to ctr_1deg_anndat ?

# 2-deg atmospherical resolution
ctr_2deg_root = '/alphadata02/born/no_backup/b40.1850.track1.2deg.003/'
path_ctr_2deg = dict(
  anndat     = ctr_2deg_root + '/annual_data/',
  nobckp     = ctr_2deg_root + '/no_backup/',                                    #?? what's in here? original runs already sep. for vars?
  nobckphist = ctr_2deg_root + '/no_backup/hist',                                #?? what's in here? orig. run not even sep. for vars?
  anndat0    = ctr_2deg_root + '/restart_yr999/')                                #?? restart important?
 

# ---------------------------------------------------------------------------------------
# Directories for Variables named by auxillary grid
def get_path2var(auxgrd_name):
    # lat: 170 equally spaced boxes from 80S to 90N | z: 60 boxes
    if auxgrd_name == 'lateq80S90N_zeq60': 
        return('vars_lateq80S90N_zeq60/')
    # lat: as in ncdat.lat_aux_grid but only every other entry | z: 60 boxes
    elif auxgrd_name == 'latMOCmodelEveryOther_zeq60':
        return('vars_latMOCmodelEveryOther_zeq60/')
