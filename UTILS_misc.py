#################################
# ----------- Toolbox -----------
# ----- miscellaneous tools -----
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import UTILS_misc as utils_misc
#################################
# contained functions:
#################################
# - ProgBar()
# - loadvar()
# - savevar()
# - checkdir()
# - primefactors()
# - LGS()
# - getnsave()
#################################
# please log your changes below
#################################
# 14-Mai-2016 - buerki@climate.unibe.ch: created this toolbox
#                                        created ProgBar()
# 20?-Mai-2016 -buerki@climate.unibe.ch: created loadvar()
#                                        created savevar()
# 25-Mai-2016 - buerki@climate.unibe.ch: added forceinit-patch to ProgBar() -> needs some improvement
#                                        changed toolbox name from utils_spec to utils_misc
# 31-Mai-2016 - buerki@climate.unibe.ch: inProgBar added auto-calculation of barlength
#                                        added checkdir()
# 23-Jun-2016 - buerki@climate.unibe.ch: changed the name of checkdir() to mkdir()
# 19-Aug-2016 - buerki@climate.unibe.ch: created LGS()
# 05-Oct-2016 - buerki@climate.unibe.ch: created LG() and GS()
#################################

import numpy as np
import sys
import xarray as xr
import pickle
import os
import UTILS_misc as utils_misc

# --- dump 3dim matrix
def getmat():
    mat = np.array([[[1,1,1],[2,2,2],[3,3,3],[4,4,4]],[[5,5,5],[6,6,6],[7,7,7],[8,8,8]]])
    mat[:,:,1] *= 10
    mat[:,:,2] *= 100
    return mat
    
    
# --- Progress Bars ---
def ProgBar(stage, step=99, nsteps=99,  maxlen = 60, forceinit=False, minbarlen=20):
    '''
    Comments:
      > improve forceinit behaviour, such that not limited to step==1.
    '''
    if stage=='step':
        # calculate barlen
        if nsteps <= maxlen:     
            barlen = nsteps
        else:
            for i in np.arange(2,maxlen+1)[::-1]:
                if nsteps%i==0:
                    barlen = i
                    break
        if barlen < minbarlen:
            barlen = minbarlen
        # draw points
        if (step%(np.round(nsteps/barlen))==0) | ((forceinit==True) & (step == 1)):
	    if (step==0) | ((forceinit==True) & (step == 1)): # initialization
	        print('  --> Progress: ['+' '*barlen+']'), 	#! do NOT delete the pending comma!!
	        print('\b'*(2+barlen)), # set the cursor back to start
	        sys.stdout.flush()
	    else:
                print('\b.'), 				        #! do NOT delete the pending comma!!
                sys.stdout.flush(),
    elif stage=='done':
        print('\b]  Done!')

def Counter(step, nsteps, stepname='Step', mod=1):
    if not np.mod(step, mod):
        print('\r{} {} / {}'.format(stepname, step, nsteps)),
        sys.stdout.flush()

    
'''
# code for testing
import time
for i in np.arange(21):
    time.sleep(.01)
    ProgBar('step', step = i, nsteps = 21)
utils_misc.ProgBar('done')
'''

# --- load variable using pickle
def loadvar(path_to_var, str_varname='__', verbose=True):
    if verbose: print(' > loading {}\n    from {}'.format(str_varname, path_to_var))
    sys.stdout.flush()    
    with open(path_to_var, 'rb') as f: 
      var = pickle.load(f)
    if verbose: print(' > success!')
    return(var)
    
# --- load variable using xarray
def loadxvar(path_to_var, str_varname='__', verbose=True):
    if verbose: print(' > loading {}\n    from {}.nc'.format(str_varname, path_to_var))
    if str_varname == None: xvar = xr.open_dataset(path_to_var+'.nc', decode_times=False)
    else:                   xvar = xr.open_dataset(path_to_var+'.nc', decode_times=False)[str_varname]
    if verbose: print(' > success!')
    return(xvar)
    
# --- save variable using pickle
def savevar(var, path_to_var, str_varname='__', verbose=True):
    if verbose: print(' > saving {} to file...'.format(str_varname))
    sys.stdout.flush()    
    with open(path_to_var, 'wb') as f:
      pickle.dump(var, f)
    if verbose: print(' > success!')

# --- save variable using xarray
def savexvar(xvar, path_to_var, str_varname='__', verbose=True):
    if verbose: print(' > saving {} to file...'.format(str_varname))
    var2save = xr.Dataset({str_varname:xvar}).to_netcdf(path_to_var+'.nc')
    if verbose: print(' > success!')

# --- try to load, else get and save variable
def LGS(fun_to_get_var, path_to_var, str_varname='__', verbose=True, noload=False, nosave=False, format='nc'):
    ''' 
    Usage:     > LGS(lambda: fun_to_get_var(args), path_to_var)
    '''
    try:
        if noload: raise ValueError()
        if verbose: print(' > trying to load {}\n    from {}...'.format(str_varname, path_to_var))
        if format=='nc':        var = utils_misc.loadxvar(path_to_var, str_varname, verbose=False)
        elif format=='pickle':  var = utils_misc.loadvar(path_to_var, verbose=False)
        if verbose: print(' > success!\n')
    except:
        if noload & verbose: print(' > no load (note that an eventually stored variable will be overwritten)\n > calculating {}...'.format(str_varname))
        elif verbose: print(' > failed to load!\n > calculating {}...'.format(str_varname))
        try: var = fun_to_get_var()
        except: var = fun_to_get_var # as "fun" might simply be any variable (which is uncallable and thus produces an error.)
        if verbose: print(' > success!\n > saving {} to file...'.format(str_varname))
        if format=='nc':        utils_misc.savexvar(var, path_to_var, str_varname, verbose=False)
        elif format=='pickle':  utils_misc.savevar(var, path_to_var, verbose=False)
        if verbose: print(' > success!\n')
    return(var)

# get and save variable (LGS with noload=True)
def GS(fun_to_get_var, path_to_var, str_varname, verbose=True):
    LGS(fun_to_get_var, path_to_var, str_varname, verbose, noload=True)

# try to load, else get variable (LGS with nosave=True)
def LG(fun_to_get_var, path_to_var, str_varname, verbose=True):
    LGS(fun_to_get_var, path_to_var, str_varname, verbose, nosave=True)

# --- add directory if unexistant

def mkdir(dirname):
    if os.path.isdir(dirname)==False:
      os.mkdir(dirname)
      print(' -> created new directory: ' + dirname)

#################################################################################################
# COLLECTION OF UNUSED OR UNTESTED FUNCTIONS 
#################################################################################################

# --- TextWrappers for indentions
import textwrap

user = "Username"
prefix = user + ":\t\t"
expanded_indent = textwrap.fill(prefix+'$', replace_whitespace=False)[:-1]
subsequent_indent = ' ' * len(expanded_indent)
wrapper = textwrap.TextWrapper(initial_indent=prefix,
                               subsequent_indent=subsequent_indent)
message = "LEFTLEFTLEFTLEFTLEFTLEFTLEFT RIGHTRIGHTRIGHT " * 3
#print wrapper.fill(message)


# --- find primefactors
''' Source:     http://datascientist.mabs.me/prime-factors-in-python/
    Comments:   NOT tested!!
'''
import math

def primefactors(n):
    primfac = []
    if (n % 2) == 0:
        lastFactor = 2
        n=n/2
        while (n % 2) == 0:
            n=n/2
    else:
        lastFactor =1
    factor = 3
    maxFactor = math.sqrt(n)
    
    while (n > 1) and (factor<=maxFactor):
        if (n % factor) == 0:
            n=n/factor
            lastFactor=factor
            while (n % factor) == 0:
                primfac.append(d)  # supposing you want multiple factors repeated
                n = n / factor
            maxFactor= math.sqrt(n)
        factor += 2
    if n == 1:
       return (lastFactor, primfac)
    else:
        return (n, primfac)
