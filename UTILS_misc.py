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
# - loadgetsave()
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
# 19-Aug-2016 - buerki@climate.unibe.ch: created loadgetsave() and getnsave()
#################################

import numpy as np
import sys
import pickle
import os
import UTILS_misc as utils_misc

# --- Progress Bars ---
def ProgBar(stage, step=99, nsteps=99,  maxlen = 60, forceinit=False, minbarlen=1):
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

# --- save variable using pickle
def savevar(var, path_to_var, str_varname='__', verbose=True):
    if verbose: print(' > saving {} to file...'.format(str_varname))
    sys.stdout.flush()    
    with open(path_to_var, 'wb') as f:
      pickle.dump(var, f)
    if verbose: print(' > success!')

# --- try to load, else create and save
def loadgetsave(fun_to_get_var, path_to_var, str_varname='__', verbose=True, noload=False):
    ''' 
    Usage:
     > import utils_misc.loadgetsave as LGS
     > LGS(lambda: fun_to_get_var(args), path_to_var)
    '''
    try:
        if noload: raise ValueError()
        if verbose: print(' > trying to load {}\n    from {}...'.format(str_varname, path_to_var))
        var = utils_misc.loadvar(path_to_var, verbose=False)
        if verbose: print(' > success!\n')
    except:
        if noload & verbose: print(' > no load (note that an eventually stored variable will be overwritten)\n > calculating {}...'.format(str_varname))
        elif verbose: print(' > failed to load!\n > calculating {}...'.format(str_varname))
        var = fun_to_get_var()
        if verbose: print(' > success!\n > saving {} to file...'.format(str_varname))
        utils_misc.savevar(var, path_to_var, verbose=False)
        if verbose: print(' > success!\n')
    return(var)

# --- create variable and save
def getnsave(fun_to_get_var, path_to_var):
    ''' 
    Usage: 
     > getnsave(lambda: fun_to_get_var(args), path_to_var)
    '''
    var = fun_to_get_var()
    utils_misc.savevar(var, path_to_var)
    return(var)
    
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
