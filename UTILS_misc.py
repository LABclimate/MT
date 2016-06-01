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
        print('\b.]  Done!')

'''
# code for testing
import time
for i in np.arange(21):
    time.sleep(.01)
    ProgBar('step', step = i, nsteps = 21)
utils_misc.ProgBar('done')
'''

# --- load variable using pickle
def loadvar(filename):
    print('> loading "' + filename + '" from file... ')
    sys.stdout.flush()    
    with open(filename, 'rb') as f: 
      var = pickle.load(f)
    print(' --> Success!')
    return(var)

# --- save variable using pickle
def savevar(var, filename):
    print('> saving "' + filename + '" to file... ')
    sys.stdout.flush()    
    with open(filename, 'wb') as f:
      pickle.dump(var, f)
    print(' --> Success!')      

# --- add directory if unexistant
def checkdir(dirname):
    if os.path.isdir(dirname)==False:
      os.mkdir(path_vars)
      print(' -> created new directory:' + path_vars)


#################################################################################################
# COLLECTION OF UNTESTED FUNCTIONS 
#################################################################################################

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
