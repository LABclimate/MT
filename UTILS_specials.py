###############################
# ------ Special Toolbox ------
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import UTILS_specials as utils_spec
#################################
# contained functions:
#################################
# - ProgBar()
#################################
# please log your changes below
#################################
# 14-Mai-2016 - buerki@climate.unibe.ch : created this toolbox
#                                         created ProgBar()
#################################

import numpy as np
import sys
import pickle

# --- Progress Bars ---
def ProgBar(stage, barlen=10, step=0, nsteps=10):
    if stage=='step':
        if step%(np.round(nsteps/barlen))==0:
	    if step==0: # initialization
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
for i in np.arange(384):
    time.sleep(.01)
    utils_spec.ProgBar('step', barlen=32, step = i, nsteps = 384)
utils_spec.ProgBar('done')
'''

# --- load variable using pickle
def loadvar(filename):
    print('> loading "' + filename + '" from file... ')
    with open(filename, 'rb') as f: 
      var = pickle.load(f)
    print(' --> Done!')
    return(var)

# --- save variable using pickle
def savevar(var, filename):
    print('> saving "' + filename + '" to file... ')	
    with open(filename, 'wb') as f:
      pickle.dump(var, f)
    print(' --> Done!')      
