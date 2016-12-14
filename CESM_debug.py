#################################
# The CESM python toolbox at KUP
# ------- Debugging Tools -------
#################################
# Usage: directly in ipython
#################################
# contained functions:
#################################
# - Count number of negative gradients in RHO or sig2 and find the most negative one.
#################################

# -- tic-toc timing
t1 = time.time()
t2 = time.time(); print(t2-t1); t1 = t2;


# --- some statistics
    print('Statistics on ___:' \
          '\n mean:   {}\n median: {}\n min:    {}\n max:    {}'.format(\
          np.nanmean(foo), np.nanmedian(foo), \
          np.nanmin(foo), np.nanmax(foo)))
          
          
# =======================================================================================
# TEST correlation functions
# =======================================================================================
''' indeed, they produce the same result.
    Caveats:    scipy.stats.pearsonr() does not support nans
                pandas.Series.corr() does not return the p-values
    The function utils_ana.nanpearsonr(a,b) solves this problem.
'''
import CESM_utils_analysis as utils_ana
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from numpy.random import rand

a, b = pd.Series(rand(2000)), pd.Series(rand(2000))
c1 = pearsonr(a,b)[0]
p1 = pearsonr(a,b)[1]
c2 = pd.Series.corr(a,b,method='pearson')
c3 = utils_ana.nanpearsonr(a,b)
p3 = utils_ana.ttest_pval(len(a), c3)
cp3 = utils_ana.xpearsonr(a, b, [0], True)

print c1, p1
print c3, p3
print cp3[0][0], cp3[1][0]

a[40] = np.nan
c1nan = pearsonr(a,b)
c2nan = pd.Series.corr(a,b,method='pearson')
cp3nan = utils_ana.xpearsonr(a, b, [0], True)
print cp3nan[0][0], cp3nan[1][0]


''' test lag-direction of utils_ana.xpearsonr() with sine and cosine function
'''
a = np.arange(40)
b = np.cos(a/40.*2*np.pi)
c = np.sin(a/40.*2*np.pi)
#b = (b-b.mean())/b.std() # does not make any difference in corr
#c = (c-c.mean())/c.std() # does not make any difference in corr
lags = np.arange(-20, 21)
corr = utils_ana.xpearsonr(b,c,lags)
plt.figure()
plt.subplot(211); plt.plot(a,b,'.-', a,c,'.-')
plt.subplot(212); plt.plot(lags,corr, 'o-'); plt.xlabel('lag'); plt.ylabel('correaltion coeficient')
plt.show()

# =======================================================================================
# TEST resampling loop
# =======================================================================================

# for one-dimensional grid and data
og = np.array([5,6,7,8,9,10])                                   # old grid
ng = np.array([2,3,5,6.6, 8,10,11])                             # new grid
od = og                                                         # old data
nd = utils_conv.resample_colwise(od, og, ng, method='wmean')    # new data
unsort = [0,1,2,4,3,5]
nd = utils_conv.resample_colwise(od[unsort], og[unsort], ng, method='wmean')

#print(diff_2/diff_total, diff_1/diff_total)
#print(odat[idxo-1], odat[idxo], ndat[idxn])
#print('\n')


# =======================================================================================
# Test densT-->densU and Test utils_ana.canonical_cumsum()
# =======================================================================================

densT = np.array([[1,2,3,4,5,6,7,8,9,10],[2,3,4,5,6,7,8,9,10,11],[10,3,10,1,20,10,10,10,10,10]], dtype=float)
densT = np.expand_dims(densT,axis=0)
densT = np.concatenate((densT, densT*2), axis=0)
densT = densT[:,:,:-3]


# =======================================================================================
# Count number of negative gradients in RHO or sig2 and find the most negative one.
# =======================================================================================

#dens = RHO
dens = sig2

mingradient = 3
numberofneggradients = 0
for j in np.arange(sig2.shape[1]):
    for i in np.arange(sig2.shape[2]):
        diff = np.diff(RHO[:,j,i])
        mingradient = np.nanmin([mingradient, np.nanmin(diff)])
        if any(diff<0):
            numberofneggradients += 1
print(numberofneggradients, mingradient)



# =======================================================================================
# Compare runtimes
# =======================================================================================

del a
del b
del c
t0 = time.clock()
a = 
t1 = time.clock()-t0
t0 = time.clock()
b = 
t2 = time.clock()-t0
t0 = time.clock()
c = 
t3 = time.clock()-t0
print(t1,t2,t3)