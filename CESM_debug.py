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




# --- some statistics
    print('Statistics on ___:' \
          '\n mean:   {}\n median: {}\n min:    {}\n max:    {}'.format(\
          np.nanmean(foo), np.nanmedian(foo), \
          np.nanmin(foo), np.nanmax(foo)))
          
          
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