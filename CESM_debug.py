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