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
