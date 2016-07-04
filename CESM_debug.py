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


# PLOT the results
plt.figure()

lat = 85
lon = 10
plt.suptitle('location: lat={:d} lon={:d}'.format(lat, lon))

plt.subplot(1,5,1)
ax = plt.gca(); ax.invert_yaxis()
plt.plot(sig2[:,lat,lon], MW_mgrd.z_w_top, '.-r', label='sigma 2')
plt.plot(RHO[:,lat,lon], MW_mgrd.z_w_top, '.-m', label='in-situ density')
plt.title('Dens against depth')
plt.legend(loc='best')
plt.plot(35*np.ones_like(MW_mgrd.z_w_top), MW_mgrd.z_w_top, 'b.')
plt.xlim(20,44)

plt.subplot(1,5,2)
ax = plt.gca(); ax.invert_yaxis()
plt.plot(MW_mgrd[:,lat,lon], MW_mgrd.z_w_top, '.-r')
plt.plot(np.zeros_like(MW_mgrd.z_w_top), MW_mgrd.z_w_top, 'b.')
plt.title('MW against depth')
#plt.text(.1,.05,'SUM = '+str(np.nansum(MW_mgrd[:,lat,lon])), horizontalalignment='left', transform=ax.transAxes)


plt.subplot(1,5,3)
ax = plt.gca(); ax.invert_yaxis()
plt.plot(MW_dens[:,lat,lon], dens_bins[:-1],'.-r')
plt.plot(np.zeros_like(dens_bins), dens_bins, 'b.')
plt.title('MW against density')
#plt.text(.1,.05,'SUM = '+str(np.nansum(MW_dens[:,lat,lon])), horizontalalignment='left', transform=ax.transAxes)

plt.subplot(1,5,4)
ax = plt.gca(); ax.invert_yaxis()
plt.plot(T[:,lat,lon], MW_mgrd.z_w_top, '.-r')
plt.plot(np.zeros_like(MW_mgrd.z_w_top), MW_mgrd.z_w_top, 'b.')
plt.title('Temp against depth')

plt.subplot(1,5,5)
ax = plt.gca(); ax.invert_yaxis()
plt.plot(T_dens[:,lat,lon], dens_bins,'.-r')
plt.plot(np.zeros_like(dens_bins), dens_bins, 'b.')
plt.title('Temp against density')

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