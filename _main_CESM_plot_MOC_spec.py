'''
CESM controls --> plotting routines for MOC calculated from W transports

@author: buerki@climate.unibe.ch
TODO: 	 add '.values' where possible to speed up code.

'''
import matplotlib.pyplot as plt
import matplotlib as ml
import CESM_utils_plt as utils_plt

plt.ion() # enable interactive mode
path_fig = '../figures/160711/'




# -----------------------------------------------------------------------------------------
# COMBINATION of MWxint_auxgrd and MOC_auxgrd
fig = plt.figure()
plt.suptitle('density binning from {} to {} in {} steps'.format(dbc.min(), dbc.max(), len(dbc)))

plt.subplot(3,1,1)
plt.title('MOC on depth axis on auxgrd')
ax = utils_plt.plot_MOC(lat_auxgrd, z_w_auxgrd, MOC_auxgrd_W, nlevels=10, plttype='pcolor+contour', to_newfigure=False)
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])

plt.subplot(3,1,2)
plt.title('MW on depth axis longitudinally integrated on auxgrd (in Sv)')
plt.ylabel('depth')
ax = utils_plt.plot_MOC(lat_auxgrd, z_w_auxgrd, MWxint_auxgrd, nlevels=10, plttype='pcolor+contour', to_newfigure=False)
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])

plt.subplot(3,1,3)
plt.plot(lat_auxgrd, np.nansum(MWxint_auxgrd,axis=0), '.-k')
plt.colorbar()
plt.xlim([-36,90])
plt.ylabel('sum over whole density-axis (in Sv)')
plt.xlabel('latitude')
 #utils_plt.print2pdf(fig, 'testfigures/dMWxint_auxgrd')


# -----------------------------------------------------------------------------------------
# COMBINATION of dMWxint_auxgrd and dMOC_auxgrd
fig = plt.figure()
plt.suptitle('density binning from {} to {} in {} steps'.format(db.min(), db.max(), len(db)))

plt.subplot(3,1,1)
plt.title('MOC on density axis on auxgrd')
ax = utils_plt.plot_MOC(lat_auxgrd, db, dMOC_auxgrd_W, nlevels=10, plttype='pcolor+contour', to_newfigure=False)
#plt.yticks(ticks_vol_reg)
#plt.gca().set_yticklabels(ticks_dens)
plt.xlim([-36,90])

plt.subplot(3,1,2)
plt.title('MW on density axis longitudinally integrated on auxgrd (in Sv)')
plt.ylabel('density')
ax = utils_plt.plot_MOC(lat_auxgrd, db, dMWxint_auxgrd, nlevels=40, plttype='pcolor', to_newfigure=False)
plt.plot(lat_mgrd, np.nanmax(np.nanmax(sig2,0),1), 'm-', label='maximal density (on mgrd)')
#plt.yticks(ticks_vol_reg)
#plt.gca().set_yticklabels(ticks_dens)
plt.xlim([-36,90])
plt.legend(loc='lower left')

plt.subplot(3,1,3)
plt.plot(lat_auxgrd, np.nansum(dMWxint_auxgrd,axis=0), '.-k')
plt.colorbar()
plt.xlim([-36,90])
plt.ylabel('sum over whole density-axis (in Sv)')
plt.xlabel('latitude')
utils_plt.print2pdf(fig, path_fig+'dMOC_tripple'+varname_binning)




# -----------------------------------------------------------------------------
# 5 subplots showing densities, Temperature and MW against depth and density axis
# -----------------------------------------------------------------------------
lat = 85
lon = 30
plt.figure()
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
plt.plot(MW_dens[:,lat,lon], db[:-1],'.-r')
plt.plot(np.zeros_like(db), db, 'b.')
plt.title('MW against density')
#plt.text(.1,.05,'SUM = '+str(np.nansum(MW_dens[:,lat,lon])), horizontalalignment='left', transform=ax.transAxes)

plt.subplot(1,5,4)
ax = plt.gca(); ax.invert_yaxis()
plt.plot(T[:,lat,lon], MW_mgrd.z_w_top, '.-r')
plt.plot(np.zeros_like(MW_mgrd.z_w_top), MW_mgrd.z_w_top, 'b.')
plt.title('Temp against depth')

plt.subplot(1,5,5)
ax = plt.gca(); ax.invert_yaxis()
plt.plot(T_dens[:,lat,lon], db,'.-r')
plt.plot(np.zeros_like(db), db, 'b.')
plt.title('Temp against density')


# -----------------------------------------------------------------------------
# Temp and MWxint against density rsp. density-bins in order to directly compare dens versus depth space
# -----------------------------------------------------------------------------
lat = 85
lon = 51
plt.figure()
plt.subplot(1,2,1)
plt.suptitle('location: lat={:d} lon={:d}'.format(lat, lon))
ax = plt.gca(); ax.invert_yaxis()
plt.plot(T_dens[:,lat,lon], db,'.-r', label='density-T against density bins')
plt.plot(T[:,lat,lon], sig2[:,lat,lon],'.-b', label='depth-T against sigma2')
plt.title('Temperature')
plt.legend(loc='best')

plt.subplot(1,2,2)
ax = plt.gca(); ax.invert_yaxis()
plt.plot(MW_dens[:,lat,lon], db[:-1],'.-r', label='density-MW against density bins')
plt.plot(MW_mgrd[:,lat,lon], sig2[:,lat,lon],'.-b', label='depth-MW against sigma2')
plt.title('MW')
plt.legend(loc='best')



# -----------------------------------------------------------------------------
# horizontal plot of Temperatures columnwise integrated over depth/density axis
# -----------------------------------------------------------------------------

# - TEMPERATURE
plt.figure()

T_int = utils_ana.integrate_along_dens(T, ncdat.dz)
#T_dens_int = utils_ana.integrate_along_dens(T_dens[1:-1], ddb)
T_dens_int = utils_ana.integrate_along_dens(T_dens, vol_dbs_col)
min_T = np.nanmin([T_int])
max_T = np.nanmax([T_int])
cmap_T = utils_plt.shiftCMap(ml.cm.seismic, midpoint = 1-max_T/(max_T-min_T), name='shifted')
min_T_dens = np.nanmin([T_dens_int])
max_T_dens = np.nanmax([T_dens_int])
cmap_T_dens = utils_plt.shiftCMap(ml.cm.seismic, midpoint = 1-max_T_dens/(max_T_dens-min_T_dens), name='shifted')

# - plotting
plt.suptitle('Values integrated over columns')
plt.subplot(1,2,1)
plt.imshow(utils_conv.rollATL(T_int), cmap=cmap_T)
plt.colorbar()
#plt.contour(utils_conv.rollATL(T_int), cmap='Reds', levels=np.linspace(-5, 30, 5))
plt.contour(utils_conv.rollATL(ncdat.HT), levels=np.linspace(0,560000,2), cmap='hot') # draw continents
plt.title('Temperature on depth axis')
ax = plt.gca(); ax.invert_yaxis()
print(np.nanmin(T_int), np.nanmax(T_int))

plt.subplot(1,2,2)
plt.imshow(utils_conv.rollATL(T_dens_int), cmap=cmap_T_dens)
plt.colorbar()
#plt.contour(utils_conv.rollATL(T_dens_int), cmap='Reds', levels=np.linspace(-5, 30, 5))
plt.contour(utils_conv.rollATL(ncdat.HT), levels=np.linspace(0,560000,2), cmap='hot') # draw continents
plt.title('Temperature on density axis')
ax = plt.gca(); ax.invert_yaxis()
print(np.nanmin(T_dens_int), np.nanmax(T_dens_int))

# - Transports
plt.figure()

MW_int = utils_ana.integrate_along_dens(MW_mgrd, ncdat.dzw)
MW_dens_int = utils_ana.integrate_along_dens(MW_dens, ddbc)
min_MW = np.nanmin([MW_int])
max_MW = np.nanmax([MW_int])
cmap_MW = utils_plt.shiftCMap(ml.cm.seismic, midpoint = 1-max_MW/(max_MW-min_MW), name='shifted')
min_MW_dens = np.nanmin([MW_dens_int])
max_MW_dens = np.nanmax([MW_dens_int])
cmap_MW_dens = utils_plt.shiftCMap(ml.cm.seismic, midpoint = 1-max_MW_dens/(max_MW_dens-min_MW_dens), name='shifted')

plt.subplot(1,2,1)
plt.imshow(utils_conv.rollATL(MW_int), cmap=cmap_MW)
plt.colorbar()
#plt.contour(utils_conv.rollATL(MW_int), cmap='Reds', levels=np.linspace(-5, 30, 5))
plt.contour(utils_conv.rollATL(ncdat.HT), levels=np.linspace(0,560000,2), cmap='hot') # draw continents
plt.title('MW on depth axis')
ax = plt.gca(); ax.invert_yaxis()
print(np.nanmin(MW_int), np.nanmax(MW_int))

plt.subplot(1,2,2)
plt.imshow(utils_conv.rollATL(MW_dens_int), cmap=cmap_MW_dens)
plt.colorbar()
#plt.contour(utils_conv.rollATL(MW_dens_int), cmap='Reds', levels=np.linspace(-5, 30, 5))
plt.contour(utils_conv.rollATL(ncdat.HT), levels=np.linspace(0,560000,2), cmap='hot') # draw continents
plt.title('MW on density axis')
ax = plt.gca(); ax.invert_yaxis()
print(np.nanmin(MW_dens_int), np.nanmax(MW_dens_int))



