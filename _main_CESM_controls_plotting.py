'''
CESM controls --> plotting routines

@author: buerki@climate.unibe.ch
TODO: 	 add '.values' where possible to speed up code.

'''
import matplotlib.pyplot as plt
import matplotlib as ml
import CESM_utils_plt as utils_plt

plt.ion() # enable interactive mode
path_fig = '../figures/post_discussion_160701/'


# =======================================================================================
#  CCSM4 representations
# =======================================================================================
# -----------------------------------------------------------------------------------------
# BSF on geographical grid calculated by model
BSF_model = utils_mask.mask_ATLANTIC(ncdat.BSF.isel(time=0), ncdat.REGION_MASK)
fig, map = utils_plt.plot_BSF(BSF_model, 'T', nlevels = 10)
plt.title('BSF model on T grid')
utils_plt.print2pdf(fig, path_fig+'BSF_model_T')
# -----------------------------------------------------------------------------------------
# MOC on geographical grid calculated by model
MOC_model = ncdat.MOC.isel(time=0, transport_reg=1, moc_comp=0)
#MOC_model = MOC_model - MOC_model[:,-1] # normalisation
fig, ax = utils_plt.plot_MOC(MOC_model.lat_aux_grid, MOC_model.moc_z, MOC_model, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HT_auxgrd_xmax) 				# plot seafloor
plt.title('MOC model')
plt.xlim([-36,90])
utils_plt.print2pdf(fig, path_fig+'MOC_model')

# =======================================================================================
#  Calculated on model grid
# =======================================================================================
# -----------------------------------------------------------------------------------------
# BSF on model grid
fig, map = utils_plt.plot_BSF(BSF_mgrd, 'U', nlevels=100)
plt.title('BSF mgrd on U grid')
utils_plt.print2pdf(fig, 'testfigures/BSF_mgrd_U')
# -----------------------------------------------------------------------------------------
# MOC_mgrd_W
fig, ax = utils_plt.plot_MOC(lat_mgrd, ncdat.z_w_top, MOC_mgrd_W, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_mgrd,HT_mgrd_xmax) 				# plot seafloor #! it's the T-grid!!!
plt.title('MOC mgrd W')
plt.xlim([-36,90])
utils_plt.print2pdf(fig, path_fig+'MOC_mgrd_W')
# -----------------------------------------------------------------------------------------
# MOC_mgrd_V
fig, ax = utils_plt.plot_MOC(lat_mgrd, ncdat.z_t, MOC_mgrd_V, nlevels=100, plttype='pcolor+contour')
plt.plot(lat_mgrd,HT_mgrd_xmax)				# plot seafloor
plt.title('MOC mgrd V')
plt.xlim([-36,90])
 #utils_plt.print2pdf(fig, path_fig+'MOC_mgrd_V')
# -----------------------------------------------------------------------------------------
# dMOC on model grid (in Sv)
fig, ax = utils_plt.plot_MOC(lat_mgrd, dens_bins_centers, dMOC_mgrd_W_norm, nlevels=10, plttype='pcolor+contour')
plt.title('dMOC mgrd W (sigma2)')
plt.suptitle('density binning from {} to {} in {} steps'.format(dens_bins_centers.min(), dens_bins_centers.max(), len(dens_bins_centers)))
plt.xlim([-36,73])
utils_plt.print2pdf(fig, path_fig+'dMOC_mgrd_W_sig2')

# =======================================================================================
#  Calculated on auxiliary (geographical) grid
# =======================================================================================
# -----------------------------------------------------------------------------------------
# MOC_auxgrd_W
fig, ax = utils_plt.plot_MOC(lat_auxgrd, z_w_top_auxgrd, MOC_auxgrd_W, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
plt.title('MOC auxgrd W')
utils_plt.print2pdf(fig, path_fig+'MOC_auxgrd_W')
# -----------------------------------------------------------------------------------------
# MOC_auxgrd_V
fig, ax = utils_plt.plot_MOC(lat_auxgrd, zT_auxgrd, MOC_auxgrd_V_norm, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HU_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
plt.title('MOC auxgrd V')
 #utils_plt.print2pdf(fig, 'testfigures/MOC_auxgrd_V')
# -----------------------------------------------------------------------------------------
# dMOC_auxgrd_W (in Sv)
fig, ax = utils_plt.plot_MOC(lat_auxgrd, dens_bins_centers, dMOC_auxgrd_W_norm, nlevels=10, plttype='pcolor+contour')
plt.title('dMOC auxgrd W (sigma2)')
plt.suptitle('density binning from {} to {} in {} steps'.format(dens_bins_centers.min(), dens_bins_centers.max(), len(dens_bins_centers)))
plt.xlim([-36,90])
utils_plt.print2pdf(fig, path_fig+'dMOC_auxgrd_W_sig2')

# -----------------------------------------------------------------------------------------
# MWxint_auxgrd
fig, ax = utils_plt.plot_MOC(lat_auxgrd, z_w_top_auxgrd, MWxint_auxgrd, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
 #utils_plt.print2pdf(fig, 'testfigures/MWxint_auxgrd')
# -----------------------------------------------------------------------------------------
# MVxint_auxgrd
fig, ax = utils_plt.plot_MOC(lat_auxgrd, zT_auxgrd, MVxint_auxgrd, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HU_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
plt.title('MVxint auxgrd')
 #utils_plt.print2pdf(fig, 'testfigures/MVxint_auxgrd')


# -----------------------------------------------------------------------------------------
# COMBINATION of MWxint_auxgrd and MOC_auxgrd
fig = plt.figure()
plt.subplot(3,1,1)
ax = utils_plt.plot_MOC(lat_auxgrd, z_w_top_auxgrd, MOC_auxgrd_W_norm, nlevels=10, plttype='pcolor+contour', to_newfigure=False)
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.title('MOC on depth axis on auxgrd')
plt.suptitle('density binning from {} to {} in {} steps'.format(dens_bins_centers.min(), dens_bins_centers.max(), len(dens_bins_centers)))
plt.xlim([-36,90])

plt.subplot(3,1,2)
ax = utils_plt.plot_MOC(lat_auxgrd, z_w_top_auxgrd, MWxint_auxgrd, nlevels=10, plttype='pcolor+contour', to_newfigure=False)
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
plt.title('MW on depth axis longitudinally integrated on auxgrd (in Sv)')
plt.ylabel('depth')

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
plt.subplot(3,1,1)
ax = utils_plt.plot_MOC(lat_auxgrd, dens_bins[:-1], dMOC_auxgrd_W_norm, nlevels=10, plttype='pcolor+contour', to_newfigure=False)
plt.title('MOC on density axis on auxgrd')
plt.suptitle('density binning from {} to {} in {} steps'.format(dens_bins_centers.min(), dens_bins_centers.max(), len(dens_bins_centers)))
plt.xlim([-36,90])

plt.subplot(3,1,2)
ax = utils_plt.plot_MOC(lat_auxgrd, dens_bins[:-1], dMWxint_auxgrd, nlevels=10, plttype='pcolor+contour', to_newfigure=False)
plt.xlim([-36,90])
plt.plot(lat_mgrd, np.nanmax(np.nanmax(sig2,0),1), 'm-', label='maximal density (on mgrd)')
plt.title('MW on density axis longitudinally integrated on auxgrd (in Sv)')
plt.ylabel('density')
plt.legend(loc='lower left')

plt.subplot(3,1,3)
plt.plot(lat_auxgrd, np.nansum(dMWxint_auxgrd,axis=0), '.-k')
plt.colorbar()
plt.xlim([-36,90])
plt.ylabel('sum over whole density-axis (in Sv)')
plt.xlabel('latitude')
 #utils_plt.print2pdf(fig, 'testfigures/dMWxint_auxgrd')




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


# -----------------------------------------------------------------------------
# Temp and MWxint against density rsp. density-bins in order to directly compare dens versus depth space
# -----------------------------------------------------------------------------
lat = 85
lon = 50
plt.figure()
plt.subplot(1,2,1)
plt.suptitle('location: lat={:d} lon={:d}'.format(lat, lon))
ax = plt.gca(); ax.invert_yaxis()
plt.plot(T_dens[:,lat,lon], dens_bins,'.-r', label='density-T against density bins')
plt.plot(T[:,lat,lon], sig2[:,lat,lon],'.-b', label='depth-T against sigma2')
plt.title('Temperature')
plt.legend(loc='best')

plt.subplot(1,2,2)
ax = plt.gca(); ax.invert_yaxis()
plt.plot(MW_dens[:,lat,lon], dens_bins[:-1],'.-r', label='density-MW against density bins')
plt.plot(MW_mgrd[:,lat,lon], sig2[:,lat,lon],'.-b', label='depth-MW against sigma2')
plt.title('MW')
plt.legend(loc='best')


# -----------------------------------------------------------------------------
# horizontal plot of Temperatures columnwise integrated over depth/density axis
# -----------------------------------------------------------------------------

# - TEMPERATURE
plt.figure()

T_int = utils_ana.integrate_along_dens(T, ncdat.dz)
T_dens_int = utils_ana.integrate_along_dens(T_dens[1:-1], ddb)
min_T = np.nanmin([T_int, T_dens_int])
max_T = np.nanmax([T_int, T_dens_int])
cmap_T = utils_plt.shiftCMap(ml.cm.seismic, midpoint = 1-max_T/(max_T-min_T), name='shifted')

# - plotting
plt.suptitle('Values integrated over columns')
plt.subplot(1,2,1)
plt.imshow(utils_conv.rollATL(T_int), cmap=cmap_T, clim = [min_T, max_T])
plt.colorbar()
#plt.contour(utils_conv.rollATL(T_int), cmap='Reds', levels=np.linspace(-5, 30, 5))
plt.contour(utils_conv.rollATL(ncdat.HT), levels=np.linspace(0,560000,2), cmap='hot') # draw continents
plt.title('Temperature on depth axis')
ax = plt.gca(); ax.invert_yaxis()
print(np.nanmin(T_int), np.nanmax(T_int))

plt.subplot(1,2,2)
plt.imshow(utils_conv.rollATL(T_dens_int), cmap=cmap_T, clim = [min_T, max_T])
plt.colorbar()
#plt.contour(utils_conv.rollATL(T_dens_int), cmap='Reds', levels=np.linspace(-5, 30, 5))
plt.contour(utils_conv.rollATL(ncdat.HT), levels=np.linspace(0,560000,2), cmap='hot') # draw continents
plt.title('Temperature on density axis')
ax = plt.gca(); ax.invert_yaxis()
print(np.nanmin(T_dens_int), np.nanmax(T_dens_int))

# - Transports
plt.figure()

MW_int = utils_ana.integrate_along_dens(MW_mgrd, ncdat.dzw)
MW_dens_int = utils_ana.integrate_along_dens(MW_dens, ddb_centers)
min_MW = np.nanmin([MW_int, MW_dens_int])
max_MW = np.nanmax([MW_int, MW_dens_int])
cmap_MW = utils_plt.shiftCMap(ml.cm.seismic, midpoint = 1-max_MW/(max_MW-min_MW), name='shifted')

plt.subplot(1,2,1)
plt.imshow(utils_conv.rollATL(MW_int), cmap=cmap_MW, clim = [min_MW, max_MW])
plt.colorbar()
#plt.contour(utils_conv.rollATL(MW_int), cmap='Reds', levels=np.linspace(-5, 30, 5))
plt.contour(utils_conv.rollATL(ncdat.HT), levels=np.linspace(0,560000,2), cmap='hot') # draw continents
plt.title('MW on depth axis')
ax = plt.gca(); ax.invert_yaxis()
print(np.nanmin(MW_int), np.nanmax(MW_int))

plt.subplot(1,2,2)
plt.imshow(utils_conv.rollATL(MW_dens_int), cmap=cmap_MW, clim = [min_MW, max_MW])
plt.colorbar()
#plt.contour(utils_conv.rollATL(MW_dens_int), cmap='Reds', levels=np.linspace(-5, 30, 5))
#plt.contour(utils_conv.rollATL(ncdat.HT), levels=np.linspace(0,560000,2), cmap='hot') # draw continents
plt.title('MW on density axis')
ax = plt.gca(); ax.invert_yaxis()
print(np.nanmin(MW_dens_int), np.nanmax(MW_dens_int))




# -----------------------------------------------------------------------------
# horizontal plot of Temperatures columnwise integrated over depth/density axis
# -----------------------------------------------------------------------------

# - TEMPERATURE
plt.figure()

T_int = utils_ana.integrate_along_dens(T, ncdat.dz)
T_dens_int = utils_ana.integrate_along_dens(T_dens[1:-1], ddb)
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
MW_dens_int = utils_ana.integrate_along_dens(MW_dens, ddb_centers)
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
















# =======================================================================================
#  Other variables
# =======================================================================================
# -----------------------------------------------------------------------------------------
# Seafloor
fig, ax = plt.contourf(ncdat.HT.roll(nlon=54), levels=np.linspace(0,560000,100))
plt.title('Depth of Seafloor')
fig, map = utils_plt.pcolor_basemap(MW_mgrd.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
fig, map = utils_plt.pcolor_basemap(MW_mgrd.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
 #utils_plt.print2pdf(fig, 'testfigures/seafloor')
# -----------------------------------------------------------------------------------------
# Temperature
fig, map = utils_plt.pcolor_basemap(T.roll(nlon=54), 'T', nlevels=100)
plt.title('Temperature')
 #utils_plt.print2pdf(fig, 'testfigures/temp')
# -----------------------------------------------------------------------------------------
# Density
fig, ax = utils_plt.plot_MOC(ncdat.TLAT[:,0], ncdat.z_t, np.nanmean(sig2, 2), nlevels=10, plttype='contour')
plt.title('Density')
 #utils_plt.print2pdf(fig, 'density/temp')
# -----------------------------------------------------------------------------------------
# MW_mgrd
fig, map = utils_plt.pcolor_basemap(MW_mgrd.roll(nlon=54).mean(dim='z_w_top'), 'T', nlevels=100,)
plt.title('MW_mgrd')
 #utils_plt.print2pdf(fig, 'testfigures/MW_mgrd')
# -----------------------------------------------------------------------------------------
# ANGLE
plt.figure()
plt.pcolor(ncdat.ANGLE*180/np.pi)
plt.contour(ncdat.REGION_MASK)

fig, ax = utils_plt.plot_MOC(ncdat.VVEL.TLAT[:,0], z_w_top_auxgrd, MV_projauxgrd[:,:,59])
fig, ax = utils_plt.plot_MOC(ncdat.VVEL.TLAT[:,0], z_w_top_auxgrd, ncdat.VVEL[0,:,:,59])