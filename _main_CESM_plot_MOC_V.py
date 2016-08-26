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


# =======================================================================================
#  Calculated on model grid
# =======================================================================================
# -----------------------------------------------------------------------------------------
# MOC_mgrd_V
fig, ax = utils_plt.plot_MOC(lat_mgrd, ncdat.z_t, MOC_mgrd_V, nlevels=100, plttype='pcolor+contour')
plt.plot(lat_mgrd,HT_mgrd_xmax)				# plot seafloor
plt.title('MOC mgrd V')
plt.xlim([-36,90])
 #utils_plt.print2pdf(fig, path_fig+'MOC_mgrd_V')
# -----------------------------------------------------------------------------------------
# dMOC_mgrid_V (in Sv)
fig, ax = utils_plt.plot_MOC(lat_mgrd, dbc, dMOC_mgrd_V, nlevels=10, plttype='pcolor+contour')
plt.title('dMOC mgrd W (sigma2)')
plt.suptitle('density binning from {} to {} in {} steps'.format(dbc.min(), dbc.max(), len(dbc)))
plt.xlim([-36,73])
plt.yticks(ticks_vol_glob)
plt.gca().set_yticklabels(ticks_dens)
#utils_plt.print2pdf(fig, path_fig+'dMOC_mgrd_V_sig2')

# =======================================================================================
#  Calculated on auxiliary (geographical) grid
# =======================================================================================
# -----------------------------------------------------------------------------------------
# MOC_auxgrd_V
fig, ax = utils_plt.plot_MOC(lat_auxgrd, zT_auxgrd, MOC_auxgrd_V, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HU_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
plt.title('MOC auxgrd V')
 #utils_plt.print2pdf(fig, 'testfigures/MOC_auxgrd_V')
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# dMOC_auxgrd_V (in Sv)
fig, ax = utils_plt.plot_MOC(lat_auxgrd, dbc, dMOC_auxgrd_V, nlevels=10, plttype='pcolor+contour')
plt.title('dMOC auxgrd V (sigma2)')
plt.suptitle('density binning from {} to {} in {} steps'.format(dbc.min(), dbc.max(), len(dbc)))
plt.xlim([-36,90])
plt.yticks(ticks_vol_reg)
plt.gca().set_yticklabels(ticks_dens)
#utils_plt.print2pdf(fig, path_fig+'dMOC_auxgrd_V_sig2')

# -----------------------------------------------------------------------------------------
# MVxint_auxgrd
fig, ax = utils_plt.plot_MOC(lat_auxgrd, zT_auxgrd, MVxint_auxgrd, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HU_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
plt.title('MVxint auxgrd')
 #utils_plt.print2pdf(fig, 'testfigures/MVxint_auxgrd')

