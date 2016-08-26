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
# MOC_mgrd_W
fig, ax = utils_plt.plot_MOC(lat_mgrd, ncdat.z_w_top, MOC_mgrd_W, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_mgrd,HT_mgrd_xmax) 				# plot seafloor #! it's the T-grid!!!
plt.title('MOC mgrd W')
plt.xlim([-36,90])
utils_plt.print2pdf(fig, path_fig+'MOC_mgrd_W')

# -----------------------------------------------------------------------------------------
# dMOC_mgrid_W(in Sv)
fig, ax = utils_plt.plot_MOC(lat_mgrd, db, dMOC_mgrd_W, nlevels=10, plttype='pcolor+contour')
plt.title('dMOC mgrd W (sigma2)')
plt.suptitle('density binning from {} to {} in {} steps'.format(dbc.min(), dbc.max(), len(dbc)))
plt.xlim([-36,73])
plt.yticks(ticks_vol)
plt.gca().set_yticklabels(ticks_dens)
#utils_plt.print2pdf(fig, path_fig+'dMOC_mgrd_W_sig2')

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
# dMOC_auxgrd_W (in Sv)
fig, ax = utils_plt.plot_MOC(lat_auxgrd, dbc, dMOC_auxgrd_W, nlevels=10, plttype='pcolor+contour')
plt.title('dMOC auxgrd W (sigma2)')
plt.suptitle('density binning from {} to {} in {} steps'.format(dbc.min(), dbc.max(), len(dbc)))
plt.xlim([-36,90])
plt.yticks(ticks_vol)
plt.gca().set_yticklabels(ticks_dens)
#utils_plt.print2pdf(fig, path_fig+'dMOC_auxgrd_W_sig2')

# -----------------------------------------------------------------------------------------
# MWxint_auxgrd
fig, ax = utils_plt.plot_MOC(lat_auxgrd, z_w_top_auxgrd, MWxint_auxgrd, nlevels=10, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HT_auxgrd_xmax)  				# plot seafloor
plt.xlim([-36,90])
 #utils_plt.print2pdf(fig, 'testfigures/MWxint_auxgrd')
