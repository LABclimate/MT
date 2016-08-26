# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 17:03:49 2016

@author: levyn
"""
import matplotlib.pyplot as plt
import matplotlib as ml
import CESM_utils_plt as utils_plt

plt.ion() # enable interactive mode
path_fig = '../figures/160826/'

# data choice
dMOC_auxgrd_V = dMOC_auxgrd_V_0
dMOC_mgrd_V = dMOC_mgrd_V_0

# y-axes
ax_vol = ax_vol_glob
ticks_vol = ticks_vol_glob
# x-axes
ax_long = np.arange(len(lat_mgrd)) # lat_mgrd

# scaling
min = -80
max = 80
nlevels_col=81
levels_cont = np.linspace(min, max, 9)
print('levels: ' +str(levels_cont))

# dMOC_mgrid_V (in Sv)
# --stretch---------------------------------------
fig, ax = utils_plt.plot_MOC(ax_long, ax_vol, dMOC_mgrd_V, min, max, nlevels_col, levels_cont, plttype='pcolor+contour', to_subplot=[2,2,1])
plt.xlim([0,73])#[-36,73])
plt.yticks(ticks_vol)
plt.gca().set_yticklabels(ticks_dens)
# --linear----------------------------------------
fig, ax = utils_plt.plot_MOC(ax_long, dbc, dMOC_mgrd_V, min, max, nlevels_col, levels_cont, plttype='pcolor+contour', to_newfig=False, to_subplot=[2,2,2])
plt.xlim([0,73])#[-36,73])
# dMOC_auxgrd_V (in Sv)
# --stretch---------------------------------------
fig, ax = utils_plt.plot_MOC(lat_auxgrd, ax_vol, dMOC_auxgrd_V, min, max, nlevels_col, levels_cont, plttype='pcolor+contour', to_newfig=False, to_subplot=[2,2,3])
plt.xlim([-36,90])
plt.yticks(ticks_vol)
plt.gca().set_yticklabels(ticks_dens)
# --linear----------------------------------------
fig, ax = utils_plt.plot_MOC(lat_auxgrd, dbc, dMOC_auxgrd_V, min, max, nlevels_col, levels_cont, plttype='pcolor+contour', to_newfig=False, to_subplot=[2,2,4])
plt.xlim([-36,90])
# ------------------------------------------------
plt.suptitle('METHOD 0 - without proj:   dMOC mgrd V  ----  dMOC auxgrd V')
#utils_plt.print2pdf(fig, path_fig+'dMOC_V_METHOD_0_nonproj')