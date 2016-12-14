# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 17:03:49 2016

@author: levyn
"""

import matplotlib.pyplot as plt
import matplotlib as ml
import CESM_utils_plt as utils_plt

plt.ion() # enable interactive mode
path_fig = '../figures/{}/'.format(dir_DB)
utils_misc.mkdir(dir_dens)

# colorbar scaling
min = -60
max = 60
nlevels_col=31

# contour scaling
levels_cont = np.linspace(-80, 80, 9)
print('levels: ' +str(levels_cont))

 #mutation = 'deepest10nan'MVxint_mgrd     = np.nansum(MV, axis=2)
mutation = 'none' #'noGulfMex3'

# Choice of data/method
method = 'B'
dMOCm = dMOC_mgrd_V_B
dMOCaux = dMOC_auxgrd_V_B
# Choice of data/method
method = '0'
dMOCm = dMOC_mgrd_V_0
dMOCaux = dMOC_auxgrd_V_0

# x-axes
lat_mgrd = ncdat.TLAT.isel(nlon=0)          # mean of LAT for each j #! very inappropriate
ax_long = np.arange(len(lat_mgrd)) # lat_mgrd

# y-axes
ax_vol = ax_vol_glob
ticks_vol = ticks_vol_glob


# dMOC_mgrd_V (in Sv)
# --stretch---------------------------------------
fig, ax = utils_plt.plot_MOC(ax_long, ax_vol, dMOCm, min, max, nlevels_col, levels_cont, plttype='pcolor+contour', to_subplot=[2,2,1])
#plt.xlim([-36,73])
plt.title('dMOC mgrd (stretch)')
plt.xlim([81,383])
plt.yticks(ticks_vol)
plt.gca().set_yticklabels(ticks_dens_rd)
# --linear----------------------------------------
fig, ax = utils_plt.plot_MOC(ax_long, DBc, dMOCm, min, max, nlevels_col, levels_cont, plttype='pcolor+contour', to_newfig=False, to_subplot=[2,2,2])
plt.title('dMOC mgrd (linear)')
plt.xlim([81,383])
plt.ylim([38,30])

# dMOCaux (in Sv)
# --stretch---------------------------------------
fig, ax = utils_plt.plot_MOC(lat_auxgrd, ax_vol, dMOCaux, min, max, nlevels_col, levels_cont, plttype='pcolor+contour', to_newfig=False, to_subplot=[2,2,3])
plt.title('dMOC auxgrd (stretch)')
plt.xlim([-36,90])
plt.yticks(ticks_vol)
plt.gca().set_yticklabels(ticks_dens_rd)
# --linear----------------------------------------
fig, ax = utils_plt.plot_MOC(lat_auxgrd, DBc, dMOCaux, min, max, nlevels_col, levels_cont, plttype='pcolor+contour', to_newfig=False, to_subplot=[2,2,4])
plt.title('dMOC auxgrd (linear)')
plt.xlim([-36,90])
plt.ylim([38,30])
# ------------------------------------------------
max_auxgrd = round(np.nanmax(dMOCaux),2)
min_auxgrd = round(np.nanmin(dMOCaux),2)
max_mgrd = round(np.nanmax(dMOCm),2)
min_mgrd = round(np.nanmin(dMOCm),2)
plt.suptitle('CESM V'+str(CESMversion)+' | Method {}'.format(method)+' | Mutation {}'.format(mutation)+\
'\n'+str2_db+\
'\nrange mgrd: ({}, {}) | range auxgrd: ({}, {})\n'.format(min_mgrd, max_mgrd, min_auxgrd, max_auxgrd)+\
'levels: {}'.format(levels_cont))

plt.savefig(path_fig+'dMOC_V_CESMV'+str(CESMversion)+'_M'+method+'_'+mutation)
#utils_plt.print2pdf(fig, path_fig+'dMOC_V_CESMV'+str(CESMversion)+'_M'+method)
