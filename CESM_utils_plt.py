#################################
# The CESM python toolbox at KUP
# ------ Plotting-Toolbox ------
#################################
# Usage (add the following to your python scripts):
# sys.path.append('/home/buerki/Documents/MT/scripts/')
# import CESM_utils_plt as utils_plt
#################################
# contained functions:
#################################
# - plot_BSF()
# - plot_MOC()
#   - emptyBasemapFig()
#   - shiftCMap()
#   - get_cmap()            # not used anymore
#     - get_viridis()       # not used anymore
# - print2pdf()
#################################
# please log your changes below
#################################
# 28-Apr-2016 - buerki@climate.unibe.ch : created this toolbox
#                                         added makeBasemapFig()
# 29-Apr-2016 - buerki@climate.unibe.ch : added get_cmap()
#                                         added get_viridis() (/alphadata02/born/CESM_utils.py)
# 30-Apr-2016 - buerki@climate.unibe.ch : renamed makeBasemapFig() to emptyBasemapFig()
#                                         added print2pdf()
#                                         added pcolor_basemap()
# 14-Jun-2016 - buerki@climate.unibe.ch : major renaming of function names
#                                         added shiftCMap()
#################################

import numpy as np
import matplotlib as ml
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.colors import from_levels_and_colors  # used for viridis colormap
from matplotlib.backends.backend_pdf import PdfPages
import CESM_utils_mask as utils_mask
import CESM_utils_plt as utils_plt
import CESM_utils_conv as utils_conv

# =======================================================================================
#  Styling (general styling for any plot)
# =======================================================================================
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
ml.style.use('bmh')

# =======================================================================================
#  Plotting
# =======================================================================================

# wrapper to draw contour(f) or pcolorplot of any lat-depth/density variable
def plot_MOC(xvar, yvar, var, nlevels=100, plttype='contour', min = [], max = [], to_newfigure=True, to_subplot=[0,0,0]):
    ''' uses utils_plt.get_cmap()'''
    # colormap
    if min == []:   min = np.nanmin(var) # minimum value of varin
    if max == []:   max = np.nanmax(var)  # maximum value of varin
     #cmap, norm = utils_plt.get_cmap(min, max, nlevels, scheme=utils_plt.get_viridis()) # viridis
    cmap = utils_plt.shiftCMap(ml.cm.seismic, midpoint=1-max/(max-min), name='shifted') # shifted blue white red
    # open new figure
    if to_newfigure == True:
        fig = plt.figure()
    # subplotting
    if np.sum(to_subplot) > 0:
        plt.subplot(to_subplot[0],to_subplot[1],to_subplot[2])
    # draw plot
    if plttype == 'contourf':
        ax = plt.contourf(xvar, yvar, var, cmap=cmap, levels=np.linspace(min, max, nlevels))
    elif plttype == 'contour':
        ax = plt.contaxur(xvar, yvar, var, cmap=cmap, levels=np.linspace(min, max, nlevels))
    elif plttype == 'pcolor':
        ax = plt.pcolor(xvar, yvar, var, cmap=cmap)
    elif plttype == 'pcolor+contour':
        ax = plt.pcolor(xvar, yvar, var, cmap=cmap)
        plt.contour(xvar, yvar, var, colors='black', levels=np.linspace(min, max, nlevels))
    
    plt.colorbar(ax)
    plt.gca().invert_yaxis()
   # plt.gca().set_yscale("log")
   ## overwrite ticks (sofar only indices) with given values
   # if xticks:
   #     plt.xticks(xticks[[int(i) for i in plt.xticks()[0][:-1]]])
   # if yticks:
   #     plt.yticks(yticks[[int(i) for i in plt.yticks()[0][:-1]]])

    plt.text(-10,600000, 'max: {}, min: {}'.format(round(np.nanmax(var),2), round(np.nanmin(var),2)))
    
    if to_newfigure == True:    return(fig, ax)
    else:                       return(ax)

# wrapper to draw pcolorplot of any lat-lon variable using Basemap
def plot_BSF(var, TorUgrid, nlevels=100, mappingtoolbox='basemap', proj='ortho', min = [], max = []):
    ''' uses utils_plt.emptyBasemapFig()'''
    # colormap
    if min == []:   min = np.floor(np.nanmin(var)) # minimum value of varin
    if max == []:   max = np.ceil(np.nanmax(var))  # maximum value of varin
    #cmap, norm = utils_plt.get_cmap(min, max, nlevels, scheme=utils_plt.get_viridis()) # viridis
    cmap = utils_plt.shiftCMap(ml.cm.seismic, midpoint=1-max/(max-min), name='shifted') # shifted blue white red

    # add cyclic border (to prevent gap in pcolor plot)
    var = utils_conv.add_cyclic(var)
    # choose U or T grid
    if TorUgrid == 'U':
        xvar = var.ULONG.values
	yvar = var.ULAT.values
    elif TorUgrid == 'T':
        xvar = var.TLONG.values
        yvar = var.TLAT.values
    # draw plot in new figure
    if mappingtoolbox == 'basemap':
        fig, map = utils_plt.emptyBasemapFig(proj)
        c1 = map.pcolor(xvar,yvar,var.values,latlon=True, cmap=cmap, rasterized=True)
        map.colorbar()
    elif mappingtoolbox == 'cartopy':
        fig, map = True, True#TODO
    return(fig, map)

# wrapper to prepare a Basemap plot in a new figure
def emptyBasemapFig(projection='ortho'):
    fig = plt.figure()    # initialization
    # projection and extent
    if projection == 'ortho':
        map = Basemap(projection='ortho',
        lat_0=40, lon_0=-20,
        lat_ts=20, resolution='c')
    elif projection == 'merc':
        map = Basemap(projection='merc',
        llcrnrlat=0, urcrnrlat=80, llcrnrlon=-110, urcrnrlon=30, 
        lat_ts=20, resolution='c')
    # add grid, coastlines, ...
    map.drawmapboundary()
    map.drawcoastlines()
    p1 = map.drawparallels(np.arange(-80.,81,10.),labels=[True,False,False,False])
    p2 = map.drawmeridians(np.arange(0.,360.,30.),labels=[False,False,True,True])
    map.fillcontinents(color='DarkGray',lake_color='LightGray')
    return(fig, map)

# =======================================================================================
#  Printing
# =======================================================================================

# print figure to pdf
def print2pdf(fig, filename):
    pp = PdfPages(filename + '.pdf')
    pp.savefig(fig)
    pp.close()

# =======================================================================================
#  Color Maps
# =======================================================================================

# shift colormap such that zero is centered
def shiftCMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    source: http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
    
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input:
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = { 'red': [], 'green': [], 'blue': [], 'alpha': [] }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = ml.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap
 