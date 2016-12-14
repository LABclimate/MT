"""
Created on Tue Dec 06 2016

@author:    buerki@climate.unibe.ch
TODO:       

plt.pcolor(b[0]+1, hatch='//', alpha = .2)


"""

# =============================================================================
#  General Preparation
# =============================================================================
import pandas as pd
from UTILS_misc import LGS, GS, LG
from CESM_utils_plt import fg, sp, sps, spsinv, spinv, cb, tt, xl, yl, tight, colyax, invyax
from CESM_utils_conv import rollATL
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, Cursor

from matplotlib.pyplot import legend as lg

# -----------------------------------------------------------------------------
# BUG-FIXES
qMOCidx = MOCidx; qMOCidx[66] = np.nan      #! Bugfix until missing 66 is fixed

# -----------------------------------------------------------------------------
# auxgrid and paths
auxgrd_name = ['lat395model', 'lat198model', 'lat170eq80S90N', 'lat340eq80S90N'][0]       # choose aux grid
path_mgrd = paths.get_path2vars('mgrd', CESMversion=CESMversion, mkdir=True)
path_auxgrd = paths.get_path2vars(auxgrd_name, CESMversion=CESMversion, mkdir=True)
path_figscorr = '../figures/corr/.'

# -----------------------------------------------------------------------------
# reduce to pval>.95, mask and roll for pcolor-plotting

mxcorrBSFidx = xr.DataArray(xcorrBSFidx).to_masked_array(copy=True)
pmxcorrBSFidx = xr.DataArray(xcorrBSFidx).where(pvalxcorrBSFidx <= .05).to_masked_array(copy=True)
mxcorrMOCidx = rollATL(xcorrMOCidx.to_masked_array(copy=True), axis=1)
pmxcorrMOCidx = rollATL(xr.DataArray(xcorrMOCidx).where(pvalxcorrMOCidx <= .05).to_masked_array(copy=True), axis=1)

mMOC = MOC.to_masked_array(copy=True)
mBSF = rollATL(BSF.to_masked_array(copy=True), axis=2)
mBSFsel = rollATL(BSFsel.to_masked_array(copy=True), axis=2)

# -----------------------------------------------------------------------------
# axes
lat_auxgrd = np.array(ncdat.MOC.lat_aux_grid)
lat_BSF = ncdat.BSF.TLAT.mean(axis=1)
long_BSF = rollATL(ncdat.BSF.TLONG.mean(axis=0), axis=0)
zT = np.array(ncdat.MOC.moc_z)

# -----------------------------------------------------------------------------
#  zonal maxima of ocean depth
HT_mgrd_xmax = LGS(lambda: utils_mask.calc_H_mgrd_xmax(ncdat, 'T'), path_mgrd+'HT_mgrd_xmax', 'HT_mgrd_xmax')
HT_auxgrd_xmax = LGS(lambda: utils_mask.calc_H_auxgrd_xmax(lat_auxgrd, ncdat, 'T'), path_auxgrd+'HT_auxgrd_xmax', 'HT_auxgrd_xmax')

# add cyclic
'''
def twice_add_cyc(var, dim1, dim2):
    return utils_conv.add_cyclic(utils_conv.add_cyclic(var, dim1), dim2)
cyc_xcorrMOCidx = twice_add_cyc(xcorrMOCidx, 'dim_0', 'dim_1')
cyc_pvalxcorrMOCidx = twice_add_cyc(pvalxcorrMOCidx, 'dim_0', 'dim_1')
cyc_BSF = twice_add_cyc(BSF, 'nlat', 'nlon')
cyc_BSFsel = twice_add_cyc(BSFsel, 'nlat', 'nlon')

lat_auxgrd = np.hstack([lat_auxgrd,lat_auxgrd[0]])
lat_BSF = np.hstack([lat_BSF,lat_BSF[0]])
long_BSF = np.hstack([long_BSF,long_BSF[0]])
zT = np.hstack([zT,zT[0]])
HT_auxgrd_xmax = np.hstack([HT_auxgrd_xmax,HT_auxgrd_xmax[0]])
'''

# -----------------------------------------------------------------------------
# unit conversions
zT /= 100
HT_mgrd_xmax /= 100
HT_auxgrd_xmax /= 100

# -----------------------------------------------------------------------------
# round value to existing values in array  || for interactive plot click-on
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

# -----------------------------------------------------------------------------
# minima, maxima
fmax = np.zeros([len(lags), 2], dtype=int)
fmin = np.zeros([len(lags), 2], dtype=int)
fmaxsc = np.zeros([len(lags), 2], dtype=int)
fminsc = np.zeros([len(lags), 2], dtype=int)
for idx, n in enumerate(lags[:]):
    fmin[idx,:] = np.where(xcorrBSFidx[:,:,idx] == xcorrBSFidx[:,:,idx].min())[::-1]
    fmax[idx,:] = np.where(xcorrBSFidx[:,:,idx] == xcorrBSFidx[:,:,idx].max())[::-1]
    fminsc[idx,0] = lat_auxgrd[fmin[idx,0]]
    fmaxsc[idx,0] = lat_auxgrd[fmax[idx,0]]
    fminsc[idx,1] = zT[fmin[idx,1]]
    fmaxsc[idx,1] = zT[fmax[idx,1]]

# styling   
col1 = 'crimson'
col2 = 'darkcyan'
prop_txt_on_floor = dict(color='FloralWhite', fontsize=20)
prop_of_min = dict(color='DeepSkyBlue', ms=12, lw=4, label='min')
prop_of_max = dict(color='LightCoral', ms=12, lw=4, label='max')



# =============================================================================
#  Plotting
# =============================================================================
# -----------------------------------------------------------------------------
# plot corrIDX
fg, ax = sps(2,1); ax01 = ax[0].twinx()
p0 = ax[0].plot(-BSFidx, col1); xl(ax[0],'time'); yl(ax[0],'negative BSFidx'); colyax(ax[0], col1)
p1 = ax01.plot(qMOCidx, col2); xl(ax01,'time'); yl(ax01,'MOCidx');   colyax(ax01, col2)
ax[1].plot(lags,corrIDX, 'o-', lw=3); 
tt(ax[1],'Correlation between MOCidx and BSFidx'); xl(ax[1],'Lag [years]'); yl(ax[1],'Correlation coef.');
utils_plt.add_lag_indicators(ax[1])
tight()
# -----------------------------------------------------------------------------
# plot point-wise corrIDX
idxlag0 = np.where(lags==0)[0][0]
fg, ax = sps(2,1);  
p0 = ax[0].pcolor(lat_auxgrd, zT, mxcorrBSFidx[:,:,idxlag0]); cb0 = cb(p0,ax[0]); yl(cb0.ax, 'xcorrBSFidx'); 
p0s = ax[0].fill(foo, fill=False, hatch='\\')
utils_plt.ptype_lat_z(fg, ax[0], lat_auxgrd, HT_auxgrd_xmax);
#---
p1 = ax[1].pcolormesh(long_BSF, lat_BSF,  mxcorrMOCidx[:,:,idxlag0]); cb1 = cb(p1,ax[1]); yl(cb1.ax, 'xcorrMOCidx'); 
xl(ax[1],'Longitude'); yl(ax[1],'Latitude');
tight()
# -----------------------------------------------------------------------------
# plot MOC Streamfunction
times = [0,-1]
fg, ax = sps(2,1);
p = [ax[idx].pcolor(lat_auxgrd, zT, mMOC[time,:,:]) for idx, time in enumerate(times)]
for idx, time in enumerate(times):
    cb0 = cb(p[idx],ax[idx]); yl(cb0.ax,'MOC');
    ax[idx].text(60,5500, 't = {}'.format(time), **prop_txt_on_floor)
    utils_plt.ptype_lat_z(fg, ax[idx], lat_auxgrd, HT_auxgrd_xmax)
    tight()
# -----------------------------------------------------------------------------
# plot BSF Streamfunction
fg, ax = sps(2,1);
p0 = ax[0].pcolor(long_BSF, lat_BSF, mBSF[0,:,:]);  cb0 = cb(p0,ax[0]); yl(cb0.ax,'BSF'); 
tt(ax[0],'t= 0');  utils_plt.ptype_lon_lat(fg, ax[0])
#---
p1 = ax[1].pcolor(long_BSF, lat_BSF, mBSF[-1,:,:]); cb1 = cb(p1,ax[1]); yl(cb1.ax,'BSF'); 
tt(ax[1],'t= -1'); utils_plt.ptype_lon_lat(fg, ax[1])
tight()

# -----------------------------------------------------------------------------
# interactive plot corrBSFidx
# draw plot
global fg, ax, p1, p2, cb, mylag, stat_pvalfilter, myMOCsel
fg, ax = sps(2,1)
fg.subplots_adjust(right=0.85)
mylag=0; myMOCsel = MOC[:,0,0]; stat_pvalfilter=True; # default values for initial plot
xax, yax = np.array(lat_auxgrd), np.array(zT)
def draw_pcolor():
    global fg, ax, p1, cb, mylag, stat_pvalfilter
    fg.delaxes(ax[0]); ax[0]=fg.add_subplot(211);
    try: cb0.remove();
    except: pass
    if stat_pvalfilter: var = pmxcorrBSFidx; varlab='xcorrBSFidx @ pval st .05'
    else: var = mxcorrBSFidx; varlab='xcorrBSFidx'
    idxmylag = np.where(lags==mylag)[0][0]
    p1 = ax[0].pcolor(xax, yax, var[:,:,idxmylag], picker=True, vmin=var.min(), vmax=var.max()); 
    utils_plt.ptype_lat_z(fg, ax[0], lat_auxgrd, HT_auxgrd_xmax)
    cb0 = cb(p1,ax[0]);  yl(cb0.ax,varlab)
    ax[0].fill_between(xax, HT_auxgrd_xmax, ax[0].get_ylim()[0], color='DarkSlateGrey')  # plot seafloor
    plt.draw()
def draw_curve():
    global fg, ax, p2, mylag, myMOCsel
    shiftedvar = pd.Series(myMOCsel).shift(-mylag)
    p2 = ax[1].plot(shiftedvar, col2); yl(ax[1],'MOC'); colyax(ax[1], col2)
    ax[1].set_ylim([np.min(myMOCsel), np.max(myMOCsel)])
    ax3 = ax[1].twinx()
    p3 = ax3.plot(BSFidx, col1); yl(ax3,'BSFidx'); xl(ax3,'Time'); colyax(ax[1],col1)
    ax3.invert_yaxis()
    plt.draw()
def redraw_curve():
    global ax, p2, mylag, myMOCsel
    shiftedvar = pd.Series(myMOCsel).shift(-mylag)
    p2[0].set_ydata(shiftedvar)
    ax[1].set_ylim([np.min(myMOCsel), np.max(myMOCsel)])
    plt.draw()
draw_pcolor()
draw_curve()

# draw widgets for interactivity
axRBlag = plt.axes([.9,.15,.15,.85])
RBlag = RadioButtons(axRBlag, ('-100','-50','-20','-10','-5','-3','-2','-1','0','1','2','3','5','10','20','50','100'), active=8) # active 5 --> 0 lag
axCBpval = plt.axes([.9, .01, .15, .1])
CBpval = CheckButtons(axCBpval, ['filter pval'], [stat_pvalfilter])
def onpick_xcorrBSFidx(event):
    global fg, ax, p1, p2, mylag, myMOCsel
    xcoord = find_nearest(xax, event.xdata); xcoordidx = np.where(xax == xcoord)[0][0]
    ycoord = find_nearest(yax, event.ydata); ycoordidx = np.where(yax == ycoord)[0][0]
    print xax.min(), xax.max(), yax.min(), yax.max()
    if ((xcoordidx<=xax.min())|(xcoordidx>=xax.max())) & ((ycoordidx<=yax.min())|(ycoordidx>=yax.max())): return() # check wether inside pcolor-plot
    idxmylag = np.where(lags==mylag)[0][0]
    print('redrawing MOC at {} / {}\t(idx: {} / {})\t(xcorridx = {})'.format(\
    str(xcoord), str(ycoord), str(xcoordidx), str(ycoordidx), mxcorrBSFidx[ycoordidx, xcoordidx, idxmylag]))
    myMOCsel = (MOC[:,ycoordidx, xcoordidx]).squeeze()
    redraw_curve()
def updateLAG(label):
    global fg, ax, p1, p2, mylag, myMOCsel
    mylag = int(label);
    print('shifting MOC to lag={}'.format(mylag))
    draw_pcolor();
    redraw_curve()
def filterpval(label):
    global stat_pvalfilter, cb
    stat_pvalfilter = not stat_pvalfilter
    draw_pcolor();

# update widgets
fg.canvas.mpl_connect('button_press_event', onpick_xcorrBSFidx)
RBlag.on_clicked(updateLAG)
CBpval.on_clicked(filterpval)


# -----------------------------------------------------------------------------
# moving min/max of BSFidx
fg, ax = spsinv(); 
p00 = ax.pcolor(lat_auxgrd, zT, mxcorrBSFidx[:,:,idxlag0])
p01 = ax.plot(fminx, fminy, 'o-', prop_of_min); 
p02 = ax.plot(fmaxx, fmaxy, 'o-', prop_of_max); 
ax.legend()

# -----------------------------------------------------------------------------
# movie of BSFidx over lags with max/min as dots
writer = utils_plt.movie_setup()
fg, ax = spsinv(); cb00 = 'dummy'
with writer.saving(fg, path_figscorr+'xcorrBSFidx_201.mp4', 100): # 100 dpi
    for idx, n in enumerate(lags):
        print n
        ax = utils_plt.movie_newframe(fg, ax, 111, cb00)
        p00 = ax.pcolor(lat_auxgrd, zT, mxcorrBSFidx[:,:,idx]); cb00 = cb(p00,ax); 
        p01 = ax.plot(fminsc[idx,0], fminsc[idx,1], 'o', **prop_of_min);
        p02 = ax.plot(fmaxsc[idx,0], fmaxsc[idx,1], 'o', **prop_of_max);
        yl(cb00.ax, 'xcorrBSFidx'); plt.text(60, 5500, 'lag = {}'.format(n), **prop_txt_on_floor);
        utils_plt.ptype_lat_z(fg, ax, lat_auxgrd, HT_auxgrd_xmax)
        writer.grab_frame()

# -----------------------------------------------------------------------------
# movie of BSF over lags with max/min as dots
times = mBSF.shape[0]
writer = utils_plt.movie_setup()
fg, ax = spsinv(); cb00 = 'dummy'
with writer.saving(fg, path_figscorr+'xcorrBSFidx_201.mp4', 100): # 100 dpi
    for idx, n in enumerate(times):
        ax = utils_plt.movie_newframe(fg, ax, 111, cb00)
        p00 = ax.pcolor(lat_auxgrd, zT, mBSF[:,:,idx]); cb00 = cb(p00,ax); 
        yl(cb00.ax, 'xcorrBSFidx'); plt.text(60, 5500, 'lag = {}'.format(n), **prop_txt_on_floor);
        utils_plt.ptype_lat_z(fg, ax, lat_auxgrd, HT_auxgrd_xmax)
        writer.grab_frame()

# -----------------------------------------------------------------------------
# interactive plot corrMOCidx
fg, ax1 = fgspinv(211)
p1 = plt.pcolor(mxcorrMOCidxSurfex20[:,:,0], picker=True); plt.colorbar(); #! round
#axRBlag = plt.axes([0,0,.1,.1])
#RBlag = RadioButtons(axRBlag, ('-10', '-5', '-3', '-2', '-1', '0', '1', '2', '3', '5', '10'), active=5) # active 5 --> 0 lag
ax2 = host_subplot(212, axes_class=AA.Axes)
p2 = ax2.plot(MOCidx, label='MOCidx'); yl(ax2,'MOCidx');
ax3 = ax2.twinx()
p3 = ax3.plot(BSF[:,0,0], label='BSF'); yl(ax3,'BSF');
xl(ax2,'Time')
plt.legend(loc='topright')
def onpick_xcorrMOCidx(event):
    xcoord = find_nearest(np.arange(BSF.shape[1]),event.xdata)
    ycoord = find_nearest(np.arange(BSF.shape[2]),event.ydata)
    print('redrawing BSF at '+str(xcoord)+' / '+str(ycoord))
    newBSFsel = BSF[:,xcoord, ycoord]
    p3[0].set_ydata(newBSFsel)
    p3[0].set_label(str(xcoord)+' / '+str(ycoord))
    ax3.set_ylim([min(newBSFsel), max(newBSFsel)])
    plt.draw()
    print('done!')
# update
fg.canvas.mpl_connect('button_press_event', onpick_xcorrMOCidx)
#RBlag.on_clicked(updateLAG)      # for window width
plt.show()


# moving min/max of MOCidx
fmaxx, fmaxy, fminx, fminy = np.zeros(len(lags)), np.zeros(len(lags)), np.zeros(len(lags)), np.zeros(len(lags))
for idx, n in enumerate(lags[:]):
    fminx[idx] = np.where(mxcorrMOCidx[:,:,idx] == mxcorrMOCidx[:,:,idx].min())[1]
    fminy[idx] = np.where(mxcorrMOCidx[:,:,idx] == mxcorrMOCidx[:,:,idx].min())[0]
    fmaxx[idx] = np.where(mxcorrMOCidx[:,:,idx] == mxcorrMOCidx[:,:,idx].max())[1]
    fmaxy[idx] = np.where(mxcorrMOCidx[:,:,idx] == mxcorrMOCidx[:,:,idx].max())[0]
fg(); pc(mxcorrMOCidx[:,:,0])
pt(fminx, fminy, 'go', ms=12, lw=4, label='min'); pt(fmaxx, fmaxy, 'ro', ms=12, lw=4, label='max'); lg();

# movie of MOCidx over lags in two dimensions
import matplotlib.animation as manimation
FFMpegWriter = manimation.writers['avconv']
writer = FFMpegWriter(fps=2)
moviefig = fg()
with writer.saving(moviefig, "xcorrMOCidx_LAG201.mp4", 100):
    for idx, n in enumerate(lags):
        print 'frame'+str(n)
        plt.clf()
        p = plt.pcolor(mxcorrMOCidx[:,:,idx])
        pt(fminx[idx], fminy[idx], 'go', ms=12, lw=4, label='min');
        pt(fmaxx[idx], fmaxy[idx], 'ro', ms=12, lw=4, label='max');
        plt.title('lag = '+str(n))
        writer.grab_frame()


fg, ax = plt.subplots()
im1 = ax.imshow(mxcorrBSFidx[:,:,0], picker=True)

def onpick4(event):
    artist = event.artist
    if isinstance(artist, AxesImage):
        im = artist
        A = im.get_array()
        print('onpick4 image', A.shape)
        
fg.canvas.mpl_connect('pick_event', onpick4)
plt.show()


# =============================================================================
#  Streamfunctions
# =============================================================================

lat_mgrd = ncdat.TLAT.isel(nlon=0)          # mean of LAT for each j #! very inappropriate

# -----------------------------------------------------------------------------
# BSF on geographical grid calculated by model
fg, map = utils_plt.plot_BSF(BSF.isel(time=0), 'T', nlevels = 10)
plt.title('BSF model on T grid')
 #utils_plt.print2pdf(fg, 'testfigures/BSF_model_T')
# -----------------------------------------------------------------------------
# MOC on geographical grid calculated by model
MOCsel = MOC.isel(time=0)
fg, ax = utils_plt.plot_MOC(MOCsel.lat_aux_grid, MOCsel.moc_z, MOCsel, nlevels=40, plttype='pcolor+contour')
plt.plot(lat_mgrd, HT_auxgrd_xmax) 				# plot seafloor
plt.title('MOC model')
plt.xlim([-36,90])
 #utils_plt.print2pdf(fg, path_figs+'MOC_model')

# =============================================================================
#  Streamfunction Indices
# =============================================================================

# -----------------------------------------------------------------------------
# BSF on geographical grid calculated by model
fg, map = utils_plt.plot_BSF(xcorrBSFidx, 'T', nlevels = 10)
plt.title('BSF model on T grid')
 #utils_plt.print2pdf(fg, 'testfigures/BSF_model_T')
# -----------------------------------------------------------------------------
# MOC on geographical grid calculated by model
fg, ax = utils_plt.plot_MOC(MOC_model.lat_aux_grid, MOC_model.moc_z, MOC_model, nlevels=40, plttype='pcolor+contour')
plt.plot(lat_auxgrd,HT_auxgrd_xmax) 				# plot seafloor
plt.title('MOC model')
plt.xlim([-36,90])
 #utils_plt.print2pdf(fg, 'testfigures/MOC_model')

# --------------
ncdat.close()
