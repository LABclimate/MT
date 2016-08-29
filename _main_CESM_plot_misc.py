'''
CESM controls --> plotting routines

@author: buerki@climate.unibe.ch
TODO: 	 add '.values' where possible to speed up code.

'''
import matplotlib.pyplot as plt
import matplotlib as ml
import CESM_utils_plt as utils_plt

plt.ion() # enable interactive mode
path_fig = '../figures/160711/'


# =======================================================================================
#  BSF / MOC
# =======================================================================================
var = MV[0,:,:]
var[190:300,-65:-45] = np.zeros_like(var[190:300,-65:-45])
fig, map = utils_plt.plot_BSF(var, 'T', nlevels = 10)

# -----------------------------------------------------------------------------------------
# BSF on geographical grid calculated by model
BSF_model = utils_mask.mask_ATLANTIC(ncdat.BSF.isel(time=0), ncdat.REGION_MASK)
fig, map = utils_plt.plot_BSF(BSF_model, 'T', nlevels = 10)
plt.title('BSF model on T grid')
utils_plt.print2pdf(fig, path_fig+'BSF_model_T')
# -----------------------------------------------------------------------------------------
# BSF on model grid
fig, map = utils_plt.plot_BSF(BSF_mmgrd, 'U', nlevels=100)
plt.title('BSF mgrd on U grid')
utils_plt.print2pdf(fig, 'testfigures/BSF_mgrd_U')
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

fig, ax = utils_plt.plot_MOC(ncdat.VVEL.TLAT[:,0], z_w_auxgrd, MV_projauxgrd[:,:,59])
fig, ax = utils_plt.plot_MOC(ncdat.VVEL.TLAT[:,0], z_w_auxgrd, ncdat.VVEL[0,:,:,59])
