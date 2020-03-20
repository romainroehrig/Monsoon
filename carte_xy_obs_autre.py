# -*- coding: utf-8 -*-

import netCDF4 as nc
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy, cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib as mp
import numpy.ma as ma
from matplotlib import cm
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib as mpl

#Correction probleme àlon0

def shiftgrid(lon0,datain,lonsin,start=False,cyclic=360.0):
    if np.fabs(lonsin[-1]-lonsin[0]-cyclic) > 1.e-4:
        # Use all data instead of raise ValueError, 'cyclic point not included'
        start_idx = 0
    else:
        # If cyclic, remove the duplicate point
        start_idx = 1
    if lon0 < lonsin[0] or lon0 > lonsin[-1]:
        raise ValueError('lon0 outside of range of lonsin')
    i0 = np.argmin(np.fabs(lonsin-lon0))
    i0_shift = len(lonsin)-i0
    if ma.isMA(datain):
        dataout  = ma.zeros(datain.shape,datain.dtype)
    else:
        dataout  = np.zeros(datain.shape,datain.dtype)
    if ma.isMA(lonsin):
        lonsout = ma.zeros(lonsin.shape,lonsin.dtype)
    else:
        lonsout = np.zeros(lonsin.shape,lonsin.dtype)
    if start:
        lonsout[0:i0_shift] = lonsin[i0:]
    else:
        lonsout[0:i0_shift] = lonsin[i0:]-cyclic
    dataout[...,0:i0_shift] = datain[...,i0:]
    if start:
        lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]+cyclic
    else:
        lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]
    dataout[...,i0_shift:] = datain[...,start_idx:i0+start_idx]
    return dataout,lonsout	
#Choix de la variable
var='t'
level='925'

#Lecture du fichier
data_nc = '/cnrm/amacs/USERS/roehrig/datasets/ERAI/monthly_0.75x0.75/netcdf/regrid_arptl127l31r_bilinear/'+str(var)+'.regrid_ARPEGE_tl127l31r_bilinear.levARP.nc'
file = nc.Dataset(data_nc)


if var=='u':
	bounds =[ i for i in np.arange(-15,16,1)]
	unit='(m/s)'
	cmap = cm.seismic
	cmaplist = [cmap(i) for i in range(cmap.N)]
	conv=1
	ext='both'
elif var=='v':
	bounds =[ i for i in np.arange(-5,6,1)]
	unit='(m/s)'
	cmap = cm.seismic
	cmaplist = [cmap(i) for i in range(cmap.N)]
	conv=1
	ext='both'
elif var=='w':
	bounds =[ i for i in np.arange(-0.15,0.16,0.02)]
	unit='(Pa/s)'
	cmapref = cm.seismic
	conv=1
	ext='both'
elif var=='r':
	bounds =[ i for i in np.arange(0,104,4)]
	unit='%'
	cmap = cm.RdBu
	cmaplist = [cmap(i) for i in range(cmap.N)]
	conv=1
	ext='neither'
elif var=='t':
	unit = 'K'
	bounds =[ i for i in np.arange(290,316,1)]
	cmap = cm.coolwarm
	cmaplist = [cmap(i) for i in range(cmap.N)]
	ext= 'both'
	level='_'
	conv=1



#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
time = file.variables['time'][:]
lev= file.variables['plev'][:]
#dates = nc.num2date(file['time'][:],units=(file['time']).units,calendar=(file['time']).calendar)
data = (shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
lon=(shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1]
data=np.average(data[5:8:12],axis=0)

#Mise en place de la carte
proj = ccrs.PlateCarree(central_longitude=0, globe=None)
cmap = cmap.from_list('Custom cmap', cmaplist[:-20], len(cmaplist)-20)
cmaplist = [cmap(i) for i in range(cmap.N)]
norm = mpl.colors.BoundaryNorm(bounds,cmap.N)
#cmaplist[0] = (1,1,1,0)







#Plot
fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(12,10),subplot_kw=dict(projection=proj))

ax.coastlines('50m')
ax.add_feature(cf.BORDERS)

ax.set_xlim(-25,35)
ax.set_ylim(-5,37)

xticks = range(-180,181,10)
yticks = range(-90,91,10)
gl = ax.gridlines(xlocs=xticks, ylocs=yticks,linestyle='--',lw=1,color='dimgrey',draw_labels=True)
gl.xlabels_top = False
gl.xformatter = LONGITUDE_FORMATTER
#gl.xlines = False
gl.ylabels_right = False
gl.yformatter = LATITUDE_FORMATTER
#gl.ylines = False

cs=ax.contourf((shiftgrid(180,file.variables[str(var)][:,:,:,:],nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1],lat, data[0,:,:], bounds, cmap=cmap, norm=norm, extend=str(ext), transform=ccrs.PlateCarree())
cbar = fig.colorbar(cs, shrink=0.9, orientation='horizontal', pad=0.1, extend=str(ext), boundaries=bounds, ticks=bounds , label=str(unit))
plt.title('xy_'+str(var)+'_'+str(level)+'_JJAS')
plt.savefig('obs/xy_'+str(var)+'_'+str(level)+'_JJAS.png')
plt.show()

