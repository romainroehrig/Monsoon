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
var='pr'
L=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0601']
	L+=[str(1998+k)+'0701']
	L+=[str(1998+k)+'0801']	
	L+=[str(1998+k)+'0901']


data_nc='/cnrm/amacs/USERS/roehrig/datasets/TRMM/3B43/Version7.regrid_T127_conserve/3B43.19980801.regrid_T127_conserve.nc'
file=nc.Dataset(data_nc)

path='/cnrm/amacs/USERS/roehrig/datasets/TRMM/3B43/Version7.regrid_T127_conserve/3B43.'
fichier=[]
for i in L:
	fichier.append(path+str(i)+'.regrid_T127_conserve.nc')

data=[]

for j in range(len(fichier)):
	data.append((shiftgrid(180,nc.Dataset(fichier[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


data=np.average(data[:],axis=0)

#Lecture des variables
lat = file.variables['latitude'][:]
lon = file.variables['longitude'][:]



#Mise en place de la carte
proj = ccrs.PlateCarree(central_longitude=0, globe=None)
cmap = cm.hsv
unit = '(mm/day)'
bounds =[ i for i in np.arange(0,25,1)]
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist[0] = (1,1,1,0)
ext = 'max'


cmap = cmap.from_list('Custom cmap', cmaplist[:-20], len(cmaplist)-20)
norm = mpl.colors.BoundaryNorm(bounds,cmap.N)



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

cs=ax.contourf((shiftgrid(180,file.variables[str(var)][:,:,:]*24,nc.Dataset(data_nc).variables['longitude'][:],start=False,cyclic=360.0))[1],lat, data[0,:,:], bounds, cmap=cmap, norm=norm, extend=str(ext), transform=ccrs.PlateCarree())
cbar = fig.colorbar(cs, shrink=0.9, orientation='horizontal', pad=0.1, extend='max', boundaries=bounds, ticks=bounds , label=str(unit))
plt.title('xy_pr_JJAS')
plt.savefig('obs/xy_pr_JJAS.png')
plt.show()

