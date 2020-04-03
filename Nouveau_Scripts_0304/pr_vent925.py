# -*- coding: utf-8 -*-


import os

import netCDF4 as nc

import numpy as np
import numpy.ma as ma

import matplotlib as mpl

import matplotlib.pyplot as plt
from matplotlib import cm

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from myfunctions import shiftgrid
from data_atlas import getfile
from seasons import seasons


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

var='pr'
sim='LC_Sym'
cont='ideal'	
seasons = {'year': range(1,13),
           'JJAS': [6,7,8,9],
           }
if sim=='LC_Sym':
    t='time'
else:
    t='time_counter'
#Lecture du fichier
data_nc = getfile(sim,var)


file = nc.Dataset(data_nc)


#Lecture des variables
conv=86400
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
#time = file.variables['time'][:]
dates = nc.num2date(file[str(t)][:],units=(file[str(t)]).units,calendar=(file[str(t)]).calendar)
data = file.variables[str(var)][:,:,:]*conv
data = (shiftgrid(180,file.variables[str(var)][:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
inds = [d.month in seasons['JJAS'] for d in dates]
data_p=np.average(data[inds],axis=0)


#Mise en place 
proj = ccrs.PlateCarree(central_longitude=0, globe=None)
cmap1 = cm.hsv
unit = 'mm'
bounds1 =[ i for i in np.arange(0,14,1)]
cmaplist = [cmap1(i) for i in range(cmap1.N)]
cmaplist[0] = (1,1,1,0)
level='_'
data_out=data[:,:]
ext = 'max'
cmap1 = cmap1.from_list('Custom cmap', cmaplist[:-20], len(cmaplist)-20)
norm = mpl.colors.BoundaryNorm(bounds1,cmap1.N)



var='ua'
data_nc=getfile(sim,var)
file = nc.Dataset(data_nc)

#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
#time = file.variables['time'][:]
unit='(m/s)'

conv=1
ext='both'
data = (shiftgrid(180,file.variables[str(var)][:,:,:,:],nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
lon=(shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1]
inds = [d.month in seasons['JJAS'] for d in dates]
U=np.average(data[inds],axis=0)
U=U[1,30:110:2,110:190:2]

var='va'
data_nc=getfile(sim,var)
file = nc.Dataset(data_nc)
bounds =[ i for i in np.arange(-5,6,1)]
unit='(m/s)'
cmap = cm.hsv
cmaplist = [cmap(i) for i in range(cmap.N)]
conv=1
ext='both'

#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
#time = file.variables['time'][:]
unit='(m/s)'
conv=1
ext='max'
data = (shiftgrid(180,file.variables[str(var)][:,:,:,:],nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
lon=(shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1]
V=np.average(data[inds],axis=0)
V=V[1,30:110:2,110:190:2]


#Plot
fig, ax = plt.subplots(figsize=(12,10),subplot_kw=dict(projection=proj))

if sim=='LC_Sym':
    lon1=-20
    lon2=65
    lat1=-50
    lat2=50
else:
    lon1=-20
    lon2=80
    lat1=-10
    lat2=50
    
ax.set_xlim(lon1,lon2)
ax.set_ylim(lat1,lat2)

if sim=='LC_Sym':
    lon1=0
    lon2=45
    lat1=-35
    lat2=35
else:
    lon1=0
    lon2=60
    lat1=5
    lat2=35
    
xticks = range(-180,181,10)
yticks = range(-90,91,10)

xxc = [lon1,lon2,lon2,lon1,lon1]
yyc = [lat1,lat1,lat2,lat2,lat1]

if cont == 'real':
    ax.coastlines('50m')
    ax.add_feature(cf.BORDERS)
elif cont == 'ideal':
    ax.plot(xxc,yyc,color='k',linewidth=2,transform=ccrs.PlateCarree())
gl = ax.gridlines(xlocs=xticks, ylocs=yticks,linestyle='--',lw=0.5,color='dimgrey',draw_labels=True)
gl.xlabels_top = False
gl.xformatter = LONGITUDE_FORMATTER
#gl.xlines = False
gl.ylabels_right = False
gl.yformatter = LATITUDE_FORMATTER
#gl.ylines = False





cs=ax.contourf((shiftgrid(180,file.variables[str(var)][:,:,:]*24,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1],lat, data_p[:,:], bounds1, cmap=cmap1, norm=norm, extend=str(ext), transform=ccrs.PlateCarree())
q = ax.quiver(lon[110:190:2],lat[30:110:2],U, V,scale=None, scale_units='inches')
q._init()
assert isinstance(q.scale, float)
ax.quiver(lon[110:190:2],lat[30:110:2],U*1, V*1, scale=q.scale, scale_units='inches')
cbar = fig.colorbar(cs, shrink=0.9, orientation='horizontal', pad=0.1, extend='max', boundaries=bounds1, ticks=bounds1 , label='(mm/day)')
plt.title('vecteur_vent_'+str(sim))
"""
plt.savefig('images/vect_vent/vecteur_vent_'+str(sim)+'.png')
"""
plt.savefig('essaie.png')
