# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 15:51:35 2020

@author: Ludovic
"""

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


#Variables

var1='ta'
var2='hus'
var3='zg'
var4='ua'
var5='va'
var6='wap'
sim='Simu_Base'
cont='real'	
seasons = {'year': range(1,13),
           'JJAS': [6,7,8,9],
           }
if sim=='Simu_Base':
    t='time'
else:
    t='time_counter'
    
#Lecture des fichiers
    
data_nc1 = getfile(sim,var1)
file1 = nc.Dataset(data_nc1)
data_nc2 = getfile(sim,var2)
file2 = nc.Dataset(data_nc2)
data_nc3 = getfile(sim,var3)
file3 = nc.Dataset(data_nc3)
data_nc4 = getfile(sim,var4)
file4 = nc.Dataset(data_nc4)
data_nc5 = getfile(sim,var5)
file5 = nc.Dataset(data_nc5)
data_nc6 = getfile(sim,var6)
file6 = nc.Dataset(data_nc6)

#Lecture des variables

conv=86400
lat = file1.variables['lat'][:]
lon = file1.variables['lon'][:]
lv=file1.variables['plev'][:]
#Lecture MSE
dates = nc.num2date(file1[str(t)][:],units=(file1[str(t)]).units,calendar=(file1[str(t)]).calendar)
data = (shiftgrid(180,4.186*file1.variables[str(var1)][:,:,:,:]+2.257*file2.variables[str(var2)][:,:,:,:]+9.8*file3.variables[str(var3)][:,:,:,:],nc.Dataset(data_nc1).variables['lon'][:],start=False,cyclic=360.0))[0]
inds = [d.month in seasons['JJAS'] for d in dates]
data = np.ma.average(data[inds],axis=0)
#Lecture variable de vent
data_u = (shiftgrid(180,file4.variables[str(var4)][:,:,:,:],nc.Dataset(data_nc1).variables['lon'][:],start=False,cyclic=360.0))[0]
data_v = (shiftgrid(180,file5.variables[str(var5)][:,:,:,:],nc.Dataset(data_nc1).variables['lon'][:],start=False,cyclic=360.0))[0]
data_w = (shiftgrid(180,file6.variables[str(var6)][:,:,:,:],nc.Dataset(data_nc1).variables['lon'][:],start=False,cyclic=360.0))[0]
data_u = np.ma.average(data_u[inds],axis=0)
data_v = np.ma.average(data_v[inds],axis=0)
data_w = np.ma.average(data_w[inds],axis=0)

#Calcul des termes d'advection

#Advection horizontale
data_p_u=np.zeros((19,128,256))
data_p_v=np.zeros((19,128,256))
data_advhor=np.zeros((19,128,256))
#Advection verticale
data_p_w=np.zeros((19,128,256))
data_adv=np.zeros((19,128,256))
for i in range(1,len(lv)-1):
    for j in range(1,len(lat)-1):
        for k in range(1,len(lon)-1):
            data_p_v[i,j,k]= np.ma.max(data_v[i,j,k],0)*((data[i,j,k]-data[i,j-1,k])/(lat[j]-lat[j-1])) + np.ma.min(data_v[i,j,k],0)*((data[i,j+1,k]-data[i,j,k])/(lat[j+1]-lat[j]))
            data_p_u[i,j,k]= np.ma.max(data_u[i,j,k],0)*((data[i,j,k]-data[i,j,k-1])/(lon[k]-lon[k-1])) + np.ma.min(data_u[i,j,k],0)*((data[i,j,k+1]-data[i,j,k])/(lon[k+1]-lon[k]))
            data_p_w[i,j,k]= np.ma.max(data_w[i,j,k],0)*((data[i,j,k]-data[i-1,j,k])/(lv[i]-lv[i-1])) + np.ma.min(data_w[i,j,k],0)*((data[i+1,j,k]-data[i,j,k])/(lv[i+1]-lv[i]))
            data_advhor[i,j,k]=data_p_u[i,j,k]+data_p_v[i,j,k]
            data_adv[i,j,k]=data_p_u[i,j,k]+data_p_v[i,j,k]+data_p_w[i,j,k]


#Intégration
data_int_uv=np.zeros((128,256))
data_int_w=np.zeros((128,256))
data_int=np.zeros((128,256))

for j in range(1,len(lat)-1):
    for k in range(1,len(lon)-1):
        data_int_uv[j,k]=((lv[18]-lv[0])/19)*((np.sum(data_advhor[0:len(lv)-1,j,k],axis=0)+np.sum(data_advhor[1:len(lv),j,k],axis=0))/2)
        data_int_w[j,k]=((lv[18]-lv[0])/19)*((np.sum(data_p_w[0:len(lv)-1,j,k],axis=0)+np.sum(data_p_w[1:len(lv),j,k],axis=0))/2)
        data_int[j,k]=((lv[18]-lv[0])/19)*((np.sum(data_adv[0:len(lv)-1,j,k],axis=0)+np.sum(data_adv[1:len(lv),j,k],axis=0))/2)


#Mise en place de la carte 
proj = ccrs.PlateCarree(central_longitude=0, globe=None)
cmap1 = cm.hsv
unit = 'J/s'
bounds1 =[ i for i in np.arange(-10000,10000,1000)]
cmaplist = [cmap1(i) for i in range(cmap1.N)]
level='_'
ext = 'max'
cmap1 = cmap1.from_list('Custom cmap', cmaplist[:-20], len(cmaplist)-20)
norm = mpl.colors.BoundaryNorm(bounds1,cmap1.N)

fig, ax = plt.subplots(figsize=(12,10),subplot_kw=dict(projection=proj))

if sim=='LC_Sym':
    lon1=-20
    lon2=45
    lat1=-50
    lat2=50
else:
    lon1=-20
    lon2=40
    lat1=-10
    lat2=50
"""   
ax.set_xlim(lon1,lon2)
ax.set_ylim(lat1,lat2)
"""
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



cs= ax.contourf((shiftgrid(180,file1.variables[str(var1)][:,:,:,:],nc.Dataset(data_nc1).variables['lon'][:],start=False,cyclic=360.0))[1] ,lat, data_int_uv[:,:], bounds1, cmap=cmap1, norm=norm, extend=str(ext), transform=ccrs.PlateCarree())
cbar = fig.colorbar(cs, shrink=0.9, orientation='horizontal', pad=0.1, extend='max', boundaries=bounds1, ticks=bounds1)

plt.savefig('MSE.png')
