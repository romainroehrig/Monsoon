# -*- coding: utf-8 -*-

import netCDF4 as nc
import numpy as np
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
var='ua'
level='925'

#Lecture du fichier
data_nc = 'C:/Users/Ludovic/Desktop/PFE/Amon_SST/'+str(var)+'/gr/files/d20180711/'+str(var)+'_Amon_CNRM-CM6-1_amip_r1i1p1f2_gr_197901-201412.nc'

file = nc.Dataset(data_nc)

if var=='pr':
	conv=86400
else:
	conv=1
#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
time = file.variables['time'][:]
dates = nc.num2date(file['time'][:],units=(file['time']).units,calendar=(file['time']).calendar)
data = file.variables[str(var)][:,:,:]*conv
data = (shiftgrid(180,file.variables[str(var)][:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
data=np.average(data[5:8:12],axis=0)


#Mise en place de la carte
proj = ccrs.PlateCarree(central_longitude=0, globe=None)
if var=='pr' :
	cmap = cm.hsv
	unit = 'mm'
	bounds =[ i for i in np.arange(0,25,1)]
	cmaplist = [cmap(i) for i in range(cmap.N)]
	cmaplist[0] = (1,1,1,0)
	level='_'
	data_out=data[:,:]
	ext = 'max'
elif var=='tas':
	cmap = cm.coolwarm
	unit = 'K'
	bounds =[ i for i in np.arange(289,313,2)]
	cmaplist = [cmap(i) for i in range(cmap.N)]
	level='_'
	data_out=data[:,:]
	ext= 'both'
	
elif var=='ua':
	cmap = cm.seismic
	unit ='m/s'
	if level=='925':
		bounds =[ i for i in np.arange(-15,16,1)]
		k=1
	elif level=='700':
		bounds =[ i for i in np.arange(-20,22,2)]
		k=3
	elif level=='200':
		bounds =[ i for i in np.arange(-15,16,1)]
		k=7
	elif level=='600':
		bounds =[ i for i in np.arange(-20,22,2)]
		k=4
	cmaplist = [cmap(i) for i in range(cmap.N)]
	data_out=data[k,:,:]
	ext='both'
elif var=='va':
	cmap = cm.seismic
	unit ='m/s'
	if level=='925':
		bounds =[ i for i in np.arange(-20,21,1)]
		k=1
	elif level=='700':
		bounds =[ i for i in np.arange(-15,16,1)]
		k=3
	elif level=='200':
		bounds =[ i for i in np.arange(-20,22,2)]
		k=9
	cmaplist = [cmap(i) for i in range(cmap.N)]
	data_out=data[k,:,:]
	ext='both'
elif var=='prw':
	cmap =cm.Blues
	unit='mm'
	bounds =[ i for i in np.arange(0,60,2)]
	cmaplist = [cmap(i) for i in range(cmap.N)]
	level='_'
	data_out=data[:,:]
	ext='max'
elif var=='clt':
	cmap=cm.Greys
	unit='%'
	bounds =[ i for i in np.arange(0,105,5)]
	cmaplist = [cmap(i) for i in range(cmap.N)]
	level='_'
	data_out=data[:,:]
	ext='neither'

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
cs = ax.contourf((shiftgrid(180,file.variables[str(var)][:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1],lat, data_out, bounds, cmap=cmap, norm=norm, extend=str(ext), transform=ccrs.PlateCarree())
cbar = fig.colorbar(cs, shrink=0.9, orientation='horizontal', pad=0.1, extend='max', boundaries=bounds, ticks=bounds , label=str(unit))
plt.title('xy_'+str(var)+'_'+str(level)+'_JJAS')
plt.savefig('Carte_de_moyenne/'+str(var)+'/xy_'+str(var)+'_'+str(level)+'_JJAS.png')
plt.show()






