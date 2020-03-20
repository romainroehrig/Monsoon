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

var='pr'	

#Lecture du fichier
data_nc = 'C:/Users/Ludovic/Desktop/PFE/simu_WA/DELANNOY_slab_cont_WA_arpsfx_monthly_'+str(var)+'_1850-1855.nc'


file = nc.Dataset(data_nc)


#Lecture des variables
conv=86400
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
#time = file.variables['time'][:]
#dates = nc.num2date(file['time'][:],units=(file['time']).units,calendar=(file['time']).calendar)
data = file.variables[str(var)][:,:,:]*conv
data = (shiftgrid(180,file.variables[str(var)][:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
data_p=np.average(data[5:8:12],axis=0)


#Mise en place 
proj = ccrs.PlateCarree(central_longitude=0, globe=None)
cmap1 = cm.hsv
unit = 'mm'
bounds1 =[ i for i in np.arange(0,25,1)]
cmaplist = [cmap1(i) for i in range(cmap1.N)]
cmaplist[0] = (1,1,1,0)
level='_'
data_out=data[:,:]
ext = 'max'
cmap1 = cmap1.from_list('Custom cmap', cmaplist[:-20], len(cmaplist)-20)
norm = mpl.colors.BoundaryNorm(bounds1,cmap1.N)



var='ua'
data_nc='C:/Users/Ludovic/Desktop/PFE/Simu_WA/DELANNOY_slab_cont_WA_arpsfx_monthly_'+str(var)+'_1850-1855.nc'
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
U=np.average(data[5:8:12],axis=0)
U=U[1,50:90,120:175]

var='va'
data_nc='C:/Users/Ludovic/Desktop/PFE/Simu_WA/DELANNOY_slab_cont_WA_arpsfx_monthly_'+str(var)+'_1850-1855.nc'
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
ext='both'
data = (shiftgrid(180,file.variables[str(var)][:,:,:,:],nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
lon=(shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1]
V=np.average(data[5:8:12],axis=0)
V=V[1,50:90,120:175]


#Plot
fig, ax = plt.subplots(figsize=(12,10),subplot_kw=dict(projection=proj))

"""
ax.coastlines('50m')
ax.add_feature(cf.BORDERS)
"""
ax.set_xlim(-20,70)
ax.set_ylim(-35,45)

xticks = range(0,65,5)
yticks = range(5,40,5)

gl = ax.gridlines(xlocs=xticks, ylocs=yticks,linestyle='--',lw=1,color='dimgrey',draw_labels=True)
gl.xlabels_top = False
gl.xformatter = LONGITUDE_FORMATTER
#gl.xlines = False
gl.ylabels_right = False
gl.yformatter = LATITUDE_FORMATTER
#gl.ylines = False

"""
xticks = range(-180,181,10)
yticks = range(-90,91,10)
"""
gl = ax.gridlines(xlocs=xticks, ylocs=yticks,linestyle='--',linewidth=3,color='black')
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER


cs=ax.contourf((shiftgrid(180,file.variables[str(var)][:,:,:]*24,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1],lat, data_p[:,:], bounds1, cmap=cmap1, norm=norm, extend=str(ext), transform=ccrs.PlateCarree())
q = ax.quiver(lon[120:175],lat[50:90],U, V,scale=None, scale_units='inches')
q._init()
assert isinstance(q.scale, float)
ax.quiver(lon[120:175],lat[50:90],U*1.5, V*1.5, scale=q.scale, scale_units='inches')
cbar = fig.colorbar(cs, shrink=0.9, orientation='horizontal', pad=0.1, extend='max', boundaries=bounds1, ticks=bounds1 , label='(mm/day)')
plt.title('pr_vent925')
plt.savefig('pr_vent925.png')


