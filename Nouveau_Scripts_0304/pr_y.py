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
import matplotlib as mpl#Choix de la variable
from data_atlas import getfile

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
seasons = {'year': range(1,13),
           'JJAS': [6,7,8,9],
           }
if sim=='LC_Sym':
    t='time'
else:
    t='time_counter'
if sim=='LC_Sym':
    lv='plev'
else:
    lv='pstd'
#Lecture du fichier
data_nc=getfile(sim,var)

file = nc.Dataset(data_nc)
print(file)
conv=86400

#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
#time = file.variables['time'][:]
dates = nc.num2date(file[str(t)][:],units=(file[str(t)]).units,calendar=(file[str(t)]).calendar)
data = file.variables[var][:,:,:]*conv
lon=(shiftgrid(180,file.variables[str(var)][:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1]
data = (shiftgrid(180,file.variables[str(var)][:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
inds = [d.month in seasons['JJAS'] for d in dates]
data_p=np.average(data[inds],axis=0)
# Average between 10W and 10E

# In case longitude goes from 0 to 360, make them between -180 and 180.
if np.max(lon) > 180.:
    lon = np.where(lon > 180, lon-360, lon)

# select longitude indices between -10 and 10. All parentheses are important
lon_inds = np.where((lon > 0) & (lon < 45))[0]

# final average
data = np.ma.average(data_p[:,lon_inds],axis=1)

unit='mm/day'

fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(12,10))
ax.set_xlim(-35,35)
plt.plot(lat,data)
plt.xlabel('Latitude')
plt.ylabel(str(var)+'('+(str(unit))+')')
plt.title('pr_y_030'+str(sim))
plt.savefig('images/pr_y/pr_y_030'+str(sim)+'.png')


