# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 12:07:49 2020

@author: Ludovic
"""

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
from data import getfile

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

var='tas'
cont='WA'
sim='WA_10a'
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

if var=='pr':
	conv=86400
elif var=='hur':
	conv=100
else:
    conv=1
	
#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
time = file.variables[str(t)][:]
dates = nc.num2date(file[str(t)][:],units=(file[str(t)]).units,calendar=(file[str(t)]).calendar)
data = file.variables[var][:,:,:]*conv
lon=(shiftgrid(180,file.variables[str(var)][:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1]
data = (shiftgrid(180,file.variables[str(var)][:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]

if cont=='All':
    lat1=-90
    lat2=90
    lon1=-180
    lon2=180
elif cont=='WA':
    lat1=5
    lat2=35
    lon1=0
    lon2=60

# select longitude indices between -10 and 10. All parentheses are important
lat_inds = np.where((lat > lat1) & (lat < lat2))[0]

# final average
data = np.ma.average(data[:,lat_inds,:],axis=1)

# select longitude indices between -10 and 10. All parentheses are important
lon_inds = np.where((lon > lon1) & (lon < lon2))[0]

# final average
data = np.ma.average(data[:,lon_inds],axis=1)



data_new=np.zeros(120)	

for i in (range(120)):
    data_new[i]=(sum(data[0:(i+1)])/(i+2))
print(data_new[0])


if var=='pr':
	unit='mm/day'
elif var=='prw':
	unit='mm'
elif var=='ts':
	unit='K'
elif var=='hur':
    unit='%'
elif var=='ua':
    unit='m/s'
elif var=='va':
    unit='m/s'

fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(12,10))


plt.plot(data_new)
plt.xlabel('Month')
plt.xticks(range(0,132,12))
plt.ylabel(str(var)+'('+(str(unit))+')')
plt.title('spinup_'+str(var)+str(cont)+'_'+str(sim))
plt.savefig('Spinup/spinup_'+str(var)+str(cont)+'_'+str(sim)+'.png')
plt.show()
