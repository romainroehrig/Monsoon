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
import matplotlib as mpl#Choix de la variable

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
data_nc_mod = '/cnrm/cmip6/CMIP/CNRM-CERFACS/CNRM-CM6-1/amip/r1i1p1f2/Amon/'+str(var)+'/gr/files/d20180711/'+str(var)+'_Amon_CNRM-CM6-1_amip_r1i1p1f2_gr_197901-201412.nc'

file = nc.Dataset(data_nc_mod)

if var=='pr':
	conv=86400
else:
	conv=1
	
#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
time = file.variables['time'][:]
dates = nc.num2date(file['time'][:],units=(file['time']).units,calendar=(file['time']).calendar)
lon=(shiftgrid(180,file.variables[str(var)][:,:,:]*conv,nc.Dataset(data_nc_mod).variables['lon'][:],start=False,cyclic=360.0))[1]
data_mod = (shiftgrid(180,file.variables[str(var)][:,:,:]*conv,nc.Dataset(data_nc_mod).variables['lon'][:],start=False,cyclic=360.0))[0]
data_mod=np.average(data_mod[5:8:12],axis=0)
# Average between 10W and 10E
# In case longitude goes from 0 to 360, make them between -180 and 180.
if np.max(lon) > 180.:
    lon = np.where(lon > 180, lon-360, lon)

# select longitude indices between -10 and 10. All parentheses are important
lon_inds = np.where((lon > -10) & (lon < 10))[0]

# final average
data_mod = np.ma.average(data_mod[:,lon_inds],axis=1)

	
if var=='pr':
	unit='(mm/day)'
elif var=='prw':
	unit='mm'
elif var=='tas':
	unit='K'

L=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0601']
	L+=[str(1998+k)+'0701']
	L+=[str(1998+k)+'0801']	
	L+=[str(1998+k)+'0901']


data_nc_obs='/cnrm/amacs/USERS/roehrig/datasets/TRMM/3B43/Version7.regrid_T127_conserve/3B43.19980801.regrid_T127_conserve.nc'
file=nc.Dataset(data_nc_obs)

path='/cnrm/amacs/USERS/roehrig/datasets/TRMM/3B43/Version7.regrid_T127_conserve/3B43.'
fichier=[]
for i in L:
	fichier.append(path+str(i)+'.regrid_T127_conserve.nc')

data_obs=[]

for j in range(len(fichier)):
	data_obs.append((shiftgrid(180,nc.Dataset(fichier[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])

lat = file.variables['latitude'][:]
lon=(shiftgrid(180,file.variables[str(var)][:,:,:]*24,nc.Dataset(data_nc_obs).variables['longitude'][:],start=False,cyclic=360.0))[1]


# Average between 10W and 10E
# In case longitude goes from 0 to 360, make them between -180 and 180.
if np.max(lon) > 180.:
    lon = np.where(lon > 180, lon-360, lon)


# select longitude indices between -10 and 10. All parentheses are important
lon_inds = np.where((lon > -10) & (lon < 10))[0]

# final average
for i in range(len(fichier)):
	data_obs[i] = np.ma.average(data_obs[i][:,:,lon_inds],axis=2)

data_obs=np.average(data_obs[:],axis=0)

fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(12,10))
ax.set_xlim(-10,30)
ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
plt.plot(lat,data_mod[:],label=str(var)+'_de_Arpege-Climat')
plt.plot(lat,data_obs[0,:],label=str(var)+'_d une_re-analyse_de_GPCP')
plt.legend()
plt.xlabel('Latitude')
plt.ylabel(str(var)+'('+(str(unit))+')')
plt.title('y_obs_10W10E_'+str(var)+'JJAS')
plt.savefig('Coupe_compa_obs/'+str(var)+'/y_obs_10W10E_'+str(var)+'_JJAS.png')
plt.show()







