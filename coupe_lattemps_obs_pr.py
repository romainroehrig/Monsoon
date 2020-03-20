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
import matplotlib as mpl #Choix de la variable

firstwhite = True 
extend = 'max'
bounds =[ i for i in np.arange(0,10.5,0.5)]

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
conv=24

data_nc='/cnrm/amacs/USERS/roehrig/datasets/TRMM/3B43/Version7.regrid_T127_conserve/3B43.19980801.regrid_T127_conserve.nc'
file=nc.Dataset(data_nc)

path='/cnrm/amacs/USERS/roehrig/datasets/TRMM/3B43/Version7.regrid_T127_conserve/3B43.'

L=[]
fichier_janv=[]
data_janv=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0101']
for i in L:
	fichier_janv.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_janv)):
	data_janv.append((shiftgrid(180,nc.Dataset(fichier_janv[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_janv[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_fev=[]
data_fev=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0201']
for i in L:
	fichier_fev.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_fev)):
	data_fev.append((shiftgrid(180,nc.Dataset(fichier_fev[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_fev[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_mars=[]
data_mars=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0301']
for i in L:
	fichier_mars.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_mars)):
	data_mars.append((shiftgrid(180,nc.Dataset(fichier_mars[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_mars[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_avr=[]
data_avr=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0401']
for i in L:
	fichier_avr.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_avr)):
	data_avr.append((shiftgrid(180,nc.Dataset(fichier_avr[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_avr[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_mai=[]
data_mai=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0501']
for i in L:
	fichier_mai.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_mai)):
	data_mai.append((shiftgrid(180,nc.Dataset(fichier_mai[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_mai[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_juin=[]
data_juin=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0601']
for i in L:
	fichier_juin.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_juin)):
	data_juin.append((shiftgrid(180,nc.Dataset(fichier_janv[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_juin[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_juil=[]
data_juil=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0701']
for i in L:
	fichier_juil.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_juil)):
	data_juil.append((shiftgrid(180,nc.Dataset(fichier_juil[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_juil[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_aou=[]
data_aou=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0801']
for i in L:
	fichier_aou.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_aou)):
	data_aou.append((shiftgrid(180,nc.Dataset(fichier_aou[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_aou[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_sept=[]
data_sept=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'0901']
for i in L:
	fichier_sept.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_sept)):
	data_sept.append((shiftgrid(180,nc.Dataset(fichier_sept[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_sept[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_oct=[]
data_oct=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'1001']
for i in L:
	fichier_oct.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_oct)):
	data_oct.append((shiftgrid(180,nc.Dataset(fichier_oct[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_oct[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_nov=[]
data_nov=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'1101']
for i in L:
	fichier_nov.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_nov)):
	data_nov.append((shiftgrid(180,nc.Dataset(fichier_nov[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_nov[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])


L=[]
fichier_dec=[]
data_dec=[]
for k in range(0,21,1):
	L+=[str(1998+k)+'1201']
for i in L:
	fichier_dec.append(path+str(i)+'.regrid_T127_conserve.nc')	
for j in range(len(fichier_dec)):
	data_dec.append((shiftgrid(180,nc.Dataset(fichier_dec[j]).variables[str(var)][:,:,:]*24,nc.Dataset(fichier_dec[j]).variables['longitude'][:],start=False,cyclic=360.0))[0])

lon=(shiftgrid(180,file.variables[str(var)][:,:]*conv,nc.Dataset(data_nc).variables['longitude'][:],start=False,cyclic=360.0))[1]
# Average between 10W and 10E
if np.max(lon) > 180.:
    lon = np.where(lon > 180, lon-360, lon)

# select longitude indices between -10 and 10. All parentheses are important
lon_inds = np.where((lon > -10) & (lon < 10))[0]


	
# final average





for i in range(len(data_janv)):
	data_janv[i] = np.ma.average(data_janv[i][:,:,lon_inds],axis=2)
	data_fev[i] = np.ma.average(data_fev[i][:,:,lon_inds],axis=2)
	data_mars[i] = np.ma.average(data_mars[i][:,:,lon_inds],axis=2)
	data_avr[i]= np.ma.average(data_avr[i][:,:,lon_inds],axis=2)
	data_mai[i] = np.ma.average(data_mai[i][:,:,lon_inds],axis=2)
	data_juin[i] = np.ma.average(data_juin[i][:,:,lon_inds],axis=2)
	data_juil[i] = np.ma.average(data_juil[i][:,:,lon_inds],axis=2)
	data_aou[i] = np.ma.average(data_aou[i][:,:,lon_inds],axis=2)
	data_sept[i] = np.ma.average(data_sept[i][:,:,lon_inds],axis=2)
	data_oct[i] = np.ma.average(data_oct[i][:,:,lon_inds],axis=2)
	data_nov[i] = np.ma.average(data_nov[i][:,:,lon_inds],axis=2)
	data_dec[i] = np.ma.average(data_dec[i][:,:,lon_inds],axis=2)

data_janv = np.ma.average(data_janv[:],axis=0)
data_fev = np.ma.average(data_fev[:],axis=0)
data_mars = np.ma.average(data_mars[:],axis=0)
data_avr = np.ma.average(data_avr[:],axis=0)
data_mai = np.ma.average(data_mai[:],axis=0)
data_juin = np.ma.average(data_juin[:],axis=0)
data_juil = np.ma.average(data_juil[:],axis=0)
data_aou = np.ma.average(data_aou[:],axis=0)
data_sept = np.ma.average(data_sept[:],axis=0)
data_oct = np.ma.average(data_oct[:],axis=0)
data_nov = np.ma.average(data_nov[:],axis=0)
data_dec = np.ma.average(data_dec[:],axis=0)

data_new=np.zeros((12,128))
data_new[0,:]=data_janv[0]
data_new[1,:]=data_fev[0]
data_new[2,:]=data_mars[0]
data_new[3,:]=data_avr[0]
data_new[4,:]=data_mai[0]
data_new[5,:]=data_juin[0]
data_new[6,:]=data_juil[0]
data_new[7,:]=data_aou[0]
data_new[8,:]=data_sept[0]
data_new[9,:]=data_oct[0]
data_new[10,:]=data_nov[0]
data_new[11,:]=data_dec[0]




# customize colormap

cmap = cm.hsv
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist[0] = (1,1,1,0)	 
if extend == 'max':
    # remove the last 20 colors, to distinguish the 'over' color from the main colorbar
    	cmap = cmap.from_list('Custom cmap', cmaplist[:-20], len(cmaplist)-20)
	cmap.set_over(cmaplist[-1])
elif extend == 'min':
    print 'Not yet coded'
    sys.exit()
elif extend == 'both':
    # remove the first/last 20 colors, to distinguish the 'under'/'over' color from the main colorbar
    cmap = cmap.from_list('Custom cmap', cmaplist[20:-20], len(cmaplist)-40)
    cmap.set_under(cmaplist[0])
    cmap.set_over(cmaplist[-1])
 
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(12,10))
plt.xticks(range(1,13,1))
yaxis=file.variables['latitude'][:]
xaxis=range(1,13,1)
ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
ax.set_ylim(-20,25)
#ax.set_ylim(47131.5,49000)


c = ax.contourf(xaxis,yaxis,data_new.transpose(), bounds, cmap=cmap, norm=norm, extend='max')

# Add colorbar
cbar = fig.colorbar(c, shrink=0.9, orientation='horizontal', pad=0.1, extend='max', boundaries=bounds, ticks=bounds,label='(mm/day)')


# save plot as a png file
plt.title('ty_10W10E_'+str(var)+')
plt.xlabel('Month')
plt.ylabel('Pr(mm)')
plt.savefig('Varia_intra/'+str(var)+'/ty_10W10E_'+str(var)+'.png')
plt.show()
