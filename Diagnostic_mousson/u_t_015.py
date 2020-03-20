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
import matplotlib as mpl #Choix de la variable

firstwhite = True 
extend = 'max'

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



var='ua'

data_nc='C:/Users/Ludovic/Desktop/PFE/Simu_WA/DELANNOY_slab_cont_WA_arpsfx_monthly_'+str(var)+'_1850-1855.nc'
file = nc.Dataset(data_nc)
bounds1 =[ i for i in np.arange(-20,22,2)]
unit='(m/s)'
conv=1
ext='both'

#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
#time = file.variables['time'][:]
lev= file.variables['pstd'][:]
#dates = nc.num2date(file['time'][:],units=(file['time']).units,calendar=(file['time']).calendar)
data = (shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
lon=(shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1]
data=np.average(data[5:8:12],axis=0)
bounds1 =[ i for i in np.arange(-20,22,2)]
unit='(m/s)'
conv=1
ext='both'
# Average between 10W and 10E
# select longitude indices between -10 and 10. All parentheses are important
lon_inds = np.ma.where((lon > 0) & (lon < 15))[0]

# final average
data_2 = np.ma.average(data[:,:,lon_inds],axis=2)




xaxis = file.variables['lat'][:]
yaxis = file.variables['pstd'][:]*0.01


var='ta'
data_nc='C:/Users/Ludovic/Desktop/PFE/Simu_WA/DELANNOY_slab_cont_WA_arpsfx_monthly_'+str(var)+'_1850-1855.nc'
file = nc.Dataset(data_nc)
bounds2 =[ i for i in np.arange(-5,6,1)]
unit='K'
cmapref = cm.coolwarm
conv=1
ext='both'


#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
#time = file.variables['time'][:]
lev= file.variables['pstd'][:]
#dates = nc.num2date(file['time'][:],units=(file['time']).units,calendar=(file['time']).calendar)
data = (shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
lon=(shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1]
data=np.average(data[5:8:12],axis=0)


# Average between 10W and 10E
# select longitude indices between -10 and 10. All parentheses are important
lon_inds = np.ma.where((lon > 0) & (lon < 15))[0]
# final average
data = np.ma.average(data[:,:,lon_inds],axis=2)

lat_inds = np.ma.where((lat > -30) & (lat < 30))[0]
data_moy= np.ma.average(data[:,lat_inds],axis=1)

xaxis = file.variables['lat'][:]
yaxis = file.variables['pstd'][:]*0.01



# customize colormap
cmap = cmapref
cmaplist = [cmap(i) for i in range(cmap.N)]
if extend == 'max':
    # remove the last 20 colors, to distinguish the 'over' color from the main colorbar
    cmap = cmap.from_list('Custom cmap', cmaplist[:-20], len(cmaplist)-20)
    cmap.set_over(cmaplist[-1])
elif extend == 'min':
    print('Not yet coded')
    sys.exit()
elif extend == 'both':
    # remove the first/last 20 colors, to distinguish the 'under'/'over' color from the main colorbar
    cmap = cmap.from_list('Custom cmap', cmaplist[20:-20], len(cmaplist)-40)
    cmap.set_under(cmaplist[0])
    cmap.set_over(cmaplist[-1])
 
norm2 = mpl.colors.BoundaryNorm(bounds2, cmap.N)
data_1=np.ma.zeros((len(lev),128))
for i in range(len(lev)):
	for j in range(len(lat)):
		data_1[i,j]=(data[i,j]-data_moy[i])




fig, ax = plt.subplots(figsize=(12,10))
plt.gca().invert_yaxis()

c1 = ax.contour(xaxis,yaxis,data_2, bounds1, levels=range(0,20,2), colors='black',extend=ext,linewidths=2.0,linestyles='solid')
c2 = ax.contour(xaxis,yaxis,data_2, bounds1, levels=range(-20,0,2),colors='black', extend=ext,linewidths=2.0,linestyles='dotted')

norm = mpl.colors.BoundaryNorm(bounds1, cmap.N)
yaxis = file.variables['pstd'][:]*0.01
c = ax.contourf(xaxis,yaxis,data_1 , bounds2, cmap=cmap, norm=norm, extend=ext)
ax.clabel(c1, fontsize=12)
ax.clabel(c2, fontsize=12)
# Make nice labels for x/y axes
xticks = range(5,40,5)
yticks = range(1000,-1,-200)

ax.set_xticks(xticks)
ax.set_yticks(yticks)

ax.grid(True,linestyle='--',lw=1,color='dimgrey')

ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)

ax.set_xlabel('Latitude')
ax.set_ylabel('Pressure (hPa)')

ax.set_xlim(-35,35)
ax.set_ylim(1000,100)

# Add colorbar

cbar = fig.colorbar(c, shrink=0.9, orientation='horizontal', pad=0.1, extend=ext, boundaries=bounds2, ticks=bounds2,label=(str(unit)))

plt.xlabel('Latitude')
plt.ylabel('Level(Hpa)')
plt.title('u_t_015')
plt.savefig('u_t_015.png')
plt.show()

# close frame
