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
from data_atlas import getfile
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



var='wap'
sim='WA_10a'
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
data_nc=getfile(sim,var)

file = nc.Dataset(data_nc)

var=='wap'
bounds =[ i for i in np.arange(-0.15,0.16,0.01)]
unit='(Pa/s)'
cmapref = cm.seismic
conv=1
ext='both'


#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
#time = file.variables['time'][:]
lev= file.variables[str(lv)][:]
dates = nc.num2date(file[str(t)][:],units=(file[str(t)]).units,calendar=(file[str(t)]).calendar)
data = (shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
lon=(shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1]
inds = [d.month in seasons['JJAS'] for d in dates]
data=np.average(data[inds],axis=0)
# Average between 10W and 10E
# select longitude indices between -10 and 10. All parentheses are important
lon_inds = np.where((lon > 0) & (lon < 30))[0]

# final average
data = np.ma.average(data[:,:,lon_inds],axis=2)




xaxis = file.variables['lat'][:]
yaxis = file.variables[str(lv)][:]*0.01



# customize colormap
cmap = cmapref
cmaplist = [cmap(i) for i in range(cmap.N)]
if extend == 'max':
    # remove the last 20 colors, to distinguish the 'over' color from the main colorbar
    cmap = cmap.from_list('Custom cmap', cmaplist[:-20], len(cmaplist)-20)
    cmap.set_over(cmaplist[-1])
elif extend == 'min':
    print ('Not yet coded')
    sys.exit()
elif extend == 'both':
    # remove the first/last 20 colors, to distinguish the 'under'/'over' color from the main colorbar
    cmap = cmap.from_list('Custom cmap', cmaplist[20:-20], len(cmaplist)-40)
    cmap.set_under(cmaplist[0])
    cmap.set_over(cmaplist[-1])
 
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots(figsize=(12,10))
plt.gca().invert_yaxis()
c = ax.contourf(xaxis,yaxis,data[:,:], bounds, cmap=cmap, norm=norm, extend=ext)

# Make nice labels for x/y axes
xticks = range(-35,40,5)
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
cbar = fig.colorbar(c, shrink=0.9, orientation='horizontal', pad=0.1, extend=ext, boundaries=bounds, ticks=bounds,label=(str(unit)))

# save plot as a png file
plt.xlabel('Latitude')
plt.ylabel('Level(Hpa)')
plt.title('w_030'+str(sim))
plt.savefig('images/w/w_030_'+str(sim)+'.png')


# close frame
