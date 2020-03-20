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

var='wap'

#Lecture du fichier
data_nc = 'C:/Users/Ludovic/Desktop/PFE/Amon_SST/'+str(var)+'/gr/files/d20180711/'+str(var)+'_Amon_CNRM-CM6-1_amip_r1i1p1f2_gr_197901-201412.nc' 
 
file = nc.Dataset(data_nc)

if var=='ua':
	bounds =[ i for i in np.arange(-20,22,2)]
	unit='(m/s)'
	cmapref = cm.seismic
	conv=1
	ext='both'
elif var=='va':
	bounds =[ i for i in np.arange(-5,6,1)]
	unit='(m/s)'
	cmapref = cm.seismic
	conv=1
	ext='both'
elif var=='wap':
	bounds =[ i for i in np.arange(-0.15,0.16,0.02)]
	unit='(Pa/s)'
	cmapref = cm.seismic
	conv=1
	ext='both'
elif var=='ta':
	bounds =[ i for i in np.arange(190,310,5)]
	unit='K'
	cmapref = cm.coolwarm
	conv=1
	ext='both'
elif var=='hur':
	bounds =[ i for i in np.arange(0,105,5)]
	unit='%'
	cmapref = cm.RdBu
	conv=1
	ext='neither'

#Lecture des variables
lat = file.variables['lat'][:]
lon = file.variables['lon'][:]
time = file.variables['time'][:]
lev= file.variables['plev'][:]
print(lev)
dates = nc.num2date(file['time'][:],units=(file['time']).units,calendar=(file['time']).calendar)
data = (shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[0]
lon=(shiftgrid(180,file.variables[str(var)][:,:,:,:]*conv,nc.Dataset(data_nc).variables['lon'][:],start=False,cyclic=360.0))[1]
data=np.average(data[5:8:12],axis=0)

# Average between 10W and 10E
# select longitude indices between -10 and 10. All parentheses are important
lon_inds = np.where((lon > -10) & (lon < 10))[0]

# final average
data = np.ma.average(data[:,:,lon_inds],axis=2)




xaxis = file.variables['lat'][:]
yaxis = file.variables['plev'][:]*0.01



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
xticks = range(-90,91,10)
yticks = range(1000,-1,-200)

ax.set_xticks(xticks)
ax.set_yticks(yticks)

ax.grid(True,linestyle='--',lw=1,color='dimgrey')

ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)

ax.set_xlabel('Latitude')
ax.set_ylabel('Pressure (hPa)')

ax.set_xlim(-20,40)
ax.set_ylim(1000,0)

# Add colorbar
cbar = fig.colorbar(c, shrink=0.9, orientation='horizontal', pad=0.1, extend=ext, boundaries=bounds, ticks=bounds,label=(str(unit)))

# save plot as a png file
plt.title('yz_10W10E_'+str(var)+'_JJAS')
plt.xlabel('Latitude')
plt.ylabel('Level(Hpa)')
plt.savefig('Coupe_latalt/'+str(var)+'/yz_10W10E_'+str(var)+'_JJAS.png')
plt.show()

# close frame
