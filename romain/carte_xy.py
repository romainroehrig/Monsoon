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
from data import getfile
from seasons import seasons
import config as cc


def carte_xy(sim,var,lev=None,rep0='./images/',
        continent=None,season='year',
        xlim=None,ylim=None,
        xticks=None,yticks=None,
        xxc=None,yyc=None):

    repout = '{0}/{1}/{2}avg/'.format(rep0,sim,season)
    if not(os.path.exists(repout)):
        os.makedirs(repout)

    if lev is None:
        varname = var
    else:
        varname = '{0}{1}'.format(var,int(lev))

    # Get plot details
    coef       = cc.xy_plotdetails[varname]['coef']
    units      = cc.xy_plotdetails[varname]['units']
    cmap       = cc.xy_plotdetails[varname]['cmap']
    bounds     = cc.xy_plotdetails[varname]['bounds']
    firstwhite = cc.xy_plotdetails[varname]['firstwhite']
    ext        = cc.xy_plotdetails[varname]['ext']
    title      = cc.xy_plotdetails[varname]['title']

    # Ouverture du fichier
    data_nc = getfile(sim,var)
    fnc = nc.Dataset(data_nc)


    # Lecture des variables
    lat = fnc.variables['lat'][:]
    lon = fnc.variables['lon'][:]

    for taxis in ['time','time_counter']:
        try:
            time = fnc.variables[taxis]
            break
        except:
            pass
    dates = nc.num2date(time[:],units=time.units,calendar=time.calendar)

    data = fnc.variables[var]
    if len(data.shape) == 4:
        levname = 'pstd'
        levaxis = fnc.variables[levname][:]
        data = np.squeeze(data[:,levaxis == lev*100.,:,:])*coef
    else:
        data = data[:,:,:]*coef

    # shifting grid from 0-360 to -180-180
    datanew, lonnew = shiftgrid(180,data,lon)

    # average over the chosen season
    inds = [d.month in seasons[season] for d in dates]
    datanew=np.average(datanew[inds],axis=0)

    # Mise en place de la carte
    proj = ccrs.PlateCarree(central_longitude=0, globe=None)

    # prepare colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    if firstwhite:
        cmaplist[0] = (1,1,1,0)
    cmap = cmap.from_list('Custom cmap', cmaplist[:-20], len(cmaplist)-20)
    norm = mpl.colors.BoundaryNorm(bounds,cmap.N)

    # Plot
    fig = plt.figure(figsize=(10,10))
    fig, ax = plt.subplots(figsize=(12,10),subplot_kw=dict(projection=proj))

    if not(xlim is None):
        ax.set_xlim(xlim)
    if not(ylim is None):
        ax.set_ylim(ylim)

    if not(xticks is None) and not(yticks is None):
        gl = ax.gridlines(xlocs=xticks, ylocs=yticks,linestyle='--',lw=1,color='dimgrey',draw_labels=True)
        gl.xlabels_top = False
        gl.xformatter = LONGITUDE_FORMATTER
        #gl.xlines = False
        gl.ylabels_right = False
        gl.yformatter = LATITUDE_FORMATTER
        #gl.ylines = False

    cs = ax.contourf(lonnew,lat, datanew, bounds, cmap=cmap, norm=norm, extend=ext, transform=ccrs.PlateCarree())

    cbar = fig.colorbar(cs, shrink=0.9, orientation='horizontal', pad=0.1, extend=ext, boundaries=bounds, ticks=bounds , label=units)

    if continent == 'real':
        ax.coastlines('50m')
        ax.add_feature(cf.BORDERS)
    elif continent == 'ideal':
        # Add idealized continent
        ax.plot(xxc,yyc,color='k',linewidth=2,transform=ccrs.PlateCarree())

    plt.title('{0} ({1}) - {2}'.format(title,units,season))

    plt.savefig('{0}/xy_{1}.png'.format(repout,varname))
    print '{0}/xy_{1}.png'.format(repout,varname)
    plt.close()
