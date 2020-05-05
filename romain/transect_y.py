# -*- coding: utf-8 -*-

import os

import netCDF4 as nc

import numpy as np
import numpy.ma as ma

import matplotlib.pyplot as plt

from cartopy.mpl.gridliner import LATITUDE_FORMATTER

from data import getfile
from seasons import seasons
from variables import varinfo
#import config as cc

def getlonname(lon):
    if lon > 0:
        return '{0}E'.format(int(lon))
    elif lon < 0:
        return '{0}W'.format(int(-lon))
    else:
        return '0'

def transect_y(sim,var,lev=None,rep0='./images/',
        lonavg=(-180,180),
        continent=None,season='year',
        xlim=None,ylim=None,
        xticks=None,yticks=None,
        clats=None,addAll=False):

    lonavg_name = '{0}-{1}'.format(getlonname(lonavg[0]),getlonname(lonavg[1]))

    repout = '{0}/{1}/{2}avg/'.format(rep0,sim,season)
    if not(os.path.exists(repout)):
        os.makedirs(repout)

    if lev is None:
        varname = var
    else:
        varname = '{0}{1}'.format(var,int(lev))

    # Get plot details
    coef       = varinfo[var]['coef']
    units      = varinfo[var]['units']
#    cmap       = cc.xy_plotdetails[varname]['cmap']
#    bounds     = cc.xy_plotdetails[varname]['bounds']
#    firstwhite = cc.xy_plotdetails[varname]['firstwhite']
#    ext        = cc.xy_plotdetails[varname]['ext']
    if lev is None:
        title = '{0}'.format(varinfo[var]['name'])
    else:
        title = '{0} at {1} hPa'.format(varinfo[var]['name'],int(lev))

    # Ouverture du fichier
    if season == 'year':
        dtype = 'timmean_1860-1869'
    else:
        dtype = '{0}mean_1860-1869'.format(season)
    data_nc = getfile(sim,var,dtype=dtype)
    fnc = nc.Dataset(data_nc)


    # Lecture des variables
    lat = fnc.variables['lat'][:]
    lon = fnc.variables['lon'][:]

#    for taxis in ['time','time_counter']:
#        try:
#            time = fnc.variables[taxis]
#            break
#        except:
#            pass
#    dates = nc.num2date(time[:],units=time.units,calendar=time.calendar)

    data = fnc.variables[var]
    if len(data.shape) == 4:
        levname = 'pstd'
        levaxis = fnc.variables[levname][:]
        data = np.squeeze(data[:,levaxis == lev*100.,:,:])*coef
    else:
        data = np.squeeze(data[:,:,:])*coef

    # average over the chosen season
#    inds = [d.month in seasons[season] for d in dates]
#    datanew=np.average(datanew[inds],axis=0)

    # shifting grid from 0-360 to -180-180
    if np.max(lon) > 180.:
        lon = np.where(lon > 180, lon-360, lon)

    # average over lonavg
    lon_inds = np.where((lon > lonavg[0]) & (lon <= lonavg[1]))[0]

    datanew= np.ma.average(data[:,lon_inds],axis=1)

    if addAll:
        dataAll = np.ma.average(data,axis=1)

    # Mise en place de la carte

    # Plot
    fig, ax = plt.subplots(figsize=(10,5))

    if not(xlim is None):
        ax.set_xlim(xlim)
    if not(ylim is None):
        ax.set_ylim(ylim)

    if not(xticks is None):
        ax.set_xticks(xticks)
        ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)

    ax.plot(lat,datanew,color='r')

    if addAll:
        ax.plot(lat,datanew,color='r',label='[{0}]'.format(lonavg_name))
        ax.plot(lat,dataAll,color='b',label='[180W-180E]')
        plt.legend(loc='best')
    else:
        ax.plot(lat,datanew,color='k')

    if continent == 'ideal' and not(clats is None):
        ylim = ax.get_ylim()
        yyc = ylim[0],ylim[0]
        ax.plot(clats,yyc,color='k',linewidth=5.)
        ax.set_ylim(ylim)

    plt.title('{0} ({1}) - {2} [{3}]'.format(title,units,season,lonavg_name))

    plt.savefig('{0}/y_{1}_{2}.png'.format(repout,lonavg_name,varname))
    print '{0}/y_{1}_{2}.png'.format(repout,lonavg_name,varname)
    plt.close()
