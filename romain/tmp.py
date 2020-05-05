# -*- coding: utf-8 -*-

import os

import netCDF4 as nc

import numpy as np
import numpy.ma as ma

import matplotlib.pyplot as plt
from matplotlib import cm

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

def transect_y(sims,var,lev=None,rep0='./',img_name='y',
        lonavg=(-180,180),
        continent=None,season='year',
        xlim=None,ylim=None,
        xticks=None,yticks=None,
        clats=None,
        linecolors=None,labels=None):

    lonavg_name = '{0}-{1}'.format(getlonname(lonavg[0]),getlonname(lonavg[1]))

    repout = '{0}/'.format(rep0)
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

    data_nc = {}
    fnc = {}
    data = {}
    lat = {}

    for sim in sims:
        data_nc[sim] = getfile(sim,var,dtype=dtype)
        fnc[sim] = nc.Dataset(data_nc[sim])


        # Lecture des variables
        lat[sim] = fnc[sim].variables['lat'][:]
        lon = fnc[sim].variables['lon'][:]

#    for taxis in ['time','time_counter']:
#        try:
#            time = fnc.variables[taxis]
#            break
#        except:
#            pass
#    dates = nc.num2date(time[:],units=time.units,calendar=time.calendar)

        data[sim] = fnc[sim].variables[var]
        if len(data[sim].shape) == 4:
            levname = 'pstd'
            levaxis = fnc[sim].variables[levname][:]
            data[sim] = np.squeeze(data[sim][:,levaxis == lev*100.,:,:])*coef
        else:
            data[sim] = np.squeeze(data[sim][:,:,:])*coef

    # average over the chosen season
#    inds = [d.month in seasons[season] for d in dates]
#    datanew=np.average(datanew[inds],axis=0)

    # shifting grid from 0-360 to -180-180
        if np.max(lon) > 180.:
            lon = np.where(lon > 180, lon-360, lon)

    # average over lonavg
        lon_inds = np.where((lon > lonavg[0]) & (lon <= lonavg[1]))[0]

        data[sim] = np.ma.average(data[sim][:,lon_inds],axis=1)

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

    if labels is None:
        labels = []
        for sim in sims:
            labels.append(sim)

    if linecolors is None:
        for i,sim in enumerate(sims):
             ax.plot(lat[sim],data[sim],label=labels[i])
    else:
        for i,sim in enumerate(sims):
             ax.plot(lat[sim],data[sim],color=linecolors[i],label=labels[i])

    plt.legend(loc='best')

    if continent == 'ideal' and not(clats is None):
        ylim = ax.get_ylim()
        yyc = ylim[0],ylim[0]
        ax.plot(clats,yyc,color='k',linewidth=5.)
        ax.set_ylim(ylim)

    plt.title('{0} ({1}) - {2} [{3}]'.format(title,units,season,lonavg_name))

    plt.savefig('{0}/{1}_{2}_{3}.png'.format(repout,img_name,lonavg_name,varname))
    print '{0}/{1}_{2}_{3}.png'.format(repout,img_name,lonavg_name,varname)
    plt.close()

simroot = 'DELANNOY_slab_cont_WA_evap'
simroot = 'DELANNOY_slab_cont_WA_qflux_evap'
simroot = 'DELANNOY_slab_cont_tracmip_evap'
simulations = ['{0}{1:0>2}'.format(simroot,i) for i in range(0,11)]

season = 'JJAS'
lonavg = (0,60)
var = 'pr'

xlim = (-40,40)
xticks = range(-40,41,10)

continent = 'ideal'
cont = 'WA'
cont = 'tracmip'
if cont == 'WA':
    clats = (5,35)
elif cont == 'tracmip':    
    clats = (-30,30)

cmap = cm.coolwarm
linecolors = []
labels = []
for i in range(0,len(simulations)):
    linecolors.append(cmap(i*int(cmap.N/len(simulations))))
    labels.append(i/10.)

transect_y(simulations,var,rep0='./',img_name='{0}_y'.format(simroot),
        lonavg=lonavg,
        continent=continent,season=season,
        xlim=xlim,ylim=None,
        xticks=xticks,yticks=None,
        clats=clats,linecolors=linecolors,labels=labels)
