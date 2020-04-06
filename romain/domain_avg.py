# -*- coding: utf-8 -*-

import os
import sys

import netCDF4 as nc

import numpy as np
import numpy.ma as ma

import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.dates as mdates

from data import getfile
from seasons import seasons
from domains import domains
from variables import varinfo

years = mdates.YearLocator()
months = mdates.MonthLocator()
years_fmt = mdates.DateFormatter('%Y')

def domain_avg(sim,var,lev=None,rep0='./images/',
        domain='global',season=None,
        ylim=None,yticks=None):

    if season is None:
        repout = '{0}/{1}/{2}_avg/'.format(rep0,sim,domain)
    else:
        repout = '{0}/{1}/{2}_avg_{3}/'.format(rep0,sim,domain,season)
    if not(os.path.exists(repout)):
        os.makedirs(repout)

    if lev is None:
        varname = var
    else:
        varname = '{0}{1}'.format(var,int(lev))

    # Get plot details
    coef       = varinfo[var]['coef']
    units      = varinfo[var]['units']
    if lev is None:
        title = varinfo[var]['name']
    else:
        title = '{0} at {1} hPa'.format(varinfo[var]['name'],lev)

    # Ouverture du fichier
    data_nc = getfile(sim,var)
    fnc = nc.Dataset(data_nc)


    # Lecture des variables
    lat = fnc.variables['lat'][:]

    # select latitude in the required domain. All parentheses are important
    latmin = domains[domain]['latmin']
    latmax = domains[domain]['latmax']
    lat_inds = np.where((lat >= latmin) & (lat <= latmax))[0]

    #print lat[lat_inds]

    lon = fnc.variables['lon'][:]

    # in case longitude goes from -180 to 180, make them between 0 and 360
    if np.min(lon) < 0.:
        lon = np.where(lon < 0., lon+360, lon)

    # select longitude in the required domain. All parentheses are important
    lonmin = domains[domain]['lonmin']
    lonmax = domains[domain]['lonmax']
    lon_inds = np.where((lon >= lonmin) & (lon <= lonmax))[0]

    #print lon[lon_inds]

    for taxis in ['time','time_counter']:
        try:
            time = fnc.variables[taxis]
            break
        except:
            pass
    dates = nc.num2date(time[:],units=time.units,calendar=time.calendar,
            only_use_cftime_datetimes=False,only_use_python_datetimes=True)
    #dd = [i.strftime("%Y-%m-%d %H:%M") for i in dates]
    #print dd

    data = fnc.variables[var]
    if len(data.shape) == 4:
        levname = 'pstd'
        levaxis = fnc.variables[levname][:]
        data = np.squeeze(data[:,levaxis == lev*100.,lat_inds,lon_inds])*coef
    else:
        data = data[:,lat_inds,lon_inds]*coef

    # Domain and season average
    if season is None:
        datanew = np.average(data,axis=(1,2))
    else:
        print 'season not yet coded:', season
        sys.exit()



    # Plot
    fig, ax = plt.subplots()

    ax.plot(dates,datanew,color='k',linewidth=2)

    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(years_fmt)
    ax.xaxis.set_minor_locator(months)

    datemin = np.datetime64(dates[0],'Y')
    datemax = np.datetime64(dates[-1],'Y')

    ax.set_xlim(datemin,datemax)

    if not(ylim is None):
        ax.set_ylim(ylim)

    if not(yticks is None):
        ax.set_yticks(yticks)    

    ax.set_ylabel(units)

    if season is None:
        ax.set_title('{0} - {2}'.format(title,units,domain))
    else:
        ax.set_title('{0} - {2}, {3)'.format(title,units,domain,season))

    fig.autofmt_xdate()

    plt.savefig('{0}/{1}.png'.format(repout,varname))
    print '{0}/{1}.png'.format(repout,varname)
    plt.close()
