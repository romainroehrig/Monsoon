# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 15:51:35 2020

@author: Ludovic
"""
 
import os

import math

import netCDF4 as nc

import numpy as np

from myfunctions import shiftgrid
#from data_atlas import getfile
from seasons import seasons

# Constants

Cp = 1004. # J kg-1 K-1
g = 9.80665 # m s-2
Lv = 2.501e6 # J kg-1 (at 0°C)

Rt = 6371000 # m

# Variables

var2read = ['ta','hus','zg','ua','va','wap','ps']

#sim='Simu_Base'
#cont='real'	
seasons = {'year': range(1,13),
           'JJAS': [6,7,8,9],
           }
#if sim=='Simu_Base':
#    t='time'
#else:
#    t='time_counter'
    
#Lecture des fichiers

data = {}

for var in var2read:
    tmp = '/Volumes/CNRM/TMP/simulations/DELANNOY/DELANNOY_slab_cont_WA_evap00/ymonmean_1860-1869/{0}.nc'.format(var) #getfile(sim,var)
    fnc = nc.Dataset(tmp)

    lat = fnc['lat'][:]
    
    if var == 'ps':
        data[var],lon = shiftgrid(180,fnc[var][:,:,:],fnc['lon'][:],start=False,cyclic=360.)
    else:
        lev = fnc['pstd'][:]
        data[var],lon = shiftgrid(180,fnc[var][:,:,:,:],fnc['lon'][:],start=False,cyclic=360.)

    dates = nc.num2date(fnc['time_counter'][:],units=fnc['time_counter'].units,calendar=fnc['time_counter'].calendar)

    inds = [d.month in seasons['JJAS'] for d in dates]
    # Average over JJAS
    data[var] = np.ma.average(data[var][inds],axis=0)

    fnc.close()


mse = Cp*data['ta']+g*data['zg']+Lv*data['hus'] # J kg-1

nlev,nlat,nlon = mse.shape

dlon = (lon[1:]-lon[:-1])*math.pi/180.*Rt
dx = mse*0. # Note we do not consider the upper boundary
for ilat in range(0,nlat):
    for ilon in range(0,nlon-1):
        dx[:,ilat,ilon] = dlon[ilon]*math.cos(lat[ilat]*math.pi/180.)

dlat = (lat[1:]-lat[:-1])*math.pi/180.*Rt
dy = mse*0.
for ilat in range(0,nlat-1):
    dy[:,ilat,:] = dlat[ilat]

dlev = lev[1:]-lev[:-1]
dp = mse*0.
for ilev in range(0,nlev-1):
    dp[ilev,:,:] = dlev[ilev]

upos = np.ma.maximum(data['ua'],0)
uneg = np.ma.minimum(data['ua'],0)

vpos = np.ma.maximum(data['va'],0)
vneg = np.ma.minimum(data['va'],0)

wpos = np.ma.maximum(data['wap'],0)
wneg = np.ma.minimum(data['wap'],0)

mse_dx = mse*0.
mse_dx[:,:,:-1] = (mse[:,:,1:]-mse[:,:,:-1])/dx[:,:,:-1]

mse_advu = mse*0.
mse_advu[:,:,1:] = upos[:,:,1:]*mse_dx[:,:,:-1]
mse_advu[:,:,:-1] = mse_advu[:,:,:-1] + uneg[:,:,:-1]*mse_dx[:,:,:-1]

mse_dy = mse*0.
mse_dy[:,:-1,:] = (mse[:,1:,:]-mse[:,:-1,:])/dy[:,:-1,:]

mse_advv = mse*0.
mse_advv[:,1:,:] = vpos[:,1:,:]*mse_dy[:,:-1,:]
mse_advv[:,:-1,:] = mse_advv[:,:-1,:] + vneg[:,:-1,:]*mse_dy[:,:-1,:]

mse_dp = mse*0.
mse_dp[:-1,:,:] = (mse[1:,:,:]-mse[:-1,:,:])/dp[:-1,:,:]

mse_advw = mse*0.
mse_advw[1:,:,:] = wpos[1:,:,:]*mse_dp[:-1,:,:] + wneg[:-1,:,:]*mse_dp[:-1,:,:]
mse_advw[:-1,:,:] = mse_advw[:-1,:,:] + wneg[:-1,:,:]*mse_dp[:-1,:,:]

# Vertical integration

int_mse_advu = np.zeros((nlat,nlon))
int_mse_advv = np.zeros((nlat,nlon))
int_mse_advw = np.zeros((nlat,nlon))

int_mse_advu[:,:] = mse_advu[0,:,:]*(data['ps'][:,:]-(lev[0]+lev[1])/2.)
int_mse_advv[:,:] = mse_advv[0,:,:]*(data['ps'][:,:]-(lev[0]+lev[1])/2.)
int_mse_advw[:,:] = mse_advw[0,:,:]*(data['ps'][:,:]-(lev[0]+lev[1])/2.)


for ilev in range(1,nlev-1):
    int_mse_advu[:,:] = int_mse_advu[:,:] + mse_advu[ilev,:,:]*(lev[ilev-1]-lev[ilev+1])/2.
    int_mse_advv[:,:] = int_mse_advv[:,:] + mse_advv[ilev,:,:]*(lev[ilev-1]-lev[ilev+1])/2.
    int_mse_advw[:,:] = int_mse_advw[:,:] + mse_advw[ilev,:,:]*(lev[ilev-1]-lev[ilev+1])/2.

int_mse_advu = int_mse_advu/data['ps']
int_mse_advv = int_mse_advv/data['ps']
int_mse_advw = int_mse_advw/data['ps']

int_mse_advhor = mse_advu + mse_advv

# Saving data
fout = nc.Dataset('/Volumes/CNRM/TMP/simulations/DELANNOY/DELANNOY_slab_cont_WA_evap00/ymonmean_1860-1869/output.nc','w')
fout.createDimension('lev', nlev)
fout.createDimension('lat', nlat)
fout.createDimension('lon', nlon)

levAxis = fout.createVariable('lev','f4',('lev',))
levAxis[:] = lev[:]
latAxis = fout.createVariable('lat','f4',('lat',))
latAxis[:] = lat[:]
lonAxis = fout.createVariable('lon','f4',('lon',))
lonAxis[:] = lon[:]


nc_mse = fout.createVariable('mse','f4',('lev','lat','lon',))
nc_mse[:,:,:] =  mse[:,:,:]

nc_mse_dx = fout.createVariable('mse_dx','f4',('lev','lat','lon',))
nc_mse_dx[:,:,:] =  mse_dx[:,:,:]

nc_mse_dy = fout.createVariable('mse_dy','f4',('lev','lat','lon',))
nc_mse_dy[:,:,:] =  mse_dy[:,:,:]

nc_mse_dp = fout.createVariable('mse_dp','f4',('lev','lat','lon',))
nc_mse_dp[:,:,:] =  mse_dp[:,:,:]


nc_mse_advu = fout.createVariable('mse_advu','f4',('lev','lat','lon',))
nc_mse_advu[:,:,:] =  mse_advu[:,:,:]

nc_int_mse_advu = fout.createVariable('int_mse_advu','f4',('lat','lon',))
nc_int_mse_advu[:,:] = int_mse_advu[:,:]

nc_mse_advv = fout.createVariable('mse_advv','f4',('lev','lat','lon',))
nc_mse_advv[:,:,:] =  mse_advv[:,:,:]

nc_int_mse_advv = fout.createVariable('int_mse_advv','f4',('lat','lon',))
nc_int_mse_advv[:,:] = int_mse_advv[:,:]

nc_mse_advw = fout.createVariable('mse_advw','f4',('lev','lat','lon',))
nc_mse_advw[:,:,:] =  mse_advw[:,:,:]

nc_int_mse_advw = fout.createVariable('int_mse_advw','f4',('lat','lon',))
nc_int_mse_advw[:,:] = int_mse_advw[:,:]

fout.close()


