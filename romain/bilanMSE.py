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
from data import getfile
#from seasons import seasons

sim = 'DELANNOY_slab_cont_WA_evap00'
season = 'JJAS'

ss = season
if season == 'year':
    ss = 'tim'

# Constants

Cp = 1004. # J kg-1 K-1
g = 9.80665 # m s-2
Lv = 2.501e6 # J kg-1 (at 0°C)

Rt = 6371000 # m

# Variables

var2read = ['ta','hus','zg','ua','va','wap','ps','rsdt','rsut','rsds','rsus','rlut','rlds','rlus','hfls','hfss']

#Lecture des fichiers

data = {}

for var in var2read:
    tmp = getfile(sim,var,dtype='{0}mean_1860-1869'.format(ss)) #'/Volumes/CNRM/TMP/simulations/DELANNOY/DELANNOY_slab_cont_WA_evap00/ymonmean_1860-1869/{0}.nc'.format(var) #getfile(sim,var)
    print var, tmp
    fnc = nc.Dataset(tmp)

    lat = fnc['lat'][:]
    
    if 'pstd' in fnc.dimensions.keys():
        lev = fnc['pstd'][:]
        data[var],lon = shiftgrid(180,np.ma.squeeze(fnc[var][:,:,:,:]),fnc['lon'][:],start=False,cyclic=360.)
    else:
        data[var],lon = shiftgrid(180,np.ma.squeeze(fnc[var][:,:,:]),fnc['lon'][:],start=False,cyclic=360.)

#    dates = nc.num2date(fnc['time_counter'][:],units=fnc['time_counter'].units,calendar=fnc['time_counter'].calendar)

#    inds = [d.month in seasons['JJAS'] for d in dates]
    # Average over JJAS
#    data[var] = np.ma.average(data[var][inds],axis=0)

    fnc.close()

tmp = tmp.split('/')[:-1]
dirout = os.path.join('/',*tmp)

# Computing MSE

mse = Cp*data['ta']+g*data['zg']+Lv*data['hus'] # J kg-1

nlev,nlat,nlon = mse.shape

# Computing MSE gradients and advections

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

data['mse'] = mse

data['mse_dx'] = mse_dx
data['mse_dy'] = mse_dy
data['mse_dp'] = mse_dp

data['mse_advu'] = mse_advu
data['mse_advv'] = mse_advv
data['mse_advw'] = mse_advw
data['mse_advh'] = mse_advu + mse_advv
data['mse_adv'] = data['mse_advh'] + mse_advw


# Vertical integration

for var in ['u','v','w']:
    varloc = 'mse_adv' + var
    varint = 'int_mse_adv' + var

    data[varint] = data[varloc][0,:,:]*(data['ps'][:,:]-(lev[0]+lev[1])/2.)

    for ilev in range(1,nlev-1):
        data[varint][:,:] = data[varint][:,:] + data[varloc][ilev,:,:]*(lev[ilev-1]-lev[ilev+1])/2./g # W m-2

data['int_mse_advh'] = data['int_mse_advu'] + data['int_mse_advv']
data['int_mse_adv'] = data['int_mse_advh'] + data['int_mse_advw']

# Other data

data['rsa'] = data['rsdt'] - data['rsut'] - (data['rsds'] - data['rsus'])
data['rla'] = 0. - data['rlut'] - (data['rlds'] - data['rlus'])
data['rna'] = data['rsa'] + data['rla']
data['hfs'] = data['hfls'] + data['hfss']
data['fnet'] = data['rna'] + data['hfs']

# Residual

data['res'] = data['fnet'] - data['int_mse_adv']

# Saving data
var2save = ['mse'] + ['mse_d'+v for v in ['x','y','p']]\
        + ['mse_adv'+v for v in ['u','v','w']]\
        + ['int_mse_adv'+v for v in ['u','v','w']]\
        + ['mse_advh','int_mse_advh','mse_adv','int_mse_adv']\
        + ['rsa','rla','rna','hfs','fnet']\
        + ['res']

for var in var2save:
    fout = nc.Dataset('/{0}/{1}.nc'.format(dirout,var),'w')

    lvar3D = len(data[var].shape) == 3

    if lvar3D:
        fout.createDimension('pstd', nlev)
        levAxis = fout.createVariable('pstd','f4',('pstd',))
        levAxis[:] = lev[:]

    fout.createDimension('lat', nlat)
    latAxis = fout.createVariable('lat','f4',('lat',))
    latAxis[:] = lat[:]

    fout.createDimension('lon', nlon)
    lonAxis = fout.createVariable('lon','f4',('lon',))
    lonAxis[:] = lon[:]

    if lvar3D:
        nc_var = fout.createVariable(var,'f4',('pstd','lat','lon',))
        nc_var[:,:,:] = data[var][:,:,:]
    else:
        nc_var = fout.createVariable(var,'f4',('lat','lon',))
        nc_var[:,:] = data[var][:,:]

    fout.close()

