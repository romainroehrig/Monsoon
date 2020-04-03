# -*- coding: utf-8 -*-

import numpy as np

from matplotlib import cm

xy_plotdetails = {
        'pr':    {'coef': 86400., 'units': r'mm day$^{-1}$', 'title': 'Precipitation', 
                  'cmap': cm.hsv,      'bounds': [ i for i in np.arange(0,15,1)],
                  'firstwhite': True,  'ext': 'max'},
        'tas':   {'coef': 1.,     'units': r'K', 'title': '2-m Temperature',
                  'cmap': cm.coolwarm, 'bounds': [ i for i in np.arange(270,301,1)],
                  'firstwhite': False, 'ext': 'both'},
        'clt':   {'coef': 1.,     'units': r'%', 'title': 'Total Cloud Fraction',
                  'cmap': cm.Greys, 'bounds': [ i for i in np.arange(0,105,5)],
                  'firstwhite': True,  'ext': 'neither'},
        'prw':   {'coef': 1.,     'units': r'kg m$^{-2}$', 'title': 'Precipitable Water',
                  'cmap': cm.Blues, 'bounds': [ i for i in np.arange(0,60,2)],
                  'firstwhite': False, 'ext': 'max'},
        'ua925': {'coef': 1.,     'units': r'm s$^{-1}$', 'title': 'Zonal Wind at 925 hPa',
                  'cmap': cm.seismic,  'bounds': [ i for i in np.arange(-15,16,1)],
                  'firstwhite': False, 'ext': 'both'},
        'ua700': {'coef': 1.,     'units': r'm s$^{-1}$', 'title': 'Zonal Wind at 700 hPa',
                  'cmap': cm.seismic,  'bounds': [ i for i in np.arange(-20,22,2)],
                  'firstwhite': False, 'ext': 'both'},
        'ua600': {'coef': 1.,     'units': r'm s$^{-1}$', 'title': 'Zonal Wind at 600 hPa',
                  'cmap': cm.seismic,  'bounds': [ i for i in np.arange(-20,22,2)],
                  'firstwhite': False, 'ext': 'both'},
        'ua200': {'coef': 1.,     'units': r'm s$^{-1}$', 'title': 'Zonal Wind at 200 hPa',
                  'cmap': cm.seismic,  'bounds': [ i for i in np.arange(-15,31,1)],
                  'firstwhite': False, 'ext': 'both'},
        'va925': {'coef': 1.,     'units': r'm s$^{-1}$', 'title': 'Meridional Wind at 925 hPa',
                  'cmap': cm.seismic,  'bounds': [ i for i in np.arange(-20,21,1)],
                  'firstwhite': False, 'ext': 'both'},
        'va200': {'coef': 1.,     'units': r'm s$^{-1}$', 'title': 'Meridional Wind at 200 hPa',
                  'cmap': cm.seismic,  'bounds': [ i for i in np.arange(-15,16,1)],
                  'firstwhite': False, 'ext': 'both'},
        'hfss': {'coef': 1.,     'units': 'W.m^(-2)', 'title': 'Sensible',
                  'cmap': cm.hsv,  'bounds': [ i for i in np.arange(0,44,4)],
                  'firstwhite': False, 'ext': 'both'},
        'hfls': {'coef': 1.,     'units': 'W.m^(-2)', 'title': 'Latent',
                  'cmap': cm.hsv,  'bounds': [ i for i in np.arange(200,420,20)],
                  'firstwhite': False, 'ext': 'both'}
}   


