# -*- coding: utf-8 -*-

import numpy as np

from matplotlib import cm

xy_plotdetails = {
        'pr':    {'cmap': cm.coolwarm_r,'bounds': [ i for i in np.arange(0,13,1)],
                  'firstwhite': True,   'ext': 'max'},
        'tas':   {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(270,301,2)],
                  'firstwhite': False,  'ext': 'both'},
        'clt':   {'cmap': cm.Greys_r,   'bounds': [ i for i in np.arange(0,105,5)],
                  'firstwhite': False,  'ext': 'neither'},
        'prw':   {'cmap': cm.Blues,     'bounds': [ i for i in np.arange(0,60,2)],
                  'firstwhite': False,  'ext': 'max'},
        'ua925': {'cmap': cm.seismic,   'bounds': [ i for i in np.arange(-20,22,2)],
                  'firstwhite': False,  'ext': 'both'},
        'ua700': {'cmap': cm.seismic,   'bounds': [ i for i in np.arange(-20,22,2)],
                  'firstwhite': False,  'ext': 'both'},
        'ua600': {'cmap': cm.seismic,   'bounds': [ i for i in np.arange(-20,22,2)],
                  'firstwhite': False,  'ext': 'both'},
        'ua200': {'cmap': cm.seismic,   'bounds': [ i for i in np.arange(-60,61,5)],
                  'firstwhite': False,  'ext': 'both'},
        'va925': {'cmap': cm.seismic,   'bounds': [ i for i in np.arange(-4,4.5,0.5)],
                  'firstwhite': False,  'ext': 'both'},
        'va700': {'cmap': cm.seismic,   'bounds': [ i for i in np.arange(-3,3.5,0.5)],
                  'firstwhite': False,  'ext': 'both'},
        'va200': {'cmap': cm.seismic,   'bounds': [ i for i in np.arange(-6,7,1)],
                  'firstwhite': False,  'ext': 'both'}
        }
