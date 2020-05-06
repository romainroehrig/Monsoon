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
        'fnet':  {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(-100,101,10)],
                  'firstwhite': False,  'ext': 'both'},
        'rsa':   {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(0,111,5)],
                  'firstwhite': False,  'ext': 'both'},
        'rla':   {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(-220,1,10)],
                  'firstwhite': False,  'ext': 'both'},
        'rna':   {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(-180,1,10)],
                  'firstwhite': False,  'ext': 'both'},
        'hfs':   {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(0,201,10)],
                  'firstwhite': False,  'ext': 'both'},
        'hfls':  {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(0,201,10)],
                  'firstwhite': False,  'ext': 'both'},
        'hfss':  {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(0,121,5)],
                  'firstwhite': False,  'ext': 'both'},
        'int_mse_advh':  {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(-240,241,20)],
                  'firstwhite': False,  'ext': 'both'},
        'int_mse_advw':  {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(-120,121,10)],
                  'firstwhite': False,  'ext': 'both'},
        'res':   {'cmap': cm.coolwarm,  'bounds': [ i for i in np.arange(-200,201,20)],
                  'firstwhite': False,  'ext': 'both'},
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
                  'firstwhite': False,  'ext': 'both'},
        'wap400': {'cmap': cm.seismic,   'bounds': [ i for i in np.arange(-160,161,20)],
                  'firstwhite': False,  'ext': 'both'},
        'wap850': {'cmap': cm.seismic,   'bounds': [ i for i in np.arange(-160,161,20)],
                  'firstwhite': False,  'ext': 'both'},
        'mse925': {'cmap': cm.coolwarm,   'bounds': [ i for i in np.arange(260,341,5)],
                  'firstwhite': False,  'ext': 'both'},
        'mse200': {'cmap': cm.coolwarm,   'bounds': [ i for i in np.arange(300,343,3)],
                  'firstwhite': False,  'ext': 'both'}        
        }
