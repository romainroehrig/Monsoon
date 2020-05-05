#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os

import argparse

from carte_xy import carte_xy
from transect_y import transect_y

repout = './images/'
if not(os.path.exists(repout)):
    os.makedirs(repout)

parser = argparse.ArgumentParser()
parser.add_argument("-sim", help="simulations to plot", type=str, required=True)
parser.add_argument("-cont", help="continent type (WA, tracmip)", type=str, required=False, default='no')
args = parser.parse_args()

#sim = 'DELANNOY_slab_cont_WA'
sim = args.sim
cont = args.cont


#### Cartes global

seasons = ['year','JJAS']

var2D_xy = ['pr','tas','clt']
var3D_xy = ['ua','va']

levels = {'ua': [925,700,200],
          'va': [925,700,200]
          }

xlim = (-180,180)
ylim = (-90,90)
xticks=range(-180,181,60)
yticks=range(-90,91,30)

continent = 'ideal'
if cont == 'WA':
    xxc = [0,60,60,0,0]
    yyc = [5,5,35,35,5]
elif cont == 'tracmip':
    xxc = [0,60,60,0,0]
    yyc = [-30,-30,30,30,-30]
elif cont == 'no':
    continent = None
    xxc = None
    yyc = None

for s in seasons:
    for var in var2D_xy:
        carte_xy(sim,var,rep0=repout,
                 continent=continent,season=s,
                 xlim=xlim,ylim=ylim,
                 xticks=xticks,yticks=yticks,
                 xxc=xxc,yyc=yyc)

    for var in var3D_xy:
        for lev in levels[var]:
            carte_xy(sim,var,lev=lev,rep0=repout,
                     continent=continent,season=s,
                     xlim=xlim,ylim=ylim,
                     xticks=xticks,yticks=yticks,
                     xxc=xxc,yyc=yyc)

stop
#### Transect continentaux    

seasons = ['JJAS']
lonavg = (0,60)

var2D_y = ['pr','tas','clt']
var3D_y = ['ua','va']

levels = {'ua': [925,700,200],
          'va': [925,700,200]
          }

xlim = (-40,40)
xticks = range(-40,41,10)

continent = 'ideal'
if cont == 'WA':
    clats = (5,35)
elif cont == 'tracmip':    
    clats = (-30,30)

for s in seasons:
    for var in var2D_y:
        transect_y(sim,var,rep0=repout,
                lonavg=lonavg,
                continent=continent,season=s,
                xlim=xlim,ylim=None,
                xticks=xticks,yticks=None,
                clats=clats,addAll=True)

    for var in var3D_y:
        for lev in levels[var]:
            transect_y(sim,var,lev=lev,rep0=repout,
                    lonavg=lonavg,
                    continent=continent,season=s,
                    xlim=xlim,ylim=None,
                    xticks=xticks,yticks=None,
                    clats=clats,addAll=True)
