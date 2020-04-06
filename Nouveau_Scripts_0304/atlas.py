# -*- coding: utf-8 -*-

import os

from carte_xy import carte_xy

repout = './images/'
if not(os.path.exists(repout)):
    os.makedirs(repout)

sim = 'WA_10a'


#### Cartes global

seasons = ['year','JJAS']

var2D_xy = ['pr','tas','clt','hfss','hfls']
var3D_xy = ['ua','va']


levels = {'ua': [925,700,200],
          'va': [925,200]
          }



xlim = (-180,180)
ylim = (-90,90)
xticks=range(-180,181,60)
yticks=range(-90,91,30)

continent = 'ideal'
xxc = [0,60,60,0,0]
yyc = [5,5,35,35,5]

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


