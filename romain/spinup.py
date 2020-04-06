#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os

import argparse

from domain_avg import domain_avg

repout = './images/'
if not(os.path.exists(repout)):
    os.makedirs(repout)

parser = argparse.ArgumentParser()
parser.add_argument("-sim", help="simulations to plot", type=str, required=True)
args = parser.parse_args()

#sim = 'DELANNOY_slab_cont_WA'
sim = args.sim


domains = ['global','WA']

var2D_xy = ['pr','ts','rlut','rsut','hfls','hfss',]
var3D_xy = ['ta','hus']

levels = {'ta':  [925,700,200,50,10],
          'hus': [925,700,200,50,10]
          }

for d in domains:
    for var in var2D_xy:
        domain_avg(sim,var,rep0=repout,
                   domain=d)

    for var in var3D_xy:
        for lev in levels[var]:
            domain_avg(sim,var,rep0=repout,
                       lev=lev,
                       domain=d)
