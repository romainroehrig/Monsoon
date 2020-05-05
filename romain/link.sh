#!/bin/sh

DIR0=`pwd`

sim_pattern=DELANNOY_slab_cont_WA_qflux_evap
sim_pattern=DELANNOY_slab_cont_tracmip_evap

var=pr

rm -rf tmp
mkdir tmp

for i in `seq -w 00 10`
do
  ln -s $DIR0/images//${sim_pattern}$i/JJASavg//xy_$var.png $DIR0/tmp/xy_${var}_evap$i.png
done
