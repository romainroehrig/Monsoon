#!/bin/sh

DIR0=`pwd`


var=pr

rm -rf tmp
mkdir tmp

for i in `seq -w 00 10`
do
  ln -s $DIR0/images//DELANNOY_slab_cont_WA_evap$i/JJASavg//xy_$var.png $DIR0/tmp/xy_${var}_evap$i.png
done
