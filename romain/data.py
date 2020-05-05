# -*- coding: utf-8 -*-

rep0 = '/Volumes/CNRM/TMP/simulations/DELANNOY/'


def getfile(sim,var,dtype = None):

    if sim == 'CNRM-CM6-1':
        f = 'C:/Users/Ludovic/Desktop/PFE/Amon_SST/{0}/gr/files/d20180711/{0}_Amon_CNRM-CM6-1_amip_r1i1p1f2_gr_197901-201412.nc'.format(var)
    elif sim in ['DELANNOY_slab_cont_WA']:
        f = '{0}/{1}/M/{1}_arpsfx_monthly_{2}_1850-1859.nc'.format(rep0,sim,var)

    elif sim in ['DELANNOY_slab_cont_WA_evap{0:0>2}'.format(i) for i in range(0,11)] + ['DELANNOY_slab_ref',]:
        f = '{0}/{1}/concat/{1}_arpsfx_monthly_{2}_1850-1869.nc'.format(rep0,sim,var)
        if not(dtype is None):
            f = '{0}/{1}/{2}/{3}.nc'.format(rep0,sim,dtype,var)

    elif sim in ['DELANNOY_slab_cont_WA_qflux_evap{0:0>2}'.format(i) for i in range(0,11)]\
              + ['DELANNOY_slab_cont_tracmip_evap{0:0>2}'.format(i) for i in range(0,11)]\
              + ['DELANNOY_slab_ref_qflux',]:
        f = '{0}/{1}/M/{1}_arpsfx_monthly_{2}_1850-1859.nc'.format(rep0,sim,var)
        if not(dtype is None):
            f = '{0}/{1}/M/{2}/{3}.nc'.format(rep0,sim,dtype,var)

    return f
