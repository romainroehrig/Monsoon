# -*- coding: utf-8 -*-




def getfile(sim,var):

    if sim == 'Simu_Base':
        f = 'C:/Users/Ludovic/Desktop/PFE/Amon_SST/{0}/gr/files/d20180711/{0}_Amon_CNRM-CM6-1_amip_r1i1p1f2_gr_197901-201412.nc'.format(var)
    if sim == 'LC_Sym':
        f = 'C:/Users/Ludovic/Desktop/PFE/Amon/{0}_Amon_CNRM-AM6-DIA-v2_LandControl_r1i1p1_196901-200812.nc'.format(var)
    if sim == 'WA_10a':
        f = 'C:/Users/Ludovic/Desktop/PFE/Simu_WA_10ans/DELANNOY_slab_cont_WA_arpsfx_monthly_{0}_1850-1859.nc'.format(var)
    if sim == 'WA_6a':
        f = 'C:/Users/Ludovic/Desktop/PFE/Simu_WA/DELANNOY_slab_cont_WA_arpsfx_monthly_{0}_1850-1855.nc'.format(var)
    return f
