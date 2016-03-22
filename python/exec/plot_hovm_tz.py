#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate time-depth Hovmuller diagrams of 3D fields out of NEMO output files...
#
# L. Brodeau, 2011
#

import sys
import os
import numpy as nmp

from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_orca as bo
import barakuda_plot as bp


venv_needed = {'ORCA','RUN','DIAG_D','NN_T','NN_S'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

cname_temp = vdic['NN_T']
cname_sali = vdic['NN_S']



#if len(sys.argv) != 4 and len(sys.argv) != 6 :
#    print 'Usage: '+sys.argv[0]+' <YYYY1> <YYYY2> <Nb. Levels> (<name temp.> <name salin.>)'
#    sys.exit(0)

#cy1 = sys.argv[1] ; cy2 = sys.argv[2] ; jy1=int(cy1); jy2=int(cy2)



path_fig=vdic['DIAG_D']+'/'

fig_type='png'

voceans_u = [cc.upper() for cc in bo.voce2treat]  ; # same list but in uppercase


jo = 0
for coce in bo.voce2treat:

    cf_temp = cname_temp+'_mean_Vprofile_'+CONFRUN+'_'+coce+'.nc' ; bt.chck4f(cf_temp)
    cf_sali = cname_sali+'_mean_Vprofile_'+CONFRUN+'_'+coce+'.nc' ; bt.chck4f(cf_sali)

    id_temp = Dataset(cf_temp)
    if jo == 0:
        vyears = id_temp.variables['time'][:]
        vdepth = id_temp.variables['deptht'][:]
    XT = id_temp.variables[cname_temp][:,:]
    id_temp.close()

    id_sali = Dataset(cf_sali)
    XS = id_sali.variables[cname_sali][:,:]
    id_sali.close()


    if jo == 0:
        vyears = nmp.trunc(vyears) + 0.5 ; # in case 1990 and not 1990.5 !!!
        yr1=float(int(min(vyears)))
        yr2=float(int(max(vyears)))


    [nby, nz] = XT.shape

    ixtics = bt.iaxe_tick(nby)

    # Number of NaN vertical points:
    visnan = nmp.isnan(XT[0,:])
    nz_nan = nmp.sum(visnan)
    nz = nz - nz_nan

    XTe = nmp.zeros((nz,nby))
    XTe[:,:] = nmp.flipud(nmp.rot90(XT[:,:nz]))

    XSe = nmp.zeros((nz,nby))
    XSe[:,:] = nmp.flipud(nmp.rot90(XS[:,:nz]))



    # Removing value for first year to all years:
    vy1 = nmp.zeros(nz) ; vy1[:] = XTe[:,0]
    for jy in range(nby): XTe[:,jy] = XTe[:,jy] - vy1[:]
    vy1 = nmp.zeros(nz) ; vy1[:] = XSe[:,0]
    for jy in range(nby): XSe[:,jy] = XSe[:,jy] - vy1[:]

    z0 = vdepth[0]
    zK = max(vdepth)

    [ rmin, rmax, rdf ] = bt.get_min_max_df(XTe,40)
    bp.plot("vert_section")(vyears[:], vdepth[:nz], XTe[:,:], XTe[:,:]*0.+1., rmin, rmax, rdf,
                            cpal='bbr2', xmin=yr1, xmax=yr2+1., dx=ixtics, lkcont=False,
                            zmin = z0, zmax = zK, l_zlog=True,
                            cfignm=path_fig+'hov_temperature_'+CONFRUN+'_'+coce, cbunit=r'$^{\circ}$C', cxunit='',
                            czunit='Depth (m)',
                            ctitle=CONFRUN+': Spatially-averaged temperature evolution, '+voceans_u[jo]+', ('+str(int(yr1))+'-'+str(int(yr2))+')',
                            cfig_type=fig_type, lforce_lim=True, i_cb_subsamp=2)

    XSe = 1000.*XSe
    [ rmin, rmax, rdf ] = bt.get_min_max_df(XSe,40)
    bp.plot("vert_section")(vyears[:], vdepth[:nz], XSe[:,:], XSe[:,:]*0.+1., rmin, rmax, rdf,
                            cpal='bbr2', xmin=yr1, xmax=yr2+1., dx=ixtics, lkcont=False,
                            zmin = z0, zmax = zK, l_zlog=True,
                            cfignm=path_fig+'hov_salinity_'+CONFRUN+'_'+coce, cbunit=r'10$^{-3}$PSU', cxunit='',
                            czunit='Depth (m)',
                            ctitle=CONFRUN+': Spatially-averaged salinity evolution, '+voceans_u[jo]+', ('+str(int(yr1))+'-'+str(int(yr2))+')',
                            cfig_type=fig_type, lforce_lim=True, i_cb_subsamp=2)


    jo = jo +1
