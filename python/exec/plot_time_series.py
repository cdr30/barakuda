#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate misc. time-series out of NEMO output files...
#
#       L. Brodeau, 2013
#

import sys
import os
import numpy as nmp

from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_ncio as bn
import barakuda_orca as bo
import barakuda_plot as bp

DEFAULT_LEGEND_LOC = 'lower left'

venv_needed = {'ORCA','RUN','NN_SST','NN_SSS','NN_SSH','NN_T','NN_S','NN_MLD','LMOCLAT','TRANSPORT_SECTION_FILE'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

narg = len(sys.argv)
if narg != 2:
    print 'Usage: {} <diag>'.format(sys.argv[0])
    sys.exit(0)
cdiag = sys.argv[1]

print '\n plot_time_series.py: diag => "'+cdiag+'"'

if cdiag == 'mean_tos':
    cvar  = vdic['NN_SST']
    idfig = 'simple'
    clnm  = 'Globally-averaged sea surface temperature'
    cyu   = r'$^{\circ}$C'
    ym    = yp = 0.

elif cdiag == 'mean_sos':
    cvar  = vdic['NN_SSS']
    idfig = 'simple'
    clnm  = 'Globally-averaged sea surface salinity'
    cyu   = r'PSU'
    ym = yp = 0.

elif cdiag == 'mean_zos':
    cvar  = vdic['NN_SSH']
    idfig = 'simple'
    clnm  = 'Globally-averaged sea surface height'
    cyu   = r'm'
    ym = yp = 0.


elif  cdiag == '3d_thetao':
    cvar  = vdic['NN_T']
    idfig = 'ts3d'
    clnm = 'Globally-averaged temperature'
    cyu  = r'$^{\circ}$C'
    #ym = 3.6 ; yp = 4.
    ym = 0. ; yp = 0.
    #ym0  = 1.5 ; yp0 = 20.
    ym0  = yp0 = 0.

elif cdiag == '3d_so':
    cvar  = vdic['NN_S']
    idfig = 'ts3d'
    clnm = 'Globally-averaged salinity'
    cyu  = r'PSU'
    #ym  = 34.6 ; yp  = 35.
    #ym0 = 34.6 ; yp0 = 35.
    ym  = yp  = 0.
    ym0 = yp0 = 0.

elif cdiag == 'amoc':
    idfig = 'amoc'
    cyu  = r'Sv'
    ym = 3.5
    yp = 24.5

elif cdiag == 'mean_mldr10_1':
    cvar  = vdic['NN_MLD']
    idfig = 'mld'
    clnm  = 'Mean mixed-layer depth, '
    cyu   = r'm'
    ym = yp = 0.


elif cdiag == 'transport_sections':
    idfig = 'transport'
    print '  Using TRANSPORT_SECTION_FILE = '+vdic['TRANSPORT_SECTION_FILE']
    list_sections = bo.get_sections_names_from_file(vdic['TRANSPORT_SECTION_FILE'])
    print 'List of sections to treat: ', list_sections



elif cdiag == 'seaice':
    idfig = 'ice'
    cyu  = r'10$^6$km$^2$'



else:
    print 'ERROR: plot_time_series.py => diagnostic '+cdiag+' unknown!'; sys.exit(0)









##########################################
# Basic temp., sali. and SSH time series #
##########################################

if idfig == 'simple':

    cf_in = 'mean_'+cvar+'_'+CONFRUN+'_global.nc' ;  bt.chck4f(cf_in, script_name='plot_time_series.py')
    id_in = Dataset(cf_in)
    vtime = id_in.variables['time'][:] ; nbm = len(vtime)
    vvar  = id_in.variables[cvar][:]
    id_in.close()

    if nbm%12 != 0:
        print 'ERROR: plot_time_series.py => '+cvar+', numberof records not a multiple of 12!'
        sys.exit(0)

    # Annual data
    VY, FY = bt.monthly_2_annual(vtime[:], vvar[:])

    ittic = bt.iaxe_tick(nbm/12)

    # Time to plot
    bp.plot("1d_mon_ann")(vtime, VY, vvar, FY, cfignm=cdiag+'_'+CONFRUN, dt_year=ittic,
                          cyunit=cyu, ctitle = CONFRUN+': '+clnm, ymin=ym, ymax=yp)




if idfig == 'ts3d':

    nb_oce = len(bo.voce2treat)

    vzrange = [ '0-bottom', '0-100'  , '100-1000',   '1000-bottom'  ] ;  nbzrange = len(vzrange)
    vlab    = [ 'AllDepth', '0m-100m', '100m-1000m', '1000m-bottom' ]

    joce = 0
    for coce in bo.voce2treat[:]:

        cf_in = '3d_'+cvar+'_'+CONFRUN+'_'+coce+'.nc' ;  bt.chck4f(cf_in, script_name='plot_time_series.py')
        id_in = Dataset(cf_in)
        vtime = id_in.variables['time'][:] ; nbm = len(vtime)
        jz = 0
        for czr in vzrange:
            if not joce and not jz:
                FM = nmp.zeros(nbm*nbzrange*nb_oce)
                FM.shape = [ nb_oce, nbzrange, nbm ]
            print '   * reading '+cvar+'_'+czr+' in '+cf_in
            FM[joce,jz,:]  = id_in.variables[cvar+'_'+czr][:]
            jz = jz + 1
        id_in.close()

        # Annual data:
        if not joce:
            nby = nbm/12
            FY = nmp.zeros(nby*4*nb_oce) ; FY.shape = [ nb_oce, 4, nby ]
        VY, FY[joce,:,:] = bt.monthly_2_annual(vtime[:], FM[joce,:,:])

        print ' *** '+coce+' done...\n'
        joce = joce + 1

    ittic = bt.iaxe_tick(nby)

    # One plot only for global:
    bp.plot("1d_mon_ann")(vtime, VY, FM[0,0,:], FY[0,0,:], cfignm=cdiag+'_'+CONFRUN, dt_year=ittic,
                          cyunit=cyu, ctitle = CONFRUN+': '+clnm, ymin=ym, ymax=yp)


    # Global for different depth:
    bp.plot("1d_multi")(vtime, FM[0,:,:], vlab[:], cfignm=cdiag+'_lev_'+CONFRUN, dt_year=ittic,
                        cyunit=cyu, ctitle = CONFRUN+': '+clnm, ymin=ym0, ymax=yp0)


    # Show each ocean (All depth):
    bp.plot("1d_multi")(vtime, FM[:,0,:], bo.voce2treat, cfignm=cdiag+'_basins_'+CONFRUN, dt_year=ittic,
                        cyunit=cyu, ctitle = CONFRUN+': '+clnm, ymin=ym0, ymax=yp0)






##########################################
# AMOC
##########################################

if idfig == 'amoc':
    clmoc = vdic['LMOCLAT']
    list_lat = clmoc.split() ; nblat = len(list_lat)
    print '\n AMOC: '+str(nblat)+' latitude bands!'

    i45 = 3 ; # position of AMOC at 45!

    jl = 0
    for clr in list_lat:
        [ c1, c2 ] = clr.split('-') ; clat_info = '+'+c1+'N+'+c2+'N'
        cf_in = 'max_moc_atl_'+clat_info+'.nc' ; bt.chck4f(cf_in, script_name='plot_time_series.py')
        id_in = Dataset(cf_in)
        if not jl:
            vtime = id_in.variables['time'][:] ; nbm = len(vtime)
            vlabels = nmp.zeros(nblat, dtype = nmp.dtype('a8'))
            Xamoc   = nmp.zeros(nbm*(nblat)) ; Xamoc.shape = [ nblat , nbm ]
        vlabels[jl] = clat_info
        Xamoc[jl,:] = id_in.variables['moc_atl'][:]
        id_in.close()

        jl = jl + 1

    if nbm%12 != 0:
        print 'ERROR: plot_time_series.py => '+cdiag+', numberof records not a multiple of 12!'
        sys.exit(0)
    VY, FY = bt.monthly_2_annual(vtime, Xamoc[i45,:])

    ittic = bt.iaxe_tick(nbm/12)

    # Time to plot
    bp.plot("1d_mon_ann")(vtime, VY, Xamoc[i45,:], FY, cfignm=cdiag+'_'+CONFRUN, dt_year=ittic,
                          cyunit=cyu, ctitle = CONFRUN+': '+r'Max. of AMOC between '+vlabels[i45],
                          ymin=ym, ymax=yp, dy=1., i_y_jump=2)

    # Annual:
    VY, FY  = bt.monthly_2_annual(vtime, Xamoc[:,:])

    # Time to plot
    bp.plot("1d_multi")(VY, FY, vlabels, cfignm=cdiag+'_'+CONFRUN+'_comp', dt_year=ittic,
                        cyunit=cyu, ctitle = CONFRUN+': '+r'Max. of AMOC', ymin=0, ymax=0,
                        loc_legend='lower left')





if idfig == 'ice':

    vlab = [ 'Arctic', 'Antarctic' ]

    # montly sea-ice volume and extent, Arctic and Antarctic...
    cf_in = 'seaice_diags.nc' ;  bt.chck4f(cf_in, script_name='plot_time_series.py')
    id_in = Dataset(cf_in)
    vtime = id_in.variables['time'][:] ; nbm = len(vtime)
    vvolu_n  = id_in.variables['volu_ne'][:]
    varea_n  = id_in.variables['area_ne'][:]
    vvolu_s  = id_in.variables['volu_se'][:]
    varea_s  = id_in.variables['area_se'][:]
    id_in.close()

    cyua = r'10$^6$km$^2$'
    cyuv = r'10$^3$km$^3$'

<<<<<<< HEAD
    if nbm%12 != 0: print 'ERROR: plot_time_series.py => '+cdiag+', numberof records not a multiple of 12!', sys.exit(0)
=======
    if nbm%12 != 0:
        print 'ERROR: plot_time_series.py => '+cdiag+', numberof records not a multiple of 12!'
        sys.exit(0)
>>>>>>> master
    nby = nbm/12

    ittic = bt.iaxe_tick(nby)

    vtime_y = nmp.zeros(nby)
    Xplt = nmp.zeros(2*nby) ; Xplt.shape = [2 , nby]

    vtime_y, FY = bt.monthly_2_annual(vtime[:], vvolu_n[:])

    # End local summer
    Xplt[0,:] = varea_n[8::12] ; # extent Arctic september
    Xplt[1,:] = varea_s[2::12] ; # extent Antarctic march
<<<<<<< HEAD
    bp.plot_1d_multi(vtime_y, Xplt, vlab, cfignm='seaice_extent_summer_'+CONFRUN, dt_year=ittic,
                     cyunit=cyua, ctitle = CONFRUN+': '+r'Sea-Ice extent, end of local summer', ymin=0., ymax=0.)

    Xplt[0,:] = vvolu_n[8::12] ; # volume Arctic september
    Xplt[1,:] = vvolu_s[2::12] ; # volume Antarctic march
    bp.plot_1d_multi(vtime_y, Xplt, vlab, cfignm='seaice_volume_summer_'+CONFRUN, dt_year=ittic,
                     cyunit=cyuv, ctitle = CONFRUN+': '+r'Sea-Ice volume, end of local summer', ymin=0., ymax=0.)
=======
    bp.plot("1d_multi")(vtime_y, Xplt, vlab, cfignm='seaice_extent_summer_'+CONFRUN, dt_year=ittic,
                        cyunit=cyua, ctitle = CONFRUN+': '+r'Sea-Ice extent, end of local summer', ymin=0., ymax=0.)

    Xplt[0,:] = vvolu_n[8::12] ; # volume Arctic september
    Xplt[1,:] = vvolu_s[2::12] ; # volume Antarctic march
    bp.plot("1d_multi")(vtime_y, Xplt, vlab, cfignm='seaice_volume_summer_'+CONFRUN, dt_year=ittic,
                        cyunit=cyuv, ctitle = CONFRUN+': '+r'Sea-Ice volume, end of local summer', ymin=0., ymax=0.)
>>>>>>> master

    # End of local winter
    Xplt[0,:] = varea_n[2::12] ; # extent Arctic march
    Xplt[1,:] = varea_s[8::12] ; # extent Antarctic september
<<<<<<< HEAD
    bp.plot_1d_multi(vtime_y, Xplt, vlab, cfignm='seaice_extent_winter_'+CONFRUN, dt_year=ittic,
                     cyunit=cyua, ctitle = CONFRUN+': '+r'Sea-Ice extent, end of local winter', ymin=0., ymax=0.)

    Xplt[0,:] = vvolu_n[2::12] ; # volume Arctic march
    Xplt[1,:] = vvolu_s[8::12] ; # volume Antarctic september
    bp.plot_1d_multi(vtime_y, Xplt, vlab, cfignm='seaice_volume_winter_'+CONFRUN, dt_year=ittic,
                     cyunit=cyuv, ctitle = CONFRUN+': '+r'Sea-Ice volume, end of local winter', ymin=0., ymax=0.)
=======
    bp.plot("1d_multi")(vtime_y, Xplt, vlab, cfignm='seaice_extent_winter_'+CONFRUN, dt_year=ittic,
                        cyunit=cyua, ctitle = CONFRUN+': '+r'Sea-Ice extent, end of local winter', ymin=0., ymax=0.)

    Xplt[0,:] = vvolu_n[2::12] ; # volume Arctic march
    Xplt[1,:] = vvolu_s[8::12] ; # volume Antarctic september
    bp.plot("1d_multi")(vtime_y, Xplt, vlab, cfignm='seaice_volume_winter_'+CONFRUN, dt_year=ittic,
                        cyunit=cyuv, ctitle = CONFRUN+': '+r'Sea-Ice volume, end of local winter', ymin=0., ymax=0.)
>>>>>>> master




if idfig == 'transport':

    js = 0
    for csec in list_sections:

        print ' * treating section '+csec

        cf_in = 'transport_sect_'+csec+'.nc' ;   bt.chck4f(cf_in, script_name='plot_time_series.py')
        id_in = Dataset(cf_in)
        if js == 0:
            vtime = id_in.variables['time'][:]
            nbm = len(vtime)
        Xtrsp   = nmp.zeros(nbm*3) ; Xtrsp.shape = [ 3 , nbm ] ; # time + 3 types of transport
        Xtrsp[0,:] = id_in.variables['trsp_volu'][:]
        Xtrsp[1,:] = id_in.variables['trsp_heat'][:]
        Xtrsp[2,:] = id_in.variables['trsp_salt'][:]
        id_in.close()


        if nbm%12 != 0: print 'ERROR: plot_time_series.py => '+cdiag+', numberof records not a multiple of 12!', sys.exit(0)
        VY, FY  = bt.monthly_2_annual(vtime, Xtrsp[:,:])

        ittic = bt.iaxe_tick(nbm/12)

        # Transport of volume:
        bp.plot("1d_mon_ann")(vtime, VY, Xtrsp[0,:], FY[0,:], cfignm='transport_vol_'+csec+'_'+CONFRUN,
                             dt_year=ittic, cyunit='Sv', ctitle = CONFRUN+': transport of volume, '+csec,
                             ymin=0, ymax=0)

        # Transport of heat:
        bp.plot("1d_mon_ann")(vtime, VY, Xtrsp[1,:], FY[1,:], cfignm='transport_heat_'+csec+'_'+CONFRUN,
                             dt_year=ittic, cyunit='PW', ctitle = CONFRUN+': transport of heat, '+csec,
                             ymin=0, ymax=0, mnth_col='g')


        js = js + 1




if idfig == 'mld':
    jbox = 0
    for cbox in bo.cname_mld_boxes:
        cf_in_m = 'mean_'+cvar+'_'+CONFRUN+'_'+cbox+'.nc'
        if os.path.exists(cf_in_m):
            print ' Opening '+cf_in_m
            vt0, vd0 = bn.read_1d_series(cf_in_m, cvar, cv_t='time', l_return_time=True)
            nbm = len(vt0)
            if nbm%12 != 0:
                print 'ERROR: plot_time_series.py => '+cvar+', number of records not a multiple of 12!'
                sys.exit(0)
            VY, FY = bt.monthly_2_annual(vt0, vd0)
            ittic = bt.iaxe_tick(nbm/12)
            bp.plot("1d_mon_ann")(vt0, VY, vd0, FY, cfignm=cdiag+'_'+CONFRUN+'_'+cbox, dt_year=ittic, cyunit=cyu,
                                  ctitle = CONFRUN+': '+clnm+bo.clgnm_mld_boxes[jbox], ymin=ym, ymax=yp, plt_m03=True, plt_m09=True)
        else:
            print 'WARNING: plot_time_series.py => MLD diag => '+cf_in_m+' not found!'
        jbox = jbox+1


print 'plot_time_series.py done...\n'
