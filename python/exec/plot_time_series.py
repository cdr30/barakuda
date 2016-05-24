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

venv_needed = {'ORCA','RUN','NN_SST','NN_SSS','NN_SSH','NN_T','NN_S','NN_MLD','LMOCLAT','TRANSPORT_SECTION_FILE','FIG_FORM'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

ff = vdic['FIG_FORM'] ; # format for figures (usually "png" or "svg")

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

elif cdiag == 'mean_fwf':
    venv_ndd = {'NN_FWF','NN_EMP','NN_RNF','NN_P'}
    vdic_fwf = bt.check_env_var(sys.argv[0], venv_ndd)
    idfig = 'fwf'
    cvar  = 'EmPmR'
    clnm  = 'Globally-averaged upward net freshwater flux (E-P-R)'
    cvr2  = 'R'
    cln2  = 'Globally-averaged continental runoffs + ice calving (R)'
    cvr3  = 'EmP'
    cln3  = 'Globally-averaged Evaporation - Precipitation (E-P)'
    cvr4  = 'P'
    cln4  = 'Globally-averaged Precipitation (P)'
    cyu   = r'Sv'
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
        print 'ERROR: plot_time_series.py => '+cvar+', number of records not a multiple of 12!'
        sys.exit(0)

    # Annual data
    VY, FY = bt.monthly_2_annual(vtime[:], vvar[:])

    ittic = bt.iaxe_tick(nbm/12)

    # Time to plot
    bp.plot("1d_mon_ann")(vtime, VY, vvar, FY, cfignm=cdiag+'_'+CONFRUN, dt_year=ittic,
                          cyunit=cyu, ctitle = CONFRUN+': '+clnm, ymin=ym, ymax=yp, cfig_type=ff)




if idfig == 'fwf':

    l_rnf = False ; l_emp = False ; l_prc = False

    cf_in = cdiag+'_'+CONFRUN+'_global.nc' ;  bt.chck4f(cf_in, script_name='plot_time_series.py')
    id_in = Dataset(cf_in)
    list_var = id_in.variables.keys()
    vtime = id_in.variables['time'][:] ; nbm = len(vtime)
    vfwf  = id_in.variables[cvar][:]
    if cvr2 in list_var[:]:
        l_rnf = True
        vrnf  = id_in.variables[cvr2][:]
    if cvr3 in list_var[:]:
        l_emp = True
        vemp  = id_in.variables[cvr3][:]
    if cvr4 in list_var[:]:
        l_prc = True
        vprc  = id_in.variables[cvr3][:]
    id_in.close()

    # Checking if there a potential file for IFS:
    l_fwf_ifs = False
    cf_IFS_in = cdiag+'_IFS_'+vdic['RUN']+'_global.nc'
    print '  *** Checking for the existence of '+cf_IFS_in
    if os.path.exists(cf_IFS_in):
        print "  *** IFS FWF files found!"
        id_IFS_in = Dataset(cf_IFS_in)
        vemp_ifs = id_IFS_in.variables['flx_emp_sv'][:]
        ve_ifs   = id_IFS_in.variables['flx_e_sv'][:]
        vp_ifs   = id_IFS_in.variables['flx_p_sv'][:]
        vemp_glb_ifs = id_IFS_in.variables['flx_emp_glb_sv'][:]
        #ve_glb_ifs   = id_IFS_in.variables['flx_e_glb_sv'][:]
        #vp_glb_ifs   = id_IFS_in.variables['flx_p_glb_sv'][:]
        vemp_land_ifs = id_IFS_in.variables['flx_emp_land_sv'][:]
        ve_land_ifs   = id_IFS_in.variables['flx_e_land_sv'][:]
        vp_land_ifs   = id_IFS_in.variables['flx_p_land_sv'][:]
        id_IFS_in.close()
        if len(vemp_ifs) != nbm:
            print 'ERROR: plot_time_series.py => length of E-P of IFS in '+cf_IFS_in+' does not agree with its NEMO counterpart!'
            print '       =>', len(vemp_ifs), nbm
            sys.exit(0)
        l_fwf_ifs = True
    else:
        print '       => Nope!\n'


    if nbm%12 != 0:
        print 'ERROR: plot_time_series.py => '+cvar+', number of records not a multiple of 12!'
        sys.exit(0)

    ittic = bt.iaxe_tick(nbm/12)

    # Annual data
    VY, FY = bt.monthly_2_annual(vtime, vfwf)
    # Time to plot
    bp.plot("1d_mon_ann")(vtime, VY, vfwf, FY, cfignm=cdiag+'_fwf_'+CONFRUN, dt_year=ittic,
                          cyunit=cyu, ctitle = CONFRUN+': '+clnm, ymin=ym, ymax=yp, cfig_type=ff)

    if l_rnf:
        VY, FY = bt.monthly_2_annual(vtime, vrnf)
        bp.plot("1d_mon_ann")(vtime, VY, vrnf, FY, cfignm=cdiag+'_rnf_'+CONFRUN, dt_year=ittic,
                              cyunit=cyu, ctitle = CONFRUN+': '+cln2, ymin=ym, ymax=yp, cfig_type=ff)

    if l_emp:
        VY, FY = bt.monthly_2_annual(vtime, vemp)
        bp.plot("1d_mon_ann")(vtime, VY, vemp, FY, cfignm=cdiag+'_emp_'+CONFRUN, dt_year=ittic,
                              cyunit=cyu, ctitle = CONFRUN+': '+cln3, ymin=ym, ymax=yp, cfig_type=ff)

    if l_prc:
        VY, FY = bt.monthly_2_annual(vtime, vprc)
        bp.plot("1d_mon_ann")(vtime, VY, vprc, FY, cfignm=cdiag+'_prc_'+CONFRUN, dt_year=ittic,
                              cyunit=cyu, ctitle = CONFRUN+': '+cln4, ymin=ym, ymax=yp, cfig_type=ff)

    if l_fwf_ifs and l_emp:
        # Only E-P for NEMO and IFS:
        Xplt = nmp.zeros((3,nbm))
        Xplt[0,:] = vemp[:]
        Xplt[1,:] = vemp_ifs[:]
        Xplt[2,:] = vemp_land_ifs[:]
        bp.plot("1d_multi")(vtime, Xplt, ['E-P NEMO','E-P IFS (oceans)','E-P IFS (land)'],
                            cfignm=cdiag+'_emp_IFS_'+CONFRUN, dt_year=ittic,
                            cyunit=cyu, ctitle = CONFRUN+': E-P (monthly)', ymin=ym, ymax=yp, cfig_type=ff)

        # Same but annual:
        nby = nbm/12
        Xplt = nmp.zeros((2,nby))
        VY, Xplt[0,:] = bt.monthly_2_annual(vtime[:], vemp[:])
        VY, Xplt[1,:] = bt.monthly_2_annual(vtime[:], vemp_ifs[:])
        bp.plot("1d_multi")(VY, Xplt, ['E-P NEMO','E-P IFS'],
                            cfignm=cdiag+'_emp_IFS_annual_'+CONFRUN, dt_year=ittic,
                            cyunit=cyu, ctitle = CONFRUN+': E-P over oceans (annual)',
                            loc_legend='upper center', ymin=ym, ymax=yp, cfig_type=ff)
        



    if l_fwf_ifs and l_rnf:
        # Runoff of NEMO compares to ( E-P global - E-P ocean ) of IFS: #lulu
        Xplt = nmp.zeros((2,nbm))
        Xplt[0,:] = vrnf[:]
        Xplt[1,:] = -vemp_land_ifs[:]
        bp.plot("1d_multi")(vtime, Xplt, ['R NEMO','-(E-P) over land IFS'], cfignm=cdiag+'_rnf_IFS_'+CONFRUN, dt_year=ittic,
                            cyunit=cyu, ctitle = CONFRUN+': Continental runoffs (monthly)', loc_legend='upper center',
                            ymin=ym, ymax=yp, cfig_type=ff)

        # Same but annual:
        nby = nbm/12
        Xplt = nmp.zeros((2,nby))
        VY, Xplt[0,:] = bt.monthly_2_annual(vtime[:], vrnf[:])
        VY, Xplt[1,:] = bt.monthly_2_annual(vtime[:], -vemp_land_ifs[:])
        #VY, Xplt[2,:] = bt.monthly_2_annual(vtime[:], -(vemp_glb_ifs[:]-vemp_ifs[:])) ; # ,'-(Global(E-P)-Oceans(E-P) IFS'
        bp.plot("1d_multi")(VY, Xplt, ['R NEMO','-(E-P) over land IFS'], cfignm=cdiag+'_rnf_IFS_annual_'+CONFRUN, dt_year=ittic,
                            cyunit=cyu, ctitle = CONFRUN+': Continental runoffs (annual)', loc_legend='upper center',
                            ymin=ym, ymax=yp, cfig_type=ff)

    if l_fwf_ifs and l_prc:
        # Only P for NEMO and IFS:
        Xplt = nmp.zeros((3,nbm))
        Xplt[0,:] = vprc[:]
        Xplt[1,:] = vp_ifs[:]
        Xplt[2,:] = vp_land_ifs[:]
        bp.plot("1d_multi")(vtime, Xplt, ['P NEMO','P IFS (oceans)','P IFS (land)'], cfignm=cdiag+'_prc_IFS_'+CONFRUN, dt_year=ittic,
                            cyunit=cyu, ctitle = CONFRUN+': Precip', ymin=ym, ymax=yp, cfig_type=ff)


        # Everything possible
        Xplt = nmp.zeros((7,nbm))
        vlab = []
        Xplt[0,:] = vemp[:]      ; vlab.append('E-P NEMO ('+vdic_fwf['NN_EMP']+')')
        Xplt[1,:] = vemp_ifs[:]  ; vlab.append('E-P IFS')
        Xplt[2,:] = vfwf[:]      ; vlab.append('E-P-R NEMO ('+vdic_fwf['NN_FWF']+')')
        Xplt[3,:] = vrnf[:]      ; vlab.append('R NEMO ('+vdic_fwf['NN_RNF']+')')
        Xplt[4,:] = ve_ifs[:]    ; vlab.append('E IFS')
        Xplt[5,:] = vprc[:]      ; vlab.append('P NEMO ('+vdic_fwf['NN_P']+')')
        Xplt[6,:] = vp_ifs[:]    ; vlab.append('P IFS')


        bp.plot("1d_multi")(vtime, Xplt, vlab,
                            cfignm=cdiag+'_emp_ALL_IFS_'+CONFRUN, dt_year=ittic,
                            loc_legend='center', cyunit=cyu,
                            ctitle = CONFRUN+': fresh-water budgets', ymin=ym, ymax=yp, cfig_type=ff)

        #lulu



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
                FM = nmp.zeros((nb_oce, nbzrange, nbm))
            print '   * reading '+cvar+'_'+czr+' in '+cf_in
            FM[joce,jz,:]  = id_in.variables[cvar+'_'+czr][:]
            jz = jz + 1
        id_in.close()

        # Annual data:
        if not joce:
            nby = nbm/12
            FY = nmp.zeros((nb_oce, 4, nby))
        VY, FY[joce,:,:] = bt.monthly_2_annual(vtime[:], FM[joce,:,:])

        print ' *** '+coce+' done...\n'
        joce = joce + 1

    ittic = bt.iaxe_tick(nby)

    # One plot only for global:
    bp.plot("1d_mon_ann")(vtime, VY, FM[0,0,:], FY[0,0,:], cfignm=cdiag+'_'+CONFRUN, dt_year=ittic,
                          cyunit=cyu, ctitle = CONFRUN+': '+clnm, ymin=ym, ymax=yp, cfig_type=ff)


    # Global for different depth:
    bp.plot("1d_multi")(vtime, FM[0,:,:], vlab[:], cfignm=cdiag+'_lev_'+CONFRUN, dt_year=ittic,
                        cyunit=cyu, ctitle = CONFRUN+': '+clnm, ymin=ym0, ymax=yp0, cfig_type=ff)


    # Show each ocean (All depth):
    bp.plot("1d_multi")(vtime, FM[:,0,:], bo.voce2treat, cfignm=cdiag+'_basins_'+CONFRUN, dt_year=ittic,
                        cyunit=cyu, ctitle = CONFRUN+': '+clnm, ymin=ym0, ymax=yp0, cfig_type=ff)






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
            Xamoc   = nmp.zeros((nblat , nbm))
        vlabels[jl] = clat_info
        Xamoc[jl,:] = id_in.variables['moc_atl'][:]
        id_in.close()

        jl = jl + 1

    if nbm%12 != 0:
        print 'ERROR: plot_time_series.py => '+cdiag+', number of records not a multiple of 12!'
        sys.exit(0)
    VY, FY = bt.monthly_2_annual(vtime, Xamoc[i45,:])

    ittic = bt.iaxe_tick(nbm/12)

    # Time to plot
    bp.plot("1d_mon_ann")(vtime, VY, Xamoc[i45,:], FY, cfignm=cdiag+'_'+CONFRUN, dt_year=ittic,
                          cyunit=cyu, ctitle = CONFRUN+': '+r'Max. of AMOC between '+vlabels[i45],
                          ymin=ym, ymax=yp, dy=1., i_y_jump=2, cfig_type=ff)

    # Annual:
    VY, FY  = bt.monthly_2_annual(vtime, Xamoc[:,:])

    # Time to plot
    bp.plot("1d_multi")(VY, FY, vlabels, cfignm=cdiag+'_'+CONFRUN+'_comp', dt_year=ittic,
                        cyunit=cyu, ctitle = CONFRUN+': '+r'Max. of AMOC', ymin=0, ymax=0,
                        loc_legend='lower left', cfig_type=ff)





if idfig == 'ice':

    vlab_sum = [ 'Arctic (Sept.)'   , 'Antarctic (March)' ]
    vlab_win = [ 'Arctic (March)'   , 'Antarctic (Sept.)' ]

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

    if nbm%12 != 0:
        print 'ERROR: plot_time_series.py => '+cdiag+', number of records not a multiple of 12!'
        sys.exit(0)
    nby = nbm/12

    ittic = bt.iaxe_tick(nby)

    vtime_y = nmp.zeros(nby)
    Xplt = nmp.zeros((2 , nby))

    vtime_y, FY = bt.monthly_2_annual(vtime[:], vvolu_n[:])

    # End local summer
    Xplt[0,:] = varea_n[8::12] ; # extent Arctic september
    Xplt[1,:] = varea_s[2::12] ; # extent Antarctic march
    bp.plot("1d_multi")(vtime_y, Xplt, vlab_sum, cfignm='seaice_extent_summer_'+CONFRUN, dt_year=ittic,
                        cyunit=cyua, ctitle = CONFRUN+': '+r'Sea-Ice extent, end of local summer', ymin=0., ymax=0., cfig_type=ff)

    Xplt[0,:] = vvolu_n[8::12] ; # volume Arctic september
    Xplt[1,:] = vvolu_s[2::12] ; # volume Antarctic march
    bp.plot("1d_multi")(vtime_y, Xplt, vlab_sum, cfignm='seaice_volume_summer_'+CONFRUN, dt_year=ittic,
                        cyunit=cyuv, ctitle = CONFRUN+': '+r'Sea-Ice volume, end of local summer', ymin=0., ymax=0., cfig_type=ff)

    # End of local winter
    Xplt[0,:] = varea_n[2::12] ; # extent Arctic march
    Xplt[1,:] = varea_s[8::12] ; # extent Antarctic september
    bp.plot("1d_multi")(vtime_y, Xplt, vlab_win, cfignm='seaice_extent_winter_'+CONFRUN, dt_year=ittic,
                        cyunit=cyua, ctitle = CONFRUN+': '+r'Sea-Ice extent, end of local winter', ymin=0., ymax=0., cfig_type=ff)

    Xplt[0,:] = vvolu_n[2::12] ; # volume Arctic march
    Xplt[1,:] = vvolu_s[8::12] ; # volume Antarctic september
    bp.plot("1d_multi")(vtime_y, Xplt, vlab_win, cfignm='seaice_volume_winter_'+CONFRUN, dt_year=ittic,
                        cyunit=cyuv, ctitle = CONFRUN+': '+r'Sea-Ice volume, end of local winter', ymin=0., ymax=0., cfig_type=ff)




if idfig == 'transport':

    js = 0
    for csec in list_sections:

        print ' * treating section '+csec

        cf_in = 'transport_sect_'+csec+'.nc' ;   bt.chck4f(cf_in, script_name='plot_time_series.py')
        id_in = Dataset(cf_in)
        if js == 0:
            vtime = id_in.variables['time'][:]
            nbm = len(vtime)
        Xtrsp   = nmp.zeros((3 , nbm)) ; # time + 3 types of transport
        Xtrsp[0,:] = id_in.variables['trsp_volu'][:]
        Xtrsp[1,:] = id_in.variables['trsp_heat'][:]
        Xtrsp[2,:] = id_in.variables['trsp_salt'][:]
        id_in.close()


        if nbm%12 != 0: print 'ERROR: plot_time_series.py => '+cdiag+', number of records not a multiple of 12!', sys.exit(0)
        VY, FY  = bt.monthly_2_annual(vtime, Xtrsp[:,:])

        ittic = bt.iaxe_tick(nbm/12)

        # Transport of volume:
        bp.plot("1d_mon_ann")(vtime, VY, Xtrsp[0,:], FY[0,:], cfignm='transport_vol_'+csec+'_'+CONFRUN,
                             dt_year=ittic, cyunit='Sv', ctitle = CONFRUN+': transport of volume, '+csec,
                             ymin=0, ymax=0, cfig_type=ff)

        # Transport of heat:
        bp.plot("1d_mon_ann")(vtime, VY, Xtrsp[1,:], FY[1,:], cfignm='transport_heat_'+csec+'_'+CONFRUN,
                             dt_year=ittic, cyunit='PW', ctitle = CONFRUN+': transport of heat, '+csec,
                             ymin=0, ymax=0, mnth_col='g', cfig_type=ff)


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
                                  ctitle = CONFRUN+': '+clnm+bo.clgnm_mld_boxes[jbox], ymin=ym, ymax=yp,
                                  plt_m03=True, plt_m09=True, cfig_type=ff)
        else:
            print 'WARNING: plot_time_series.py => MLD diag => '+cf_in_m+' not found!'
        jbox = jbox+1


print 'plot_time_series.py done...\n'
