#!/usr/bin/env python

# L. Brodeau, november 2013

import sys
import os
import numpy as nmp

from netCDF4 import Dataset

import barakuda_ncio as bn
import barakuda_orca as bo
import barakuda_plot as bp
import barakuda_tool as bt

DEFAULT_LEGEND_LOC = 'lower right'

iamoc  = 1
i2dfl  = 1
i3dfl  = 1
imld   = 1
iice   = 1
itrsp = 1


LIST_RUNS = os.getenv('LIST_RUNS')
if LIST_RUNS == None: print 'The LIST_RUNS environement variable is no set'; sys.exit(0)

DIAG_DIR = os.getenv('DIAG_DIR')
if DIAG_DIR == None: print 'The DIAG_DIR environement variable is no set'; sys.exit(0)

CONF = os.getenv('CONF')
if CONF == None: print 'The CONF environement variable is no set'; sys.exit(0)

FIG_FORMAT = os.getenv('FIG_FORMAT')
if FIG_FORMAT == None: print 'The FIG_FORMAT environement variable is no set'; sys.exit(0)

NN_SST = os.getenv('NN_SST')
if NN_SST == None: print 'The NN_SST environement variable is no set'; sys.exit(0)
NN_SSS = os.getenv('NN_SSS')
if NN_SSS == None: print 'The NN_SSS environement variable is no set'; sys.exit(0)
NN_SSH = os.getenv('NN_SSH')
if NN_SSH == None: print 'The NN_SSH environement variable is no set'; sys.exit(0)
NN_T = os.getenv('NN_T')
if NN_T == None: print 'The NN_T environement variable is no set'; sys.exit(0)
NN_S = os.getenv('NN_S')
if NN_S == None: print 'The NN_S environement variable is no set'; sys.exit(0)
NN_MLD = os.getenv('NN_MLD')
if NN_MLD == None: print 'The NN_MLD environement variable is no set'; sys.exit(0)



narg = len(sys.argv)
if narg != 3: print 'Usage: '+sys.argv[0]+' <first_year> <last_year>'; sys.exit(0)
cy1 = sys.argv[1] ; y1 = int(cy1)
cy2 = sys.argv[2] ; y2 = int(cy2)

nb_years = y2 - y1 + 1


clist_runs = LIST_RUNS.split()
clist_confruns = []

for crun in clist_runs:
    clist_confruns.append(CONF+'-'+crun)

print sys.argv[0]+': will compare following runs: '; print clist_confruns
print ' ... saved into '+DIAG_DIR+'\n'

nbrun = len(clist_confruns)



ittic = bt.iaxe_tick(nb_years)

vtime = nmp.zeros(nb_years)
Xf = nmp.zeros(nbrun*nb_years) ; Xf.shape = [ nbrun, nb_years ]












# Only one column to read:
##########################

if i2dfl == 1:

    vvar  = [  NN_SSH, NN_SSS,          NN_SST    ]
    vname = [ 'SSH'     ,  'SSS'      , 'SST'     ]
    vunit = [ r'm'    ,  r'PSU'   , r'$^{\circ}$C']

    
    jvar=0
    for cvar in vvar:
        cdiag = 'mean_'+cvar
        print '\n Treating '+cdiag
    
        for cocean in bo.voce2treat:
            Xf[:,:] = 0. ; jrun = 0
            for confrun in clist_confruns:

                cf_in = DIAG_DIR+'/'+confrun+'/'+cdiag+'_'+confrun+'_'+cocean+'.nc'
                vt0, vd0 = bn.read_1d_series(cf_in, cvar, cv_t='time', l_return_time=True)
                nbm = len(vt0) ; print ' *** nmb =', nbm
                if nbm%12 != 0: print 'ERROR: compare_time_series.py => PROBLEM#1. diag ='+cdiag; sys.exit(0)
                if nbm/12 > nb_years: print 'ERROR: compare_time_series.py => PROBLEM#2. diag ='+cdiag, nbm/12, nb_years ; sys.exit(0)
                VY, FY = bt.monthly_2_annual(vt0, vd0)
                vtime[:nbm/12]   = VY[:]
                Xf[jrun,:nbm/12] = FY[:] ; Xf[jrun,nbm/12:] = -999.
                jrun = jrun + 1
    
            bp.plot("1d_multi")(vtime[:], Xf[:,:], clist_confruns, cfig_type=FIG_FORMAT,
                                cfignm=cdiag+'_comparison_'+cocean, dt_year=ittic, loc_legend=DEFAULT_LEGEND_LOC,
                                cyunit=vunit[jvar], ctitle = vname[jvar]+', '+cocean, ymin=0, ymax=0)
    
        jvar = jvar+1
    
    









if imld == 1:

    cvar  = NN_MLD
    cdiag = 'mean_'+cvar
    vname = [ r'Mixed layer depth' ]
    lplot = True
    
    print '\n Treating '+cdiag
    Xf[:,:] = 0
    jbox = 0
    for cbox in bo.cname_mld_boxes:
        jrun = 0
        for confrun in clist_confruns:
            cf_in = DIAG_DIR+'/'+confrun+'/'+cdiag+'_'+confrun+'_'+cbox+'.nc'
            if os.path.exists(cf_in):
                vt0, vd0 = bn.read_1d_series(cf_in, cvar, cv_t='time', l_return_time=True)
                nbm = len(vt0)
                if nbm%12 != 0: print 'ERROR: compare_time_series.py => PROBLEM#1. diag ='+cdiag; sys.exit(0)
                if nbm/12 > nb_years: print 'ERROR: compare_time_series.py => PROBLEM#2. diag ='+cdiag; sys.exit(0)
                VY, FY = bt.monthly_2_annual(vt0, vd0)
                vtime[:nbm/12]   = VY[:]
                Xf[jrun,:nbm/12] = FY[:] ; Xf[jrun,nbm/12:] = -999.
                lplot = lplot and lplot
            jrun = jrun + 1
    
        if lplot:
            bp.plot("1d_multi")(vtime[:], Xf[:,:], clist_confruns, cfig_type=FIG_FORMAT,
                                cfignm=cdiag+'_'+cbox+'_comparison', dt_year=ittic, loc_legend=DEFAULT_LEGEND_LOC,
                                cyunit='m', ctitle = 'Mixed layer depth, '+bo.clgnm_mld_boxes[jbox], ymin=0, ymax=0)
        jbox = jbox+1
    
    
    
    
    
    
    
    
    
    
    
# Several columns to read
#=========================

if i3dfl == 1:

    vvar  = [ NN_S       ,      NN_T ]
    vname = [ 'Salinity' , 'Potential Temperature' ]
    vunit = [ r'PSU'     ,  r'$^{\circ}$C' ]
    
    vdepth_infil = [ '0-bottom', '0-100', '100-1000', '1000-bottom' ] ; # for 3d_thetao and 3d_so
    vdepth_range = [ 'All', '0m-100m', '100-1000m', '1000m-bottom' ] ; # for 3d_thetao and 3d_so
    
    jdiag=0
    for cvar in vvar:
        cdiag = '3d_'+cvar
        print '\n Treating '+cdiag
    
        for cocean in bo.voce2treat:
            
            # along the 4 columns of temperature
            idepth=0
            for cdepth in vdepth_range:
                cdif = vdepth_infil[idepth]
                Xf[:,:] = 0.
                jrun = 0
                for confrun in clist_confruns:
                    
                    cf_in = DIAG_DIR+'/'+confrun+'/'+cdiag+'_'+confrun+'_'+cocean+'.nc'
                    vt0, vd0 = bn.read_1d_series(cf_in, cvar+'_'+cdif, cv_t='time', l_return_time=True)
                    nbm = len(vt0)
                    if nbm%12 != 0: print 'ERROR: compare_time_series.py => PROBLEM#1. diag ='+cdiag; sys.exit(0)
                    if nbm/12 > nb_years: print 'ERROR: compare_time_series.py => PROBLEM#2. diag ='+cdiag; sys.exit(0)
                    VY, FY = bt.monthly_2_annual(vt0, vd0)
                    vtime[:nbm/12]   = VY[:]
                    Xf[jrun,:nbm/12] = FY[:]  ; Xf[jrun,nbm/12:] = -999.
                    jrun = jrun + 1
    
                bp.plot("1d_multi")(vtime[:], Xf[:,:], clist_confruns, cfig_type=FIG_FORMAT,
                                    cfignm=cdiag+'_comparison_'+cocean+'_'+cdepth, dt_year=ittic, loc_legend=DEFAULT_LEGEND_LOC,
                                    cyunit=vunit[jdiag], ctitle = vname[jdiag]+', '+cocean+', depth range = '+cdepth, ymin=0, ymax=0)
    
                idepth = idepth + 1
    
    
        jdiag = jdiag+1
    
    
    
    


# Sea-ice
#########

if iice == 1:

    vvar  = [ 'volu'           ,     'area'       ]
    vvnm  = [ 'Sea-ice Volume' , 'Sea-ice Extent' ]
    vunit = [ r'$10^3km^3$'    , r'$10^6$km$^2$'  ]    
    
    vpole = [ 'Arctic', 'Antarctic' ]
    vlab  = [ 'ne'    , 'se'        ]

    jvar = -1
    for cvar in vvar:
        jvar = jvar + 1
        
        vmnth = [            2          ,              8             ]
        vname = [ vvnm[jvar]+' in March', vvnm[jvar]+' in September' ]
    
        for jdiag in range(len(vmnth)):
            print '\n Treating '+vname[jdiag]
    
            # along the 2 columns of Arcit/Antarctic
            ipole = 0
            for cpole in vpole:
                Xf[:,:] = 0.
                jrun = 0
                for confrun in clist_confruns:
                    cf_in = DIAG_DIR+'/'+confrun+'/seaice_diags.nc'                
                    vt0, vd0 = bn.read_1d_series(cf_in, cvar+'_'+vlab[ipole], cv_t='time', l_return_time=True)
                    nbm = len(vt0)
                    if nbm%12 != 0: print 'ERROR: compare_time_series.py => PROBLEM#1. diag ='+cdiag; sys.exit(0)
                    if nbm/12 > nb_years: print 'ERROR: compare_time_series.py => PROBLEM#2. diag ='+cdiag, nbm/12, nb_years ; sys.exit(0)
                    vtime[:nbm/12]   = vt0[vmnth[jdiag]::12]
                    Xf[jrun,:nbm/12] = vd0[vmnth[jdiag]::12] ; Xf[jrun,nbm/12:] = -999.
                    jrun = jrun + 1

                cdiag = 'seaice_'+cvar
                cmnth = '%2.2i'%(vmnth[jdiag]+1)
                bp.plot("1d_multi")(vtime, Xf, clist_confruns, cfig_type=FIG_FORMAT,
                                    cfignm=cdiag+'_m'+str(cmnth)+'_comparison_'+cpole, dt_year=ittic, loc_legend='upper center',
                                    cyunit=vunit[jdiag], ctitle = vname[jdiag]+', '+cpole, ymin=0, ymax=0)
    
                ipole = ipole + 1
    
    
    
    
    
    
    
    
    

# Transport through sections
############################

if itrsp == 1:

    TRANSPORT_SECTION_FILE = os.getenv('TRANSPORT_SECTION_FILE')
    if TRANSPORT_SECTION_FILE == None: print 'The TRANSPORT_SECTION_FILE environement variable is no set'; sys.exit(0)

    print '\nUsing TRANSPORT_SECTION_FILE = '+TRANSPORT_SECTION_FILE
    list_sections = bo.get_sections_names_from_file(TRANSPORT_SECTION_FILE)
    print 'List of sections to treat: ', list_sections
    nbsect = len(list_sections)
        
    vstuff = [ 'volume', 'heat' , 'salt' ]
    vunit  = [ 'Sv'    , 'PW'   , 'kt/s' ]


    jrun = 0
    for confrun in clist_confruns:
    
        jsect=0
        for csect in list_sections:
            print '\n Treating transports through '+csect

            cf_in = DIAG_DIR+'/'+confrun+'/transport_sect_'+csect+'.nc' ; bt.chck4f(cf_in, script_name='compare_time_series.py')
            id_in = Dataset(cf_in)
            if jsect == 0:
                if jrun == 0:
                    vtime = nmp.zeros(nb_years*12)
                    vyear = nmp.zeros(nb_years)
                    Xtrsp   = nmp.zeros(nbrun*nb_years*3*nbsect) ; Xtrsp.shape = [ nbrun, nbsect , 3, nb_years ]
                
                vtime_t = id_in.variables['time'][:] ; nbm = len(vtime_t) ; nby = nbm/12
                if nbm%12 != 0: print 'ERROR: compare_time_series.py => PROBLEM#1. diag ='+cdiag; sys.exit(0)
                if nby > nb_years: print 'ERROR: compare_time_series.py => PROBLEM#2. diag ='+cdiag, nby, nb_years ; sys.exit(0)
                

                if nby == nb_years: vtime[:] = vtime_t[:]

            Xtrsp[jrun,jsect,:,:] = -999.
            vyear[:nby], Xtrsp[jrun,jsect,0,:nby] = bt.monthly_2_annual(vtime_t, id_in.variables['trsp_volu'][:nbm])
            vyear[:nby], Xtrsp[jrun,jsect,1,:nby] = bt.monthly_2_annual(vtime_t, id_in.variables['trsp_heat'][:nbm])
            vyear[:nby], Xtrsp[jrun,jsect,2,:nby] = bt.monthly_2_annual(vtime_t, id_in.variables['trsp_salt'][:nbm])
            
            id_in.close()
    
            jsect = jsect + 1    
        jrun = jrun + 1

    # All data read!
    
    jsect=0
    for csect in list_sections:
        jstuff = 0
        for cstuff in vstuff:

            bp.plot("1d_multi")(vyear[:], Xtrsp[:,jsect,jstuff,:], clist_confruns, cfig_type=FIG_FORMAT,
                                cfignm='transport_'+cstuff+'_'+csect+'_comparison', dt_year=ittic, loc_legend='upper left',
                                cyunit=vunit[jstuff], ctitle = 'Transport of '+cstuff+' through section '+csect,
                                ymin=0, ymax=0)
            
            jstuff = jstuff + 1
        jsect = jsect+1
    



    


# AMOC
if iamoc == 1:

    LMOCLAT = os.getenv('LMOCLAT')
    if LMOCLAT == None: print 'The LMOCLAT environement variable is no set'; sys.exit(0)
    
    list_lat = LMOCLAT.split() ; nblat = len(list_lat)
    print '\n AMOC: '+str(nblat)+' latitude bands!'

    nbm_prev = 0
    jrun = 0
    for confrun in clist_confruns:
    
        jl = 0
        for clr in list_lat:
            [ c1, c2 ] = clr.split('-') ; clat_info = '+'+c1+'N+'+c2+'N'
            cf_in = DIAG_DIR+'/'+confrun+'/max_moc_atl_'+clat_info+'.nc' ; bt.chck4f(cf_in, script_name='compare_time_series.py')
            
            id_in = Dataset(cf_in)
            if jl == 0:
                if jrun == 0:
                    vtime = nmp.zeros(nb_years*12)
                    vyear = nmp.zeros(nb_years)
                    Xamoc = nmp.zeros(nbrun*nb_years*nblat) ; Xamoc.shape = [ nbrun, nblat, nb_years ]
                vtime_t = id_in.variables['time'][:] ; nbm = len(vtime_t) ; nby = nbm/12
                if nbm%12 != 0: print 'ERROR: compare_time_series.py => PROBLEM#1. diag ='+cdiag; sys.exit(0)
                if nby > nb_years: print 'ERROR: compare_time_series.py => PROBLEM#2. diag ='+cdiag, nby, nb_years ; sys.exit(0)

                if nby == nb_years: vtime[:] = vtime_t[:]
    
            Xamoc[jrun,jl,:] = -999.
            vyear[:nby], Xamoc[jrun,jl,:nby] = bt.monthly_2_annual(vtime_t[:nbm], id_in.variables['moc_atl'][:nbm])
            id_in.close()
    
            jl = jl + 1    
        jrun = jrun + 1
    
    
    jl = 0
    for clr in list_lat:

        bp.plot("1d_multi")(vyear[:], Xamoc[:,jl,:], clist_confruns, cfig_type=FIG_FORMAT,
                            cfignm='AMOC_'+clr+'_comparison', loc_legend='lower left',
                            dt_year=ittic, cyunit='Sv', ctitle = 'AMOC ('+clr+')', ymin=0, ymax=0)
    
        jl = jl + 1










print  '\n'+sys.argv[0]+' done...\n'

