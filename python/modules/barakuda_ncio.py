import numpy as nmp
import sys

from netCDF4 import Dataset

from os import path



def wrt_1d_series(vt, vd, cvar, cinfo,
                  cu_t='unknown', cu_d='unknown', cln_d='unknown', nsmooth=0,
                  vd2=[], vd3=[], vd4=[], vd5=[],
                  cvar2='', cvar3='', cvar4='', cvar5='',
                  cln_d2='', cln_d3='', cln_d4='', cln_d5=''):

    cf_o  = cvar+'_'+cinfo+'.nc'

    lsmooth = False

    if nsmooth > 0:
        import barakuda_stat as bs
        lsmooth = True

        if nsmooth == 11:
            vd_sm = bs.running_mean_11(vd, l_fill_bounds=False)
        elif nsmooth == 5:
            vd_sm = bs.running_mean_5(vd, l_fill_bounds=False)
        else:
            print 'ERROR: wrt_1d_series.barakuda_ncio => smoothing with nsmooth='+str(nsmooth)+' not supported!'; sys.exit(0)


    f_o = Dataset(cf_o, 'w', format='NETCDF3_CLASSIC')

    nt = len(vt)
    if len(vd) != nt:  print 'ERROR: wrt_1d_series.barakuda_ncio => data & time have different lengths!'; sys.exit(0)

    l_do_v2=False ; l_do_v3=False ; l_do_v4=False ; l_do_v5=False
    if len(vd2) == nt: l_do_v2=True
    if len(vd3) == nt: l_do_v3=True
    if len(vd4) == nt: l_do_v4=True
    if len(vd5) == nt: l_do_v5=True


    f_o.createDimension('time', None)
    id_t = f_o.createVariable('time','f4',('time',)) ;  id_t.units = cu_t

    id_d = f_o.createVariable(cvar,'f4',('time',))
    id_d.units = cu_d ;  id_d.long_name = cln_d

    if l_do_v2: id_d2 = f_o.createVariable(cvar2,'f4',('time',)); id_d2.units = cu_d; id_d2.long_name = cln_d2
    if l_do_v3: id_d3 = f_o.createVariable(cvar3,'f4',('time',)); id_d3.units = cu_d; id_d3.long_name = cln_d3
    if l_do_v4: id_d4 = f_o.createVariable(cvar4,'f4',('time',)); id_d4.units = cu_d; id_d4.long_name = cln_d4
    if l_do_v5: id_d5 = f_o.createVariable(cvar5,'f4',('time',)); id_d5.units = cu_d; id_d5.long_name = cln_d5


    if lsmooth:
        id_sm = f_o.createVariable(cvar+'_'+str(nsmooth)+'yrm','f4',('time',))
        id_sm.units = cu_d ;  id_sm.long_name = str(nsmooth)+'-year running mean of '+cln_d


    for jt in range(nt):
        id_t[jt]   = vt[jt]
        id_d[jt]   = vd[jt]
        if lsmooth: id_sm[jt] = vd_sm[jt]
        if l_do_v2: id_d2[jt] = vd2[jt]
        if l_do_v3: id_d3[jt] = vd3[jt]
        if l_do_v4: id_d4[jt] = vd4[jt]
        if l_do_v5: id_d5[jt] = vd5[jt]

    f_o.Author = 'L. Brodeau (barakuda_ncio.py of Barakuda)'
    f_o.close()
    print ' * wrt_1d_series => '+cf_o+' written!\n'

    return 0



def read_1d_series(cf_i, cv_i, cv_t='time', l_return_time=True):

    if not path.exists(cf_i): print 'ERROR: read_1d_series.barakuda_ncio => '+cf_i+' not found!'; sys.exit(0)

    id_i = Dataset(cf_i)
    if l_return_time: vt = id_i.variables[cv_t][:]
    vd = id_i.variables[cv_i][:]
    id_i.close()

    if l_return_time:
        print '  * read_1d_series => just read '+cv_t+' and '+cv_i+' into '+cf_i+'\n'
        return vt, vd
    else:
        print '  * read_1d_series => just read '+cv_i+' into '+cf_i+'\n'
        return vd






def read_1d_series_2(cf_i, cv_i, cv_t='time'):

    if not path.exists(cf_i): print 'ERROR: read_1d_series_2.barakuda_ncio => '+cf_i+' not found!'; sys.exit(0)

    id_i = Dataset(cf_i)
    vt   = id_i.variables[cv_t][:]
    nt   = len(vt)
    xout = nmp.zeros((2,nt))
    xout[0,:] = vt[:]

    xout[1,:] = id_i.variables[cv_i][:]

    id_i.close()

    print '  * read_1d_series_2 => just read '+cv_t+' and '+cv_i+' into '+cf_i+'\n'
    return xout





def wrt_appnd_1d_series(vt, vd, cf, cvar1,  cu_t='year', cu_d='unknown', cln_d='unknown',
                        vd2=[], cvar2='', cln_d2='', cun2='unknown',
                        vd3=[], cvar3='', cln_d3='', cun3='unknown',
                        vd4=[], cvar4='', cln_d4='', cun4='unknown',
                        vd5=[], cvar5='', cln_d5='', cun5='unknown',
                        vd6=[], cvar6='', cln_d6='', cun6='unknown',
                        vd7=[], cvar7='', cln_d7='', cun7='unknown',
                        vd8=[], cvar8='', cln_d8='', cun8='unknown',
                        vd9=[], cvar9='', cln_d9='', cun9='unknown',
                        vd10=[],cvar10='',cln_d10='',cun10='unknown' ):

    l_nc_is_new = not path.exists(cf)

    Nt = len(vt)
    if len(vd) != Nt:
        print 'ERROR: wrt_appnd_1d_series.barakuda_ncio => data & time have different lengths!'
        sys.exit(0)

    l_do_v2=False ; l_do_v3=False ; l_do_v4=False ; l_do_v5=False
    l_do_v6=False ; l_do_v7=False ; l_do_v8=False ; l_do_v9=False ; l_do_v10=False
    if len(vd2) == Nt: l_do_v2=True
    if len(vd3) == Nt: l_do_v3=True
    if len(vd4) == Nt: l_do_v4=True
    if len(vd5) == Nt: l_do_v5=True
    if len(vd6) == Nt: l_do_v6=True
    if len(vd7) == Nt: l_do_v7=True
    if len(vd8) == Nt: l_do_v8=True
    if len(vd9) == Nt: l_do_v9=True
    if len(vd10)== Nt: l_do_v10=True


    if l_nc_is_new:
        f_out = Dataset(cf, 'w', format='NETCDF3_CLASSIC')
    else:
        f_out = Dataset(cf, 'a', format='NETCDF3_CLASSIC')

    if l_nc_is_new:
        jrec2write = 0
        f_out.createDimension('time', None)
        id_t = f_out.createVariable('time','f4',('time',)) ;      id_t.units = 'year'

        id_v01   = f_out.createVariable(cvar1 ,'f4',('time',))
        id_v01.units     = cu_d
        id_v01.long_name = cln_d
        if l_do_v2:
            id_v02   = f_out.createVariable(cvar2 ,'f4',('time',))
            id_v02.units = cu_d
            if not cun2=='unknown': id_v02.units = cun2
            id_v02.long_name = cln_d2
        if l_do_v3:
            id_v03   = f_out.createVariable(cvar3 ,'f4',('time',))
            id_v03.units     = cu_d
            if not cun3=='unknown': id_v03.units = cun3
            id_v03.long_name = cln_d3
        if l_do_v4:
            id_v04   = f_out.createVariable(cvar4 ,'f4',('time',))
            id_v04.units     = cu_d
            if not cun4=='unknown': id_v04.units = cun4
            id_v04.long_name = cln_d4
        if l_do_v5:
            id_v05   = f_out.createVariable(cvar5 ,'f4',('time',))
            id_v05.units = cu_d
            if not cun5=='unknown': id_v05.units = cun5
            id_v05.long_name = cln_d5
        if l_do_v6:
            id_v06   = f_out.createVariable(cvar6 ,'f4',('time',))
            id_v06.units     = cu_d
            if not cun6=='unknown': id_v06.units = cun6
            id_v06.long_name = cln_d6
        if l_do_v7:
            id_v07   = f_out.createVariable(cvar7 ,'f4',('time',))
            id_v07.units     = cu_d
            if not cun7=='unknown': id_v07.units = cun7
            id_v07.long_name = cln_d7
        if l_do_v8:
            id_v08   = f_out.createVariable(cvar8 ,'f4',('time',))
            id_v08.units = cu_d
            if not cun8=='unknown': id_v08.units = cun8
            id_v08.long_name = cln_d8
        if l_do_v9:
            id_v09   = f_out.createVariable(cvar9 ,'f4',('time',))
            id_v09.units     = cu_d
            if not cun9=='unknown': id_v09.units = cun9
            id_v09.long_name = cln_d9
        if l_do_v10:
            id_v10   = f_out.createVariable(cvar10 ,'f4',('time',))
            id_v10.units     = cu_d
            if not cun7=='unknown': id_v10.units = cun10
            id_v10.long_name = cln_d10




        jrw = 0
        for jt in range(Nt):
            jrw = jrec2write + jt
            id_t[jrw]   = vt[jt]
            id_v01[jrw] = vd[jt]
            if l_do_v2: id_v02[jrw] = vd2[jt]
            if l_do_v3: id_v03[jrw] = vd3[jt]
            if l_do_v4: id_v04[jrw] = vd4[jt]
            if l_do_v5: id_v05[jrw] = vd5[jt]
            if l_do_v6: id_v06[jrw] = vd6[jt]
            if l_do_v7: id_v07[jrw] = vd7[jt]
            if l_do_v8: id_v08[jrw] = vd8[jt]
            if l_do_v9: id_v09[jrw] = vd9[jt]
            if l_do_v10:id_v10[jrw] = vd10[jt]
        f_out.about = 'Diagnostics created with BaraKuda (https://github.com/brodeau/barakuda)'
    else:
        vtr = f_out.variables['time']
        jrec2write = len(vtr)
        v01 = f_out.variables[cvar1]
        if l_do_v2: v02 = f_out.variables[cvar2]
        if l_do_v3: v03 = f_out.variables[cvar3]
        if l_do_v4: v04 = f_out.variables[cvar4]
        if l_do_v5: v05 = f_out.variables[cvar5]
        if l_do_v6: v06 = f_out.variables[cvar6]
        if l_do_v7: v07 = f_out.variables[cvar7]
        if l_do_v8: v08 = f_out.variables[cvar8]
        if l_do_v9: v09 = f_out.variables[cvar9]
        if l_do_v10:v10 = f_out.variables[cvar10]
        jrw = 0
        for jt in range(Nt):
            jrw = jrec2write + jt
            vtr[jrw]  = vt[jt]
            v01[jrw] = vd[jt]
            if l_do_v2: v02[jrw] = vd2[jt]
            if l_do_v3: v03[jrw] = vd3[jt]
            if l_do_v4: v04[jrw] = vd4[jt]
            if l_do_v5: v05[jrw] = vd5[jt]
            if l_do_v6: v06[jrw] = vd6[jt]
            if l_do_v7: v07[jrw] = vd7[jt]
            if l_do_v8: v08[jrw] = vd8[jt]
            if l_do_v9: v09[jrw] = vd9[jt]
            if l_do_v10:v10[jrw] = vd10[jt]
    f_out.close()

    print '   *** '+str(Nt)+' snapshots of fields saved into '+cf+' !\n'


