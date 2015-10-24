import numpy as nmp
import sys

from netCDF4 import Dataset

from os import path



def wrt_1d_series(vt, vd, cvar, cinfo,                   
                  cu_t='unknown', cu_d='unknown', cln_d='unknown', nsmooth=0,
                  vd2=[], vd3=[], vd4=[], vd5=[],
                  cvar2='', cvar3='', cvar4='', cvar5='',
                  cln_d2='', cln_d3='', cln_d4='', cln_d5='',):

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
