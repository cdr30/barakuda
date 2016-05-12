#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate misc. spatial averaging out of NEMO output files...
#
#       L. Brodeau, november 2013

import sys
import os
import numpy as nmp

from netCDF4 import Dataset
from string  import replace

import barakuda_tool as bt
import barakuda_orca as bo
import barakuda_ncio as bnc

# Box nino 3.4:
lon1_nino = 360. - 170.  ; # east
lat1_nino = -5.
lon2_nino = 360. - 120.  ; # east
lat2_nino = 5.

venv_needed = {'ORCA','RUN','DIAG_D','MM_FILE','BM_FILE','NEMO_SAVED_FILES','FILE_FLX_SUFFIX','NN_FWF','NN_EMP','NN_P','NN_RNF','NN_SST','NN_SSS','NN_SSH','NN_T','NN_S','NN_MLD'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

if len(sys.argv) != 3:
    print 'Usage : sys.argv[1] <ORCA1_RUN_grid_T.nc> <year>'
    sys.exit(0)

cnexec = sys.argv[0]
cf_T_in  = sys.argv[1]
cyear  = sys.argv[2] ; jyear = int(cyear); cyear = '%4.4i'%jyear

print 'Current year is '+cyear+' !\n'



#lolo:
vtime = nmp.zeros(12)
for jt in range(12):
    vtime[jt] = float(jyear) + (float(jt)+0.5)*1./12.




# Checking if the land-sea mask file is here:
for cf in [vdic['MM_FILE'], vdic['BM_FILE']]:
    if not os.path.exists(cf):
        print 'Mask file '+cf+' not found'; sys.exit(0)

# Reading the grid metrics:
id_mm = Dataset(vdic['MM_FILE'])
list_variables = id_mm.variables.keys()
rmask  = id_mm.variables['tmask'][0,:,:,:]
rlon   = id_mm.variables['glamt'][0,:,:]
rlat   = id_mm.variables['gphit'][0,:,:]
re1t   = id_mm.variables['e1t'][0,:,:]
re2t   = id_mm.variables['e2t'][0,:,:]
if 'e3t_0' in list_variables[:]:
    Xe3t = id_mm.variables['e3t_0'][0,:,:,:] # we need the 3D field becaus partial steps!!!
else:
    print 'ERROR: '+cnexec+' => how do we retrieve 3D e3t???'; sys.exit(0)
id_mm.close()

[ nk, nj, ni ] = rmask.shape


Xarea_t = nmp.zeros((nj, ni))
Xarea_t[:,:] = re1t[:,:]*re2t[:,:]*rmask[0,:,:]
Socean = nmp.sum( Xarea_t[:,:] )
print '\n  *** Surface of the ocean = ', Socean* 1.E-12, '  [10^6 km^2]\n'



cfe_sflx = vdic['FILE_FLX_SUFFIX']
l_fwf = False
if cfe_sflx in vdic['NEMO_SAVED_FILES']:
    l_fwf = True
    cf_F_in = replace(cf_T_in, 'grid_T', cfe_sflx)

Xe1t = nmp.zeros((nk, nj, ni))
Xe2t = nmp.zeros((nk, nj, ni))
for jk in range(nk):
    Xe1t[jk,:,:] = re1t[:,:]
    Xe2t[jk,:,:] = re2t[:,:]
    
del re1t, re2t

print 'Opening different basin masks in file '+vdic['BM_FILE']
id_bm = Dataset(vdic['BM_FILE'])
mask_atl = id_bm.variables['tmaskatl'][:,:]
mask_pac = id_bm.variables['tmaskpac'][:,:]
mask_ind = id_bm.variables['tmaskind'][:,:]
id_bm.close()

mask = nmp.zeros((4,nk,nj,ni))

mask[0,:,:,:] = rmask[:,:,:] ; # global
for jk in range(nk):
    mask[1,jk,:,:] = mask_atl[:,:]*rmask[jk,:,:]
    mask[2,jk,:,:] = mask_pac[:,:]*rmask[jk,:,:]
    mask[3,jk,:,:] = mask_ind[:,:]*rmask[jk,:,:]
del rmask, mask_atl, mask_pac, mask_ind







##############################################################
# Time-series of globally averaged surface freshwater fluxes #
##############################################################

if l_fwf:

    cv_fwf = vdic['NN_FWF']
    cv_emp = vdic['NN_EMP']
    cv_prc = vdic['NN_P']
    cv_rnf = vdic['NN_RNF']

    id_in = Dataset(cf_F_in)
    list_variables = id_in.variables.keys()
    FWF_m = id_in.variables[cv_fwf][:,:,:]
    print '   *** E-P-R ('+cv_fwf+') read!'

    l_emp = False
    if  cv_emp in list_variables[:]:
        l_emp = True
        EMP_m = id_in.variables[cv_emp][:,:,:]
        print '   *** E-P ('+cv_emp+') read!'
             
    l_prc = False
    if  cv_prc in list_variables[:]:
        l_prc = True
        PRC_m = id_in.variables[cv_prc][:,:,:]
        print '   *** P ('+cv_prc+') read!'

    l_rnf = False
    if  cv_rnf in list_variables[:]:
        l_rnf = True
        RNF_m = id_in.variables[cv_rnf][:,:,:]
        print '   *** P ('+cv_rnf+') read!'

    id_in.close()

               
    [ nt, nj0, ni0 ] = FWF_m.shape
    
    if l_emp and not l_rnf:
        l_rnf = True
        RNF_m = nmp.zeros((nj0,ni0))
        RNF_m = - ( FWF_m - EMP_m )

    vtime = nmp.zeros(nt)
    for jt in range(nt): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./12.

    vfwf = nmp.zeros(nt)
    
    vemp = [] ; vrnf = [] ; vprc = []
    if l_emp: vemp = nmp.zeros(nt)
    if l_rnf: vrnf = nmp.zeros(nt)
    if l_prc: vprc = nmp.zeros(nt)


    for jt in range(nt):
        vfwf[jt]           = nmp.sum( FWF_m[jt,:,:]*Xarea_t ) * 1.E-9 ;  # to Sv
        if l_emp: vemp[jt] = nmp.sum( EMP_m[jt,:,:]*Xarea_t ) * 1.E-9 ;  # to Sv
        if l_rnf: vrnf[jt] = nmp.sum( RNF_m[jt,:,:]*Xarea_t ) * 1.E-9 ;  # to Sv
        if l_prc: vprc[jt] = nmp.sum( PRC_m[jt,:,:]*Xarea_t ) * 1.E-9 ;  # to Sv

    cf_out   = vdic['DIAG_D']+'/mean_fwf_'+CONFRUN+'_global.nc'

    bnc.wrt_appnd_1d_series(vtime, vfwf, cf_out, 'EmPmR',
                            cu_t='year', cu_d='Sv', cln_d ='Globally averaged net freshwater flux (nemo:'+cv_fwf+')',
                            vd2=vemp, cvar2='EmP',  cln_d2='Globally averaged Evap - Precip (nemo:'+cv_emp+')',
                            vd3=vrnf, cvar3='R',    cln_d3='Globally averaged continental runoffs',
                            vd4=vprc, cvar4='P',    cln_d4='Globally averaged total precip (nemo:'+cv_prc+')',)








####################################
# MLD time serie in different boxes:
####################################
l_mld = False
print '\nSpatially-averaged MLD in different boxes'

cvar = vdic['NN_MLD']
id_in = Dataset(cf_T_in)
list_variables = id_in.variables.keys()
if cvar in list_variables[:]: # check if MLD variable is present!
    MLD_m = id_in.variables[cvar][:,:,:]
    print '   *** MLD ('+cvar+') found and read!'
    l_mld = True
else:
    print '   *** OOPS! MLD ('+cvar+') not found, skipping MLD time series diag...'
id_in.close()

if l_mld:

    [ nt, nj0, ni0 ] = MLD_m.shape

    if nt != 12: print 'ERROR: '+cnexec+' => only treating monthly data so far...'; sys.exit(0)

    if [ nj0, ni0 ] != [ nj, ni ]: print 'ERROR: '+cnexec+' => Field and metrics do not agree in size!'; sys.exit(0)

    vtime = nmp.zeros(nt)
    for jt in range(nt): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./12.

    mask2d = nmp.zeros((nj,ni))


    # Reading boxes definitions into barakuda_orca.py:
    cname_b = bo.cname_mld_boxes
    nb_boxes = len(cname_b)

    for ib in range(nb_boxes):

        cbox = cname_b[ib] ; print '    *** treating '+cvar+' for '+cbox+', ('+bo.clgnm_mld_boxes[ib]+')'

        i1 = 0 ; j1 = 0 ; i2 = ni-1 ; j2 = nj-1

        rx1 = bo.r_lon_p1_mld[ib] ; rx2 = bo.r_lon_p2_mld[ib] ; ry1 = bo.r_lat_p1_mld[ib] ; ry2 = bo.r_lat_p2_mld[ib]

        # Need to itterate because ORCA grid disytorded in the North...
        vold = [ -999, -999, -999, -999 ] ;  itt = 0
        while [ i1, i2, j1, j2 ] != vold and itt < 10 :
            itt = itt+1
            #print ' .... itt =', itt
            vold = [ i1, i2, j1, j2 ]
            #print 'seraching for rx1, rx2, ry1, ry2 = ', rx1, rx2, ry1, ry2
            if rx1 > -900.: i1 = bt.find_index_from_value( rx1, rlon[j1,:] )
            if rx2 > -900.: i2 = bt.find_index_from_value( rx2, rlon[j2,:] )
            if ry1 > -900.: j1 = bt.find_index_from_value( ry1, rlat[:,i1] )
            if ry2 > -900.: j2 = bt.find_index_from_value( ry2, rlat[:,i2] )
            #print '   => i1, i2, j1, j2 =>', i1, i2, j1, j2, '\n'



        mask2d[:,:] = 0.
        mask2d[j1:j2,i1:i2] = mask[0,0,j1:j2,i1:i2]

        Vts = bo.mean_2d(MLD_m, mask2d[:,:], Xe1t[0,:,:], Xe2t[0,:,:])

        # NETCDF:
        cf_out   = vdic['DIAG_D']+'/mean_'+cvar+'_'+CONFRUN+'_'+cbox+'.nc' ;  cv1 = cvar

        bnc.wrt_appnd_1d_series(vtime, Vts, cf_out, cv1,
                                cu_t='year', cu_d='m', cln_d='2D-average of '+cvar+' on rectangular box '+cbox)








#############################################
# 2D (surface) averaging for temperature and salinity #
#############################################

jvar = 0

for cvar in [ vdic['NN_SST'], vdic['NN_SSS'], vdic['NN_SSH'] ]:

    # DATA:
    print '  *** reading '+cvar+' into '+cf_T_in
    id_in = Dataset(cf_T_in)
    if cvar == 'thetao' or cvar == 'so':
        Xs_m = id_in.variables[cvar][:,0,:,:]
    else:
        Xs_m = id_in.variables[cvar][:,:,:]
    id_in.close()
    print '  ...read!'


    [ nt, nj0, ni0 ] = Xs_m.shape

    if nt != 12: print 'ERROR: '+cnexec+' => only treating monthly data so far...'; sys.exit(0)

    if [ nj0, ni0 ] != [ nj, ni ]: print 'ERROR: '+cnexec+' => Field and metrics do not agree in size!'; sys.exit(0)

    if jvar == 0:
        vtime = nmp.zeros(nt)
        for jt in range(nt): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./12.
        print ' * Montly calendar: ', vtime[:]


    joce = 0

    for cocean in bo.voce2treat[:]:

        print 'Treating '+cvar+' for '+cocean

        Vts = bo.mean_2d(Xs_m, mask[joce,0,:,:], Xe1t[0,:,:], Xe2t[0,:,:])


        #if 'cf_out' in locals() or 'cf_out' in globals():
        #    f = open(cf_out, 'a'); # 'w' would erase...
        #    f.write('#      Year       Mean ('+cocean+')\n')
        #    for jt in range(nt):
        #        f.write('%13.6f'%vtime[jt])
        #        f.write('  '+'%13.6f'%round(Vts[jt],6))
        #        f.write('\n')
        #    f.close()
        #    print cf_out+' written!'


        # NETCDF:
        cf_out   = vdic['DIAG_D']+'/mean_'+cvar+'_'+CONFRUN+'_'+cocean+'.nc' ;  cv1 = cvar
        bnc.wrt_appnd_1d_series(vtime, Vts, cf_out, cv1,
                                cu_t='year', cu_d='m', cln_d='2D-average of '+cvar+' on ocean '+cocean)

        joce = joce + 1

    jvar = jvar + 1



print '\n'






###################
# El nino box 3.4 #
###################

[i1, j1] = bo.find_ij(lon1_nino, lat1_nino, rlon, rlat, 'c')
[i2, j2] = bo.find_ij(lon2_nino, lat2_nino, rlon, rlat, 'c')
print ' Nino box 3.4, longitude: '+str(rlon[10,i1])+' => '+str(rlon[10,i2])+' \ latitude: '+str(rlat[j1,50])+' => '+str(rlat[j2,50])

id_in = Dataset(cf_T_in)
if vdic['NN_SST'] == 'thetao':
    Xs_m = id_in.variables[vdic['NN_SST']][:,0,:,:]
else:
    Xs_m = id_in.variables[vdic['NN_SST']][:,:,:]
id_in.close()

Vts = bo.mean_2d(Xs_m[:,j1:j2+1,i1:i2+1], mask[0,0,j1:j2+1,i1:i2+1], Xe1t[0,j1:j2+1,i1:i2+1], Xe2t[0,j1:j2+1,i1:i2+1])


#if 'cf_out' in locals() or 'cf_out' in globals():
#    f = open(cf_out, 'a'); # 'w' would erase...
#    f.write('#      Year       Mean ()\n')
#    for jt in range(nt):
#        f.write('%13.6f'%vtime[jt])
#        f.write('  '+'%13.6f'%round(Vts[jt],6))
#        f.write('\n')
#    f.close()
#    print cf_out+' written!'

# NETCDF:
cf_out   = vdic['DIAG_D']+'/Nino34_'+CONFRUN+'.nc' ;  cv1 = vdic['NN_SST']
bnc.wrt_appnd_1d_series(vtime, Vts, cf_out, cv1,
                        cu_t='year', cu_d='K', cln_d='2D-average of SST Nino box 3.4')







#############################################
# 3D averaging for temperature and salinity #
#############################################

jvar = 0

for cvar in [ vdic['NN_T'] , vdic['NN_S'] ]:



    # DATA:
    id_in = Dataset(cf_T_in)
    vdepth = id_in.variables['deptht'][:]
    Xd_m = id_in.variables[cvar][:,:,:,:]
    id_in.close()


    if jvar == 0:
        j100m  = bt.find_index_from_value(100.  , vdepth) ; print 'j100m  = ', j100m,  '=> ', vdepth[j100m]
        j1000m = bt.find_index_from_value(1000. , vdepth) ; print 'j1000m = ', j1000m, '=> ', vdepth[j1000m]

    [ nt, nk0, nj0, ni0 ] = Xd_m.shape

    if nt != 12: print 'ERROR: '+cnexec+' => only treating monthly data so far...'; sys.exit(0)

    if [ nk0, nj0, ni0 ] != [ nk, nj, ni ]: print 'ERROR: '+cnexec+' => Field and metrics do not agree in size!'; sys.exit(0)

    if jvar == 0:
        vtime = nmp.zeros(nt)
        for jt in range(nt): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./12.
        print ' * Montly calendar: ', vtime[:]



    # Annual mean array for current year:
    Xd_y = nmp.zeros((1, nk, nj, ni))

    Xd_y[0,:,:,:] = nmp.mean(Xd_m, axis=0)



    joce = 0

    for cocean in bo.voce2treat[:]:


        print 'Treating '+cvar+' for '+cocean



        # I) Montly mean for diffrent depth ranges
        # ========================================

        Vts_tot   = bo.mean_3d(Xd_m, mask[joce,:,:,:], Xe1t, Xe2t, Xe3t) ; # Top to bottom
        Vts_0_100 = bo.mean_3d(Xd_m[:,:j100m,:,:], mask[joce,:j100m,:,:], Xe1t[:j100m,:,:], Xe2t[:j100m,:,:], Xe3t[:j100m,:,:])
        Vts_100_1000 = bo.mean_3d(Xd_m[:,j100m:j1000m,:,:], mask[joce,j100m:j1000m,:,:], Xe1t[j100m:j1000m,:,:], Xe2t[j100m:j1000m,:,:], Xe3t[j100m:j1000m,:,:])
        Vts_1000_bot = bo.mean_3d(Xd_m[:,j1000m:,:,:], mask[joce,j1000m:,:,:], Xe1t[j1000m:,:,:], Xe2t[j1000m:,:,:], Xe3t[j1000m:,:,:])

        cf_out = vdic['DIAG_D']+'/3d_'+cvar+'_'+CONFRUN+'_'+cocean+'.nc'
        cv1 = cvar+'_0-bottom'
        cv2 = cvar+'_0-100'
        cv3 = cvar+'_100-1000'
        cv4 = cvar+'_1000-bottom'


        bnc.wrt_appnd_1d_series(vtime, Vts_tot, cf_out, cv1,
                                cu_t='year', cu_d='Unknown', cln_d ='3D-average of '+cvar+': surface to bottom, '+cocean,
                                vd2=Vts_0_100,    cvar2=cv2, cln_d2='3D-average of '+cvar+': surface to 100m, '+cocean,
                                vd3=Vts_100_1000, cvar3=cv3, cln_d3='3D-average of '+cvar+': 100m to 1000m, '+cocean,
                                vd4=Vts_1000_bot, cvar4=cv4, cln_d4='3D-average of '+cvar+': 1000m to bottom, '+cocean)



        # II) Annual mean vertical profile
        # ================================

        Vf = nmp.zeros(nk)

        for jk in range(nk):

            [ rf ] = bo.mean_2d(Xd_y[:,jk,:,:], mask[joce,jk,:,:], Xe1t[jk,:,:], Xe2t[jk,:,:])

            Vf[jk] = rf



        # NETCDF:
        cf_out = vdic['DIAG_D']+'/'+cvar+'_mean_Vprofile_'+CONFRUN+'_'+cocean+'.nc'
        l_nc_is_new = not os.path.exists(cf_out)
        #
        # Creating/Opening output Netcdf file:
        if l_nc_is_new:
            f_out = Dataset(cf_out, 'w', format='NETCDF3_CLASSIC')
        else:
            f_out = Dataset(cf_out, 'a', format='NETCDF3_CLASSIC')

        if l_nc_is_new:
            jrec2write = 0

            # Creating Dimensions:
            f_out.createDimension('time', None)
            f_out.createDimension('deptht', nk)

            # Creating variables:
            id_t = f_out.createVariable('time','f4',('time',)) ;      id_t.units = 'year'
            id_z = f_out.createVariable('deptht','f4',('deptht',)) ;  id_z.units = 'm'
            id_v01   = f_out.createVariable(cvar ,'f4',('time','deptht',))
            id_v01.long_name = 'Horizontally-averaged '+cvar+': '+cocean
            # Writing depth vector
            id_z[:] = vdepth[:]
            id_t[jrec2write] = float(jyear)+0.5
            id_v01[jrec2write,:] = Vf[:]
            f_out.Author = 'L. Brodeau ('+cnexec+' of Barakuda)'

        else:
            vt = f_out.variables['time']
            jrec2write = len(vt)
            v01 = f_out.variables[cvar]
            vt[jrec2write] = float(jyear)+0.5
            v01[jrec2write,:] = Vf[:]

        f_out.close()
        print cf_out+' written!'





        print ''

        joce = joce + 1


    jvar = jvar + 1
    print '\n'









print '\n'



