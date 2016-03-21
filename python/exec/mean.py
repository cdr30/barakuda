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

import barakuda_tool as bt
import barakuda_orca as bo

# Box nino 3.4:
lon1_nino = 360. - 170.  ; # east
lat1_nino = -5.
lon2_nino = 360. - 120.  ; # east
lat2_nino = 5.

venv_needed = {'ORCA','RUN','DIAG_D','MM_FILE','BM_FILE','NN_SST','NN_SSS','NN_SSH','NN_T','NN_S','NN_MLD'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

if len(sys.argv) != 3:
    print 'Usage : sys.argv[1] <ORCA1_RUN_grid_T.nc> <year>'
    sys.exit(0)

cnexec = sys.argv[0]
cf_in  = sys.argv[1]
cyear  = sys.argv[2] ; jyear = int(cyear); cyear = '%4.4i'%jyear

print 'Current year is '+cyear+' !\n'

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





####################################
# MLD time serie in different boxes:
####################################
l_mld = False
print '\nSpatially-averaged MLD in different boxes'

cvar = vdic['NN_MLD']

id_in = Dataset(cf_in)
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
    print ' * Montly calendar: ', vtime[:], '\n'

    mask2d = nmp.zeros((nj,ni))
    

    # Reading boxes definitions into barakuda_orca.py:
    cname_b = bo.cname_mld_boxes
    nb_boxes = len(cname_b)

    for ib in range(nb_boxes):

        cbox = cname_b[ib] ; print '    *** treating '+cvar+' for '+cbox+', ('+bo.clgnm_mld_boxes[ib]+')'

        i1 = 0 ; j1 = 0 ; i2 = ni-1 ; j2 = nj-1

        rx1 = bo.r_lon_p1_mld[ib] ; rx2 = bo.r_lon_p2_mld[ib] ; ry1 = bo.r_lat_p1_mld[ib] ; ry2 = bo.r_lat_p2_mld[ib]

        # Need to itterate because ORCA grid disytorded in the North...
        vold = [ 0, 0, 0, 0 ] ;  itt = 0
        while [ i1, i2, j1, j2 ] != vold and itt < 10 :
            itt = itt+1
            vold = [ i1, i2, j1, j2 ]
            if rx1 > -900.: i1 = bt.find_index_from_value( rx1, rlon[j1,:] )
            if rx2 > -900.: i2 = bt.find_index_from_value( rx2, rlon[j2,:] )
            if ry1 > -900.: j1 = bt.find_index_from_value( ry1, rlat[:,i1] )
            if ry2 > -900.: j2 = bt.find_index_from_value( ry2, rlat[:,i2] )
            
    
        mask2d[:,:] = 0.
        mask2d[j1:j2,i1:i2] = mask[0,0,j1:j2,i1:i2]

        Vts = bo.mean_2d(MLD_m, mask2d[:,:], Xe1t[0,:,:], Xe2t[0,:,:])
        
        # NETCDF:
        cf_out   = vdic['DIAG_D']+'/mean_'+cvar+'_'+CONFRUN+'_'+cbox+'.nc' ;  cv1 = cvar        
        l_nc_is_new = not os.path.exists(cf_out)
        if l_nc_is_new:
            f_out = Dataset(cf_out, 'w', format='NETCDF3_CLASSIC')
        else:
            f_out = Dataset(cf_out, 'a', format='NETCDF3_CLASSIC')
        if l_nc_is_new:
            jrec2write = 0
            f_out.createDimension('time', None)
            id_t = f_out.createVariable('time','f4',('time',)) ;      id_t.units = 'year'
            id_v01   = f_out.createVariable(cv1 ,'f4',('time',))
            id_v01.long_name = '2D-average of '+cvar+' on rectangular box '+cbox
            jrw = 0
            for jt in range(nt):
                jrw = jrec2write + jt
                id_t[jrw]   = float(jyear) + 1./12.*(float(jt)+0.5)
                id_v01[jrw] = Vts[jt]
            f_out.box_def = cbox+' => ji:'+str(i1)+'->'+str(i2)+' jj:'+str(j1)+'->'+str(j2)
            f_out.Author = 'L. Brodeau ('+cnexec+' of Barakuda)'
        else:
            vt = f_out.variables['time']
            jrec2write = len(vt)
            v01 = f_out.variables[cv1]
            jrw = 0
            for jt in range(nt):
                jrw = jrec2write + jt
                vt[jrw]  = float(jyear) + 1./12.*(float(jt)+0.5)
                v01[jrw] = Vts[jt]
        f_out.close()
        print cf_out+' written!'








#############################################
# 2D (surface) averaging for temperature and salinity #
#############################################

jvar = 0

for cvar in [ vdic['NN_SST'], vdic['NN_SSS'], vdic['NN_SSH'] ]:

    # DATA:
    print '  *** reading '+cvar+' into '+cf_in
    id_in = Dataset(cf_in)
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


        if 'cf_out' in locals() or 'cf_out' in globals():  
            f = open(cf_out, 'a'); # 'w' would erase...
            f.write('#      Year       Mean ('+cocean+')\n')
            for jt in range(nt):
                f.write('%13.6f'%vtime[jt])
                f.write('  '+'%13.6f'%round(Vts[jt],6))
                f.write('\n')
            f.close()
            print cf_out+' written!'


        # NETCDF:
        cf_out   = vdic['DIAG_D']+'/mean_'+cvar+'_'+CONFRUN+'_'+cocean+'.nc' ;  cv1 = cvar
        l_nc_is_new = not os.path.exists(cf_out)
        if l_nc_is_new:
            f_out = Dataset(cf_out, 'w', format='NETCDF3_CLASSIC')
        else:
            f_out = Dataset(cf_out, 'a', format='NETCDF3_CLASSIC')
        if l_nc_is_new:
            jrec2write = 0
            f_out.createDimension('time', None)
            id_t = f_out.createVariable('time','f4',('time',)) ;  id_t.units = 'year'
            id_v01   = f_out.createVariable(cv1 ,'f4',('time',))
            id_v01.long_name = '2D-average of '+cvar+' on ocean '+cocean
            jrw = 0
            for jt in range(nt):
                jrw = jrec2write + jt
                id_t[jrw]   = float(jyear) + 1./12.*(float(jt)+0.5) ; id_v01[jrw] = Vts[jt]
            f_out.Author = 'L. Brodeau ('+cnexec+' of Barakuda)'
        else:
            vt = f_out.variables['time']
            jrec2write = len(vt)
            v01 = f_out.variables[cv1]
            jrw = 0
            for jt in range(nt):
                jrw = jrec2write + jt
                vt[jrw]  = float(jyear) + 1./12.*(float(jt)+0.5) ; v01[jrw] = Vts[jt]
        f_out.close()
        print cf_out+' written!'

        joce = joce + 1

    jvar = jvar + 1



print '\n'






###################
# El nino box 3.4 #
###################

[i1, j1] = bo.find_ij(lon1_nino, lat1_nino, rlon, rlat, 'c')
[i2, j2] = bo.find_ij(lon2_nino, lat2_nino, rlon, rlat, 'c')
print ' Nino box 3.4, longitude: '+str(rlon[10,i1])+' => '+str(rlon[10,i2])+' \ latitude: '+str(rlat[j1,50])+' => '+str(rlat[j2,50])

id_in = Dataset(cf_in)
if vdic['NN_SST'] == 'thetao':
    Xs_m = id_in.variables[vdic['NN_SST']][:,0,:,:]
else:
    Xs_m = id_in.variables[vdic['NN_SST']][:,:,:]
id_in.close()

Vts = bo.mean_2d(Xs_m[:,j1:j2+1,i1:i2+1], mask[0,0,j1:j2+1,i1:i2+1], Xe1t[0,j1:j2+1,i1:i2+1], Xe2t[0,j1:j2+1,i1:i2+1])


if 'cf_out' in locals() or 'cf_out' in globals():  
    f = open(cf_out, 'a'); # 'w' would erase...
    f.write('#      Year       Mean ()\n')
    for jt in range(nt):
        f.write('%13.6f'%vtime[jt])
        f.write('  '+'%13.6f'%round(Vts[jt],6))
        f.write('\n')
    f.close()
    print cf_out+' written!'

# NETCDF:
cf_out   = vdic['DIAG_D']+'/Nino34_'+CONFRUN+'.nc' ;  cv1 = vdic['NN_SST']
l_nc_is_new = not os.path.exists(cf_out)
if l_nc_is_new:
    f_out = Dataset(cf_out, 'w', format='NETCDF3_CLASSIC')
else:
    f_out = Dataset(cf_out, 'a', format='NETCDF3_CLASSIC')
if l_nc_is_new:
    jrec2write = 0
    f_out.createDimension('time', None)
    id_t = f_out.createVariable('time','f4',('time',)) ;  id_t.units = 'year'
    id_v01   = f_out.createVariable(cv1 ,'f4',('time',)) ; id_v01.long_name = '2D-average of SST Nino box 3.4'
    jrw = 0
    for jt in range(nt):
        jrw = jrec2write + jt
        id_t[jrw]   = float(jyear) + 1./12.*(float(jt)+0.5) ; id_v01[jrw] = Vts[jt]
    f_out.Author = 'L. Brodeau ('+cnexec+' of Barakuda)'
else:
    vt = f_out.variables['time']
    jrec2write = len(vt)
    v01 = f_out.variables[cv1]
    jrw = 0
    for jt in range(nt):
        jrw = jrec2write + jt
        vt[jrw]  = float(jyear) + 1./12.*(float(jt)+0.5) ; v01[jrw] = Vts[jt]
f_out.close()
print cf_out+' written!'












#############################################
# 3D averaging for temperature and salinity #
#############################################

jvar = 0

for cvar in [ vdic['NN_T'] , vdic['NN_S'] ]:



    # DATA:
    id_in = Dataset(cf_in)
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
        
        l_nc_is_new = not os.path.exists(cf_out)
        if l_nc_is_new:
            f_out = Dataset(cf_out, 'w', format='NETCDF3_CLASSIC')
        else:
            f_out = Dataset(cf_out, 'a', format='NETCDF3_CLASSIC')

        if l_nc_is_new:
            jrec2write = 0
            f_out.createDimension('time', None)
            id_t = f_out.createVariable('time','f4',('time',)) ;      id_t.units = 'year'
            
            id_v01   = f_out.createVariable(cv1 ,'f4',('time',))
            id_v01.long_name = '3D-average of '+cvar+': surface to bottom, '+cocean

            id_v02   = f_out.createVariable(cv2 ,'f4',('time',))
            id_v02.long_name = '3D-average of '+cvar+': surface to 100m, '+cocean

            id_v03   = f_out.createVariable(cv3 ,'f4',('time',))
            id_v03.long_name = '3D-average of '+cvar+': 100m to 1000m, '+cocean

            id_v04   = f_out.createVariable(cv4 ,'f4',('time',))
            id_v04.long_name = '3D-average of '+cvar+': 1000m to bottom, '+cocean

            jrw = 0
            for jt in range(nt):
                jrw = jrec2write + jt
                id_t[jrw]   = float(jyear) + 1./12.*(float(jt)+0.5)
                id_v01[jrw] = Vts_tot[jt]
                id_v02[jrw] = Vts_0_100[jt]
                id_v03[jrw] = Vts_100_1000[jt]
                id_v04[jrw] = Vts_1000_bot[jt]
            f_out.Author = 'L. Brodeau ('+cnexec+' of Barakuda)'

        else:
            vt = f_out.variables['time']
            jrec2write = len(vt)
            v01 = f_out.variables[cv1]
            v02 = f_out.variables[cv2]
            v03 = f_out.variables[cv3]
            v04 = f_out.variables[cv4]
            
            jrw = 0
            for jt in range(nt):
                jrw = jrec2write + jt
                vt[jrw]  = float(jyear) + 1./12.*(float(jt)+0.5)
                v01[jrw] = Vts_tot[jt]
                v02[jrw] = Vts_0_100[jt]
                v03[jrw] = Vts_100_1000[jt]
                v04[jrw] = Vts_1000_bot[jt]

        f_out.close()
        print cf_out+' written!'




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



