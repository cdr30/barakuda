#!/usr/bin/env python

# L. Brodeau 2014

#####################################
import sys
import os
import numpy as nmp
from netCDF4 import Dataset

import barakuda_orca as bo
import barakuda_tool as bt

#####################################


DIAG_D = os.getenv('DIAG_D')
if DIAG_D == None: print 'The DIAG_D environement variable is no set'; sys.exit(0)

ORCA = os.getenv('ORCA')
if ORCA == None: print 'The ORCA environement variable is no set'; sys.exit(0)

RUN = os.getenv('RUN')
if RUN == None: print 'The RUN environement variable is no set'; sys.exit(0)

NN_SST = os.getenv('NN_SST')
if NN_SST == None: print 'The NN_SST environement variable is no set'; sys.exit(0)
NN_SSS = os.getenv('NN_SSS')
if NN_SSS == None: print 'The NN_SSS environement variable is no set'; sys.exit(0)

FILE_DEF_BOXES = os.getenv('FILE_DEF_BOXES')
if FILE_DEF_BOXES == None: print 'The FILE_DEF_BOXES environement variable is no set'; sys.exit(0)





print '\n '+sys.argv[0]+':'
print ' ORCA = '+ORCA;
print ' RUN = '+RUN

CONFRUN = ORCA+'-'+RUN




cnexec = sys.argv[0]

na = len(sys.argv)
if na != 3 and na != 5 :
    print 'Usage : '+cnexec+' <ORCA1_RUN_grid_T.nc> <year> (<name_sst> <name_sss>) '
    sys.exit(0)

cf_in  = sys.argv[1]
cyear  = sys.argv[2] ; jyear = int(cyear); cyear = '%4.4i'%jyear


cv_sst = NN_SST
cv_sss = NN_SSS

if na == 5:
    cv_sst = sys.argv[3]
    cv_sss = sys.argv[4]


print 'Current year is '+cyear+' !\n'



cf_mesh_mask = './mesh_mask.nc'
cf_basin     = './new_maskglo.nc'







bt.chck4f(cf_mesh_mask)
id_mm = Dataset(cf_mesh_mask)

ni = id_mm.dimensions
rmsk = id_mm.variables['tmask'][0,0,:,:]
rlon = id_mm.variables['glamt'][0,:,:]
rlat = id_mm.variables['gphit'][0,:,:]
Xe1t = id_mm.variables['e1t'][0,:,:]
Xe2t = id_mm.variables['e2t'][0,:,:]
id_mm.close()

[ nj, ni ] = rmsk.shape



bt.chck4f(cf_in)
id_in = Dataset(cf_in)

if cv_sst == 'thetao':
    XSST  = id_in.variables[cv_sst][:,0,:,:]
else:
    XSST  = id_in.variables[cv_sst][:,:,:]
    
if cv_sss == 'so':
    XSSS  = id_in.variables[cv_sss][:,0,:,:]
else:
    XSSS  = id_in.variables[cv_sss][:,:,:]
    
id_in.close()



# Need Atlantic basin mask:
#bt.chck4f(cf_basin)
#id_mm = Dataset(cf_basin)
#rmsk_atl = id_mm.variables['tmaskatl'][0,:,:]
#id_mm.close()




[ nt, nj0, ni0 ] = XSST.shape

if [ nj0, ni0 ] != [ nj, ni ]: print 'ERROR: ssx_boxes.py => mask and field disagree in shape!'; sys.exit(0)

print 'nt, nj, ni =', nt, nj, ni



vtime = nmp.zeros(nt)
for jt in range(nt): vtime[jt] = float(jyear) + (float(jt) + 0.5)/float(nt)






# First will read name and coordinates of rectangular boxes to treat into file FILE_DEF_BOXES
##############################################################################################
vboxes, vi1, vj1, vi2, vj2 = bt.read_box_coordinates_in_ascii(FILE_DEF_BOXES)
nbb = len(vboxes)
print ''

rmean_sst = nmp.zeros(nt*nbb) ; rmean_sst.shape = [nt,nbb]
rmean_sss = nmp.zeros(nt*nbb) ; rmean_sss.shape = [nt,nbb]



for jb in range(nbb):

    cbox = vboxes[jb] ; print '\n   *** Focus on '+cbox+' box'

    i1 = vi1[jb]
    j1 = vj1[jb]
    i2 = vi2[jb]+1
    j2 = vj2[jb]+1
    
    nx_b = i2 - i1
    ny_b = j2 - j1

    Xbox = nmp.zeros(ny_b*nx_b) ; Xbox.shape = [ ny_b, nx_b ]
    

    idx_msk = nmp.where( rmsk[j1:j2,i1:i2] < 0.5 )



    for jt in range(nt):
        rmean_sst[jt,jb] = nmp.sum(XSST[jt,j1:j2,i1:i2]*Xe1t[j1:j2,i1:i2]*Xe2t[j1:j2,i1:i2]*rmsk[j1:j2,i1:i2]) \
                           / nmp.sum(Xe1t[j1:j2,i1:i2]*Xe2t[j1:j2,i1:i2]*rmsk[j1:j2,i1:i2])
        rmean_sss[jt,jb] = nmp.sum(XSSS[jt,j1:j2,i1:i2]*Xe1t[j1:j2,i1:i2]*Xe2t[j1:j2,i1:i2]*rmsk[j1:j2,i1:i2]) \
                           / nmp.sum(Xe1t[j1:j2,i1:i2]*Xe2t[j1:j2,i1:i2]*rmsk[j1:j2,i1:i2])
    


    # NETCDF:
    # Time to save both montly spatially-averaged time-series and 2D boxes
    # in netcdf files...

    cf_out =  DIAG_D+'/'+'SSX_'+cbox+'_'+CONFRUN+'.nc'
    l_nc_is_new = not os.path.exists(cf_out)


    # Creating/Opening output Netcdf file:
    if l_nc_is_new:
        f_out = Dataset(cf_out, 'w', format='NETCDF3_CLASSIC')
    else:
        f_out = Dataset(cf_out, 'a', format='NETCDF3_CLASSIC')

    if l_nc_is_new:
        jrec2write = 0
        
        # Creating Dimensions:
        f_out.createDimension('time', None)
        f_out.createDimension('y', ny_b)
        f_out.createDimension('x', nx_b)

        # Creating variables:
        id_t = f_out.createVariable('time','f4',('time',)) ;  id_t.units = 'year'
        id_lat = f_out.createVariable('nav_lat','f4',('y','x',)) ;  id_lat.units = 'deg.N'
        id_lon = f_out.createVariable('nav_lon','f4',('y','x',)) ;  id_lon.units = 'deg.E'
        id_v01   = f_out.createVariable(cv_sst+'_sa' ,'f4',('time',))
        id_v02   = f_out.createVariable(cv_sss+'_sa' ,'f4',('time',))
        id_x01   = f_out.createVariable(cv_sst ,'f4',('time','y','x',), fill_value=-9999.)
        id_x02   = f_out.createVariable(cv_sss ,'f4',('time','y','x',), fill_value=-9999.)

        # Writing depth vector
        id_lat[:,:] = rlat[j1:j2,i1:i2]
        id_lon[:,:] = rlon[j1:j2,i1:i2]

        for jt in range(nt):
            id_t[jt]   = vtime[jt]
            id_v01[jt] = rmean_sst[jt,jb]
            id_v02[jt] = rmean_sss[jt,jb]
            
            Xbox[:,:] = XSST[jt,j1:j2,i1:i2]
            Xbox[idx_msk] = -9999.
            id_x01[jt,:,:] = Xbox[:,:]

            Xbox[:,:] = XSSS[jt,j1:j2,i1:i2]
            Xbox[idx_msk] = -9999.
            id_x02[jt,:,:] = Xbox[:,:]

            
        f_out.Author = 'L. Brodeau (ssx_boxes.py of Barakuda)'

    else:
        vt = f_out.variables['time']
        jrec2write = len(vt)
        v01 = f_out.variables[cv_sst+'_sa']
        x01 = f_out.variables[cv_sst]
        v02 = f_out.variables[cv_sss+'_sa']
        x02 = f_out.variables[cv_sss]
        for jt in range(nt):
            vt [jrec2write+jt] = vtime[jt]
            v01[jrec2write+jt] = rmean_sst[jt,jb]
            v02[jrec2write+jt] = rmean_sss[jt,jb]

            Xbox[:,:] = XSST[jt,j1:j2,i1:i2] ; Xbox[idx_msk] = -9999.
            x01[jrec2write+jt,:,:] = Xbox[:,:]

            Xbox[:,:] = XSSS[jt,j1:j2,i1:i2] ; Xbox[idx_msk] = -9999.
            x02[jrec2write+jt,:,:] = Xbox[:,:]

    f_out.close()
    print cf_out+' written!'


