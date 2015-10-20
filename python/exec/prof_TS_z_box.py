#!/usr/bin/python

# L. Brodeau, June 2014

import sys
import os
import numpy as nmp
from netCDF4 import Dataset
from os.path import basename

# Laurent's:
import barakuda_physics as bphy
import barakuda_plot as brkdp
import barakuda_tool as bt


#l_plot_debug = True
l_plot_debug = False


FILE_DEF_BOXES = os.getenv('FILE_DEF_BOXES')
if FILE_DEF_BOXES == None: print 'The FILE_DEF_BOXES environement variable is no set'; sys.exit(0)

CONF = os.getenv('CONF')
if CONF == None: print 'The CONF environement variable is no set'; sys.exit(0)

RUN = os.getenv('RUN')
if RUN == None: print 'The RUN environement variable is no set'; sys.exit(0)
CONFRUN = CONF+'-'+RUN

CPREF = os.getenv('CPREF')
if CPREF == None: print 'The CPREF environement variable is no set'; sys.exit(0)

MM_FILE = os.getenv('MM_FILE')
if MM_FILE == None: print 'The MM_FILE environement variable is no set'; sys.exit(0)

DIAG_D = os.getenv('DIAG_D')
if DIAG_D == None: print 'The DIAG_D environement variable is no set'; sys.exit(0)


# NEMO variable names:
NN_T = os.getenv('NN_T')
if NN_T == None: print 'The NN_T environement variable is no set'; sys.exit(0)

NN_S = os.getenv('NN_S')
if NN_S == None: print 'The NN_S environement variable is no set'; sys.exit(0)



cname_script = basename(sys.argv[0])

print '\n'+cname_script


narg = len(sys.argv)
if narg < 2: print 'Usage: '+cname_script+' <year>'; sys.exit(0)
cy = sys.argv[1] ; jy=int(cy)





print '\n '+cname_script+':'
print ' CONF = '+CONF;
print ' RUN = '+RUN;
CONFRUN=CONF+'-'+RUN
print ' CONFRUN = '+CONFRUN
print ' CPREF = '+CPREF






# First will read name and coordinates of rectangular boxes to treat into file FILE_DEF_BOXES
##############################################################################################
vboxes, vi1, vj1, vi2, vj2 = bt.read_box_coordinates_in_ascii(FILE_DEF_BOXES)
nbb = len(vboxes)
print ''



# Opening mesh-mask and T-grid file of NEMO:
############################################

bt.chck4f(MM_FILE)
id_mm = Dataset(MM_FILE)
Xmask = id_mm.variables['tmask'][0,:,:,:]
ze1t  = id_mm.variables['e1t']    [0,:,:]
ze2t  = id_mm.variables['e2t']    [0,:,:]
id_mm.close()


cf_in = CPREF+cy+'0101_'+cy+'1231_grid_T.nc' ; bt.chck4f(cf_in)
id_in = Dataset(cf_in)
Vtime  = id_in.variables['time_counter'][:]
Vdepth = id_in.variables['deptht'][:]
Xlon   = id_in.variables['nav_lon'][:,:]
Xlat   = id_in.variables['nav_lat'][:,:]
Xtheta = id_in.variables[NN_T][:,:,:,:]
Xsali  = id_in.variables[NN_S][:,:,:,:]
print '(has ',Xtheta.shape[0],' time snapshots)\n'
id_in.close()


[ Nt, nk, nj, ni ] = Xtheta.shape

print 'Nt, nk, nj, ni =', Nt, nk, nj, ni


if Nt == 12:
    for jt in range(Nt): Vtime[jt] = float(jy) + (float(jt) + 0.5)/float(Nt)




# Loop along boxes:
for jb in range(nbb):

    cbox = vboxes[jb]

    i1 = vi1[jb]
    j1 = vj1[jb]
    i2 = vi2[jb]+1
    j2 = vj2[jb]+1
    
    print '\n *** Treating box '+cbox+' => ', i1, j1, i2-1, j2-1



    # Filling box arrays:
    # ~~~~~~~~~~~~~~~~~~~

    nx_b = i2 - i1
    ny_b = j2 - j1

    if l_plot_debug: print "nx_b , ny_b => ", nx_b , ny_b

    shape_array = [ Nt, nk, ny_b, nx_b ]

    # Temporary arrays:
    Ztmp = nmp.zeros(ny_b*nx_b) ;     Ztmp.shape = [ ny_b, nx_b ]
    Zar2 = nmp.zeros(ny_b*nx_b) ;     Zar2.shape = [ ny_b, nx_b ]    
    
    # Time series for each level of the mean sigma0 on the box:
    Zm1 = nmp.zeros(Nt*nk) ; Zm1.shape = [ Nt, nk ]
    Tm1 = nmp.zeros(Nt*nk) ; Tm1.shape = [ Nt, nk ]
    Sm1 = nmp.zeros(Nt*nk) ; Sm1.shape = [ Nt, nk ]


    
    Ztmp[:,:] = ze1t[j1:j2, i1:i2]*ze2t[j1:j2, i1:i2]
        
    for jk in range(nk):
        
        rarea = nmp.sum(Ztmp[:,:]*Xmask[jk, j1:j2, i1:i2])
        
        if rarea > 0.:
            Zar2[:,:] = Ztmp[:,:]*Xmask[jk, j1:j2, i1:i2]
            for jt in range(Nt):
                # Sigma0
                Zm1[jt,jk] = nmp.sum(bphy.sigma0(Xtheta[jt,jk, j1:j2, i1:i2],Xsali[jt,jk, j1:j2, i1:i2])*Zar2[:,:]) / rarea
                # Theta
                Tm1[jt,jk] = nmp.sum(Xtheta[jt,jk, j1:j2, i1:i2]*Zar2[:,:]) / rarea
                # Salinity
                Sm1[jt,jk] = nmp.sum( Xsali[jt,jk, j1:j2, i1:i2]*Zar2[:,:]) / rarea
        else:
            Zm1[:,jk] = nmp.nan ; Tm1[:,jk] = nmp.nan ; Sm1[:,jk] = nmp.nan




















    # Writing in output file
    ########################
    
    cf_out =  DIAG_D+'/TS_z_box_'+cbox+'_'+CONFRUN+'.nc'

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
        f_out.createDimension('deptht', nk)
    
        # Creating variables:
        id_t = f_out.createVariable('time','f4',('time',)) ;      id_t.units = 'year'
        id_z = f_out.createVariable('deptht','f4',('deptht',)) ;   id_z.units = 'm'

        #lili
        id_v01 = f_out.createVariable('sigma0',  'f4',('time','deptht',))
        id_v01.unit = ''; id_v01.long_name = 'sigma0 density averaged on box '+cbox
    
        id_v02 = f_out.createVariable('theta','f4',('time','deptht',))
        id_v02.unit = 'deg.C'; id_v02.long_name = 'potential temperature on box '+cbox

        id_v03 = f_out.createVariable('S','f4',('time','deptht',))
        id_v03.unit = 'PSU'; id_v03.long_name = 'salinity on box '+cbox
    

        

        id_z[:]    = Vdepth[:]
    
        for jm in range(Nt):
            id_t[jrec2write+jm]     = Vtime[jm]
            id_v01[jrec2write+jm,:] = Zm1[jm,:]
            id_v02[jrec2write+jm,:] = Tm1[jm,:]
            id_v03[jrec2write+jm,:] = Sm1[jm,:]


        f_out.box_coordinates = cbox+' => '+str(i1)+','+str(j1)+' -> '+str(i2-1)+','+str(j2-1)
        f_out.box_file        = FILE_DEF_BOXES
        f_out.Author          = 'L. Brodeau ('+cname_script+' of Barakuda)'
    
    else:
        vt  = f_out.variables['time']
        jrec2write = len(vt)
        v01 = f_out.variables['sigma0']
        v02 = f_out.variables['theta']
        v03 = f_out.variables['S']
    
        
        for jm in range(Nt):
            vt [jrec2write+jm]   = Vtime[jm]
            v01[jrec2write+jm,:] = Zm1[jm,:]
            v02[jrec2write+jm,:] = Tm1[jm,:]
            v03[jrec2write+jm,:] = Sm1[jm,:]
            
    f_out.close()
    
    print cf_out+' written!\n'
