#!/usr/bin/python

# L. Brodeau 2015

#####################################
import sys
import os
import numpy as nmp
from netCDF4 import Dataset
from os.path import basename

# Laurent's:
import barakuda_plot as bp
import barakuda_tool as bt

#####################################


#l_plot_debug = True
l_plot_debug = False

rdiv = 1000. ; # => result in 10^3 km^3   (DMV divided 4 times by rdiv)

FILE_DEF_BOXES = os.getenv('FILE_DEF_BOXES')
if FILE_DEF_BOXES == None: print 'The FILE_DEF_BOXES environement variable is no set'; sys.exit(0)

CONF = os.getenv('CONF')
if CONF == None: print 'The CONF environement variable is no set'; sys.exit(0)

RUN = os.getenv('RUN')
if RUN == None: print 'The RUN environement variable is no set'; sys.exit(0)
CONFRUN = CONF+'-'+RUN

CPREF = os.getenv('CPREF')
if CPREF == None: print 'The CPREF environement variable is no set'; sys.exit(0)

MLD_CRIT = os.getenv('MLD_CRIT')
if MLD_CRIT == None:
    print 'WARNING: the MLD_CRIT environement variable is no set!';
    print '   => setting to 1000m!'; MLD_CRIT='1000,500'

NN_MLD = os.getenv('NN_MLD')
if NN_MLD == None:
    print 'The NN_MLD environement variable is no set!'; sys.exit(0)

DIAG_D = os.getenv('DIAG_D')
if DIAG_D == None: print 'The DIAG_D environement variable is no set'; sys.exit(0)

I_MLD_FIG = os.getenv('I_MLD_FIG')
if I_MLD_FIG == None:
    print 'WARNING: The I_MLD_FIG environement variable was no set, forcing to 0 ...\n'
    i_figures = 0
else:
    i_figures = int(I_MLD_FIG)


print '\n '+sys.argv[0]+':'
print ' CONF = '+CONF;
print ' RUN = '+RUN; print ' CONFRUN = '+CONFRUN; print ' CPREF = '+CPREF
print ' i_figures = ', i_figures
print ' NN_MLD = ', NN_MLD,'\n'

path_fig='./'

# Image type? eps, png, jpg...
FIG_FORM = os.getenv('FIG_FORM')
if FIG_FORM == None: FIG_FORM='png'


cname_script = basename(sys.argv[0])
print '\n'+cname_script


narg = len(sys.argv)
if narg < 2: print 'Usage: '+sys.argv[0]+' <year>'; sys.exit(0)
cy = sys.argv[1] ; jy=int(cy)

cf_in = CPREF+cy+'0101_'+cy+'1231_grid_T.nc'

cf_mesh_mask = './mesh_mask.nc'



# Vector containing the different z_crit:
vMLD_crit = []
vv = MLD_CRIT.split(',')
for cv in vv: vMLD_crit.append(float(cv))
print "\n All the z_crit to use:", vMLD_crit[:]







# First will read name and coordinates of rectangular boxes to treat into file FILE_DEF_BOXES
##############################################################################################
vboxes, vi1, vj1, vi2, vj2 = bt.read_box_coordinates_in_ascii(FILE_DEF_BOXES)
nbb = len(vboxes)
print ''


bt.chck4f(cf_mesh_mask)
id_mm = Dataset(cf_mesh_mask)
zmask_orca = id_mm.variables['tmask'][0,0,:,:]
ze1t_orca  = id_mm.variables['e1t']    [0,:,:]
ze2t_orca  = id_mm.variables['e2t']    [0,:,:]
id_mm.close()


bt.chck4f(cf_in)
id_in = Dataset(cf_in)
Xmld_orca  = id_in.variables[NN_MLD][:,:,:]
print '(has ',Xmld_orca.shape[0],' time snapshots)\n'
#vdepth = id_in.variables['deptht'][:]
xlat   = id_in.variables['nav_lat'][:,:]
id_in.close()


[ nt, nj, ni ] = Xmld_orca.shape
print 'nt, nj, ni =', nt, nj, ni





# Loop along depth crtiteria:

for rMLD_crit in vMLD_crit:

    czcrit = str(int(rMLD_crit))

    # Loop along convection boxes:
    for jb in range(nbb):
    
        cbox = vboxes[jb]
    
        i1 = vi1[jb]
        j1 = vj1[jb]
        i2 = vi2[jb]+1
        j2 = vj2[jb]+1
    
        nx_b = i2 - i1
        ny_b = j2 - j1
        
        print '\n *** Treating box '+cbox+' => ', i1, j1, i2-1, j2-1


        
        
        # NH
        icold  = 2 ; # march
        ccold  = 'march'
        ccold_ln = 'March'
        vinter = [0, 1, 2]
        cvinter = 'JFM'
    
        if xlat[j1,i1] < 0:
            # SH
            icold  = 8 ; # september
            ccold  = 'sept'
            ccold_ln = 'September'
            vinter = [6, 7, 8]
            cvinter = 'JAS'
    
    
        # Filling Convection regions arrays:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        xmld     = nmp.zeros((nt, ny_b, nx_b))
        msk_deep = nmp.zeros((nt, ny_b, nx_b))
        xe1t     = nmp.zeros((ny_b, nx_b))
        xe2t     = nmp.zeros((ny_b, nx_b))
        xtmp     = nmp.zeros((ny_b, nx_b))
        VDMV     = nmp.zeros(3)
        
        # MLD:
        xmld[:,:,:] = Xmld_orca[:,j1:j2,i1:i2]
    
        
        #E1Tb & E2T:
        xe1t[:,:] = ze1t_orca[j1:j2,i1:i2] / rdiv
        xe2t[:,:] = ze2t_orca[j1:j2,i1:i2] / rdiv
        
    
        # Building deep ML mask for the 3 first months:
        for jm in vinter:
            xtmp[:,:] = 0.
            xtmp[:,:] = zmask_orca[j1:j2,i1:i2]
            
            # Excluding points where MLD < rMLD_crit
            idx1 = nmp.where(xmld[jm,:,:] < rMLD_crit)
            xtmp[idx1] = 0
            msk_deep[jm,:,:] = xtmp[:,:]
        
        if l_plot_debug:
            bp.check_with_fig_2(xmld[icold,:,:], zmask_orca[j1:j2,i1:i2], cbox+'_xmld')
            bp.check_with_fig_2(xmld[icold,:,:], msk_deep[icold,:,:], cbox+'_xmld_m')
        
    
        xtmp[:,:] = xe1t[:,:]*xe2t[:,:]
    
    
        # Deepest ML in march:
        rML_max       = nmp.max(xmld[icold,:,:])
    
        # Mean ML in march where ML > zcrit:
        rd = nmp.sum(xtmp[:,:]*msk_deep[icold,:,:])
        if rd > 0:
            rML_deep_mean = nmp.sum(xmld[icold,:,:]*xtmp[:,:]*msk_deep[icold,:,:])/rd
        else:
            rML_deep_mean = 0.
    
    
        # DMV
        #####

        jc=-1
        for jm in vinter:
            jc=jc+1
            VDMV[jc] = nmp.sum((xmld[jm,:,:]-rMLD_crit)*xtmp[:,:]*msk_deep[jm,:,:])/(rdiv*rdiv)
            
        
        rc_WINT = nmp.mean(VDMV)
        
        
        if l_plot_debug: print "VDMV[2] = ", VDMV[2]
    
    
    
        ########################
        # Writing in output file
        ########################
    
        cv_dmv_m   = 'DMV_'+czcrit+'_'+ccold
        cv_dmv_jfm = 'DMV_'+czcrit+'_'+cvinter
        
        cf_out =  DIAG_D+'/DMV_'+czcrit+'_box_'+cbox+'_'+CONFRUN+'.nc'
    
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
        
            # Creating variables:
            id_t = f_out.createVariable('time','f4',('time',)) ;      id_t.units = 'year'
    
            id_v01 = f_out.createVariable(cv_dmv_m,  'f4',('time',))
            id_v01.unit = '10^3 km^3'; id_v01.long_name = 'Deep Mixed Volume (crit = '+czcrit+'m) for '+ccold_ln+' on box '+cbox
            id_v02 = f_out.createVariable(cv_dmv_jfm,  'f4',('time',))
            id_v02.unit = '10^3 km^3'; id_v02.long_name = 'Deep Mixed Volume (crit = '+czcrit+'m) for '+cvinter+' on box '+cbox
    
            id_v03 = f_out.createVariable('ML_max',  'f4',('time',))
            id_v03.unit = '10^3 km^3'; id_v03.long_name = 'Deepest ML point in '+ccold_ln+' on box '+cbox
    
            id_v04 = f_out.createVariable('ML_deep_mean',  'f4',('time',))
            id_v04.unit = '10^3 km^3'; id_v04.long_name = 'Mean MLD in '+ccold_ln+' where MLD > '+czcrit+'m on box '+cbox
    
            id_t[jrec2write]   = float(jy)
            id_v01[jrec2write] = VDMV[2]
            id_v02[jrec2write] = rc_WINT        
            id_v03[jrec2write] = rML_max
            id_v04[jrec2write] = rML_deep_mean
            
            f_out.box_coordinates = cbox+' => '+str(i1)+','+str(j1)+' -> '+str(i2-1)+','+str(j2-1)
            f_out.box_file        = FILE_DEF_BOXES
            f_out.Author          = 'L. Brodeau ('+cname_script+' of Barakuda)'
        
        else:
            vt  = f_out.variables['time']
            jrec2write = len(vt)
            v01 = f_out.variables[cv_dmv_m]
            v02 = f_out.variables[cv_dmv_jfm]
            v03 = f_out.variables['ML_max']
            v04 = f_out.variables['ML_deep_mean']
            
            vt [jrec2write] = float(jy)
            v01[jrec2write] = VDMV[2]
            v02[jrec2write] = rc_WINT
            v03[jrec2write] = rML_max
            v04[jrec2write] = rML_deep_mean
                
        f_out.close()
        
        print cf_out+' written!\n'

    print '\n Z_crit '+str(int(czcrit))+'m done!\n\n'


print '\nBye!\n'

