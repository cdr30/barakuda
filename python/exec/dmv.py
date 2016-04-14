#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate 2D plots and maps of the Mixed layer depth
#
#       L. Brodeau, 2014
#
# Computes the DMV (Deep Mixed Volume) in user-defined deep-convective mixing
# regions (rectangle defined by 2 grid points) according to the method by
# Brodeau & Koenigk 2015.
#
# Rectangular boxes must be defined in a file which name is known as the global
# environment variable "FILE_DEF_BOXES"
#  => template: data/def_boxes_convection_ORCA1.txt
#
### Reference:
#  L. Brodeau and T. Koenigk. Extinction of the northern oceanic deep convection in
#  an ensemble of climate model simulations of the 20th and 21st centuries. Climate
#  Dynamics, Online, 2015.
#  DOI: 10.1007/s00382-015-2736-5
#
#
# List of environement variables that should be "known/set" when launching this script:
#
#  * ORCA   : global ORCA grid you use (ex: "ORCA1.L75")
#  * RUN    : name of your NEMO experiment
#  * DIAG_D : full path to the directory where the diagnostics (here a netcdf file) are saved
#  * FILE_DEF_BOXES: full path to the ASCII file containing definition of convection boxes
#                    => template: data/def_boxes_convection_ORCA1.txt
#  * MM_FILE : full path to the "mesh_mask.nc" file relevent to your ORCA !
#  * NN_MLD  : name of the Mixed-Layer Depth variable inside the *_grid_T.nc files of your NEMO experiment
#  * MLD_CRIT:  contains a list of all the "depth criterion to use" (wil be used for each box)
#              => ex:  MLD_CRIT="2000,1500,1000,725,500"
#
###########################################################################################

import barakuda_tool as bt
import barakuda_ncio as bnc


rdiv = 1000. ; # => result in 10^3 km^3   (DMV divided 4 times by rdiv)


# Getting all required environment variables needed inside dictionary vdic:
venv_needed = {'ORCA','RUN','DIAG_D','FILE_DEF_BOXES','MM_FILE','NN_MLD','MLD_CRIT'}
vdic = bt.check_env_var(sys.argv[0], venv_needed)


CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

cname_script = basename(sys.argv[0])

if len(sys.argv) != 3:
    print 'Usage : '+cname_script+' <ORCA1_RUN_grid_T.nc> <year>'
    sys.exit(0)
cf_in  = sys.argv[1]
cyear  = sys.argv[2] ; jyear = int(cyear); cyear = '%4.4i'%jyear

print 'Current year is '+cyear+' !\n'

# Vector containing the different z_crit:
vMLD_crit = []
vv = vdic['MLD_CRIT'].split(',')
for cv in vv: vMLD_crit.append(float(cv))
print "\n All the z_crit to use:", vMLD_crit[:]




# First will read name and coordinates of rectangular boxes to treat into file FILE_DEF_BOXES
##############################################################################################
vboxes, vi1, vj1, vi2, vj2 = bt.read_box_coordinates_in_ascii(vdic['FILE_DEF_BOXES'])
nbb = len(vboxes)
print ''


bt.chck4f(vdic['MM_FILE'])
id_mm = Dataset(vdic['MM_FILE'])
zmask_orca = id_mm.variables['tmask'][0,0,:,:]
ze1t_orca  = id_mm.variables['e1t']    [0,:,:]
ze2t_orca  = id_mm.variables['e2t']    [0,:,:]
id_mm.close()


bt.chck4f(cf_in)
id_in = Dataset(cf_in)
Xmld_orca  = id_in.variables[vdic['NN_MLD']][:,:,:]
print '(has ',Xmld_orca.shape[0],' time snapshots)\n'
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
        
        
    
    
        ##########################################
        # Writing/Appending in output netcdf file
        ##########################################

        # Appending only 1 record for 1 year into the netcdf file!
        
        cf_out =  vdic['DIAG_D']+'/DMV_'+czcrit+'_box_'+cbox+'_'+CONFRUN+'.nc'

        cv_dmv_m   = 'DMV_'+czcrit+'_'+ccold
        cv_dmv_jfm = 'DMV_'+czcrit+'_'+cvinter

        long_name1 = 'Deep Mixed Volume (crit = '+czcrit+'m) for '+ccold_ln+' on box '+cbox
        long_name2 = 'Deep Mixed Volume (crit = '+czcrit+'m) for '+cvinter+' on box '+cbox
        long_name3 = 'Deepest ML point in '+ccold_ln+' on box '+cbox
        long_name4 = 'Mean MLD in '+ccold_ln+' where MLD > '+czcrit+'m on box '+cbox


        bnc.wrt_appnd_1d_series([float(jy)], [VDMV[2]], cf_out, cv_dmv_m,  cu_t='year', cu_d='10^3 km^3', cln_d=long_name1,
                                vd2=[rc_WINT],       cvar2=cv_dmv_jfm,     cln_d2=long_name2,
                                vd3=[rML_max],       cvar3='ML_max',       cln_d3=long_name3,
                                vd4=[rML_deep_mean], cvar4='ML_deep_mean', cln_d4=long_name4)




    print '\n Z_crit '+str(int(czcrit))+'m done!\n\n'


print '\nBye!\n'

