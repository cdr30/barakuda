#!/usr/bin/env python

# L. Brodeau 2015

#####################################
import sys
import os
import numpy as nmp
from netCDF4 import Dataset
from os.path import basename

# Laurent's:
import barakuda_tool    as bt
import barakuda_physics as bph

#####################################

rmiss = -999.0

l_overflows = True


# 2 points at overflow for Denmark Strait:
PO_ds1 = [ 261, 244, 22 ] ; PO_ds2 = [ 261, 245, 22 ]

# 2 points at overflow for Faroe Bamks Channel:
PO_fb1 = [ 277, 236, 24 ] ; PO_fb2 = [ 277, 237, 24 ]


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

NN_SST = os.getenv('NN_SST')
if NN_SST == None:
    print 'The NN_SST environement variable is no set!'; sys.exit(0)

NN_SSS = os.getenv('NN_SSS')
if NN_SSS == None:
    print 'The NN_SSS environement variable is no set!'; sys.exit(0)

NN_MLD = os.getenv('NN_MLD')
if NN_MLD == None:
    print 'The NN_MLD environement variable is no set!'; sys.exit(0)

NN_S = os.getenv('NN_S')
if NN_S == None:
    print 'The NN_S environement variable is no set!'; sys.exit(0)

NN_T = os.getenv('NN_T')
if NN_T == None:
    print 'The NN_T environement variable is no set!'; sys.exit(0)

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




# First will read name and coordinates of rectangular boxes to treat into file FILE_DEF_BOXES
vboxes, vi1, vj1, vi2, vj2 = bt.read_box_coordinates_in_ascii(FILE_DEF_BOXES)
nbb = len(vboxes)
print ''






#[ nt, nj, ni ] = Xmld_orca.shape
#print 'nt, nj, ni =', nt, nj, ni




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
    

    # Filling Convection regions arrays:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    xsurf    = nmp.zeros(ny_b*nx_b) ;     xsurf.shape = [ ny_b, nx_b ]
    xtmp0    = nmp.zeros(ny_b*nx_b) ;    xtmp0.shape = [ ny_b, nx_b ]    
    xtmp1    = nmp.zeros(ny_b*nx_b) ;    xtmp1.shape = [ ny_b, nx_b ]
    xtmp2    = nmp.zeros(ny_b*nx_b) ;    xtmp2.shape = [ ny_b, nx_b ]


    bt.chck4f(cf_in)
    print ' Reading SST, SSS and MLD in box '+cbox+' in file '+cf_in
    id_in = Dataset(cf_in)
    if NN_SST == 'thetao':
        xsst  = id_in.variables[NN_SST][:,0,j1:j2,i1:i2]
    else:
        xsst  = id_in.variables[NN_SST][:,j1:j2,i1:i2]
    if NN_SSS == 'so':
        xsss  = id_in.variables[NN_SSS][:,0,j1:j2,i1:i2]
    else:
        xsss  = id_in.variables[NN_SSS][:,j1:j2,i1:i2]
        
    xmld  = id_in.variables[NN_MLD][:,j1:j2,i1:i2]
    xtpot = id_in.variables[NN_T][:,:,j1:j2,i1:i2]
    xsali = id_in.variables[NN_S][:,:,j1:j2,i1:i2]
    id_in.close()
    print ''

    if jb == 0: [ Nt, nj0, ni0 ] = xmld.shape

    bt.chck4f(cf_mesh_mask)
    id_mm = Dataset(cf_mesh_mask)
    Xmask = id_mm.variables['tmask'][0,:,j1:j2,i1:i2]
    ze1t  = id_mm.variables['e1t']    [0,j1:j2,i1:i2]
    ze2t  = id_mm.variables['e2t']    [0,j1:j2,i1:i2]
    if jb == 0:
        zgdept     = id_mm.variables['gdept_1d'][0,:]
        zgdepw     = id_mm.variables['gdepw_1d'][0,:]
        nk = len(zgdept)
    id_mm.close()

    xsurf[:,:] = ze1t[:,:]*ze2t[:,:]


    # Surface Sigma0
    xsg0 = bph.sigma0(xsst, xsss)

    # 3D sigma0
    Xsig0 = bph.sigma0(xtpot, xsali)

    # Depth of interest (between 300. and 2000.):
    #if jb == 0:
    #    jk1 = bt.find_index_from_value(300., zgdepw)
    #    jk2 = bt.find_index_from_value(2000., zgdepw)
    #    nb_zcrit = len(zgdepw[jk1:jk2])

    nb_zcrit = nk ; jk1 = 0 ; jk2 = nk-1

    rmean_sss0_deep_jfm  = nmp.zeros(nb_zcrit)
    rmean_sss0_deep_m03  = nmp.zeros(nb_zcrit)
    nbp_deeper_zcrit = nmp.zeros(nb_zcrit)


    # zcrit loop of w-depth:
    jk  = jk1 - 1
    for jcrit in range(nb_zcrit):

        jk = jk + 1
        
        zcrit = zgdepw[jk]

        #print '\n TESTING FOR zcrit=', str(zcrit)
        
        nb_month_to_count_jfm = 0
        nb_month_to_count_m03 = 0
        
        for jm in [ 0, 1, 2]:
    
            #print ' *** Month = ', jm+1
            #print ' Mean surface sigma0 in box '+cbox+':', nmp.mean(xsg0[jm,:,:])
        
            xtmp1[:,:] = xmld[jm,:,:]
            xtmp2[:,:] = xsg0[jm,:,:]
            
            ideep     = nmp.where(xtmp1[:,:] >= zcrit)
    
            [ n0, nbdeep ] = nmp.shape(ideep)
            #print '  Number of points where MLD > '+str(zcrit)+' =>', nbdeep
            nbp_deeper_zcrit[jcrit] = nbp_deeper_zcrit[jcrit] + float(nbdeep)

            # At least 1 point with deep mld:
            if nbdeep >= 1:
                nb_month_to_count_jfm = nb_month_to_count_jfm + 1
                
                xtmp0[:,:] = xsurf[:,:]*Xmask[0,:,:]
                xtmp1[:,:] = xtmp2[:,:]*xtmp0[:,:]                
                rmean_s0_deep = nmp.sum(xtmp1[ideep])/nmp.sum(xtmp0[ideep])
                
                rmean_sss0_deep_jfm[jcrit] = rmean_sss0_deep_jfm[jcrit] + rmean_s0_deep
                if jm == 2:
                    nb_month_to_count_m03 = 1
                    rmean_sss0_deep_m03[jcrit] = rmean_s0_deep
            else:
                rmean_s0_deep = 0.
            
            #print ' rmean_s0_deep => ', rmean_s0_deep
            #print ''
    


        #print str(nb_month_to_count_jfm)+' to count!!!'
        if nb_month_to_count_jfm > 0:
            rmean_sss0_deep_jfm[jcrit]  = rmean_sss0_deep_jfm[jcrit]/nb_month_to_count_jfm
            nbp_deeper_zcrit[jcrit] = nbp_deeper_zcrit[jcrit]/nb_month_to_count_jfm
        else:
            rmean_sss0_deep_jfm[jcrit] = rmiss

        # No convection in march:
        if nb_month_to_count_m03 == 0:
            rmean_sss0_deep_m03[jcrit] = rmiss

    # END / for jcrit in range(nb_zcrit)

    # Summary:    
    #for jk in range(jk1,jk2): print zgdepw[jk], rmean_sss0_deep_jfm[jk-jk1], nbp_deeper_zcrit[jk-jk1]





    # Mean annual value of the vertical profile of sigma0 in the box:
    ##################################################################

    vprof_sig0_ann = nmp.zeros(nk)
    vprof_sig0_jfm = nmp.zeros(nk)
    vprof_sig0_m03 = nmp.zeros(nk)    

    
    for jk in range(nk):

        xtmp0[:,:] = xsurf[:,:]*Xmask[jk,:,:]

        rsum = nmp.sum(xtmp0[:,:])

        if rsum > 0.: 
            vprof_sig0_ann[jk] = nmp.sum(nmp.mean(Xsig0[:,jk,:,:], axis=0)*xtmp0[:,:]) / rsum
            vprof_sig0_jfm[jk] = nmp.sum(nmp.mean(Xsig0[:3,jk,:,:],axis=0)*xtmp0[:,:]) / rsum
            vprof_sig0_m03[jk] = nmp.sum(         Xsig0[2,jk,:,:]         *xtmp0[:,:]) / rsum
        else:
            vprof_sig0_ann[jk] = rmiss
            vprof_sig0_jfm[jk] = rmiss
            vprof_sig0_m03[jk] = rmiss
        












    ########################
    # Writing in output file
    ########################

    cf_out =  DIAG_D+'/zcrit_box_'+cbox+'_'+CONFRUN+'.nc'

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
        f_out.createDimension('depthw', nk)
        f_out.createDimension('deptht', nk)
        
    
        # Creating variables:
        id_t = f_out.createVariable('time','f4',('time',)) ;  id_t.units = 'year'
        id_zw = f_out.createVariable('depthw','f4',('depthw',)) ;   id_zw.units = 'm'
        id_zt = f_out.createVariable('deptht','f4',('deptht',)) ;   id_zt.units = 'm'

        id_v01 = f_out.createVariable('SSsig0_jfm',  'f4',('time','depthw',), fill_value=rmiss)
        id_v01.unit = ''; id_v01.long_name = 'Mean winter sea surface sigma0 on the area where MLD>=depthw(depthw) ('+cbox+')'
        
        id_v02 = f_out.createVariable('SSsig0_m03',  'f4',('time','depthw',), fill_value=rmiss)
        id_v02.unit = ''; id_v02.long_name = 'Mean March sea surface sigma0 on the area where MLD>=depthw(depthw) ('+cbox+')'
        
        id_v03 = f_out.createVariable('Nbp_w_deep',  'f4',('time','depthw',), fill_value=rmiss)
        id_v03.unit = ''; id_v03.long_name = 'Mean winter average of the number of points where MLD>=depthw(depthw) ('+cbox+')'

        id_v04 = f_out.createVariable('sig0_ann',  'f4',('time','deptht',), fill_value=rmiss)
        id_v04.unit = ''; id_v04.long_name = 'Mean annual Sigma0 horizontally averaged on box '+cbox
        id_v05 = f_out.createVariable('sig0_jfm',  'f4',('time','deptht',), fill_value=rmiss)
        id_v05.unit = ''; id_v05.long_name = 'Mean JFM Sigma0 horizontally averaged on box '+cbox
        id_v06 = f_out.createVariable('sig0_m03',  'f4',('time','deptht',), fill_value=rmiss)
        id_v06.unit = ''; id_v06.long_name = 'Mean March Sigma0 horizontally averaged on box '+cbox





        id_t[jrec2write]     = float(jy)
        id_zw[:]              = zgdepw[:]
        id_zt[:]              = zgdept[:]        
        id_v01[jrec2write,:] = rmean_sss0_deep_jfm[:]
        id_v02[jrec2write,:] = rmean_sss0_deep_m03[:]
        id_v03[jrec2write,:] = nbp_deeper_zcrit[:]
        id_v04[jrec2write,:] = vprof_sig0_ann[:]
        id_v05[jrec2write,:] = vprof_sig0_jfm[:]
        id_v06[jrec2write,:] = vprof_sig0_m03[:]        
        
        f_out.box_coordinates = cbox+' => '+str(i1)+','+str(j1)+' -> '+str(i2-1)+','+str(j2-1)
        f_out.box_file        = FILE_DEF_BOXES
        f_out.Author          = 'L. Brodeau ('+cname_script+' of Barakuda)'
    
    else:
        vt  = f_out.variables['time']
        jrec2write = len(vt)
        v01 = f_out.variables['SSsig0_jfm']
        v02 = f_out.variables['SSsig0_m03']
        v03 = f_out.variables['Nbp_w_deep']
        v04 = f_out.variables['sig0_ann']
        v05 = f_out.variables['sig0_jfm']
        v06 = f_out.variables['sig0_m03']
            
        vt [jrec2write]   = float(jy)
        v01[jrec2write,:] = rmean_sss0_deep_jfm[:]
        v02[jrec2write,:] = rmean_sss0_deep_m03[:]
        v03[jrec2write,:] = nbp_deeper_zcrit[:]
        v04[jrec2write,:] = vprof_sig0_ann[:]
        v05[jrec2write,:] = vprof_sig0_jfm[:]
        v06[jrec2write,:] = vprof_sig0_m03[:]        
    f_out.close()
    
    print cf_out+' written!\n'






if l_overflows:


    
    id_mm = Dataset(cf_mesh_mask)
    Xmsk = id_mm.variables['tmask'][0,:,:,:]
    id_mm.close()


    # Water properties at overflows:

    id_in = Dataset(cf_in)

    ip1 = PO_ds1[0] ; jp1 = PO_ds1[1] ; kp1 = PO_ds1[2]
    ip2 = PO_ds2[0] ; jp2 = PO_ds2[1] ; kp2 = PO_ds2[2]
    print '\n Denmark Strait:', Xmsk[kp1,jp1,ip1], Xmsk[kp2,jp2,ip2], zgdept[kp1], zgdepw[kp1+1]

    vt_p1 = id_in.variables[NN_T][:,kp1,jp1,ip1] ; vs_p1 = id_in.variables[NN_S][:,kp1,jp1,ip1]
    vz_p1 = bph.sigma0(vt_p1,vs_p1)
    t_p1 = nmp.mean(vt_p1) ; s_p1 = nmp.mean(vs_p1) ; z_p1 = nmp.mean(vz_p1)

    vt_p2 = id_in.variables[NN_T][:,kp2,jp2,ip2] ; vs_p2 = id_in.variables[NN_S][:,kp2,jp2,ip2]
    vz_p2 = bph.sigma0(vt_p2,vs_p2)
    t_p2 = nmp.mean(vt_p2) ; s_p2 = nmp.mean(vs_p2) ; z_p2 = nmp.mean(vz_p2)

    t_of_DS = 0.5*(t_p1 + t_p2) ; s_of_DS = 0.5*(s_p1 + s_p2) ; z_of_DS = 0.5*(z_p1 + z_p2)


    ip1 = PO_fb1[0] ; jp1 = PO_fb1[1] ; kp1 = PO_fb1[2]
    ip2 = PO_fb2[0] ; jp2 = PO_fb2[1] ; kp2 = PO_fb2[2]
    print '\n Faroe Banks Channel:', Xmsk[kp1,jp1,ip1], Xmsk[kp2,jp2,ip2], zgdept[kp1], zgdepw[kp1+1]
    
    vt_p1 = id_in.variables[NN_T][:,kp1,jp1,ip1] ; vs_p1 = id_in.variables[NN_S][:,kp1,jp1,ip1]
    vz_p1 = bph.sigma0(vt_p1,vs_p1)
    t_p1 = nmp.mean(vt_p1) ; s_p1 = nmp.mean(vs_p1) ; z_p1 = nmp.mean(vz_p1)

    vt_p2 = id_in.variables[NN_T][:,kp2,jp2,ip2] ; vs_p2 = id_in.variables[NN_S][:,kp2,jp2,ip2]
    vz_p2 = bph.sigma0(vt_p2,vs_p2)
    t_p2 = nmp.mean(vt_p2) ; s_p2 = nmp.mean(vs_p2) ; z_p2 = nmp.mean(vz_p2)

    t_of_FB = 0.5*(t_p1 + t_p2) ; s_of_FB = 0.5*(s_p1 + s_p2) ; z_of_FB = 0.5*(z_p1 + z_p2)

    print ''

    id_in.close()







    ########################
    # Writing in output file
    ########################

    cf_out =  DIAG_D+'/overflows_properties_'+CONFRUN+'.nc'

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
        id_t = f_out.createVariable('time','f4',('time',)) ;  id_t.units = 'year'

        id_v05 = f_out.createVariable('t_of_DS',  'f4',('time',), fill_value=rmiss)
        id_v05.unit = ''; id_v05.long_name = 'Temperature at Denmark Strait overflow '
        id_v06 = f_out.createVariable('s_of_DS',  'f4',('time',), fill_value=rmiss)
        id_v06.unit = ''; id_v06.long_name = 'Salinity at Denmark Strait overflow '
        id_v07 = f_out.createVariable('z_of_DS',  'f4',('time',), fill_value=rmiss)
        id_v07.unit = ''; id_v07.long_name = 'Sigma0 at Denmark Strait overflow '

        id_v15 = f_out.createVariable('t_of_FB',  'f4',('time',), fill_value=rmiss)
        id_v15.unit = ''; id_v15.long_name = 'Temperature at Faroe Banks Channel overflow '
        id_v16 = f_out.createVariable('s_of_FB',  'f4',('time',), fill_value=rmiss)
        id_v16.unit = ''; id_v16.long_name = 'Salinity at Faroe Banks Channel overflow '
        id_v17 = f_out.createVariable('z_of_FB',  'f4',('time',), fill_value=rmiss)
        id_v17.unit = ''; id_v17.long_name = 'Sigma0 at Faroe Banks Channel overflow '


        id_t[jrec2write]     = float(jy)
        id_v05[jrec2write] = t_of_DS ; id_v06[jrec2write] = s_of_DS ; id_v07[jrec2write] = z_of_DS
        id_v15[jrec2write] = t_of_FB ; id_v16[jrec2write] = s_of_FB ; id_v17[jrec2write] = z_of_FB
                    
        f_out.Author          = 'L. Brodeau ('+cname_script+' of Barakuda)'
    
    else:
        vt  = f_out.variables['time']
        jrec2write = len(vt)
        r05 = f_out.variables['t_of_DS'] ; r06 = f_out.variables['s_of_DS'] ; r07 = f_out.variables['z_of_DS']
        r15 = f_out.variables['t_of_FB'] ; r16 = f_out.variables['s_of_FB'] ; r17 = f_out.variables['z_of_FB']
            
        vt [jrec2write]   = float(jy)
        r05[jrec2write] = t_of_DS ; r06[jrec2write] = s_of_DS ; r07[jrec2write] = z_of_DS
        r15[jrec2write] = t_of_FB ; r16[jrec2write] = s_of_FB ; r17[jrec2write] = z_of_FB 
        
    f_out.close()
    
    print cf_out+' written!\n'













print '\nBye!\n'












