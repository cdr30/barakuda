#!/home/x_laubr/bin/python

# L. Brodeau, december 2011

import sys
import os
import numpy as nmp
from netCDF4 import Dataset

import lb_util_orca as luo
import lb_manip  as lbm


ORCA = os.getenv('ORCA')
if ORCA == None: print 'The ORCA environement variable is no set'; sys.exit(0)
CASE = os.getenv('CASE')
if CASE == None: print 'The CASE environement variable is no set'; sys.exit(0)
NAME = os.getenv('NAME')
if NAME == None: print 'The NAME environement variable is no set'; sys.exit(0)
STORE_DIR = os.getenv('STORE_DIR')
if STORE_DIR == None: print 'The STORE_DIR environement variable is no set'; sys.exit(0)

print '\n surf_fluxes.py:'
print ' ORCA = '+ORCA;
print ' CASE = '+CASE; print ' NAME = '+NAME; print ' STORE_DIR = '+STORE_DIR


if 'ORCA2' in ORCA:
    ji_lat0 = 132
elif 'ORCA1' in ORCA:
    ji_lat0 = 265
else:
    print 'FIX ME!!! ssh.py => dont know ji_lat0 for conf '+ORCA+' !!!'; sys.exit(0)




path_fig='./'
fig_type='png'

# Mesh-mask file:
cf_mesh_mask = os.getenv('MM_FILE')
if cf_mesh_mask == None: print 'The MM_FILE environement variable (mesh_mask) is no set'; sys.exit(0)
print '\n Mesh-Mask file is:\n', cf_mesh_mask, '\n'



narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <year1> <year2>'; sys.exit(0)
cy1 = sys.argv[1] ; cy2=sys.argv[2]; jy1=int(cy1); jy2=int(cy2)

print ' First and last year to treat:', jy1, jy2, '\n'





#     C L I M A T O L O G Y
#     *********************
# NOCS2.0 climatology: monthly on ORCA1 grid...
f_nocs_sw  = STORE_DIR+'/CLIM-LB/Qsw_1degx1deg-'+ORCA+'_NOCS2.0.nc'   ; cv_nocs_sw ='Qsw'
f_nocs_lw  = STORE_DIR+'/CLIM-LB/Qlw_1degx1deg-'+ORCA+'_NOCS2.0.nc'   ; cv_nocs_lw ='Qlw'
f_nocs_la  = STORE_DIR+'/CLIM-LB/Qlat_1degx1deg-'+ORCA+'_NOCS2.0.nc'  ; cv_nocs_la ='Qlat'
f_nocs_se  = STORE_DIR+'/CLIM-LB/Qsen_1degx1deg-'+ORCA+'_NOCS2.0.nc'  ; cv_nocs_se ='Qsen'


lbm.chck4f(f_nocs_sw)
id_nocs_sw = Dataset(f_nocs_sw)
XSW_OBS    = id_nocs_sw.variables[cv_nocs_sw][:,:,:]; print ' has ', XSW_OBS.shape[0],' time snapshots\n'
id_nocs_sw.close()

lbm.chck4f(f_nocs_lw)
id_nocs_lw = Dataset(f_nocs_lw)
XLW_OBS    = id_nocs_lw.variables[cv_nocs_lw][:,:,:]; print ' has ', XLW_OBS.shape[:],' time snapshots\n'
id_nocs_lw.close()

lbm.chck4f(f_nocs_la)
id_nocs_la = Dataset(f_nocs_la)
XLA_OBS    = id_nocs_la.variables[cv_nocs_la][:,:,:]; print ' has ', XLA_OBS.shape[:],' time snapshots\n'
id_nocs_la.close()

lbm.chck4f(f_nocs_se)
id_nocs_se = Dataset(f_nocs_se)
XSE_OBS    = id_nocs_se.variables[cv_nocs_se][:,:,:]; print ' has ', XSE_OBS.shape[:],' time snapshots\n'
id_nocs_se.close()

[ nmn , nj0 , ni0 ] = XSW_OBS.shape


# Getting coordinates:
cf_coor = STORE_DIR+'/'+ORCA+'-I/coordinates.nc' ; print 'Coordinates =', cf_coor
lbm.chck4f(cf_coor) ; id_coor = Dataset(cf_coor)
xlon   = id_coor.variables['nav_lon'][:,:] ; xlat   = id_coor.variables['nav_lat'][:,:]
id_coor.close()




print '\n\n\n Treating run', NAME ; print '============================\n'

cd_str = STORE_DIR+'/diagnoce/'+CASE

# Getting land-sea mask:
lbm.chck4f(cf_mesh_mask)
id_mask = Dataset(cf_mesh_mask) ; imask  = id_mask.variables['tmask'][0,:,:,:]; id_mask.close()





# Getting NEMO mean seasonal climatology of temperature and salinity:
# ------------------------------------------------------------------
cf_nemo = cd_str+'/clim/mclim_'+NAME+'_'+cy1+'-'+cy2+'_grid_T.nc4' ; lbm.chck4f(cf_nemo)

id_mod = Dataset(cf_nemo)
XSW_MOD = id_mod.variables['soshfldo'][:,:,:] ; print '(has ',XSW_MOD.shape[0],' time snapshots)\n'
XLW_MOD = id_mod.variables['qlw_oce'] [:,:,:]
XLA_MOD = id_mod.variables['qla_oce'] [:,:,:]
XSE_MOD = id_mod.variables['qsb_oce'] [:,:,:]
id_mod.close()


# Testing if grids agree...
[ nt, nj, ni ] = XSW_OBS.shape
#if nj != nj0 or ni != ni0:
#    print 'ERROR: 3D clim and NEMO file do no agree in shape!'
#    print '       clim => '+str(ni0)+', '+str(nj0)+', '+str(nk0)
#    print '       NEMO => '+str(ni)+', '+str(nj)+', '+str(nk)



    
# Creating Annual Seasonal climatologies
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
XOBS = nmp.zeros(5*3*nj*ni) ; XOBS.shape = [ 5, 3, nj, ni ]
XMOD = nmp.zeros(5*3*nj*ni) ; XMOD.shape = [ 5, 3, nj, ni ]

# Qnet => [0,* ; Qsw => [1,* ; Qlw => [2,* ; Qlat => [3,* ; Qsen => [4,*


# Annual => index 0
XOBS[1,0,:,:] = nmp.mean(XSW_OBS, axis=0)
XOBS[2,0,:,:] = nmp.mean(XLW_OBS, axis=0 )
XOBS[3,0,:,:] = nmp.mean(XLA_OBS, axis=0 )
XOBS[4,0,:,:] = nmp.mean(XSE_OBS, axis=0 )

XMOD[1,0,:,:] = nmp.mean(XSW_MOD, axis=0)
XMOD[2,0,:,:] = nmp.mean(XLW_MOD, axis=0 )
XMOD[3,0,:,:] = nmp.mean(XLA_MOD, axis=0 )
XMOD[4,0,:,:] = nmp.mean(XSE_MOD, axis=0 )


# DJF => index 1
XOBS[1,1,:,:] = 1./3. * ( XSW_OBS[0,:,:] + XSW_OBS[1,:,:] + XSW_OBS[11,:,:] )
XOBS[2,1,:,:] = 1./3. * ( XLW_OBS[0,:,:] + XLW_OBS[1,:,:] + XLW_OBS[11,:,:] )
XOBS[3,1,:,:] = 1./3. * ( XLA_OBS[0,:,:] + XLA_OBS[1,:,:] + XLA_OBS[11,:,:] )
XOBS[4,1,:,:] = 1./3. * ( XSE_OBS[0,:,:] + XSE_OBS[1,:,:] + XSE_OBS[11,:,:] )

XMOD[1,1,:,:] = 1./3. * ( XSW_MOD[0,:,:] + XSW_MOD[1,:,:] + XSW_MOD[11,:,:] )
XMOD[2,1,:,:] = 1./3. * ( XLW_MOD[0,:,:] + XLW_MOD[1,:,:] + XLW_MOD[11,:,:] )
XMOD[3,1,:,:] = 1./3. * ( XLA_MOD[0,:,:] + XLA_MOD[1,:,:] + XLA_MOD[11,:,:] )
XMOD[4,1,:,:] = 1./3. * ( XSE_MOD[0,:,:] + XSE_MOD[1,:,:] + XSE_MOD[11,:,:] )

# JJA => index 2
XOBS[1,2,:,:] = 1./3. * ( XSW_OBS[5,:,:] + XSW_OBS[6,:,:] + XSW_OBS[7,:,:] )
XOBS[2,2,:,:] = 1./3. * ( XLW_OBS[5,:,:] + XLW_OBS[6,:,:] + XLW_OBS[7,:,:] )
XOBS[3,2,:,:] = 1./3. * ( XLA_OBS[5,:,:] + XLA_OBS[6,:,:] + XLA_OBS[7,:,:] )
XOBS[4,2,:,:] = 1./3. * ( XSE_OBS[5,:,:] + XSE_OBS[6,:,:] + XSE_OBS[7,:,:] )

XMOD[1,2,:,:] = 1./3. * ( XSW_MOD[5,:,:] + XSW_MOD[6,:,:] + XSW_MOD[7,:,:] )
XMOD[2,2,:,:] = 1./3. * ( XLW_MOD[5,:,:] + XLW_MOD[6,:,:] + XLW_MOD[7,:,:] )
XMOD[3,2,:,:] = 1./3. * ( XLA_MOD[5,:,:] + XLA_MOD[6,:,:] + XLA_MOD[7,:,:] )
XMOD[4,2,:,:] = 1./3. * ( XSE_MOD[5,:,:] + XSE_MOD[6,:,:] + XSE_MOD[7,:,:] )


# Creating 1D long. and lat.:
vlon = nmp.zeros(ni) ;  vlon.shape = [ ni ] ; vlon[:] = xlon[0,:]
vlat = nmp.zeros(nj) ;  vlat.shape = [ nj ] ; vlat[:] = xlat[:,ji_lat0]




#          Qnet   Qsw    Qlw     Qla     Qse
var_max = [ 0.,    320.,   -30.,     0.,    20. ]
var_min = [ 0.,      0.,   -80.,  -300.,  -100. ]
var_cnt = [ 0.,     10.,    2.5,    10.,     5. ]
var_pal = [ '',   'sst','sst_inv','sst_inv','sst_inv' ]

jv=1
for cv in [ 'Qsw' , 'Qlw' , 'Qla' , 'Qse' ]:
    print ''; print ''; print ''; print cv; print '  => # ', jv; print ''


    # A N N U A L
    ##############

    # Ploting annual obs. for Qsw:
    luo.plot_2d(vlon, vlat, XOBS[jv,0,:,:], imask[0,:,:], var_min[jv], var_max[jv], var_cnt[jv],
               corca=ORCA, lkcont=False, cpal=var_pal[jv], cfignm=path_fig+cv+'_obs_annual_'+NAME,
               cbunit='(W/m^2)', cfig_type=fig_type, lat_min=-58., lat_max=65., 
               ctitle='Annual mean obs. : '+cv+' (NOCS2.0)',
               lforce_lim=True)

    # Ploting annual NEMO for Qsw:
    luo.plot_2d(vlon, vlat, XMOD[jv,0,:,:], imask[0,:,:], var_min[jv], var_max[jv], var_cnt[jv],
               corca=ORCA, lkcont=False, cpal=var_pal[jv], cfignm=path_fig+cv+'_mod_annual_'+NAME,
               cbunit='(W/m^2)', cfig_type=fig_type, lat_min=-58., lat_max=65., 
               ctitle='Annual mean NEMO : '+cv,
               lforce_lim=True)
    
    # Ploting diff NEMO - obs.:
    luo.plot_2d(vlon, vlat, XMOD[jv,0,:,:]-XOBS[jv,0,:,:], imask[0,:,:], -30., 30., 2.,
               corca=ORCA, lkcont=False, cpal='bbr', cfignm=path_fig+cv+'_diff_annual_'+NAME,
               cbunit='(W/m^2)', cfig_type=fig_type, lat_min=-58., lat_max=65.,  
               ctitle='Annual mean NEMO - obs : '+cv,
               lforce_lim=True)


    # W I N T E R
    #############

    # Ploting DJF obs. for Qsw:
    luo.plot_2d(vlon, vlat, XOBS[jv,1,:,:], imask[0,:,:], var_min[jv], var_max[jv], var_cnt[jv],
               corca=ORCA, lkcont=False, cpal=var_pal[jv], cfignm=path_fig+cv+'_obs_DJF_'+NAME,
               cbunit='(W/m^2)', cfig_type=fig_type, lat_min=-58., lat_max=65., 
               ctitle='DJF obs. : '+cv+' (NOCS2.0)',
               lforce_lim=True)

    # Ploting DJF NEMO for Qsw:
    luo.plot_2d(vlon, vlat, XMOD[jv,1,:,:], imask[0,:,:], var_min[jv], var_max[jv], var_cnt[jv],
               corca=ORCA, lkcont=False, cpal=var_pal[jv], cfignm=path_fig+cv+'_mod_DJF_'+NAME,
               cbunit='(W/m^2)', cfig_type=fig_type, lat_min=-58., lat_max=65., 
               ctitle='DJF NEMO : '+cv,
               lforce_lim=True)
    
    # Ploting diff NEMO - obs.:
    luo.plot_2d(vlon, vlat, XMOD[jv,1,:,:]-XOBS[jv,1,:,:], imask[0,:,:], -30., 30., 2.,
               corca=ORCA, lkcont=False, cpal='bbr', cfignm=path_fig+cv+'_diff_DJF_'+NAME,
               cbunit='(W/m^2)', cfig_type=fig_type, lat_min=-58., lat_max=65.,  
               ctitle='DJF NEMO - obs : '+cv,
               lforce_lim=True)
    



    # S U M M E R
    ##############

    # Ploting JJA obs. for Qsw:
    luo.plot_2d(vlon, vlat, XOBS[jv,2,:,:], imask[0,:,:], var_min[jv], var_max[jv], var_cnt[jv],
               corca=ORCA, lkcont=False, cpal=var_pal[jv], cfignm=path_fig+cv+'_obs_JJA_'+NAME,
               cbunit='(W/m^2)', cfig_type=fig_type, lat_min=-58., lat_max=65., 
               ctitle='JJA obs. : '+cv+' (NOCS2.0)',
               lforce_lim=True)
    
    # Ploting JJA NEMO for Qsw:
    luo.plot_2d(vlon, vlat, XMOD[jv,2,:,:], imask[0,:,:], var_min[jv], var_max[jv], var_cnt[jv],
               corca=ORCA, lkcont=False, cpal=var_pal[jv], cfignm=path_fig+cv+'_mod_JJA_'+NAME,
               cbunit='(W/m^2)', cfig_type=fig_type, lat_min=-58., lat_max=65., 
               ctitle='JJA NEMO : '+cv,
               lforce_lim=True)

   # Ploting diff NEMO - obs.:
    luo.plot_2d(vlon, vlat, XMOD[jv,2,:,:]-XOBS[jv,2,:,:], imask[0,:,:], -30., 30., 2.,
               corca=ORCA, lkcont=False, cpal='bbr', cfignm=path_fig+cv+'_diff_JJA_'+NAME,
               cbunit='(W/m^2)', cfig_type=fig_type, lat_min=-58., lat_max=65.,  
               ctitle='JJA NEMO - obs : '+cv,
               lforce_lim=True)


    jv = jv + 1


print '\n Bye!'
#
