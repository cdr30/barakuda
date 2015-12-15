# L. Brodeau, november 2013

import sys
import os
import numpy as nmp
from netCDF4 import Dataset

import barakuda_orca as brkdo
import barakuda_plot as brkdp
import barakuda_tool as brkdt




fig_type='png'

ORCA = os.getenv('ORCA')
if ORCA == None: print 'The ORCA environement variable is no set'; sys.exit(0)

MM_FILE = os.getenv('MM_FILE')
if MM_FILE == None: print 'The MM_FILE environement variable is no set'; sys.exit(0)

BM_FILE = os.getenv('BM_FILE')
if BM_FILE == None: print 'The BM_FILE environement variable is no set'; sys.exit(0)

DIAG_DIR = os.getenv('DIAG_DIR')
if DIAG_DIR == None: print 'The DIAG_DIR environement variable is no set'; sys.exit(0)

print '\n '+sys.argv[0]+':'


narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <RUN1> <RUN2>'; sys.exit(0)
crun1 = sys.argv[1] ; crun2=sys.argv[2]

CONFRUN1 = ORCA+'-'+crun1 ; print CONFRUN1
CONFRUN2 = ORCA+'-'+crun2 ; print CONFRUN2


# Getting coordinates and mask
cf_mask = MM_FILE ; print 'Reading coordinates into '+cf_mask
brkdt.chck4f(cf_mask)
id_mask = Dataset(cf_mask)
mask  = id_mask.variables['tmask'][0,:,:,:]
xlon   = id_mask.variables['nav_lon'][:,:] ; xlat   = id_mask.variables['nav_lat'][:,:]
id_mask.close()


[ nk, nj, ni ] = nmp.shape(mask)


# Getting basin mask:
cf_basin_mask = BM_FILE
brkdt.chck4f(cf_basin_mask)
id_bm = Dataset(cf_basin_mask)
Xmask_atl = id_bm.variables['tmaskatl'][:,:,:]
id_bm.close()


ji_lat_mid_atlantic = brkdt.find_index_from_value( -28., xlon[0,:] )

# Creating 1D long. and lat.:
vlon = nmp.zeros(ni) ;  vlon.shape = [ ni ] ; vlon[:] = xlon[0,:]
vlat = nmp.zeros(nj) ;  vlat.shape = [ nj ] ; vlat[:] = xlat[:,ji_lat_mid_atlantic]


# Getting first and last year of NEMO clim from "last_clim" files:
ccf1 = DIAG_DIR+'/'+CONFRUN1+'/last_clim' ; brkdt.chck4f(ccf1)
ccf2 = DIAG_DIR+'/'+CONFRUN2+'/last_clim' ; brkdt.chck4f(ccf2)
f = open(ccf1, 'r'); cread_lines = f.readlines() ; f.close()
cy1_1, cy2_1 = cread_lines[0].split("-") ; cy2_1, ca = cy2_1.split("\n")
f = open(ccf2, 'r'); cread_lines = f.readlines() ; f.close()
cy1_2, cy2_2 = cread_lines[0].split("-") ; cy2_2, ca = cy2_2.split("\n")
#print cy1_1, cy2_1, cy1_2, cy2_2






# Getting AMOC clim:
# -------------------

cf_nemo_moc1 = DIAG_DIR+'/'+CONFRUN1+'/clim/aclim_'+CONFRUN1+'_'+cy1_1+'-'+cy2_1+'_MOC.nc4' ; brkdt.chck4f(cf_nemo_moc1)
cf_nemo_moc2 = DIAG_DIR+'/'+CONFRUN2+'/clim/aclim_'+CONFRUN2+'_'+cy1_2+'-'+cy2_2+'_MOC.nc4' ; brkdt.chck4f(cf_nemo_moc2)


id_nemo_moc1 = Dataset(cf_nemo_moc1)
Amoc1   = id_nemo_moc1.variables[brkdo.cv_amoc][:,:,:]
vz = id_nemo_moc1.variables['depthw'][:]
id_nemo_moc1.close()

[ nt, nk0, nj0 ] = Amoc1.shape
if [ nk0, nj0 ] != [ nk, nj ]: print 'ERROR: AMOC clim and mesh_mask do no agree in shape! #1'; sys.exit(0)

id_nemo_moc2 = Dataset(cf_nemo_moc2)
Amoc2   = id_nemo_moc2.variables[brkdo.cv_amoc][:,:,:]
id_nemo_moc2.close()

[ nt0, nk0, nj0 ] = Amoc2.shape
if [ nt0, nk0, nj0 ] != [ nt, nk, nj ]: print 'ERROR: AMOC clim and mesh_mask do no agree in shape! #2'; sys.exit(0)

if nt != 1: print 'ERROR: AMOC clim supposed to have only 1 time record'; sys.exit(0)

# Building the vertical mask:
msk_vert = nmp.zeros(nk*nj) ; msk_vert.shape = [ nk, nj ]
msk_vert[:,:] = nmp.sum(Xmask_atl[:,:,:],axis=2)
idxm = nmp.where(msk_vert[:,:] > 0.);
msk_vert[idxm] = 1.


brkdp.plot_amoc_lat_depth(vlat[:], -vz[:], Amoc1[0,:,:], msk_vert[:,:], -3.5, 22., 0.5, \
                          cfig_type=fig_type, lkcont=True, cpal='amoc', ymin=0., ymax=70.,
                          cfignm='AMOC_'+ORCA+'-'+crun1, cbunit='Sv',
                          cxunit=r'Latitude ($^{\circ}$N)', zmin = 5000., zmax = 0., l_zlog=False,
                          czunit='Depth (m)', ctitle='AMOC, ('+cy1_1+'-'+cy2_1+'), '+ORCA+'-'+crun1, lforce_lim=True)

brkdp.plot_amoc_lat_depth(vlat[:], -vz[:], Amoc2[0,:,:], msk_vert[:,:], -3.5, 22., 0.5, \
                          cfig_type=fig_type, lkcont=True, cpal='amoc', ymin=0., ymax=70.,
                          cfignm='AMOC_'+ORCA+'-'+crun2, cbunit='Sv',
                          cxunit=r'Latitude ($^{\circ}$N)', zmin = 5000., zmax = 0., l_zlog=False,
                          czunit='Depth (m)', ctitle='AMOC, ('+cy1_2+'-'+cy2_2+'), '+ORCA+'-'+crun1, lforce_lim=True)


# Difference of AMOC:
brkdp.plot_amoc_lat_depth(vlat[:], -vz[:], Amoc2[0,:,:]-Amoc1[0,:,:], msk_vert[:,:], -0.36, 0.36, 0.02, \
                          cfig_type=fig_type, lkcont=True, cpal='bbr2', ymin=0., ymax=70.,
                          cfignm='AMOC_'+ORCA+'_'+crun2+'-'+crun1, cbunit='Sv',
                          cxunit=r'Latitude ($^{\circ}$N)', zmin = 5000., zmax = 0., l_zlog=False,
                          czunit='Depth (m)', ctitle='AMOC, ('+cy1_1+'-'+cy2_1+'), '+ORCA+' '+crun2+'-'+crun1, lforce_lim=True)

print '  Diff. of AMOC done!\n'






# Getting NEMO mean monthly climatology on grid-T file:
# -----------------------------------------------------

cf_nemo_mnmc1 = DIAG_DIR+'/'+CONFRUN1+'/clim/mclim_'+CONFRUN1+'_'+cy1_1+'-'+cy2_1+'_grid_T.nc4'; brkdt.chck4f(cf_nemo_mnmc1)
cf_nemo_mnmc2 = DIAG_DIR+'/'+CONFRUN2+'/clim/mclim_'+CONFRUN2+'_'+cy1_2+'-'+cy2_2+'_grid_T.nc4'; brkdt.chck4f(cf_nemo_mnmc2)

id_nemo_mnmc1 = Dataset(cf_nemo_mnmc1)
Tnemo1  = id_nemo_mnmc1.variables[brkdo.cv_temp][:,:,:,:]
Snemo1  = id_nemo_mnmc1.variables[brkdo.cv_sali][:,:,:,:]
vdepth  = id_nemo_mnmc1.variables['deptht'][:]
id_nemo_mnmc1.close()
[ nt, nk0, nj0, ni0 ] = Tnemo1.shape
if nk != nk0 or nj != nj0 or ni != ni0:
    print 'ERROR: 3D clim and NEMO file do no agree in shape! #1'; sys.exit(0)

id_nemo_mnmc2 = Dataset(cf_nemo_mnmc2)
Tnemo2  = id_nemo_mnmc2.variables[brkdo.cv_temp][:,:,:,:]
Snemo2  = id_nemo_mnmc2.variables[brkdo.cv_sali][:,:,:,:]
id_nemo_mnmc2.close()
[ nt0, nk0, nj0, ni0 ] = Tnemo1.shape
if nk != nk0 or nj != nj0 or ni != ni0 or nt != nt0:
    print 'ERROR: 3D clim and NEMO file do no agree in shape! #2'; sys.exit(0)



ji30W = brkdt.find_index_from_value( -30., xlon[0,:] ) ; lon30W = xlon[0,ji30W]; print ' ji30W =', ji30W, lon30W

lat_max_north = 40.



brkdp.plot_vert_section(xlat[:,ji30W], vdepth, nmp.mean(Tnemo1[:,:,:,ji30W],axis=0),
                        mask[:,:,ji30W], -1., 15., 0.5, cpal='mld', xmin=-70., xmax=lat_max_north,
                        cfignm='section_temp_'+ORCA+'-'+crun1, cbunit=r'$^{\circ}$C', cxunit=r'Latitude ($^{\circ}$N)',
                        czunit='Depth (m)', ctitle='Temperature, ('+cy1_1+'-'+cy2_1+'), '+ORCA+'-'+crun1+', lon = 30W',
                        cfig_type=fig_type, lforce_lim=True, i_sub_samp=2)

brkdp.plot_vert_section(xlat[:,ji30W], vdepth, nmp.mean(Tnemo2[:,:,:,ji30W],axis=0),
                        mask[:,:,ji30W], -1., 15., 0.5, cpal='mld', xmin=-70., xmax=lat_max_north,
                        cfignm='section_temp_'+ORCA+'-'+crun2, cbunit=r'$^{\circ}$C', cxunit=r'Latitude ($^{\circ}$N)',
                        czunit='Depth (m)', ctitle='Temperature, ('+cy1_1+'-'+cy2_1+'), '+ORCA+'-'+crun2+', lon = 30W',
                        cfig_type=fig_type, lforce_lim=True, i_sub_samp=2)

brkdp.plot_vert_section(xlat[:,ji30W], vdepth, nmp.mean(Tnemo2[:,:,:,ji30W],axis=0) - nmp.mean(Tnemo1[:,:,:,ji30W],axis=0),
                        mask[:,:,ji30W], -0.1, 0.1, 0.005, cpal='bbr2', xmin=-70., xmax=lat_max_north,
                        cfignm='section_temp_'+ORCA+'_'+crun2+'-'+crun1, cbunit=r'$^{\circ}$C', cxunit=r'Latitude ($^{\circ}$N)',
                        czunit='Depth (m)', ctitle='Temperature, ('+cy1_1+'-'+cy2_1+'), '+ORCA+' '+crun2+'-'+crun1+', lon = 30W',
                        cfig_type=fig_type, lforce_lim=True, i_sub_samp=2)







brkdp.plot_vert_section(xlat[:,ji30W], vdepth, nmp.mean(Snemo2[:,:,:,ji30W],axis=0) - nmp.mean(Snemo1[:,:,:,ji30W],axis=0),
                        mask[:,:,ji30W], -0.02, 0.02, 0.002, cpal='bbr2', xmin=-70., xmax=lat_max_north,
                        cfignm='section_sali_'+ORCA+'_'+crun2+'-'+crun1, cbunit=r'PSU', cxunit=r'Latitude ($^{\circ}$N)',
                        czunit='Depth (m)', ctitle='Salinity, ('+cy1_1+'-'+cy2_1+'), '+ORCA+' '+crun2+'-'+crun1+', lon = 30W',
                        cfig_type=fig_type, lforce_lim=True, i_sub_samp=1)

