#!/home/x_laubr/bin/python

# L. Brodeau, november 2009

import sys
import os
import numpy as nmp
from netCDF4 import Dataset

import barakuda_plot as brkdp
import barakuda_tool as brkdt


ldebug = False

ORCA = os.getenv('ORCA')
if ORCA == None: print 'The ORCA environement variable is no set'; sys.exit(0)
RUN = os.getenv('RUN')
if RUN == None: print 'The RUN environement variable is no set'; sys.exit(0)
DIAG_D = os.getenv('DIAG_D')
if DIAG_D == None: print 'The DIAG_D environement variable is no set'; sys.exit(0)

#COMP2D = os.getenv('COMP2D')
#if COMP2D == None: print 'The COMP2D environement variable is no set'; sys.exit(0)


#ICE_CLIM_12 = os.getenv('ICE_CLIM_12')
#if ICE_CLIM_12 == None: print 'The ICE_CLIM_12 environement variable is no set'; sys.exit(0)

#F_T_CLIM_3D_12 = os.getenv('F_T_CLIM_3D_12')
#if F_T_CLIM_3D_12 == None: print 'The F_T_CLIM_3D_12 environement variable is no set\n'
#F_S_CLIM_3D_12 = os.getenv('F_S_CLIM_3D_12')
#if F_S_CLIM_3D_12 == None: print 'The F_S_CLIM_3D_12 environement variable is no set\n'


print ' ORCA = '+ORCA
print ' RUN = '+RUN
print ' DIAG_D = '+DIAG_D
#print ' COMP2D = '+COMP2D
#print ' ICE_CLIM_12 = '+ICE_CLIM_12


CONFRUN = ORCA+'-'+RUN


cv_moc = 'zomsfatl'
path_fig='./'
fig_type='png'

# Mesh-mask file:
cf_mesh_mask = os.getenv('MM_FILE')
if cf_mesh_mask == None: print 'The MM_FILE environement variable (mesh_mask) is no set'; sys.exit(0)
print '\n Mesh-Mask file is:\n', cf_mesh_mask, '\n'

# Basin mask file:
cf_basin_mask = os.getenv('BM_FILE')
if cf_basin_mask == None: print 'The BM_FILE environement variable (basin_mask) is no set'; sys.exit(0)
print '\n Basin mask file is:\n', cf_basin_mask, '\n'



narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <year1> <year2>'; sys.exit(0)
cy1 = sys.argv[1] ; cy2=sys.argv[2]; jy1=int(cy1); jy2=int(cy2)


jy1_clim = jy1 ; jy2_clim = jy2

print ' => mean on the clim : ', jy1_clim, jy2_clim, '\n'



store_dir = os.getenv('DIAG_D')
if store_dir == '': print 'The DIAG_D environement variable is no set'; sys.exit(0)


# Getting coordinates:
brkdt.chck4f(cf_mesh_mask)
id_mm = Dataset(cf_mesh_mask)
xlat = id_mm.variables['gphit'][0,:,:]
xlon = id_mm.variables['glamt'][0,:,:]
Xmask = id_mm.variables['tmask'][0,:,:,:]
#vlev  = id_mm.variables['gdept_0'][0,:]
id_mm.close()


[ nk, nj, ni ] = nmp.shape(Xmask)


# Getting basin mask:
brkdt.chck4f(cf_basin_mask)
id_bm = Dataset(cf_basin_mask)
Xmask_atl = id_bm.variables['tmaskatl'][:,:]
id_bm.close()

for jk in range(nk): Xmask[jk,:,:] = Xmask[jk,:,:] * Xmask_atl[:,:]


#  Getting NEMO mean monthly climatology of MLD coverage:
cf_nemo_moc  = DIAG_D+'/clim/aclim_'+CONFRUN+'_'+cy1+'-'+cy2+'_MOC.nc4'


brkdt.chck4f(cf_nemo_moc)
id_nemo = Dataset(cf_nemo_moc)
vz = id_nemo.variables['depthw'][:]
amoc   = id_nemo.variables[cv_moc][0,:,:]
id_nemo.close()

[ nk, nj ] = amoc.shape ; print ' Shape of AMOC :', nk, nj, '\n'


# Building a latitude vector:
vlat = nmp.zeros(nj)
ji_lat_mid_atlantic = brkdt.find_index_from_value( -28., xlon[0,:] )
vlat[:] = xlat[:,ji_lat_mid_atlantic]


# Building the vertical mask:
msk_vert = nmp.zeros(nk*nj) ; msk_vert.shape = [ nk, nj ]
msk_vert[:,:] = nmp.sum(Xmask[:,:,:],axis=2)
idxm = nmp.where(msk_vert[:,:] > 0.);
msk_vert[idxm] = 1.


brkdp.plot_amoc_lat_depth(vlat[:], -vz[:], amoc[:,:], msk_vert[:,:], -3., 22., 1., \
                         cfig_type=fig_type, lkcont=True, cpal='amoc', ymin=0., ymax=70.,
                         cfignm='AMOC_annual_'+CONFRUN, cbunit='Sv',
                         cxunit=r'Latitude ($^{\circ}$N)', zmin = 5000., zmax = 0., l_zlog=False,
                         czunit='Depth (m)', ctitle='AMOC, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True)



print '\n Bye!'

