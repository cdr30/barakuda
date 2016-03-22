#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate lat-depth plot of the Atlantic Meridional Overturning Circulation
#
#       L. Brodeau, 2009

import sys
import numpy as nmp
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp

ldebug = False

venv_needed = {'ORCA','RUN','DIAG_D','MM_FILE','BM_FILE'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

cv_moc = 'zomsfatl'
path_fig='./'
fig_type='png'


narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <year1> <year2>'; sys.exit(0)
cy1 = sys.argv[1] ; cy2=sys.argv[2]; jy1=int(cy1); jy2=int(cy2)


jy1_clim = jy1 ; jy2_clim = jy2

print ' => mean on the clim : ', jy1_clim, jy2_clim, '\n'


# Getting coordinates:
bt.chck4f(vdic['MM_FILE'])
id_mm = Dataset(vdic['MM_FILE'])
xlat = id_mm.variables['gphit'][0,:,:]
xlon = id_mm.variables['glamt'][0,:,:]
Xmask = id_mm.variables['tmask'][0,:,:,:]
id_mm.close()


[ nk, nj, ni ] = nmp.shape(Xmask)


# Getting basin mask:
bt.chck4f(vdic['BM_FILE'])
id_bm = Dataset(vdic['BM_FILE'])
Xmask_atl = id_bm.variables['tmaskatl'][:,:]
id_bm.close()

for jk in range(nk): Xmask[jk,:,:] = Xmask[jk,:,:] * Xmask_atl[:,:]


#  Getting NEMO mean monthly climatology of MLD coverage:
cf_nemo_moc  = vdic['DIAG_D']+'/clim/aclim_'+CONFRUN+'_'+cy1+'-'+cy2+'_MOC.nc4'


bt.chck4f(cf_nemo_moc)
id_nemo = Dataset(cf_nemo_moc)
vz = id_nemo.variables['depthw'][:]
amoc   = id_nemo.variables[cv_moc][0,:,:]
id_nemo.close()

[ nk, nj ] = amoc.shape ; print ' Shape of AMOC :', nk, nj, '\n'


# Building a latitude vector:
vlat = nmp.zeros(nj)
ji_lat_mid_atlantic = bt.find_index_from_value( -28., xlon[0,:] )
vlat[:] = xlat[:,ji_lat_mid_atlantic]


# Building the vertical mask:
msk_vert = nmp.zeros((nk,nj))
msk_vert[:,:] = nmp.sum(Xmask[:,:,:],axis=2)
idxm = nmp.where(msk_vert[:,:] > 0.);
msk_vert[idxm] = 1.


bp.plot("amoc_lat_depth")(vlat[:], -vz[:], amoc[:,:], msk_vert[:,:], -3., 20., 1., \
                          cfig_type=fig_type, lkcont=True, cpal='amoc', ymin=0., ymax=70.,
                          cfignm='AMOC_annual_'+CONFRUN, cbunit='Sv',
                          cxunit=r'Latitude ($^{\circ}$N)', zmin = 5000., zmax = 0., l_zlog=False,
                          czunit='Depth (m)', ctitle='AMOC, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                          lforce_lim=True, i_cb_subsamp=1)



print '\n Bye!'

