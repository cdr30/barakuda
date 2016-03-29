#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate global plot of sea surface height
#
#       L. Brodeau, 2009

import sys
import numpy as nmp
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp
import barakuda_physics as bphys

venv_needed = {'ORCA','RUN','DIAG_D','MM_FILE','NN_SSH'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

path_fig='./'
fig_type='png'


# What ji point to use to extract latitude vector:
if 'ORCA2' in vdic['ORCA']:
    ji_lat0 = 132
elif 'ORCA1' in vdic['ORCA']:
    ji_lat0 = 265
else:
    print 'FIX ME!!! ssh.py => dont know ji_lat0 for conf '+vdic['ORCA']+' !!!'; sys.exit(0)


narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <year1> <year2>'; sys.exit(0)
cy1 = sys.argv[1] ; cy2=sys.argv[2]; jy1=int(cy1); jy2=int(cy2)


jy1_clim = jy1 ; jy2_clim = jy2
print ' => mean on the clim : ', jy1_clim, jy2_clim, '\n'


# Getting coordinates:
bt.chck4f(vdic['MM_FILE'])
id_mm = Dataset(vdic['MM_FILE'])
xlon   = id_mm.variables['glamt'][0,:,:] ; xlat = id_mm.variables['gphit'][0,:,:]
Xmask = id_mm.variables['tmask'][0,0,:,:]
Xe1t = id_mm.variables['e1t'][0,:,:]
Xe2t = id_mm.variables['e2t'][0,:,:]
id_mm.close()


#  Getting NEMO mean monthly climatology of SSH coverage:
cf_nemo_mnmc = vdic['DIAG_D']+'/clim/mclim_'+CONFRUN+'_'+cy1+'-'+cy2+'_grid_T.nc4'

bt.chck4f(cf_nemo_mnmc)
id_nemo = Dataset(cf_nemo_mnmc)
ssh   = id_nemo.variables[vdic['NN_SSH']][:,:,:]
id_nemo.close()

[ nt, nj, ni ] = ssh.shape ; print ' Shape of SSH :', nt, nj, ni, '\n'

ssh_plot = nmp.zeros((nj,ni))

ssh_plot[:,:] = nmp.mean(ssh[:,:,:],axis=0)


ztot = nmp.sum(ssh_plot*Xmask*Xe1t*Xe2t)/nmp.sum(Xmask*Xe1t*Xe2t)
print 'ztot =', ztot

ssh_plot = ssh_plot - ztot

bp.plot("2d")(xlon[0,:], xlat[:,ji_lat0], ssh_plot[:,:], Xmask, -2., 2., 0.1,
              corca=vdic['ORCA'], lkcont=True, cpal='BrBG_r',
              cfignm=path_fig+'ssh_mean_'+CONFRUN, cbunit=r'$(m)$',
              ctitle='Mean SSH (corrected about z=0), '+CONFRUN+' ('+cy1+'-'+cy2+')',
              lforce_lim=True, i_cb_subsamp=2,
              cfig_type=fig_type, lat_min=-77., lat_max=75., lpix=False, vcont_spec = [ 0. ])

