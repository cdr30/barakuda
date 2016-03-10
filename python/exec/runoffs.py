#!/usr/bin/env python

# L. Brodeau, 2016

import sys
import numpy as nmp
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp


ldebug = False

rmult = 1.E3
zmax_rnf_atl = 0.1 ; dz_rnf = 0.005

venv_needed = {'ORCA','RUN','DIAG_D','NN_RNF','MM_FILE'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)


if 'ORCA2' in vdic['ORCA']:
    ji_lat0 = 132
elif 'ORCA1' in vdic['ORCA']:
    #ji_lat0 = 265
    ji_lat0 = 100
else:
    print 'FIX ME!!! '+sys.argv[0]+' => dont know ji_lat0 for conf '+vdic['ORCA']+' !!!'; sys.exit(0)




CONFRUN = vdic['ORCA']+'-'+vdic['RUN']


path_fig='./'
fig_type='png'



narg = len(sys.argv)
if narg < 3: 
    print 'Usage: '+sys.argv[0]+' <year1> <year2>'; sys.exit(0)
cy1 = sys.argv[1] ; cy2=sys.argv[2]; jy1=int(cy1); jy2=int(cy2)


jy1_clim = jy1 ; jy2_clim = jy2

print ' => mean on the clim : ', jy1_clim, jy2_clim, '\n'




# Getting coordinates:
bt.chck4f(vdic['MM_FILE'])
id_mm = Dataset(vdic['MM_FILE'])
xlon   = id_mm.variables['glamt'][0,:,:] ; xlat = id_mm.variables['gphit'][0,:,:]
Xmask = id_mm.variables['tmask'][0,:,:,:]
vlev  = id_mm.variables['gdept_1d'][0,:]
id_mm.close()

nk = len(vlev)



#  Getting NEMO mean monthly climatology of RNF coverage:
cf_nemo_mnmc = vdic['DIAG_D']+'/clim/mclim_'+CONFRUN+'_'+cy1+'-'+cy2+'_SBC.nc4'

bt.chck4f(cf_nemo_mnmc)
id_nemo = Dataset(cf_nemo_mnmc)
rnf   = rmult*id_nemo.variables[vdic['NN_RNF']][:,:,:]
id_nemo.close()

[ nt, nj, ni ] = rnf.shape ; print ' Shape of Runoffs :', nt, nj, ni, '\n'





rnf_plot = nmp.zeros(nj*ni) ; rnf_plot.shape = [ nj , ni ]


rnf_plot[:,:] = nmp.mean(rnf[:,:,:],axis=0)
#rnf_plot[:,:] = nmp.log(rnf_plot[:,:]+1.0)

# With lat-lon axis:
#bp.plot_2d(xlon[0,:], xlat[:,ji_lat0], rnf_plot[:,:], Xmask[0,:,:], 0., zmax_rnf_atl, dz_rnf,
#              corca=vdic['ORCA'], lkcont=False, cpal='sst',
#              cfignm=path_fig+'runoffs_mean_'+CONFRUN, cbunit=r'10$^{-3}$mm/day',
#              ctitle='Mean Runoffs, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True, i_cb_subsamp=2,
#              cfig_type=fig_type, lat_min=-79., lat_max=85., lpix=True)

# Without:
bp.plot("2d")([0], [0], rnf_plot[:,:], Xmask[0,:,:], 0., zmax_rnf_atl, dz_rnf,
           corca=vdic['ORCA'], lkcont=False, cpal='sst',
           cfignm=path_fig+'runoffs_mean_'+CONFRUN, cbunit=r'10$^{-3}$mm/day',
           ctitle='Mean Runoffs, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True, i_cb_subsamp=2,
           cfig_type=fig_type, lpix=True)


print '\n Bye!'

