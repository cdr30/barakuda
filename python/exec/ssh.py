#!/home/x_laubr/bin/python

# L. Brodeau, november 2009

import sys
import os
import numpy as nmp
from netCDF4 import Dataset

import barakuda_orca    as brkdo
import barakuda_plot    as brkdp
import barakuda_physics as bphys
import barakuda_tool    as brkdt


ldebug = False

ORCA = os.getenv('ORCA')
if ORCA == None: print 'The ORCA environement variable is no set'; sys.exit(0)
RUN = os.getenv('RUN')
if RUN == None: print 'The RUN environement variable is no set'; sys.exit(0)
DIAG_D = os.getenv('DIAG_D')
if DIAG_D == None: print 'The DIAG_D environement variable is no set'; sys.exit(0)
NN_SSH = os.getenv('NN_SSH')
if NN_SSH == None: print 'The NN_SSH environement variable is no set'; sys.exit(0)

#COMP2D = os.getenv('COMP2D')
#if COMP2D == None: print 'The COMP2D environement variable is no set'; sys.exit(0)


print ' ORCA = '+ORCA
print ' RUN = '+RUN
print ' DIAG_D = '+DIAG_D
#print ' COMP2D = '+COMP2D







CONFRUN = ORCA+'-'+RUN


path_fig='./'
fig_type='png'


# What ji point to use to extract latitude vector:

if 'ORCA2' in ORCA:
    ji_lat0 = 132
elif 'ORCA1' in ORCA:
    ji_lat0 = 265
else:
    print 'FIX ME!!! ssh.py => dont know ji_lat0 for conf '+ORCA+' !!!'; sys.exit(0)



# Mesh-mask file:
cf_mesh_mask = os.getenv('MM_FILE')
if cf_mesh_mask == None: print 'The MM_FILE environement variable (mesh_mask) is no set'; sys.exit(0)
print '\n Mesh-Mask file is:\n', cf_mesh_mask, '\n'



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
xlon   = id_mm.variables['glamt'][0,:,:] ; xlat = id_mm.variables['gphit'][0,:,:]
Xmask = id_mm.variables['tmask'][0,:,:,:]
id_mm.close()


#  Getting NEMO mean monthly climatology of SSH coverage:
cf_nemo_mnmc = DIAG_D+'/clim/mclim_'+CONFRUN+'_'+cy1+'-'+cy2+'_grid_T.nc4'

brkdt.chck4f(cf_nemo_mnmc)
id_nemo = Dataset(cf_nemo_mnmc)
ssh   = id_nemo.variables[NN_SSH][:,:,:]
id_nemo.close()

[ nt, nj, ni ] = ssh.shape ; print ' Shape of SSH :', nt, nj, ni, '\n'



ssh_plot = nmp.zeros(nj*ni) ; ssh_plot.shape = [ nj , ni ]


ssh_plot[:,:] = nmp.mean(ssh[:,:,:],axis=0)



brkdp.plot_2d(xlon[0,:], xlat[:,ji_lat0], ssh_plot[:,:], Xmask[0,:,:], -3.5, 0.6, 0.1,
              corca=ORCA, lkcont=True, cpal='jet',
              cfignm=path_fig+'ssh_mean_'+CONFRUN, cbunit='m',
              ctitle='Mean SSH, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True, i_sub_samp=2,
              cfig_type=fig_type, lat_min=-77., lat_max=75., lpix=False, vcont_spec = [ 0. ])

