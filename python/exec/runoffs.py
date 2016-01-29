#!/home/x_laubr/bin/python

# L. Brodeau, 2016

import sys
import os
import numpy as nmp
from netCDF4 import Dataset

import barakuda_orca    as brkdo
import barakuda_plot    as brkdp
import barakuda_physics as bphys
import barakuda_tool    as brkdt

ldebug = False

zmax_rnf_atl = 1.E-4 ; dz_rnf = 1.E-3

ORCA = os.getenv('ORCA')
if ORCA == None: print 'The ORCA environement variable is no set'; sys.exit(0)
RUN = os.getenv('RUN')
if RUN == None: print 'The RUN environement variable is no set'; sys.exit(0)
DIAG_D = os.getenv('DIAG_D')
if DIAG_D == None: print 'The DIAG_D environement variable is no set'; sys.exit(0)

NN_RNF = os.getenv('NN_RNF')
if NN_RNF == None: print 'The NN_RNF environement variable is no set'; sys.exit(0)


print ' ORCA = '+ORCA
print ' RUN = '+RUN
print ' DIAG_D = '+DIAG_D



if 'ORCA2' in ORCA:
    ji_lat0 = 132
elif 'ORCA1' in ORCA:
    ji_lat0 = 265
else:
    print 'FIX ME!!! '+sys.argv[0]+' => dont know ji_lat0 for conf '+ORCA+' !!!'; sys.exit(0)




CONFRUN = ORCA+'-'+RUN


path_fig='./'
fig_type='png'

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
vlev  = id_mm.variables['gdept_1d'][0,:]
id_mm.close()

nk = len(vlev)



#  Getting NEMO mean monthly climatology of RNF coverage:
cf_nemo_mnmc = DIAG_D+'/clim/mclim_'+CONFRUN+'_'+cy1+'-'+cy2+'_SBC.nc4'

brkdt.chck4f(cf_nemo_mnmc)
id_nemo = Dataset(cf_nemo_mnmc)
rnf   = 1.E3*id_nemo.variables[NN_RNF][:,:,:]
id_nemo.close()

[ nt, nj, ni ] = rnf.shape ; print ' Shape of RNF :', nt, nj, ni, '\n'





rnf_plot = nmp.zeros(nj*ni) ; rnf_plot.shape = [ nj , ni ]


rnf_plot[:,:] = nmp.mean(rnf[:,:,:],axis=0)


#lolo zmax_rnf_atl, dz_rnf,
brkdp.plot_2d(xlon[0,:], xlat[:,ji_lat0], rnf_plot[:,:], Xmask[0,:,:], 0., 0.1, 0.005,
              corca=ORCA, lkcont=False, cpal='jet',
              cfignm=path_fig+'rnf_mean_'+CONFRUN, cbunit='mm/day',
              ctitle='Mean RNF, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True, i_sub_samp=2,
              cfig_type=fig_type, lat_min=-82., lat_max=80., lpix=True)

#, vcont_spec = [ 0. ])

sys.exit(0)





# FIGURES MARCH #
#################


brkdp.plot_nproj('nseas', 200., zmax_mld_atl, dz_mld, xlon, xlat, mldr10[imnth,:,:],
                cfignm=path_fig+'mld_NSeas_march_'+CONFRUN, cpal='sst0', cbunit='m',
                ctitle='MLD, March, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                lkcont=True, cfig_type=fig_type,
                lforce_lim=True)

brkdp.plot_nproj('spstere', 50., 200., 10., xlon, xlat, mldr10[imnth,:,:],
                cfignm=path_fig+'mld_ACC_march_'+CONFRUN, cpal='sst0', cbunit='m',
                ctitle='MLD, March, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                lkcont=True, cfig_type=fig_type,
                lforce_lim=True)

brkdp.plot_2d(xlon[0,:], xlat[:,ji_lat0], mldr10[imnth,:,:], Xmask[0,:,:], 0., 600., 20.,
             corca=ORCA, lkcont=True, cpal='sst0',
             cfignm=path_fig+'mld_Global_march_'+CONFRUN, cbunit='m',
             ctitle='MLD, March, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True, i_sub_samp=1,
             cfig_type=fig_type, lat_min=-80., lat_max=75., lpix=False)

if l_obs_mld:
    brkdp.plot_nproj('nseas', 200., zmax_mld_atl, dz_mld, xlon, xlat, Xmld_obs[:,:],
                    cfignm=path_fig+'mld_obs_001_NSeas_march_'+CONFRUN, cpal='sst0', cbunit='m',
                    ctitle='MLD (obs., 0.01 crit.), March (Levitus 1980-1999)',
                    lkcont=True, cfig_type=fig_type, lforce_lim=True)

    brkdp.plot_2d(xlon[0,:], xlat[:,ji_lat0], Xmld_obs[:,:], Xmask[0,:,:], 0., 600., 20.,
                 corca=ORCA, lkcont=True, cpal='sst0',
                 cfignm=path_fig+'mld_obs_001_Global_march_'+CONFRUN, cbunit='m',
                 ctitle='MLD (obs., 0.01 crit.), March (Levitus 1980-1999)', lforce_lim=True, i_sub_samp=1,
                 cfig_type=fig_type, lat_min=-80., lat_max=75., lpix=False)








#####################
# S E P T E M B E R #
#####################

imnth = 8 ; # september


if l_obs_mld:
    # Computing sigma0 3D field:
    Sigma0 = bphys.sigma0(Tclim[imnth,:,:,:], Sclim[imnth,:,:,:])*Xmask[:,:,:]

    mmask[:,:] = 1.
    for jk in range(nk-1):
        zz = vlev[jk]
        Sigma0[jk,:,:]   = Sigma0[jk,:,:]*mmask[:,:]
        Sigma0[jk+1,:,:] = Sigma0[jk+1,:,:]*mmask[:,:]
        ijloc = nmp.where( Sigma0[jk+1,:,:] - 0.01 > Sigma0[jk,:,:] )
        #lolo:Xmld_obs[ijloc] = zz ;# Problem lolo!!!
        mmask[ijloc] = 0. ; # these points won't be checked again only first occurence of the criterion matters!



# Figures september:
brkdp.plot_nproj('spstere', 100., 2000., dz_mld, xlon, xlat, mldr10[imnth,:,:],
                cfignm=path_fig+'mld_ACC_september_'+CONFRUN, cpal='sst0', cbunit='m',
                ctitle='MLD, September, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                lkcont=True, cfig_type=fig_type,
                lforce_lim=True)

brkdp.plot_2d(xlon[0,:], xlat[:,ji_lat0], mldr10[imnth,:,:], Xmask[0,:,:], 0., 600., 20.,
             corca=ORCA, lkcont=True, cpal='sst0',
             cfignm=path_fig+'mld_Global_september_'+CONFRUN, cbunit='m',
             ctitle='MLD, September, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True, i_sub_samp=1,
             cfig_type=fig_type, lat_min=-80., lat_max=75., lpix=False)


if l_obs_mld:
    brkdp.plot_2d(xlon[0,:], xlat[:,ji_lat0], Xmld_obs[:,:], Xmask[0,:,:], 0., 600., 20.,
                 corca=ORCA, lkcont=True, cpal='sst0',
                 cfignm=path_fig+'mld_obs_001_Global_september_'+CONFRUN, cbunit='m',
                 ctitle='MLD (obs, 0.01 crit.), March (Levitus 1980-1999)', lforce_lim=True, i_sub_samp=1,
                 cfig_type=fig_type, lat_min=-80., lat_max=75., lpix=False)










print '\n Bye!'

