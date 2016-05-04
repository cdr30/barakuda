#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate 2D plots and maps of the Mixed layer depth
#
#       L. Brodeau, 2009

import sys
import numpy as nmp
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp
import barakuda_physics as bphys

venv_needed = {'ORCA','RUN','DIAG_D','COMP2D','MM_FILE','NN_MLD','NN_S_CLIM','NN_T_CLIM','F_T_CLIM_3D_12','F_S_CLIM_3D_12'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

corca = vdic['ORCA']

CONFRUN = corca+'-'+vdic['RUN']
ldebug = False

zmax_mld_atl = 1600. ; dz_mld = 100.




if 'ORCA2' in corca:
    ji_lat0 = 132
elif 'eORCA1' in vdic['ORCA']:
    ji_lat0 = 265
elif 'ORCA1' in vdic['ORCA']:
    ji_lat0 = 265
else:
    print 'FIX ME!!! ssh.py => dont know ji_lat0 for conf '+corca+' !!!'; sys.exit(0)


l_obs_mld = False
if (not vdic['F_T_CLIM_3D_12'] == None) and (not vdic['F_T_CLIM_3D_12'] == None):
    l_obs_mld = True ; print 'Since obs Temp and Sali here will compute observed MLD!'


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
xlon   = id_mm.variables['glamt'][0,:,:] ; xlat = id_mm.variables['gphit'][0,:,:]
Xmask = id_mm.variables['tmask'][0,:,:,:]
vlev  = id_mm.variables['gdept_1d'][0,:]
id_mm.close()

nk = len(vlev)



#  Getting NEMO mean monthly climatology of MLD coverage:
cf_nemo_mnmc = vdic['DIAG_D']+'/clim/mclim_'+CONFRUN+'_'+cy1+'-'+cy2+'_grid_T.nc4'

bt.chck4f(cf_nemo_mnmc)
id_nemo = Dataset(cf_nemo_mnmc)
mldr10   = id_nemo.variables[vdic['NN_MLD']][:,:,:]
id_nemo.close()

[ nt, nj, ni ] = mldr10.shape ; print ' Shape of MLD :', nt, nj, ni, '\n'




if l_obs_mld:

    # Getting 3D+T 12-month climatology of T and S
    # --------------------------------------------

    # Temperature
    bt.chck4f(vdic['F_T_CLIM_3D_12']) ; id_clim = Dataset(vdic['F_T_CLIM_3D_12'])
    Tclim  = id_clim.variables[vdic['NN_T_CLIM']][:,:,:,:]; print '(has ',Tclim.shape[0],' time snapshots)\n'
    id_clim.close()

    # Salinity
    bt.chck4f(vdic['F_S_CLIM_3D_12']) ; id_clim = Dataset(vdic['F_S_CLIM_3D_12'])
    Sclim  = id_clim.variables[vdic['NN_S_CLIM']][:,:,:,:]; print '(has ',Sclim.shape[0],' time snapshots)\n'
    id_clim.close()

    [ nmn , nk0, nj0 , ni0 ] = Tclim.shape

    if nj != nj0 or ni != ni0 or nk != nk0:
        print 'ERROR: 3D clim and NEMO file do no agree in shape!'
        print '       clim => '+str(ni0)+', '+str(nj0)+', '+str(nk0)+' ('+vdic['F_T_CLIM_3D_12']+')'
        print '       NEMO => '+str(ni)+', '+str(nj)+', '+str(nk)
        sys.exit(0)





    #############
    # M A R C H #
    #############

    imnth = 2 ; # march

    # Computing sigma0 3D field:
    Sigma0 = bphys.sigma0(Tclim[imnth,:,:,:], Sclim[imnth,:,:,:])*Xmask[:,:,:]

    if ldebug: Sigma0_nemo = bphys.sigma0(Tnemo[imnth:,:,:], Snemo[imnth,:,:,:])*Xmask[:,:,:]

    Xmld_obs = nmp.zeros((nj,ni))
    mmask = nmp.zeros((nj,ni))

    mmask[:,:] = 1.
    for jk in range(nk-1):
        zz = vlev[jk]
        Sigma0[jk,:,:]   = Sigma0[jk,:,:]*mmask[:,:]
        Sigma0[jk+1,:,:] = Sigma0[jk+1,:,:]*mmask[:,:]
        ijloc = nmp.where( Sigma0[jk+1,:,:] - 0.01 > Sigma0[jk,:,:] )
        mmask[ijloc] = 0. ; # these points won't be checked again only first occurence of the criterion matters!


    if ldebug:
        # Testing my MLD method built of T and S from NEMO (to check the obs. MLD I build the same way...)
        Xmld_obs[:,:] = 0.
        mmask[:,:] = 1.
        for jk in range(nk-1):
            zz = vlev[jk]
            Sigma0_nemo[jk,:,:]   = Sigma0_nemo[jk,:,:]*mmask[:,:]
            Sigma0_nemo[jk+1,:,:] = Sigma0_nemo[jk+1,:,:]*mmask[:,:]
            ijloc = nmp.where( Sigma0_nemo[jk+1,:,:] - 0.01 > Sigma0_nemo[jk,:,:] )
            Xmld_obs[ijloc] = zz
            mmask[ijloc] = 0. ; # these points won't be checked again only first occurence of the criterion matters!
        bp.plot("nproj")('nseas', 200., zmax_mld_atl, dz_mld, xlon, xlat, Xmld_obs[:,:],
                         cfignm=path_fig+'mld_NEMO_001_NSeas_march_'+CONFRUN+'_vs_'+vdic['COMP2D'], cpal='sst0', cbunit='m',
                         ctitle='MLD NEMO (0.01 crit.), March, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                         lkcont=True, cfig_type=fig_type, lforce_lim=True)



# FIGURES MARCH #
#################


bp.plot("nproj")('nseas', 200., zmax_mld_atl, dz_mld, xlon, xlat, mldr10[imnth,:,:],
                 cfignm=path_fig+'mld_NSeas_march_'+CONFRUN, cpal='sst0', cbunit='m',
                 ctitle='MLD, March, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                 lkcont=True, cfig_type=fig_type,
                 lforce_lim=True)

bp.plot("nproj")('spstere', 50., 200., 10., xlon, xlat, mldr10[imnth,:,:],
                 cfignm=path_fig+'mld_ACC_march_'+CONFRUN, cpal='sst0', cbunit='m',
                 ctitle='MLD, March, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                 lkcont=True, cfig_type=fig_type,
                 lforce_lim=True)

bp.plot("2d")(xlon[0,:], xlat[:,ji_lat0], mldr10[imnth,:,:], Xmask[0,:,:], 0., 600., 20.,
              corca=vdic['ORCA'], lkcont=True, cpal='sst0',
              cfignm=path_fig+'mld_Global_march_'+CONFRUN, cbunit='m',
              ctitle='MLD, March, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True, i_cb_subsamp=1,
              cfig_type=fig_type, lat_min=-80., lat_max=75., lpix=False)

if l_obs_mld:
    bp.plot("nproj")('nseas', 200., zmax_mld_atl, dz_mld, xlon, xlat, Xmld_obs[:,:],
                     cfignm=path_fig+'mld_obs_001_NSeas_march_'+CONFRUN, cpal='sst0', cbunit='m',
                     ctitle='MLD (obs., 0.01 crit.), March (Levitus 1980-1999)',
                     lkcont=True, cfig_type=fig_type, lforce_lim=True)

    bp.plot("2d")(xlon[0,:], xlat[:,ji_lat0], Xmld_obs[:,:], Xmask[0,:,:], 0., 600., 20.,
                  corca=vdic['ORCA'], lkcont=True, cpal='sst0',
                  cfignm=path_fig+'mld_obs_001_Global_march_'+CONFRUN, cbunit='m',
                  ctitle='MLD (obs., 0.01 crit.), March (Levitus 1980-1999)', lforce_lim=True, i_cb_subsamp=1,
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
bp.plot("nproj")('spstere', 100., 2000., dz_mld, xlon, xlat, mldr10[imnth,:,:],
                 cfignm=path_fig+'mld_ACC_september_'+CONFRUN, cpal='sst0', cbunit='m',
                 ctitle='MLD, September, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                 lkcont=True, cfig_type=fig_type,
                 lforce_lim=True)

bp.plot("2d")(xlon[0,:], xlat[:,ji_lat0], mldr10[imnth,:,:], Xmask[0,:,:], 0., 600., 20.,
              corca=vdic['ORCA'], lkcont=True, cpal='sst0',
              cfignm=path_fig+'mld_Global_september_'+CONFRUN, cbunit='m',
              ctitle='MLD, September, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True, i_cb_subsamp=1,
              cfig_type=fig_type, lat_min=-80., lat_max=75., lpix=False)


if l_obs_mld:
    bp.plot("2d")(xlon[0,:], xlat[:,ji_lat0], Xmld_obs[:,:], Xmask[0,:,:], 0., 600., 20.,
                  corca=vdic['ORCA'], lkcont=True, cpal='sst0',
                  cfignm=path_fig+'mld_obs_001_Global_september_'+CONFRUN, cbunit='m',
                  ctitle='MLD (obs, 0.01 crit.), March (Levitus 1980-1999)', lforce_lim=True, i_cb_subsamp=1,
                  cfig_type=fig_type, lat_min=-80., lat_max=75., lpix=False)










print '\n Bye!'

