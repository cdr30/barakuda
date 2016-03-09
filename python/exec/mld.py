#!/usr/bin/env python

# L. Brodeau, november 2009

import sys
import os
import numpy as nmp
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_orca as bo
import barakuda_plot as bp
import barakuda_physics as bphys

ldebug = False

zmax_mld_atl = 1600. ; dz_mld = 100.

ORCA = os.getenv('ORCA')
if ORCA == None: print 'The ORCA environement variable is no set'; sys.exit(0)
RUN = os.getenv('RUN')
if RUN == None: print 'The RUN environement variable is no set'; sys.exit(0)
DIAG_D = os.getenv('DIAG_D')
if DIAG_D == None: print 'The DIAG_D environement variable is no set'; sys.exit(0)

COMP2D = os.getenv('COMP2D')
if COMP2D == None: print 'The COMP2D environement variable is no set'; sys.exit(0)

NN_MLD = os.getenv('NN_MLD')
if NN_MLD == None: print 'The NN_MLD environement variable is no set'; sys.exit(0)
NN_S = os.getenv('NN_S')
if NN_S == None: print 'The NN_S environement variable is no set'; sys.exit(0)
NN_T = os.getenv('NN_T')
if NN_T == None: print 'The NN_T environement variable is no set'; sys.exit(0)

NN_S_CLIM = os.getenv('NN_S_CLIM')
if NN_S_CLIM == None: print 'The NN_S_CLIM environement variable is no set'; sys.exit(0)
NN_T_CLIM = os.getenv('NN_T_CLIM')
if NN_T_CLIM == None: print 'The NN_T_CLIM environement variable is no set'; sys.exit(0)



F_T_CLIM_3D_12 = os.getenv('F_T_CLIM_3D_12')
if F_T_CLIM_3D_12 == None: print 'The F_T_CLIM_3D_12 environement variable is no set\n'
F_S_CLIM_3D_12 = os.getenv('F_S_CLIM_3D_12')
if F_S_CLIM_3D_12 == None: print 'The F_S_CLIM_3D_12 environement variable is no set\n'


print ' ORCA = '+ORCA
print ' RUN = '+RUN
print ' DIAG_D = '+DIAG_D
print ' COMP2D = '+COMP2D



if 'ORCA2' in ORCA:
    ji_lat0 = 132
elif 'ORCA1' in ORCA:
    ji_lat0 = 265
else:
    print 'FIX ME!!! ssh.py => dont know ji_lat0 for conf '+ORCA+' !!!'; sys.exit(0)




l_obs_mld = False
if (not F_T_CLIM_3D_12 == None) and (not F_T_CLIM_3D_12 == None):
    l_obs_mld = True ; print 'Since obs Temp and Sali here will compute observed MLD!'




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
bt.chck4f(cf_mesh_mask)
id_mm = Dataset(cf_mesh_mask)
xlon   = id_mm.variables['glamt'][0,:,:] ; xlat = id_mm.variables['gphit'][0,:,:]
Xmask = id_mm.variables['tmask'][0,:,:,:]
vlev  = id_mm.variables['gdept_1d'][0,:]
id_mm.close()

nk = len(vlev)



#  Getting NEMO mean monthly climatology of MLD coverage:
cf_nemo_mnmc = DIAG_D+'/clim/mclim_'+CONFRUN+'_'+cy1+'-'+cy2+'_grid_T.nc4'

bt.chck4f(cf_nemo_mnmc)
id_nemo = Dataset(cf_nemo_mnmc)
mldr10   = id_nemo.variables[NN_MLD][:,:,:]
id_nemo.close()

[ nt, nj, ni ] = mldr10.shape ; print ' Shape of MLD :', nt, nj, ni, '\n'




if l_obs_mld:

    # Getting 3D+T 12-month climatology of T and S
    # --------------------------------------------

    # Temperature
    bt.chck4f(F_T_CLIM_3D_12) ; id_clim = Dataset(F_T_CLIM_3D_12)
    Tclim  = id_clim.variables[NN_T_CLIM][:,:,:,:]; print '(has ',Tclim.shape[0],' time snapshots)\n'
    id_clim.close()

    # Salinity
    bt.chck4f(F_S_CLIM_3D_12) ; id_clim = Dataset(F_S_CLIM_3D_12)
    Sclim  = id_clim.variables[NN_S_CLIM][:,:,:,:]; print '(has ',Sclim.shape[0],' time snapshots)\n'
    id_clim.close()

    [ nmn , nk0, nj0 , ni0 ] = Tclim.shape

    if nj != nj0 or ni != ni0 or nk != nk0:
        print 'ERROR: 3D clim and NEMO file do no agree in shape!'
        print '       clim => '+str(ni0)+', '+str(nj0)+', '+str(nk0)+' ('+F_T_CLIM_3D_12+')'
        print '       NEMO => '+str(ni)+', '+str(nj)+', '+str(nk)
        sys.exit(0)





    #############
    # M A R C H #
    #############

    imnth = 2 ; # march

    # Computing sigma0 3D field:
    Sigma0 = bphys.sigma0(Tclim[imnth,:,:,:], Sclim[imnth,:,:,:])*Xmask[:,:,:]

    if ldebug: Sigma0_nemo = bphys.sigma0(Tnemo[imnth:,:,:], Snemo[imnth,:,:,:])*Xmask[:,:,:]



    Xmld_obs = nmp.zeros(nj*ni) ; Xmld_obs.shape = [nj,ni]
    mmask = nmp.zeros(nj*ni)    ; mmask.shape = [nj,ni]



    mmask[:,:] = 1.
    for jk in range(nk-1):
        zz = vlev[jk]
        Sigma0[jk,:,:]   = Sigma0[jk,:,:]*mmask[:,:]
        Sigma0[jk+1,:,:] = Sigma0[jk+1,:,:]*mmask[:,:]
        ijloc = nmp.where( Sigma0[jk+1,:,:] - 0.01 > Sigma0[jk,:,:] )
        #lolo Xmld_obs[ijloc] = zz ; # lolo, gives bug!
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
                         cfignm=path_fig+'mld_NEMO_001_NSeas_march_'+CONFRUN+'_vs_'+COMP2D, cpal='sst0', cbunit='m',
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
              corca=ORCA, lkcont=True, cpal='sst0',
              cfignm=path_fig+'mld_Global_march_'+CONFRUN, cbunit='m',
              ctitle='MLD, March, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True, i_cb_subsamp=1,
              cfig_type=fig_type, lat_min=-80., lat_max=75., lpix=False)

if l_obs_mld:
    bp.plot("nproj")('nseas', 200., zmax_mld_atl, dz_mld, xlon, xlat, Xmld_obs[:,:],
                     cfignm=path_fig+'mld_obs_001_NSeas_march_'+CONFRUN, cpal='sst0', cbunit='m',
                     ctitle='MLD (obs., 0.01 crit.), March (Levitus 1980-1999)',
                     lkcont=True, cfig_type=fig_type, lforce_lim=True)

    bp.plot("2d")(xlon[0,:], xlat[:,ji_lat0], Xmld_obs[:,:], Xmask[0,:,:], 0., 600., 20.,
                  corca=ORCA, lkcont=True, cpal='sst0',
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
              corca=ORCA, lkcont=True, cpal='sst0',
              cfignm=path_fig+'mld_Global_september_'+CONFRUN, cbunit='m',
              ctitle='MLD, September, '+CONFRUN+' ('+cy1+'-'+cy2+')', lforce_lim=True, i_cb_subsamp=1,
              cfig_type=fig_type, lat_min=-80., lat_max=75., lpix=False)


if l_obs_mld:
    bp.plot("2d")(xlon[0,:], xlat[:,ji_lat0], Xmld_obs[:,:], Xmask[0,:,:], 0., 600., 20.,
                  corca=ORCA, lkcont=True, cpal='sst0',
                  cfignm=path_fig+'mld_obs_001_Global_september_'+CONFRUN, cbunit='m',
                  ctitle='MLD (obs, 0.01 crit.), March (Levitus 1980-1999)', lforce_lim=True, i_cb_subsamp=1,
                  cfig_type=fig_type, lat_min=-80., lat_max=75., lpix=False)










print '\n Bye!'

