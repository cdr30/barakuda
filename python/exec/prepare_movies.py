#!/usr/bin/env python

#       B a r a K u d a
#
#  Prepare 2D maps (monthly) that will later become a GIF animation!
#  NEMO output and observations needed
#
#    L. Brodeau, november 2016

import sys
import os
import numpy as nmp

from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp


venv_needed = {'ORCA','RUN','DIAG_D','MM_FILE','NN_SST','NN_T','NN_S','NN_ICEF',
               'F_T_CLIM_3D_12','F_S_CLIM_3D_12','SST_CLIM_12','NN_SST_CLIM','NN_T_CLIM','NN_S_CLIM'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

tmin=-4.  ;  tmax=-tmin ;  dtemp = 0.25
smin=-1. ;  smax=-smin ;  dsali = 0.05

fig_type='png'



narg = len(sys.argv)
if narg < 4:
    print 'Usage: '+sys.argv[0]+' <NEMO file (1 year, monthyly)> <year> <var>'
    print '          with var one of "sst", "sss", "ice"'
    sys.exit(0)

cf_in = sys.argv[1]
cy    = sys.argv[2] ; jy=int(cy)
cvar  = sys.argv[3]

if not cvar in ['sst','sss','ice']:
    print 'ERROR (prepare_movies.py): variable '+cvar+' not supported yet!'
    sys.exit(0)

path_fig = 'movies'

os.system("mkdir -p "+path_fig)


# 3D climatology :
# ------------

# Temperature
#    bt.chck4f(vdic['F_T_CLIM_3D_12'])
#    id_clim = Dataset(vdic['F_T_CLIM_3D_12'])
#    Tclim  = id_clim.variables[vdic['NN_T_CLIM']][:,0,:,:]; print '(has ',Tclim.shape[0],' time snapshots)\n'
#    id_clim.close()
#    [ nmn , nk0 , nj0 , ni0 ] = Tclim.shape


# Salinity
if cvar == 'sss':
    bt.chck4f(vdic['F_S_CLIM_3D_12'])
    id_clim = Dataset(vdic['F_S_CLIM_3D_12'])
    Vclim  = id_clim.variables[vdic['NN_S_CLIM']][:,0,:,:]; print '(has ',Vclim.shape[0],' time snapshots)\n'
    id_clim.close()

# 2D SST obs :
if cvar == 'sst':
    bt.chck4f(vdic['SST_CLIM_12'])
    id_clim_sst = Dataset(vdic['SST_CLIM_12'])
    Vclim  = id_clim_sst.variables[vdic['NN_SST_CLIM']][:,:,:]; print '(has ',Vclim.shape[0],' time snapshots)\n'
    id_clim_sst.close()

# Sea-ice concentration :
# => no clim used!
if cvar == 'sss' or cvar == 'sst': [ nmn , nj0 , ni0 ] = Vclim.shape



# Getting land-sea mask and coordinates:
bt.chck4f(vdic['MM_FILE'])
id_mask = Dataset(vdic['MM_FILE'])
xlon  = id_mask.variables['nav_lon'][:,:]
xlat  = id_mask.variables['nav_lat'][:,:]
imask = id_mask.variables['tmask'][0,0,:,:]
id_mask.close()



# Getting NEMO variables:
# -----------------------

bt.chck4f(cf_in)

id_in = Dataset(cf_in)

if cvar == 'sst':
    if vdic['NN_SST'] == 'thetao' or vdic['NN_SST'] == 'votemper' :   #lolo:bad !!! should check shape!!!
        Vnemo = id_in.variables[vdic['NN_SST']][:,0,:,:]
    else:
        Vnemo = id_in.variables[vdic['NN_SST']][:,:,:]
    cv = 'dsst'

if cvar == 'sss':
    Vnemo = id_in.variables[vdic['NN_S']][:,0,:,:]
    cv = 'dsss'

if cvar == 'ice':
    Vnemo = id_in.variables[vdic['NN_ICEF']][:,:,:]

id_in.close()




[ nt, nj, ni ] = Vnemo.shape

if nt != 12:
    print 'ERROR (prepare_movies.py): we expect 12 montly records in NEMO grid_T file!'
    sys.exit(0)

if cvar == 'sss' or cvar == 'sst':
    if nj != nj0 or ni != ni0:
        print 'ERROR (prepare_movies.py): NEMO file and clim do no agree in shape!'
        print '       clim => '+str(ni0)+', '+str(nj0)+', '+str(nk0),' ('+vdic['F_T_CLIM_3D_12']+')'
        print '       NEMO => '+str(ni)+', '+str(nj)+', '+str(nk)
        sys.exit(0)
    # Creating 1D long. and lat.:
    ji_lat0 = nmp.argmax(xlat[nj-1])  ; #lolo
    vlon = nmp.zeros(ni) ; vlon[:] = xlon[20,:]
    vlat = nmp.zeros(nj) ; vlat[:] = xlat[:,ji_lat0]


if cvar == 'ice':
    # Extraoplating sea values on continents:
    bt.drown(Vnemo[:,:,:], imask, k_ew=2, nb_max_inc=10, nb_smooth=10)



for jt in range(nt):

    cm = "%02d" % (jt+1)
    cdate = cy+cm

    if cvar == 'sst':
        bp.plot("2d")(vlon, vlat, Vnemo[jt,:,:] - Vclim[jt,:,:],
                      imask[:,:],  tmin, tmax, dtemp,
                      corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r',
                      cfignm=path_fig+'/'+cv+'_'+cdate,
                      cbunit='K', cfig_type=fig_type, lat_min=-65., lat_max=75.,
                      ctitle='SST (NEMO - obs) '+CONFRUN+' ('+cdate+')',
                      lforce_lim=True, i_cb_subsamp=2)

    if cvar == 'sss':
        bp.plot("2d")(vlon, vlat, Vnemo[jt,:,:] - Vclim[jt,:,:],
                      imask[:,:],  smin, smax, dsali,
                      corca=vdic['ORCA'], lkcont=False, cpal='PiYG_r',
                      cfignm=path_fig+'/'+cv+'_'+cdate,
                      cbunit='PSU', cfig_type=fig_type, lat_min=-65., lat_max=75.,
                      ctitle='SSS (NEMO - obs) '+CONFRUN+' ('+cdate+')',
                      lforce_lim=True, i_cb_subsamp=2)

    
    if cvar == 'ice':

        # Extraoplating sea values on continents:
        bt.drown(Vnemo[jt,:,:], imask, k_ew=2, nb_max_inc=10, nb_smooth=10)

        # ICE north:
        cv = "icen"
        bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, Vnemo[jt,:,:],
                         cfignm=path_fig+'/'+cv+'_'+cdate, cpal='ice', cbunit='(frac.)',
                         ctitle='Ice frac. '+CONFRUN+' ('+cdate+')',
                         lkcont=True, cfig_type=fig_type, lforce_lim=True)

        cv = "ices"
        bp.plot("nproj")('spstere', 0., 1., 0.1, xlon, xlat, Vnemo[jt,:,:],
                         cfignm=path_fig+'/'+cv+'_'+cdate, cpal='ice', cbunit='(frac.)',
                         ctitle='Ice frac. '+CONFRUN+' ('+cdate+')',
                         lkcont=True, cfig_type=fig_type, lforce_lim=True)
