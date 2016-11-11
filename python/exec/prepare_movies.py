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


#lfig0 = True
#lfig1 = True
#lfig2 = True


venv_needed = {'ORCA','RUN','DIAG_D','MM_FILE','NN_SST','NN_T','NN_S',
               'F_T_CLIM_3D_12','F_S_CLIM_3D_12','SST_CLIM_12','NN_SST_CLIM','NN_T_CLIM','NN_S_CLIM'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

tmin=-4.  ;  tmax=-tmin ;  dtemp = 0.25
smin=-1. ;  smax=-smin ;  dsali = 0.05

fig_type='png'



narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <NEMO grid_T file (1 year, monthyly)> <year>'; sys.exit(0)
cf_in = sys.argv[1]
cy=sys.argv[2] ; jy=int(cy)

path_fig = 'movies'

os.system("mkdir -p "+path_fig)


# 3D climatology :
# ------------

# Temperature
bt.chck4f(vdic['F_T_CLIM_3D_12'])
id_clim = Dataset(vdic['F_T_CLIM_3D_12'])
Tclim  = id_clim.variables[vdic['NN_T_CLIM']][:,:,:,:]; print '(has ',Tclim.shape[0],' time snapshots)\n'
id_clim.close()
[ nmn , nk0 , nj0 , ni0 ] = Tclim.shape

# Salinity
bt.chck4f(vdic['F_S_CLIM_3D_12'])
id_clim = Dataset(vdic['F_S_CLIM_3D_12'])
#Sclim  = id_clim.variables[vdic['NN_S_CLIM']][:,:,:,:]; print '(has ',Sclim.shape[0],' time snapshots)\n'
SSSclim  = id_clim.variables[vdic['NN_S_CLIM']][:,0,:,:]; print '(has ',SSSclim.shape[0],' time snapshots)\n'
id_clim.close()

print [ nmn , nk0 , nj0 , ni0 ]


# 2D SST obs :
print 'We use the following SST climatology:'; print vdic['SST_CLIM_12']
bt.chck4f(vdic['SST_CLIM_12'])
id_clim_sst = Dataset(vdic['SST_CLIM_12'])
SSTclim  = id_clim_sst.variables[vdic['NN_SST_CLIM']][:,:,:]; print '(has ',SSTclim.shape[0],' time snapshots)\n'
id_clim_sst.close()



# Getting land-sea mask and coordinates:
bt.chck4f(vdic['MM_FILE'])
id_mask = Dataset(vdic['MM_FILE'])
xlon  = id_mask.variables['nav_lon'][:,:]
xlat  = id_mask.variables['nav_lat'][:,:]
imask = id_mask.variables['tmask'][0,:,:,:]
id_mask.close()



# Getting SST, THETA and S in NEMO monthly grid_T file:
# -----------------------------------------------------

bt.chck4f(cf_in) ; id_in = Dataset(cf_in)

if vdic['NN_SST'] == 'thetao':
    SSTnemo = id_in.variables[vdic['NN_SST']][:,0,:,:]
else:
    SSTnemo = id_in.variables[vdic['NN_SST']][:,:,:]

SSSnemo = id_in.variables[vdic['NN_S']][:,0,:,:]

#Tnemo  = id_in.variables[vdic['NN_T']][:,:,:,:]
#print '(has ',Tnemo.shape[0],' time snapshots)\n'
#Snemo  = id_in.variables[vdic['NN_S']][:,:,:,:]
vdepth = id_in.variables['deptht'][:]
id_in.close()

#[ nt, nk, nj, ni ] = Tnemo.shape
[ nt, nj, ni ] = SSTnemo.shape
#if nk != nk0 or nj != nj0 or ni != ni0:
if nj != nj0 or ni != ni0:
    print 'ERROR: NEMO file do no agree in shape!'
    print '       clim => '+str(ni0)+', '+str(nj0)+', '+str(nk0),' ('+vdic['F_T_CLIM_3D_12']+')'
    print '       NEMO => '+str(ni)+', '+str(nj)+', '+str(nk)
    sys.exit(0)
if nt != 12:
    print 'ERROR: we expect 12 montly records in NEMO grid_T file!'
    sys.exit(0)

# Creating 1D long. and lat.:
ji_lat0 = nmp.argmax(xlat[nj-1])  ; #lolo
vlon = nmp.zeros(ni) ; vlon[:] = xlon[20,:]
vlat = nmp.zeros(nj) ; vlat[:] = xlat[:,ji_lat0]


###cy = '1990' ; #lulu

cv_dsst = 'dsst'
cv_dsss = 'dsss'

for jt in range(nt):

    cm = "%02d" % (jt+1)
    cdate = cy+cm

    # SST:
    bp.plot("2d")(vlon, vlat, SSTnemo[jt,:,:] - SSTclim[jt,:,:],
                  imask[0,:,:],  tmin, tmax, dtemp,
                  corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r',
                  cfignm=path_fig+'/'+cv_dsst+'_'+cdate,
                  cbunit='K', cfig_type=fig_type, lat_min=-65., lat_max=75.,
                  ctitle='SST (NEMO - obs) '+CONFRUN+' ('+cdate+')',
                  lforce_lim=True, i_cb_subsamp=2)

    # SSS:
    bp.plot("2d")(vlon, vlat, SSSnemo[jt,:,:] - SSSclim[jt,:,:],
                  imask[0,:,:],  smin, smax, dsali,
                  corca=vdic['ORCA'], lkcont=False, cpal='PiYG_r',
                  cfignm=path_fig+'/'+cv_dsss+'_'+cdate,
                  cbunit='PSU', cfig_type=fig_type, lat_min=-65., lat_max=75.,
                  ctitle='SSS (NEMO - obs) '+CONFRUN+' ('+cdate+')',
                  lforce_lim=True, i_cb_subsamp=2)



#for cv in [ cv_dsst, cv_dsss ]:
#    cfig_out = path_fig+'/'+cv+'_'+CONFRUN+'_'+cy+'.gif'
#    os.system("convert -delay 100 -loop 0 "+cv+"*.png "+cfig_out+" > out_conv_"+cv+".out")
