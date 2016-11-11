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
import barakuda_orca as bo
import barakuda_plot as bp


#lfig0 = True
#lfig1 = True
#lfig2 = True


venv_needed = {'ORCA','RUN','DIAG_D','MM_FILE','NN_SST','NN_T','NN_S',
               'F_T_CLIM_3D_12','F_S_CLIM_3D_12','SST_CLIM_12','NN_SST_CLIM','NN_T_CLIM','NN_S_CLIM'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

tmin=-4.  ;  tmax=-tmin ;  dtemp = 0.25
smin=-0.5 ;  smax=-smin ;  dsali = 0.025

fig_type='png'



narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <NEMO grid_T file (1 year, monthyly)> <year>'; sys.exit(0)
cf_in = sys.argv[1]
cy=sys.argv[2] ; jy=int(cy)

path_fig = 'flicks'

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
vlon = nmp.zeros(ni) ; vlon[:] = xlon[0,:]
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



for cv in [ cv_dsst, cv_dsss ]:
    cfig_out = path_fig+'/'+cv+'_'+CONFRUN+'_'+cy+'.gif'
    os.system("convert -delay 100 -loop 0 "+cv+"*.png "+cfig_out+" > out_conv_"+cv+".out")


sys.exit(0)


# Saving some array to avoid to call 'nmp.mean' all the time:

#Annual:
#Tnemo_annual = nmp.zeros((nk,nj,ni))
#Tnemo_annual[:,:,:] = nmp.mean(Tnemo[:,:,:,:], axis=0)
#Snemo_annual = nmp.zeros((nk,nj,ni))
#Snemo_annual[:,:,:] = nmp.mean(Snemo[:,:,:,:], axis=0)
#SSTnemo_annual = nmp.zeros((nj,ni))
#SSTnemo_annual[:,:] = nmp.mean(SSTnemo[:,:,:], axis=0)

#Tclim_annual = nmp.zeros((nk,nj,ni))
#Tclim_annual[:,:,:] = nmp.mean(Tclim[:,:,:,:], axis=0)
#Sclim_annual = nmp.zeros((nk,nj,ni))
#Sclim_annual[:,:,:] = nmp.mean(Sclim[:,:,:,:], axis=0)
#SSTclim_annual = nmp.zeros((nj,ni))
#SSTclim_annual[:,:] = nmp.mean(SSTclim[:,:,:], axis=0)


jk100  = bt.find_index_from_value(100.  , vdepth) ; print 'jk100  = ', jk100,  '=> ', vdepth[jk100]
jk1000 = bt.find_index_from_value(1000. , vdepth) ; print 'jk1000 = ', jk1000, '=> ', vdepth[jk1000]
jk3000 = bt.find_index_from_value(3000. , vdepth) ; print 'jk3000 = ', jk3000, '=> ', vdepth[jk3000]

tdj = [ jk100,   jk1000, jk3000  ]

tdd_true = [ str(int(round(vdepth[jk100])))+'m' , str(int(round(vdepth[jk1000])))+'m' , str(int(round(vdepth[jk3000])))+'m' ]
tdd      = [ '100m', '1000m', '3000m' ]

print '\n', tdd_true[:], '\n'













# Time for figures:
# -----------------



if lfig1:

    #                    SST / Reynolds
    # JFM
    bp.plot("2d")(vlon, vlat, SSTnemo_JFM[:,:] - SSTclim_JFM[:,:],
                  imask[0,:,:], tmin, tmax, dtemp,
                  corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dsst_JFM_'+CONFRUN+'_-_'+CC,
                  cbunit='K', cfig_type=fig_type,
                  ctitle='SST difference to '+CC+', JFM, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)
    # JAS
    bp.plot("2d")(vlon, vlat, SSTnemo_JAS[:,:] - SSTclim_JAS[:,:],
                  imask[0,:,:], tmin, tmax, dtemp,
                  corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dsst_JAS_'+CONFRUN+'_-_'+CC,
                  cbunit='K', cfig_type=fig_type,
                  ctitle='SST difference to '+CC+', JAS, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)

    # Annual
    bp.plot("2d")(vlon, vlat, SSTnemo_annual[:,:] - SSTclim_annual[:,:],
                  imask[0,:,:],  tmin, tmax, dtemp,
                  corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dsst_annual_'+CONFRUN+'_-_'+CC,
                  cbunit='K', cfig_type=fig_type,
                  ctitle='SST difference to '+CC+', '+CONFRUN+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)





    # Temperature 100m, 1000m... / climatology

    for jd in range(nmp.size(tdj)):
        jdepth = tdj[jd] ; cdepth = tdd[jd] ; cdepth_true = tdd_true[jd]

        print '\n Treating depth '+str(vdepth[jdepth])+' !!!'


        if jd < 1:
            # JFM
            bp.plot("2d")(vlon, vlat, Tnemo_JFM[jdepth,:,:] - Tclim_JFM[jdepth,:,:],
                          imask[jdepth,:,:], tmin, tmax, dtemp,
                          corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dT_JFM_'+cdepth+'_'+CONFRUN+'_-_'+CC,
                          cbunit='K', cfig_type=fig_type,
                          ctitle='Temperature diff. to '+CC+' at '+cdepth_true+', JFM, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                          lforce_lim=True)
            # JAS
            bp.plot("2d")(vlon, vlat, Tnemo_JAS[jdepth,:,:] - Tclim_JAS[jdepth,:,:],
                          imask[jdepth,:,:], tmin, tmax, dtemp,
                          corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dT_JAS_'+cdepth+'_'+CONFRUN+'_-_'+CC,
                          cbunit='K', cfig_type=fig_type,
                          ctitle='Temperature diff. to '+CC+' at '+cdepth_true+', JAS, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                          lforce_lim=True)

        # Annual
        bp.plot("2d")(vlon, vlat, Tnemo_annual[jdepth,:,:] - Tclim_annual[jdepth,:,:],
                      imask[jdepth,:,:], tmin, tmax, dtemp,
                      corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dT_annual_'+cdepth+'_'+CONFRUN+'_-_'+CC,
                      cbunit='K', cfig_type=fig_type,
                      ctitle='Temperature diff. to '+CC+' at '+cdepth_true+', '+CONFRUN+' ('+cy1+'-'+cy2+')',
                      lforce_lim=True)






    #                   S S S
    # JFM
    bp.plot("2d")(vlon, vlat, Snemo_JFM[0,:,:] - Sclim_JFM[0,:,:],
                  imask[0,:,:], smin, smax, dsali, cpal='PiYG_r',
                  corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dsss_JFM_'+CONFRUN+'_-_'+CC,
                  cbunit='PSU', cfig_type=fig_type,
                  ctitle='SSS difference to '+CC+', JFM, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)
    # JAS
    bp.plot("2d")(vlon, vlat, Snemo_JAS[0,:,:] - Sclim_JAS[0,:,:],
                  imask[0,:,:], smin, smax, dsali, cpal='PiYG_r',
                  corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dsss_JAS_'+CONFRUN+'_-_'+CC,
                  cbunit='PSU', cfig_type=fig_type,
                  ctitle='SSS difference to '+CC+', JAS, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)
    # Annual
    bp.plot("2d")(vlon, vlat, Snemo_annual[0,:,:] - Sclim_annual[0,:,:],
                  imask[0,:,:], smin, smax, dsali, cpal='PiYG_r',
                  corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dsss_annual_'+CONFRUN+'_-_'+CC,
                  cbunit='PSU', cfig_type=fig_type,
                  ctitle='SSS difference to '+CC+', '+CONFRUN+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)






    #                   Salinity 100m / climatology
    for jd in range(nmp.size(tdj)):
        jdepth = tdj[jd] ; cdepth = tdd[jd] ; cdepth_true = tdd_true[jd]

        if jd < 1:
            # JFM
            bp.plot("2d")(vlon, vlat, Snemo_JFM[jdepth,:,:] - Sclim_JFM[jdepth,:,:],
                          imask[jdepth,:,:], smin, smax, dsali, cpal='PiYG_r',
                          corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dS_JFM_'+cdepth+'_'+CONFRUN+'_-_'+CC,
                          cbunit='PSU', cfig_type=fig_type,
                          ctitle='Salinity diff. to '+CC+' at '+cdepth_true+', JFM, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                          lforce_lim=True)
            # JAS
            bp.plot("2d")(vlon, vlat, Tnemo_JAS[jdepth,:,:] - Tclim_JAS[jdepth,:,:],
                          imask[jdepth,:,:], smin, smax, dsali, cpal='PiYG_r',
                          corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dS_JAS_'+cdepth+'_'+CONFRUN+'_-_'+CC,
                          cbunit='PSU', cfig_type=fig_type,
                          ctitle='Salinity diff. to '+CC+' at '+cdepth_true+', JAS, '+CONFRUN+' ('+cy1+'-'+cy2+')',
                          lforce_lim=True)

        # Annual
        bp.plot("2d")(vlon, vlat, Tnemo_annual[jdepth,:,:] - Tclim_annual[jdepth,:,:],
                      imask[jdepth,:,:], smin, smax, dsali, cpal='PiYG_r',
                      corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dS_annual_'+cdepth+'_'+CONFRUN+'_-_'+CC,
                      cbunit='PSU', cfig_type=fig_type,
                      ctitle='Salinity diff. to '+CC+' at '+cdepth_true+', '+CONFRUN+' ('+cy1+'-'+cy2+')',
                      lforce_lim=True)





if lfig2: # Temperature and salinity vertical sections:

    vcoup_ind = bo.coor2ind('atl_vert', xlon, xlat); jic = vcoup_ind[0]     # Atlantic transect


    #  N E M O
    #  ~~~~~~~
    bp.plot("vert_section")(xlat[:,jic], vdepth, Tnemo_annual[:,:,jic],
                            imask[:,:,jic], -1., 25., 1., cpal='mld', xmin=-75., xmax=65., dx=15.,
                            cfignm=path_fig+'section_temp_'+CONFRUN, cbunit=r'$^{\circ}$C', cxunit=r'Latitude ($^{\circ}$N)',
                            czunit='Depth (m)', ctitle='Temperature, ('+cy1+'-'+cy2+'), '+CONFRUN+', lon = 30W',
                            cfig_type=fig_type, lforce_lim=True)

    bp.plot("vert_section")(xlat[:,jic], vdepth, Snemo_annual[:,:,jic],
                            imask[:,:,jic], 33.9, 35.9, 0.1, cpal='mld', xmin=-75., xmax=65., dx=15.,
                            cfignm=path_fig+'section_sali_'+CONFRUN, cbunit='PSU', cxunit=r'Latitude ($^{\circ}$N)',
                            czunit='Depth (m)', ctitle='Salinity, ('+cy1+'-'+cy2+'), '+CONFRUN+', lon = 30W',
                            cfig_type=fig_type, lforce_lim=True)


    #
    #  L E V I T U S
    #  ~~~~~~~~~~~~~
    bp.plot("vert_section")(xlat[:,jic], vdepth, Tclim_annual[:,:,jic],
                            imask[:,:,jic], -1., 25., 1., cpal='mld', xmin=-75., xmax=65., dx=15.,
                            cfignm=path_fig+'section_temp_'+CC, cbunit=r'$^{\circ}$C',
                            cxunit=r'Latitude ($^{\circ}$N)',
                            czunit='Depth (m)', ctitle='Temperature, '+CC+', lon = 30W',
                            cfig_type=fig_type, lforce_lim=True)
    #
    bp.plot("vert_section")(xlat[:,jic], vdepth, Sclim_annual[:,:,jic],
                            imask[:,:,jic], 33.9, 35.9, 0.1, cpal='mld', xmin=-75., xmax=65., dx=15.,
                            cfignm=path_fig+'section_sali_'+CC, cbunit='PSU',
                            cxunit=r'Latitude ($^{\circ}$N)',
                            czunit='Depth (m)', ctitle='Salinity, '+CC+', lon = 30W',
                            cfig_type=fig_type, lforce_lim=True)



print '\n Bye!'
