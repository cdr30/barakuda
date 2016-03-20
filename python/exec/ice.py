#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate stereographic plot of sea-ice extent at both poles
#
#       L. Brodeau, 2009

import sys
import numpy as nmp
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp

venv_needed = {'ORCA','RUN','DIAG_D','MM_FILE','COMP2D','ICE_CLIM_12','NN_ICEF_CLIM','NN_ICEF','NN_ICET','FILE_ICE_SUFFIX'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

crun = vdic['RUN']
CONFRUN = vdic['ORCA']+'-'+crun

path_fig='./'
fig_type='png'



narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <year1> <year2>'; sys.exit(0)
cy1 = sys.argv[1] ; cy2=sys.argv[2]; jy1=int(cy1); jy2=int(cy2)


jy1_clim = jy1 ; jy2_clim = jy2

print ' => mean on the clim : ', jy1_clim, jy2_clim, '\n'


# Getting coordinates and mask :
bt.chck4f(vdic['MM_FILE'])
id_mm = Dataset(vdic['MM_FILE'])
xlon    = id_mm.variables['nav_lon'][:,:]
xlat    = id_mm.variables['nav_lat'][:,:]
xmask   = id_mm.variables['tmask'][0,0,:,:]
id_mm.close()


# Getting obs:
bt.chck4f(vdic['ICE_CLIM_12'])
id_clim = Dataset(vdic['ICE_CLIM_12'])
xclim03 = id_clim.variables[vdic['NN_ICEF_CLIM']][2,:,:]
xclim09 = id_clim.variables[vdic['NN_ICEF_CLIM']][8,:,:]
id_clim.close()


#  Getting NEMO mean monthly climatology of sea-ice coverage:
cf_nemo_mnmc = vdic['DIAG_D']+'/clim/mclim_'+CONFRUN+'_'+cy1+'-'+cy2+'_'+vdic['FILE_ICE_SUFFIX']+'.nc4'

bt.chck4f(cf_nemo_mnmc)
id_ice = Dataset(cf_nemo_mnmc)
xnemo_frac_03   = id_ice.variables[vdic['NN_ICEF']][2,:,:]
xnemo_frac_09   = id_ice.variables[vdic['NN_ICEF']][8,:,:]
xnemo_thic_03   = id_ice.variables[vdic['NN_ICET']][2,:,:]
xnemo_thic_09   = id_ice.variables[vdic['NN_ICET']][8,:,:]
id_ice.close()

[ nj, ni ] = xnemo_frac_03.shape ; print ' Shape of sea-ice :', nj, ni, '\n'


# Extraoplating sea values on continents:
bt.drown(xnemo_frac_03, xmask, k_ew=2, nb_max_inc=10, nb_smooth=10)
bt.drown(xnemo_frac_09, xmask, k_ew=2, nb_max_inc=10, nb_smooth=10)
bt.drown(xnemo_thic_03, xmask, k_ew=2, nb_max_inc=10, nb_smooth=10)
bt.drown(xnemo_thic_09, xmask, k_ew=2, nb_max_inc=10, nb_smooth=10)
bt.drown(xclim03, xmask, k_ew=2, nb_max_inc=10, nb_smooth=10)
bt.drown(xclim09, xmask, k_ew=2, nb_max_inc=10, nb_smooth=10)




# Time for figures:       
# -----------------
#
# Extending to 90S: (from 78 to 90):
#
js_ext = 12; nje = nj + js_ext
xlat0       = nmp.zeros((nje,ni))
xlon0       = nmp.zeros((nje,ni))
xnemo_frac_030    = nmp.zeros((nje,ni))
xnemo_frac_090    = nmp.zeros((nje,ni))
xnemo_thic_030    = nmp.zeros((nje,ni))
xnemo_thic_090    = nmp.zeros((nje,ni))
xclim030    = nmp.zeros((nje,ni))
xclim090    = nmp.zeros((nje,ni))

#
for jj in range(js_ext):
    xlat0[jj,:] = -90. + float(jj)
    xlon0[jj,:] = xlon[0,:]
    xnemo_frac_030[jj,:] = xnemo_frac_03[1,:] #persistence
    xnemo_frac_090[jj,:] = xnemo_frac_09[1,:] #persistence
    xnemo_thic_030[jj,:] = xnemo_thic_03[1,:] #persistence
    xnemo_thic_090[jj,:] = xnemo_thic_09[1,:] #persistence
    xclim030[jj,:] = xclim03[1,:] #persistence
    xclim090[jj,:] = xclim09[1,:] #persistence

xlat0[js_ext:nje,:] = xlat[:,:]
xlon0[js_ext:nje,:] = xlon[:,:]
xnemo_frac_030[js_ext:nje,:] = xnemo_frac_03[:,:]
xnemo_frac_090[js_ext:nje,:] = xnemo_frac_09[:,:]
xnemo_thic_030[js_ext:nje,:] = xnemo_thic_03[:,:]
xnemo_thic_090[js_ext:nje,:] = xnemo_thic_09[:,:]
xclim030[js_ext:nje,:] = xclim03[:,:]
xclim090[js_ext:nje,:] = xclim09[:,:]



ratio = 1.

#if vdic['COMP2D'] == 'CLIM': ratio = 100.
if xclim03.max()>90.:
    ratio = 100.


#DEBUG:
if False:
    bp.plot("nproj")('spstere', 0., 1., 0.1, xlon0, xlat0, xclim090[:,:]/ratio,
                     cfignm=path_fig+'sea-ice_SP_sept_obs', cpal='ice', cbunit='(frac.)',
                     ctitle='Ice fraction, Sept., obs.',
                     lkcont=True, cfig_type=fig_type, 
                     lforce_lim=True)
    
    bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, xclim03[:,:]/ratio,
                     cfignm=path_fig+'sea-ice_NP_march_obs', cpal='ice', cbunit='(frac.)',
                     ctitle='Ice fraction, March, obs.',
                     lkcont=True, cfig_type=fig_type, 
                     lforce_lim=True)

    sys.exit(0)
#DEBUG.



# September
# ~~~~~~~~~
if vdic['COMP2D'] == 'CLIM':
    ctit_clim = 'Ice fraction, Sept., obs.'
else:
    ctit_clim = 'Ice fraction, Sept., '+vdic['COMP2D']



# Nordic Seas:
bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, xnemo_frac_09[:,:],
                 cfignm=path_fig+'sea-ice_NP_sept_'+CONFRUN, cpal='ice', cbunit='(frac.)',
                 ctitle='Ice fraction, Sept., '+crun+' ('+cy1+'-'+cy2+')',
                 lkcont=True, cfig_type=fig_type, 
                 lforce_lim=True)

bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, xclim09[:,:]/ratio,
                 cfignm=path_fig+'sea-ice_NP_sept_'+vdic['COMP2D'], cpal='ice', cbunit='(frac.)',
                 ctitle=ctit_clim,
                 lkcont=True, cfig_type=fig_type,
                 lforce_lim=True)
    
bp.plot("nproj")('npol2', 0., 10., 0.25, xlon, xlat, xnemo_thic_09[:,:],
                 cfignm=path_fig+'sea-ice_thickness_NP_sept_'+CONFRUN, cpal='cubehelix_r', cbunit='(m)',
                 ctitle='Ice thickness, Sept., '+crun+' ('+cy1+'-'+cy2+')',
                 lkcont=True, cfig_type=fig_type, i_cb_subsamp=2,
                 lforce_lim=True)



    
# September, big South:
bp.plot("nproj")('spstere', 0., 1., 0.1, xlon0, xlat0, xnemo_frac_090[:,:],
                 cfignm=path_fig+'sea-ice_SP_sept_'+CONFRUN, cpal='ice', cbunit='(frac.)',
                 ctitle='Ice fraction, Sept., '+crun+' ('+cy1+'-'+cy2+')',
                 lkcont=True, cfig_type=fig_type, 
                 lforce_lim=True)


bp.plot("nproj")('spstere', 0., 1., 0.1, xlon0, xlat0, xclim090[:,:]/ratio,
                 cfignm=path_fig+'sea-ice_SP_sept_'+vdic['COMP2D'], cpal='ice', cbunit='(frac.)',
                 ctitle=ctit_clim,
                 lkcont=True, cfig_type=fig_type, 
                 lforce_lim=True)

bp.plot("nproj")('spstere', 0., 5., 0.1, xlon0, xlat0, xnemo_thic_090[:,:],
                 cfignm=path_fig+'sea-ice_thickness_SP_sept_'+CONFRUN, cpal='cubehelix_r', cbunit='(m)',
                 ctitle='Ice thickness, Sept., '+crun+' ('+cy1+'-'+cy2+')',
                 lkcont=True, cfig_type=fig_type, i_cb_subsamp=2,
                 lforce_lim=True)










# March:
# ~~~~~~

if vdic['COMP2D'] == 'CLIM':
    ctit_clim = 'Ice fraction, March, obs.'
else:
    ctit_clim = 'Ice fraction, March, '+vdic['COMP2D']



# Nordic Seas:

#print "  -- doing stereographic projection of sea-ice fraction / March / N.E."

bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, xnemo_frac_03[:,:],
              cfignm=path_fig+'sea-ice_NP_march_'+CONFRUN, cpal='ice', cbunit='(frac.)',
              ctitle='Ice fraction, March, '+crun+' ('+cy1+'-'+cy2+')',
              lkcont=True, cfig_type=fig_type, 
              lforce_lim=True)

bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, xclim03[:,:]/ratio,
                cfignm=path_fig+'sea-ice_NP_march_'+vdic['COMP2D'], cpal='ice', cbunit='(frac.)',
                ctitle=ctit_clim,
                lkcont=True, cfig_type=fig_type, 
                lforce_lim=True)

#print "  -- doing stereographic projection of sea-ice thickness / March / N.E."

bp.plot("nproj")('npol2', 0., 10., 0.25, xlon, xlat, xnemo_thic_03[:,:],
              cfignm=path_fig+'sea-ice_thickness_NP_march_'+CONFRUN, cpal='cubehelix_r', cbunit='(m)',
              ctitle='Ice thickness, March, '+crun+' ('+cy1+'-'+cy2+')',
              lkcont=True, cfig_type=fig_type, i_cb_subsamp=2,
              lforce_lim=True)




# Big south:

bp.plot("nproj")('spstere', 0., 1., 0.1, xlon0, xlat0, xnemo_frac_030[:,:],
              cfignm=path_fig+'sea-ice_SP_march_'+CONFRUN, cpal='ice', cbunit='(frac.)',
              ctitle='Ice fraction, March, '+crun+' ('+cy1+'-'+cy2+')',
              lkcont=True, cfig_type=fig_type, 
              lforce_lim=True)

bp.plot("nproj")('spstere', 0., 1., 0.1, xlon0, xlat0, xclim030[:,:]/ratio,
                cfignm=path_fig+'sea-ice_SP_march_'+vdic['COMP2D'], cpal='ice', cbunit='(frac.)',
                ctitle=ctit_clim,
                lkcont=True, cfig_type=fig_type, 
                lforce_lim=True)

bp.plot("nproj")('spstere', 0., 5., 0.1, xlon0, xlat0, xnemo_thic_030[:,:],
              cfignm=path_fig+'sea-ice_thickness_SP_march_'+CONFRUN, cpal='cubehelix_r', cbunit='(m)',
              ctitle='Ice thickness, March, '+crun+' ('+cy1+'-'+cy2+')',
              lkcont=True, cfig_type=fig_type, i_cb_subsamp=2,
              lforce_lim=True)


print '\n Bye!'

