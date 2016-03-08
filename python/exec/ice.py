#!/usr/bin/env python

# L. Brodeau, november 2009

import sys
import os
import numpy as nmp
from netCDF4 import Dataset

import barakuda_plot as bp
import barakuda_tool as bt



ORCA = os.getenv('ORCA')
if ORCA == None: print 'The ORCA environement variable is no set'; sys.exit(0)
RUN = os.getenv('RUN')
if RUN == None: print 'The RUN environement variable is no set'; sys.exit(0)
DIAG_D = os.getenv('DIAG_D')
if DIAG_D == None: print 'The DIAG_D environement variable is no set'; sys.exit(0)
COMP2D = os.getenv('COMP2D')
if COMP2D == None: print 'The COMP2D environement variable is no set'; sys.exit(0)

ICE_CLIM_12 = os.getenv('ICE_CLIM_12')
if ICE_CLIM_12 == None: print 'The ICE_CLIM_12 environement variable is no set'; sys.exit(0)
NN_ICEF_CLIM = os.getenv('NN_ICEF_CLIM')
if NN_ICEF_CLIM == None: print 'The NN_ICEF_CLIM environement variable is no set'; sys.exit(0)

NN_ICEF = os.getenv('NN_ICEF')
if NN_ICEF == None: print 'The NN_ICEF environement variable is no set'; sys.exit(0)

FILE_ICE_SUFFIX = os.getenv('FILE_ICE_SUFFIX')
if FILE_ICE_SUFFIX == None: print 'The FILE_ICE_SUFFIX environement variable is no set'; sys.exit(0)

print ' ORCA = '+ORCA
print ' RUN = '+RUN
print ' DIAG_D = '+DIAG_D
print ' COMP2D = '+COMP2D
print ' Obs: ICE_CLIM_12 = '+ICE_CLIM_12
print ' In obs: NN_ICEF_CLIM = '+NN_ICEF_CLIM
print ' In NEMO: FILE_ICE_SUFFIX = '+FILE_ICE_SUFFIX
print ' In NEMO: NN_ICEF = '+NN_ICEF



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
bt.chck4f(cf_mesh_mask) ; id_coor = Dataset(cf_mesh_mask)
xlon   = id_coor.variables['nav_lon'][:,:] ; xlat   = id_coor.variables['nav_lat'][:,:]
id_coor.close()


# Getting obs:
bt.chck4f(ICE_CLIM_12) ; id_clim = Dataset(ICE_CLIM_12)
xclim03 = id_clim.variables[NN_ICEF_CLIM][2,:,:]
xclim09 = id_clim.variables[NN_ICEF_CLIM][8,:,:]
id_clim.close()



# Getting land-sea mask:
#-----------------------
xmask = 0
id_mask = Dataset(cf_mesh_mask)
xmask = id_mask.variables['tmask'][0,0,:,:]
id_mask.close()




#  Getting NEMO mean monthly climatology of sea-ice coverage:
cf_nemo_mnmc = DIAG_D+'/clim/mclim_'+CONFRUN+'_'+cy1+'-'+cy2+'_'+FILE_ICE_SUFFIX+'.nc4'

bt.chck4f(cf_nemo_mnmc)
id_ice = Dataset(cf_nemo_mnmc)
xnemo03   = id_ice.variables[NN_ICEF][2,:,:]
xnemo09   = id_ice.variables[NN_ICEF][8,:,:]
id_ice.close()

[ nj, ni ] = xnemo03.shape ; print ' Shape of sea-ice :', nj, ni, '\n'


# Extraoplating sea values on continents:
bt.drown(xnemo03, xmask, k_ew=2, nb_max_inc=10, nb_smooth=10)
bt.drown(xnemo09, xmask, k_ew=2, nb_max_inc=10, nb_smooth=10)
bt.drown(xclim03, xmask, k_ew=2, nb_max_inc=10, nb_smooth=10)
bt.drown(xclim09, xmask, k_ew=2, nb_max_inc=10, nb_smooth=10)




# Time for figures:       
# -----------------
#
# Extending to 90S: (from 78 to 90):
#
js_ext = 12; nje = nj + js_ext
xlat0     = nmp.zeros(nje*ni); xlat0.shape     = [ nje, ni ]
xlon0     = nmp.zeros(nje*ni); xlon0.shape     = [ nje, ni ]
xnemo030    = nmp.zeros(nje*ni); xnemo030.shape    = [ nje, ni ]
xnemo090    = nmp.zeros(nje*ni); xnemo090.shape    = [ nje, ni ]
xclim030    = nmp.zeros(nje*ni); xclim030.shape    = [ nje, ni ]
xclim090    = nmp.zeros(nje*ni); xclim090.shape    = [ nje, ni ]

#
for jj in range(js_ext):
    xlat0[jj,:] = -90. + float(jj)
    xlon0[jj,:] = xlon[0,:]
    xnemo030[jj,:] = xnemo03[1,:] #persistence
    xnemo090[jj,:] = xnemo09[1,:] #persistence
    xclim030[jj,:] = xclim03[1,:] #persistence
    xclim090[jj,:] = xclim09[1,:] #persistence

xlat0[js_ext:nje,:] = xlat[:,:]
xlon0[js_ext:nje,:] = xlon[:,:]
xnemo030[js_ext:nje,:] = xnemo03[:,:]
xnemo090[js_ext:nje,:] = xnemo09[:,:]
xclim030[js_ext:nje,:] = xclim03[:,:]
xclim090[js_ext:nje,:] = xclim09[:,:]



ratio = 1.

#if COMP2D == 'CLIM': ratio = 100.
if xclim03.max()>90.:
    ratio = 100.


#DEBUG:
if False:
    bp.plot("nproj")('spstere', 0., 1., 0.1, xlon0, xlat0, xclim090[:,:]/ratio,
                     cfignm=path_fig+'sea-ice_SP_sept_obs', cpal='ice', cbunit='frac.',
                     ctitle='Sea-Ice, Sept., obs.',
                     lkcont=True, cfig_type=fig_type, 
                     lforce_lim=True)
    
    bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, xclim03[:,:]/ratio,
                     cfignm=path_fig+'sea-ice_NP_march_obs', cpal='ice', cbunit='frac.',
                     ctitle='Sea-Ice, March, obs.',
                     lkcont=True, cfig_type=fig_type, 
                     lforce_lim=True)

    sys.exit(0)
#DEBUG.



# September
# ~~~~~~~~~
if COMP2D == 'CLIM':
    ctit_clim = 'Sea-Ice, Sept., obs.'
else:
    ctit_clim = 'Sea-Ice, Sept., '+COMP2D



# Nordic Seas:
bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, xnemo09[:,:],
                 cfignm=path_fig+'sea-ice_NP_sept_'+CONFRUN, cpal='ice', cbunit='frac.',
                 ctitle='Sea-Ice, Sept., '+CONFRUN+' ('+cy1+'-'+cy2+')',
                 lkcont=True, cfig_type=fig_type, 
                 lforce_lim=True)

bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, xclim09[:,:]/ratio,
                 cfignm=path_fig+'sea-ice_NP_sept_'+COMP2D, cpal='ice', cbunit='frac.',
                 ctitle=ctit_clim,
                 lkcont=True, cfig_type=fig_type,
                 lforce_lim=True)
    
    
# September, big South:
bp.plot("nproj")('spstere', 0., 1., 0.1, xlon0, xlat0, xnemo090[:,:],
                 cfignm=path_fig+'sea-ice_SP_sept_'+CONFRUN, cpal='ice', cbunit='frac.',
                 ctitle='Sea-Ice, Sept., '+CONFRUN+' ('+cy1+'-'+cy2+')',
                 lkcont=True, cfig_type=fig_type, 
                 lforce_lim=True)


bp.plot("nproj")('spstere', 0., 1., 0.1, xlon0, xlat0, xclim090[:,:]/ratio,
                 cfignm=path_fig+'sea-ice_SP_sept_'+COMP2D, cpal='ice', cbunit='frac.',
                 ctitle=ctit_clim,
                 lkcont=True, cfig_type=fig_type, 
                 lforce_lim=True)











# March:
# ~~~~~~

if COMP2D == 'CLIM':
    ctit_clim = 'Sea-Ice, March, obs.'
else:
    ctit_clim = 'Sea-Ice, March, '+COMP2D



# Nordic Seas:
#
bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, xnemo03[:,:],
              cfignm=path_fig+'sea-ice_NP_march_'+CONFRUN, cpal='ice', cbunit='frac.',
              ctitle='Sea-Ice, March, '+CONFRUN+' ('+cy1+'-'+cy2+')',
              lkcont=True, cfig_type=fig_type, 
              lforce_lim=True)

bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, xclim03[:,:]/ratio,
                cfignm=path_fig+'sea-ice_NP_march_'+COMP2D, cpal='ice', cbunit='frac.',
                ctitle=ctit_clim,
                lkcont=True, cfig_type=fig_type, 
                lforce_lim=True)
#
#
#
# Big south:

bp.plot("nproj")('spstere', 0., 1., 0.1, xlon0, xlat0, xnemo030[:,:],
              cfignm=path_fig+'sea-ice_SP_march_'+CONFRUN, cpal='ice', cbunit='frac.',
              ctitle='Sea-Ice, March, '+CONFRUN+' ('+cy1+'-'+cy2+')',
              lkcont=True, cfig_type=fig_type, 
              lforce_lim=True)

bp.plot("nproj")('spstere', 0., 1., 0.1, xlon0, xlat0, xclim030[:,:]/ratio,
                cfignm=path_fig+'sea-ice_SP_march_'+COMP2D, cpal='ice', cbunit='frac.',
                ctitle=ctit_clim,
                lkcont=True, cfig_type=fig_type, 
                lforce_lim=True)



print '\n Bye!'

