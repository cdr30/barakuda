#!/usr/bin/env python
#
# L. Brodeau, July 2011
#

import sys
import os
#import math
import numpy as nmp
from netCDF4 import Dataset

import barakuda_orca as bo
import barakuda_plot as bp
import barakuda_tool as bt









ORCA = os.getenv('ORCA')
if ORCA == None: print 'The ORCA environement variable is no set'; sys.exit(0)

RUN = os.getenv('RUN')
if RUN == None: print 'The RUN environement variable is no set'; sys.exit(0)

DIAG_D = os.getenv('DIAG_D')
if DIAG_D == None: print 'The DIAG_D environement variable is no set'; sys.exit(0)


CONFRUN = ORCA+'-'+RUN





#if len(sys.argv) != 4:
#    print 'Usage: '+sys.argv[0]+' <YYYY1> <YYYY2> <Nb. lat. points>'
#    sys.exit(0)

#cy1  = sys.argv[1] ; cy2 = sys.argv[2] ; jy1=int(cy1); jy2=int(cy2)
#nj   = int(sys.argv[3])


path_fig=DIAG_D+'/'

fig_type='png'

cf_in = DIAG_D+'/merid_transport_T_S_'+CONFRUN+'.nc'


if not os.path.exists(cf_in):
    print ' ERROR: plot_hovm_merid_trsp.py => old ascii file system not supported anymore!'
    print '        => need file '+cf_in+' !!!' ; sys.exit(0)




nbasins = len(bo.voce2treat)

id_in = Dataset(cf_in)

vyear = id_in.variables['time'][:]   ; Nby = len(vyear) ; ittic = bt.iaxe_tick(Nby)
vyear = vyear - 0.5
vlat  = id_in.variables['lat'][:]    ; Nlat = len(vlat)

for jb in range(nbasins):

    cbasin = bo.voce2treat[jb] ; # long name of basin
    cbas   = cbasin[:3] ;           # name as in cf_in ...
    
    if jb == 0:
        Xheat = nmp.zeros(nbasins*Nby*Nlat) ; Xheat.shape = [ nbasins, Nby, Nlat ]
        Xsalt = nmp.zeros(nbasins*Nby*Nlat) ; Xsalt.shape = [ nbasins, Nby, Nlat ]

    Xheat[jb,:,:] = id_in.variables['zomht_'+cbas][:,:]
    Xsalt[jb,:,:] = id_in.variables['zomst_'+cbas][:,:]
    print ' *** zomht_'+cbas+' and zomst_'+cbas+' sucessfully read into '+cf_in

id_in.close()

print ''




imask  = nmp.zeros(Nlat*Nby); imask.shape = [ Nby, Nlat ]

for jb in range(nbasins):

    if jb == 0:
        vyear = nmp.trunc(vyear) + 0.5 ; # in case 1990 and not 1990.5 !!!
        yr1=float(int(min(vyear)))
        yr2=float(int(max(vyear)))


    cbasin = bo.voce2treat[jb]; print '\n *** Basin: '+cbasin

    imask[:,:] = 0
    Lfinite = nmp.isfinite(Xheat[jb,:,:]) ; idx_good = nmp.where(Lfinite)
    imask[idx_good] = 1
            
    [ rmin, rmax, rdf ] = bt.get_min_max_df(Xheat[jb,5:,:],40)

    bp.plot("vert_section")(vyear[:], vlat[:], nmp.flipud(nmp.rot90(Xheat[jb,:,:])), nmp.flipud(nmp.rot90(imask[:,:])),
                            rmin, rmax, rdf,
                            cpal='RdBu', xmin=yr1, xmax=yr2+1., dx=ittic, lkcont=False,
                            zmin = vlat[0], zmax = vlat[Nlat-1], l_zlog=False, 
                            cfignm=path_fig+'MHT_'+CONFRUN+'_'+cbasin, cbunit='PW', cxunit='',
                            czunit=r'Latitude ($^{\circ}$N)',
                            ctitle=CONFRUN+': Northward advective meridional heat transport, '+cbasin,
                            cfig_type=fig_type, lforce_lim=False, i_cb_subsamp=2, l_z_increase=True)



    # Salt transport

    [ rmin, rmax, rdf ] = bt.get_min_max_df(Xsalt[jb,5:,:],40)

    bp.plot("vert_section")(vyear[:], vlat[:], nmp.flipud(nmp.rot90(Xsalt[jb,:,:])), nmp.flipud(nmp.rot90(imask[:,:])),
                            rmin, rmax, rdf,
                            cpal='RdBu', xmin=yr1, xmax=yr2+1., dx=ittic, lkcont=False,
                            zmin = vlat[0], zmax = vlat[Nlat-1], l_zlog=False, 
                            cfignm=path_fig+'MST_'+CONFRUN+'_'+cbasin, cbunit=r'10$^3$ tons/s', cxunit='',
                            czunit=r'Latitude ($^{\circ}$N)',
                            ctitle=CONFRUN+': Northward advective meridional salt transport, '+cbasin,
                            cfig_type=fig_type, lforce_lim=False, i_cb_subsamp=2, l_z_increase=True)


