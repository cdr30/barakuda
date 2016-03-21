#!/usr/bin/env python

# Misc :
import os
import sys
from netCDF4 import Dataset
import numpy as nmp

from string import replace

# Mine:
import barakuda_tool as bt
import barakuda_plot as bp


#FIF='eps'
#FIF='svg'
FIF='png'


na = len(sys.argv)

if len(sys.argv) < 3:
    print 'Usage : '+sys.argv[0]+' <NEMO_FILE.nc4> <2D_variable_name>'
    sys.exit(0)

cf_in = sys.argv[1]
cv_in = sys.argv[2]

if cv_in == 'siconc' or 'ileadfra':
    rmax = 1. ; rmin = 0. ; dr = 0.1
    colmap = 'ice' ; cunit = 'frac.'

elif cv_in == 'sithic' or 'iicethic':
    rmax = 6. ; rmin = 0. ; dr = 0.25
    colmap = 'jet' ; cunit = 'm'

else:
    print 'ERROR!!! variable '+cv_in+' is unknown!!!' ; sys.exit(0)


cfig_suff = replace(os.path.basename(cf_in), '.nc4', '')


bt.chck4f(cf_in)
id_in = Dataset(cf_in)
xlon = id_in.variables['nav_lon'][:,:]
xlat = id_in.variables['nav_lat'][:,:]
XF2D = id_in.variables[cv_in][:,:,:]
id_in.close()


[Nt, nj, ni] = nmp.shape(XF2D)

print 'Shape of Ice =', [Nt, nj, ni]


#Nt = 24

jm = 0
for jt in range(Nt):

    jy = jt/12 + 1
    if jt % 12 == 0: jm = 0    
    jm = jm + 1


    #print 'jt, jm, jy =', jt, jm, jy

    ct = str(jt).zfill(4)
    cm = str(jm).zfill(2)
    cy = str(jy).zfill(4)



    cfig = cv_in+'_'+cfig_suff+'_'+ct

    print '  *** will create fig '+cfig+' (year = '+cy+', month = '+cm+')'
    
    bp.plot("nproj")('spstere', rmin, rmax, dr, xlon, xlat, XF2D[jt,:,:],
                     cfignm=cfig, cpal=colmap, cbunit=cunit,
                     ctitle='Sea-Ice, year = '+cy+', month = '+cm,
                     lkcont=True, cfig_type=FIF,
                     lforce_lim=True)



cfig_out = cv_in+'_'+cfig_suff+'.gif'

os.system("convert -delay 100 -loop 0 "+cv_in+"*.png "+cfig_out+" > out_conv.out")
