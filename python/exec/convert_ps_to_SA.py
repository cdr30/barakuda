#!/usr/bin/env python


# L. Brodeau, 2015


# Practical salinity to absolute salinity (TEOS 10)

import sys
import os
import numpy as nmp
from netCDF4 import Dataset
from string import replace

import gsw

#l_accurate = True
l_accurate = False

SSO = 35.16504


if len(sys.argv) < 3:
    print 'Usage: '+sys.argv[0]+' <Salinity_file_to_convert> <salinity_name> (2d)'
    sys.exit(0)


cf_sal  = sys.argv[1]
cv_sal  = sys.argv[2]

l2d = False

if len(sys.argv) == 4:
    cv_2d  = sys.argv[3]
    if cv_2d != '2d': print 'Usage: '+sys.argv[0]+' <Salinity_file_to_convert> <salinity_name> (2d)'
    l2d=True


cf_out = replace(cf_sal, cf_sal, 'absolute_salinity_'+cf_sal)


os.system('rm -f '+cf_out)
os.system('cp '+cf_sal+' '+cf_out)







print '\n'

# Opening the Netcdf file:
f_out = Dataset(cf_out, 'r+')     # r+ => can read and write in the file... )
print 'File ', cf_out, 'is open...\n'

# Extracting tmask at surface level:
if l2d:
    xsal  = f_out.variables[cv_sal][:,:,:]
else:
    xsal  = f_out.variables[cv_sal][:,:,:,:]

if l_accurate and not l2d:
    
    vz    = f_out.variables['deptht'][:]
    
    [nt,nk,nj,ni] = nmp.shape(xsal)

    xdepth = nmp.zeros((nk,nj,ni))

    # building 3d arrays of depth, lon and lat:
    for jk in range(nk): xdepth[jk,:,:] = vz[jk]

    # pressure should be in dbar and it's the same as the depth in metre actually:
    for jt in range(nt):
        print ' jt =', jt
        f_out.variables[cv_sal][jt,:,:,:] = gsw.SA_from_SP(xsal[jt,:,:,:], xdepth, -140., 0.)

else:
    # Fabien says it's enough:
    if l2d:
        f_out.variables[cv_sal][:,:,:]   = xsal[:,:,:]*SSO/35.
    else:
        f_out.variables[cv_sal][:,:,:,:] = xsal[:,:,:,:]*SSO/35. 


f_out.variables[cv_sal].long_name = 'Absolute Salinity (TEOS10) build from practical salinity (*35.16504/35)'

f_out.close()




print cf_out+' sucessfully created!'

