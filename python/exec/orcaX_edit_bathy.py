#!/usr/bin/env python

import sys
import numpy as nmp
from os      import system
from netCDF4 import Dataset
from string  import replace

if len(sys.argv) != 2:
    print 'Usage: '+sys.argv[0]+' <nemo_bathy.nc>'
    sys.exit(0)

cf_old = sys.argv[1]
cf_new = replace(cf_old, '.nc', '_new.nc')

cv_bathy = 'Bathymetry'


# First, creating a copy:
system('rm -f '+cf_new)
system('cp '+cf_old+' '+cf_new)





print '\n'

# Opening the Netcdf file:
f_new = Dataset(cf_new,  'r+')
print 'File ', cf_new, 'is open...\n'

Xbathy = f_new.variables[cv_bathy][:,:]


# Edit zone:

# Removing caspian sea in ORCA1:
Xbathy[200:238,332:346] = 0.

f_new.variables[cv_bathy][:,:] = Xbathy[:,:]


f_new.close()






print cf_new+' sucessfully created!'

