#!/usr/bin/env python


# L. Brodeau, 2015


# Potential temperature to conservative temperature (TEOS 10)

import sys
import os
import numpy as nmp
from netCDF4 import Dataset
from string import replace

import gsw

#SSO = 35.16504


if len(sys.argv) != 5:
    print 'Usage: '+sys.argv[0]+' <Temperature_file_to_convert> <temperature_name> <Absolute_salinity_file> <salinity_name>'
    sys.exit(0)


cf_temp  = sys.argv[1]
cv_temp  = sys.argv[2]
cf_sal   = sys.argv[3]
cv_sal   = sys.argv[4]

cf_out = replace(cf_temp, cf_temp, 'conservative_temperature_'+cf_temp)

os.system('rm -f '+cf_out)
os.system('cp '+cf_temp+' '+cf_out)







print '\n'


f_sal = Dataset(cf_sal)     # r+ => can read and write in the file... )
xsal  = f_sal.variables[cv_sal][:,:,:,:]
f_sal.close()


print '\n'


# Opening the Netcdf file:
f_out = Dataset(cf_out, 'r+')     # r+ => can read and write in the file... )
print 'File ', cf_out, 'is open...\n'

# Extracting tmask at surface level:
xtemp  = f_out.variables[cv_temp][:,:,:,:]

#xtemp[:,:,:,:] = xtemp[:,:,:,:]*2.


#gsw.CT_from_pt(SA, pt)

f_out.variables[cv_temp][:,:,:,:] = gsw.CT_from_pt(xsal, xtemp)


f_out.variables[cv_temp].long_name = 'Conservative Temperature (TEOS10) built from potential temperature'

f_out.close()




print cf_out+' sucessfully created!'

