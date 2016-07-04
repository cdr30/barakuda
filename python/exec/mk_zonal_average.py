#!/usr/bin/env python

# L. Brodeau, June 2012

import sys
from os.path import splitext
import numpy as nmp
from netCDF4 import Dataset
from string import replace as rplc

import barakuda_tool as bt

rmv = -9999.

if len(sys.argv) != 3:
    print 'Usage: '+sys.argv[0]+' <FILE_lat-lon.nc> <variable>'
    sys.exit(0)

cf_in  = sys.argv[1]
cv_in  = sys.argv[2]

cn_file, cn_ext = splitext(cf_in)

print cn_ext

cn_file = rplc(cn_file, cv_in+'_', '')
print cn_file

print "\n"



bt.chck4f(cf_in) ; f_in = Dataset(cf_in)

# Extracting the longitude and 1D array:
vlon     = f_in.variables['lon'][:]
cunt_lon = f_in.variables['lon'].units
print 'LONGITUDE: ', cunt_lon

# Extracting the longitude 1D array:
vlat     = f_in.variables['lat'][:]
cunt_lat = f_in.variables['lat'].units
print 'LATITUDE: ', cunt_lat

# Extracting time 1D array:
vtime     = f_in.variables['time'][:] ; cunt_time = f_in.variables['time'].units
print 'TIME: ', cunt_time, '\n'

# Field !!!
#rmv    = f_in.variables[cv_in]._FillValue
xfield = f_in.variables[cv_in][:,:,:]

#print 'Missing value for '+cv_in+' is : ', rmv, '\n'

f_in.close()


nt = len(vtime)
print ' nt = '+str(nt)




# Checking dimensions
# ~~~~~~~~~~~~~~~~~~~
[ nt, nj, ni ] = xfield.shape
print ' DIMENSION =>  ni, nj, nt = ', ni, nj, nt


VZ = nmp.zeros((nt,nj))

for jt in range(nt):

    cjt  = '%3.3d' %(jt+1)

    # Zonally-averaging:
    for jj in range(nj):
    
        cpt = 0
    
        for ji in range(ni):
            val = xfield[jt,jj,ji]
            if val != rmv:
                cpt = cpt + 1
                VZ[jt,jj] = VZ[jt,jj] + val

        if cpt >= 1 : VZ[jt,jj] = VZ[jt,jj]/cpt




print '\n Creating ascii file!'

cf_out = 'zonal_'+cv_in+'_'+cn_file+'.dat'

f = open(cf_out, 'w')

f.write('# created with '+sys.argv[0]+' from file '+cf_in+'\n')


for jj in range(nj):

    mean_val = nmp.mean(VZ[:,jj])
    

    # Writing latitudes:
    f.write(str(vlat[jj]))

    # time-averaged column first if relevant
    if nt > 1:
        f.write('   '+str(mean_val))

    # snapshot column
    for jt in range(nt):
        f.write('   '+str(VZ[jt,jj]) )

    f.write('\n')            


print ' ASCII file '+cf_out+' created!\n' ; print '\n'
        



print '\n Creating netcdf file!'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cf_out = 'zonal_'+cv_in+'_'+cn_file+'.nc'

f_out = Dataset(cf_out, 'w',format='NETCDF3_CLASSIC')

# Dimensions:
f_out.createDimension('lat', nj)
f_out.createDimension('time', None)

# Variables
id_lat = f_out.createVariable('lat','f4',('lat',))
id_tim = f_out.createVariable('time','f4',('time',))
id_f1  = f_out.createVariable(cv_in,'f4',('time','lat',))
id_f2  = f_out.createVariable(cv_in+'_anom','f4',('time','lat',))

# Attributes
#id_tim.units = cunt_time
#
#id_lat.long_name     = clnm_lat
#id_lat.units         = cunt_lat
#id_lat.standard_name = csnm_lat
#
#id_tim.units         = cunt_time
#
#
#id_f1.long_name = clnm_flx
#id_f1.units = cunit
#id_f1.code  = cvin_code
#id_f1.table = cvin_table
#
f_out.about = 'Diagnostics created with BaraKuda (https://github.com/brodeau/barakuda)'

# Filling variables:
id_lat[:] = vlat[:]


for jt in range(nt):
    id_tim[jt] = vtime[jt]
    id_f1[jt,:] = VZ[jt,:]
    id_f2[jt,:] = VZ[jt,:] - nmp.mean(VZ[:,:], axis=0)

f_out.close()
        
