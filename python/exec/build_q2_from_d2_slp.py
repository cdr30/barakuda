#!/usr/bin/env python
#
# L. Brodeau, Feb.2001

import sys
import numpy as nmp
from netCDF4 import Dataset
import string

from os.path import basename

import barakuda_tool as bt
import barakuda_thermo as bthermo

cv_d2='D2M'
cv_p0='MSL'
cv_q2='Q2M' ; # OUTPUT !


if len(sys.argv) != 2:
    print 'Usage: '+sys.argv[0]+' <IN_FILE_D2.nc>'
    sys.exit(0)

cf_d2   = sys.argv[1]
cf_p0 = string.replace(cf_d2, cv_d2, cv_p0)
cf_q2 = basename(string.replace(cf_d2, cv_d2, cv_q2))



# First need time length:

bt.chck4f(cf_d2) ; f_d2_in = Dataset(cf_d2)
vlon     = f_d2_in.variables['lon'][:]
cunt_lon = f_d2_in.variables['lon'].units
print 'LONGITUDE: ', cunt_lon
# Extracting the longitude 1D array:
vlat     = f_d2_in.variables['lat'][:]
cunt_lat = f_d2_in.variables['lat'].units
print 'LATITUDE: ', cunt_lat

# Extracting time 1D array:
vtime     = f_d2_in.variables['time'][:] ; cunt_time = f_d2_in.variables['time'].units
print 'TIME: ', cunt_time, '\n'
f_d2_in.close()





Nt = len(vtime)


print 'Nt = ', Nt

for jt in range(Nt):

    print ' *** jt = ', jt
    
        
    # D2M
    # ~~~
    if jt == 0:
        bt.chck4f(cf_d2)
        f_d2_in = Dataset(cf_d2)
        cunt_d2 = f_d2_in.variables[cv_d2].units
    xd2     = f_d2_in.variables[cv_d2][jt,:,:]
    if jt == Nt-1: f_d2_in.close()
    
    
    # MSL
    # ~~~
    if jt == 0:
        bt.chck4f(cf_p0)
        f_p0_in = Dataset(cf_p0)
        cunt_p0 = f_p0_in.variables[cv_p0].units
    xp0     = f_p0_in.variables[cv_p0][jt,:,:]
    if jt == Nt-1: f_p0_in.close()
    
    
    
    
    
    # Checking dimensions
    # ~~~~~~~~~~~~~~~~~~~
    if jt == 0:
        dim_d2 = xd2.shape ; dim_p0 = xp0.shape
        if dim_d2 != dim_p0:
            print 'Shape problem!!!'; print dim_d2 , dim_p0
        print '\n'
        [ nj, ni ] = dim_d2
        print 'ni, nj, nt = ', ni, nj, Nt
        xq2 = nmp.zeros(nj*ni) ; xq2.shape = dim_d2
    
    
    # Building q2
    # ~~~~~~~~~~~
    xq2 = bthermo.qa_e_p(bthermo.e_sat(xd2), xp0)
    
    
    
    
    
    # Creating output file
    # ~~~~~~~~~~~~~~~~~~~~
    if jt == 0:
        f_out = Dataset(cf_q2, 'w', format='NETCDF3_CLASSIC')
    
        # Dimensions:
        f_out.createDimension('lon', ni)
        f_out.createDimension('lat', nj)
        f_out.createDimension('time', None)
    
        # Variables
        id_lon = f_out.createVariable('lon','f4',('lon',))
        id_lat = f_out.createVariable('lat','f4',('lat',))
        id_tim = f_out.createVariable('time','f4',('time',))
        id_q2  = f_out.createVariable(cv_q2,'f4',('time','lat','lon',))
    
        # Attributes
        id_tim.units = cunt_time
    
        #id_lat.long_name     = clnm_lat
        id_lat.units         = cunt_lat
        #id_lat.standard_name = csnm_lat
    
        #id_lon.long_name     = clnm_lon
        id_lon.units         = cunt_lon
        #id_lon.standard_name = csnm_lon
    
        id_tim.units         = cunt_time
    
        id_q2.long_name = 'Surface specific humidity at 2m, built from D2M and MSL'
        id_q2.units = 'kg/kg'
        id_q2.code  = '133'
        id_q2.table = '128'
    
        f_out.About = 'Created by L. Brodeau using MSL and D2M corresponding fields'
    
        # Filling variables:
        id_lat[:] = vlat[:]
        id_lon[:] = vlon[:]
        
    id_tim[jt]     = vtime[jt]
    id_q2[jt,:,:]  = xq2[:,:] 
    
    if jt == Nt-1: f_out.close()
    
    
    
    
    
    
        
print 'Bye!'
