#!/usr/bin/env python
#
# L. Brodeau, Dec. 2015

import sys
import numpy as nmp
from netCDF4 import Dataset
import string

from os.path import basename

import barakuda_tool as bt
import barakuda_thermo as bthermo

cv_u='U10M'
cv_v='V10M'

cv_wm ='WIND_MOD_10M' ; # OUTPUT !
cv_ke ='EKE'

if len(sys.argv) < 2 and len(sys.argv) > 3:
    print 'Usage: '+sys.argv[0]+' <IN_FILE_U.nc> (Year)'
    sys.exit(0)
cf_u  = sys.argv[1]
cf_v  = string.replace(cf_u, cv_u, cv_v)
cf_wm = basename(string.replace(cf_u, cv_u, cv_wm))


if len(sys.argv) == 3:
    cvy = sys.argv[2]
    iyear = int(cvy)

#print iyear ; sys.exit(0)


# First need time length:

bt.chck4f(cf_u) ; f_u_in = Dataset(cf_u)
vlon     = f_u_in.variables['lon'][:]
cunt_lon = f_u_in.variables['lon'].units
print 'LONGITUDE: ', cunt_lon
# Extracting the longitude 1D array:
vlat     = f_u_in.variables['lat'][:]
cunt_lat = f_u_in.variables['lat'].units
print 'LATITUDE: ', cunt_lat

# Extracting time 1D array:
vtime     = f_u_in.variables['time'][:] ; cunt_time = f_u_in.variables['time'].units
print 'TIME: ', cunt_time, '\n'
f_u_in.close()




Nt = len(vtime)

print 'Nt = ', Nt

#ii = Nt/31

for jt in range(Nt):


    #vtime[jt] = float(iyear) + 365./float(Nt) * float(jt)+0.5


    print ' *** jt = ', jt
    
    
    #  U10M
    #  ~~~~
    if jt == 0:
        bt.chck4f(cf_u)
        f_u_in = Dataset(cf_u)
        cunt_u = f_u_in.variables[cv_u].units
    xu     = f_u_in.variables[cv_u][jt,:,:]
    if jt == Nt-1: f_u_in.close()
    
    
    
    # V10M
    # ~~~
    if jt == 0:
        bt.chck4f(cf_v)
        f_v_in = Dataset(cf_v)
        cunt_v = f_v_in.variables[cv_v].units
    xv     = f_v_in.variables[cv_v][jt,:,:]
    if jt == Nt-1: f_v_in.close()
    
    

    
    
    # Checking dimensions
    # ~~~~~~~~~~~~~~~~~~~
    if jt == 0:
        dim_u = xu.shape ; dim_v = xv.shape
        if dim_u != dim_v:
            print 'Shape problem!!!'; print dim_u , dim_v
        print '\n'
        [ nj, ni ] = dim_u
        print 'ni, nj, nt = ', ni, nj, Nt
        xwm = nmp.zeros(nj*ni) ; xwm.shape = dim_u
        xke = nmp.zeros(nj*ni) ; xke.shape = dim_u
    
    
    # Building wm
    # ~~~~~~~~~~~
    xke = xu*xu + xv*xv
    xwm = nmp.sqrt(xke)
    xke = 0.5*xke
    

    
    # Creating output file
    # ~~~~~~~~~~~~~~~~~~~~
    if jt == 0:
        f_out = Dataset(cf_wm, 'w', format='NETCDF3_CLASSIC')
    
        # Dimensions:
        f_out.createDimension('lon', ni)
        f_out.createDimension('lat', nj)
        f_out.createDimension('time', None)
    
        # Variables
        id_lon = f_out.createVariable('lon','f4',('lon',))
        id_lat = f_out.createVariable('lat','f4',('lat',))
        id_tim = f_out.createVariable('time','f4',('time',))
        id_wm  = f_out.createVariable(cv_wm,'f4',('time','lat','lon',))
        id_ke  = f_out.createVariable(cv_ke,'f4',('time','lat','lon',))
    
        # Attributes
        id_tim.units = cunt_time
    
        #id_lat.long_name     = clnm_lat
        id_lat.units         = cunt_lat
        #id_lat.standard_name = csnm_lat
    
        #id_lon.long_name     = clnm_lon
        id_lon.units         = cunt_lon
        #id_lon.standard_name = csnm_lon
    
        id_tim.units         = cunt_time
    
        #id_wm.long_name = 'Surface specific humidity at 2m, built from U10M, V10M and MSL'
        #id_wm.units = 'kg/kg'
        #id_wm.code  = '133'
        #id_wm.table = '128'
    
        f_out.About = 'Created by L. Brodeau using MSL, U10M and V10M corresponding fields'
    
        # Filling variables:
        id_lat[:] = vlat[:]
        id_lon[:] = vlon[:]
        
    id_tim[jt]     = vtime[jt]
    id_wm[jt,:,:]  = xwm[:,:]
    id_ke[jt,:,:]  = xke[:,:] 
    
    if jt == Nt-1: f_out.close()
    
    
    
    
    
print ' *** file '+cf_wm+' created !'    
        
print 'Bye!'
