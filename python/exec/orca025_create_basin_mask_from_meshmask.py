#!/usr/bin/env python

# Petteri Uotila, 2015

import sys
import numpy as nmp
from netCDF4 import Dataset
from string import replace

if len(sys.argv) != 2:
    print 'Usage: '+sys.argv[0]+' <mesh_mask_ORCA025_file.nc>'
    sys.exit(0)

cf_mm  = sys.argv[1]
cf_out = replace(cf_mm, 'mesh_mask', 'basin_mask')


print '\n'

# Opening the Netcdf file:
f_mm = Dataset(cf_mm)
print 'File ', cf_mm, 'is open...\n'

# Extracting the longitude 2D array:
xlon = f_mm.variables['nav_lon'][:,:]

# Extracting the longitude 2D array:
xlat = f_mm.variables['nav_lat'][:,:]

# Extracting tmask at surface level:
tmask  = f_mm.variables['tmask'][0,0,:,:]

f_mm.close()



# Info on the shape of t:
[ nj, ni ] = tmask.shape

print 'Dimension = ', ni, nj, '\n'

mask_atl = nmp.zeros((nj,ni))
mask_pac = nmp.zeros((nj,ni))
mask_ind = nmp.zeros((nj,ni))
mask_soc = nmp.zeros((nj,ni))

f_mo025 = Dataset('new_maskglo_ORCA025.nc')
# ATL for ORCA025
# ~~~~~~~~~~~~~
mask_atl[:,:] = f_mo025.variables['tmaskatl'][:]
mask_atl[835:,:] = tmask[835:,:]

# PAC for ORCA1
# ~~~~~~~~~~~~~
mask_pac[:,:] = f_mo025.variables['tmaskpac'][:]
mask_pac[:357,:] = 0
mask_pac[355:376,225:266] = 0 # mask the Great Australian Bight

# IND for ORCA1
# ~~~~~~~~~~~~~
mask_ind[:,:] = f_mo025.variables['tmaskind'][:]

# Indo-Pacific
# ~~~~~~~~~~~~
idx = nmp.where((mask_pac==1)|(mask_ind==1))
mask_inp[idx] = 1

# Southern Ocean
mask_soc[:,:] = f_mo025.variables['tmaskant'][:]
mask_soc[255:357,:] = tmask[255:357,:]

f_mo025.close()

# Creating output file:
f_out = Dataset(cf_out, 'w',format='NETCDF3_CLASSIC')

# Dimensions:
f_out.createDimension('x', ni)
f_out.createDimension('y', nj)

# Variables
id_lon  = f_out.createVariable('nav_lon','f4',('y','x',))
id_lat  = f_out.createVariable('nav_lat','f4',('y','x',))

id_atl  = f_out.createVariable('tmaskatl' ,'f4',('y','x',)) ; id_atl.long_name = 'Atlantic Basin'
id_pac  = f_out.createVariable('tmaskpac' ,'f4',('y','x',)) ; id_pac.long_name = 'Pacific Basin'
id_ind  = f_out.createVariable('tmaskind' ,'f4',('y','x',)) ; id_ind.long_name = 'Indian Basin'
id_soc  = f_out.createVariable('tmasksoc' ,'f4',('y','x',)) ; id_soc.long_name = 'Southern Basin'
id_inp  = f_out.createVariable('tmaskinp' ,'f4',('y','x',)) ;  id_inp.long_name = 'Indo-Pacific Basin'


# Filling variables:
id_lat[:,:]    =   xlat[:,:]
id_lon[:,:]    =   xlon[:,:]

id_atl[:,:]  =  mask_atl[:,:]
id_pac[:,:]  =  mask_pac[:,:]
id_ind[:,:]  =  mask_ind[:,:]
id_soc[:,:]  =  mask_soc[:,:]

id_inp[:,:]  =  mask_inp[:,:]


f_out.About  = 'ORCA1 main oceanic basin land-sea mask created from '+cf_mm
f_out.Author = ' Generated with "orca025_create_basin_mask_from_meshmask.py" of BaraKuda (https://github.com/brodeau/barakuda)'

f_out.close()


print cf_out+' sucessfully created!'

