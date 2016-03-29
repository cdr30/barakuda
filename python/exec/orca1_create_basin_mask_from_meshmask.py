#!/usr/bin/env python

# L. Brodeau, april 2010

import sys
import numpy as nmp
from netCDF4 import Dataset
from string import replace

if len(sys.argv) != 2:
    print 'Usage: '+sys.argv[0]+' <mesh_mask_ORCA1_file.nc>'
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

mask_ip1 = nmp.zeros((nj,ni))
mask_inp = nmp.zeros((nj,ni))




# ATL for ORCA1
# ~~~~~~~~~~~~~
mask_atl[:,:] = tmask[:,:]

# Removing Southern Ocean:
mask_atl[:95,:] = 0

# Removing Pacific and Indian
mask_atl[0:246,0:190] = 0 # 246 => to keep Pacific side of the arctic basin...
mask_atl[0:168,0:223] = 0 ; mask_atl[0:255,310:] = 0
mask_atl[165:177,190:204] = 0; mask_atl[165:180,190:198] = 0; mask_atl[165:170,200:206] = 0
mask_atl[188:209,282:] = 0; mask_atl[209:215,288:] = 0



# REMOVING INDONESIA + AUSTRALIA
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mask_ip1[:,:] = tmask[:,:]

mask_ip1[114:122,53:75] = 0
mask_ip1[119:126,68:74] = 0
mask_ip1[124:143,44:59] = 0
mask_ip1[128:159,33:42] = 0
mask_ip1[120:142,52:61] = 0
mask_ip1[124:136,41:70] = 0
mask_ip1[127:128,37:42] = 0
mask_ip1[120:126,60:70] = 0
mask_ip1[141:158,30:33] = 0
mask_ip1[152:162,26:30] = 0






# PAC for ORCA1
# ~~~~~~~~~~~~~
mask_pac[:,:] = tmask[:,:]

# Removing Southern Ocean until souther Australia:
mask_pac[:95,:] = 0


# Removing Indonesian side
mask_pac[:,:45] = 0
mask_pac[88:145,45:61] = 0
mask_pac[112:125,59:70] = 0
mask_pac[123:136,60:67] = 0
mask_pac[88:99,60:71] = 0 # bottom Australia


# V2
#mask_pac[:,:26] = 0



# Removing Atlantic
idxatl = nmp.where(mask_atl == 1.0)
mask_pac[idxatl] = 0

# Removing atlantic bottom and the rest (Indian)
mask_pac[83:,224:] = 0


# IND for ORCA1
# ~~~~~~~~~~~~~
mask_ind[:,:] = tmask[:,:]

# Removing Southern Ocean until southern Australia:
mask_ind[:95,:] = 0

# Removing Atl and Pac
mask_ind[:,:] = mask_ind[:,:] - mask_atl[:,:] - mask_pac[:,:]

mask_ind[93:100,46:68] = 0 # australia bottom

# Removing Mediterranean+Caspian sea:
mask_ind[192:228,279:329] = 0
mask_ind[198:242,328:344] = 0







# Indo-Pacific
# ~~~~~~~~~~~~
mask_inp[:,:] = tmask[:,:]
mask_inp[:95,:] = 0

# Removing Atlantic
idxatl = nmp.where(mask_atl == 1.0)
mask_inp[idxatl] = 0

mask_inp[93:100,46:68] = 0 # australia bottom

# Removing Mediterranean sea:
mask_inp[192:228,279:329] = 0
mask_inp[198:242,328:344] = 0

# Removing indonesia
#mask_inp[:,:] = mask_inp[:,:] * mask_ip1[:,:]





# Souther Ocean

mask_soc[:,:] = tmask[:,:]

idxatl = nmp.where(mask_atl+mask_pac+mask_ind > 0.5)
mask_soc[idxatl] = 0
mask_soc[122:,:] = 0






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
f_out.Author = ' Generated with "orca1_create_basin_mask_from_meshmask.py" of BaraKuda (https://github.com/brodeau/barakuda)'

f_out.close()


print cf_out+' sucessfully created!'

