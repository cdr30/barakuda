#!/home/x_laubr/bin/python

# L. Brodeau, April 2011

import sys
import numpy as nmp
import string
import os
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot  as bp

# Some constants:
Pi = 3.141592654 ; to_rad = Pi/180. ; RE = 6.36 ; # (10^6 m)


NN_SST = os.getenv('NN_SST')
if NN_SST == None: print 'The NN_SST environement variable is no set'; sys.exit(0)


if len(sys.argv) != 2:
    print 'Usage: '+sys.argv[0]+' <nino_file.nc>'
    sys.exit(0)

cf_in  = sys.argv[1]


cname = string.replace(os.path.basename(cf_in), '.nc', '')

print "\n"

fig_type='png'


bt.chck4f(cf_in)
id_in = Dataset(cf_in)
vtime = id_in.variables['time'][:] ; nbm = len(vtime)
vsst  = id_in.variables[NN_SST][:]
id_in.close()


nt = len(vsst)
print ' => '+str(nt)+' months  => '+str(nt/12)+' years\n'


ittic = bt.iaxe_tick(nt/12)


# Array to contain nino series:
xnino = nmp.zeros(nt*4) ; xnino.shape = [ nt, 4 ]


for jt in nmp.arange(nt): xnino[jt,0] = vsst[jt]


# 5-month running mean:

for jt in nmp.arange(2,nt-2):
    xnino[jt,1] = (xnino[jt-2,0] + xnino[jt-1,0] + xnino[jt,0] + xnino[jt+1,0] + xnino[jt+2,0]) / 5.

xnino[0:2,1] = xnino[2,1] ; xnino[nt-2:nt,1] = xnino[nt-3,1]

print '\n'

print 'mean value for sst mean = ', nmp.sum(xnino[:,0])/nt
print 'mean value for sst 5-m-r mean = ', nmp.sum(xnino[:,1])/nt

# least-square curve for 5-month running mean:
sumx  = nmp.sum(vtime[:]) ; sumy  = nmp.sum(xnino[:,1])
sumxx = nmp.sum(vtime[:]*vtime[:])
sumxy = nmp.sum(vtime[:]*xnino[:,1])
a = ( sumx*sumy - nt*sumxy ) / ( sumx*sumx - nt*sumxx )
b = ( sumy - a*sumx ) / nt
print 'a, b =', a, b

# least-square linear trend:
xnino[:,2] = a*vtime[:] + b
print 'mean value for least-square linear trend = ', nmp.sum(xnino[:,2])/nt

# anomaly
xnino[:,3] = xnino[:,1] - xnino[:,2] ; # anomaly for 5-month running mean
print 'mean value for anomaly = ', nmp.sum(xnino[:,3])/nt



# save serie into a file:
#cf_out = cname+'_anomaly.nc'
#f = open(cf_out, 'w')
#f.write('# created with '+sys.argv[0]+' from file '+cf_in+'\n')
#f.write('# time, mean SST,  5-month RM,  linear trend,    anomaly     \n')
#for jt in range(nt):
#    f.write(str(round(vtime[jt],3))+'   '+str(round(xnino[jt,0],5))+'   '\
#            +str(round(xnino[jt,1],5))+'   '+str(round(xnino[jt,2],5))+'   '+str(round(xnino[jt,3],5))+'\n')
#print '\n'
#print ' ascii file '+cf_out+' created!\n'



bp.plot_enso( vtime, xnino[:,0], cfignm=cname, dt_year=ittic )


