import math
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import sys
# setup stereographic basemap.
# lat_ts is latitude of true scale.
# lon_0,lat_0 is central point.

#m = Basemap(width=12000000,height=8000000,
#            resolution='l',projection='stere',\
#            lat_ts=50,lat_0=50,lon_0=-107.)


m = Basemap(width=12000000,height=8000000,
            resolution='l',projection='stere',\
            lat_ts=-50,lat_0=-90,lon_0=-107.)



m.drawcoastlines()


ni = 360
nj = 180

vlon = np.arange(0,ni)  + 0.5
vlat = np.arange(-90,nj/2) + 0.5


xlon = np.zeros(nj*ni) ; xlon.shape = [nj, ni]
xlat = np.zeros(nj*ni) ; xlat.shape = [nj, ni]
XF = np.zeros(nj*ni) ; XF.shape = [nj, ni]


for ji in range(ni): xlat[:,ji] = vlat[:]
for jj in range(nj):
    xlon[jj,:] = vlon[:]
    XF[jj,:] = math.cos(vlat[jj]/15.)






x0,y0 = m(xlon,xlat)


print vlon[:]
print ''
print vlat[:]




m.fillcontinents(color='coral',lake_color='aqua')


# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='aqua')





m.contourf(x0, y0, XF, 10)


plt.title("Stereographic Projection")

plt.show()


sys.exit(0)


# draw tissot's indicatrix to show distortion.
ax = plt.gca()
for y in np.linspace(m.ymax/20,19*m.ymax/20,9):
    for x in np.linspace(m.xmax/20,19*m.xmax/20,12):
        lon, lat = m(x,y,inverse=True)
        poly = m.tissot(lon,lat,1.5,100,\
                        facecolor='green',zorder=10,alpha=0.5)

