
import sys
import numpy as nmp

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from math import trunc

# Mine:
import barakuda_orca as brkdo


# Old savefig option:
# savefig(cfignm+'.'+cfig_type, dpi=100, facecolor='w', edgecolor='w', orientation='portrait'


# Time-series:
WDTH_TS     = 13.2
FIG_SIZE_TS = (WDTH_TS,4.2)
DPI_TS      = 120
AXES_TS     = [0.08, 0.05, 0.89, 0.89]




# For projections :
# =================
# lcc = Lambert conformal conic
#
#         ['natarct', 'lcc', -60., 40., 80., 72.,            55., -32., 10., 'l' ],   # NATL + Arctic

# my zone PROJ llcrnrlon llcrnrlat urcrnrlon urcrnrlat   lat1  lon0 mer/par continent-res
projection_def = [
         ['nseas',   'lcc',  -55., 40., 55., 75.,           60., -20., 10., 'l' ],   # Nordic seas
         ['natarct', 'lcc', -52., 38., 110., 75.,            60., -15., 10., 'l' ],   # NATL + Arctic
         ['labir',   'lcc',  -62., 48., -10., 75.,          50., -30.,  5., 'l' ],
         ['labsp',   'lcc',  -60., 48., 50., 75.5,          50., -30., 10., 'l' ],
         ['npol',    'stere', -75., 45., 100., 60.,          80., -30., 10., 'l' ],
         ['npol2',   'stere', -55., 40., 145., 40.,          80.,  -5., 10., 'l' ],   # North Pole
         ['spstere', 'stere',  0.,  0.,  0., 0.,            -48.,  90., 10., 'l' ],   # South Pole Default matplotlib!
         ['matl' ,   'cyl',  -82.,-21.,  12., 79.,          30., -30., 15., 'l' ],   # Nordic seas
         ['atmed',   'lcc',  -18., 33.,  -2., 42.,          30., -10.,  5., 'h' ],
         ['kav7' ,   'kav',    0.,  0.,   0.,  0.,           0.,   0.,  0., 'l' ] ] # global map-monde
         
#
#         ['spol2',   'stere', -45.,-35., 130., -35.,        -88.,   5., 10., 'l' ],   # South Pole
#         ['spol2','stere', -45.,-38., 130.,-35.,         -90.,  -5., 10., 'l' ],   # South Pole

#=================================================================================
#          ['natarct','lcc', -65., 35., 88., 70.,     50., -30., 10., 'h' ],   # NATL + Arctic










def plot_nproj_extra(czone, rmin, rmax, dc, xlon, xlat, XF, XI,
                     cfignm='fig', lkcont=False, cpal='jet', cbunit=' ',
                     cfig_type='pdf', ctitle=' ', lforce_lim=False,
                     cb_orient='vertical', i_cb_subsamp=1, dpi_fig=140,
                     xlon_o=[ 0. ], xlat_o=[ 0. ], XI_o=[ 0. ], rice_crit=0.15):

    # Plot projection with basemap...

    #===================================================================================
    # INPUT:
    #          xlon and xlat can be 1D or 2D !!!
    #
    #===================================================================================


    
    from mpl_toolkits.basemap import Basemap
    from mpl_toolkits.basemap import shiftgrid
    import barakuda_colmap

    font_ttl, font_ylb, font_clb = __font_unity__()

    # For projections :    
    vp = __give_proj__(czone) ; # projection information

    l_ice_obs = False
    if len(nmp.shape(XI_o)) > 1: l_ice_obs = True



    # must work with XFtmp rather than XF, because sometimes XF is overwrited...
    [ny, nx] = nmp.shape(XF)
    XFtmp = nmp.zeros(ny*nx) ; XFtmp.shape = [ny, nx]
    XFtmp[:,:] = XF[:,:]

    XItmp = nmp.zeros(ny*nx) ; XItmp.shape = [ny, nx]
    XItmp[:,:] = XI[:,:]
    
    if l_ice_obs:
        [ny2, nx2] = nmp.shape(XI_o)
        XIotmp = nmp.zeros(ny2*nx2) ; XIotmp.shape = [ny2, nx2]
        XIotmp[:,:] = XI_o[:,:]



    if len(nmp.shape(xlat)) == 1 and len(nmp.shape(xlon)) == 1:
        #if czone == 'kav7' and xlon[0] >= 0.:
        #    # Shifting data and longitude to be consistent with map projection
        #    XItmp, xlon = shiftgrid(180.+xlon[0], XItmp, xlon, start=False, cyclic=360.0)            
        #    XFtmp, xlon = shiftgrid(180.+xlon[0], XFtmp, xlon, start=False, cyclic=360.0)
        LON_2D, LAT_2D = nmp.meshgrid(xlon,xlat)
    else:
        LAT_2D = nmp.zeros(ny*nx) ; LAT_2D.shape = [ny, nx] ; LAT_2D[:,:] = xlat[:,:]
        LON_2D = nmp.zeros(ny*nx) ; LON_2D.shape = [ny, nx] ; LON_2D[:,:] = xlon[:,:]



    if lforce_lim: __force_min_and_max__(rmin, rmax, XFtmp)

    vc  = __vcontour__(rmin, rmax, dc)
    vci = __vcontour__(0., 1., 0.1)

    # Colorbar position/size if horizontal
    vcbar = [0.1, 0.08, 0.86, 0.03]

    # Figure/canvas size:
    if cb_orient == 'horizontal':
        if czone == 'natarct':
            vfig_size = [ 5.8, 5.6 ]; vsporg = [0.08, 0.1, 0.9,  0.92]
            vcbar = [0.05, 0.08, 0.9, 0.03]
        if czone == 'npol2':
            vfig_size = [ 4.4, 5.6 ];  vsporg = [0.01, 0.15, 1., 0.8]
            vcbar = [0.05, 0.065, 0.92, 0.03]
        if czone == 'kav7':
            vfig_size = [ 8.1, 5.6 ];  vsporg = [0.001, 0.15, 1., 0.8]
            vcbar = [0.04, 0.08, 0.92, 0.03]
            
    else:
        # Vertical color bar on the rhs
        vfig_size = [ 7., 7. ]; vsporg = [0.1, 0.1, 0.85, 0.85]
        if czone == 'nseas':   vfig_size = [ 7., 5.4 ]; vsporg = [0.085,  0.03, 0.9, 0.94]
        if czone == 'natarct': vfig_size = [ 7.7, 7. ]; vsporg = [0.07,  0.033, 0.92, 0.93]
        if czone == 'spstere': vfig_size = [ 7., 5.8 ]; vsporg = [0.075, 0.035, 0.93, 0.93]
        if czone == 'npol2':   vfig_size = [ 7., 7.1 ]; vsporg = [0.085, 0.03, 0.91, 0.94]
        #if czone == 'kav7':    vfig_size = [ 7., 5.  ]; vsporg = [0.085, 0.03, 0.91, 0.94]
        

        
    fig = plt.figure(num = 1, figsize=(vfig_size), dpi=None, facecolor='w', edgecolor='k')
    ax  = plt.axes(vsporg, axisbg = 'w')


    ## Palette:
    palette = barakuda_colmap.chose_palette(cpal)    
    pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)
    mpl.rcParams['contour.negative_linestyle'] = 'solid'; plt.contour.negative_linestyle='solid'

    palettei  = barakuda_colmap.chose_palette('blanc')
    #pal_normi = colors.Normalize(vmin = 0., vmax = 1., clip = True)




    


    if vp[1] == 'lcc' or vp[1] == 'cyl' :
        carte = Basemap(llcrnrlon=vp[2],llcrnrlat=vp[3],urcrnrlon=vp[4],urcrnrlat=vp[5],\
                        resolution=vp[9],area_thresh=1000.,projection=vp[1],\
                        lat_1=vp[6],lon_0=vp[7])

    elif vp[1] == 'stere' :
        if vp[0] == 'spstere' or vp[0] == 'npstere':
            carte = Basemap(projection=vp[0], boundinglat=vp[6], lon_0=vp[7], resolution=vp[9])
        else:
            carte = Basemap(llcrnrlon=vp[2],llcrnrlat=vp[3],urcrnrlon=vp[4],urcrnrlat=vp[5],\
                          resolution=vp[9],area_thresh=1000.,projection='stere',\
                          lat_0=vp[6],lon_0=vp[7])
    elif vp[1] == 'kav' :
        print ' *** plot_nproj.barakuda_plot => Projection '+vp[0]+' / '+str(vp[7])+' / '+vp[9]
        carte = Basemap(projection=vp[0],lon_0=vp[7],resolution=vp[9])
            
    else:
        print 'ERROR: barakuda_plot.py => proj type '+vp[1]+' unknown!!!'; sys.exit(0)

    x0,y0 = carte(LON_2D,LAT_2D)

    if l_ice_obs:
        x2,y2 = carte(xlon_o,xlat_o)

    

    cf = carte.contourf(x0, y0, XFtmp, vc, cmap = palette, norm = pal_norm)
    # Black contours if needed :
    if lkcont:
        ckf = carte.contour(x0, y0, XFtmp, vc, colors='k', linewidths=0.5)
        if cpal != 'ice':
            for c in cf.collections: c.set_zorder(0.5)   # Changing zorder so black cont. on top
        for c in ckf.collections: c.set_zorder(0.5) # of filled cont. and under continents (zorder 1)


    # Adding Sea-ice:
    cfi = carte.contourf(x0, y0, XItmp, [rice_crit, 1.], cmap = palettei)
    cf2 = carte.contour( x0, y0, XItmp, [rice_crit], colors='k', linewidths=0.5)
    for c in cfi.collections: c.set_zorder(0.75)
    for c in cf2.collections: c.set_zorder(1)


    if l_ice_obs:
        cf3 = carte.contour( x2, y2, XIotmp, [rice_crit], colors='r', linewidths=2.)
        for c in cf3.collections: c.set_zorder(1)


    carte.drawcoastlines() ; carte.fillcontinents(color='grey') ; carte.drawmapboundary()



    

    if vp[1] == 'lcc' or vp[1] == 'cyl':
        carte.drawmeridians(nmp.arange(-360,360,vp[8]), labels=[0,0,0,1])
        carte.drawparallels(nmp.arange( -90, 90,vp[8]), labels=[1,0,1,0])

    if vp[1] == 'stere':
        carte.drawmeridians(nmp.arange(-180,180,20), labels=[0,0,0,1])
        carte.drawparallels(nmp.arange( -90, 90,10), labels=[1,0,0,0])


    #plt.title(ctitle, **font_ttl)


    # ADDING COLORBAR
    # ===============

    if cb_orient == 'horizontal':
        clbax = fig.add_axes(vcbar) # axes for colorbar
        clb   = plt.colorbar(cf, cax=clbax, ticks=vc, drawedges=lkcont, orientation='horizontal')
        for t in clb.ax.get_xticklabels(): t.set_font.size(10)
    else:
        clb = plt.colorbar(cf, ticks=vc, drawedges=lkcont)
        for t in clb.ax.get_yticklabels(): t.set_fontsize(16)


    if i_cb_subsamp > 1:
        cn_clb = [] ; jcpt = 0
        for rtck in vc:
            if jcpt % i_cb_subsamp == 0:                
                if float(int(rtck)) == round(rtck,0):
                    cn_clb.append(str(int(rtck))) ; # we can drop the ".0"
                else:
                    cn_clb.append(str(rtck)) ; # keeping the decimals...
            else:
                cn_clb.append(' ')
            jcpt = jcpt + 1
        clb.ax.set_xticklabels(cn_clb)

    clb.set_label('('+cbunit+')', **font_clb)

    font_lab = { 'fontname':'Arial', 'fontweight':'normal', 'fontsize':24 }
    ax.annotate(ctitle, xy=(0.4, 0.95), xycoords='axes fraction', **font_lab)



    plt.savefig(cfignm+'.'+cfig_type, dpi=dpi_fig, orientation='portrait', transparent=True) ; #, transparent=True, acecolor='w', edgecolor='w',

    plt.close(1)

    del LON_2D, LAT_2D, XFtmp

    return










# LOCAL functions
# ===============




def __message__(ccr):

    # Find the CONF from CONF-RUN and exit if CONF does not exist!
    i = 0 ; conf = ''
    while i < len(ccr) and ccr[i] != '-' : conf = conf+ccr[i]; i=i+1
    print 'conf =', conf, '\n'
    return conf



def __get_mat__(cf):

    f1 = open(cf, 'r') # for reading
    lines1=f1.readlines()
    f1.close()

    zm   = []
    jy   = 0

    for l in lines1:
        if l[0] != '#':
            jy = jy + 1
            ls = l.split()
            zm.append([])
            for c in ls:
                zm[jy-1].append(float(c))

    zxm = array(zm)
    print 'Shape zxm = ',nmp.shape(zxm), '\n'
    return zxm



def __vcontour__(zmin, zmax, zdc):
    #
    #
    lngt = zmax - zmin
    #
    ncont = lngt/zdc
    #
    vcont = nmp.arange(zmin, zmax + zdc, zdc)
    #
    #lat_min
    return vcont



def __name_coor_ticks__(lon_min=0., lon_max=360., dlon=30., lat_min=-90., lat_max=90., dlat=15., lon_ext=0):
    #
    # Builds nice ticks for X and Y (lon, lat) axes!
    #
    # Arrange longitude axis !
    VX = nmp.arange(lon_min, lon_max+lon_ext+dlon, dlon); VX0 = nmp.arange(lon_min, lon_max+lon_ext+dlon, dlon);
    ivf = nmp.where(VX>180); VX0[ivf] = VX[ivf] - 360
    cn_lon = []
    for rlon in VX0:
        jlon = int(rlon)
        if jlon < 0:
            cn_lon.append(str(-jlon)+r'$^{\circ}$W')
        else:
            if jlon == 0:
                cn_lon.append(str(jlon)+r'$^{\circ}$')
            else:
                cn_lon.append(str(jlon)+r'$^{\circ}$E')
    #
    # Arrange latitude axis !
    VY = nmp.arange(lat_min, lat_max+dlat, dlat)
    cn_lat = []
    for rlat in VY:
        jlat = int(rlat)
        if jlat < 0:
            cn_lat.append(str(-jlat)+r'$^{\circ}$S')
        else:
            if jlat == 0:
                cn_lat.append(str(jlat)+r'$^{\circ}$')
            else:
                cn_lat.append(str(jlat)+r'$^{\circ}$N')
    #
    return VX, VY, cn_lon, cn_lat



def __give_proj__(cname):
    #
    #
    nb =nmp.shape(projection_def)[0] ; #print 'nb =', nb
    #
    # Initializing :
    vproj = [ 'NC', 'NC', 0.,  0.,  0.,  0.,  0.,  0., 'NC' ]
    #
    #
    jb = 0
    while jb < nb :
        if projection_def[jb][0] == cname:
            break
        else :
            jb = jb + 1
    #
    if jb == nb :
        print 'Zone "'+cname+'" does not exist!\n'
        print 'so far choice is :'
        for jb in range(nb): print projection_def[jb][0]
        sys.exit(0)
        #
        #
    vproj = projection_def[jb][:]
    #
    #print 'For ', projection_def[jb][0], ' we have vproj =', vproj, '\n'
    #
    return vproj





def __font_unity__():
    #
    params = {'font.family':'Arial','font.size':16,'xtick.labelsize':16,'ytick.labelsize':16,'axes.labelsize':16}
    mpl.rcParams.update(params)    
    big_fixed_fonts = { 'fontname':'Arial',       'fontweight':'normal', 'fontsize':16 }
    label_fonts     = { 'fontname':'Arial', 'fontweight':'normal', 'fontsize':16 }
    colorbar_fonts  = { 'fontname':'Arial',       'fontweight':'normal', 'fontsize':16 }    
    return big_fixed_fonts, label_fonts, colorbar_fonts




def __force_min_and_max__(rm, rp, Xin):
    idx_bad  = nmp.where(nmp.logical_not(nmp.isfinite(Xin)))
    Xin[idx_bad] = 0.
    idx1 = nmp.where(Xin <= rm); Xin[idx1] = rm + abs(rp-rm)*1.E-4
    idx2 = nmp.where(Xin >= rp); Xin[idx2] = rp - abs(rp-rm)*1.E-4
    Xin[idx_bad] = nmp.nan


def __font_unity_big__():
    #
    params = {'font.family':'Trebuchet MS','text.fontsize':20,'xtick.labelsize':16,'ytick.labelsize': 16,'axes.labelsize':18}
    mpl.rcParams.update(params)
    title_fonts     = { 'fontname':'Trebuchet MS', 'fontweight':'normal', 'fontsize':20 }
    big_fixed_fonts = { 'fontname':'monaco',       'fontweight':'normal', 'fontsize':20 }
    label_fonts     = { 'fontname':'Trebuchet MS', 'fontweight':'normal', 'fontsize':16 }
    colorbar_fonts  = { 'fontname':'Tahoma',       'fontweight':'normal', 'fontsize':14 }    
    return title_fonts, big_fixed_fonts, label_fonts, colorbar_fonts
    



#    params = { 'font.family': 'Ubuntu Mono',
#               'legend.fontsize': 14,
#               'text.fontsize':   14,
#               'xtick.labelsize': 12,
#               'ytick.labelsize': 12,
#               'axes.labelsize':  14}
#    mpl.rcParams.update(params)
#
#    font_ttl = { 'fontname':'Bitstream Vera Sans Mono', 'fontweight':'normal', 'fontsize':14 }
#    font_ylb = { 'fontname':'Tahoma', 'fontweight':'normal', 'fontsize':12 }
