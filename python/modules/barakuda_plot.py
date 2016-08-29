
########################################################################
#
#   B a r a K u d a
#
#   Ploting functions and utilities
#
## Authors:
#  --------
# 2010-2015: Laurent Brodeau (original primitive code)
#      2016: Saeed Falahat   (update to fancy grown-up coding!) :D
#
#######################################################################

import os
import sys
import numpy as nmp

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import barakuda_orca as bo


# Some defaults:
WDTH_DEF     = 10.
HGHT_DEF     =  4.
FIG_SIZE_DEF = ( WDTH_DEF , HGHT_DEF )
DPI_DEF      = 120
AXES_DEF     = [0.1, 0.082, 0.87, 0.84]

# Colors for line:    (https://www.daftlogic.com/projects-hex-colour-tester.htm)
b_blu = '#2C558A'
b_red = '#AD0000'
b_gre = '#3A783E' ; # b_gre = '#42BD82'
b_prp = '#8C008C'
b_org = '#ED7C4C'

v_dflt_colors = [b_blu, b_red, b_gre, b_org, b_prp, 'pink', '0.5', 'b', 'g', 'brown', 'orange' ]  ; # extend if more...

# Some projections to use with BaseMap:
#
#         zone PROJ llcrnrlon llcrnrlat urcrnrlon urcrnrlat   lat1  lon0 mer/par continent-res
#                               (lcc = Lambert conformal conic)
projection_def = [
         ['nseas',   'lcc',  -55., 40., 55., 75.,           60., -20., 10., 'l' ],   # Nordic seas
         ['natarct', 'lcc', -60., 40., 80., 72.,            55., -32., 10., 'l' ],   # NATL + Arctic
         ['labir',   'lcc',  -62., 48., -10., 75.,          50., -30.,  5., 'l' ],
         ['labsp',   'lcc',  -60., 48., 50., 75.5,          50., -30., 10., 'l' ],
         ['npol',    'stere', -75., 45., 100., 60.,          80., -30., 10., 'l' ],
         ['npol2',   'stere', -55., 40., 145., 40.,          80.,  -5., 10., 'l' ],   # North Pole
         ['spstere', 'stere',  0.,  0.,  0., 0.,            -48.,  90., 10., 'l' ],   # South Pole Default matplotlib!
         ['matl' ,   'cyl',  -82.,-21.,  12., 79.,          30., -30., 15., 'l' ],   # Nordic seas
         ['atmed',   'lcc',  -18., 33.,  -2., 42.,          30., -10.,  5., 'h' ],
         ['kav7' ,   'kav',    0.,  0.,   0.,  0.,           0.,   0.,  0., 'l' ] ] # global map-monde







class plot :
    ''' This class encapsulates all the plot routines
    In order to use it you need to type as follows:

    plot(function name without the prefix __) (all the arguments of the function)
    for example for __vert_section we do as follows:

    plot("vert_section")(VX, VZ, XF, XMSK, rmin, rmax, dc, lkcont=True, cpal='jet',
        xmin=-80., xmax=85., dx=5, cfignm='fig', cbunit='', cxunit=' ',
        zmin = 0., zmax = 5000., l_zlog=False, cfig_type='png',
        czunit=' ', ctitle=' ', lforce_lim=False, i_cb_subsamp=1, l_z_increase=False )

    The reason to prefix all the function with double underscore __ is that all these function become private members of
    the plot class and they can not be accessed directly outside of the class. That is, you need to call these functios through the call wrapper
    as seen below in the function __call__.
    '''

    def __init__(self,splot) :

        self.splot = splot


    def __call__(self,*args, **kw) :

        if "_"+self.__class__.__name__+ "__" + self.splot in self.__class__.__dict__.keys() :

            self.__class__.__dict__["_"+self.__class__.__name__+ "__" + self.splot](self,*args, **kw)

        else :
            print "function " + "__" + self.splot + " does not exist"
            sys.exit()



    # Functions:


    def __vert_section(self,VX, VZ, XF, XMSK, rmin, rmax, dc, lkcont=True, cpal='jet',
                       xmin=-80., xmax=85., dx=5, cfignm='fig', cbunit='', cxunit=' ',
                       zmin = 0., zmax = 5000., l_zlog=False, cfig_type='png',
                       czunit=' ', ctitle=' ', lforce_lim=False, i_cb_subsamp=1, l_z_increase=False ):

        import barakuda_colmap as bcm

        font_ttl, font_xylb, font_clb = __font_unity__()

        zVZ = __prepare_z_log_axis__(l_zlog, VZ)

        if lforce_lim: __force_min_and_max__(rmin, rmax, XF)

        # Masking where mask is zero!
        XF = nmp.ma.masked_where(XMSK == 0, XF)

        fig = plt.figure(num = 1, figsize=(WDTH_DEF,5.), dpi=None, facecolor='w', edgecolor='k') ; #vert_section
        ax  = plt.axes([0.085, 0.06, 0.98, 0.88], axisbg = 'gray')
        vc  = __vcontour__(rmin, rmax, dc)

        # Colormap:
        palette = bcm.chose_palette(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        cf = plt.contourf(VX, zVZ, XF, vc, cmap = palette, norm = pal_norm)
        plt.hold(True)
        if lkcont: plt.contour(VX, zVZ, XF, vc, colors='k', linewidths=0.2)

        # Colorbar:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, cunit=cbunit, cfont=font_clb, fontsize=10)

        # X-axis:
        __nice_x_axis__(ax, plt, xmin, xmax, dx, cunit=cxunit, cfont=font_xylb)

        # Y-axis:
        plt.ylabel(czunit, **font_xylb)

        # Fixing z ticks:
        __fix_z_axis__(ax, plt, zmin, zmax, l_log=l_zlog, l_z_inc=l_z_increase)

        plt.title(ctitle, **font_ttl)
        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=False) ; #vert_section
        print cfignm+'.'+cfig_type+' created!\n'
        plt.close(1)

        return


    def __2d(self,VX, VY, XF, XMSK, rmin, rmax, dc, corca='ORCA1', lkcont=True, cpal='jet',
             cfignm='fig', cbunit='', ctitle=' ', lforce_lim=False, i_cb_subsamp=1, cb_orient='vertical',
             cfig_type='pdf', lat_min=-75., lat_max=75., lpix=False, vcont_spec = []):

        #
        # Plot nicely a field given on ORCA coordinates on 2D world map without using any projection
        #
        # if VX = [0] and VY = [0] => ignoring lon and lat...


        import barakuda_tool   as bt
        import barakuda_colmap as bcm


        font_ttl, font_xylb, font_clb = __font_unity__()

        i_lat_lon = 1
        if len(VX) == 1 or len(VY) == 1: i_lat_lon = 0 ; # no long. and lat. provided !
        if corca[:5] == 'eORCA':
            # don't know how to lat-lon 2d plot eORCA
            # so just plot without projections
            i_lat_lon = 0



        # Don't want to modify XF array, working with XFtmp:
        [ny, nx] = nmp.shape(XF)
        XFtmp = nmp.zeros(ny*nx) ; XFtmp.shape = [ny, nx]
        XFtmp[:,:] = XF[:,:]

        # First drowning the field:
        bt.drown(XFtmp, XMSK, k_ew=2, nb_max_inc=20, nb_smooth=10)

        rlon_ext = 32.

        if lforce_lim: __force_min_and_max__(rmin, rmax, XFtmp)

        XMSK0 = bo.lon_reorg_orca(XMSK,corca, rlon_ext)
        XF0   = bo.lon_reorg_orca(XFtmp,  corca, rlon_ext)

        if i_lat_lon == 1:
            [ny, nx] = nmp.shape(XF0)
            dlong  = abs(VX[11] - VX[10])
            #VX0 = nmp.arange(0.,nx,dlong) + dlong/2. ; #lolo:
            VX0 = nmp.arange(0.,nx)
            VX0 = VX0*dlong + dlong/2.

        if i_lat_lon == 1:
            vert_rat = (lat_max - lat_min)/(75. + 75.)
            fig_size = (WDTH_DEF,4.8*vert_rat)    ; #2d  lulu
        else:
            fig_size = (WDTH_DEF , float(nx)/float(ny)*5.) ; #2d


        # FIGURE
        # ~~~~~~
        fig = plt.figure(num = 1, figsize=fig_size, dpi=None, facecolor='w', edgecolor='k')
        ax  = plt.axes([0.05, 0.06, 1., 0.86], axisbg = 'white')

        vc = __vcontour__(rmin, rmax, dc); #print vc, '\n'

        # Palette:
        palette = bcm.chose_palette(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)


        if lpix:
            # Pixelized plot:
            if i_lat_lon == 1:
                cf = plt.pcolor(VX0, VY, XF0, cmap = palette, norm = pal_norm)
            else:
                cf = plt.pcolor(         XF0, cmap = palette, norm = pal_norm)

        else:
            # Contour fill plot:
            if i_lat_lon == 1:
                cf = plt.contourf(VX0, VY, XF0, vc, cmap = palette, norm = pal_norm)
            else:
                cf = plt.contourf(XF0, vc, cmap = palette, norm = pal_norm)

            for c in cf.collections: c.set_zorder(0.15)

            if lkcont:
                if i_lat_lon == 1:
                    cfk = plt.contour(VX0, VY, XF0, vc, colors='k', linewidths = 0.2)
                else:
                    cfk = plt.contour(XF0, vc, colors='k', linewidths = 0.2)

                for c in cfk.collections: c.set_zorder(0.25)

            # contour for specific values on the ploted field:
            if len(vcont_spec) >= 1:
                if i_lat_lon == 1:
                    cfs = plt.contour(VX0, VY, XF0, vcont_spec, colors='black', linewidths = 1.5)
                else:
                    cfs = plt.contour(XF0, vcont_spec, colors='black', linewidths = 1.5)

                plt.clabel(cfs, inline=1, fmt='%4.1f', fontsize=10)
                for c in cfs.collections: c.set_zorder(0.35)

            ## Contours of continents:
            ##cfm = plt.contour(VX0, VY, XMSK0, [ 0.7 ], colors='k', linewidths = 0.4)
            ##for c in cfm.collections: c.set_zorder(1.)



        # Putting land-sea mask on top of current plot, cleaner than initial masking...
        # because won't influence contours since they are done
        # field needs to be DROWNED prior to this though!!!
        idx_land = nmp.where(XMSK0[:,:] < 0.5)
        XF0 = nmp.ma.masked_where(XMSK0[:,:] > 0.5, XF0)
        XF0[idx_land] = 1000.
        if i_lat_lon == 1:
            cf0 = plt.pcolor(VX0, VY, XF0, cmap = bcm.chose_palette("mask"))
        else:
            cf0 = plt.pcolor(XF0, cmap = bcm.chose_palette("mask"))

        # Colorbar:
        ifsize = 14
        if i_lat_lon == 1: ifsize = int(ifsize*vert_rat); ifsize=max(ifsize,6)
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, lkc=lkcont, cb_or=cb_orient, cunit=cbunit, cfont=font_clb, fontsize=ifsize)

        # X and Y nice ticks:
        if i_lat_lon == 1:
            [vvx, vvy, clon, clat] = __name_coor_ticks__(lon_ext=rlon_ext);
            plt.yticks(vvy,clat) ; plt.xticks(vvx,clon)
            plt.axis([ 0., 360.+rlon_ext-2., lat_min, lat_max])
        else:
            #ax.set_xlim(0., 360.+rlon_ext-2.)
            plt.axis([ 0., float(nx)+rlon_ext-2., 0, float(ny)])

        plt.title(ctitle, **font_ttl)

        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=False) ; #2d

        print cfignm+'.'+cfig_type+' created!\n'
        plt.close(1)

        del XFtmp, XF0

        return


    def __2d_reg(self,VX, VY, XF, XMSK, rmin, rmax, dc, lkcont=False, cpal='jet',
                 cfignm='fig', cfig_type='pdf', cbunit=' ', ctitle='',
                 cb_orient='vertical', lat_min=-77., lat_max=77., i_cb_subsamp=1,
                 lpix=False, l_continent_pixel=True, colorbar_fs=14,
                 col_min='k', col_max='k', vcont_spec = []):
        
        import barakuda_tool   as bt
        import barakuda_colmap as bcm

        font_ttl, font_xylb, font_clb = __font_unity__()


        # Don't want to modify XF array, working with XFtmp:
        [ny, nx] = nmp.shape(XF)
        XFtmp = nmp.zeros(ny*nx) ; XFtmp.shape = [ny, nx]
        XFtmp[:,:] = XF[:,:]

        # First drowning the field:
        bt.drown(XFtmp, XMSK, k_ew=0, nb_max_inc=20, nb_smooth=10)


        iskp = 28 ; iext = 32

        # Extending / longitude:
        VXe   = bt.extend_domain(VX,    iext, skp_west_deg=iskp) ; nxe = len(VXe)
        XFe   = bt.extend_domain(XFtmp, iext, skp_west_deg=iskp)
        XMSKe = bt.extend_domain(XMSK,  iext, skp_west_deg=iskp)



        # FIGURE

        rat_vert = 1. / ( ( 77. + 77. ) / ( lat_max - lat_min ) )


        if cb_orient == 'horizontal':
            # Horizontal colorbar!
            if ctitle == '':
                fig = plt.figure(num = 1, figsize=(12.4,7.*rat_vert), dpi=None, facecolor='w', edgecolor='k')
                ax = plt.axes([0.05, -0.01, 0.93, 1.], axisbg = 'white')
            else:
                fig = plt.figure(num = 1, figsize=(12.4,7.4*rat_vert), dpi=None, facecolor='w', edgecolor='k')
                ax = plt.axes([0.05, -0.01, 0.93, 0.96], axisbg = 'white')
        else:
            # Vertical colorbar!
            fig = plt.figure(num = 1, figsize=(12.4,6.*rat_vert), dpi=None, facecolor='w', edgecolor='k')
            ax = plt.axes([0.046, 0.06, 1.02, 0.88], axisbg = 'white')

        vc = __vcontour__(rmin, rmax, dc)

        # Palette:
        #palette = bcm.chose_palette(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)
        mpl.rcParams['contour.negative_linestyle'] = 'solid'
        plt.contour.negative_linestyle='solid'

        if lpix:
            cf = plt.pcolormesh(VXe, VY, XFe, cmap = cpal, norm = pal_norm)
        else:
            cf = plt.contourf(VXe, VY, XFe, vc, cmap = cpal, norm = pal_norm, extend="both")
            for c in cf.collections: c.set_zorder(0.15)
            cf.cmap.set_under(col_min)
            cf.cmap.set_over(col_max)

        # contour for specific values on the ploted field:
        if len(vcont_spec) >= 1:
            cfs = plt.contour(VXe, VY, XFe, vcont_spec, colors='w', linewidths = 1.)
            #plt.clabel(cfs, inline=1, fmt='%4.1f', fontsize=12)


        if lkcont:
            cfk = plt.contour(VXe, VY, XFe, vc, colors='k', linewidths = 0.2)
            for c in cfk.collections: c.set_zorder(0.25)




        # Putting land-sea mask on top of current plot, cleaner than initial masking...
        # because won't influence contours since they are done
        # field needs to be DROWNED prior to this though!!!

        if l_continent_pixel:
            idx_land = nmp.where(XMSKe[:,:] < 0.5)
            XFe = nmp.ma.masked_where(XMSKe[:,:] > 0.5, XFe)
            XFe[idx_land] = 1000.
            cf0 = plt.pcolor(VXe, VY, XFe, cmap = bcm.chose_palette("mask"))
        else:
            # Masking with contour rather than pixel:
            cf0 = plt.contourf(VXe, VY, XMSKe, [ 0., 0.1 ], cmap = bcm.chose_palette("mask"))
            for c in cf0.collections: c.set_zorder(5)
            plt.contour(VXe, VY, XMSKe, [ 0.25 ], colors='k', linewidths = 1.)


        # COLOR BAR:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, lkc=lkcont, cb_or=cb_orient, cunit=cbunit, cfont=font_clb, fontsize=colorbar_fs)

        # X and Y nice ticks:
        print "VXe[0], VXe[nxe-1] =>", VXe[0], VXe[nxe-1]
        rlon_min = round(VXe[0],0) ; rlon_max = round(VXe[nxe-1],0)
        print "rlon_min, rlon_max =>", rlon_min, rlon_max

        [vvx, vvy, clon, clat] = __name_coor_ticks__(lon_min=rlon_min, lon_max=rlon_max, dlon=30., lon_ext=iext-iskp)
        plt.yticks(vvy,clat) ; plt.xticks(vvx,clon)

        plt.axis([ rlon_min, rlon_max, lat_min, lat_max])


        if ctitle != ' ': plt.title(ctitle, **font_ttl)

        # Prevents from using scientific notations in axess ticks numbering:
        #ax.get_xaxis().get_major_formatter().set_useOffset(False)

        plt.savefig(cfignm+'.'+cfig_type, dpi=110, orientation='portrait', transparent=False)

        print cfignm+'.'+cfig_type+' created!\n'
        plt.close(1)

        return


    def __2d_box(self,XF, XMSK, rmin, rmax, dc, lkcont=True,
                 cpal='jet', cfignm='fig', cbunit='', ctitle=' ', lforce_lim=False,
                 i_cb_subsamp=1, cfig_type='pdf', lcontours=True,
                 x_offset=0., y_offset=0., vcont_spec = [], lcont_mask=False):

        import barakuda_colmap as bcm



        if lforce_lim: __force_min_and_max__(rmin, rmax, XF)

        [ ny , nx ] = XF.shape
        vert_rat = float(ny)/float(nx)
        print "Vert. ratio, nx, ny =", vert_rat, nx, ny

        # Masking field:
        if lcontours:
            idxm = nmp.where(XMSK[:,:] == 0); XF[idxm] = -9999.9  # c'est NaN qui merde!!!
        else:
            XF = nmp.ma.masked_where(XMSK == 0, XF)


        font_ttl, font_xylb, font_clb = __font_unity__()


        # FIGURE
        # ~~~~~~
        fig = plt.figure(num = 1, figsize=(7.,6.*vert_rat), dpi=None, facecolor='w', edgecolor='k')

        ax = plt.axes([0.07, 0.05, 0.9, 0.9], axisbg = 'gray')

        vc = __vcontour__(rmin, rmax, dc); #print vc, '\n'

        # Palette:
        palette = bcm.chose_palette(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        plt.hold(True)


        if lcontours:
            cf = plt.contourf(XF, vc, cmap = palette, norm = pal_norm)
            for c in cf.collections: c.set_zorder(0.5)
        else:
            cf = plt.pcolor(XF, cmap = palette, norm = pal_norm)




        # contour for specific values on the ploted field:
        if len(vcont_spec) >= 1:
            cfs = plt.contour(XF, vcont_spec, colors='white', linewidths = 1.)
            plt.clabel(cfs, inline=1, fmt='%4.1f', fontsize=10)



        if lkcont:
            cfk = plt.contour(XF, vc, colors='k', linewidths = 0.1)
            for c in cfk.collections: c.set_zorder(0.75)



        # contour for continents:
        if lcontours and lcont_mask:
            cfm = plt.contour(XMSK, [ 0.7 ], colors='k', linewidths = 1.)
            for c in cfm.collections: c.set_zorder(1.)

        # Colorbar:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, lkc=lkcont, cunit=cbunit, cfont=font_clb)

        if x_offset > 0 or y_offset > 0 :  __add_xy_offset__(plt, x_offset, y_offset)

        plt.axis([ 0., nx-1, 0, ny-1])

        plt.title(ctitle, **font_ttl)

        # Prevents from using scientific notations in axess ticks numbering:
        ax.get_xaxis().get_major_formatter().set_useOffset(False)

        plt.savefig(cfignm+'.'+cfig_type, dpi=100, orientation='portrait', transparent=True)

        print cfignm+'.'+cfig_type+' created!\n'
        plt.close(1)

        return


    def __zonal(self,VY, VZn, VZ1=[0.], VZ2=[0.], VZ3=[0.],
                cfignm='fig_zonal', zmin=-100., zmax=100., dz=25., i_z_jump=1,
                xmin=-90., xmax=90., cfig_type='png', cxunit=r'Latitude ($^{\circ}$N)',
                czunit='', ctitle='', lab='', lab1='', lab2='', lab3='', box_legend=(0.6, 0.75),
                lnarrow=False):

        font_ttl, font_xylb, font_clb = __font_unity__()

        ny = len(VY)
        if len(VZn) != ny: print 'ERROR: plot_zonal.barakuda_plot => VY and VZn do not agree in size'; sys.exit(0)

        lp1=False ; lp2=False ; lp3=False
        if len(VZ1) == ny: lp1=True
        if len(VZ2) == ny: lp2=True
        if len(VZ3) == ny: lp3=True

        if lnarrow:
            fig = plt.figure(num = 1, figsize=(WDTH_DEF/2.,4.), dpi=None)  ; #zonal
            ax  = plt.axes([0.14 , 0.12, 0.85, 0.82])   #, axisbg = 'gray')
        else:
            fig = plt.figure(num = 1, figsize=(WDTH_DEF,6.), dpi=None)  ; #zonal
            ax = plt.axes([0.08, 0.11, 0.9, 0.82])   #, axisbg = 'gray')

        plt.plot(VY, VZn*0.0, 'k', linewidth=1) ; plt.hold(True)

        plt.plot(VY, VZn, 'k', linewidth=3., label=lab)
        if lp1: plt.plot(VY, VZ1, color=b_blu, linewidth=2., label=lab1)
        if lp2: plt.plot(VY, VZ2, color=b_red, linewidth=2., label=lab2)
        if lp3: plt.plot(VY, VZ3, color=b_gre, linewidth=2., label=lab3)

        plt.legend(bbox_to_anchor=box_legend, shadow=False, fancybox=True)

        # X-axis
        __nice_x_axis__(ax, plt, xmin, xmax, 15., cunit=cxunit, cfont=font_xylb)

        # Z-axis:
        __nice_z_axis__(ax, plt, zmin, zmax, dz, i_sbsmp=i_z_jump, cunit=czunit, cfont=font_xylb)

        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.title(ctitle, **font_ttl)

        # Prevents from using scientific notations in axess ticks numbering:
        #lolo: ax.get_xaxis().get_major_formatter().set_useOffset(False)

        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=False)  ; #zonal

        plt.close(1)

        return







    def __nproj(self,czone, rmin, rmax, dc, xlon, xlat, XF,
                cfignm='fig', lkcont=False, cpal='jet', cbunit=' ',
                cfig_type='pdf', ctitle=' ', lforce_lim=False,
                cb_orient='vertical', i_cb_subsamp=1, dpi_fig=DPI_DEF, lpcont=True):

        # Plot projection with basemap...

        #===================================================================================
        # INPUT:
        #          xlon and xlat can be 1D or 2D !!!
        #
        #   lpcont=True  => do contourf
        #   lpcont=False => do pcolor
        #
        #===================================================================================



        from mpl_toolkits.basemap import Basemap
        from mpl_toolkits.basemap import shiftgrid
        import barakuda_colmap as bcm

        font_ttl, font_xylb, font_clb = __font_unity__()

        # For projections :
        vp = __give_proj__(czone) ; # projection information


        # must work with XFtmp rather than XF, because sometimes XF is overwrited...
        [ny, nx] = nmp.shape(XF)
        XFtmp = nmp.zeros(ny*nx) ; XFtmp.shape = [ny, nx]
        XFtmp[:,:] = XF[:,:]


        if len(nmp.shape(xlat)) == 1 and len(nmp.shape(xlon)) == 1:
            if czone == 'kav7' and xlon[0] >= 0.:
                # Shifting data and longitude to be consistent with map projection
                XFtmp, xlon = shiftgrid(180.+xlon[0], XFtmp, xlon, start=False, cyclic=360.0)
            LON_2D, LAT_2D = nmp.meshgrid(xlon,xlat)
        else:
            LAT_2D = nmp.zeros(ny*nx) ; LAT_2D.shape = [ny, nx] ; LAT_2D[:,:] = xlat[:,:]
            LON_2D = nmp.zeros(ny*nx) ; LON_2D.shape = [ny, nx] ; LON_2D[:,:] = xlon[:,:]


        if lforce_lim: __force_min_and_max__(rmin, rmax, XFtmp)

        vc = __vcontour__(rmin, rmax, dc)

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
            rw = WDTH_DEF/2
            vfig_size                        = [ rw, rw ];        vsporg = [0.1, 0.1, 0.85, 0.85]
            if czone == 'nseas':   vfig_size = [ rw , rw ];       vsporg = [0.085,  0.03, 0.9, 0.94]
            if czone == 'natarct': vfig_size = [ rw , rw ];       vsporg = [0.065,  0.04, 0.95, 0.92]
            if czone == 'spstere': vfig_size = [ rw , 0.81*rw ];  vsporg = [0.05, 0.05, 0.94, 0.89]
            if czone == 'npol2':   vfig_size = [ rw , 0.96*rw ] ; vsporg = [0.1,  0.05, 0.86, 0.9]
            #if czone == 'kav7':    vfig_size = [ 7., 5.  ]; vsporg = [0.085, 0.03, 0.91, 0.94]
            #lulu


        fig = plt.figure(num = 1, figsize=(vfig_size), dpi=None, facecolor='w', edgecolor='k')
        ax  = plt.axes(vsporg, axisbg = 'w')


        ## Palette:
        palette = bcm.chose_palette(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)
        mpl.rcParams['contour.negative_linestyle'] = 'solid'; plt.contour.negative_linestyle='solid'






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

        if lpcont:
            cf = carte.contourf(x0, y0, XFtmp, vc, cmap = palette, norm = pal_norm)
            # Black contours if needed :
            if lkcont:
                ckf = carte.contour(x0, y0, XFtmp, vc, colors='k', linewidths=0.5)
                if cpal != 'ice':
                    for c in cf.collections: c.set_zorder(0.5)   # Changing zorder so black cont. on top
                for c in ckf.collections: c.set_zorder(1.) # of filled cont. and under continents (zorder 1)

        else:
            cf = carte.pcolor(x0, y0, XFtmp, cmap = palette, norm = pal_norm)




        carte.drawcoastlines() ; carte.fillcontinents(color='grey') ; carte.drawmapboundary()





        if vp[1] == 'lcc' or vp[1] == 'cyl':
            carte.drawmeridians(nmp.arange(-360,360,vp[8]), labels=[0,0,0,1])
            carte.drawparallels(nmp.arange( -90, 90,vp[8]), labels=[1,0,0,0])

        if vp[1] == 'stere':
            carte.drawmeridians(nmp.arange(-180,180,20), labels=[0,0,0,1])
            carte.drawparallels(nmp.arange( -90, 90,10), labels=[1,0,0,0])

        plt.title(ctitle, **font_ttl)

        # Colorbar:
        if cb_orient == 'horizontal':
            clbax = fig.add_axes(vcbar) # new axes for colorbar!
            __nice_colorbar__(cf, plt, vc, cax_other=clbax, i_sbsmp=i_cb_subsamp, lkc=(lkcont and lpcont), cb_or='horizontal', cunit=cbunit, cfont=font_clb, fontsize=10)
        else:
            __nice_colorbar__(cf, plt, vc,                  i_sbsmp=i_cb_subsamp, lkc=(lkcont and lpcont),                     cunit=cbunit, cfont=font_clb, fontsize=12)

        plt.savefig(cfignm+'.'+cfig_type, dpi=dpi_fig, orientation='portrait', transparent=False) ; #, transparent=True, acecolor='w', edgecolor='w',trans

        plt.close(1)

        print ' *** created figure '+cfignm+'.'+cfig_type+'\n'

        del LON_2D, LAT_2D, XFtmp

        return





    def __sig_transport(self,cstck, ccr, cd_eps):

        #_____________________________________________________________
        #
        #  cstck  : root directory containing monitoring output
        #  crr    : CONF-RUN  (ex: ORCA2-GXXX)
        #  cd_eps : directory to save eps output image
        #_____________________________________________________________

        ccnf = __message__(ccr)

        cf = cstck+'/'+ccnf+'/'+ccr+'-'+'MONITOR/'+ccr+'_TRPSIG.mtl'
        print 'File =', cf, '\n'

        xmat = __get_mat__(cf)
        [nl,nc] =nmp.shape(xmat)
        print nl,nc

        vsigma = xmat[0,1:nc]   ; print 'vsigma =', vsigma, '\n' ; nbs = len(vsigma)
        vyears = xmat[1:nl:2,0] ; print 'vyears =', vyears, '\n' ; nby = len(vyears)

        print 'Number of years and sigma levels =', nby, 'and', nbs, '\n'

        X_ds = -xmat[1:nl:2,1:nc]   # Denmark Straight
        X_fb = -xmat[2:nl:2,1:nc]   # Faro Banks channel

        X_ds = nmp.ma.masked_where(X_ds  == 0, X_ds)
        X_fb = nmp.ma.masked_where(X_fb  == 0, X_fb)

        print 'Size of X_ds =',nmp.shape(X_ds)
        print 'Size of X_fb =',nmp.shape(X_fb)

        cf_eps = cd_eps+'/tsc_DS_'+ccr+'.png'
        xmt=transpose(X_ds)

        __LB_2D__(1, 0., 1.5, 30, 0.1, vyears, vsigma, xmt,
                  'Transport by sigma class, Denmark Straight', cf_eps, cbunit='(Sv)')

        cf_eps = cd_eps+'/tsc_FB_'+ccr+'.eps'

        xmt=transpose(X_fb)
        __LB_2D__(1, 0., 1.5, 40, 0.1, vyears, vsigma, xmt,
                  'Transport by sigma class, Denmark Straight', cf_eps,  cbunit='(Sv)')
        #
        #
        return





    def __amoc_lat_depth(self,VY, VZ, Xamoc, XMSK, rmin, rmax, dc, lkcont=True, cpal='jet', ymin=-80., ymax=85.,
                            cfignm='fig', cbunit='', cxunit=' ', zmin = 0., zmax = 5000., l_zlog=False,
                            cfig_type='pdf', czunit=' ', ctitle=' ', lforce_lim=False, i_cb_subsamp=1):

        import matplotlib.colors as colors   # palette and co.
        import barakuda_colmap as bcm

        zVZ = __prepare_z_log_axis__(l_zlog, VZ)

        Xamoc = nmp.ma.masked_where(XMSK == 0, Xamoc)

        if lforce_lim: __force_min_and_max__(rmin, rmax, Xamoc)

        font_ttl, font_xylb, font_clb = __font_unity__()

        fig = plt.figure(num = 1, figsize=(WDTH_DEF,6.), facecolor='w', edgecolor='k')  ; #amoc_lat_depth
        ax  = plt.axes([0.1,  0.1,  0.94, 0.85], axisbg = 'gray')

        vc = __vcontour__(rmin, rmax, dc)

        # Colormap:
        palette = bcm.chose_palette(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        # Plot:
        cf = plt.contourf(VY, zVZ, Xamoc, vc, cmap = palette, norm = pal_norm)
        if lkcont: plt.contour(VY, zVZ, Xamoc, vc, colors='k', linewidths=0.2)

        # Colorbar:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, cunit=cbunit, cfont=font_clb)  #, fontsize=14)

        plt.axis([ ymin, ymax, zmin, zmax])
        plt.xlabel(cxunit, **font_xylb); plt.ylabel(czunit, **font_xylb)

        # Fixing z ticks:
        __fix_z_axis__(ax, plt, zmin, zmax, l_log=l_zlog)

        plt.title(ctitle, **font_ttl)

        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=False)  ; #amoc_lat_depth
        print cfignm+'.'+cfig_type+' created!\n'

        plt.close(1)

        return







    def __2d_box_2f(self,XF1, XF2, XMSK, rmin, rmax, dc, vcont_spec2, corca='ORCA1', lkcont=True,
                    cpal='jet', cfignm='fig', cbunit='', ctitle=' ', lforce_lim=False,
                    i_cb_subsamp=1, cfig_type='pdf', lcontours=True,
                    x_offset=0., y_offset=0., vcont_spec1 = []):

        # Take 2 fields as imput and shows contours of second field (vcont_spec2) on top of field 1

        import matplotlib.colors as colors   # palette and co.
        import barakuda_colmap as bcm


        if nmp.shape(XF1) != nmp.shape(XF2):
            print 'ERROR barakuda_plot.plot_2d_box_2f: fields F1 and F2 dont have the same shape!'
            sys.exit(0)



        font_ttl, font_xylb, font_clb = __font_unity__()


        if lforce_lim: __force_min_and_max__(rmin, rmax, XF1)

        [ ny , nx ] = XF1.shape
        vert_rat = float(ny)/float(nx)
        print "Vert. ratio, nx, ny =", vert_rat, nx, ny

        # Masking field:
        if lcontours:
            idxm = nmp.where(XMSK[:,:] == 0); XF1[idxm] = -9999.9  # c'est NaN qui merde!!!
        else:
            XF1 = nmp.ma.masked_where(XMSK == 0, XF1)



        # FIGURE
        # ~~~~~~
        fig = plt.figure(num = 1, figsize=(7.,6.*vert_rat), dpi=None, facecolor='w', edgecolor='k')

        ax = plt.axes([0.07, 0.05, 0.9, 0.9], axisbg = 'gray')

        vc = __vcontour__(rmin, rmax, dc); #print vc, '\n'

        # Palette:
        palette = bcm.chose_palette(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        plt.hold(True)


        if lcontours:
            cf = plt.contourf(XF1, vc, cmap = palette, norm = pal_norm)
            for c in cf.collections: c.set_zorder(0.5)
        else:
            cf = plt.pcolor(XF1, cmap = palette, norm = pal_norm)

        # contour for specific values on the ploted field:
        if len(vcont_spec1) >= 1:
            cfs1 = plt.contour(XF1, vcont_spec1, colors='white', linewidths = 1.)
            plt.clabel(cfs1, inline=1, fmt='%4.1f', fontsize=10)

        # Contours of field F2:
        cfs2 = plt.contour(XF2, vcont_spec2, colors=b_red, linewidths = 1.3)
        #plt.clabel(cfs1, inline=1, fmt='%4.1f', fontsize=10)


        if lkcont:
            cfk = plt.contour(XF1, vc, colors='k', linewidths = 0.1)
            for c in cfk.collections: c.set_zorder(0.75)





        # contour for continents:
        if lcontours:
            cfm = plt.contour(XMSK, [ 0.7 ], colors='k', linewidths = 0.4)
            for c in cfm.collections: c.set_zorder(1.)

        # Colorbar:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, lkc=lkcont, cunit=cbunit, cfont=font_clb)

        if x_offset > 0 or y_offset > 0 :  __add_xy_offset__(plt, x_offset, y_offset)

        plt.axis([ 0., nx-1, 0, ny-1])

        plt.title(ctitle, **font_ttl)

        plt.savefig(cfignm+'.'+cfig_type, dpi=100, orientation='portrait', transparent=True)

        print cfignm+'.'+cfig_type+' created!\n'
        plt.close(1)

        del Xtmp

        return






    def __trsp_sig_class(self,VT, vsigma_bounds, XF, rmin, rmax, dc, dsig,
                         lkcont=True, cpal='sigtr', dt_year=5., cfignm='fig',
                         cfig_type='pdf', ctitle='', lforce_lim=False, vcont_spec1 = [],
                         i_cb_subsamp=2):

        # Plot transport by sigma class...

        import matplotlib.colors as colors   # palette and co.
        import barakuda_colmap as bcm

        font_ttl, font_xylb, font_clb = __font_unity__()

        if lforce_lim: __force_min_and_max__(rmin, rmax, XF)

        fig = plt.figure(num = 1, figsize=(WDTH_DEF,6.), dpi=None, facecolor='w', edgecolor='k') ; #trsp_sig_class
        ax = plt.axes([0.055,  -0.025, 0.92, 0.98], axisbg = 'white')


        vc = __vcontour__(rmin, rmax, dc); #print 'plot_time_depth_hovm: contours =>\n', vc, '\n'

        nbins = len(vsigma_bounds) - 1

        # Palette:
        palette = bcm.chose_palette(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)
        mpl.rcParams['contour.negative_linestyle'] = 'solid'
        plt.contour.negative_linestyle='solid'

        cf = plt.contourf(VT, vsigma_bounds[:-1], XF, vc, cmap = palette, norm = pal_norm)
        #cf = plt.pcolor(VT, vsigma_bounds[:-1], XF, cmap = palette, norm = pal_norm)
        if lkcont:
            cfc = plt.contour(VT, vsigma_bounds[:-1], XF, nmp.arange(-3.,3.,0.5), colors='k', linewidths=0.4)

        # contour for specific values on the ploted field:
        if len(vcont_spec1) >= 1:
            cfs1 = plt.contour(VT, vsigma_bounds[:-1], XF, vcont_spec1, colors='white', linewidths = 1.)
            plt.clabel(cfs1, inline=1, fmt='%4.1f', fontsize=11, manual=[(2080,2.)] )


        # COLOR BAR:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, lkc=True, cb_or='horizontal', cunit='Sv', cfont=font_clb, fontsize=12)

        # AXES:
        y1 = int(min(VT))  ; y2 = int(max(VT))+1
        plt.axis([y1, y2, vsigma_bounds[nbins], vsigma_bounds[0]])

        __nice_x_axis__(ax, plt, y1, y2, dt_year, cfont=font_xylb)

        plt.yticks( nmp.flipud(vsigma_bounds) )

        label_big = { 'fontname':'Trebuchet MS', 'fontweight':'normal', 'fontsize':18 }
        plt.ylabel(r'$\sigma_0$', **label_big)

        plt.title(ctitle, **font_ttl)
        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=True) ; #trsp_sig_class
        print cfignm+'.'+cfig_type+' created!\n'
        plt.close(1)

        return












    def __vert_section_extra(self,VX, VZ, XF, XMSK, Vcurve, rmin, rmax, dc, lkcont=True, cpal='jet', xmin=-80., xmax=85.,
                             cfignm='fig', cbunit='', cxunit=' ', zmin = 0., zmax = 5000., l_zlog=False,
                             cfig_type='pdf', czunit=' ', ctitle=' ', lforce_lim=False, fig_size=(8.,8.) ):

        import matplotlib.colors as colors   # palette and co.
        import barakuda_colmap as bcm

        zVZ = __prepare_z_log_axis__(l_zlog, VZ)

        XF = nmp.ma.masked_where(XMSK == 0, XF)

        if lforce_lim: __force_min_and_max__(rmin, rmax, XF)

        font_ttl, font_xylb, font_clb = __font_unity__()


        fig = plt.figure(num = 1, figsize=fig_size, dpi=None, facecolor='w', edgecolor='k')
        ax = plt.axes([0.1,  0.065,   0.92,       0.89], axisbg = 'gray')
        vc = __vcontour__(rmin, rmax, dc)

        # Palette:
        palette = bcm.chose_palette(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        cf = plt.contourf(VX, zVZ, XF, vc, cmap = palette, norm = pal_norm)
        plt.hold(True)
        if lkcont: plt.contour(VX, zVZ, XF, vc, colors='k', linewidths=0.2)

        # Colorbar:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, cunit=cbunit, cfont=font_clb, fontsize=10)

        # X-axis:
        __nice_x_axis__(ax, plt, xmin, xmax, dx, cunit=cxunit, cfont=font_xylb)

        plt.ylabel(czunit, **font_xylb)

        plt.plot(VX,Vcurve, 'w', linewidth=2)

        for zz in zVZ[:]:
            plt.plot(VX,VX*0.+zz, 'k', linewidth=0.3)

        # Fixing z ticks:
        __fix_z_axis__(ax, plt, zmin, zmax, l_log=l_zlog, l_z_inc=False)

        plt.title(ctitle, **font_ttl)
        plt.savefig(cfignm+'.'+cfig_type, dpi=100, orientation='portrait', transparent=True)
        print cfignm+'.'+cfig_type+' created!\n'
        plt.close(1)
        #
        return





    def __time_depth_hovm(self,VT, VZ, XF, XMSK, rmin, rmax, dc, lkcont=True, cpal='jet',
                          tmin=0., tmax=100., dt=5.,
                          cfignm='fig', cbunit='', cxunit=' ', zmin = 0., zmax = 5000., l_zlog=False,
                          cfig_type='pdf', czunit=' ', ctitle=' ', lforce_lim=False, i_cb_subsamp=1,
                          vcont_spec1 = [], col_cont_spec1='w'):

        import matplotlib.colors as colors   # palette and co.
        import barakuda_colmap as bcm

        font_ttl, font_xylb, font_clb = __font_unity__()

        zVZ = __prepare_z_log_axis__(l_zlog, VZ)

        XF = nmp.ma.masked_where(XMSK == 0, XF)

        if lforce_lim: __force_min_and_max__(rmin, rmax, XF)


        fig = plt.figure(num = 1, figsize=(WDTH_DEF,6.), dpi=None, facecolor='w', edgecolor='k') ; #time_depth_hovm
        ax = plt.axes([0.07,  0.08,   1., 0.82], axisbg = 'gray')


        vc = __vcontour__(rmin, rmax, dc); #print 'plot_time_depth_hovm: contours =>\n', vc, '\n'

        # Palette:
        palette = bcm.chose_palette(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        cf = plt.contourf(VT, zVZ, XF, vc, cmap = palette, norm = pal_norm)
        plt.hold(True)
        if lkcont: plt.contour(VT, zVZ, XF, vc, colors='k', linewidths=0.2)

        # contour for specific values on the ploted field:
        if len(vcont_spec1) >= 1:
            cfs1 = plt.contour(VT, zVZ, XF, vcont_spec1, colors=col_cont_spec1, linewidths = 2.)
            plt.clabel(cfs1, inline=1, fmt='%4.1f', fontsize=11, manual=[(2080,2.)] )

        ax.grid(color='k', linestyle='-', linewidth=0.1)

        # Colorbar:
        clb = plt.colorbar(cf, ticks=vc)
        if i_cb_subsamp > 1: __subsample_colorbar__(i_cb_subsamp, vc, clb, cb_or=cb_orient)
        for t in clb.ax.get_yticklabels(): t.set_fontsize(13)
        if cbunit != '': clb.set_label('('+cbunit+')', **font_clb)

        plt.axis([ tmin, tmax, zmax, zmin])
        plt.ylabel(czunit, **font_xylb)

        # X-axis:
        __nice_x_axis__(ax, plt, tmin, tmax, dt, cfont=font_xylb)

        # Fixing z ticks:
        __fix_z_axis__(ax, plt, zmin, zmax, l_log=l_zlog, l_z_inc=False)

        plt.title(ctitle, **font_ttl)
        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=False) ; #time_depth_hovm
        print cfignm+'.'+cfig_type+' created!\n'
        plt.close(1)

        return









    def __enso(self,VT, VSST, cfignm='fig', dt_year=5):

        font_ttl, font_xylb, font_clb = __font_unity__()

        Nt = len(VT)

        if len(VT) != len(VSST): print 'ERROR: plot_enso.barakuda_plot => VT and VSST do not agree in size'; sys.exit(0)

        print ' Nt =>', Nt

        # Array to contain nino series:
        xnino = nmp.zeros(Nt*4) ; xnino.shape = [ Nt, 4 ]

        xnino[:,0] = VSST[:]

        # 5-month running mean:
        for jt in nmp.arange(2,Nt-2):
            xnino[jt,1] = (xnino[jt-2,0] + xnino[jt-1,0] + xnino[jt,0] + xnino[jt+1,0] + xnino[jt+2,0]) / 5.

        xnino[0:2,1] = xnino[2,1] ; xnino[Nt-2:Nt,1] = xnino[Nt-3,1]

        print '\n'

        print 'mean value for sst mean = ', nmp.sum(xnino[:,0])/Nt
        print 'mean value for sst 5-m-r mean = ', nmp.sum(xnino[:,1])/Nt


        # least-square curve for 5-month running mean:
        sumx  = nmp.sum(VT[:]) ; sumy  = nmp.sum(xnino[:,1])
        sumxx = nmp.sum(VT[:]*VT[:])
        sumxy = nmp.sum(VT[:]*xnino[:,1])
        a = ( sumx*sumy - Nt*sumxy ) / ( sumx*sumx - Nt*sumxx )
        b = ( sumy - a*sumx ) / Nt
        print 'a, b =', a, b

        # least-square linear trend:
        xnino[:,2] = a*VT[:] + b
        print 'mean value for least-square linear trend = ', nmp.sum(xnino[:,2])/Nt

        # anomaly
        xnino[:,3] = xnino[:,1] - xnino[:,2] ; # anomaly for 5-month running mean
        print 'mean value for anomaly = ', nmp.sum(xnino[:,3])/Nt

        # FIGURE ENSO
        #############

        vsst_plus  = nmp.zeros(Nt) ; vsst_plus.shape = [ Nt ]
        vsst_minus = nmp.zeros(Nt) ; vsst_minus.shape = [ Nt ]

        vsst_plus[:]  = xnino[:,3]
        vsst_minus[:] = xnino[:,3]

        vsst_plus[nmp.where(xnino[:,3] < 0. )] = 0.
        vsst_minus[nmp.where(xnino[:,3] > 0. )] = 0.

        vsst_plus[0] = 0. ; vsst_minus[0] = 0.
        vsst_plus[Nt-1] = 0. ; vsst_minus[Nt-1] = 0.

        y1 = int(min(VT))
        y2 = int(max(VT)+0.25)

        fig = plt.figure(num = 2, figsize=FIG_SIZE_DEF, facecolor='w', edgecolor='k') ; #enso
        ax = plt.axes(AXES_DEF)

        xnino[:,0] =  0.4 ; plt.plot(VT, xnino[:,0], 'r--', linewidth=1.5)
        xnino[:,0] = -0.4 ; plt.plot(VT, xnino[:,0], 'b--', linewidth=1.5)

        plt.fill(VT, vsst_plus, b_red, VT, vsst_minus, b_blu, linewidth=0)
        plt.plot(VT, xnino[:,3], 'k', linewidth=0.7)
        xnino[:,3] = 0.0
        plt.plot(VT, xnino[:,3], 'k', linewidth=0.7)
        plt.axis([min(VT), max(VT), -2.5, 2.5])

        __nice_x_axis__(ax, plt, y1, y2, dt_year, cfont=font_xylb)

        plt.yticks( nmp.arange(-2.5,2.501,0.5) )

        ax.grid(color='k', linestyle='-', linewidth=0.2)
        plt.ylabel(r'SST anomaly ($^{\circ}$C)', **font_xylb)
        plt.title('SST anomaly on Nino region 3.4', **font_ttl)
        cf_fig = cfignm+'.png'
        plt.savefig(cf_fig, dpi=DPI_DEF, orientation='portrait', transparent=True) ; #enso

        plt.close(2)


        del xnino

        return









    def __1d_mon_ann(self,VTm, VTy, VDm, VDy, cfignm='fig', dt_year=5, cyunit='', ctitle='',
                        ymin=0, ymax=0, dy=0, i_y_jump=1, mnth_col='b', plt_m03=False, plt_m09=False,
                        cfig_type='png', l_tranparent_bg=True, fig_size=FIG_SIZE_DEF):

        # if you specify ymin and ymax you can also specify y increment (for y grid) as dy
        #
        # plt_m03 => plot march values on top in green
        # plt_m09 => plot september values on top in green

        font_ttl, font_xylb, font_clb = __font_unity__()

        Nt1 = len(VTm) ; Nt2 = len(VTy)

        if len(VTm) != len(VDm): print 'ERROR: plot_1d_mon_ann.barakuda_plot => VTm and VDm do not agree in size'; sys.exit(0)
        if len(VTy) != len(VDy): print 'ERROR: plot_1d_mon_ann.barakuda_plot => VTy and VDy do not agree in size'; sys.exit(0)

        fig = plt.figure(num = 1, figsize=fig_size, facecolor='w', edgecolor='k')

        ax = plt.axes(AXES_DEF) ; #1d_mon_ann

        if mnth_col == 'g': mnth_col = b_gre
        if mnth_col == 'b': mnth_col = b_blu

        plt.plot(VTm, VDm, mnth_col, label=r'monthly', linewidth=1)
        plt.plot(VTy, VDy, b_red, label=r'annual', linewidth=2)

        if plt_m03: plt.plot(VTm[2:Nt1:12], VDm[2:Nt1:12], 'orange', label=r'March',     linewidth=2)
        if plt_m09: plt.plot(VTm[8:Nt1:12], VDm[8:Nt1:12], 'orange', label=r'September', linewidth=2)

        if plt_m03 or plt_m09: plt.legend(loc='lower center', shadow=True, fancybox=True)


        y1 = int(min(VTy)-0.5)
        y2 = int(max(VTy)+0.5)

        mean_val = nmp.mean(VDy)
        df = max( abs(min(VDm)-mean_val), abs(max(VDm)-mean_val) )


        y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)

        if ymin==0 and ymax==0:
            plt.axis( [y1, y2, min(VDm)-0.2*df, max(VDm)+0.2*df] )

        else:
            plt.axis([y1, y2, ymin, ymax])

            if dy != 0:
                plt.yticks( nmp.arange(float(int(ymin+0.5)), float(int(ymax))+dy, dy) )
                locs, labels = plt.yticks() #lolo?
                if i_y_jump > 1: __subsample_axis__(plt, 'y', i_y_jump)

        __nice_x_axis__(ax, plt, y1, y2, dt_year, cfont=font_xylb)

        ax.grid(color='k', linestyle='-', linewidth=0.2)

        plt.ylabel('('+cyunit+')', **font_xylb)

        plt.title(ctitle, **font_ttl)

        cf_fig = cfignm+'.'+cfig_type

        plt.savefig(cf_fig, dpi=DPI_DEF, orientation='portrait', transparent=l_tranparent_bg) ; #1d_mon_ann

        plt.close(1)








    def __1d_multi(self,vt, XD, vlabels, cfignm='fig', dt_year=5, i_t_jump=1, cyunit='', ctitle='',
                   cfig_type='png', ymin=0, ymax=0, lzonal=False, xmin=0, xmax=0,
                   loc_legend='lower center', line_styles=[], fig_size=FIG_SIZE_DEF,
                   l_tranparent_bg=True, cxunit='', lmask=True):

        # lzonal => zonally averaged curves...


        if lzonal:
            font_ttl, font_big_fixed, font_xylb, font_clb = __font_unity__(size='big')
        else:
            font_ttl, font_xylb, font_clb = __font_unity__()

        # Number of lines to plot:
        [ nb_plt, nbt ] = XD.shape

        if len(vt) != nbt: print 'ERROR: plot_1d_multi.barakuda_plot.py => vt and XD do not agree in shape!'; sys.exit(0)
        if len(vlabels) != nb_plt: print 'ERROR: plot_1d_multi.barakuda_plot.py => wrong number of labels...'; sys.exit(0)

        n0 = len(line_styles)
        if n0 > 0 and n0 != nb_plt: print 'ERROR: plot_1d_multi.barakuda_plot.py => wrong number line styles!!!'; sys.exit(0)



        # Masking the time-series shorter than others (masked with -999.)
        if lmask: XD = nmp.ma.masked_where(XD < -900., XD)

        if lzonal:
            fig = plt.figure(num = 1, figsize=fig_size, facecolor='w', edgecolor='k')
            ax = plt.axes([0.08, 0.11, 0.88, 0.83])
        else:
            fig = plt.figure(num = 1, figsize=fig_size, facecolor='w', edgecolor='k')
            ax = plt.axes(AXES_DEF) ; #1d_multi

        if lzonal: plt.plot(vt[:], XD[0,:]*0., 'k', linewidth=1)

        for jp in range(nb_plt):
            if n0 > 0:
                plt.plot(vt[:], XD[jp,:], line_styles[jp],   label=vlabels[jp], linewidth=2)
            else:
                plt.plot(vt[:], XD[jp,:], v_dflt_colors[jp], label=vlabels[jp], linewidth=2)


        if loc_legend != '0':  plt.legend(loc=loc_legend, ncol=int(nb_plt/4+1), shadow=True, fancybox=True)

        # Prevents from using scientific notations in axess ticks numbering:
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)

        if lzonal:
            dt_year = 15. ; # x-axis increment (latitude!)
            if xmin == 0 and xmax == 0:
                y1 = -90. ; y2 = 90.
            else:
                y1 = xmin ;  y2 = xmax
        else:
            if xmin == 0 and xmax == 0:
                y1 = int(vt[0])
                y2 = int(round(vt[len(vt)-1]+0.4))
            else:
                y1 = xmin ; y2 = xmax

        if ymin==0 and ymax==0:
            mean_val = nmp.mean(XD[:,:])
            df = max( abs(nmp.min(XD[:,:])-mean_val), abs(nmp.max(XD[:,:])-mean_val) )
            plt.axis( [y1, y2, nmp.min(XD[:,:])-0.2*df, nmp.max(XD[:,:])+0.2*df] )
        else:
            plt.axis([y1, y2, ymin,     ymax])

        if lzonal:
            __nice_x_axis__(ax, plt, y1, y2, 10., cfont=font_xylb)
        else:
            __nice_x_axis__(ax, plt, y1, y2, dt_year, i_sbsmp=i_t_jump, cfont=font_xylb)

        ax.grid(color='k', linestyle='-', linewidth=0.2)

        if lzonal: plt.xlabel(r'Latitude ($^{\circ}$N)', **font_xylb)

        if cyunit != '': plt.ylabel('('+cyunit+')', **font_xylb)
        if cxunit != '': plt.xlabel('('+cxunit+')', **font_xylb)

        plt.title(ctitle, **font_ttl)

        cf_fig = cfignm+'.'+cfig_type

        plt.savefig(cf_fig, dpi=DPI_DEF, orientation='portrait', transparent=l_tranparent_bg) ; #1d_multi

        plt.close(1)
        print '   => Multi figure "'+cf_fig+'" created!'










    def __1d(self,vt, VF, cfignm='fig', dt_year=5, i_t_jump=1, cyunit='', ctitle='',
                cfig_type='png', ymin=0, ymax=0, xmin=0, xmax=0,
                loc_legend='lower center', line_styles='-', fig_size=FIG_SIZE_DEF,
                l_tranparent_bg=False, cxunit='', lmask=True):

        font_ttl, font_xylb, font_clb = __font_unity__()

        # Number of lines to plot:
        nbt = len(VF)

        if len(vt) != nbt: print 'ERROR: plot_1d.barakuda_plot.py => vt and VF do not agree in shape!'; sys.exit(0)


        # Masking the time-series shorter than others (masked with -999.)
        if lmask: VF = nmp.ma.masked_where(VF < -900., VF)

        fig = plt.figure(num = 1, figsize=fig_size, facecolor='w', edgecolor='k')
        ax = plt.axes(AXES_DEF) ; #1d

        plt.plot(vt[:], VF[:], line_styles, linewidth=2)


        # Prevents from using scientific notations in axess ticks numbering:
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)

        if xmin == 0 and xmax == 0:
            y1 = int(vt[0])
            y2 = int(round(vt[len(vt)-1]+0.4))
        else:
            y1 = xmin ; y2 = xmax

        if ymin==0 and ymax==0:
            mean_val = nmp.mean(VF[:])
            df = max( abs(nmp.min(VF[:])-mean_val), abs(nmp.max(VF[:])-mean_val) )
            plt.axis( [y1, y2, nmp.min(VF[:])-0.2*df, nmp.max(VF[:])+0.2*df] )
        else:
            plt.axis([y1, y2, ymin,     ymax])

        print nmp.arange(y1, y2+dt_year, dt_year)

        __nice_x_axis__(ax, plt, y1, y2, dt_year, i_sbsmp=i_t_jump, cfont=font_xylb)

        ax.grid(color='k', linestyle='-', linewidth=0.2)

        if cyunit != '': plt.ylabel('('+cyunit+')', **font_xylb)
        if cxunit != '': plt.xlabel('('+cxunit+')', **font_xylb)

        plt.title(ctitle, **font_ttl)

        cf_fig = cfignm+'.'+cfig_type

        plt.savefig(cf_fig, dpi=DPI_DEF, orientation='portrait', transparent=l_tranparent_bg) ; #1d

        plt.close(1)
        print '   => Multi figure "'+cf_fig+'" created!'











    def __spectrum(self,vfrq, Vspec, cfignm='fig', cyunit='', log_x=True, log_y=False,
                      year_min=3., year_max = 50., rmax_amp = 10., rmin_amp = 0.,
                      cfig_type='png', vnoise=[ 0 ], vrci95=[ 0 ], lab95_xpos=0.5, lplot_1onF=False,
                      cnoise='White noise', lplot_freq_ax=True):

        l_do_ci95 = False ; l_do_ci95m = False

        nnoise = len(vnoise); nl95 = len(vrci95)

        if nnoise != 1 and nl95 != 1:
            if nl95 != len(Vspec) or nl95 != nnoise:
                print "ERROR: plot_spectrum.barakuda_plot.py => length of 95 CI array and/or noise array doesnt match spectrum length!"
                sys.exit(0)
            l_do_ci95 = True
            l_do_ci95m = True

        font_ttl, font_xylb, font_clb = __font_unity__()


        print "avant:", rmin_amp, rmax_amp
        if log_y:
            if rmin_amp <= 0.: rmin_amp = 0.01
            rmin_amp = 20.*nmp.log10(rmin_amp); rmax_amp = 20.*nmp.log10(rmax_amp)
        print "apres:", rmin_amp, rmax_amp


        # Spectral axis:
        x_min = 1./year_max ; x_max = 1./year_min ; # min and max in frequency!
        clbnd = str(int(round(year_min)))

        if log_x:
            cvee = [ '50','45','40','35','30','25','22','20','17','15','13','12','11','10','9','8','7','6','5','4','3' ]
            if year_max == 40.: cvee = cvee[2:]
            if year_max == 35.: cvee = cvee[3:]
            if year_min ==  5.: cvee = cvee[:-2]
        else:
            cvee = [ '50','30','20','15','12','10','9','8','7','6','5','4', '3' ]


        lvee = []
        civee = []
        for ce in cvee:
            lvee.append(float(ce))
            civee.append(str(round(1./float(ce),3)))
        vee  = nmp.asarray(lvee)

        print civee[:]


        rnoise = nmp.mean(vnoise[5:20])
        rrci95 = nmp.mean(vrci95[5:20])


        fig = plt.figure(num = 1, figsize=(8.,4.), facecolor='w', edgecolor='k')

        if lplot_freq_ax:
            ax = plt.axes([0.069, 0.13, 0.9, 0.8])
        else:
            ax = plt.axes([0.08, 0.13, 0.9, 0.83])

        if log_x:
            ifr1 = 1
            vfl = nmp.log10(vfrq[ifr1:])

            if not log_y:
                if l_do_ci95:
                    plt.plot(vfl, vnoise[ifr1:],               '--k'  , linewidth=1.8, label=cnoise)
                    plt.plot(vfl, vnoise[ifr1:]+vrci95[ifr1:], '0.4', linewidth=1.8, label='95% CI')
                    plt.plot(vfl, vnoise[ifr1:]-vrci95[ifr1:], '0.4', linewidth=1.8)
                plt.plot(vfl, Vspec[ifr1:],   '*-k', linewidth=2.)
                #if lplot_1onF: plt.plot(vfl, 1./vfl,   b_red, linewidth=2)
            else:
                if l_do_ci95:
                    plt.plot(vfl, 20.*nmp.log10(vnoise[ifr1:]),               '--k'  , linewidth=1.8, label=cnoise)
                    plt.plot(vfl, 20.*nmp.log10(vnoise[ifr1:]+vrci95[ifr1:]), '0.4', linewidth=1.8, label='95% CI')
                    plt.plot(vfl, 20.*nmp.log10(vnoise[ifr1:]-vrci95[ifr1:]), '0.4', linewidth=1.8)
                plt.plot(vfl, 20.*nmp.log10(Vspec[ifr1:]),   '*-k', linewidth=2.)

        else:
            if not log_y:
                if l_do_ci95:
                    plt.plot(vfrq, vnoise,        '--k'  , linewidth=1.8)
                    plt.plot(vfrq, vnoise+vrci95, '0.4', linewidth=1.8)
                    plt.plot(vfrq, vnoise-vrci95, '0.4', linewidth=1.8)
                plt.plot(vfrq, Vspec,   '*-k', linewidth=2)
                if lplot_1onF: plt.plot(vfrq[1:], 0.03*1./vfrq[1:],   b_red, linewidth=2)
            else:
                if l_do_ci95:
                    plt.plot(vfrq, 20.*nmp.log10(vnoise),        '--k'  , linewidth=1.8)
                    plt.plot(vfrq, 20.*nmp.log10(vnoise+vrci95), '0.4', linewidth=1.8)
                    plt.plot(vfrq, 20.*nmp.log10(vnoise-vrci95), '0.4', linewidth=1.8)
                plt.plot(vfrq, 20.*nmp.log10(Vspec),   '*-k', linewidth=2)

        plt.ylabel('Amplitude Spectrum ('+cyunit+')', color='k', **font_xylb)

        plt.xlabel('Period (years)', color='k', **font_xylb)

        if log_x:
            x1=nmp.log10(x_max) ; x2=nmp.log10(x_min)
            plt.axis([x1, x2, rmin_amp, rmax_amp])
            if lplot_freq_ax:
                plt.xticks(nmp.log10(1./vee[:]),cvee[:])
            else:
                print ''
                vee_n = nmp.arange(vee[0], vee[len(vee)-1]-1, -1.)
                print vee_n[:]
                cvee_n = []
                for rr in vee_n:
                    cr = str(int(rr))
                    if cr in cvee:
                        cvee_n.append(cr)
                    else:
                        cvee_n.append('')
                print 'cvee =>', cvee[:]
                print 'cvee_n =>', cvee_n[:]
                plt.xticks(nmp.log10(1./vee_n[:]),cvee_n[:])
        else:
            x1=x_max; x2=x_min
            plt.axis([x1, x2, rmin_amp, rmax_amp])
            plt.xticks(1./vee[:],cvee[:], color='k')

        ax.grid(color='0.4', linestyle='-', linewidth=0.3)

        plt.legend(loc='upper left', shadow=False, fancybox=True)

        if lplot_freq_ax:
            ax2 = ax.twiny()
            ax2.set_xlabel('Frequency (cpy)', color='k', **font_xylb)
            if log_x:
                plt.axis([x1, x2, rmin_amp, rmax_amp])
                for jp in range(1,18,2): civee[jp] = ''
                plt.xticks(nmp.log10(1./vee[:]),civee)
                #ax2.xaxis.set_ticks(nmp.log10(1./vee[:]))
            else:
                plt.axis([x1, x2, rmin_amp, rmax_amp])
                plt.xticks(1./vee[:],civee)

            for t in ax2.get_xticklabels(): t.set_fontsize(14)
            ax2.xaxis.labelpad = 12 ; # move label upwards a bit...

        plt.savefig(cfignm+'.'+cfig_type, dpi=100, facecolor='w', edgecolor='w', orientation='portrait', transparent=False)
        plt.close(1)

        def __del__(self) :

            plot.__counter -=  1





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





def __font_unity__(size='normal'):


    rat = 1.
    if size == 'big': rat = 1.25

    params = { 'font.family':'Trebuchet MS',
               'font.size':       int(14.*rat),
               'legend.fontsize': int(14.*rat),
               'xtick.labelsize': int(14.*rat),
               'ytick.labelsize': int(14.*rat),
               'axes.labelsize':  int(14.*rat) }

    mpl.rcParams.update(params)

    title_fonts    = { 'fontname':'Trebuchet MS', 'fontweight':'normal', 'fontsize':int(16.*rat) }
    label_fonts    = { 'fontname':'Arial'       , 'fontweight':'normal', 'fontsize':int(15.*rat) }
    colorbar_fonts = { 'fontname':'Arial'       , 'fontweight':'normal', 'fontsize':int(14.*rat) }

    return title_fonts, label_fonts, colorbar_fonts




def __force_min_and_max__(rm, rp, Xin):
    idx_bad  = nmp.where(nmp.logical_not(nmp.isfinite(Xin)))
    Xin[idx_bad] = 0.
    idx1 = nmp.where(Xin <= rm); Xin[idx1] = rm + abs(rp-rm)*1.E-4
    idx2 = nmp.where(Xin >= rp); Xin[idx2] = rp - abs(rp-rm)*1.E-4
    Xin[idx_bad] = nmp.nan


def __subsample_colorbar__(i_sbsmp, vcc, clb_hndl, cb_or='vertical'):
    cb_labs = []

    # First checking if vcc countains integers or not...
    lcint = False
    vc = vcc.astype(nmp.int64)  ; # integer version of vcc
    if nmp.max(nmp.abs(vcc))>5 and nmp.sum(vcc-vc) == 0. : lcint=True

    cpt = 0
    if int(vcc[0]) % 2 != 0: cpt = 1   # not an even number !

    for rr in vcc:
        if cpt % i_sbsmp == 0:
            if lcint:
                cr = str(int(rr))
            else:
                cr = str(round(float(rr),6))
            cb_labs.append(cr)
        else:
            cb_labs.append(' ')
        cpt = cpt + 1
    if cb_or == 'horizontal':
        clb_hndl.ax.set_xticklabels(cb_labs)
    else:
        clb_hndl.ax.set_yticklabels(cb_labs)
    del cb_labs, vc



def __nice_colorbar__(fig_hndl, plt_hndl, vcc,
                      cax_other=None, i_sbsmp=1, lkc=False, cb_or='vertical', cunit=None, cfont=None, fontsize=0):

    if cb_or not in {'horizontal','vertical'}:
        print "ERROR: only 'vertical' or 'horizontal' can be specified for the colorbar orientation!"
        cb_or = 'vertical'

    if cb_or == 'horizontal':
        if cax_other is not None:
            clb = plt_hndl.colorbar(fig_hndl, cax=cax_other, ticks=vcc, drawedges=lkc, orientation='horizontal', pad=0.07, shrink=1., aspect=40)
        else:
            clb = plt_hndl.colorbar(fig_hndl,                ticks=vcc, drawedges=lkc, orientation='horizontal', pad=0.07, shrink=1., aspect=40)
    else:
        if cax_other is not None:
            clb = plt_hndl.colorbar(fig_hndl, cax=cax_other, ticks=vcc, drawedges=lkc, pad=0.03)
        else:
            clb = plt_hndl.colorbar(fig_hndl,                ticks=vcc, drawedges=lkc, pad=0.03)

    if i_sbsmp > 1: __subsample_colorbar__(i_sbsmp, vcc, clb, cb_or=cb_or)

    if not cunit is None:
        if cfont is None:
            clb.set_label(cunit)
        else:
            clb.set_label(cunit, **cfont)

    if fontsize > 0:
        if cb_or == 'horizontal':
            for t in clb.ax.get_xticklabels(): t.set_fontsize(fontsize) # Font size for colorbar ticks!
        else:
            for t in clb.ax.get_yticklabels(): t.set_fontsize(fontsize) # Font size for colorbar ticks!



def _add_xy_offset__(plt_hndl, ixo, iyo):
    if ( ixo != 0. ):
        locs, labels = plt_hndl.xticks() ; jl=0
        vlabs = []
        for ll in locs:
            clab = str(int(locs[jl])+int(ixo))
            vlabs.append(clab); jl=jl+1
        plt_hndl.xticks(locs,vlabs)
    if ( y_offset != 0. ):
        locs, labels = plt_hndl.yticks() ; jl=0; vlabs = []
        for ll in locs:
            clab = str(int(locs[jl])+int(iyo))
            vlabs.append(clab); jl=jl+1
        plt_hndl.yticks(locs,vlabs)
    del vlabs

def __subsample_axis__(plt_hndl, cax, i_sbsmp, icpt=1):
    ax_lab = []
    if   cax == 'x':
        locs, labels = plt_hndl.xticks()
    elif cax == 'y':
        locs, labels = plt_hndl.yticks()
    else:
        print ' Error: __subsample_axis__.barakuda_plot => only "x" or "y" please'; sys.exit(0)
    cpt = icpt # with ipct = 1: tick priting will start at y1+dt_year on x axis rather than y1
    for rr in locs:
        if cpt % i_sbsmp == 0:
            if rr%1.0 == 0.:
                cr = str(int(rr))  # it's something like 22.0, we want 22 !!!
            else:
                cr = str(rr)
            ax_lab.append(cr)
        else:
            ax_lab.append(' ')
        cpt = cpt + 1
    if cax == 'x': plt_hndl.xticks(locs,ax_lab)
    if cax == 'y': plt_hndl.yticks(locs,ax_lab)
    del ax_lab



def __nice_x_axis__(ax_hndl, plt_hndl, x_0, x_L, dx, i_sbsmp=1, cunit=None, cfont=None):
    plt_hndl.xticks( nmp.arange(x_0, x_L+dx, dx) )
    locs, labels = plt_hndl.xticks()
    #jl=0
    #xlabs = []
    #for ll in locs:
    #    xlabs.append(str(int(locs[jl])))
    #    jl=jl+1
    #plt_hndl.xticks(locs,xlabs)

    if i_sbsmp > 1: __subsample_axis__( plt, 'x', i_sbsmp)

    ax_hndl.set_xlim(x_0,x_L)

    if not cunit is None:
        if cfont is None:
            plt_hndl.xlabel(cunit)
        else:
            plt_hndl.xlabel(cunit, **cfont)

    # Prevents from using scientific notations in axess ticks numbering:
    #ax_hndl.get_xaxis().get_major_formatter().set_useOffset(False)

    #del xlabs


def __nice_z_axis__(ax_hndl, plt_hndl, z0, zK, dz, i_sbsmp=1, cunit=None, cfont=None):
    iia = 1
    if zK <= 0. and z0 < 0.:
        iia = 0
    vv = nmp.arange(z0, zK+dz, dz)
    plt_hndl.yticks(vv)
    if i_sbsmp > 1: __subsample_axis__(plt, 'y', i_sbsmp, icpt=0)
    ax_hndl.set_ylim(z0-dz/2., zK+dz/2.)
    if not cunit is None:
        if cfont is None:
            plt_hndl.ylabel(cunit)
        else:
            plt_hndl.ylabel(cunit, **cfont)








def __prepare_z_log_axis__(l_log, vz):

    import math

    nk  = len(vz)
    zvz = nmp.zeros(nk)

    if l_log:
        for jk in range(nk):
            zvz[jk] = math.log10(vz[jk])
    else:
        zvz= vz

    return zvz


def __fix_z_axis__(ax_hndl, plt_hndl, z0, zK, l_log=False, l_z_inc=True):

    if l_log:
        y_log_ofs = 10.
        vyview_list = [ 3. , 10. , 25., 50. , 100. , 250. , 500. , 1000. , 2500.,  5000. ]
        nd = len(vyview_list)
        vyview      = nmp.zeros(nd)
        for jn in range(nd): vyview[jn] = vyview_list[jn]
        vyview_log = nmp.log10(vyview + y_log_ofs)
        ylab = []
        for rr in vyview_list: ylab.append(str(int(rr)))
        z0 = nmp.log10(z0+y_log_ofs)
        zK = nmp.log10(zK+y_log_ofs)
        ax_hndl.set_yticks(vyview_log)
        ax_hndl.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        ax_hndl.set_yticklabels(ylab)

    if l_z_inc:
        ax_hndl.set_ylim(z0,zK)
    else:
        ax_hndl.set_ylim(zK+(zK-z0)/50. , z0)

    ax_hndl.grid(color='k', linestyle='-', linewidth=0.5)
