# My own colormaps...

# Check http://matplotlib.org/examples/color/colormaps_reference.html !!!


# Colors:  http://www.pitt.edu/~nisg/cis/web/cgi/rgb.html

# Last Updated: L. Brodeau, January 2013

import sys
import numpy as nmp
import matplotlib
matplotlib.use('Agg') ; # important so no DISPLAY is needed!!!
from matplotlib.pylab import cm


ctrl = 0.2 ; # for logarythmic scale in 'pal_eke'



class pal :
    '''This class encapsulates all the pal routines                                                                                                                                                              
    
    In order to use it you need to type as follows:                                                                                                                                                  
    
    pal(pal name without the prefix __) (all the arguments of the function)
    
    for example for __eke we do as follows:
    
    pal("eke")()
        
    Functions of the pal class and they can not be accessed directly outside of the class. That is, you need to call these functios through the call wrapper                                                                  
    as seen below in the function __call__.                                                                                                                                                                        
    '''

    def __init__(self,spal) :

        self.spal = spal


    def __call__(self,*args, **kw) :

        if "_"+self.__class__.__name__+ "__" + self.spal in self.__class__.__dict__.keys() :

            self.__class__.__dict__["_"+self.__class__.__name__+ "__" + self.spal](self,*args, **kw)

        else :
            print "pal " + "__" + self.spal + " does not exist"
            sys.exit()


    def _blk(self):
        M = nmp.array( [
            [ 0. , 0., 0. ], # black
            [ 0. , 0., 0. ]  # black
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap




    # ['rgb(255,255,204)','rgb(161,218,180)','rgb(65,182,196)','rgb(44,127,184)','rgb(37,52,148)']

    
    # http://colorbrewer2.org/?type=sequential&scheme=YlGnBu&n=5
    def __cb1(self):
        M = nmp.array( [
            [ 255,255,204  ],
            [ 161,218,180  ],
            [ 65,182,196  ],
            [ 44,127,184  ],
            [ 37,52,148  ]
        ] ) / 255.
        my_cmap = __build_colormap__(M, log_ctrl=ctrl)
        return my_cmap



        
    def __eke(self):
        M = nmp.array( [
            [ 0.  , 0.0 , 0.2  ], # black
            [ 0.1 , 0.5 , 1.0  ], # blue
            [ 0.2 , 1.0 , 0.0  ], # green
            [ 1.  , 1.0 , 0.0  ], # yellow
            [ 1.  , 0.0 , 0.0  ], # red
            [0.2  , 0.27, 0.07 ] # brown
        ] )
        my_cmap = __build_colormap__(M, log_ctrl=ctrl)
        return my_cmap


    def __bathy(self):
        M = nmp.array( [
            [ 0.0 , 0.0 , 0.4 ], # dark blue
            [ 0.1 , 0.5 , 1.0  ], # blue
            [ 0.2 , 1.0 , 0.0  ], # green
            [ 1.  , 1.0 , 0.0  ], # yellow
            [ 1.  , 0.0 , 0.0  ], # red
            [0.2  , 0.27, 0.07 ] # brown
        ] )
        my_cmap = __build_colormap__(M, log_ctrl=ctrl)
        return my_cmap


    #        [ 0.4 , 0.0 , 0.6 ], # violet
    def __mld(self):
        M = nmp.array( [
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 0.1 , 0.5 , 1.0 ], # light blue
            [ 0.13, 0.54, 0.13], # dark green
            [ 0.2 , 1.0 , 0.0 ], # light green
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 0.2 , 0.3 , 0.1 ] # dark redish brown
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    def __mld_r(self):
        M = nmp.array( [
            [ 0.4 , 0.0 , 0.6 ], # violet
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 0.1 , 0.5 , 1.0 ], # light blue
            [ 0.13, 0.54, 0.13], # dark green
            [ 0.2 , 1.0 , 0.0 ], # light green
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 0.2 , 0.3 , 0.1 ] # dark redish brown
        ] )
        my_cmap = __build_colormap__(M[::-1,:])
        return my_cmap



    def __jetblanc(self):
        M = nmp.array( [
            [ 0.6 , 0.0 , 0.8 ], # violet
            [ 0.0 , 0.0 , 0.4 ], # dark blue
            [ 0.1 , 0.5 , 1.0 ], # light blue
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 0.2 , 0.3 , 0.1 ]  # dark redish brown
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap


    def __jetblanc_r(self):
        M = nmp.array( [
            [ 0.6 , 0.0 , 0.8 ], # violet
            [ 0.0 , 0.0 , 0.4 ], # dark blue
            [ 0.1 , 0.5 , 1.0 ], # light blue
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 0.2 , 0.3 , 0.1 ] # dark redish brown
        ] )
        my_cmap = __build_colormap__(M[::-1,:])
        return my_cmap



    def __amoc(self):
        M = nmp.array( [
            [ 0.4 , 0.0 , 0.6 ], # violet
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 1.0 , 1.0 , 1.0 ], # white
            [0.68 , 0.98, 0.98], # light blue
            [ 0.0 , 0.0 , 0.95], # dark blue
            [ 0.2 , 1.0 , 0.0 ], # green
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 0.2 , 0.3 , 0.1 ] # dark read
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap


    #        [ 0.2 , 1.0 , 0.0 ], # green

    def __sst(self):
        M = nmp.array( [
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 0.4 , 0.0 , 0.6 ], # violet
            [ 0. , 0.2 , 0.99], # dark blue
            [0.68 , 0.98, 0.98], # light blue
            [ 0.13, 0.54, 0.13], # dark green
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 0.2 , 0.3 , 0.1 ] # dark read
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap


    def __sst_r(self):
        M = nmp.array( [
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 0.4 , 0.0 , 0.6 ], # violet
            [ 0. , 0.2 , 0.99], # dark blue
            [0.68 , 0.98, 0.98], # light blue
            [ 0.13, 0.54, 0.13], # dark green
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 0.2 , 0.3 , 0.1 ] # dark read
        ] )
        my_cmap = __build_colormap__(M[::-1,:])
        return my_cmap


    def __sst0(self):
        M = nmp.array( [
        [ 1.0 , 1.0 , 1.0 ], # white
        [ 0.4 , 0.0 , 0.6 ], # violet
        [ 0.0 , 0.0 , 0.95], # dark blue
        [0.68 , 0.98, 0.98], # light blue
        [46./255., 203./255., 35./255.], # green
        [ 1.0 , 1.0 , 0.0 ], # yellow
        [ 1.0 , 0.0 , 0.0 ], # red
        [ 0.2 , 0.3 , 0.1 ] # dark read
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    #        [ 253./255., 238./255., 1. ], # really pale pink
    #        [ 0.2 , 1.0 , 0.0 ], # green

    def __sst0_r(self):
        M = nmp.array( [
        [ 1.0 , 1.0 , 1.0 ], # white
        [ 0.4 , 0.0 , 0.6 ], # violet
        [ 0.0 , 0.0 , 0.95], # dark blue
        [0.68 , 0.98, 0.98], # light blue
        [46./255., 203./255., 35./255.], # green
        [ 1.0 , 1.0 , 0.0 ], # yellow
        [ 1.0 , 0.0 , 0.0 ], # red
        [ 0.2 , 0.3 , 0.1 ] # dark read
        ] )
        my_cmap = __build_colormap__(M[::-1,:])
        return my_cmap





    def __std(self):
        M = nmp.array( [
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 0.8 , 0.8 , 0.8 ], # grey
            [ 0.4 , 0.0 , 0.6 ], # violet
            [ 0.0 , 0.0 , 0.95], # dark blue
            [0.68 , 0.98, 0.98], # light blue
            [ 0.13, 0.54, 0.13], # dark green
            [ 0.2 , 1.0 , 0.0 ], # green
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 0.2 , 0.3 , 0.1 ] # dark read
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    #        [  0. , 0. ,  1.0 ], # blue
    def __ice(self):
        M = nmp.array( [
            [  0. , 0.  , 0.3 ], # dark blue
            [ 0.6 , 0.6 , 0.8 ], # light grey
            [ 0.95 , 0.95 , 0.95 ],  # white
            [ 1.0 , 1.0 , 1.0 ]  # white
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    def __blanc(self):
        M = nmp.array( [
            [ 1.0 , 1.0 , 1.0 ],  # white
            [ 1.0 , 1.0 , 1.0 ]  # white
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap


    def __rms(self):
        M = nmp.array( [
            [ 1.0 , 1.0 , 1.0 ],
            [ 0.1 , 0.5 , 1.0 ],
            [ 0.2 , 1.0 , 0.0 ],
            [ 1.0 , 1.0 , 0.0 ],
            [ 1.0 , 0.0 , 0.0 ],
            [ 0.2 , 0.3 , 0.1 ]
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap



    def __sigtr(self):
        #[ 1.0 , 0.4 , 1.0 ], # violet pinkish
        #[ 1.0 , 1.0 , 1.0 ], # white
        M = nmp.array( [
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 0.0 , 0.8 , 1.0 ], #light blue
            [ 0.1 , 0.5 , 1.0 ], #light blue
            [ 0.0 , 0.0 , 0.4 ], # blue
            [ 0.0 , 0.4 , 0.0 ], # dark green
            [ 0.1 , 1.0 , 0.0 ], # green
            [ 0.4 , 1.0 , 0.0 ], # vert pomme
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 1.0 , 0.4 , 0.0 ], # orange
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 0.6 , 0.0 , 0.0 ], # red
            [ 0.2 , 0.3 , 0.1 ]  # dark red
        ] )
        my_cmap = __build_colormap__(M[::-1,:])
        return my_cmap


    def __sigtr_r(self):
        #[ 1.0 , 0.4 , 1.0 ], # violet pinkish
        #[ 1.0 , 1.0 , 1.0 ], # white
        M = nmp.array( [
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 0.0 , 0.8 , 1.0 ], #light blue
            [ 0.1 , 0.5 , 1.0 ], #light blue
            [ 0.0 , 0.0 , 0.4 ], # blue
            [ 0.0 , 0.4 , 0.0 ], # dark green
            [ 0.1 , 1.0 , 0.0 ], # green
            [ 0.4 , 1.0 , 0.0 ], # vert pomme
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 1.0 , 0.4 , 0.0 ], # orange
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 0.6 , 0.0 , 0.0 ], # red
            [ 0.2 , 0.3 , 0.1 ]  # dark red
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap


    def __bbr(self):
        M = nmp.array( [
            [ 0.  , 0. , 0.2 ],
            [ 0.  , 0. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 0. , 0.  ],
            [ 0.6 , 0. , 0.  ]
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    def __bbr_r(self):
        M = nmp.array( [
            [ 0.  , 0. , 0.2 ],
            [ 0.  , 0. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 0. , 0.  ],
            [ 0.6 , 0. , 0.  ]
        ] )
        my_cmap = __build_colormap__(M[::-1,:])
        return my_cmap




    def __bbr2(self):
        M = nmp.array( [
            [ 0.  , 1. , 1.  ], # cyan
            [ 0.  , 0. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 0. , 0.  ],
            [ 1.  , 1. , 0.  ]  # jaune
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    def __bbr2_r(self):
        M = nmp.array( [
            [ 0.  , 1. , 1.  ], # cyan
            [ 0.  , 0. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 0. , 0.  ],
            [ 1.  , 1. , 0.  ]  # jaune
        ] )
        my_cmap = __build_colormap__(M[::-1,:])
        return my_cmap


    def __bbr0(self):
        M = nmp.array( [
            [ 0.  , 1. , 1.  ], # cyan
            [ 0.  , 0. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 0. , 0.  ],
            [ 1.  , 1. , 0.  ]  # jaune
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap



    def __bbr_cold(self):
        M = nmp.array( [
            [ 0.  , 1. , 1.  ], # cyan
            [ 0.  , 0. , 1.  ],
            [ 19./255.  , 7./255. , 129./255  ], # dark blue
            [ .1  , .1 , .9  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 0.7  , 0. , 0.  ] # Dark red        
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap


    # light red:         [ 1.  , 0.3 , 0.3  ]

    def __bbr_warm(self):
        M = nmp.array( [
            #        [ 0.3  , 0.3 , 1.  ],
            #        [ 0.6  , 0.6 , 1.  ],
            [ 19./255.  , 7./255. , 129./255  ], # dark blue
            [ 1.  , 1. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 0.9  , 0.1 , 0.1  ],
            [ 0.7  , 0. , 0.  ], # Dark red
            #[ 1.  , 0.2 , 0.  ],
            [ 1.  , 1. , 0.  ],  # jaune
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    #        [ 19./255.  , 7./255. , 129./255  ], # dark blue
    #         [ 0.  , 1. , 1.  ], # cyan
    #


    def __bbr0_r(self):
        M = nmp.array( [
            [ 0.  , 1. , 1.  ], # cyan
            [ 0.  , 0. , 1.  ],
            [ 1.  , 1. , 1.  ],
            [ 1.  , 0. , 0.  ],
            [ 1.  , 1. , 0.  ]  # jaune
        ] )
        my_cmap = __build_colormap__(M[::-1,:])
        return my_cmap

    def __cold0(self):
        M = nmp.array( [
            [ 177./255.  , 250./255. , 122./255. ],   # greenish
            [ 0.  , 1. , 1.  ], # cyan
            [ 7./255.  , 11./255. , 122./255. ], # dark blue
            [ 0.  , 0. , 1.  ], # true blue
            [ 177./255.  , 189./255. , 250./255. ], # light blue
            [ 1.  , 1. , 1.  ],
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap
    # #        [ 177./255.  , 250./255. , 122./255. ],   # greenish



    def __warm0(self):
        M = nmp.array( [
            [ 1.  , 1. , 1.  ],
            [ 255./255.  , 254./255. , 198./255.  ], # very light yellow
            [ 1.  , 1. , 0.  ],  # yellow
            [ 244./255.  , 78./255. , 255./255.  ], # pink
            [ 1.  , 0. , 0.  ], # true red
            [ 139./255.  , 5./255. , 5./255.  ] # dark red
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap
    #        [ 247./255.  , 150./255. , 176./255.  ], # light red

    def __graylb(self):
        M = nmp.array( [
            [ 1.  , 1. , 1. ],
            [ 0.1  , 0.1 , 0.1 ]
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    def __graylb_r(self):
        M = nmp.array( [
            [ 0.  , 0. , 0. ],
            [ 1.  , 1. , 1. ]
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    def __graylb2(self):
        M = nmp.array( [
            [ 0.6  , 0.6 , 0.6 ],
            [ 1.  , 1. , 1. ]
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap



    def __sigma(self):
        M = nmp.array( [
            [ 1.0 , 1.0 , 1.0 ], # white
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 0.2 , 1.0 , 0.0 ], # green
            [ 0.1 , 0.5 , 1.0 ], # light blue
            [ 0.0 , 0.0 , 0.4 ], # dark blue
            [ 0.6 , 0.0 , 0.8 ]  # violet
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    def __sigma0(self):
        M = nmp.array( [
            [ 0.2 , 0.3 , 0.1 ], # dark redish brown
            [ 1.0 , 0.0 , 0.0 ], # red
            [ 1.0 , 1.0 , 0.0 ], # yellow
            [ 0.2 , 1.0 , 0.0 ], # green
            [ 0.1 , 0.5 , 1.0 ], # light blue
            [ 0.0 , 0.0 , 0.4 ], # dark blue
            [ 0.6 , 0.0 , 0.8 ], # violet
            [ 1.0 , 1.0 , 1.0 ]  # white
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    def __mask(self):
        M = nmp.array( [
            [ 0.5 , 0.5 , 0.5 ], # gray
            [ 0.5 , 0.5 , 0.5 ]  # gray
        ] )
        my_cmap = __build_colormap__(M)
        return my_cmap

    def __jet(self):
        
        return cm.jet


# ===== local functions ======
def __build_colormap__(MC, log_ctrl=0):

    [ nc, n3 ] = nmp.shape(MC)

    # Make x vector :
    x =[]
    for i in range(nc): x.append(255.*float(i)/((nc-1)*255.0))
    x = nmp.array(x)
    if log_ctrl > 0: x = nmp.log(x + log_ctrl)
    rr = x[nc-1] ; x  = x/rr

    y =nmp.zeros(nc)
    for i in range(nc): y[i] = x[nc-1-i]

    x = 1 - y ; rr = x[nc-1] ; x  = x/rr

    red  = [] ; blue = [] ; green = []

    for i in range(nc):
        red.append  ([x[i],MC[i,0],MC[i,0]])
        green.append([x[i],MC[i,1],MC[i,1]])
        blue.append ([x[i],MC[i,2],MC[i,2]])

    cdict = {'red':red, 'green':green, 'blue':blue}
    my_cm = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

    return my_cm
