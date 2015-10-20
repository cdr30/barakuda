#
# Laurent Brodeau
# 
# November 2013

import sys
import numpy as nmp
import math




########################
# SIGMA DENSITY STUFFS #
########################

rsigma_dense = 27.8

zrau0 = 1000.
rt0   = 273.16
grav  = 9.8          # gravity



def sigma0(XT, XS):

    #--------------------------------------------------------------------
    #                  ***  FUNCTION sigma0  ***
    #
    # ** Purpose : Compute the in situ density (ratio rho/rau0) and the
    #              potential volumic mass (Kg/m3) from potential temperature 
    #              and salinity fields using an equation of state defined 
    #              through the namelist parameter neos.
    #
    # ** Method  : Jackett and McDougall (1994) equation of state.
    #              The in situ density is computed directly as a function of
    #              potential temperature relative to the surface (the opa t
    #              variable), salt and pressure (assuming no pressure variation
    #              along geopotential surfaces, i.e. the pressure p in decibars
    #              is approximated by the depth in meters.
    #              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    #              rhop(t,s)  = rho(t,s,0)
    #              with pressure                      p        decibars
    #                   potential temperature         t        deg celsius
    #                   salinity                      s        psu
    #                   reference volumic mass        rau0     kg/m**3
    #                   in situ volumic mass          rho      kg/m**3
    #                   in situ density anomalie      prd      no units
    #
    #----------------------------------------------------------------------

    Vshape = nmp.shape(XT)
    nsize  = nmp.size(XT) ; # number of elements in array
    nbdim = len(Vshape)   ; # number of dimensions...
    
    if nmp.shape(XS) != Vshape:
        print 'ERROR in sigma0.lb_ocean_orca1: XT and XS dont have the same shape!'
        sys.exit(0)

    zr1 = nmp.zeros(nsize); zr1.shape = Vshape
    zr2 = nmp.zeros(nsize); zr2.shape = Vshape
    zr3 = nmp.zeros(nsize); zr3.shape = Vshape    
    zr4 = nmp.zeros(nsize); zr4.shape = Vshape
            
    # compute volumic mass pure water at atm pressure
    zr1 = ( ( ( ( 6.536332e-9*XT -1.120083e-6 )*XT + 1.001685e-4)*XT - 9.095290e-3 )*XT + 6.793952e-2 )*XT + 999.842594

    # seawater volumic mass atm pressure
    zr2 = ( ( ( 5.3875e-9*XT-8.2467e-7 )*XT + 7.6438e-5 )*XT - 4.0899e-3 ) *XT + 0.824493
    zr3 = ( -1.6546e-6*XT+1.0227e-4 )*XT - 5.72466e-3
    zr4 = 4.8314e-4
    
    # potential volumic mass (reference to the surface)
    return ( zr4*XS + zr3*nmp.sqrt(nmp.abs(XS)) + zr2 ) *XS + zr1 - zrau0

