# Misc :

import sys
import numpy as nmp
import math

rt0 = 273.16


grav  = 9.8          # gravity
Rgas  = 287.04     
Patm  = 101000.    
ctv   = 0.608        # for virtual temperature
eps   = 0.62197      # humidity constant
cte   = 0.622     
kappa = 0.4          # Von Karman's constant
Cp    = 1000.5    
Pi    = 3.141592654 
eps_w = 0.987        # emissivity of water
sigma = 5.67E-8      # Stefan Boltzman constamt
alfa  = 0.066        # Surface albedo over ocean


sensit = 0.1



def Lvap(zsst):
    #
    # INPUT  : zsst => water temperature in [K]
    # OUTPUT : Lvap => Latent Heat of Vaporization [J/Kg]
    return ( 2.501 - 0.00237*(zsst - rt0) )*1.E6

def e_sat(rt):
    # vapour pressure at saturation  [Pa]
    # rt      ! temperature (K)
    zrtrt0 = rt/rt0

    return 100*( nmp.power(10.,(10.79574*(1. - rt0/rt) - 5.028*nmp.log10(zrtrt0)     
                 + 1.50475*0.0001*(1. - nmp.power(10.,(-8.2969*(zrtrt0 - 1.))) )
                 + 0.42873*0.001 *(nmp.power(10.,(4.76955*(1. - rt0/rt))) - 1.) + 0.78614 ) ) )




def e_air(q_air, zslp):
    #
    #--------------------------------------------------------------------
    #                  **** Function e_air ****
    #
    # Gives vapour pressure of air from pressure and specific humidity
    #
    #--------------------------------------------------------------------
    #

    diff  = 1.E8
    e_old = q_air*zslp/eps

    while diff > 1.:
        #print "Again... diff = ", diff
        ee = q_air/eps*(zslp - (1. - eps)*e_old)
        diff  = nmp.sum(abs( ee - e_old ))
        e_old = ee

    return ee






#def rh_air(q_air, t_air, zslp)
#    #
#    REAL             ::    rh_air      #: relative humidity             [%]
#    REAL, INTENT(in) ::          &
#         &                 q_air,   &  #: specific humidity of air      [kg/kg]
#         &                 t_air,   &  #: air temperature               [K]
#         &                 zslp        #: atmospheric pressure          [Pa]
#    #
#    REAL             :: ea, es
#    #
#    #
#    ea = e_air(q_air, zslp)
#    es = e_sat(t_air)
#    #
#    rh_air = ea/es
#    #
#  END FUNCTION rh_air


  
#def q_air_rh(rha, ta, zslp)
#    # Specific humidity from RH 
#    REAL, DIMENSION(ni,nj) :: q_air_rh
#    INTEGER, INTENT(in)    :: ni, nj
#    
#    REAL, DIMENSION(ni,nj), INTENT(in) :: &
#         &     rha,     &   !: relative humidity      [fraction, not %#!]
#         &     ta,      &   !: air temperature        [K]
#         &     zslp         !: atmospheric pressure          [Pa]
#    
#    REAL, DIMENSION(ni,nj) :: ea
#    
#    ea = rha*e_sat(ni,nj, ta)
##    
 #   q_air_rh = ea*eps/(zslp - (1. - eps)*ea)

#  END FUNCTION q_air_rh




#def q_air_dp(da, zslp)
    #
    # Air specific humidity from dew point temperature
    #
#    INTEGER, INTENT(in) :: ni, nj
#    REAL, DIMENSION(ni,nj) :: q_air_dp  !: kg/kg
    #
#    REAL, DIMENSION(ni,nj), INTENT(in) :: &
#         &     da,     &    !: dew-point temperature   [K]
#         &     zslp         !: atmospheric pressure    [Pa]
    #
#    q_air_dp = e_sat(da)*eps/(zslp - (1. - eps)*e_sat(da))
    #


  #
  #
  #
  # Humidity :
  # ----------
  #            - ea is the water vapour pressure  (h.Pa)
  #            - qa is the specific hymidity      (g/kg)
  #              rqa = rqa/1000.     ! puts specific humidity in kg/kg instead of g/kg    
  #             rea = (rqa*rpa)/(0.378*rqa + cte)
  #
  #    Virtual temperature :
  #
  # 
  # Tv = T*(1 + 0.608*q)    
  #
  # eps = 0.622        --> 0.608 = (1 - eps) / eps
  #
  


  
def rho_air(zt, zq, zP):
    #
    #INTEGER, INTENT(in)    :: ni, nj
    #REAL, DIMENSION(ni,nj) ::   rho_air      !: density of air [kg/m^3] 
    #REAL, DIMENSION(ni,nj), INTENT(in) ::  &
    #     &      zt,       &     !: air temperature in (K)
    #     &      zq,       &     !: air spec. hum. (kg/kg)
    #     &      zP              !: pressure in       (Pa)
    #
    rho_air = zP/(Rgas*zt*(1. + ctv*zq))
    return rho_air




def q_sat(zsst, zslp):
    #
    #REAL, DIMENSION(ni,nj) :: q_sat
    #INTEGER, INTENT(in) :: ni, nj
    #REAL, DIMENSION(ni,nj), INTENT(in) ::  &
    #     &                  zsst,  &   !: sea surface temperature         [K]  
    #     &                  zslp       !: sea level atmospheric pressure  [Pa]
    #
    #
    # Local :
    # -------
    #REAL, DIMENSION(ni,nj) ::  &
    #     &    e_s
    #
    #
    #
    # Specific humidity at saturation
    # -------------------------------
    #
    # Vapour pressure at saturation :
    e_s = 100*(10.^(10.79574*(1-rt0/zsst)-5.028*math.log10(zsst/rt0)  \
                   + 1.50475*10.^(-4)*(1 - 10.^(-8.2969*(zsst/rt0 - 1)) )    \
                   + 0.42873*10.^(-3)*(10.^(4.76955*(1 - rt0/zsst)) - 1) + 0.78614) )
    #
    return eps*e_s/(zslp - (1. - eps)*e_s)






#def e_sat_ice(ni,nj, zrt)
    #
#    INTEGER, INTENT(in) :: ni, nj
#    REAL, DIMENSION(ni,nj) :: e_sat_ice !: vapour pressure at saturation in presence of ice [Pa]
#    REAL, DIMENSION(ni,nj), INTENT(in) :: zrt
#    #
#    e_sat_ice = 100.*(10.^( -9.09718*(273.16/zrt - 1.) - 3.56654*math.log10(273.16/zrt) &
#         &                  + 0.876793*(1. - zrt/273.16) + math.log10(6.1071) ) )
#    #





#def q_sat_simple(zsst)
#    REAL, DIMENSION(ni,nj) :: q_sat_simple
#    INTEGER, INTENT(in)    :: ni, nj
#    REAL, DIMENSION(ni,nj), INTENT(in) ::  &
#         &                  zsst     !: sea surface temperature         [K]
#    
#    q_sat_simple = 640380./1.22 * exp(-5107.4/zsst)
#    

  
#def q_sat_simple_with_rho(zsst, zslp)
#    REAL, DIMENSION(ni,nj) :: q_sat_simple_with_rho
#    INTEGER, INTENT(in)    :: ni, nj
#    REAL, DIMENSION(ni,nj), INTENT(in) ::  &
#         &                  zsst, &     !: sea surface temperature         [K]
#         &                  zslp        !: sea level atmospheric pressure  [Pa]
#    REAL, DIMENSION(ni,nj) :: ztmp, ztmp2
#    ! we need to know specific humidity to get a good estimate of density:
#    ztmp2 = 0.99 ! RH! air is saturated #!
#    ztmp  = 0.98 * q_air_rh(ztmp2, zsst, zslp)
#    
#    ztmp2 = 0.0
#    
#    ztmp2 = rho_air(zsst, ztmp, zslp) ! rho_air
#
#    q_sat_simple_with_rho = 640380./ztmp2 * exp(-5107.4/zsst)




#def e_sat(rt):
#    #  Vapour pressure at saturation for a given temperature [Pa]
#    #
#    #  * rt      ! temperature (K)
#    #
#    rt0 = 273.16
#    #
#    e_sat = 100*( 10**(10.79574*(1 - rt0/rt) - 5.028*log10(rt/rt0)
#                       + 1.50475*10**(-4)*(1 - 10**(-8.2969*(rt/rt0 - 1)) )
#                       + 0.42873*10**(-3)*(10**(4.76955*(1 - rt0/rt)) - 1) + 0.78614) )
#    return e_sat



def qa_e_p(res, rp):
    #  Specific humidity from pressure and vapour pressure at saturation
    #
    #  * res    : vapour pressure at saturation [Pa]
    #  * rp     : atmospheric pressure [Pa]
    #
    reps = 0.62197
    #
    qa_e_p = reps*res / ( rp - (1. - reps)*res )
    #
    return qa_e_p








