#!/usr/bin/env python

#       B a r a K u d a
#
#     Budget and other stuffs on rectangular boxes!
#
#       L. Brodeau, 2013

import sys
import numpy as nmp
from netCDF4 import Dataset
from os.path import basename

import barakuda_tool as bt
import barakuda_ncio as bnc
from barakuda_physics import sigma0


#next = 10
next = 10


# The density of seawater is about 1025 kg/m^3 and the specific heat is about 3850 J/(kg C)
#rho0 = 1025.  ; # kg/m^3
#rCp  = 3850.  ; # 3850 J/kg/deg.C
rho0 = 1000.  ; # kg/m^3           => to stay consistent with cdftransportiz.f90 of CDFTOOLS...
rCp  = 4000.  ; # 3850 J/kg/deg.C  => to stay consistent with cdftransportiz.f90 of CDFTOOLS...

#l_plot_debug = True
l_plot_debug = False




venv_needed = {'ORCA','RUN','CPREF','DIAG_D','MM_FILE','FILE_DEF_BOXES','NN_T','NN_S', \
               'NN_SST','NN_SSS','NN_SSH'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

corca = vdic['ORCA']

CONFRUN = corca+'-'+vdic['RUN']




cname_script = basename(sys.argv[0])

print '\n'+cname_script

narg = len(sys.argv)
if narg < 3 or narg > 4:
    print 'Usage: '+cname_script+' <year> <depth for surface properties> (uv)'
    print '      by specifying "uv" as 3rd argument, budget will be extended to'
    print '      some variables found in grid_U and grid_V files such as wind stress\n'
    sys.exit(0)
cy = sys.argv[1] ; jy=int(cy)
czs= sys.argv[2] ; zs = float(int(czs))


# Shall we need U and V files???
luv = False
if narg == 4 and sys.argv[3] == 'uv':
    luv = True
    venv_uv = {'NN_TAUX','NN_TAUY'}
    vdic_uv = bt.check_env_var(sys.argv[0], venv_uv)



path_fig = vdic['DIAG_D']+'/'

# Image type? eps, png, jpg...
#FIG_FORM = 'pdf'
FIG_FORM = 'png'



# First will read name and coordinates of rectangular boxes to treat into file FILE_DEF_BOXES
##############################################################################################
vboxes, vi1, vj1, vi2, vj2 = bt.read_box_coordinates_in_ascii(vdic['FILE_DEF_BOXES'])
nbb = len(vboxes)
print ''





# Checking presence of NEMO files:
cfroot  = vdic['CPREF']+cy+'0101_'+cy+'1231'
cf_in_T = cfroot+'_grid_T.nc'; bt.chck4f(cf_in_T, script_name=cname_script)
if luv:
    cf_in_U = cfroot+'_grid_U.nc'; bt.chck4f(cf_in_U, script_name=cname_script)
    cf_in_V = cfroot+'_grid_V.nc'; bt.chck4f(cf_in_V, script_name=cname_script)
    

# Coordinates, mask and metrics:
bt.chck4f(vdic['MM_FILE'], script_name=cname_script)
id_mm = Dataset(vdic['MM_FILE'])
Xmask = id_mm.variables['tmask']   [0,:,:,:]
ze1t  = id_mm.variables['e1t']     [0,:,:]
ze2t  = id_mm.variables['e2t']     [0,:,:]
ve3t  = id_mm.variables['e3t_1d']  [0,:]
zlon  = id_mm.variables['nav_lon'] [:,:]
zlat  = id_mm.variables['nav_lat'] [:,:]
vzt   = id_mm.variables['gdept_1d'][0,:]
vzw   = id_mm.variables['gdepw_1d'][0,:]
id_mm.close()




lqnet = False ;  lqsw = False ;  lpme = False; ltau = False


# NEMO output, Grid T
# ~~~~~~~~~~~~~~~~~~~
id_in_T = Dataset(cf_in_T)

list_variables = id_in_T.variables.keys()

Vtime = id_in_T.variables['time_counter'][:]

if vdic['NN_SST'] == 'thetao':
    Zsst  = id_in_T.variables[vdic['NN_SST']][:,0,:,:]
else:
    Zsst  = id_in_T.variables[vdic['NN_SST']][:,:,:]
    
if vdic['NN_SSS'] == 'so':
    Zsss  = id_in_T.variables[vdic['NN_SSS']][:,0,:,:]
else:
    Zsss  = id_in_T.variables[vdic['NN_SSS']][:,:,:]
    
Zssh  = id_in_T.variables[vdic['NN_SSH']][:,:,:]

Xtemp = id_in_T.variables[vdic['NN_T']][:,:,:,:]
Xsali = id_in_T.variables[vdic['NN_S']][:,:,:,:]

if 'sohefldo' in list_variables[:]:
    Zqnet = id_in_T.variables['sohefldo']  [:,:,:] ; # Net Downward Heat Flux
    lqnet = True
if 'tohfls' in list_variables[:]:
    Zqnet = id_in_T.variables['tohfls']  [:,:,:] ; # Net Downward Heat Flux
    lqnet = True
    
if 'soshfldo' in list_variables[:]:
    Zqsw = id_in_T.variables['soshfldo']  [:,:,:] ; # Shortwave Radiation
    lqsw = True
if 'rsntds' in list_variables[:]:
    Zqsw = id_in_T.variables['rsntds']  [:,:,:] ; # Shortwave Radiation
    lqsw = True

# Want PmE (positive when ocean gains FW), in NEMO files its the oposite EmP...
if 'wfo' in list_variables[:]:
    Zpme = -id_in_T.variables['wfo']  [:,:,:] ; # wfo same as below = EmP, > 0 when ocean losing water
    lpme = True
if 'sowaflup' in list_variables[:]:
    Zpme = -id_in_T.variables['sowaflup']  [:,:,:] ; # Net Downward Heat Flux
    #                                                # sowaflup = EmP (>0 if more evaporation than P)
    lpme = True
    
print '(has ',Xtemp.shape[0],' time snapshots)\n'
id_in_T.close()

[ Nt, nk, nj, ni ] = Xtemp.shape ; print 'Nt, nk, nj, ni =', Nt, nk, nj, ni


Zss0 = sigma0( Zsst, Zsss )


print len(ve3t[:])
if len(ve3t[:]) != nk: print 'Problem with nk!!!'; sys.exit(0)





# NEMO output, Grid U and V
# ~~~~~~~~~~~~~~~~~~~~~~~~~
if luv:
    id_in_U = Dataset(cf_in_U)
    list_variables = id_in_U.variables.keys()
    if vdic_uv['NN_TAUX'] in list_variables[:]:
        Ztaux = id_in_U.variables[vdic_uv['NN_TAUX']][:,:,:] ; # Net Downward Heat Flux
        ltau = True
        print vdic_uv['NN_TAUX']+' found in '+cf_in_U
    else:
        print vdic_uv['NN_TAUX']+' NOT found in '+cf_in_U
    id_in_U.close()
    id_in_V = Dataset(cf_in_V)
    list_variables = id_in_V.variables.keys()
    if ltau and vdic_uv['NN_TAUY'] in list_variables[:]:
        Ztauy = id_in_V.variables[vdic_uv['NN_TAUY']][:,:,:] ; # Net Downward Heat Flux
        print vdic_uv['NN_TAUY']+' found in '+cf_in_V+'\n'
    else:
        print vdic_uv['NN_TAUY']+' NOT found in '+cf_in_V
        ltau = False
    id_in_V.close()



if ltau:

    # Must interpolate Taux and Tauy on T-grid:
    Ztau = nmp.zeros(Nt*nj*ni) ; Ztau.shape = [ Nt, nj, ni ]
    xtmp1 = nmp.zeros(nj*ni) ; xtmp1.shape = [ nj, ni ]
    xtmp2 = nmp.zeros(nj*ni) ; xtmp2.shape = [ nj, ni ]

    for jm in range(Nt):
        xtmp1[:,1:]  = 0.5*(Ztaux[jm,:,1:]+Ztaux[jm,:,:ni-1]) ; # u on Tgrid
        xtmp2[1:,:]  = 0.5*(Ztauy[jm,1:,:]+Ztauy[jm,:nj-1,:]) ; # v on Tgrid
        Ztau[jm,:,:] = nmp.sqrt(xtmp1*xtmp1 + xtmp2*xtmp2)


    #print Ztau[3,100,:]
    #del Ztaux, Ztaux, xtmp1, xtmp2





jks = 0
for jk in range(nk-1):
    if  zs >= vzw[jk] and zs < vzw[jk+1]: jks = jk+1

czs = str(int(round(vzw[jks],0))) ; print czs

print '  *** for depth '+czs+': jks = '+str(jks)+', depthw = '+str(vzw[jks])+' => '+czs+'m'
print '   => will average from jk=0 to jk='+str(jks)+'-1 on T-points => X[:jks-1]'
jks = jks - 1
print '   => that\'s on T-points from jk=0 to jk='+str(jks)+' (depth of deepest T-point used ='+str(vzt[jks])+'m)'


# Loop along boxes:
for jb in range(nbb):

    cbox = vboxes[jb]

    i1 = vi1[jb]
    j1 = vj1[jb]
    i2 = vi2[jb]+1
    j2 = vj2[jb]+1
    
    print '\n *** Treating box '+cbox+' => ', i1, j1, i2-1, j2-1



    # Filling box arrays:
    # ~~~~~~~~~~~~~~~~~~~
    nx_b = i2 - i1
    ny_b = j2 - j1

    shape_array = [ Nt, nk, ny_b, nx_b ]



    XVolu = nmp.zeros(nk*ny_b*nx_b) ;  XVolu.shape = [ nk, ny_b, nx_b ]
    ZArea = nmp.zeros(   ny_b*nx_b) ;  ZArea.shape = [     ny_b, nx_b ]
    Ztmp  = nmp.zeros(   ny_b*nx_b) ;  Ztmp.shape  = [     ny_b, nx_b ]


    Xs0 = nmp.zeros(nk*ny_b*nx_b) ;  Xs0.shape = [ nk, ny_b, nx_b ]

    ssh_m   = nmp.zeros(Nt)    ;  ssh_m.shape  = [ Nt ]    
    sst_m   = nmp.zeros(Nt)    ;  sst_m.shape  = [ Nt ]
    sss_m   = nmp.zeros(Nt)    ;  sss_m.shape  = [ Nt ]
    ss0_m   = nmp.zeros(Nt)    ;  ss0_m.shape  = [ Nt ]

    surf_T_m   = nmp.zeros(Nt)    ;  surf_T_m.shape  = [ Nt ]
    surf_S_m   = nmp.zeros(Nt)    ;  surf_S_m.shape  = [ Nt ]
    surf_s0_m   = nmp.zeros(Nt)    ;  surf_s0_m.shape  = [ Nt ]

    T_m     = nmp.zeros(Nt)    ;  T_m.shape    = [ Nt ]
    Tau_m  = nmp.zeros(Nt) ;  Tau_m.shape = [ Nt ]
    Qnet_m  = nmp.zeros(Nt) ;  Qnet_m.shape = [ Nt ]
    Qnet_x_S_m  = nmp.zeros(Nt) ;  Qnet_x_S_m.shape = [ Nt ]
    Qsw_m   = nmp.zeros(Nt) ;  Qsw_m.shape   = [ Nt ]
    Qsw_x_S_m   = nmp.zeros(Nt) ;  Qsw_x_S_m.shape   = [ Nt ]
    PmE_m   = nmp.zeros(Nt) ; PmE_m.shape   = [ Nt ]
    H_m     = nmp.zeros(Nt)    ;  H_m.shape    = [ Nt ] ; # Heat coNtent
    S_m     = nmp.zeros(Nt)    ;  S_m.shape    = [ Nt ]
    Vol_m   = nmp.zeros(Nt)    ;  Vol_m.shape  = [ Nt ] ; # Volume derived from SSH!
    
    
    
    
    # On the sea of the box only:
    
    ZArea[:,:] = ze1t[j1:j2, i1:i2]*ze2t[j1:j2, i1:i2]

    
    for jk in range(nk): XVolu[jk,:,:] = Xmask[jk, j1:j2, i1:i2]*ZArea[:,:]*ve3t[jk]

    ZArea[:,:] = Xmask[0, j1:j2, i1:i2]*ZArea[:,:]

    #if l_plot_debug: bp.check_with_fig_2(ZArea, ZArea*0.+1., 'ZArea', fig_type=FIG_FORM)


    
    Tot_area = nmp.sum(ZArea) ; print 'Total area of '+cbox+' (m^2): ', Tot_area
    
    Tot_vol  = nmp.sum(XVolu)
    
    Tot_vol_jks  = nmp.sum(XVolu[:jks,:,:])
    
    
    for jm in range(Nt):
    
    
        # 3D sigma0 density for current month
        Xs0[:,:,:] = 0.    
        Xs0[:,:,:] = sigma0( Xtemp[jm,:, j1:j2, i1:i2], Xsali[jm,:, j1:j2, i1:i2] )


        # Mean SSH
        ssh_m[jm] = nmp.sum(Zssh[jm, j1:j2, i1:i2]*ZArea) / Tot_area

        # Mean SST
        sst_m[jm] = nmp.sum(Zsst[jm, j1:j2, i1:i2]*ZArea) / Tot_area

        # Mean SSS
        sss_m[jm] = nmp.sum(Zsss[jm, j1:j2, i1:i2]*ZArea) / Tot_area

        # Mean SS0
        ss0_m[jm] = nmp.sum(Zss0[jm, j1:j2, i1:i2]*ZArea) / Tot_area

    
        # Mean surface temp (first Xm)
        surf_T_m[jm] = nmp.sum(Xtemp[jm, :jks, j1:j2, i1:i2]*XVolu[:jks,:,:]) / Tot_vol_jks
        
        # Mean temperature
        T_m[jm]      = nmp.sum(Xtemp[jm,  :  , j1:j2, i1:i2]*XVolu[:,:,:]) / Tot_vol
        
        # Heat content in Peta Joules:
        H_m[jm] = nmp.sum(Xtemp[jm,:, j1:j2, i1:i2]*(XVolu/1.E6))*rho0*rCp * 1.E-9 ; # => PJ (E15)



        # Mean surface salinity (first Xm)
        surf_S_m[jm] = nmp.sum(Xsali[jm,:jks, j1:j2, i1:i2]*XVolu[:jks,:,:]) / Tot_vol_jks        
        # Mean salinity
        S_m[jm] = nmp.sum(Xsali[jm,:, j1:j2, i1:i2]*XVolu) / Tot_vol
    
    
        # Mean surface sigma0 density (first Xm)
        surf_s0_m[jm] = nmp.sum(Xs0[:jks,:,:]*XVolu[:jks,:,:]) / Tot_vol_jks    
    
        
        # Sea-ice area
        #Aice_m[jm] = nmp.sum(Zicec[jm,:,:]*ZArea) * 1.E-12; # Million km^2
    
    
        # For open-sea:
        #Ztmp[:,:] = 1. - Zicec[jm,:,:] ; # 1 => 100% open sea / 0 => 100% ice
        Ztmp[:,:] = 1.
    
    
        # ZArea is in m^2

        if ltau:
            # Surface wind stress:
            Tau_m[jm] = nmp.sum(Ztau[jm, j1:j2, i1:i2]*ZArea) / Tot_area


        
        # Surface heat flux
        if lqnet:
            rr = nmp.sum(Zqnet[jm, j1:j2, i1:i2]*ZArea)
            Qnet_m[jm]     = rr / Tot_area    # W/m^2
            Qnet_x_S_m[jm] = rr * 1.E-15      # PW
        
        # Shortwave heat flux
        if lqsw:
            rr = nmp.sum(Zqsw[jm, j1:j2, i1:i2]*ZArea)
            Qsw_m[jm]      = rr / Tot_area    # W/m^2
            Qsw_x_S_m[jm]  = rr * 1.E-15      # PW
        
        # PmE
        if lpme:  PmE_m[jm]  = nmp.sum( Zpme[jm, j1:j2, i1:i2]*ZArea) * 1.E-9; # (Sv) 1 kg/m^2/s x S => 1E-3 m^3/s
        #                                                         #     1 Sv = 1E6 m^3
        
        # Volume associated with SSH
        Vol_m[jm] = nmp.sum(Zssh[jm, j1:j2, i1:i2]*ZArea) * 1.E-9; # km^3
    
    
    
        if l_plot_debug:
            VMN = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
            print ' *** month ', str(jm+1), '***'
            print 'Mean sst = ', sst_m[jm]
            print 'Mean ssh = ', ssh_m[jm]
            print 'Mean sss (0-'+czs+'m) = ', sss_m[jm]
            print 'Mean ss0 (0-'+czs+'m) = ', ss0_m[jm]
            print 'Mean T   = ', T_m[jm]
            print 'Heat content (PJ) / ref T=0C   => ', H_m[jm]
            if jm>0: print 'Volume heat flux (PW) / ref T=0C   => ', (H_m[jm] - H_m[jm-1]) / (VMN[jm-1]*24.*3600.)
            print 'Mean S   = ', S_m[jm]
            #print 'Sea-ice (10^6 km^2) = ', Aice_m[jm]
            print 'Shortwave Radiation  (PW) = ', Qsw_m[jm]
            print 'Net surface heat flux (PW) = ', Qnet_m[jm], '\n'
            print 'Surface freshwater flux (Sv) = ', PmE_m[jm], '\n'        
            print 'Volume associated with SSH (km^3) = ', Vol_m[jm], '\n'
    


    Vtime = nmp.zeros(Nt)
    for jm in range(Nt): Vtime[jm] = float(jy) + (float(jm)+0.5)/12.
    
    cc = 'Box-averaged '
    cf_out = vdic['DIAG_D']+'/budget_'+CONFRUN+'_box_'+cbox+'.nc'
    
    bnc.wrt_appnd_1d_series(Vtime, ssh_m, cf_out, 'ssh',  cu_t='year', cu_d='m', cln_d=cc+'sea surface height',
                            vd2=sst_m,    cvar2='sst',      cln_d2=cc+'sea surface temperature',      cun2='deg.C',
                            vd3=sss_m,    cvar3='sss',      cln_d3=cc+'sea surface salinity',         cun3='PSU',
                            vd4=surf_T_m, cvar4='surf_T',   cln_d4=cc+'Temperature (first '+czs+'m)', cun4='deg.C',
                            vd5=T_m, cvar5='theta',    cln_d5=cc+'potential temperature',        cun5='deg.C',
                            vd6=H_m, cvar6='HC',       cln_d6=cc+'heat content',                 cun6='PJ',
                            vd7=surf_S_m, cvar7='surf_S',   cln_d7=cc+'salinity (first '+czs+'m)',    cun7='PSU',
                            vd8=S_m, cvar8='S',        cln_d8=cc+'salinity',                     cun8='PSU',                    
                            vd9=ss0_m, cvar9='SSs0',     cln_d9=cc+'sea surface sigma0 (sst&sss)', cun9='',
                            vd10=surf_s0_m,cvar10='surf_s0', cln_d10=cc+'surface sigma0 (first '+czs+'m)', cun10=''
                            )
    

    cf_out = vdic['DIAG_D']+'/budget_srf_flx_'+CONFRUN+'_box_'+cbox+'.nc'
    
    bnc.wrt_appnd_1d_series(Vtime, Qnet_m, cf_out, 'Qnet',  cu_t='year', cu_d='W/m^2', cln_d=cc+'net heat flux',
                            vd2=Qnet_x_S_m,    cvar2='Qnet_x_S',      cln_d2=cc+'net heat flux x Surface',      cun2='PW',
                            vd3=Qsw_m,    cvar3='Qsw',      cln_d3=cc+'solar radiation',         cun3='W/m^2',
                            vd4=Qsw_x_S_m, cvar4='Qsw_x_S',   cln_d4=cc+'solar radiation x Surface', cun4='PW',
                            vd5=PmE_m, cvar5='PmE',    cln_d5=cc+'net freshwater flux',        cun5='Sv',
                            vd6=Tau_m, cvar6='Tau',       cln_d6=cc+'wind stress module',                 cun6='N/m^2'
                            )

