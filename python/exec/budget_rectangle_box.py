#!/home/x_laubr/bin/python

# L. Brodeau, May 2014

# TO DO: test for the presence of variables like fluxes to know if treat them
#   put real sst and sss
# rename SST and SSS and ... to surf_theta, surf_s, ....


import sys
import os
import numpy as nmp
from netCDF4 import Dataset
from os.path import basename


# Laurent's:
import barakuda_tool as bt
import barakuda_plot as bp
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


CONF = os.getenv('CONF')
if CONF == None: print 'The CONF environement variable is no set'; sys.exit(0)

RUN = os.getenv('RUN')
if RUN == None: print 'The RUN environement variable is no set'; sys.exit(0)
CONFRUN = CONF+'-'+RUN

CPREF = os.getenv('CPREF')
if CPREF == None: print 'The CPREF environement variable is no set'; sys.exit(0)

MM_FILE = os.getenv('MM_FILE')
if MM_FILE == None: print 'The MM_FILE environement variable is no set'; sys.exit(0)

DIAG_D = os.getenv('DIAG_D')
if DIAG_D == None: print 'The DIAG_D environement variable is no set'; sys.exit(0)

FILE_DEF_BOXES = os.getenv('FILE_DEF_BOXES')
if FILE_DEF_BOXES == None: print 'The FILE_DEF_BOXES environement variable is no set'; sys.exit(0)

# Name for salinity and pot. temperature:
NN_T = os.getenv('NN_T')
if NN_T == None: print 'The NN_T environement variable is no set'; sys.exit(0)
NN_S = os.getenv('NN_S')
if NN_S == None: print 'The NN_S environement variable is no set'; sys.exit(0)

NN_SST = os.getenv('NN_SST')
if NN_SST == None: print 'The NN_SST environement variable is no set'; sys.exit(0)
NN_SSS = os.getenv('NN_SSS')
if NN_SSS == None: print 'The NN_SSS environement variable is no set'; sys.exit(0)
NN_SSH = os.getenv('NN_SSH')
if NN_SSH == None: print 'The NN_SSH environement variable is no set'; sys.exit(0)





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
    NN_TAUX = os.getenv('NN_TAUX')
    if NN_TAUX == None: print 'The NN_TAUX environement variable is no set'; sys.exit(0)
    NN_TAUY = os.getenv('NN_TAUY')
    if NN_TAUY == None: print 'The NN_TAUY environement variable is no set'; sys.exit(0)





print ' RUN = '+RUN; print ' CONFRUN = '+CONFRUN
print ' CPREF = '+CPREF
print ' DIAG_D = '+DIAG_D

path_fig = DIAG_D+'/'

# Image type? eps, png, jpg...
#FIG_FORM = 'pdf'
FIG_FORM = 'png'



# First will read name and coordinates of rectangular boxes to treat into file FILE_DEF_BOXES
##############################################################################################
vboxes, vi1, vj1, vi2, vj2 = bt.read_box_coordinates_in_ascii(FILE_DEF_BOXES)
nbb = len(vboxes)
print ''





# Checking presence of NEMO files:
cfroot = CPREF+cy+'0101_'+cy+'1231'
cf_in_T = cfroot+'_grid_T.nc'; bt.chck4f(cf_in_T, script_name=cname_script)
if luv:
    cf_in_U = cfroot+'_grid_U.nc'; bt.chck4f(cf_in_U, script_name=cname_script)
    cf_in_V = cfroot+'_grid_V.nc'; bt.chck4f(cf_in_V, script_name=cname_script)
    

# Coordinates, mask and metrics:
bt.chck4f(MM_FILE, script_name=cname_script)
id_mm = Dataset(MM_FILE)
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

if NN_SST == 'thetao':
    Zsst  = id_in_T.variables[NN_SST][:,0,:,:]
else:
    Zsst  = id_in_T.variables[NN_SST][:,:,:]
    
if NN_SSS == 'so':
    Zsss  = id_in_T.variables[NN_SSS][:,0,:,:]
else:
    Zsss  = id_in_T.variables[NN_SSS][:,:,:]
    
Zssh  = id_in_T.variables[NN_SSH][:,:,:]

Xtemp = id_in_T.variables[NN_T][:,:,:,:]
Xsali = id_in_T.variables[NN_S][:,:,:,:]

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
    if NN_TAUX in list_variables[:]:
        Ztaux = id_in_U.variables[NN_TAUX][:,:,:] ; # Net Downward Heat Flux
        ltau = True
        print NN_TAUX+' found in '+cf_in_U
    else:
        print NN_TAUX+' NOT found in '+cf_in_U
    id_in_U.close()
    id_in_V = Dataset(cf_in_V)
    list_variables = id_in_V.variables.keys()
    if ltau and NN_TAUY in list_variables[:]:
        Ztauy = id_in_V.variables[NN_TAUY][:,:,:] ; # Net Downward Heat Flux
        print NN_TAUY+' found in '+cf_in_V+'\n'
    else:
        print NN_TAUY+' NOT found in '+cf_in_V
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
    if l_plot_debug: bp.check_with_fig_2(ZArea, ZArea*0.+1., 'ZArea', fig_type=FIG_FORM)


    
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
    
    
    
    
    
    
    
    # NETCDF:
    
    cf_out = DIAG_D+'/budget_'+CONFRUN+'_box_'+cbox+'.nc'
    
    l_nc_is_new = not os.path.exists(cf_out)
    
    # Creating/Opening output Netcdf file:
    if l_nc_is_new:
        f_out = Dataset(cf_out, 'w', format='NETCDF3_CLASSIC')
    else:
        f_out = Dataset(cf_out, 'a', format='NETCDF3_CLASSIC')
    
    if l_nc_is_new:
        jrec2write = 0
        
        # Creating Dimensions:
        f_out.createDimension('time', None)
    
        # Creating variables:
        id_t = f_out.createVariable('time','f4',('time',)) ;      id_t.units = 'year'


        id_v01 = f_out.createVariable('ssh',    'f4',('time',))
        id_v01.unit = 'm'; id_v01.long_name = 'Box-averaged sea surface height'
        
        id_v02 = f_out.createVariable('sst',    'f4',('time',))
        id_v02.unit = 'deg.C'; id_v02.long_name = 'Box-averaged sea surface temperature'
        
        id_v03 = f_out.createVariable('sss',    'f4',('time',))
        id_v03.unit = 'PSU'; id_v03.long_name = 'Box-averaged sea surface salinity'
        
    
        id_v04 = f_out.createVariable('surf_T',  'f4',('time',))
        id_v04.unit = 'deg.C'; id_v04.long_name = 'Box-averaged Temperature (first '+czs+'m)'
    
        id_v05 = f_out.createVariable('theta','f4',('time',))
        id_v05.unit = 'deg.C'; id_v05.long_name = 'Box-averaged potential temperature'
    
        id_v06 = f_out.createVariable('H',    'f4',('time',))
        id_v06.unit = 'PJ'; id_v06.long_name = 'Box-averaged heat content'


        id_v07 = f_out.createVariable('surf_S',  'f4',('time',))
        id_v07.unit = 'PSU'; id_v07.long_name = 'Box-averaged salinity (first '+czs+'m)'
    
        id_v08 = f_out.createVariable('s',    'f4',('time',))
        id_v08.unit = 'PSU'; id_v08.long_name = 'Box-averaged salinity'
        

        id_v09 = f_out.createVariable('SSsigma0',  'f4',('time',))
        id_v09.unit = ''; id_v09.long_name = 'Box-averaged sea surface sigma0 (from sst and sss)'

        id_v09b = f_out.createVariable('surf_sigma0',  'f4',('time',))
        id_v09b.unit = ''; id_v09b.long_name = 'Box-averaged surface sigma0 (first '+czs+'m)'

        if lqnet:
            id_v10 = f_out.createVariable('Qnet',    'f4',('time',))
            id_v10.unit = 'W/m^2'; id_v10.long_name = 'Box-averaged net heat flux'
            id_v10b = f_out.createVariable('Qnet_x_S',    'f4',('time',))
            id_v10b.unit = 'PW'; id_v10b.long_name = 'Box-averaged net heat flux x Surface'
        
        if lqsw:
            id_v11 = f_out.createVariable('Qsw',    'f4',('time',))
            id_v11.unit = 'W/m^2'; id_v11.long_name = 'Box-averaged solar radiation'        
            id_v11b = f_out.createVariable('Qsw_x_S',    'f4',('time',))
            id_v11b.unit = 'PW'; id_v11b.long_name = 'Box-averaged solar radiation x Surface'
        
        if lpme:
            id_v12 = f_out.createVariable('PmE',    'f4',('time',))
            id_v12.unit = 'Sv'; id_v12.long_name = 'Box-averaged net freshwater flux'
        
        if ltau:
            id_v13 = f_out.createVariable('Tau',    'f4',('time',))
            id_v13.unit = 'N/m^2'; id_v13.long_name = 'Box-averaged wind stress module'
    
    
    
        for jm in range(Nt):
            id_t[jrec2write+jm] = float(jy) + (float(jm)+0.5)/12.
            id_v01[jrec2write+jm] =  ssh_m[jm]
            id_v02[jrec2write+jm] =  sst_m[jm]
            id_v03[jrec2write+jm] =  sss_m[jm]            
            
            id_v04[jrec2write+jm] =  surf_T_m[jm]
            id_v05[jrec2write+jm] =    T_m[jm]
            id_v06[jrec2write+jm] =    H_m[jm]

            id_v07[jrec2write+jm] =  surf_S_m[jm]
            id_v08[jrec2write+jm] =    S_m[jm]
            
            id_v09[jrec2write+jm] =  ss0_m[jm]
            id_v09b[jrec2write+jm] =  surf_s0_m[jm]

            if lqnet:
                id_v10[jrec2write+jm] = Qnet_m[jm]
                id_v10b[jrec2write+jm] = Qnet_x_S_m[jm]            
            if lqsw:
                id_v11[jrec2write+jm] =  Qsw_m[jm]
                id_v11b[jrec2write+jm] =  Qsw_x_S_m[jm]
                
            if lpme:  id_v12[jrec2write+jm] = PmE_m[jm]
            if ltau: id_v13[jrec2write+jm] = Tau_m[jm]            


        f_out.box_coordinates = cbox+' => '+str(i1)+','+str(j1)+' -> '+str(i2-1)+','+str(j2-1)
        f_out.box_area = str(Tot_area*1.E-6)+' (km^2) = (10^6 m^2)'
        f_out.box_file        = FILE_DEF_BOXES            
        f_out.Author = 'L. Brodeau ('+cname_script+' of Barakuda)'
    
    else:
        vt = f_out.variables['time']
        jrec2write = len(vt)
        v01 = f_out.variables['ssh']
        v02 = f_out.variables['sst']
        v03 = f_out.variables['sss']
        
        v04 = f_out.variables['surf_T']
        v05 = f_out.variables['theta']
        v06 = f_out.variables['H']
        
        v07 = f_out.variables['surf_S']
        v08 = f_out.variables['s']        

        v09 = f_out.variables['SSsigma0']
        v09b = f_out.variables['surf_sigma0']

        if lqnet:
            v10 = f_out.variables['Qnet']
            v10b = f_out.variables['Qnet_x_S']
        if lqsw:
            v11 = f_out.variables['Qsw']
            v11b = f_out.variables['Qsw_x_S']
            
        if lpme:  v12 = f_out.variables['PmE']
        if ltau: v13 = f_out.variables['Tau']

        
        for jm in range(Nt):
            vt [jrec2write+jm] = float(jy) + (float(jm)+0.5)/12.
            v01[jrec2write+jm] =  ssh_m[jm]            
            v02[jrec2write+jm] =  sst_m[jm]
            v03[jrec2write+jm] =  sss_m[jm]
            
            v04[jrec2write+jm] =  surf_T_m[jm]
            v05[jrec2write+jm] =    T_m[jm]
            v06[jrec2write+jm] =    H_m[jm]

            v07[jrec2write+jm] =  surf_S_m[jm]            
            v08[jrec2write+jm] =    S_m[jm]

            v09[jrec2write+jm] =  ss0_m[jm]
            v09b[jrec2write+jm] =  surf_s0_m[jm]

            if lqnet:
                v10[jrec2write+jm] = Qnet_m[jm]
                v10b[jrec2write+jm] = Qnet_x_S_m[jm]
            if lqsw:
                v11[jrec2write+jm] =  Qsw_m[jm]
                v11b[jrec2write+jm] =  Qsw_x_S_m[jm]
                
            if lpme:  v12[jrec2write+jm] = PmE_m[jm]
            if ltau: v13[jrec2write+jm] = Tau_m[jm]
            
            
    f_out.close()
    
    print cf_out+' written!'
    
