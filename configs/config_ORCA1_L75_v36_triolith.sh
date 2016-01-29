#!/bin/bash

#==========================================================
#
#         Configuration file for
#
# OCEAN MONITORING for NEMO v3.6 ORCA1 on 75 levels
#
#            HPC: triolith
#
#        L. Brodeau, 2015
#
#===========================================================

l_clim_diag=true ; # should we try to perform climatology-related diagnostics? (clim must be built!)

export CONF=ORCA1.L75 ; # horizontal global configuration
export NBL=75     ; # number of levels

# Root directory where NEMO output files are stored:
export STORE_DIR="/proj/bolinc/users/x_laubr"

# List of suffixed of files that have been saved by NEMO:
export NEMO_SAVED_FILES="grid_T grid_U grid_V icemod SBC"


# Directory structure in which to find NEMO output file (use <ORCA> and <RUN>):
export NEMO_OUT_STRCT="${STORE_DIR}/<ORCA>/<ORCA>-<RUN>-S"

export TSTAMP="1m"   ; # output time-frequency stamp as in NEMO output files...

# How does the nemo files prefix looks like
# Everything before "<year_related_info>_grid_<X>" or "<year_related_info>_icemod"
# use <ORCA>, <RUN> and <TSTAMP>=>  Ex: export NEMO_FILE_PREFIX="<ORCA>-<RUN>_<TSTAMP>_"
export NEMO_FILE_PREFIX="<ORCA>-<RUN>_<TSTAMP>_"
# => should get rid of TSTAMP actually...

# Temporary file system (scratch) on which to perform the job you can use <JOB_ID> if scracth depends on JOB ID:
export SCRATCH="/scratch/local/<JOB_ID>"

export CANOPY_PATH=${HOME}/opt/Canopy_64bit/User

# If variables names in NEMO files are not the default...
export NN_SST="tos"
export NN_SSS="sos"
export NN_SSH="zos"
export NN_T="thetao"
export NN_S="so"
export NN_MLD="mldr10_1"
export NN_U="uo"
export NN_V="vo"
#export NN_U_EIV="vozoeivu"
#export NN_V_EIV="vomeeivv"
export NN_U_EIV="0" ; # ignore
export NN_V_EIV="0" ; # ignore
export NN_TAUX="tauuo"
export NN_TAUY="tauvo"

export NN_ICEF="siconc" ; # name of ice fraction in "FILE_ICE_SUFFIX" file...
export NN_ICET="sithic" ; # ice thickness but 'sit' is only in icemod file !!!


export L_CONV2NC3=false ; # Set to true if your NEMO output is in Netcdf4 and your NCO does not support netcdf4!

export L_RENAME=false ; # set to true if your ORCA output has old name convention (ex: votemper instead of thetao)


export JTITLE="NEMO v3.6 ${CONF} (L${NBL}) - LIM3 / ocean-only experiment"

# Land-sea mask and basins files:
export MM_FILE="/proj/bolinc/users/x_laubr/${CONF}/mesh_mask_${CONF}_20150929.nc"
export BM_FILE="/proj/bolinc/users/x_laubr/${CONF}/basin_mask_${CONF}_20150929.nc"

# 3D monthly climatologies of potential temperature and salinity (can be those you used for the NEMO run):
export F_T_CLIM_3D_12=${STORE_DIR}/${CONF}/${CONF}-I/thetao_1degx1deg_WOA2009_monthly_${CONF}_cut.nc
export F_S_CLIM_3D_12=${STORE_DIR}/${CONF}/${CONF}-I/so_1degx1deg_WOA2009_monthly_${CONF}_cut.nc
export SST_CLIM_12=${STORE_DIR}/ORCA1.L75/ORCA1.L75-I/tos_180x360-ORCA1_Reynolds_monthly_mean1982-2005.nc
export NN_T_CLIM="thetao"
export NN_S_CLIM="so"
export NN_SST_CLIM="tos"

export ICE_CLIM_12=${STORE_DIR}/ORCA1.L75/ORCA1.L75-I/ice_cover_180x360-ORCA1_Hurrell_monthly_mean1980-1999.nc4
export NN_ICEF_CLIM="ice_cover"


# A text file where the vertical hydraugraphical sections of interest are defined :
#export TRANSPORT_SECTION_FILE="${BARAKUDA_ROOT}/data/transportiz_${CONF}_light.dat"
export TRANSPORT_SECTION_FILE="${BARAKUDA_ROOT}/data/transportiz_ORCA1.dat"

# For transport by sigma-class:
export DENSITY_SECTION_FILE="${BARAKUDA_ROOT}/data/dens_section_ORCA1.dat"


# In what directory of the local machine to save the diagnostics:
export DIAG_DIR="/proj/bolinc/users/x_laubr/tmp/barakuda/${CONF}_v36"


# Files with the list of rectangular boxes to look at more closely:
export FILE_DEF_BOXES="${BARAKUDA_ROOT}/data/def_boxes_convection_ORCA1.txt"



# About remote HOST to install HTML pages to:
ihttp=1 ; # do we export on a remote http server (1) or keep on the local machine (0)
RHOST=misu228.misu.su.se ; # remote host to send diagnostic page to///
RUSER=laurent ; # username associated to remote host (for file export)
RWWWD=/data/www/barakuda/${CONF} ; # directory of the local or remote host to send the diagnostic page to




#########################
# Diags to be performed #
#########################





# Basic 3D and surface averages:
i_do_mean=1

# AMOC:
i_do_amoc=1
export LMOCLAT="20-23 30-33 40-43 45-48 50-53" ; # List of latitude bands to look in for max of AMOC


# Transport of mass, heat and salt through specified sections (into TRANSPORT_SECTION_FILE):
i_do_trsp=2  ; # transport of mass, heat and salt through specified sections
#              # i_do_trsp=2 => treat also different depths range!
z1_trsp=100  ; # first  depth: i_do_trsp must be set to 2
z2_trsp=1000 ; # second depth: i_do_trsp must be set to 2


# meridional heat/salt transport (advective)
i_do_mht=1

# Transport by sigma class
i_do_sigt=1

# sea-ice diags
i_do_ice=1  ; # Sea-ice diags
export FILE_ICE_SUFFIX="icemod" ; # in what file to find ice fraction NN_ICEF? => "icemod" or "grid_T"


i_do_bb=1   ; # Budget and other stuffs on a given rectangular box!
#             # => needs file FILE_DEF_BOXES !!!
# => produces time-series f(t)  (mean of 2D fields)


i_do_ssx_box=0 ; # zoom on given boxes (+spatially-averaged values) for surface properties
#                # boxes defined into barakuda_orca.py ...


# Vertical profiles on of box-averaged as a function of time...
i_do_box_TS_z=1 ; # do sigma vert. profiles on given boxes... # 1 => no figures, 2 => figures
#                 # => needs file FILE_DEF_BOXES !!!
# => produces time-series f(t,z)

#
# Deep Mixed volume in prescribed boxes:
i_do_dmv=1
export MLD_CRIT="1000,725,500"


# Some nerdy stuffs about the critical depth in prescribed boxes: 
i_do_zcrit=0






# BETA / TESTING:

# Fresh-water transport associated to sea-ice transport
#  => must compile cdficeflux.x but depends on more recent CDFTOOLS module...
i_do_icet=0 ; # treat sea-ice volume transport!
export TRANSPORT_ICE_SECTION_FILE="${BARAKUDA_ROOT}/data/transportiz_ORCA1_ARCTIC.dat"




i_do_flx=0  ; # surface fluxes diags




i_do_amo=0 ;  # buit a SST time serie usable to build Atlantic Multidecadal Oscilation index


i_do_sect=0 ; # do sigma vert. profiles on given boxes...
VSECT_NM=( "Indian_77p5_E" "Atlantic_21p5_W" )
VSECT_JI=(      "5,5"          "266,266"     ) ; # X range in C convention
VSECT_JJ=(    "25,170"          "7,291"      ) ; # Y range in C convention
