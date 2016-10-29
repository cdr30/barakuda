#!/bin/bash

#==========================================================
#
#         Configuration file for
#
# OCEAN MONITORING for NEMO v3.6 ORCA2 on 31 levels
#
#            HPC: triolith
#
#        L. Brodeau, 2015
#
#===========================================================

export CONF=ORCA2 ; # horizontal global configuration
export NBL=31     ; # number of levels

# File system where NEMO config and output files are stored:
export STORE_DIR="/proj/bolinc/users/x_laubr"


# Is it an ec-earth run?
export ece_run=0 ; # means that NEMO files in something like ${STORE_DIR}/<RUN>/output/nemo/<YYY>
#                  # where YYY starts from '001' to
export Y_INI_EC=1990 ;    # initial year if ec-earth run...

# List of suffixed of files that have been saved by NEMO and that are needed for the diags:
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

export PYTHON_HOME=${HOME}/opt/Canopy_64bit/User

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

export FILE_ICE_SUFFIX="icemod" ; # in what file to find ice fraction and volume?
export NN_ICEF="siconc" ; # name of ice fraction in "FILE_ICE_SUFFIX" file...
export NN_ICET="sivolu" ; # ice thickness or rather volume...

export FILE_FLX_SUFFIX="SBC" ; # in what file to find surface fluxes ?
export NN_FWF="wfo"       ; # name of net freshwater flux (E-P-R) in "FILE_FLX_SUFFIX" file...
export NN_EMP="emp_oce"   ; # name of E-P in "FILE_FLX_SUFFIX" file...
export NN_P="precip"   ; # name of P in "FILE_FLX_SUFFIX" file...
export NN_RNF="XXX"          ; # name of continental runoffs in "FILE_FLX_SUFFIX" file...
export NN_CLV="XXX"  ; # calving from icebergs in "FILE_FLX_SUFFIX" file...
export NN_E="XXX"           ; # evaporation in "FILE_FLX_SUFFIX" file...


export L_CONV2NC3=false ; # Set to true if your NEMO output is in Netcdf4 and your NCO does not support netcdf4!

export L_RENAME=false ; # set to true if your ORCA output has old name convention (ex: votemper instead of thetao)


export JTITLE="NEMO v3.6 ${CONF} (L${NBL}) - LIM3 / ocean-only experiment"

# Land-sea mask and basins files:
export MM_FILE="/proj/bolinc/users/x_laubr/${CONF}/mesh_mask_${CONF}_20151014.nc"
export BM_FILE="/proj/bolinc/users/x_laubr/${CONF}/basin_mask_${CONF}.nc"

# 3D monthly climatologies of potential temperature and salinity (can be those you used for the NEMO run):
export F_T_CLIM_3D_12=${STORE_DIR}/${CONF}/${CONF}-I/data_1m_potential_temperature_nomask_${CONF}.nc
export F_S_CLIM_3D_12=${STORE_DIR}/${CONF}/${CONF}-I/data_1m_salinity_nomask_${CONF}.nc
export SST_CLIM_12=${STORE_DIR}/${CONF}/${CONF}-I/sst_data_${CONF}.nc
export NN_T_CLIM="votemper"
export NN_S_CLIM="vosaline"
export NN_SST_CLIM="sst"

export ICE_CLIM_12="${STORE_DIR}/${CONF}/${CONF}-I/ice_cover_180x360-${CONF}_Hurrell_monthly_mean1980-1999.nc4"
export NN_ICEF_CLIM="ice_cover"


# A text file where the vertical hydraugraphical sections of interest are defined :
#export TRANSPORT_SECTION_FILE="${BARAKUDA_ROOT}/data/transportiz_${CONF}_light.dat"
export TRANSPORT_SECTION_FILE="${BARAKUDA_ROOT}/data/transportiz_${CONF}_all.dat"

# For transport by sigma-class:
export DENSITY_SECTION_FILE="${BARAKUDA_ROOT}/data/dens_section_${CONF}.dat"


# In what directory of the local machine to save the diagnostics:
export DIAG_DIR="/proj/bolinc/users/x_laubr/tmp/barakuda/${CONF}_v36"


# Files with the list of rectangular boxes to look at more closely:
export FILE_DEF_BOXES="${BARAKUDA_ROOT}/data/def_boxes_convection_${CONF}.txt"



# About remote HOST to install HTML pages to:
ihttp=1 ; # do we export on a remote http server (1) or keep on the local machine (0)
RHOST=misu228.misu.su.se ; # remote host to send diagnostic page to///
RUSER=laurent ; # username associated to remote host (for file export)
RWWWD=/data/www/barakuda/${CONF} ; # directory of the local or remote host to send the diagnostic page to




#########################
# Diags to be performed #
#########################

# In what format should figures be produced:
export FIG_FORM="png"



# Basic 3D and surface averages:
i_do_mean=1

# FreshWater fluxes at the surface spatially averaged over the ocean, E-P-R, E-P, R, P, ...
#i_do_fwf=1

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
i_do_dmv=0
export MLD_CRIT="1000,725,500"


# Some nerdy stuffs about the critical depth in prescribed boxes:
i_do_zcrit=0






# BETA / TESTING:

# Fresh-water transport associated to sea-ice transport
#  => must compile cdficeflux.x but depends on more recent CDFTOOLS module...
i_do_icet=0 ; # treat sea-ice volume transport!
export TRANSPORT_ICE_SECTION_FILE="${BARAKUDA_ROOT}/data/transportiz_${CONF}_ARCTIC.dat"




i_do_flx=0  ; # surface fluxes diags




i_do_amo=0 ;  # buit a SST time serie usable to build Atlantic Multidecadal Oscilation index


i_do_sect=0 ; # do sigma vert. profiles on given boxes...
VSECT_NM=( "Indian_77p5_E" "Atlantic_21p5_W" )
VSECT_JI=(      "5,5"          "266,266"     ) ; # X range in C convention
VSECT_JJ=(    "25,170"          "7,291"      ) ; # Y range in C convention
