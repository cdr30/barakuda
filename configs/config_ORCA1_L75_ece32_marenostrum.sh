#!/bin/bash

#==========================================================
#
#         Configuration file for
#
# OCEAN MONITORING for NEMO v3.6 of EC-Earth 3.2 beta tunning on 75 levels
#
#            HPC: triolith
#
#        L. Brodeau, 2015
#
#===========================================================

export HOST=MARENOSTRUM

export CONF=ORCA1.L75 ; # horizontal global configuration
export NBL=75     ; # number of levels

# Path to directory containing NEMO output files:
export STORE_DIR="/gpfs/scratch/bsc32/bsc32325"

# Path to directory containing some 2D and 3D climatologies on the relevant ORCA grid:
export CONF_INI_DIR="/gpfs/projects/bsc32/bsc32325/ORCA1/ORCA1-I/barakuda_clim"

# In what directory of the local machine to save the diagnostics:
export DIAG_DIR="${STORE_DIR}/barakuda/${CONF}_ece32"

module add NCO/4.2.3
module add PYTHON/2.7.3

# Is it an ec-earth run?
export ece_run=2 ; # 0 => not an EC-Earth run, it's a "pure" ocean-only NEMO run done from traditional NEMO setup
#                  # 1 => it's an OCEAN-ONLY EC-Earth run done from a EC-Earth setup
#                  # 2 => it's a  COUPLED  EC-Earth run
#                  #      Both 1 and 2 imply that NEMO files are stored in something like
#                  #       ${STORE_DIR}/<RUN>/output/nemo/<YYY>
#                  #       where YYY starts from '001' to
#                  #   If you select '2', make sure 'cdo' is available and working!!!
#
export Y_INI_EC=1990 ;    # initial year if ec-earth run...


# List of suffixed of files that have been saved by NEMO and that are needed for the diags:
export NEMO_SAVED_FILES="grid_T grid_U grid_V icemod SBC"

# Directory structure in which to find NEMO output file (use <ORCA> and <RUN>):
export NEMO_OUT_STRCT="${STORE_DIR}/run_ece/<RUN>/output/nemo"

export TSTAMP="1m"   ; # output time-frequency stamp as in NEMO output files...

# How does the nemo files prefix looks like
# Everything before "<year_related_info>_grid_<X>" or "<year_related_info>_icemod"
# use <ORCA>, <RUN> and <TSTAMP>=>  Ex: export NEMO_FILE_PREFIX="<ORCA>-<RUN>_<TSTAMP>_"
export NEMO_FILE_PREFIX="<RUN>_<TSTAMP>_"
# => should get rid of TSTAMP actually...

# Temporary file system (scratch) on which to perform the job you can use <JOB_ID> if scracth depends on JOB ID:
#export SCRATCH="/scratch/local/<JOB_ID>" ; # triolith
export SCRATCH="/scratch/tmp"

export PYTHON_HOME="/apps/PYTHON/2.7.3"

# If variables names in NEMO files are not the default...
export NN_SST="sosstsst"
export NN_SSS="sosaline"
export NN_SSH="sossheig"
export NN_T="votemper"
export NN_S="vosaline"
export NN_MLD="mldr10_1"
export NN_U="vozocrtx"
export NN_V="vomecrty"
#export NN_U_EIV="vozoeivu"
#export NN_V_EIV="vomeeivv"
export NN_U_EIV="0" ; # ignore
export NN_V_EIV="0" ; # ignore
export NN_TAUX="sozotaux"
export NN_TAUY="sometauy"

export FILE_ICE_SUFFIX="icemod" ; # in what file to find ice fraction and volume?
export NN_ICEF="siconc" ; # name of ice fraction in "FILE_ICE_SUFFIX" file...
export NN_ICET="sivolu" ; # ice thickness or rather volume...

export FILE_FLX_SUFFIX="SBC" ; # in what file to find surface fluxes ?
export NN_FWF="wfo"          ; # name of net freshwater flux (E-P-R) in "FILE_FLX_SUFFIX" file...
export NN_EMP="emp_oce"      ; # name of E-P in "FILE_FLX_SUFFIX" file...
export NN_P="precip"         ; # name of P in "FILE_FLX_SUFFIX" file...
export NN_RNF="runoffs"      ; # name of continental runoffs in "FILE_FLX_SUFFIX" file...
export NN_CLV="calving_cea"  ; # calving from icebergs in "FILE_FLX_SUFFIX" file...
export NN_E="evap"           ; # evaporation in "FILE_FLX_SUFFIX" file...


export L_CONV2NC3=false ; # Set to true if your NEMO output is in Netcdf4 and your NCO does not support netcdf4!

export L_RENAME=false ; # set to true if your ORCA output has old name convention (ex: votemper instead of thetao)


export JTITLE="LIM3, NEMO 3.6 (EC-Earth 3.2b_tuning)"

# Land-sea mask and basins files:
export MM_FILE=/gpfs/projects/bsc32/bsc32325/ORCA1/ec-earth3.2/mesh_mask.nc4
export BM_FILE=/gpfs/projects/bsc32/bsc32325/ORCA1/ec-earth3.2/basin_mask.nc4

# 3D monthly climatologies of potential temperature and salinity (can be those you used for the NEMO run):
export F_T_CLIM_3D_12=${CONF_INI_DIR}/thetao_1degx1deg-ORCA1.L75_WOA2009_monthly_LB_20160223.nc4
export F_S_CLIM_3D_12=${CONF_INI_DIR}/so_1degx1deg-ORCA1.L75_WOA2009_monthly_LB_20160223.nc4
export SST_CLIM_12=${CONF_INI_DIR}/tos_180x360-ORCA1_Reynolds_monthly_mean1982-2005.nc4
export NN_T_CLIM="thetao"
export NN_S_CLIM="so"
export NN_SST_CLIM="tos"

export ICE_CLIM_12=${CONF_INI_DIR}/ice_cover_180x360-ORCA1_Hurrell_monthly_mean1980-1999.nc4
export NN_ICEF_CLIM="ice_cover"


# A text file where the vertical hydraugraphical sections of interest are defined :
export TRANSPORT_SECTION_FILE="${BARAKUDA_ROOT}/data/transportiz_ORCA1.dat"

# For transport by sigma-class:
export DENSITY_SECTION_FILE="${BARAKUDA_ROOT}/data/dens_section_ORCA1.dat"



# Files with the list of rectangular boxes to look at more closely:
export FILE_DEF_BOXES="${BARAKUDA_ROOT}/data/def_boxes_convection_ORCA1.txt"
export FILE_DMV_BOXES="${BARAKUDA_ROOT}/data/def_boxes_convection_ORCA1.txt"



# About remote HOST to install HTML pages to:
ihttp=0 ; # do we export on a remote http server (1) or keep on the local machine (0)
RHOST=misu228.misu.su.se ; # remote host to send diagnostic page to///
RUSER=laurent ; # username associated to remote host (for file export)
RWWWD=/data/www/barakuda/ec-earth_3.2b ; # directory of the local or remote host to send the diagnostic page to




#########################
# Diags to be performed #
#########################

# In what format should figures be produced:
export FIG_FORM="png"


# Movies of SST and SSS compared to OBS:
export i_do_movi=1

# Basic 3D and surface averages:
export i_do_mean=1

# IFS freshWater fluxes at the surface spatially averaged over the ocean, E-P-R, E-P, R, P, ...
export i_do_fwf=1 ; # only relevant when ece_run=2...

# AMOC:
export i_do_amoc=1
export LMOCLAT="20-23 30-33 40-43 45-48 50-53" ; # List of latitude bands to look in for max of AMOC


# Transport of mass, heat and salt through specified sections (into TRANSPORT_SECTION_FILE):
export i_do_trsp=2  ; # transport of mass, heat and salt through specified sections
#              # i_do_trsp=2 => treat also different depths range!
z1_trsp=100  ; # first  depth: i_do_trsp must be set to 2
z2_trsp=1000 ; # second depth: i_do_trsp must be set to 2


# meridional heat/salt transport (advective)
export i_do_mht=1

# Transport by sigma class
export i_do_sigt=1

# sea-ice diags
export i_do_ice=1  ; # Sea-ice diags


export i_do_bb=1   ; # Budget and other stuffs on a given rectangular box!
#             # => needs file FILE_DEF_BOXES !!!
# => produces time-series f(t)  (mean of 2D fields)


export i_do_ssx_box=0 ; # zoom on given boxes (+spatially-averaged values) for surface properties
#                # boxes defined into barakuda_orca.py ...


# Vertical profiles on of box-averaged as a function of time...
export i_do_box_TS_z=1 ; # do sigma vert. profiles on given boxes... # 1 => no figures, 2 => figures
#                 # => needs file FILE_DEF_BOXES !!!
# => produces time-series f(t,z)

#
# Deep Mixed volume in prescribed boxes:
export i_do_dmv=0
export MLD_CRIT="1000,725,500"


# Some nerdy stuffs about the critical depth in prescribed boxes:
export i_do_zcrit=0






# BETA / TESTING:

# Fresh-water transport associated to sea-ice transport
#  => must compile cdficeflux.x but depends on more recent CDFTOOLS module...
export i_do_icet=0 ; # treat sea-ice volume transport!
export TRANSPORT_ICE_SECTION_FILE="${BARAKUDA_ROOT}/data/transportiz_ORCA1_ARCTIC.dat"




export i_do_flx=0  ; # surface fluxes diags




export i_do_amo=0 ;  # buit a SST time serie usable to build Atlantic Multidecadal Oscilation index


export i_do_sect=0 ; # do sigma vert. profiles on given boxes...
VSECT_NM=( "Indian_77p5_E" "Atlantic_21p5_W" )
VSECT_JI=(      "5,5"          "266,266"     ) ; # X range in C convention
VSECT_JJ=(    "25,170"          "7,291"      ) ; # Y range in C convention


#========================== Marenostrum ================================================================
### Shouldn't be needed elsewhere than MareNostrum, where it's a hello to have CDO working...
## => Only if you specified ece_run=2 and i_do_fwf
export MOD_CDO="gcc/4.7.2 intel/13.0.1 openmpi/1.8.1 NETCDF/4.1.3 HDF5/1.8.10 UDUNITS/2.1.24 CDO/1.7.0"
#=======================================================================================================



