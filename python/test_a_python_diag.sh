#!/bin/bash

issh=0
irnf=0
iice=1
iemp=0
icmip5=0
ihov=0



export STORE_DIR="/proj/bolinc/users/x_laubr"


HERE=`pwd`

if [ ${ihov} -eq 1 ]; then
    export RUN=cp70
    export ORCA=ORCA1.L75
    export DIAG_D=/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${RUN}
    export MM_FILE=/proj/bolinc/users/x_laubr/klaus/mesh_mask.nc
    export BM_FILE=/proj/bolinc/users/x_laubr/klaus/basin_mask.nc
    export NN_T="thetao"
    export NN_S="so"
    #
    cd ${DIAG_D}/
    python /home/x_laubr/DEV/barakuda/python/exec/plot_hovm_tz.py 1996 2000
    
    mv -f hov_*_ORCA1.L75-${RUN}*.png ${HERE}/
    #
fi





if [ ${iemp} -eq 1 ]; then
    
    export ORCA="ORCA1.L75"
    #export RUN="32bI"
    export RUN="cp00"
    export TSTAMP="1m"
    export DIAG_D="/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${RUN}"
    #export NN_RNF="runoffs"
    export MM_FILE="/proj/bolinc/users/x_laubr/tmp/barakuda/test/mesh_mask.nc"
    export TRANSPORT_SECTION_FILE="boo"
    export LMOCLAT="boo" ; export NN_SSH="boo" ; export NN_SSS="boo" ; export NN_S="boo"
    export NN_MLD="boo" ; export NN_SST="boo" ; export NN_T="boo"
    export NN_FWF="wfo"       ; # name of net freshwater flux (E-P-R) in "FILE_FLX_SUFFIX" file...               
    export NN_EMP="emp_oce"   ; # name of E-P in "FILE_FLX_SUFFIX" file...                                       
    export NN_P="precip"   ; # name of P in "FILE_FLX_SUFFIX" file...                                            
    export NN_RNF="XXX"   ; # name of continental runoffs in "FILE_FLX_SUFFIX" file...



    export FILE_DEF_BOXES="/home/x_laubr/DEV/barakuda/data/def_boxes_convection_ORCA1.txt"

    cd ${DIAG_D}/

    python /home/x_laubr/DEV/barakuda/python/exec/plot_time_series.py mean_fwf



    

fi





if [ ${icmip5} -eq 1 ]; then
    export RUN="SPIN"
    export CPREF="ORCA1-${RUN}_MM_"
    export ORCA=ORCA1.L42 ; # horizontal global configuration                             
    export NBL=42         ; # number of levels                                            
    export STORE_DIR="/proj/bolinc/users/x_laubr"
    export TSTAMP="MM"
    export DIAG_D="."
    export MM_FILE="${STORE_DIR}/${ORCA}/${ORCA}-I/mesh_mask_ORCA1_ecearth2_42l_NoCaspian.nc4"
    export NN_SST="sosstsst"
    export NN_SSS="sosaline"
    export NN_SSH="sossheig"
    export NN_T="votemper"
    export NN_S="vosaline"    

    export NN_TAUX="sozotaux" 
    export NN_TAUY="sometauy"

    export FILE_DEF_BOXES="/home/x_laubr/DEV/barakuda/data/def_boxes_convection_ORCA1.txt"

    python /home/x_laubr/DEV/barakuda/python/exec/budget_rectangle_box.py 2250 100 uv

fi 






if [ ${irnf} -eq 1 ]; then

    export ORCA="ORCA1.L75"
    export RUN="LB03"
    export DIAG_D="/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${RUN}"
    export NN_RNF="runoffs"
    export MM_FILE="/proj/bolinc/users/x_laubr/tmp/barakuda/test/mesh_mask.nc"

    python exec/runoffs.py 1997 1999

fi

if [ ${iice} -eq 1 ]; then
    export RUN=cp70
    export ORCA=ORCA1.L75
    export DIAG_D=/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${RUN}
    export MM_FILE=/proj/bolinc/users/x_laubr/klaus/mesh_mask.nc
    export BM_FILE=/proj/bolinc/users/x_laubr/klaus/basin_mask.nc

    export COMP2D="CLIM"
    export FILE_ICE_SUFFIX="icemod"
    export NN_ICEF="siconc" ; # name of ice fraction in "FILE_ICE_SUFFIX" file...                                
    export NN_ICET="sivolu" ; # ice thickness or rather volume...
    #export NN_ICEF="iiceconc" ; # name of ice fraction in "FILE_ICE_SUFFIX" file...
    #export NN_ICET="iicethic" ; # ice thickness but 'sit' is only in icemod file !!!

    export ICE_CLIM_12=${STORE_DIR}/ORCA1.L75/ORCA1.L75-I/ice_cover_180x360-ORCA1_Hurrell_monthly_mean1980-1999.nc4
    export NN_ICEF_CLIM="ice_cover"

    python exec/ice.py 1996 2000

fi


if [ ${issh} -eq 1 ]; then
    
    #export CONF="ORCA1.L46"
    #export ORCA="ORCA1.L46"
    #export RUN="LCM2"
    #export DIAG_D="/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L46_ece31/ORCA1.L46-LCM2"
    #export NN_SSH="sossheig"
    #export MM_FILE="/proj/bolinc/users/x_laubr/${CONF}/mesh_mask.nc"
    #python exec/ssh.py 1024 1028
    
    export CONF="ORCA1.L75"
    export ORCA="ORCA1.L75"
    #export RUN="LCM2"
    #export DIAG_D="/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece31/ORCA1.L75-LCM2"
    #export NN_SSH="sossheig"
    #export MM_FILE="/proj/bolinc/users/x_laubr/${CONF}/mesh_mask.nc"
    #python exec/ssh.py 1024 1028

fi

