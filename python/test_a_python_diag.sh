#!/bin/bash

#  B E T A  ! ! !

# Diag to test:
ifwf=0
imov=0
issh=0
its=1
imld=0
irnf=0
iice=0
iemp=0
icmip5=0
ihov=0

CONFIG="ORCA1_L75"
ARCH="ece32_marenostrum"

export RUN="LB01"

jyear=1990 ; # year to test on (if relevant)


export BARAKUDA_ROOT=`pwd | sed -e "s|/python||g"`




. ${BARAKUDA_ROOT}/src/bash/bash_functions.bash
. ${BARAKUDA_ROOT}/configs/config_${CONFIG}_${ARCH}.sh

ORCA_LIST="ORCA1.L75 ORCA1.L46 ORCA1.L42 ORCA2 ORCA2_L46"

for og in ${ORCA_LIST}; do
    ca=""; ca=`echo ${CONFIG} | grep ${og}` ; if [ "${ca}" != "" ]; then export ORCA=${og}; fi
done
if [ "${ORCA}" = "" ]; then echo "ORCA grid of config ${CONFIG} not supported yet"; exit; fi

export CONFRUN=${ORCA}-${RUN}
export DIAG_D=${DIAG_DIR}/${CONFRUN}

HERE=`pwd`

finfoclim=${DIAG_D}/clim/last_clim

y1_clim=`cat ${finfoclim} | cut -d - -f1`
y2_clim=`cat ${finfoclim} | cut -d - -f2`

export COMP2D="CLIM"

# To know the name of NEMO output files:
export NEMO_OUT_D=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g"`
if [ ! -d ${NEMO_OUT_D} ]; then echo "Unfortunately we could not find ${NEMO_OUT_D}"; exit; fi
YEAR_INI=1990 ; YEAR_INI_F=1990
export cyear=`printf "%04d" ${jyear}`
if [ ${ece_run} -gt 0 ]; then
    iy=$((${jyear}-${YEAR_INI}+1+${YEAR_INI}-${YEAR_INI_F}))
    dir_ece="`printf "%03d" ${iy}`/"
fi
CPREF=`echo ${NEMO_FILE_PREFIX} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g" -e "s|<TSTAMP>|${TSTAMP}|g"`
ft=${NEMO_OUT_D}/${dir_ece}${CPREF}${cyear}0101_${cyear}1231_grid_T.nc4
check_if_file ${ft}
fj=${NEMO_OUT_D}/${dir_ece}${CPREF}${cyear}0101_${cyear}1231_icemod.nc4
check_if_file ${fj}



export PYTHONPATH=${PYTHON_HOME}/lib/python2.7/site-packages:${BARAKUDA_ROOT}/python/modules

echo ; echo " *** DIAG_D=${DIAG_D} !"; echo


rm -f *.png *.nc

# Time for diags:

if [ ${its} -eq 1 ]; then
    #diag=3d_thetao
    diag=mean_fwf
    ln -sf ${DIAG_D}/${diag}*.nc .
    CMD="python exec/plot_time_series.py ${diag}"

fi



if [ ${ifwf} -eq 1 ]; then
    CMD="${BARAKUDA_ROOT}/src/scripts/do_fwf_series_ifs.sh"
fi


if [ ${imov} -eq 1 ]; then
    CMD="python exec/prepare_movies.py ${ft} ${jyear} sss"
    #CMD="python exec/prepare_movies.py ${fj} ${jyear} ice"
fi

if [ ${issh} -eq 1 ]; then
    CMD="python exec/ssh.py ${y1_clim} ${y2_clim}"
fi

if [ ${imld} -eq 1 ]; then
    CMD="python exec/mld.py ${y1_clim} ${y2_clim}"
fi



echo
echo "DOING: ${CMD}"
${CMD}


# Add other diags here:






exit
# BELOW = OLD STUFFS, fix!










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



