#!/usr/bin/env bash

#==============================================================
#
#                    B A R A K U D A
#
#    An OCEAN MONITORING python environment for NEMO
#
#             L. Brodeau, 2009-2016
#
#===============================================================

iuv=0      ; # Do a climatology for current...
ivt=1      ; # Do a climatology for VT
iamoc=1    ; # Do a climatology for 2D lat-depth AMOC?
ibpsi=0    ; # Do a climatology for barotropic stream function

export BARAKUDA_ROOT=`pwd`

# Checking available configs
list_conf=`\ls configs/config_*.sh` ; list_conf=`echo ${list_conf} | sed -e s/'configs\/config_'/''/g -e s/'.sh'/''/g`

# Important bash functions:
. ${BARAKUDA_ROOT}/src/bash/bash_functions.bash

barakuda_init

while getopts R:f:i:e:C:h option ; do
    case $option in
        R) export RUN=${OPTARG};;
        f) export IFREQ_SAV_YEARS=${OPTARG} ;;
        i) export Y1=${OPTARG} ;;
        e) export Y2=${OPTARG} ;;
        C) export CONFIG=${OPTARG};;
        h)  usage;;
        \?) usage ;;
    esac
done

barakuda_check

# sourcing configuration file
fconfig=${BARAKUDA_ROOT}/configs/config_${CONFIG}.sh

if [ -f ${fconfig} ]; then
    echo "Sourcing ${fconfig} !"
    . ${fconfig}
else
    echo "PROBLEM: cannot find file ${fconfig} !"; exit
fi

barakuda_setup

echo
echo " SETTINGS: "
echo "   *** DIAG_D   = ${DIAG_D} "
echo "   *** CLIM_DIR = ${CLIM_DIR} "
echo "   *** TMP_DIR  = ${TMP_DIR} "
echo "   *** Y1       = ${Y1} "
echo "   *** Y2       = ${Y2} "
echo "   *** CONFIG   = ${CONFIG} "
echo "   *** GRID     = ${ORCA} "
echo "   *** RUN      = ${RUN} "
echo "   *** CPREF    = ${CPREF} "
echo "   *** IFREQ_SAV_YEARS = ${IFREQ_SAV_YEARS} "
echo

CP2NC4="${NCDF_DIR}/bin/nccopy -k 4 -d 9"

Y1=$((${Y1}+0))
Y2=$((${Y2}+0))
CY1=`printf "%04d" ${Y1}`
CY2=`printf "%04d" ${Y2}`


mkdir -p ${CLIM_DIR}


# Variables to extract:
V2E="${NN_SST},${NN_SSS},${NN_SSH},${NN_T},${NN_S},${NN_MLD}"
#,,qns,wfo,taum,ice_cover,qsb_oce,qla_oce,qlw_oce"
#V2E_lev="thetao,so"

C2ET="nav_lon,nav_lat,deptht" #,time_counter"
C2EU="nav_lon,nav_lat,depthu" #,time_counter"
C2EV="nav_lon,nav_lat,depthv" #,time_counter"
C2EW="nav_lon,nav_lat,depthw" #,time_counter"

GRID_IMP="grid_T"

if [ ${ivt} -eq 1 -o ${ibpsi} -eq 1 -o ${iuv} -eq 1 ]; then
    GRID_IMP="${GRID_IMP} grid_U"
fi

if [ ${iamoc} -eq 1 -o ${ivt} -eq 1 -o ${ibpsi} -eq 1 -o ${iuv} -eq 1 ]; then
    GRID_IMP="${GRID_IMP} grid_V"
fi

if [ `contains_string ${FILE_ICE_SUFFIX} ${NEMO_SAVED_FILES}` -eq 1 ]; then
    GRID_IMP="${GRID_IMP} ${FILE_ICE_SUFFIX}"
fi

if [ `contains_string SBC ${NEMO_SAVED_FILES}` -eq 1 ]; then
    GRID_IMP="${GRID_IMP} SBC"
fi

echo; echo " GRID_IMP = ${GRID_IMP}"; echo



# Checking what files we have / plan to use:
if [ "${NEMO_SAVED_FILES}" = "" ]; then
    echo "Please specify which NEMO files are saved (file suffixes, grid_T, ..., icemod) ?"
    echo " => set the variable NEMO_SAVED_FILES in your config_${CONFIG}.sh file!"; exit
fi
VAF=( "grid_T" "grid_U" "grid_V" "icemod" "SBC" )
js=0 ; gimp_new=""
for sf in ${VAF[*]}; do
    echo "Checking ${sf}..."
    ca=`echo "${NEMO_SAVED_FILES}" | grep ${sf}`; #if [ "${ca}" = "" ]; then VOK[${js}]=0; fi
    cb=`echo "${GRID_IMP}"         | grep ${sf}`
    if [ "${ca}" = "" ]; then
        if [ "${cb}" != "" ]; then
            echo "PROBLEM! The diags you specified say you need ${sf} files"
            echo "     => but you have not specified ${sf} in NEMO_SAVED_FILES !"; exit
        fi
    else
        gimp_new="${sf} ${gimp_new}"
    fi
    ((js++))
done
GRID_IMP=${gimp_new}
echo; echo "File types to import: ${GRID_IMP}"; echo; echo


VCM=( "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" )

cd ${TMP_DIR}/

echo; echo "In:"; pwd


barakuda_import_mesh_mask


ls -l ; echo; echo


nby=$((${Y2}-${Y1}+1))

if [ ! $((${nby}%${IFREQ_SAV_YEARS})) -eq 0 ]; then
    echo " Number of years should be a multiple of ${IFREQ_SAV_YEARS}!"; exit
fi



if [ ${ece_run} -gt 0 ]; then
    if [ ! -d ${NEMO_OUT_D}/001 ]; then echo "ERROR: since ece_run=${ece_run}, there should be a directory 001 in:"; echo " ${NEMO_OUT_D}"; fi
fi


ffirsty="${DIAG_D}/first_year.info"
if [ ! -f  ${ffirsty} ]; then echo "ERROR: file ${ffirsty} not found!!!"; exit; fi
export YEAR_INI_F=`cat ${ffirsty}`

export jyear=${Y1}

while [ ${jyear} -le ${Y2} ]; do

    export cyear=`printf "%04d" ${jyear}` ; echo ; echo "Year = ${cyear}"

    cpf=""
    if [ ${ece_run} -gt 0 ]; then
        iy=$((${jyear}-${Y1}+1+${Y1}-${YEAR_INI_F}))
        dir_ece=`printf "%03d" ${iy}`
        echo " *** ${cyear} => dir_ece = ${dir_ece}"
        cpf="${dir_ece}/"
    fi

    TTAG_ann=${cyear}0101_${cyear}1231

    i_get_file=0
    if [ $((${jyear}%${IFREQ_SAV_YEARS})) -eq 0 ]; then
        barakuda_check_year_is_complete  ; # lcontinue might be updated to false!
    fi

    CRTM=${CPREF}${TTAG}
    CRT1=${CPREF}${TTAG_ann}

    barakuda_import_files

    # Testing if ALL required files are present now:
    for gt in ${GRID_IMP}; do
        ftt="./${CRT1}_${gt}.nc" ;  check_if_file ${ftt}
    done
    echo; echo "All required files are in `pwd` for year ${cyear} !"; echo

    # Files to work with for current year:
    ft=${CRT1}_grid_T.nc
    fu=${CRT1}_grid_U.nc
    fv=${CRT1}_grid_V.nc
    #fg=${CRT1}_${FILE_ICE_SUFFIX}.nc ; # can be icemod or grid_T ....
    #fvt=${CRT1}_VT.nc

    echo

    if [ ${ivt} -eq 1 ]; then
        # Creating VT files:
        echo "Doing: ${BARAKUDA_ROOT}/cdftools_light/bin/cdfvT.x ${CPREF}${TTAG_ann} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV}"
        ${BARAKUDA_ROOT}/cdftools_light/bin/cdfvT.x ${CPREF}${TTAG_ann} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} &
        echo
    fi

    echo

    if [ ${iamoc} -eq 1 ]; then
        rm -f moc.nc
        echo " *** doing: ./cdfmoc.x ${fv} ${NN_V} ${NN_V_EIV}"
        ${BARAKUDA_ROOT}/cdftools_light/bin/cdfmoc.x ${fv} ${NN_V} ${NN_V_EIV} &
        echo
    fi

    if [ ${ibpsi} -eq 1 ]; then
        rm -f psi.nc
        echo " *** doing: ./cdfpsi.x ${fu} ${fv} V"
        ${BARAKUDA_ROOT}/cdftools_light/bin/cdfpsi.x ${fu} ${fv} V &
        echo
    fi

    wait

    ncwa -O -a x moc.nc ${CPREF}${TTAG_ann}_MOC.nc ; # removing degenerate x record...
    rm -f moc.nc

    if [ ${ibpsi} -eq 1 ]; then mv -f psi.nc ${CPREF}${TTAG_ann}_PSI.nc; fi

    echo "After year ${jyear}:"; ls -l *.nc*
    echo

    ((jyear++))
    export jyear=${jyear}

done



if [ ${ivt} -eq 1 ]; then
    fo=aclim_${CONFRUN}_${CY1}-${CY2}_VT.nc
    # Averaged VT file:
    ncra -O ${CPREF}*_VT.nc -o ${fo}
    rm -f ${CPREF}*_VT.nc

    # Converting to netcdf4 with maximum compression level:
    ${CP2NC4} ${fo} ${fo}4 &

fi

if [ ${iamoc} -eq 1 ]; then
    fo=aclim_${CONFRUN}_${CY1}-${CY2}_MOC.nc
    # Averaged MOC file:
    ncra -O ${CPREF}*_MOC.nc -o ${fo}
    # Converting to netcdf4 with maximum compression level:
    ${CP2NC4} ${fo} ${fo}4 &
fi


if [ ${ibpsi} -eq 1 ]; then
    fo=aclim_${CONFRUN}_${CY1}-${CY2}_PSI.nc
    # Averaged PSI file:
    ncra -O ${CPREF}*_PSI.nc -o ${fo}
    # Converting to netcdf4 with maximum compression level:
    ${CP2NC4} ${fo} ${fo}4 &
fi


echo;echo;echo;






echo "Phase 2:"; ls ; echo



# Mean monthly climatology

for suff in grid_T grid_U grid_V icemod SBC VT MOC PSI; do


    if [ -f ./${CPREF}${CY1}0101_${CY1}1231_${suff}.nc ]; then

        echo ; echo ; echo ; echo " Treating ${suff} files!"; echo

        f2c=mclim_${CONFRUN}_${CY1}-${CY2}_${suff}.nc
        f2c_reg=mclim_${CONFRUN}_${CY1}-${CY2}_${REGG}.nc

        rm -f ${CLIM_DIR}/${f2c}* ${CLIM_DIR}/${f2c_reg}*

        echo

        jm=0
        for cm in ${VCM[*]}; do

            ((jm++))

            if [ -f ./${CPREF}${CY1}0101_${CY1}1231_${suff}.nc ]; then
                echo; ls ; echo
                echo "ncra -F -O -d time_counter,${jm},,12 ${CPREF}*0101_*1231_${suff}.nc -o mean_m${cm}_${suff}.nc"
                ncra -F -O -d time_counter,${jm},,12 ${CPREF}*0101_*1231_${suff}.nc -o mean_m${cm}_${suff}.nc
                echo
            fi

        done

        rm -f ${CPREF}*0101_*1231_${suff}.nc

        echo; ls ; echo
        echo "ncrcat -O  mean_m*_${suff}.nc out_${suff}.nc"
        ncrcat -O  mean_m*_${suff}.nc out_${suff}.nc
        rm mean_m*_${suff}.nc
        echo

        echo "mv -f out_${suff}.nc ${f2c}"
        mv -f out_${suff}.nc ${f2c}
        echo

        # Converting to netcdf4 with maximum compression level:
        echo; ls ; echo
        echo "${CP2NC4} ${f2c} ${f2c}4"
        ${CP2NC4} ${f2c} ${f2c}4 &
        echo

    else
        echo ; echo ; echo ; echo " Ignoring monthly ${suff} files!"; echo
    fi

done ; # loop along files suffixes

wait
wait

for cl in aclim mclim; do
    echo "mv -f ${cl}_${CONFRUN}*.nc4 ${CLIM_DIR}/"
    mv -f ${cl}_${CONFRUN}*.nc4 ${CLIM_DIR}/
    echo
done

echo;echo
echo "${CY1}-${CY2}" > ${CLIM_DIR}/last_clim
echo "Climatology saved into: ${CLIM_DIR}/"
echo;echo

cd /tmp/
rm -rf ${TMP_DIR} 2>/dev/null
exit
