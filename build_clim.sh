#!/bin/bash

#==============================================================
#
#                    B A R A K U D A
#
#    An OCEAN MONITORING python environment for NEMO
#
#             L. Brodeau, August 2009-2015
#
#===============================================================
### Specific header for you batch manager:
# voima.fmi.fi definitions
#PBS -N build_clim 
#PBS -q workq
#PBS -l mppwidth=20
#PBS -l mppnppn=20
#PBS -l mppdepth=1
#PBS -l walltime=05:00:00
#PBS -j oe
# You may need to comment this on other clusters than voima.
set -x
cd $PBS_O_WORKDIR
# -C ORCA1_L46_v36_voima -R ECN-D501 -i 1983 -e 2012
export CONFIG=eORCA1_L75_v36_voima
#export CONFIG=ORCA1_L75_v36_voima
#export CONFIG=ORCA1_L46_v36_voima
#export CONFIG=ORCA1_L75_v36_LIM2_voima
export RUN=eO1L7501
export Y1=1961
export Y2=1971

#SBATCH -A snic2014-10-3
#SBATCH --reservation=dcs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J CLIM
#SBATCH -t 05:59:00
#SBATCH -o out_clim_%J.out
#SBATCH -e err_clim_%J.err
###

export BARAKUDA_ROOT=`pwd`

iremap=0 ; REGG="360x180"; # remap to regular lat-lon grid?
iuv=0      ; # Do a climatology for current...
ivt=1      ; # Do a climatology for VT
iamoc=1    ; # Do a climatology for 2D lat-depth AMOC?
ibpsi=0    ; # Do a climatology for barotropic stream function


# Supported ORCA grids:
ORCA_LIST="eORCA1.L75 ORCA1.L75 ORCA1.L46 ORCA1.L42 ORCA2 ORCA2_L46"

# Checking available configs
list_conf=`\ls configs/config_*.sh` ; list_conf=`echo ${list_conf} | sed -e s/'configs\/config_'/''/g -e s/'.sh'/''/g`

# Important bash functions:
. ${BARAKUDA_ROOT}/configs/bash_functions.bash

usage()
{
    echo
    echo "USAGE: ${0} -C <config> -R <run>  (options)"
    echo
    echo "     Available configs are:"
    echo "             => ${list_conf}"
    echo
    echo "   OPTIONS:"
    echo "      -f <years> => how many years per NEMO file? default = 1"
    echo "      -i <YYYY>  => first year, default=1984"
    echo "      -e <YYYY>  => last  year, default=2006"
    echo "      -k   4     => input NEMO files in netcdf 4 (.nc4)"
    echo "      -h         => print this message"
    echo
    exit
}

# Defaults:
NCE="nc"
IFREQ_SAV_YEARS=1

while getopts R:f:i:e:k:C:h option ; do
    case $option in
        R) RUN=${OPTARG};;
        f) IFREQ_SAV_YEARS=${OPTARG} ;;
        i) Y1=${OPTARG} ;;
        e) Y2=${OPTARG} ;;
        k) NCE="nc${OPTARG}" ;;
        C) CONFIG=${OPTARG};;
        h)  usage;;
        \?) usage ;;
    esac
done

if [ "${CONFIG}" = "" -o "${RUN}" = "" ]; then usage ; exit ; fi

for og in ${ORCA_LIST}; do
    ca=""; ca=`echo ${CONFIG} | grep ^${og}` ; if [ "${ca}" != "" ]; then ORCA=${og}; fi
done

if [ "${ORCA}" = "" ]; then echo "ORCA grid of config ${CONFIG} not supported yet"; exit; fi

echo

# sourcing configuration file
fconfig=${BARAKUDA_ROOT}/configs/config_${CONFIG}.sh
echo "Sourcing ${fconfig} !"
if [ -f ${fconfig} ]; then
    echo "${fconfig} found!" ; . ${fconfig}
else
    echo "PROBLEM: cannot find file ${fconfig} !"; exit
fi
echo


if [ ! "${ORCA}" = "${CONF}" ]; then echo "ERROR: ORCA and CONF disagree! => ${ORCA} ${CONF}"; exit; fi
export ORCA=${CONF}



export CONFRUN=${ORCA}-${RUN}


# file prefix (before "_<calendar_info>_grid_X.nc"):
export CPREF=`echo ${NEMO_FILE_PREFIX} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g" -e "s|<TSTAMP>|${TSTAMP}|g"`


CPR=`echo ${CPREF} | sed -e s/"_${TSTAMP}_"/""/g` ; # very begining of file without time stamp and underscore


echo
echo " SETTINGS: "
echo "   *** Y1    = ${Y1} "
echo "   *** Y2    = ${Y2} "
echo "   *** CONFIG  = ${CONFIG} "
echo "   *** GRID  = ${ORCA} "
echo "   *** RUN  = ${RUN} "
echo "   *** CPREF  = ${CPREF} "
echo "   *** IFREQ_SAV_YEARS  = ${IFREQ_SAV_YEARS} "
echo



# Need to be in accordance with the netcdf installation upon whicg cdftools_light is compiled:
export NCDF_DIR=`cat cdftools_light/make.macro | grep ^NCDF_DIR | cut -d = -f2 | sed -e s/' '//g`
echo ; echo "NCDF_DIR = ${NCDF_DIR}"; echo
export LD_LIBRARY_PATH=${NCDF_DIR}/lib:${LD_LIBRARY_PATH}




Y1=`expr ${Y1} + 0`
Y2=`expr ${Y2} + 0`
CY1=`printf "%04d" ${Y1}`
CY2=`printf "%04d" ${Y2}`


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
1            echo "     => but you have not specified ${sf} in NEMO_SAVED_FILES !"; exit
        fi
    else
        gimp_new="${sf} ${gimp_new}"
    fi
    js=`expr ${js} + 1`
done
GRID_IMP=${gimp_new}
echo; echo "File types to import: ${GRID_IMP}"; echo; echo








# We need a scratch/temporary directory to copy these files to and gunzip them:
if [ ! "${SLURM_JOBID}" = "" ]; then
    SCRATCH=`echo ${SCRATCH} | sed -e "s|<JOB_ID>|${SLURM_JOBID}|g"`
    export TMP_DIR=${SCRATCH}
else
        # Likely to be running interactively
    #export SCRATCH=/proj/bolinc/users/x_laubr/tmp
    export TMP_DIR=${SCRATCH}/${RUN}_tmp
fi





VCM=( "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" )



export NEMO_OUT_D=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g"`
if [ ! -d ${NEMO_OUT_D} ]; then echo "Unfortunately we could not find ${NEMO_OUT_D}"; exit; fi



export DIAG_D=${DIAG_DIR}/${CONFRUN} ; mkdir -p ${DIAG_D}/clim

echo; echo "Finding NEMO files into:"; echo "${NEMO_OUT_D}"; echo
echo "Will store files into:"; echo "${DIAG_D}"; echo



mkdir -p ${TMP_DIR} ; cd ${TMP_DIR}/

echo; echo "In:"; pwd




#If doing AMOC, need mesh_mask and basin mask
if [ ${iamoc} -eq 1 -o ${ibpsi} -eq 1 -o ${ivt} -eq 1 ]; then
    if [ "`cext ${MM_FILE}`" = ".gz" ]; then
        cp ${MM_FILE} mesh_mask.nc.gz    ; gunzip -f mesh_mask.nc.gz
    else
        cp ${MM_FILE} mesh_mask.nc
    fi

    #Fix, in case old nemo (prior version 3.6) must rename some metrics param:
    ca=""; ca=`${NCDF_BIN}/ncdump -h mesh_mask.nc  | grep 'e3t('`
    if [ ! "${ca}" = "" ]; then
        echo "Renaming some metrics into mesh_mask.nc !!!"
        ncrename -v e3t_0,e3t_1d -v e3w_0,e3w_1d -v gdept_0,gdept_1d -v gdepw_0,gdepw_1d  mesh_mask.nc
        ncrename -v e3t,e3t_0 -v e3u,e3u_0 -v e3v,e3v_0 -v e3w,e3w_0                      mesh_mask.nc
        echo
    fi

    if [ "`cext ${BM_FILE}`" = ".gz" ]; then
        cp ${BM_FILE} new_maskglo.nc.gz  ; gunzip -f new_maskglo.nc.gz
    else
        cp ${BM_FILE} new_maskglo.nc
    fi
    ln -sf mesh_mask.nc mesh_hgr.nc ; ln -sf mesh_mask.nc mesh_zgr.nc ; ln -sf mesh_mask.nc mask.nc
fi

ls -l ; echo; echo




nby=`expr ${Y2} - ${Y1} + 1`

if [ ! `expr ${nby} % ${IFREQ_SAV_YEARS}` -eq 0 ]; then
    echo " Number of years should be a multiple of ${IFREQ_SAV_YEARS}!"; exit
fi



if [ ${ece_run} -eq 1 ]; then
    if [ ! -d ${NEMO_OUT_D}/001 ]; then echo "ERROR: since ece_run=${ece_run}, there should be a directory 001 in:"; echo " ${NEMO_OUT_D}"; fi
fi


export jyear=${Y1}

while [ ${jyear} -le ${Y2} ]; do


    cyear=`printf "%04d" ${jyear}`
    echo; echo "Year = ${cyear}:"


    cpf=""
    if [ ${ece_run} -eq 1 ]; then
        # Need initial year in 001:
        fin=`\ls ${NEMO_OUT_D}/001/${CPREF}*_grid_T.${NCE}` ; fin=`basename ${fin}`
        YEAR_INI=`echo ${fin} | sed -e "s|${CPREF}||g" | cut -c1-4`
        echo " *** YEAR_INI = ${YEAR_INI}"
        iy=`expr ${jyear} - ${YEAR_INI} + 1` ; dir_ece=`printf "%03d" ${iy}`
        echo " *** ${cyear} => dir_ece = ${dir_ece}"
        cpf="${dir_ece}/"
        echo
    fi
    
    #exit ; #lolo

    TTAG_ann=${cyear}0101_${cyear}1231

    i_get_file=0


    if [ ${IFREQ_SAV_YEARS} -gt 1 ]; then

        FPREF="BIG_"
        if [ $((${jyear}%${IFREQ_SAV_YEARS})) -eq 0 ]; then
            jy1=${jyear} ; jy2=$((${jyear}+${IFREQ_SAV_YEARS}-1))
            cy1=`printf "%04d" ${jy1}` ; cy2=`printf "%04d" ${jy2}`
            i_get_file=1
            TTAG=${cy1}0101_${cy2}1231
        fi

    else
        
        jy1=${jyear} ; jy2=${jyear}
        cy1=`printf "%04d" ${jy1}` ; cy2=`printf "%04d" ${jy2}`
        FPREF=""
        i_get_file=1
        TTAG=${cy1}0101_${cy2}1231
        
    fi

    echo " Year ${cyear} => using files ${cy1} to ${cy2} => ${TTAG}" ; echo
    

    CRTM=${CPREF}${TTAG}
    CRT1=${CPREF}${TTAG_ann}


    if [ ${i_get_file} -eq 1 ]; then

        # Importing required files to tmp dir and gunzipping:
        for gt in ${GRID_IMP}; do
            f2i=${CRTM}_${gt}.${NCE} ; sgz=""
            if [ -f ${NEMO_OUT_D}/${cpf}${f2i}.gz ]; then sgz=".gz"; fi

            check_if_file ${NEMO_OUT_D}/${cpf}${f2i}${sgz}
            if [ ! -f ./${FPREF}${f2i} ]; then
                echo "Importing ${f2i}${sgz} ..."
                cp -L ${NEMO_OUT_D}/${cpf}${f2i}${sgz} ./${FPREF}${f2i}${sgz}
                if [ "${sgz}" = ".gz" ]; then gunzip -f ./${FPREF}${f2i}.gz ; fi
                #
                if [ ! "${NCE}" = "nc4" ]; then
                    # Converting to netcd4 compressed:
                    echo "  + converting ${FPREF}${f2i} to compressed netcdf4..."
                    mv -f ${FPREF}${f2i} ${FPREF}${f2i}3
                    #${NCDF_DIR}/bin/nccopy -k 4 -d 9 ${FPREF}${f2i}3 ${FPREF}${f2i} ; rm -f ${FPREF}${f2i}3
                    ${NCDF_BIN}/nccopy -k 4 -d 9 ${FPREF}${f2i}3 ${FPREF}${f2i} ; rm -f ${FPREF}${f2i}3
                fi
                check_if_file ${FPREF}${f2i}
                echo " ... done!"
            fi
        done


        # Need to create annual files if needed
        if [ ${IFREQ_SAV_YEARS} -gt 1 ]; then
            for gt in grid_T grid_U grid_V ${FILE_ICE_SUFFIX} SBC; do
                ftd=./${FPREF}${CRTM}_${gt}.${NCE} ; # file to divide!
                if [ -f ${ftd} ]; then
                    jy=0
                    while [ ${jy} -lt ${IFREQ_SAV_YEARS} ]; do
                        im1=$((${jy}*12+1)) ;  im2=$((${im1}+11))
                        jytc=$((${jy}+${jy1})) ; cjytc=`printf "%04d" ${jytc}`
                        ftc="./${CPREF}${cjytc}0101_${cjytc}1231_${gt}.${NCE}" ; # file to create!
                        if [ ! -f ${ftc} ]; then
                            echo "Extracting file ${ftc} from ${ftd}, month records: ${im1}=>${im2}"
                            ncks -a -O -F -d time_counter,${im1},${im2} ${ftd} -o ${ftc}
                            echo
                        fi
                        jy=$((${jy}+1))
                    done
                    #debug:
                    rm -f ${ftd}
                fi
            done
        fi

    fi  # if [ ${i_get_file} -eq 1 ]



    # Testing if required files are present now:
    for gt in ${GRID_IMP}; do
        ftt="./${CRT1}_${gt}.${NCE}"
        check_if_file ${ftt}
    done


    echo; echo "All required files are in `pwd` for year ${cyear} !"; echo


    fu=${CPREF}${TTAG_ann}_grid_U.${NCE}
    fv=${CPREF}${TTAG_ann}_grid_V.${NCE}


    echo

    if [ ${ivt} -eq 1 ]; then
        # Creating VT files:
        echo "Doing: ${BARAKUDA_ROOT}/cdftools_light/bin/cdfvT.x ${CPREF}${TTAG_ann} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV}"
        ${BARAKUDA_ROOT}/cdftools_light/bin/cdfvT.x ${CPREF}${TTAG_ann} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV}
        echo "After:"; ls *VT.nc*
        echo
    fi

    echo

    if [ ${iamoc} -eq 1 ]; then
        rm -f moc.nc
        echo " *** doing: ./cdfmoc.x ${fv} ${NN_V} ${NN_V_EIV}"
        ${BARAKUDA_ROOT}/cdftools_light/bin/cdfmoc.x ${fv} ${NN_V} ${NN_V_EIV}
        echo "Done!"; echo; echo; echo
        ncwa -O -a x moc.nc ${CPREF}${TTAG_ann}_MOC.nc ; # removing degenerate x record...
        rm -f moc.nc
        echo "After:"; ls *MOC.nc*
        echo
    fi

    if [ ${ibpsi} -eq 1 ]; then
        rm -f psi.nc
        echo " *** doing: ./cdfpsi.x ${fu} ${fv} V"
        ${BARAKUDA_ROOT}/cdftools_light/bin/cdfpsi.x ${fu} ${fv} V
        echo "Done!"; echo; echo; echo
        mv -f psi.nc ${CPREF}${TTAG_ann}_PSI.nc
        echo "After:"; ls *PSI.nc*
        echo
    fi

    export jyear=`expr ${jyear} + 1`

done



if [ ${ivt} -eq 1 ]; then
    fo=aclim_${CONFRUN}_${CY1}-${CY2}_VT.nc
    # Averaged VT file:
    ncra -O ${CPREF}*_VT.nc -o ${fo}
    rm -f ${CPREF}*_VT.nc

    # Converting to netcdf4 with maximum compression level:
    ${NCDF_BIN}/nccopy -k 4 -d 9 ${fo} ${fo}4 ;  rm -f ${fo}
    mv -f ${fo}4 ${DIAG_D}/clim/
fi

if [ ${iamoc} -eq 1 ]; then
    fo=aclim_${CONFRUN}_${CY1}-${CY2}_MOC.nc
    # Averaged MOC file:
    ncra -O ${CPREF}*_MOC.nc -o ${fo}
    # Converting to netcdf4 with maximum compression level:
    ${NCDF_BIN}/nccopy -k 4 -d 9 ${fo} ${fo}4 ;  rm -f ${fo}
    mv -f ${fo}4 ${DIAG_D}/clim/
fi


if [ ${ibpsi} -eq 1 ]; then
    fo=aclim_${CONFRUN}_${CY1}-${CY2}_PSI.nc
    # Averaged PSI file:
    ncra -O ${CPREF}*_PSI.nc -o ${fo}
    # Converting to netcdf4 with maximum compression level:
    ${NCDF_BIN}/nccopy -k 4 -d 9 ${fo} ${fo}4 ;  rm -f ${fo}
    mv -f ${fo}4 ${DIAG_D}/clim/
fi

echo;echo;echo;






echo "Phase 2:"; ls ; echo



# Mean monthly climatology

for suff in grid_T grid_U grid_V icemod SBC VT MOC PSI; do


    if [ -f ./${CPREF}${CY1}0101_${CY1}1231_${suff}.${NCE} ]; then
        
        echo ; echo ; echo ; echo " Treating ${suff} files!"; echo

        f2c=mclim_${CONFRUN}_${CY1}-${CY2}_${suff}.nc
        f2c_reg=mclim_${CONFRUN}_${CY1}-${CY2}_${REGG}.nc

        rm -f ${DIAG_D}/clim/${f2c}* ${DIAG_D}/clim/${f2c_reg}*

        echo

        jm=0
        for cm in ${VCM[*]}; do

            jm=`expr ${jm} + 1`

            if [ -f ./${CPREF}${CY1}0101_${CY1}1231_${suff}.${NCE} ]; then
                echo; ls ; echo
                echo "ncra -F -O -d time_counter,${jm},,12 ${CPREF}*0101_*1231_${suff}.${NCE} -o mean_m${cm}_${suff}.nc"
                ncra -F -O -d time_counter,${jm},,12 ${CPREF}*0101_*1231_${suff}.${NCE} -o mean_m${cm}_${suff}.nc
                echo
            fi

        done

        rm -f ${CPREF}*0101_*1231_${suff}.${NCE}

        echo; ls ; echo
        echo "ncrcat -O  mean_m*_${suff}.nc out_${suff}.nc"
        ncrcat -O  mean_m*_${suff}.nc out_${suff}.nc
        rm mean_m*_${suff}.nc
        echo

        mv -f out_${suff}.nc ${f2c}

        if [ ${iremap} -eq 1 ]; then
            echo; echo "cdo remapbil,r${REGG} ${f2c} ${f2c_reg}"
            cdo remapbil,r${REGG} ${f2c} ${f2c_reg}
            echo
        fi

        # Converting to netcdf4 with maximum compression level:
        echo; ls ; echo
        echo "${NCDF_BIN}/nccopy -k 4 -d 9 ${f2c} ${f2c}4 ;  rm -f ${f2c}"
        ${NCDF_BIN}/nccopy -k 4 -d 9 ${f2c} ${f2c}4 ;  rm -f ${f2c}
        echo
        mv -f ${f2c}4    ${DIAG_D}/clim/

        if [ ${iremap} -eq 1 ]; then mv -f ${f2c_reg} ${DIAG_D}/clim/ ; fi


    else
        echo ; echo ; echo ; echo " Ignoring monthly ${suff} files!"; echo
    fi

done ; # loop along files suffixes


echo;echo
echo "${CY1}-${CY2}" > ${DIAG_D}/clim/last_clim
echo "Climatology saved into: ${DIAG_D}/clim/"
echo;echo


#debug:
rm -rf ${TMP_DIR} 2>/dev/null

