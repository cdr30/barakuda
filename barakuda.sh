#!/usr/bin/env bash

#==============================================================
#
#                    B A R A K U D A
#
#    An OCEAN MONITORING python environment for NEMO
#
#             L. Brodeau, 2009-2015
#
#===============================================================

export BARAKUDA_ROOT=`pwd`

# Supported ORCA grids:
ORCA_LIST="ORCA1.L75 ORCA1.L46 ORCA1.L42 ORCA2 ORCA2_L46"

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
    echo "      -y <YYYY> => force initial year to YYYY"
    echo
    echo "      -z <YYYY> => force BaraKuda to stop at year YYYY"
    echo
    echo "      -f        => forces creation of diagnostic page eventhough treatment of output files is not finished"
    echo
    echo "      -e        => create the HTML diagnostics page on local or remote server"
    echo
    echo "      -E        => same as '-e' but also create the 2D plots / observations"
    echo "                   you need to have built a climatology of your run with 'build_clim.sh' first!"
    echo
    echo "      -c <run>  => when '-E' specified, 2D comparison diagnostics are performed "
    echo "                   against the climatology of another run <run> rather than observations"
    echo
    echo "      -h        => print this message"
    echo
    exit
}


# Some defaults:
YEAR0=""
LFORCE_YEAR1=false
LFORCE_END=false
export RUNREF=""
ISTAGE=1 ; # 1 => generation of data diagnostic files
#          # 2 => creation of figures and diagnostic HTML page
LFORCEDIAG=false
l_clim_diag=false

while getopts C:R:y:z:c:feEh option ; do
    case $option in
        C) CONFIG=${OPTARG};;
        R) RUN=${OPTARG};;
        y) YEAR0=${OPTARG} ; LFORCE_YEAR1=true ;;
        z) YEARN=${OPTARG} ; LFORCE_END=true ;;
        c) export RUNREF=${OPTARG} ;;
        f) LFORCEDIAG=true;;
        e) ISTAGE=2;;
        E) ISTAGE=2 ; l_clim_diag=true ;;
        h)  usage;;
        \?) usage ;;
    esac
done


if [ "${CONFIG}" = "" -o "${RUN}" = "" ]; then usage ; exit ; fi

if [ "${RUNREF}" != "" -a ${ISTAGE} -eq 1 ]; then    
    echo; echo " WARNING: option '-c' only makes sense when '-e' or '-E' are specified !"
    sleep 2; echo
fi

for og in ${ORCA_LIST}; do
    echo " ${og} / ${CONFIG}"
    ca=""; ca=`echo ${CONFIG} | grep ${og}` ; if [ "${ca}" != "" ]; then ORCA=${og}; fi
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

echo
if [ "${CANOPY_PATH}" = "" ]; then echo "ERROR: CANOPY_PATH is not set! => add it to config file"; exit; fi
PYTH="${CANOPY_PATH}/bin/python -W ignore" ; # which Python installation to use
export PYTHONPATH=${CANOPY_PATH}/lib/python2.7:${CANOPY_PATH}/lib/python2.7/site-packages:${BARAKUDA_ROOT}/python/modules ; # PATH to python barakuda modules
PYBRKD_EXEC_PATH=${BARAKUDA_ROOT}/python/exec         ; # PATH to python barakuda executable

echo " CANOPY_PATH => "${CANOPY_PATH} ; echo

echo "  TRANSPORT_SECTION_FILE => ${TRANSPORT_SECTION_FILE} !" ; echo


if ${l_clim_diag} ; then
    echo
    echo " Files containing climatologies to be used:"
    echo " T 3D => ${F_T_CLIM_3D_12} => ${NN_T_CLIM}"
    echo " S 3D => ${F_S_CLIM_3D_12} => ${NN_S_CLIM}"
    echo " SST  => ${SST_CLIM_12} => ${NN_SST_CLIM}"
    echo
    for ff in ${F_T_CLIM_3D_12} ${F_S_CLIM_3D_12} ${SST_CLIM_12}; do
        if [ ! -f ${F_T_CLIM_3D_12} ]; then echo "ERRO: ${ff} is missing!"; exit; fi
    done
fi


# Testing if NCO is installed:
which ncks 1>out 2>/dev/null; ca=`cat out`; rm -f out
if [ "${ca}" = "" ]; then echo "Install NCO!!!"; echo; exit; fi


# Names for temperature, salinity, u- and v-current...
if [ "${NN_T}" = "" -o "${NN_S}" = "" -o "${NN_U}" = "" -o "${NN_V}" = "" ]; then
    echo "NN_T, NN_S, NN_U and NN_V are NOT given a value into"
    echo " in ${fconfig} "
    echo "  => using default names: thetao, so, uo, vo" ; echo
    NN_T="thetao"; NN_S="so"; NN_U="uo"; NN_V="vo"
fi

echo ; echo " *** NN_T=${NN_T}, NN_S=${NN_S}, NN_U=${NN_U} and NN_V=${NN_V} "; echo



# What grid-type files to work with:
GRID_IMP=""

if [ ${i_do_mean} -gt 0 ]; then GRID_IMP="${GRID_IMP} grid_T"; fi
if [ ${i_do_trsp} -gt 0 ]; then GRID_IMP="${GRID_IMP} grid_T grid_U grid_V"; fi
if [ ${i_do_mht}  -gt 0 ]; then GRID_IMP="${GRID_IMP} grid_T grid_U grid_V"; fi
if [ ${i_do_sigt} -gt 0 ]; then GRID_IMP="${GRID_IMP} grid_T grid_U grid_V"; fi
if [ ${i_do_amoc} -gt 0 ]; then GRID_IMP="${GRID_IMP} grid_V" ; fi
if [ ${i_do_ice}  -gt 0 ]; then GRID_IMP="${GRID_IMP} ${FILE_ICE_SUFFIX}"; fi
if [ ${i_do_ssx_box} -gt 0 ]; then GRID_IMP="${GRID_IMP} grid_T"; fi
if [ ${i_do_bb} -gt 0 ]; then GRID_IMP="${GRID_IMP} grid_T grid_U grid_V"; fi
if [ ${i_do_box_TS_z} -gt 0 ]; then GRID_IMP="${GRID_IMP} grid_T"; fi
if [ ${i_do_dmv} -gt 0 ];      then GRID_IMP="${GRID_IMP} grid_T"; fi
if [ ${i_do_zcrit} -gt 0 ];    then GRID_IMP="${GRID_IMP} grid_T"; fi

# Not fully supported yet:
ca=" => diagnostic totally beta and not fully supported yet!"
if [ ${i_do_sect} -gt 0 ]; then echo " *** i_do_sect ${ca}"; exit; fi
if [ ${i_do_amo}  -gt 0 ]; then echo " *** i_do_amo  ${ca}"; exit; fi
if [ ${i_do_icet} -gt 0 ]; then echo " *** i_do_icet ${ca}"; exit; fi
if [ ${i_do_flx}  -gt 0 ]; then echo " *** i_do_flx  ${ca}"; exit; fi


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
    js=`expr ${js} + 1`
done
GRID_IMP=${gimp_new}

echo; echo "File types to import: ${GRID_IMP}"; echo; echo

boo="EOF"


if [ ${ISTAGE} -eq 1 ]; then
    # List of CDFTOOLS executables needed for the diagnostics:
    L_EXEC="cdfmaxmoc.x cdfmoc.x cdfvT.x cdftransportiz.x cdficediags.x cdfmhst.x cdfsigtrp.x"
    for ex in ${L_EXEC}; do check_if_file cdftools_light/bin/${ex} "Compile CDFTOOLS executables!"; done
fi



# Need to be in accordance with the netcdf installation upon whicg cdftools_light is compiled:
export NCDF_DIR=`cat cdftools_light/make.macro | grep ^NCDF_DIR | cut -d = -f2 | sed -e s/' '//g`
echo ; echo "NCDF_DIR = ${NCDF_DIR}"; echo
export LD_LIBRARY_PATH=${NCDF_DIR}/lib:${LD_LIBRARY_PATH}

# Exporting some variables needed by the python scripts:
export RUN=${RUN}
export CONFRUN=${ORCA}-${RUN}
export DIAG_D=${DIAG_DIR}/${CONFRUN}


if [ ${ISTAGE} -eq 1 ]; then
    # We need a scratch/temporary directory to copy these files to and gunzip them:
    if [ "${SLURM_JOBID}" = "" -a `hostname` = "triolith1" ]; then
        # Likely to be running interactively on triolith         
        export SCRATCH=${HOME}/tmp
        export TMP_DIR=${SCRATCH}/${RUN}_tmp
    else
        # Normal case:
        SCRATCH=`echo ${SCRATCH} | sed -e "s|<JOB_ID>|${SLURM_JOBID}|g"` 
        export TMP_DIR=${SCRATCH}/${RUN}
    fi
    echo " IMPORTANT the SCRATCH work directory is set to:" ; echo " ${SCRATCH}"
else
    export TMP_DIR=${DIAG_D}/tmp
fi

echo " IMPORTANT the TMP_DIR work directory is set to:" ; echo " ${TMP_DIR}"; echo ; sleep 2

rm -rf ${TMP_DIR}
mkdir -p ${DIAG_D} ${TMP_DIR}

export NEMO_OUT_D=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g"`
if [ ! -d ${NEMO_OUT_D} ]; then echo "Unfortunately we could not find ${NEMO_OUT_D}"; exit; fi


echo; echo " * Config to be used: ${CONFIG} => ORCA grid is ${ORCA}"
echo " * Run is ${RUN} "; echo " * Files are stored into ${NEMO_OUT_D}"; echo; sleep 2







#############################
# Okay really starting now! #
#############################

# file prefix (before "_<calendar_info>_grid_X.nc"):
export CPREF=`echo ${NEMO_FILE_PREFIX} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g" -e "s|<TSTAMP>|${TSTAMP}|g"`

# Will have a look at NEMO files in the directory where they are stored
cd ${NEMO_OUT_D}/

if [ ${ISTAGE} -eq 1 ]; then
    
    if [ ${ece_run} -eq 1 ]; then
        if [ ! -d 001 ]; then echo "ERROR: since ece_run=${ece_run}, there should be a directory 001 in:"; echo " ${NEMO_OUT_D}"; exit ; fi
        nby_ece=`ls -d */ | wc -l` ; echo " ${nby_ece} years have been completed..."
        cd 001/
    fi

    if ! ${LFORCE_YEAR1}; then
    # Try to guess the first year from stored "grid_T" files:
        export YEAR_INI=`\ls ${CPREF}*${ctest}* | sed -e s/"${CPREF}"/""/g | head -1 | cut -c1-4`
        echo ${YEAR_INI} |  grep "[^0-9]" >/dev/null ;   # Checking if it's an integer:
        if [ ! "$?" -eq 1 ]; then
            echo "ERROR: it was imposible to guess initial year from your input files"
            echo "       maybe the directory contains non-related files..."
            echo "      => use the -y <YEAR> switch to force the initial year!"; exit
        fi
        echo " Initial year guessed from stored files => ${YEAR_INI}"; echo
        YEAR_INI=`expr ${YEAR_INI} + 0`  ; # example: 1 instead of 0001...
        IFREQ_SAV_YEARS=1
        Y2=`\ls ${CPREF}*${ctest}* | sed -e s/"${CPREF}"/""/g | head -1 | cut -c10-13`
        YIr=`expr ${YEAR_INI} + 0`; Y2r=`expr ${Y2} + 0`
        if [ ${Y2r} -gt ${YIr} ]; then IFREQ_SAV_YEARS=$((${Y2r}-${YIr}+1)); fi
        echo "1 NEMO file contains ${IFREQ_SAV_YEARS} year(s)"; echo
    else
    #echo " FIX me !!! With IFREQ_SAV_YEARS !!!"; exit
        export YEAR_INI=${YEAR0}
        IFREQ_SAV_YEARS=10
        echo " Initial year forced to ${YEAR_INI}"; echo
    fi


    cd ${NEMO_OUT_D}/

    if [ ${ece_run} -eq 1 ]; then
        dir_end=`printf "%03d" ${nby_ece}`
        if [ ! -d ${dir_end} ]; then echo "ERROR: since ece_run=${ece_run}, there should be a directory ${dir_end} in:"; echo " ${NEMO_OUT_D}"; exit ; fi
        YEAR_END=`expr ${YEAR_INI} + ${nby_ece}`
    else
        export YEAR_END=`\ls ${CPREF}*${ctest}* | sed -e s/"${CPREF}"/''/g | tail -1 | cut -c1-4`
        echo ${YEAR_END} |  grep "[^0-9]" >/dev/null; # Checking if it's an integer
        if [ ! "$?" -eq 1 ]; then
            echo "ERROR: it was imposible to guess the year coresponding to the last saved year!"
            echo "       => check your NEMO output directory and file naming..."; exit
        fi
        YEAR_END=`expr ${YEAR_END} + ${IFREQ_SAV_YEARS} - 1`
    fi
    echo " Last year guessed from stored files => ${YEAR_END}"; echo


    
    echo ${IFREQ_SAV_YEARS} > ${DIAG_D}/numb_year_per_file.info
    echo ${YEAR_INI}        > ${DIAG_D}/first_year.info
    
else
    
    for fc in "first_year" "numb_year_per_file" "last_year_done"; do
        ff=${DIAG_D}/${fc}.info
        if [ ! -f ${ff} ]; then echo "ERROR: file ${ff} is missing!"; exit; fi
    done

    # Only doing plots...
    export YEAR_INI=`cat ${DIAG_D}/first_year.info`
    export YEAR_END=`cat ${DIAG_D}/last_year_done.info`
    export IFREQ_SAV_YEARS=`cat ${DIAG_D}/numb_year_per_file.info`
    
fi


cyear_ini=`printf "%04d" ${YEAR_INI}` ; cyear_end=`printf "%04d" ${YEAR_END}`

if ${LFORCE_END}; then
    if [ ${YEARN} -le ${YEAR_INI} ]; then echo "ERROR: forced stop year is before first year!"; exit; fi
fi



jyear=${YEAR_INI}

fcompletion=${DIAG_D}/last_year_done.info
if [ -f ${fcompletion} ]; then jyear=`cat ${fcompletion}`; jyear=`expr ${jyear} + 1`; fi

cd ${TMP_DIR}/


# Importing mesh_mask files:
check_if_file ${MM_FILE}
if [ ! -f ./mesh_mask.nc ]; then
    if [ "`cext ${MM_FILE}`" = ".gz" ]; then
        echo "Gunzipping mesh_mask.nc.gz..."
        cp -L ${MM_FILE} mesh_mask.nc.gz ; gunzip -f mesh_mask.nc.gz
    else
        cp -L ${MM_FILE} mesh_mask.nc
    fi
fi

#Fix, in case old nemo (prior version 3.6) must rename some metrics param:
ca=""; ca=`${NCDF_DIR}/bin/ncdump -h mesh_mask.nc  | grep 'e3t('`
if [ ! "${ca}" = "" ]; then
    echo "Renaming some metrics into mesh_mask.nc !!!"
    ncrename -v e3t_0,e3t_1d -v e3w_0,e3w_1d -v gdept_0,gdept_1d -v gdepw_0,gdepw_1d  mesh_mask.nc
    ncrename -v e3t,e3t_0 -v e3u,e3u_0 -v e3v,e3v_0 -v e3w,e3w_0                      mesh_mask.nc
    echo
fi

check_if_file ${BM_FILE}
if [ "`cext ${BM_FILE}`" = ".gz" ]; then
    cp -L ${BM_FILE} new_maskglo.nc.gz ; gunzip -f new_maskglo.nc.gz
else
    cp -L ${BM_FILE} new_maskglo.nc
fi

# Setting MM BM files to the imported and unzipped version:
export MM_FILE=${TMP_DIR}/mesh_mask.nc
export BM_FILE=${TMP_DIR}/new_maskglo.nc




if [ ${ISTAGE} -eq 1 ]; then
    # Importing cdftools executables:
    for ex in ${L_EXEC}; do cp -L ${BARAKUDA_ROOT}/cdftools_light/bin/${ex} . ; done
fi

sgz=""





# L O O P   A L O N G   Y E A R S
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lcontinue=true

if ${LFORCEDIAG}; then lcontinue=false; fi

while ${lcontinue}; do

    cyear=`printf "%04d" ${jyear}`

    #lolo
    cpf=""
    if [ ${ece_run} -eq 1 ]; then
        iy=`expr ${jyear} - ${YEAR_INI} + 1` ; dir_ece=`printf "%03d" ${iy}`
        echo " *** ${cyear} => dir_ece = ${dir_ece}"
        cpf="${dir_ece}/"
    fi

    TTAG_ann=${cyear}0101_${cyear}1231

    i_get_file=0


    if [ $((${jyear}%${IFREQ_SAV_YEARS})) -eq 0 ]; then

        jy1=${jyear} ; jy2=$((${jyear}+${IFREQ_SAV_YEARS}-1))

        cy1=`printf "%04d" ${jy1}` ; cy2=`printf "%04d" ${jy2}`
        cy1m=`printf "%04d" $((${jy1}-${IFREQ_SAV_YEARS}))` ; cy2m=`printf "%04d" $((${jy2}-${IFREQ_SAV_YEARS}))`
        i_get_file=1

        echo " ${cyear} => from ${cy1} ${cy2}"

        TTAG=${cy1}0101_${cy2}1231 # calendar-related part of the file name

        # Testing if the current year-group has been done
        for ft in ${GRID_IMP}; do            
            ftst=${NEMO_OUT_D}/${cpf}${CPREF}${TTAG}_${ft} ;  cfxt="0"
            for ca in "nc" "nc.gz" "nc4"; do
                if [ -f ${ftst}.${ca} ]; then cfxt="${ca}"; fi
            done
            if [ ${cfxt} = "0" ]; then
                echo "Year(s) ${cy1}-${cy2} is not completed yet:"; echo " => ${ftst}(?) is missing"; echo
                lcontinue=false
            fi
        done
        
    fi ; # if [ $((${jyear}%${IFREQ_SAV_YEARS})) -eq 0 ]


    CRTM=${CPREF}${TTAG}
    CRT1=${CPREF}${TTAG_ann}


    if ${lcontinue}; then

        if [ ${ISTAGE} -eq 2 ]; then
            echo; echo "You cannot create figures and HTML pages yet!"
            echo " => finish treating the results first by launching barakuda.sh without the '-e' switch."
            exit
        fi

        echo; echo; echo "Run ${RUN}: Generating diagnostic data for ${cyear}..."; echo


        # On what file type to test file presence:
        cgrid_test=`echo ${GRID_IMP} | cut -d ' ' -f2`
        echo " *** testing on files \"${cgrid_test}\" !"; echo




        l_happy=false

        while ! ${l_happy} ; do

            if [ ${IFREQ_SAV_YEARS} -eq 1 ]; then l_happy=true; fi

            rm -f *.tmp

            if [ ${i_get_file} -eq 1 ]; then

                echo " => gonna get ${CRTM}_* files..."

                # Importing required files to tmp dir and unzipping:
                for gt in ${GRID_IMP}; do

                    f2i=${CRTM}_${gt}.nc ;   sgz=""
                    
                    for ca in ".gz" "4"; do
                        if [ -f ${NEMO_OUT_D}/${cpf}${f2i}${ca} ]; then sgz="${ca}"; fi
                    done
                    
                    check_if_file ${NEMO_OUT_D}/${cpf}${f2i}${sgz}
                    
                    if [ ! -f ./${f2i} ]; then
                        echo "Importing ${f2i}${sgz} ..."
                        echo "rsync -L ${NEMO_OUT_D}/${cpf}${f2i}${sgz} `pwd`/"
                        rsync -L ${NEMO_OUT_D}/${cpf}${f2i}${sgz} ./
                        if [ "${sgz}" = ".gz" ]; then gunzip -f ./${f2i}.gz ; fi
                        if [ "${sgz}" = "4"   ]; then
                            if ${L_CONV2NC3}; then
                                echo "Need to convert back to netcdf3!"; echo "nccopy -k classic ${f2i}4 ${f2i}"
                                ${NCDF_DIR}/bin/nccopy -k classic ${f2i}4 ${f2i} ;  rm ${f2i}4
                            else
                                echo "mv ./${f2i}4 ./${f2i}"
                                mv ./${f2i}4 ./${f2i}
                            fi
                        fi
                        check_if_file ${f2i}
                        echo " ... done!"; echo
                    else
                        echo " ${f2i}${sgz} was already in `pwd`"
                    fi
                done


                # Need to create annual files if more than 1 year in 1 once NEMO file
                if [ ${IFREQ_SAV_YEARS} -gt 1 ]; then
                    for gt in ${GRID_IMP}; do
                        ftd=./${CRTM}_${gt}.nc ; # file to divide!
                        if [ -f ${ftd} ]; then
                            jy=0
                            while [ ${jy} -lt ${IFREQ_SAV_YEARS} ]; do
                                im1=$((${jy}*12+1)) ;  im2=$((${im1}+11))
                                jytc=$((${jy}+${jy1})) ; cjytc=`printf "%04d" ${jytc}`
                                ftc="./${CPREF}${cjytc}0101_${cjytc}1231_${gt}.nc" ; # file to create!
                                if [ ! -f ${ftc} ]; then
                                    echo "Extracting file ${ftc} from ${ftd}, month records: ${im1}=>${im2}"
                                    echo "=> ncks -a -O -F -d time_counter,${im1},${im2} ${ftd} -o ${ftc}"
                                    ncks -a -O -F -d time_counter,${im1},${im2} ${ftd} -o ${ftc}
                                    echo
                                fi
                                jy=$((${jy}+1))
                            done
                            rm -f ${ftd}
                        else
                            echo " ${ftd} is not here!!!"
                        fi
                    done
                fi

            fi  # if [ ${i_get_file} -eq 1 ]

            echo

            echo "   *** jyear = ${jyear}, CRT1 = ${CRT1}"



            # In case the job crashed testing only on ${cgrid_test} file:
            if [ ${IFREQ_SAV_YEARS} -gt 1 ]; then
                if [ ! -f ./${CRT1}_${cgrid_test}.nc ]; then
                    echo "${CRT1}_${cgrid_test}.nc is missing !!!"
                    jy1=$((${jyear}-${jyear}%${IFREQ_SAV_YEARS})) ; jy2=$((${jy1}+${IFREQ_SAV_YEARS}-1))
                    cy1=`printf "%04d" ${jy1}` ; cy2=`printf "%04d" ${jy2}`; TTAG=${cy1}0101_${cy2}1231
                    CRTM=${CPREF}${TTAG}
                    echo "Should re-import files ${CRTM}* cause something went wrong...."; echo
                    i_get_file=1
                else
                    l_happy=true
                fi
            fi

        done ; #while ! ${l_happy}; do



        # Testing if ALL required files are present now:
        for gt in ${GRID_IMP}; do
            ftt="./${CRT1}_${gt}.nc" ;  check_if_file ${ftt}
        done


        echo; echo "All required files are in `pwd` for year ${cyear} !"; echo


        # Files to work with for current year:
        ft=${CRT1}_grid_T.nc
        fu=${CRT1}_grid_U.nc
        fv=${CRT1}_grid_V.nc
        fg=${CRT1}_${FILE_ICE_SUFFIX}.nc ; # can be icemod or grid_T ....
        fvt=${CRT1}_VT.nc


        # ~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Creating VT file if needed
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_trsp} -gt 0 -o ${i_do_mht} -eq 1 ]; then
            if [ ! -f ${fvt} ]; then
                echo; echo; echo " *** doing: ./cdfvT.x ${CPREF}${TTAG_ann} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV}"
                ./cdfvT.x ${CPREF}${TTAG_ann} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV}
                echo "Done!"; echo; echo
            fi
        fi
        

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Computing time-series of spatially-averaged variables
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_mean} -eq 1 ]; then
            echo; echo; echo "Global monthly values"
            echo "CALLING: mean.py ${ft} ${jyear}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/mean.py ${ft} ${jyear}
        fi



        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Computing time-series of spatially-averaged variables
        # on boxes (saving the variable on 2D box too...
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_ssx_box} -eq 1 ]; then
            echo; echo; echo "Box monthly values"
            echo "CALLING: ssx_boxes ${ft} ${jyear} ${NN_SST} ${NN_SSS}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/ssx_boxes.py ${ft} ${jyear} ${NN_SST} ${NN_SSS}
        fi
        

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  AMOC  ( Max of Atlantic MOC for several latitude bands )
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_amoc} -eq 1 ]; then
            echo; echo; echo "AMOC"
            rm -f moc.nc *.tmp

            if [ "${LMOCLAT}" = "" ]; then
                echo "AMOC => specify latitude bands with variable LMOCLAT into the config file!!!"; exit
            fi

            echo " *** doing: ./cdfmoc.x ${fv} ${NN_V} ${NN_V_EIV}"
            ./cdfmoc.x ${fv} ${NN_V} ${NN_V_EIV}
            echo "Done!"; echo; echo; echo

            for clat in ${LMOCLAT}; do
                cslat=`echo ${clat} | sed -e s/'-'/' '/g`
                echo "  *** ./cdfmaxmoc.x moc.nc atl ${cslat} 500 1500 ${jyear} ${DIAG_D}"
                ./cdfmaxmoc.x moc.nc atl ${cslat} 500 1500 ${jyear} ${DIAG_D}
                echo "Done!"; echo
            done

        fi




        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # VOLUME, HEAT and SALT transports through specified section
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if [ ${i_do_trsp} -gt 0 ]; then

            echo; echo; echo "Transports of volume, heat and salt through different sections"

            if [ "${TRANSPORT_SECTION_FILE}" = "" ]; then
                echo "Please specify which TRANSPORT_SECTION_FILE to use into the config file!" ; exit
            fi

            if [ ! -f ./transportiz.dat ]; then
                check_if_file ${TRANSPORT_SECTION_FILE}
                cp ${TRANSPORT_SECTION_FILE} ./transportiz.dat
            fi

            if [ ${i_do_trsp} -eq 1 ]; then z1_trsp="" ; z2_trsp=""; fi


            if [ ! -f ${fvt} ]; then
                echo "PROBLEM: file ${fvt} is not here, skipping transport section!"
            else
                echo " *** doing: ./cdftransportiz.x ${CPREF}${TTAG_ann} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} ${jyear} ${DIAG_D} ${z1_trsp} ${z2_trsp}"
                ./cdftransportiz.x ${CPREF}${TTAG_ann} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} ${jyear} ${DIAG_D} ${z1_trsp} ${z2_trsp}
                echo "Done!"; echo; echo
                
                rm -f *.tmp broken_line_*
            fi
        fi   ; # ${i_do_trsp} -gt 0
        
        





        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Budget and other stuffs on a given rectangular box!
        # It provides time-series depending only on time (not depth)
        # budget_rectangle_box.py
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if [ ${i_do_bb} -gt 0 ]; then

            echo; echo; echo "Budget and other stuffs on rectangular boxes!"

            if [ "${FILE_DEF_BOXES}" = "" ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi
            echo " *** doing: ${PYTH} ${PYBRKD_EXEC_PATH}/budget_rectangle_box.py ${cyear} 100 uv"
            ${PYTH} ${PYBRKD_EXEC_PATH}/budget_rectangle_box.py ${cyear} 100 uv
            echo "Done!"; echo; echo

        fi




        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Meridional heat and salt transport
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_mht} -eq 1 ]; then

            echo; echo; echo "Meridional transport of heat and salt"

            fo=${DIAG_D}/merid_transport_T_S_${CONFRUN}.nc

            if [ ! -f ${fvt} ]; then
                echo "PROBLEM: file ${fvt} is not here, skipping meridional transports section"
            else
                rm -f merid_heat_trp.dat merid_salt_trp.dat
                echo " *** doing: ./cdfmhst.x ${fvt} ${fo} ${jyear}"
                ./cdfmhst.x ${fvt} ${fo} ${jyear}
                echo "Done!"; echo; echo; echo
            fi

        fi













        #~~~~~~~~~~~~~~~~~~~~~~~
        #  AMO SST time series
        #~~~~~~~~~~~~~~~~~~~~~~~

        if [ ${i_do_amo} -eq 1 ]; then

            diro=${DIAG_D}/amo ; mkdir -p ${diro}

            echo; echo; echo "AMO SST time series "
            echo; echo
            echo "CALLING: amo.py ${cyear}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/amo.py ${cyear}
            echo;echo;
            mv -f  AMO_SST_Atl_${CONFRUN}_${cyear}.nc ${diro}/
            echo
        fi
        # End SST time series



        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  Transport by sigma-class
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if [ ${i_do_sigt} -eq 1 ]; then

            echo; echo
            if [ ! -f ./dens_section.dat ]; then
                if [ -f ${DENSITY_SECTION_FILE} ]; then
                    echo "Copying ${DENSITY_SECTION_FILE} to here: `pwd` !"; cp ${DENSITY_SECTION_FILE} ./dens_section.dat
                else
                    echo; echo "WARNING: Can't do Transport by sigma-class: ${DENSITY_SECTION_FILE} is missing!!!"
                fi
            fi

            echo " *** doing: ./cdfsigtrp.x ${ft} ${fu} ${fv} 24.8 28.6 19 ${jyear} ${DIAG_D}  ${NN_T} ${NN_S} ${NN_U} ${NN_V}"; echo
            ./cdfsigtrp.x ${ft} ${fu} ${fv} 24.8 28.6 19 ${jyear} ${DIAG_D} ${NN_T} ${NN_S} ${NN_U} ${NN_V}
            echo "Done!"; echo

        fi
        # End Sigma-Class








        # Vertical meridional or zonal sections:
        if [ ${i_do_sect} -eq 1 ]; then
            diro=${DIAG_D}/sections ; mkdir -p ${diro}
            if [ ${VSECT_NM} = "" -o ${VSECT_JI} = "" -o ${VSECT_JJ} = "" ]; then
                echo "VSECT_NM, VSECT_JI and VSECT_JJ must be defined in:"
                echo "${fconfig}" ; exit
            fi
            js=0
            for cs in ${VSECT_NM[*]}; do
                cc=`echo ${ft} | sed -e s/".nc"/"_${cs}.nc"/g`
                fo="section_${cc}"
                ncks -h -O -v nav_lat,nav_lon,thetao,so \
                    -d x,${VSECT_JI[${js}]} \
                    -d y,${VSECT_JJ[${js}]} \
                    ${ft} -o ${fo}
                mv -f ${fo} ${diro}/
                js=`expr ${js} + 1`
            done
            echo
        fi











        # SEA-ICE
        # ~~~~~~~
        if [ ${i_do_ice} -eq 1 ]; then

            echo; echo; echo "Sea-ice extent and volume..."

            rm -f tmp_ice.nc

            echo "ncks  -A -v ${NN_ICEF} ${fg} -o tmp_ice.nc"
            ncks  -A -v ${NN_ICEF} ${fg} -o tmp_ice.nc
            ncrename -v ${NN_ICEF},ice_frac tmp_ice.nc

            coic=""
            if [ "${NN_ICET}" = "" ]; then
                coic="oic" ; # means only ice concentration available!
            else
                echo "ncks  -A -v ${NN_ICET} ${fg} -o tmp_ice.nc"
                ncks  -A -v ${NN_ICET} ${fg} -o tmp_ice.nc
                ncrename -v ${NN_ICET},ice_thic tmp_ice.nc
            fi
            
            echo " *** doing: ./cdficediags.x tmp_ice.nc ${jyear} ${DIAG_D} ${coic}"
            ./cdficediags.x tmp_ice.nc ${jyear} ${DIAG_D} ${coic}
            echo "Done!"; echo; echo; echo
            rm -f tmp_ice.nc

        fi







        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  Vertical profiles of T,S and density on a given box
        #     => it provides time-series depending on time and depth
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if [ ${i_do_box_TS_z} -gt 0 ]; then

            if [ "${FILE_DEF_BOXES}" = "" ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi

            echo "CALLING: prof_TS_z_box.py ${cyear}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/prof_TS_z_box.py ${cyear}
            echo;echo

        fi






        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  Deep Mixed Volume (DMV) on a given box
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if [ ! -z ${i_do_dmv} -a ${i_do_dmv} -gt 0 ]; then

            if [ "${FILE_DEF_BOXES}" = "" ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi

            echo "CALLING: dmv.py ${cyear}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/dmv.py ${cyear}
            echo;echo

        fi


        if [ ${i_do_zcrit} -gt 0 ]; then

            if [ "${FILE_DEF_BOXES}" = "" ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi

            echo "CALLING: zcrit_conv.py ${cyear}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/zcrit_conv.py ${cyear}
            echo;echo

        fi











        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Solid freshwater transport associated with sea-ice drift
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_icet} -eq 1 ]; then

            echo; echo; echo "Transports of solid freshwater (sea-ice) through different sections"

            if [ ! "${FILE_ICE_SUFFIX}" = "icemod" ]; then
                echo "ERROR: cannot compute ice transport if ice file is set to ${FILE_ICE_SUFFIX} !"; exit
            fi

            check_if_file ${TRANSPORT_ICE_SECTION_FILE}

            cp ${TRANSPORT_ICE_SECTION_FILE} ./transport_ice.dat

            diro=${DIAG_D}/transport_sections ; mkdir -p ${diro}

            echo "CALLING: /home/x_laubr/DEV/CDFTOOLS/bin/cdficeflux ${fg}"
            /home/x_laubr/DEV/CDFTOOLS/bin/cdficeflux ${fg}

            list_ice=`cat transport_ice.dat | grep '-'`

            for sect in ${list_ice}; do
                fo=${diro}/transport_solid_FW_${sect}_${CONFRUN}.dat
                mv -f section_ice-trp_${sect}.dat  ${sect}.tmp
                echo "# Time       VolTrans(Sv)     (${cyear})" >> ${fo}
                cat ${sect}.tmp | grep -v '\#' | awk -v "Y=${jyear}" '{print Y+($1-0.5)*1./12.," ",$2}'>> ${fo}
                cat ${sect}.tmp | grep -v '\#' | awk '{print $2}'> ${sect}_1.tmp
                # Mean val for the current year:
                fo=${diro}/transport_solid_FW_${sect}_${CONFRUN}_annual.dat
                mean_val1=`cat ${sect}_1.tmp | awk '{ SUM += $1} END { printf("%.15g\n", SUM/12) }'`
                ymid=`echo ${jyear} | awk  '{print $1+0.5}'`
                if [ ${jyear} -eq ${YEAR_INI} ]; then
                    echo "# Annually-averaged transports"  >> ${fo}
                    echo "# time       VolTrans(Sv)" >> ${fo}
                fi
                echo "${ymid}   ${mean_val1}" >> ${fo}
            done
            echo; echo; echo
        fi












        # DIAGS ARE DONE !!!

        rm -f *.tmp

        #debug:
        rm -f ${CRT1}_*.nc

        echo "${cyear}" > ${fcompletion}


    fi ; # if ${lcontinue}; then


    if ${LFORCE_END}; then
        if [ ${jyear} -eq ${YEARN} ]; then lcontinue=false; fi
    fi

    jyear=`expr ${jyear} + 1`



# end loop years...
done ; # while ${lcontinue}; do







# Combining files that need to be combined!


#yr1=`cat ${DIAG_D}/${CPREF}progression.dat | grep -v '\#' | head -1 | cut -c1-4`
#yr2=`cat ${DIAG_D}/${CPREF}progression.dat | grep -v '\#' | tail -1 | cut -c1-4`
#echo; echo; echo "yr1, yr2 = ${yr1} and ${yr2}"; echo


#if [ ${i_do_amo} -eq 1 ]; then
#    # Combining all the sigma files into 1 big file!
#    diro=${DIAG_D}/amo
#    l0=`\ls ${diro}/AMO_SST_Atl_${CONFRUN}_*.nc 2>/dev/null`
#    if [ ! "${l0}" = "" ]; then
#        rm -f ${DIAG_D}/AMO_SST_Atl_*.nc
#        echo
#        echo "ncrcat -h -O ${diro}/AMO_SST_Atl_${CONFRUN}_*.nc -o \
#            ${DIAG_D}/AMO_SST_Atl_${CONFRUN}_${yr1}0101-${yr2}1231.nc"
#        ncrcat -h -O ${diro}/AMO_SST_Atl_${CONFRUN}_*.nc -o \
#            ${DIAG_D}/AMO_SST_Atl_${CONFRUN}_${yr1}0101-${yr2}1231.nc
#        echo; echo
#        rm -f ${diro}/AMO_SST_Atl_${CONFRUN}_*.nc
#    fi
#fi





#if [ ${i_do_sect} -eq 1 ]; then
#    # Combining all the section files into 1 big file!
#    diro=${DIAG_D}/sections
#    for cs in ${VSECT_NM[*]}; do
#        l0=`\ls ${diro}/section_*_grid_T_${cs}.nc 2>/dev/null`
#        if [ ! "${l0}" = "" ]; then
#            rm -f ${DIAG_D}/section_*_grid_T_${cs}.nc
#            echo "ncrcat -h -O ${diro}/section_*_grid_T_${cs}.nc -o \
#                ${DIAG_D}/section_${CONFRUN}_${yr1}0101-${yr2}1231_grid_T_${cs}.nc"
#            ncrcat -h -O ${diro}/section_*_grid_T_${cs}.nc -o \
#                ${DIAG_D}/section_${CONFRUN}_${yr1}0101-${yr2}1231_grid_T_${cs}.nc
#            rm -f ${diro}/section_*_grid_T_${cs}.nc
#            echo
#        fi
#    done
#    echo
#
#fi



l_pclim=false


# PREPARING HTML PAGE
# ~~~~~~~~~~~~~~~~~~~

if [ ${ISTAGE} -eq 2 ]; then

    rm -rf ${DIAG_D}/${RUN}
    
    if [ ${IFREQ_SAV_YEARS} -gt 1 ]; then
        fnamelist=namelist.${cy1m}-${cy2m}
    else
        fnamelist=namelist.${cy2}
    fi


    # Agreement between last year from output files and 'fcompletion' file:
    ydum=`cat ${fcompletion}`
    if [ ! ${ydum} -eq ${YEAR_END} ]; then
        echo;
        echo "###################################################################"
        echo "PROBLEM: in ${fcompletion} last_year = ${ydum}"
        echo "         and from stored files files last_year = ${YEAR_END} !"
        echo "         (maybe you're calling $0 with the '-f' swicth...)"
        echo "        => forcing YEAR_END to ${ydum}"
        echo "###################################################################"
        echo
        sleep 4
        export YEAR_END=${ydum}
    fi


    cd ${BARAKUDA_ROOT}/

    echo; echo; echo "RUN ${RUN}: creating plots"; echo



#    if false ; then


    # Getting list of hydrographic sections:
    #list=`cat ${BARAKUDA_ROOT}/data/transportiz.dat | grep '-'`


    # 1D plots to perform
    # ~~~~~~~~~~~~~~~~~~~

    DIAG_1D_LIST=""

    if [ ${i_do_mean} -eq 1 ]; then
        DIAG_1D_LIST="${DIAG_1D_LIST} 3d_so mean_sos 3d_thetao mean_tos mean_zos mean_mldr10_1"
    fi
    if [ ${i_do_amoc} -eq 1 ]; then DIAG_1D_LIST="${DIAG_1D_LIST} amoc";        fi
    if [ ${i_do_trsp} -gt 0 ]; then
        DIAG_1D_LIST="${DIAG_1D_LIST} transport_sections"
    fi
    if [ ${i_do_ice}  -eq 1 ]; then DIAG_1D_LIST="${DIAG_1D_LIST} seaice";        fi

    dy=`expr ${YEAR_END} - ${YEAR_INI} + 1` ; export YF2=`expr ${YEAR_END} + 1`





    # Doing 1D plots
    # ~~~~~~~~~~~~~~

    cd ${DIAG_D}/


    echo ; echo; echo "Going to perform the following 1D plots:"
    echo "    => ${DIAG_1D_LIST}"; echo

    for fd in ${DIAG_1D_LIST}; do
        ${PYTH} ${PYBRKD_EXEC_PATH}/plot_time_series.py ${fd} ; echo
    done
    echo ; echo

    if [ ${i_do_mean} -eq 1 ]; then

         # 5-month-running mean SST anomaly on Nino region 3.4 graph:
        echo "CALLING: enso.py Nino34_${CONFRUN}.dat"
        ${PYTH} ${PYBRKD_EXEC_PATH}/plot_enso.py Nino34_${CONFRUN}.nc
        echo; echo; echo
        
        # Hovmuller of temperature and salinity
        echo "CALLING: hovm_tz.py ${YEAR_INI} ${YEAR_END} ${NBL}"
        ${PYTH} ${PYBRKD_EXEC_PATH}/plot_hovm_tz.py
        echo; echo; echo
        #
    fi


    if [ ${i_do_sigt} -eq 1 ]; then
        # Transport by sigma-class
        # ~~~~~~~~~~~~~~~~~~~~~~~~
        echo "CALLING: plot_trsp_sigma.py"
        ${PYTH} ${PYBRKD_EXEC_PATH}/plot_trsp_sigma.py
        echo; echo; echo
    fi




    if [ ${i_do_mht} -eq 1 ]; then
        #
        # Hovmullers of advective meridional heat/salt transport
        echo; echo
        echo "CALLING: plot_hovm_merid_trsp.py"
        ${PYTH} ${PYBRKD_EXEC_PATH}/plot_hovm_merid_trsp.py
        echo; echo; echo
        #
    fi

    echo; echo




    if ${l_clim_diag} ; then

        ###########################################################################
        # Climatology over X years (12-month file average of X consecutive years)
        #   => has to be built with the 'build_clim.sh' script
        ###########################################################################

        echo; echo; echo "Checking for presence of ${DIAG_D}/clim/last_clim..."
        if [ -f ${DIAG_D}/clim/last_clim ]; then
            cat ${DIAG_D}/clim/last_clim
            export CLIM_PER=`cat ${DIAG_D}/clim/last_clim`
            ftcli=${DIAG_D}/clim/mclim_${CONFRUN}_${CLIM_PER}_grid_T.nc4
            ficli=${DIAG_D}/clim/mclim_${CONFRUN}_${CLIM_PER}_${FILE_ICE_SUFFIX}.nc4
            fclvt=${DIAG_D}/clim/aclim_${CONFRUN}_${CLIM_PER}_VT.nc4
            fcmoc=${DIAG_D}/clim/aclim_${CONFRUN}_${CLIM_PER}_MOC.nc4
            fcpsi=${DIAG_D}/clim/mclim_${CONFRUN}_${CLIM_PER}_PSI.nc4
            iclyear=`echo ${CLIM_PER} | sed -e s/'-'/' '/g`
        else
            echo; echo "PROBLEM! => you set l_clim_diag to true but no file 'last_clim' was found in:"
            echo "            ${DIAG_D}/clim/"; echo
            exit
        fi

        echo; echo; echo "Checking for presence of ${ftcli}..."
        if [ -f ${ftcli} ]; then
            echo; echo;
            echo "~~~~~~~~~~~~~~~~~~~~~"
            echo "*   CLIMATO FOUND !!!"
            echo "~~~~~~~~~~~~~~~~~~~~~"
            echo
            echo "  => for years ${CLIM_PER}" ; echo "  => using ${ftcli}"

            list_comp_2d="CLIM"
            l_pclim=true
            lcomp_to_run=false

            if [ ! "${RUNREF}" = "" ]; then
                lcomp_to_run=true
                list_comp_2d="CLIM ${RUNREF}"
                # Must check if climatology for run ${RUNREF} is there:
                fclim_ref=`echo "${ftcli}" | sed -e "s|${RUN}|${RUNREF}|g"`
                check_if_file ${fclim_ref}
                echo "Going to compare also against run ${fclim_ref}!"
                echo
            fi


            echo; echo
            check_if_file ${F_T_CLIM_3D_12} "name:F_T_CLIM_3D_12"
            check_if_file ${F_S_CLIM_3D_12} "name:F_S_CLIM_3D_12"
            check_if_file ${SST_CLIM_12}    "name:SST_CLIM_12"            
            if [ ${i_do_ice}  -gt 0 ]; then check_if_file ${ICE_CLIM_12}    "name:ICE_CLIM_12" ; fi
            echo; echo


        #######################################
        # Diags that don't imply a comparison #
        #######################################

            export COMP2D="CLIM"



        # Lat-Depth AMOC
        # ~~~~~~~~~~~~~~
            if [ -f ${fcmoc} ]; then
                echo; echo
                echo " Ploting lat-depth MOC !"
                cd ${DIAG_D}/
                DIRS_2_EXP="${DIRS_2_EXP} moc"
                rm -rf moc; mkdir moc; cd moc/
                echo; echo; echo "CALLING: moc.py ${iclyear}"; echo
                ${PYTH} ${PYBRKD_EXEC_PATH}/moc.py ${iclyear}
                cd ../
                echo
            fi


            # March Mixed layer depth in Nordic Seas
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if [ `ipresent_var_in_ncf ${ftcli} ${NN_MLD}` -eq 1 ]; then
                echo; echo
                echo " Performing 2D mapping of March Mixed layer depth in Nordic Seas"
                cd ${DIAG_D}/
                DIRS_2_EXP="${DIRS_2_EXP} mld"
                rm -rf mld; mkdir mld; cd mld/
                echo; echo; echo "CALLING: mld.py ${iclyear}"; echo
                ${PYTH} ${PYBRKD_EXEC_PATH}/mld.py ${iclyear}
                cd ../
                echo
            fi

            # Sea-ice extent stereographic polar projection South and North
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if [ ${i_do_ice}  -gt 0 -a `ipresent_var_in_ncf ${ficli} ${NN_ICEF}` -eq 1 ]; then
                echo; echo
                echo " Performing 2D Sea-ice extent stereographic polar projection South and North"
                cd ${DIAG_D}/
                DIRS_2_EXP="${DIRS_2_EXP} sea_ice"
                rm -rf sea_ice; mkdir sea_ice; cd sea_ice/
                echo; echo; echo "CALLING: ice.py ${iclyear}"; echo
                ${PYTH} ${PYBRKD_EXEC_PATH}/ice.py ${iclyear}
                cd ../
                echo
            fi



        # Sea-surface height
        # ~~~~~~~~~~~~~~~~~~
            if [ `ipresent_var_in_ncf ${ftcli} ${NN_SSH}` -eq 1 ]; then
                echo; echo
                echo " SSH map"
                DIRS_2_EXP="${DIRS_2_EXP} ssh"
                cd ${DIAG_D}/
                rm -rf ssh; mkdir ssh; cd ssh/
                echo; echo; echo "CALLING: ssh.py ${iclyear}"; echo
                ${PYTH} ${PYBRKD_EXEC_PATH}/ssh.py ${iclyear}
                cd ../
                echo
            else
                echo; echo "WARNING: did not find ${NN_SSH} into ${ftcli} !!!!"; echo
            fi




        ##################################################
        # Diags that imply a comparison against "COMP2D" #
        ##################################################

            cd ${DIAG_D}/
            rm -rf temp_sal surf_fluxes

            for COMP2D in ${list_comp_2d}; do

                export COMP2D=${COMP2D}

                echo; echo; echo "Clim. comparisons against ${COMP2D}"

                if [ "${COMP2D}" = "${RUNREF}" ]; then
                    export F_T_CLIM_3D_12=${fclim_ref}; check_if_file ${F_T_CLIM_3D_12} "name:F_T_CLIM_3D_12"
                    export F_S_CLIM_3D_12=${fclim_ref}; check_if_file ${F_S_CLIM_3D_12} "name:F_S_CLIM_3D_12"
                    export SST_CLIM_12=${fclim_ref}   ; check_if_file ${SST_CLIM_12}    "name:SST_CLIM_12"
                    if [ ${i_do_ice}  -gt 0 ]; then export ICE_CLIM_12=${fclim_ref}   ; check_if_file ${ICE_CLIM_12}    "name:ICE_CLIM_12"; fi
                fi


            # Temperature and Salinity
            # ~~~~~~~~~~~~~~~~~~~~~~~~
                echo; echo
                echo " Performing 2D Temperature and Salinity maps and sections"
                cd ${DIAG_D}/
                DIRS_2_EXP="${DIRS_2_EXP} temp_sal"
                DIRS_2_EXP_RREF="${DIRS_2_EXP_RREF} temp_sal"
                mkdir -p temp_sal; cd temp_sal/
                echo; echo; echo "CALLING: temp_sal.py ${iclyear}"; echo
                ${PYTH} ${PYBRKD_EXEC_PATH}/temp_sal.py ${iclyear}
                cd ../
                echo


            # Surface fluxes diagnostics
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~
            #cd ${DIAG_D}/
            #if [ ${i_do_flx} -eq 1 ]; then
            #    COMFLUX1='' ; COMFLUX2='' ; DIRS_2_EXP="${DIRS_2_EXP} surf_fluxes"
            #    mkdir -p surf_fluxes; cd surf_fluxes/
            #    echo; echo; echo "CALLING: surf_fluxes.py ${iclyear}"; echo
            #    ${PYTH} ${PYBRKD_EXEC_PATH}/surf_fluxes.py ${iclyear}
            #    echo; echo; echo
            #    cd ../
            #fi
            #echo



            done ; # for COMP2D in ${list_comp_2d}; do




        #if [ ${i_do_mht} -eq 1 ]; then
        ## Meridional transports of heat and salt
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #    if [ -f ${fclvt} ]; then
        #        cd ${TMP_DIR}/
        #        echo "in `pwd` !"; ls ; echo
        #        echo "${BARAKUDA_ROOT}/cdftools_light/bin/cdfmhst.x ${fclvt}"
        #        ${BARAKUDA_ROOT}/cdftools_light/bin/cdfmhst.x ${fclvt}
        #            # mhst.nc, merid_heat_trp.dat and merid_salt_trp.dat created...
        #        rm -f mhst.nc
        #        mv -f  merid_heat_trp.dat ${DIAG_D}/mht_clim_${CONFRUN}_${CLIM_PER}.dat
        #        mv -f  merid_salt_trp.dat ${DIAG_D}/mst_clim_${CONFRUN}_${CLIM_PER}.dat
        #            #
        #        cd ${BARAKUDA_ROOT}/
        #        ./scripts/mk_plots.sh mht_clim
        #        ./scripts/mk_plots.sh mst_clim
        #    else
        #        echo "PROBLEM: ${fclvt} not found"
        #        echo "   => skipping meridional heat and salt transport diag..."
        #        echo
        #    fi
        #
        #fi



        else
            echo; echo
            echo " NO CLIMATO FOUND ...";
            echo "   => you can use 'build_clim.sh' to build a climato of your run"
            echo; echo
        fi

    fi ; # if ${l_clim_diag}


    cd ${DIAG_D}/






    
    
    
    
    echo; echo; echo ; echo "Creating HTML file!"
    
    cd ${DIAG_D}/
    
    inmlst=0
    if [ -f ${NEMO_OUT_D}/${fnamelist} ]; then inmlst=1; fi
    
    rm -f index.php
    
    if [ "${JTITLE}" = "" ]; then
        echo "Problem, variable JTITLE is not set!" ; exit
    else
        TITLE="Ocean diagnostics, run ${RUN}, conf ${ORCA}, ${JTITLE}"
    fi
    

    sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${CONFRUN}|g" -e "s|{DATE}|`date`|g" -e "s|{HOST}|`hostname`|g" \
        ${BARAKUDA_ROOT}/scripts/html/index_skel_header.html > index.php

    # Namelist section
    if [ ${inmlst} -eq 1 ]; then
        cat >> index.php <<EOF
        <big> <a href="./namelist.html"> Last namelist </a> </big>
        <br><br><br>
EOF
    fi
    
    # Climato section
    if ${l_pclim}; then
        cat >> index.php <<EOF
        <br><big><big> Diags from climatology (${CLIM_PER}) </big></big><br><br>
        <big> <a href="./temp_sal/index.php"> Temperature and Salinity vs CLIM</a> </big>             <br><br><br>
        <big> <a href="./ssh/index.php">  Sea Surface Height </a> </big>                              <br><br><br>
        <big> <a href="./sea_ice/index.php">  Arctic and Antarctic sea-ice extent vs CLIM </a> </big> <br><br><br>
        <big> <a href="./mld/index.php">  Mixed Layer depth in relevent regions </a> </big>           <br><br><br>
        <big> <a href="./moc/index.php">  Meridional Overturning Circulation </a> </big>              <br><br><br>
EOF
        if ${lcomp_to_run}; then
            cat >> index.php <<EOF
        <br><big><big> Comparison with run ${RUNREF}, climatology (2004-2007) </big></big><br><br>
        <big> <a href="./temp_sal/index_${RUNREF}.html"> Temperature and Salinity vs ${RUNREF}</a> </big>             <br><br><br>
<!--        <big> <a href="./ssh/index_${RUNREF}.html">  Sea Surface Height </a> </big>                              <br><br><br>
        <big> <a href="./sea_ice/index_${RUNREF}.html">  Arctic and Antarctic sea-ice extent vs ${RUNREF} </a> </big> <br><br><br>
        <big> <a href="./mld/index_${RUNREF}.html">  Mixed Layer depth in relevent regions </a> </big>           <br><br><br>
        <big> <a href="./moc/index_${RUNREF}.html">  Meridional Overturning Circulation </a> </big>              <br><br><br>
-->
EOF
        fi
    fi

    # Temperature section     
    cat >> index.php <<EOF
        <br><br><br><big><big> Temperature time-series </big></big><br><br>
        <img style="border: 0px solid" alt="" src="3d_thetao_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="mean_tos_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="3d_thetao_lev_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="3d_thetao_basins_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="Nino34_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="hov_temperature_${CONFRUN}_global.png"> <br><br>
        <img style="border: 0px solid" alt="" src="hov_temperature_${CONFRUN}_atlantic.png"> <br><br>
        <img style="border: 0px solid" alt="" src="hov_temperature_${CONFRUN}_pacific.png"> <br><br>
        <img style="border: 0px solid" alt="" src="hov_temperature_${CONFRUN}_indian.png"> <br><br><br><br>
EOF

    # Salinity section     
    cat >> index.php <<EOF
        <br><br><br><big><big> Salinity time-series </big></big><br><br>
        <img style="border: 0px solid" alt="" src="3d_so_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="mean_sos_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="3d_so_lev_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="3d_so_basins_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="hov_salinity_${CONFRUN}_global.png"> <br><br>
        <img style="border: 0px solid" alt="" src="hov_salinity_${CONFRUN}_atlantic.png"> <br><br>
        <img style="border: 0px solid" alt="" src="hov_salinity_${CONFRUN}_pacific.png"> <br><br>
        <img style="border: 0px solid" alt="" src="hov_salinity_${CONFRUN}_indian.png"> <br><br> <br><br>
EOF

    # MISC section     
    cat >> index.php <<EOF
        <br><br><br><big><big> Misc. time-series </big></big><br><br>
        <img style="border: 0px solid" alt="" src="mean_zos_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="amoc_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="amoc_${CONFRUN}_comp.png"> <br><br> <br><br> 
EOF

    # Sea-ice section     
    if [ ${i_do_ice}  -gt 0 ]; then
        cat >> index.php <<EOF
        <br><br><br><big><big> Arctic/Antarctic sea-ice time-series</big></big><br><br>
        <img style="border: 0px solid" alt="" src="seaice_extent_winter_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="seaice_extent_summer_${CONFRUN}.png"> <br><br><br>
        <img style="border: 0px solid" alt="" src="seaice_volume_winter_${CONFRUN}.png"> <br><br>
        <img style="border: 0px solid" alt="" src="seaice_volume_summer_${CONFRUN}.png"> <br><br><br><br>
EOF
    fi

    if [ ${i_do_trsp} -gt 0 ]; then
        # Adding transport section part:
        echo "<br><br><br><big><big> Transport through sections</big></big><br><br>" >> index.php
        list_section=`cat ${TRANSPORT_SECTION_FILE} | grep '-'`
        for cs in ${list_section}; do
            echo ${cs}
            echo "<img style=\"border: 0px solid\" alt=\"\" src=\"transport_vol_${cs}_${CONFRUN}.png\"> <br>"  >> index.php
            echo "<img style=\"border: 0px solid\" alt=\"\" src=\"transport_heat_${cs}_${CONFRUN}.png\"> <br><br>" >> index.php
            echo "<br><br>" >> index.php
        done
    fi
        
    # Checking if figures with time-series of MLD in specified boxes are here and adding them:
    if [ ${i_do_mean} -eq 1 ]; then
        list_mld_figs=`\ls mean_mldr10_1_${CONFRUN}*.png`
        if [ ! "${list_mld_figs}" = "" ]; then
            echo "<br><br><br><big><big> Horizontally-averaged Mixed-Layer Depth in different regions</big></big><br><br>" >> index.php
            for fmld in ${list_mld_figs}; do
                echo "<img style=\"border: 0px solid\" alt=\"\" src=\"${fmld}\"> <br>"  >> index.php
                echo "<br><br>" >> index.php
            done
        fi
    fi

    if [ ${i_do_sigt} -eq 1 ]; then
        # Adding transport by sigma class section part:
        echo "<br><br><br><big><big> Transport by sigma class at Nordic sills</big></big><br><br>" >> index.php
        list_section=`cat ${DENSITY_SECTION_FILE} | grep '_'`
        for cs in ${list_section}; do
            echo ${cs}
            echo "<img style=\"border: 0px solid\" alt=\"\" src=\"transport_sigma_class_${cs}_${CONFRUN}.png\"> <br>"  >> index.php
        done
        echo "<img style=\"border: 0px solid\" alt=\"\" src=\"tr_sigma_gt278_${CONFRUN}.png\"> <br>"  >> index.php
        echo "<br><br>" >> index.php
    fi
    
    if [ ${i_do_mht} -eq 1 ]; then
        # Adding meridional heat transport:
        echo "<br><br><br><big><big> Meridional transports</big></big><br><br>"  >> index.php
        for coce in "global" "atlantic" "pacific" "indian"; do
            echo "<img style=\"border: 0px solid\" alt=\"\" src=\"MHT_${CONFRUN}_${coce}.png\"> <br>"     >> index.php
            echo "<img style=\"border: 0px solid\" alt=\"\" src=\"MST_${CONFRUN}_${coce}.png\"> <br><br>" >> index.php
        done
        echo "<br><br>" >> index.php
    fi

    cat ${BARAKUDA_ROOT}/scripts/html/index_skel_footer.html >> index.php ; # Closing HTML file...



    # If climatology built, sub 2D html pages
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ${l_pclim}; then

       # T, S, SSH and ice HTML page:
        for cdiag in ${DIRS_2_EXP}; do
            sed -e "s|{CONFRUN}|${CONFRUN}|g" -e "s|{COMP2D}|CLIM|g" \
                ${BARAKUDA_ROOT}/scripts/html/${cdiag}/index_${cdiag}.html > ${cdiag}/index.php
        done
        for var in "sst" "sss" "sections_ts" "ts_100m" "ts_1000m" "ts_3000m"; do
            sed -e "s|{CONFRUN}|${CONFRUN}|g" -e "s|{COMP2D}|CLIM|g" \
                ${BARAKUDA_ROOT}/scripts/html/temp_sal/${var}.html > temp_sal/${var}_CLIM.php
        done
        
        if ${lcomp_to_run}; then
            for cdiag in ${DIRS_2_EXP_RREF}; do
                sed -e "s|{CONFRUN}|${CONFRUN}|g" -e "s|{COMP2D}|${RUNREF}|g" \
                    ${BARAKUDA_ROOT}/scripts/html/${cdiag}/index_${cdiag}.html > ${cdiag}/index_${RUNREF}.php
            done
            for var in "sst" "sss" "sections_ts" "ts_100m" "ts_1000m" "ts_3000m"; do
                sed -e "s|{CONFRUN}|${CONFRUN}|g" -e "s|{COMP2D}|${RUNREF}|g" \
                    ${BARAKUDA_ROOT}/scripts/html/temp_sal/${var}.html > temp_sal/${var}_${RUNREF}.php
            done
        fi
                
        # Surface fluxes HTML page:
        if [ ${i_do_flx} -eq 1 ]; then
            sed -e "s|{CONFRUN}|${CONFRUN}|g" -e "s|{COMP2D}|${COMP2D}|g" \
                ${BARAKUDA_ROOT}/scripts/html/surf_fluxes/index_fluxes.html > surf_fluxes/index.php
            for flx in "Qsw" "Qlw" "Qla" "Qse"; do
                sed -e "s|{CONFRUN}|${CONFRUN}|g" -e "s|{COMP2D}|${COMP2D}|g" \
                    ${BARAKUDA_ROOT}/scripts/html/surf_fluxes/${flx}.html > surf_fluxes/${flx}.php
            done
        fi

    fi

    echo "Done!"; echo; echo

    echo; echo; echo

    if [ ${ihttp} -eq 1 ]; then
        echo "Preparing to export to remote host!"; echo
        mkdir ${RUN}

        cp -r ${BARAKUDA_ROOT}/scripts/html/conf_*.html ${RUN}/
        mv -f index.php ${RUN}/
        mv -f *.png      ${RUN}/ >/dev/null 2>/dev/null
        mv -f ./merid_transport/*.png ${RUN}/ >/dev/null 2>/dev/null
        mv -f *.svg      ${RUN}/ >/dev/null 2>/dev/null

        if [ ${inmlst} -eq 1 ]; then cp ${NEMO_OUT_D}/${fnamelist} ${RUN}/; fi

        cp -r ${DIRS_2_EXP} ${RUN}/ >/dev/null 2>/dev/null

        tar cvf ${RUN}.tar ${RUN}
        ssh ${RUSER}@${RHOST} "mkdir -p ${RWWWD}"
        scp ${RUN}.tar ${RUSER}@${RHOST}:${RWWWD}/
        ssh ${RUSER}@${RHOST} "cd ${RWWWD}/; rm -rf ${RUN}; tar xf ${RUN}.tar 2>/dev/null; rm -f ${RUN}.tar; \
            chmod -R a+r ${RUN}; cd ${RUN}/; source-highlight -i ${fnamelist} -s fortran -o namelist.html"
        echo; echo
        echo "Diagnostic page installed on  http://${RHOST}${RWWWD}/${RUN}/ !"
        echo "( Also browsable on local host in ${DIAG_D}/${RUN} )"
        rm -f ${RUN}.tar

    else

        if [ ${ihttp} -eq 0 ]; then
            echo "Diagnostic page installed in ${DIAG_D}/"
            echo " => view this directory with a web browser..."
        else
            echo "Error: \"ihttp\" is either 0 or 1 !"
        fi

    fi

    echo; echo

    rm -rf *.eps


else
    echo; echo "Diagnostics built! Run \"${0} -e\" to create figure and HTML page..."; echo


fi


#debug:
rm -rf ${TMP_DIR} 2>/dev/null

echo

