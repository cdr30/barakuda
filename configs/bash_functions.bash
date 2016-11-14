#!/usr/bin/env bash

function barakuda_usage()
{
    echo
    echo "USAGE: ${0} -C <config> -R <run>  (options)"
    echo
    echo "     Available configs are:"
    for cc in ${list_conf}; do
        echo "         * ${cc}"
    done
    echo
    echo
    echo "   OPTIONS:"
    echo
    echo "      -f <years> => how many years per NEMO file? default = 1"
    echo
    echo "      -y <YYYY> => force initial year to YYYY"
    echo
    echo "      -F        => forces creation of diagnostic page even if treatment of output files is not finished"
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
}


function barakuda_init()
{
    # Supported ORCA grids:
    export ORCA_LIST="ORCA1.L75 ORCA1.L46 ORCA1.L42 ORCA2 ORCA2_L46"

    # Some defaults:
    export LFORCE_YINI=false
    export LFORCE_YEND=false
    export RUNREF=""
    export ISTAGE=1 ; # 1 => generation of data diagnostic files
    #          # 2 => creation of figures and diagnostic HTML page
    export LFORCEDIAG=false
    export l_clim_diag=false
    export IFREQ_SAV_YEARS=1
    #
}

function barakuda_check()
{
    script=`basename $0 | sed -e s/'.sh'/''/g`
    if [ "${CONFIG}" = "" ] || [ "${RUN}" = "" ]; then ${script}_usage ; exit ; fi

    if [ "${RUNREF}" != "" ] && [ ${ISTAGE} -eq 1 ]; then
        echo; echo " WARNING: option '-c' only makes sense when '-e' or '-E' are specified !"
        sleep 2; echo
    fi

    for og in ${ORCA_LIST}; do
        echo " ${og} / ${CONFIG}"
        ca=""; ca=`echo ${CONFIG} | grep ${og}` ; if [ "${ca}" != "" ]; then ORCA=${og}; fi
    done

    if [ "${ORCA}" = "" ]; then echo "ORCA grid of config ${CONFIG} not supported yet"; exit; fi
    echo
    
    if [ "${script}" = "build_clim" ]; then
        echo Boo
        if [ -z ${Y1} ] || [ -z ${Y2} ]; then
            ${script}_usage
            echo; echo " ==> Please specify both first and last year for climatology (-i and -e ) !!!"; echo
            exit
        fi
    fi
}


function barakuda_setup()
{
    script=`basename $0 | sed -e s/'.sh'/''/g`
    echo
    if [ ! "${ORCA}" = "${CONF}" ]; then echo "ERROR: ORCA and CONF disagree! => ${ORCA} ${CONF}"; exit; fi
    export ORCA=${CONF}
    echo

    if [ "${PYTHON_HOME}" = "" ]; then echo "ERROR: PYTHON_HOME is not set! => add it to config file"; exit; fi
    export PYTH="${PYTHON_HOME}/bin/python -W ignore" ; # which Python installation to use
    export PYTHONPATH=${PYTHON_HOME}/lib/python2.7/site-packages:${BARAKUDA_ROOT}/python/modules ; # PATH to python barakuda modules
    export PYBRKD_EXEC_PATH=${BARAKUDA_ROOT}/python/exec         ; # PATH to python barakuda executable

    echo " PYTHON_HOME => "${PYTHON_HOME} ; echo
    echo "  TRANSPORT_SECTION_FILE => ${TRANSPORT_SECTION_FILE} !" ; echo

    if ${l_clim_diag} ; then
        echo
        echo " Files containing climatologies to be used:"
        echo " T 3D => ${F_T_CLIM_3D_12} => ${NN_T_CLIM}"
        echo " S 3D => ${F_S_CLIM_3D_12} => ${NN_S_CLIM}"
        echo " SST  => ${SST_CLIM_12} => ${NN_SST_CLIM}"
        echo
        for ff in ${F_T_CLIM_3D_12} ${F_S_CLIM_3D_12} ${SST_CLIM_12}; do
            if [ ! -f ${F_T_CLIM_3D_12} ]; then echo "ERROR: ${ff} is missing!"; exit; fi
        done
    fi

# Names for temperature, salinity, u- and v-current...
    if [ "${NN_T}" = "" ] || [ "${NN_S}" = "" ] || [ "${NN_U}" = "" ] || [ "${NN_V}" = "" ]; then
        echo "NN_T, NN_S, NN_U and NN_V are NOT given a value into"
        echo " in ${fconfig} "
        echo "  => using default names: thetao, so, uo, vo" ; echo
        NN_T="thetao"; NN_S="so"; NN_U="uo"; NN_V="vo"
    fi

    echo ; echo " *** NN_T=${NN_T}, NN_S=${NN_S}, NN_U=${NN_U} and NN_V=${NN_V} "; echo

    # Checking what files we have / plan to use:
    if [ "${NEMO_SAVED_FILES}" = "" ]; then
        echo "Please specify which NEMO files are saved (file suffixes, grid_T, ..., icemod) ?"
        echo " => set the variable NEMO_SAVED_FILES in your config_${CONFIG}.sh file!"; exit
    fi
    echo; echo "File types to import (NEMO_SAVED_FILES) : ${NEMO_SAVED_FILES}"; echo; echo

    # Need to be consistent with the netcdf installation upon which cdftools_light was compiled:
    export NCDF_DIR=`cat cdftools_light/make.macro | grep ^NCDF_DIR | cut -d = -f2 | sed -e s/' '//g`
    echo ; echo "NCDF_DIR = ${NCDF_DIR}"; echo
    export LD_LIBRARY_PATH=${NCDF_DIR}/lib:${LD_LIBRARY_PATH}

    # Exporting some variables needed by the python scripts:
    export RUN=${RUN}
    export CONFRUN=${ORCA}-${RUN}
    export DIAG_D=${DIAG_DIR}/${CONFRUN}
    export CLIM_DIR=${DIAG_D}/clim
    
    if [ ${ISTAGE} -eq 1 ]; then
        # We need a scratch/temporary directory to copy these files to and gunzip them:
        if [ "${SLURM_JOBID}" = "" -a `hostname` = "triolith1" ]; then
        # Likely to be running interactively on triolith
            export SCRATCH=${HOME}/tmp/${script}_tmp
            export TMP_DIR=${SCRATCH}/${RUN}_${script}_tmp
        else
            # Normal case:
            SCRATCH=`echo ${SCRATCH} | sed -e "s|<JOB_ID>|${SLURM_JOBID}|g"`
            export TMP_DIR=${SCRATCH}/${RUN}
        fi
        echo " IMPORTANT the SCRATCH work directory is set to:" ; echo " ${SCRATCH}"
    else
        export TMP_DIR=${DIAG_D}/${script}_tmp
    fi

    echo " IMPORTANT the TMP_DIR work directory is set to:" ; echo " ${TMP_DIR}"; echo ; sleep 2

    rm -rf ${TMP_DIR}
    mkdir -p ${DIAG_D} ${TMP_DIR}

    export NEMO_OUT_D=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g"`
    if [ ! -d ${NEMO_OUT_D} ]; then echo "Unfortunately we could not find ${NEMO_OUT_D}"; exit; fi

    echo; echo " * Config to be used: ${CONFIG} => ORCA grid is ${ORCA}"
    echo " * Run is ${RUN} "; echo " * Files are stored into ${NEMO_OUT_D}"; echo; sleep 2

    export CPREF=`echo ${NEMO_FILE_PREFIX} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g" -e "s|<TSTAMP>|${TSTAMP}|g"`
    echo " NEMO files prefix = ${CPREF} "

    # only neede for barakuda.sh :
    if [ "${script}" = "barakuda" ]; then
        # Testing if NCO is installed:
        which ncks 1>out 2>/dev/null; ca=`cat out`; rm -f out
        if [ "${ca}" = "" ]; then echo "Install NCO!!!"; echo; exit; fi

        # Not fully supported yet:
        ca=" => diagnostic totally beta and not fully supported yet!"
        if [ ${i_do_sect} -gt 0 ]; then echo " *** i_do_sect ${ca}"; exit; fi
        if [ ${i_do_amo}  -gt 0 ]; then echo " *** i_do_amo  ${ca}"; exit; fi
        if [ ${i_do_icet} -gt 0 ]; then echo " *** i_do_icet ${ca}"; exit; fi
        if [ ${i_do_flx}  -gt 0 ]; then echo " *** i_do_flx  ${ca}"; exit; fi

        if [ ${ISTAGE} -eq 1 ]; then
            # List of CDFTOOLS executables needed for the diagnostics:
            L_EXEC="cdfmaxmoc.x cdfmoc.x cdfvT.x cdftransportiz.x cdficediags.x cdfmhst.x cdfsigtrp.x"
            for ex in ${L_EXEC}; do check_if_file cdftools_light/bin/${ex} "Compile CDFTOOLS executables!"; done
        fi
    fi
    
    # only needed for build_clim.sh :
    if [ "${script}" = "build_clim" ]; then
        mkdir -p ${CLIM_DIR}
    fi
    
    echo
}

function barakuda_first_last_years()
{
    cd ${NEMO_OUT_D}/
    if [ ${ece_run} -gt 0 ]; then
        if [ ! -d 001 ]; then echo "ERROR: since ece_run=${ece_run}, there should be a directory 001 in:"; echo " ${NEMO_OUT_D}"; exit ; fi
        nby_ece=`ls -d */ | wc -l` ; echo " ${nby_ece} years have been completed..."
        cd 001/
    fi

    # Try to guess the first year from stored "grid_T" files:
    YEAR_INI=`\ls ${CPREF}*${ctest}* | sed -e s/"${CPREF}"/""/g | head -1 | cut -c1-4`
    echo ${YEAR_INI} |  grep "[^0-9]" >/dev/null ;   # Checking if it's an integer:
    if [ ! "$?" -eq 1 ]; then
        echo "ERROR: it was imposible to guess initial year from your input files"
        echo "       maybe the directory contains non-related files..."
        exit
    fi
    export YEAR_INI=$((${YEAR_INI}+0))  ; # example: 1 instead of 0001...
    export YEAR_INI_F=${YEAR_INI} ; # saving the year deduced from first file

    if ${LFORCE_YINI}; then
        if [ ${YEAR0} -lt ${YEAR_INI_F} ]; then echo "ERROR: forced initial year is before first year!"; exit; fi
        export YEAR_INI=${YEAR0}
        echo " Initial year forced to ${YEAR_INI} !"
    fi

    cd ${NEMO_OUT_D}/

    if [ ${ece_run} -gt 0 ]; then
        dir_end=`printf "%03d" ${nby_ece}`
        if [ ! -d ${dir_end} ]; then echo "ERROR: since ece_run=${ece_run}, there should be a directory ${dir_end} in:"; echo " ${NEMO_OUT_D}"; exit ; fi
        export YEAR_END=$((${YEAR_INI}+${nby_ece}))
    else
        export YEAR_END=`\ls ${CPREF}*${ctest}* | sed -e s/"${CPREF}"/''/g | tail -1 | cut -c1-4`
        echo ${YEAR_END} |  grep "[^0-9]" >/dev/null; # Checking if it's an integer
        if [ ! "$?" -eq 1 ]; then
            echo "ERROR: it was imposible to guess the year coresponding to the last saved year!"
            echo "       => check your NEMO output directory and file naming..."; exit
        fi
        export YEAR_END=$((${YEAR_END}+${IFREQ_SAV_YEARS}-1))
    fi
    echo
    echo " *** Initial year set to ${YEAR_INI}"
    echo " ***   Last  year set to ${YEAR_END}"
    echo
}

function barakuda_init_plot()
{
    if [ ! -z ${YEAR0} ]; then
        echo
        echo "WARNING: using the -y switch for the creation of the plots (switch -e)"
        echo "         will have no impact. For the time slice to use, the plots are based"
        echo "         on what is found in ${DIAG_D}/ !"
        echo "         => the diagnostic netdcf files"
        echo "         => which should also be equivalent to what's found in files"
        echo "            first_year.info and last_year_done.info"
        echo ""
        sleep 5
    fi

    for fc in "first_year" "numb_year_per_file" "last_year_done"; do
        ff=${DIAG_D}/${fc}.info
        if [ ! -f ${ff} ]; then echo "ERROR: file ${ff} is missing!"; exit; fi
    done

    # Only doing plots...
    export YEAR_INI=`cat ${DIAG_D}/first_year.info`
    export YEAR_END=`cat ${DIAG_D}/last_year_done.info`
    export IFREQ_SAV_YEARS=`cat ${DIAG_D}/numb_year_per_file.info`

    echo
    echo " For time-series to be ploted:"
    echo "  * initial year = ${YEAR_INI}"
    echo "  *  last   year = ${YEAR_END}"
    echo
}

function barakuda_import_mesh_mask()
{
    script=`basename $0 | sed -e s/'.sh'/''/g`
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
    echo " *** mesh_mask to be used: ${MM_FILE} "
    echo " *** basin mask to be used: ${BM_FILE} "
    echo
    
    if [ "${script}" = "build_clim" ]; then
        ln -sf mesh_mask.nc mesh_hgr.nc ; ln -sf mesh_mask.nc mesh_zgr.nc ; ln -sf mesh_mask.nc mask.nc
    fi
}

function barakuda_check_year_is_complete()
{
    jy1=${jyear} ; jy2=$((${jyear}+${IFREQ_SAV_YEARS}-1))
    cy1=`printf "%04d" ${jy1}` ; cy2=`printf "%04d" ${jy2}`
    cy1m=`printf "%04d" $((${jy1}-${IFREQ_SAV_YEARS}))` ; cy2m=`printf "%04d" $((${jy2}-${IFREQ_SAV_YEARS}))`
    export i_get_file=1
    echo " *** (${cyear} => from ${cy1} to ${cy2})"
    export TTAG=${cy1}0101_${cy2}1231 # calendar-related part of the file name
        # Testing if the current year-group has been done
    for ft in ${NEMO_SAVED_FILES}; do
        if ${lcontinue}; then
            ftst=${NEMO_OUT_D}/${cpf}${CPREF}${TTAG}_${ft} ;  cfxt="0"
            for ca in "nc" "nc.gz" "nc4"; do
                if [ -f ${ftst}.${ca} ]; then cfxt="${ca}"; fi
            done
            if [ ${cfxt} = "0" ]; then
                echo "Year(s) ${cy1}-${cy2} is not completed yet:"; echo " => ${ftst}(?) is missing"; echo
                export lcontinue=false
            fi
        fi
    done
    if ${lcontinue}; then echo " *** All files for ${TTAG} are there!"; fi
    echo
}

function barakuda_import_files()
{
    # On what file type to test file presence:
    cgrid_test=`echo ${NEMO_SAVED_FILES} | cut -d ' ' -f2`
    echo " *** testing on files \"${cgrid_test}\" !"; echo
    l_happy=false
    while ! ${l_happy} ; do
        if [ ${IFREQ_SAV_YEARS} -eq 1 ]; then l_happy=true; fi
        rm -f *.tmp
        if [ ${i_get_file} -eq 1 ]; then
            echo " => gonna get ${CRTM}_* files..."
            # Importing required files to tmp dir and unzipping:
            for gt in ${NEMO_SAVED_FILES}; do
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
            # Need to create annual files if more than 1 year in 1 one NEMO file
            if [ ${IFREQ_SAV_YEARS} -gt 1 ]; then
                for gt in ${NEMO_SAVED_FILES}; do
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
        fi
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
    done
    echo
}

function build_clim_usage()
{
    echo
    echo "USAGE: ${0} -C <config> -R <run> -i <first_year> -e <last_year> (options)"
    echo
    echo "     Available configs are:"
    for cc in ${list_conf}; do
        echo "         * ${cc}"
    done
    echo
    echo "   OPTIONS:"
    echo "      -f <years> => how many years per NEMO file? default = 1"
    echo "      -h         => print this message"
    echo
}





# ================ Misc. Functions ===============



# Now a file extension:
function cext()
{
    nbc=`echo ${1} | wc -c`; nbb=`expr ${nbc} - 3`; echo "`echo ${1} | cut -c${nbb}-${nbc}`"
    # More elegant way would be following but there can be some "." into the ORCA grid...
    #echo ".`echo $1 | cut -d'.' -f2`
}





lb_is_leap()
{
    if [ "$1" = "" ]; then echo "USAGE: lb_is_leap <YEAR>"; exit; fi
    #
    i_mod_400=`expr ${1} % 400`
    i_mod_100=`expr ${1} % 100`
    i_mod_4=`expr ${1} % 4`
    #
    if [ ${i_mod_400} -eq 0 -o ${i_mod_4} -eq 0 -a ! ${i_mod_100} -eq 0 ]; then
        echo "1"
    else
        echo "0"
    fi
}

lb_num_leap()
{
    #
    # number of leap years comprised in <YEAR1> (included) and <YEAR2> (excluded)...
    #
    if [ "$2" = "" ]; then echo "USAGE: num_leap <YEAR1> <YEAR2>"; exit; fi
    icpt=0
    jy=${1}
    while [ ${jy} -lt ${2} ]; do
        if [ `lb_is_leap ${jy}` -eq 1 ]; then icpt=`expr ${icpt} + 1`; fi
        jy=`expr ${jy} + 1`
    done
    echo ${icpt}
}


function lb_leap_day()
{
    # We need 2 different methods to know the current date:
    # The input argument is the file date.file
    #
    if [ "$1" = "" ]; then echo "USAGE: lb_check_leap <date.file>"; exit; fi
    DATF=$1
    #
    if [ ! -f ${DATF} ]; then echo "There should be a ${DATF} !"; exit; fi
    #
    dtg1=`paste ${DATF} | cut -c 11-18`
    dtg2=`paste ${DATF} | cut -c 20-27`
    #
    yyyy=`echo ${dtg2} | cut -c1-4`
    mmdd1=`echo ${dtg1} | cut -c5-8`
    mmdd2=`echo ${dtg2} | cut -c5-8`
    #
    ileap=`expr ${yyyy} - 2000`; ileap=`expr ${ileap} % 4`
    #
    # Last DTG done:
    export CLDTG=${dtg2}
    #
    #if [ ! ${ileap} -eq 0 ]; then echo "Not a leap year !"; fi
    #
    #
    # We were at a DTBR=
    # a regular year would make a 19900225_19900229 file
    #  (should be 19900225_19900301 actually but NEMO seems to keep the same month for a file name!)
    # so when we restart a run for a leap year and the last dtg was "19900220_19900224", we force
    # DTBR=6 !
    #
    export i_leap_day=0
    #
    if [ ${ileap} -eq 0 -a "${mmdd2}" = "0224" ]; then
        #
        # Checking if the frequency was 5 days:
        if [ ! "${mmdd1}" = "0220" ]; then echo "Problem (lb_leap_day), expecting different DTG"; exit; fi
        #
        export i_leap_day=1
    fi
    #
}


function clock2hour()
{
    # convert time like hh:mm:ss
    # to decimal hour of the day
    # ex: 13:30:00  => 13.5
    #
    # Replacing ":" by " ":
    ca=`echo $1 | sed -e s/:/' '/g`
    #
    hdec=`echo ${ca} | awk '{print $1+($2*60+$3)/3600}'`
    echo ${hdec}
}



check_if_file()
{
    lspeak=true
    cmesg="$2"
    if [ "$2" = "s" ]; then
        # Silent!
        lspeak=false
        cmesg="$3"
    fi
    if [ ! -f $1 ]; then
        echo; echo "PROBLEM: file $1 is missing!!!"; echo "${cmesg}"
        exit
    else
        if ${lspeak}; then echo " ... good, file $1 is here..."; echo; fi
    fi
}




set_xtics()
{
    if [ "${1}" = "" ]; then echo "ERROR (xtics.sh): YEAR_INI is not defined!!!"; exit; fi
    if [ "${2}" = "" ]; then echo "ERROR (xtics.sh): YEAR_END is not defined!!!"; exit; fi
    #
    dy=`expr ${2} - ${1} + 1` ; export YF2=`expr ${2} + 1`
    #
    export XTICS=1
    if [ ${dy} -gt 14   -a ${dy} -le 30   ]; then export XTICS=2; fi
    if [ ${dy} -gt 30   -a ${dy} -le 120  ]; then export XTICS=5; fi
    if [ ${dy} -gt 120  -a ${dy} -le 200  ]; then export XTICS=10; fi
    if [ ${dy} -gt 200  -a ${dy} -le 400  ]; then export XTICS=20; fi
    if [ ${dy} -gt 400  -a ${dy} -le 1200 ]; then export XTICS=50; fi
    if [ ${dy} -gt 1200                   ]; then export XTICS=100; fi
}




epstopng()
{
    #
    CMD="convert -density 120x120 -quality 100" ; # Transparent background
    #CMD="convert -density 120x120 -quality 100 -background white -flatten" ; # White background :
    #
    if [ "$1" = "" ]; then
        list=`\ls *.eps`
    else
        list="$1"
    fi
    #
    if [ "`which convert`" = "" ]; then
        echo; echo "PROBLEM: convert from ImageMagick was not found in your PATH..."
        echo; exit
    fi
    #
    for feps in $list; do
        fpng=`echo ${feps} | sed -e s/'\.eps'/'\.png'/g`
        if [ ! -f ${fpng} ]; then
            echo "Creating ${fpng}"
            ${CMD} ${feps} ${fpng}
        else
            echo "${fpng} already exists!"
        fi
    done
    #
}



ipresent_var_in_ncf()
{
    ipv=0
    ca=`${NCDF_DIR}/bin/ncdump -h $1 | grep "${2}(time_counter" | grep float`
    if [ ! "${ca}" = "" ]; then
        #echo "   variable ${2} is present in file $1"
        ipv=1
    #else
    #    echo "   variable ${2} is NOT present in file $1"
    fi
    echo "${ipv}"
}



function contains_string()
{
    # Tells if string "s" (argument 1) belongs to a list (argument 2)

    nbarg="$#" ; length_list=`expr ${nbarg} - 1`
    list_all=($*)

    s=${list_all[0]}  ; # the string
    list=${list_all[@]:1:${length_list}}  ; # the list

    #echo "the string = ${s}"
    #echo "the list   = ${list_all[@]:1:${length_list}}"

    [[ ${list} =~ ${s} ]] && echo "1" || echo "0"

}


#lgzipped_file()
#{
#    lgz=false
#    nc=`echo $1 | wc -c`
#    ip=`expr ${nc} - 3`
#    fc=`echo $1 | cut -c${ip}-`
#    if [ "${fc}" = ".gz" ]; then lgz=true; fi
#    echo ${lgz}
#}

