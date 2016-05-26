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

export BARAKUDA_ROOT=`pwd`

CANOPY_PATH=${HOME}/opt/Canopy_64bit/User
PYTH="${CANOPY_PATH}/bin/python -W ignore" ; # which Python installation to use
export PYTHONPATH=${CANOPY_PATH}/lib/python2.7/site-packages:${BARAKUDA_ROOT}/python/modules ; # PATH to python barakuda modules
PYBRKD_EXEC_PATH=${BARAKUDA_ROOT}/python/exec         ; # PATH to python barakuda executable

#export FIG_FORMAT='svg'
export FIG_FORMAT='png'

# Supported ORCA grids:
ORCA_LIST="ORCA1.L75 ORCA1.L46 ORCA1.L42 ORCA2 ORCA2_L46"

# Checking available configs
list_conf=`\ls configs/config_*.sh` ; list_conf=`echo ${list_conf} | sed -e s/'configs\/config_'/''/g -e s/'.sh'/''/g`

# Important bash functions:
. ${BARAKUDA_ROOT}/configs/bash_functions.bash

usage()
{
    echo
    echo "USAGE: ${0} -C <config> -R <run1,run2,...,runN>  (options)"
    echo
    echo "     Available configs are:"
    echo "             => ${list_conf}"
    echo
    echo "   OPTIONS:"
    echo "      -y <YYYY> => force initial year to YYYY"
    echo
#    echo "      -c <run>  => 2D comparison diagnostics are performed against run <run>"
#    echo "                   instead of a climatology"
    echo
    echo "      -f        => forces creation of diagnostic page eventhough treatment of output files is not finished"
    echo
    echo "      -e        => create the HTML diagnostics page on local or remote server"
    echo
    echo "      -h        => print this message"
    echo
    exit
}


YEAR0="" ; iforcey0=0

IPREPHTML=0
IFORCENEW=0

while getopts C:R:y:feh option ; do
    case $option in
        C) CONFIG=${OPTARG};;
        R) CRUNS=${OPTARG};;
        y) YEAR0=${OPTARG} ; iforcey0=1 ;;
        f) IFORCENEW=1;;
        e) IPREPHTML=1;;
        h)  usage;;
        \?) usage ;; 
    esac
done





if [ "${CONFIG}" = "" -o "${CRUNS}" = "" ]; then usage ; exit ; fi

for og in ${ORCA_LIST}; do
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














LRUNS=`echo ${CRUNS} | sed -e s/'\,'/'\ '/g -e s/'\, '/'\ '/g`
echo; echo "Runs to be treated: ${LRUNS}"
nbr=`echo ${LRUNS} | wc -w` ; echo "  => number of runs to compare = ${nbr}"
echo " NEMO grid = ${ORCA}";
echo " reading config into: ${fconfig}"
echo; echo




export CONF=${CONF}
export LIST_RUNS=${LRUNS}





NRUNS="${ORCA}-`echo ${LRUNS} | sed -e 's/\ /_/g'`"
echo " Label to be used: ${NRUNS}" ; echo

BASE_NAME="comp_${NRUNS}"

DIAG_COMP_DIR=${DIAG_DIR}/comparisons/${BASE_NAME} ; rm -rf ${DIAG_COMP_DIR} ; mkdir -p ${DIAG_COMP_DIR}




YEAR_INI=4000
YEAR_END=0

# just that they become righ arrays...
VRUNS=( ${LRUNS} ) ;  VCONFRUNS=( ${LRUNS} ) ; VDIAGS=( ${LRUNS} ) ; 

jr=0

for run in ${LRUNS}; do
    
    echo; echo " RUN ${run} "
    RUN="${run}"
    
    echo "   => ORCA = ${ORCA}"
    
    CONFRUN=${ORCA}-${RUN}
    echo "   => CONFRUN = ${CONFRUN}"
    VCONFRUNS[${jr}]=${CONFRUN}
    

    DIAG_D="${DIAG_DIR}/${CONFRUN}"
    VDIAGS[${jr}]=${DIAG_D}
    echo "   => DIAG_D = ${DIAG_D} "; echo ; echo

    
    if [ ! -d ${DIAG_D} ]; then
        echo "PROBLEM: ${DIAG_D} does not exist!"
        echo "    =>  you must run barakuda for ${run} prior to comparison!"
        exit
    fi


    # Guessing initial and last year:
    
    check_if_file ${DIAG_D}/first_year.info
    iy=`cat ${DIAG_D}/first_year.info`
    if [ ${iy} -lt ${YEAR_INI} ]; then export YEAR_INI=${iy}; fi

    check_if_file ${DIAG_D}/last_year_done.info
    iy=`cat ${DIAG_D}/last_year_done.info`
    if [ ${iy} -gt ${YEAR_END} ]; then export YEAR_END=${iy}; fi
    
    jr=`expr ${jr} + 1`
done


echo
echo " Global YEAR_INI = ${YEAR_INI}"
echo " Gloabl YEAR_END = ${YEAR_END}"
echo



cd ${DIAG_COMP_DIR}/

${PYTH} ${PYBRKD_EXEC_PATH}/compare_time_series.py ${YEAR_INI} ${YEAR_END}



# Configuring HTML display file:
sed -e "s|{TITLE}|Ocean, ${JTITLE}: comparing ${VRUNS[*]}|g" \
    -e "s|{NRUNS}|${NRUNS}|g" -e "s|{DATE}|`date`|g" -e "s|{HOST}|`hostname`|g" \
    ${BARAKUDA_ROOT}/scripts/html/index_comp_skel.html > index.php


list_figs=`\ls *.${FIG_FORMAT}`

for ff in ${list_figs}; do
    echo "<br><br><big> `echo ${ff} | sed -e s/.${FIG_FORMAT}//g -e s/_comparison//g` </big><br>" >> index.php
    echo "  <img style=\"border: 0px solid\" alt=\"\" src=\"${ff}\"> <br><br><br>" >> index.php

done

cat ${BARAKUDA_ROOT}/scripts/html/index_skel_footer.html >> index.php ; # Closing HTML file...

cp ${BARAKUDA_ROOT}/scripts/html/conf_*.html .
cp ${BARAKUDA_ROOT}/scripts/html/logo.png    .
echo; echo


if [ ${ihttp} -eq 0 ]; then
    
    echo "Diagnostic page installed in `pwd`"
    echo " => view this directory with a web browser (index.php)..."
    
else
    
    RWWWD=${RWWWD}/COMPARISONS_time_series
    ssh ${RUSER}@${RHOST} "mkdir -p ${RWWWD}"
    
    if [ ${ihttp} -eq 1 ]; then
        
        echo "Preparing to export to remote host!"; echo
        
        cd ../

        tar cvf ${BASE_NAME}.tar ${BASE_NAME}
        scp ${BASE_NAME}.tar ${RUSER}@${RHOST}:${RWWWD}/
        ssh ${RUSER}@${RHOST} "cd ${RWWWD}/; rm -rf ${BASE_NAME}; tar xf ${BASE_NAME}.tar 2>/dev/null; rm -f ${BASE_NAME}.tar; chmod -R a+r ${BASE_NAME}"
        rm -f ${BASE_NAME}.tar
        
        echo; echo
        echo "Diagnostic page installed on remote host ${RHOST} in ${RWWWD}/${BASE_NAME}!"
        echo "( Also browsable on local host in `pwd` )"
        
        else

        echo "Error: \"ihttp\" is either 0 or 1 !"
        
    fi

fi

rm -rf ${COMP_DIR}

echo; echo
