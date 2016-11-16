#!/bin/bash

if [ "${3}" = "" ]; then
    echo "USAGE: `basename $0` <CPREF> <DIR saved NEMO files>" ; exit
fi

CPREF=${1}
DSN=${2}

HERE=`pwd`

gt=grid_T

osuff=_${gt}.nc.gz

echo; echo

cd ${DSN}/
# --time-style="+%m-%d %H:%M" should give something like : "06-21 03:19"
list=`\ls -l --time-style="+%m-%d_%H:%M" *${CPREF}*${gt}*.nc* | tr -s " " | cut -d" " -f6-9 | sed -e s/" ${CPREF}"/'.'/g | cut -c 1-16`

cd ${HERE}/

rm -f tmp_progression_${CPREF}.dat


# Getting last year done:
if [ "${DIAG_D}" = "" ]; then
    echo "ERROR: print_progression.sh does not know DIAG_D !!!"; exit
fi
fly=${DIAG_D}/last_year_done.info
if [ -f ${fly} ]; then
    jly=`cat ${fly}`
else
    echo "# year   date of completion" > tmp_progression_${CPREF}.dat
    jly=0000
fi

echo "Last year diagnosed is ${jly}"


# list of something that should lok like this: "07-20_17:23.1850"

for cl in ${list}; do
    
    ctime=`echo ${cl} | cut -d . -f1`
    model_date=`echo ${cl} | cut -d . -f2`
    echo "Model date = ${model_date} => completed ${ctime}"
    
    if [ ${model_date} -gt ${jly} ]; then
        echo "${model_date}  ${ctime}" >> tmp_progression_${CPREF}.dat
    fi
    
done

echo; echo; echo