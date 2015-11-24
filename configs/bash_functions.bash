#!/bin/sh

function lb_wait()
{
    if [ "$2" = "" ]; then
        echo "USAGE: lb_wait <number of nodes used> <run name (CASE)>"
        exit
    fi
    # $1 is the number of nodes used, useful to know which JOB...
    echo
    clrun=`spq -w | grep ${USER} | grep ${1}K | grep ${2}`; sleep 4
    clrun2=`spq -w | grep ${USER} | grep ${1}K | grep ${2}`
    #
    if [ "${clrun}" = "" -a "${clrun2}" = "" ]; then
        echo "No simulations are running!"
    else
        while [ "${clrun}" != "" ]; do
            echo "At least a simulation is running! Waiting...   `date`"
            sleep 300
            clrun=`spq -w | grep ${USER} | grep ${1}K | grep ${2}`
        done
    fi
    echo; echo
}

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
    ca=`${NCDF_BIN}/ncdump -h $1 | grep "${2}(time_counter" | grep float`
    if [ ! "${ca}" = "" ]; then
        #echo "   variable ${2} is present in file $1"
        ipv=1
    #else
    #    echo "   variable ${2} is NOT present in file $1"
    fi
    echo "${ipv}"
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

