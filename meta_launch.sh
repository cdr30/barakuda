#!/bin/bash -ue
#
# K. Wyser, July 2016
#

while getopts C:R:i:e: option ; do
    case $option in
        C) conf=${OPTARG} ;;
        R) expname=${OPTARG} ;;
        i) y1=${OPTARG} ;;
        e) y2=${OPTARG} ;;
    esac
done

job_name=b-$expname

launch_cmd='sbatch -A snic2014-10-1 --reservation=dcs -n 1 -t 4:0:0'

ll=$( $launch_cmd -J $job_name-1 barakuda.sh -C $conf -R $expname )
echo $ll
job1_id=$(echo $ll | awk '{print $4}')

ll=$( $launch_cmd -J $job_name-2 build_clim.sh -C $conf -R $expname -i $y1 -e $y2)
echo $ll
job2_id=$(echo $ll | awk '{print $4}')

$launch_cmd -d afterok:$job1_id:$job2_id -J $job_name-3 barakuda.sh -C $conf -R $expname -E

exit 0
