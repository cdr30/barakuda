#!/bin/sh


list=`squeue | grep ${USER} | grep DgnOce2 | cut -c1-9`
list="${list} `squeue | grep ${USER} | grep ORCLIM | cut -c1-9`"

#list=`echo ${list} | wc -c1-9`


for rr in ${list}; do
    echo "Cancelling batch job with ID ${rr}"
    scancel ${rr}
    sleep 0.5
done