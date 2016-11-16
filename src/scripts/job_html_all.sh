#!/bin/sh

#runs="SLB3501 SLB3502 SLB3503 SLB3504 SLB3505 SLB3506 SLBCOARE"
#CONF=ORCA1_v35
#for rr in ${runs}; do
#    ./barakuda.sh -R ${rr} -C ${CONF} -e -c SLB3500
#done

runs="SFRO251 SFRO252 SFRO253 SFRO254"
CONF=ORCA2_L46_v35
for rr in ${runs}; do
    ./barakuda.sh -R ${rr} -C ${CONF} -e
done


