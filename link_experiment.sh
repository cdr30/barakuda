#!/bin/sh
# Link experiment files for barakuda postprocessing
# Run in the directory where files are linked to

EPATH=/lustre/tmp/uotilap/ECN-D501
GRID=ORCA1.L46

for GTYPE in grid_U grid_V grid_T icemod
do
    for f in ${EPATH}/output/nemo/???/*_${GTYPE}.nc
        do 
            b=`basename $f`
            ln -s $f ${GRID}-$b
        done
done
