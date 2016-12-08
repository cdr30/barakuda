#!/bin/sh

#SBATCH -A snic2014-10-3
#SBATCH --reservation=dcs
#SBATCH -w n1581
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J 'xxMKYRxx'
#SBATCH -t 23:58:00
#SBATCH -o out_MKYR_%J.out
#SBATCH -e err_MKYR_%J.err

NCDIR=/software/apps/netcdf/4.3.2/i1214-hdf5-1.8.12-AVX-off

#RUN="32-L" ; ILDIR=216
RUN="32-M" ; ILDIR=240
#RUN="32-Z" ; ILDIR=216

cext=".L75"

ROOT="/nobackup/rossby15/sm_uflad/run/${RUN}/output/nemo"


DOUT=/proj/bolinc/users/x_laubr/ORCA1${cext}/ORCA1${cext}-${RUN}-S ; mkdir -p ${DOUT}

yyyy=1990
idi=1
jm=1

dir_job_temp=/scratch/local/${SLURM_JOBID} ; mkdir -p ${dir_job_temp}

cd ${dir_job_temp}/

while [ ${idi} -le ${ILDIR} ]; do


    cdi=`printf "%03d" "${idi}"`
    cm=`printf "%02d" "${jm}"`

    echo " ${cdi} , ${yyyy} , ${cm} "
    
    for gt in "grid_T" "grid_U" "grid_V" "icemod"; do
        lf=`\ls ${ROOT}/${cdi}/${RUN}_1d_${yyyy}${cm}01_${yyyy}${cm}*_${gt}.nc`
        fin=`basename ${lf}`
        echo "rsync -avP ${ROOT}/${cdi}/${fin} ."
        rsync -avP ${ROOT}/${cdi}/${fin} .
        echo
        
        echo "ncra -O ${fin} -o m${cm}_${yyyy}_${gt}.nc"
        ncra -O ${fin} -o m${cm}_${yyyy}_${gt}.nc
        echo
        rm -f ${fin}

        echo "${NCDIR}/bin/nccopy -k 4 -d 9 m${cm}_${yyyy}_${gt}.nc m${cm}_${yyyy}_${gt}.nc4"
        ${NCDIR}/bin/nccopy -k 4 -d 9 m${cm}_${yyyy}_${gt}.nc  m${cm}_${yyyy}_${gt}.nc4
        echo
        rm -f m${cm}_${yyyy}_${gt}.nc
    done
    



    if [ $((${idi} % 12)) -eq 0 ]; then

        for gt in "grid_T" "grid_U" "grid_V" "icemod"; do
            fo=ORCA1${cext}-${RUN}_1m_${yyyy}0101_${yyyy}1231_${gt}.nc4
            echo "ncrcat -O m*_${gt}.nc4 -o ${fo}"
            ncrcat -O m*_${yyyy}_${gt}.nc4 -o ${fo}
            echo
            rsync -avP ${fo} ${DOUT}/
            rsync -avP ${fo} ${DOUT}/
            rm -f m*_${yyyy}_${gt}.nc4
            rm -f ${fo}
        done

        yyyy=$(($yyyy+1))
        jm=0

    fi

    idi=$(($idi+1))
    jm=$(($jm+1))

done

