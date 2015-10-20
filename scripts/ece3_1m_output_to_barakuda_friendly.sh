#!/bin/sh

#SBATCH -A snic2014-10-3
#SBATCH --reservation=dcs
###SBATCH -w n1581
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J 'xxMKYRxx'
#SBATCH -t 00:58:00
#SBATCH -o out_MKYR_%J.out
#SBATCH -e err_MKYR_%J.err

NCDIR=/software/apps/netcdf/4.3.2/i1214-hdf5-1.8.12-AVX-off


#RUN="32-Z" ; ILDIR=10 ; cext=".L75" ; yyyy=1990 ; idi=1 ; cimp='cp'
#RUN="32-Y" ; ILDIR=10 ; cext=".L75" ; yyyy=1990 ; idi=1 ; cimp='cp'
#RUN="32-P" ; ILDIR=5 ; cext=".L75" ; yyyy=1990 ; idi=1 ; cimp='cp'
#RUN="32-Q" ; ILDIR=9 ; cext=".L75" ; yyyy=1990 ; idi=1 ; cimp='cp'
#RUN="1cat" ; ILDIR=5 ; cext=".L75" ; yyyy=1990 ; idi=1 ; cimp='cp'
#RUN="5ca1" ; ILDIR=5 ; cext=".L75" ; yyyy=1990 ; idi=1 ; cimp='cp'
#RUN="LB00" ; ILDIR=10 ; cext=".L75" ; yyyy=1990 ; idi=1 ; cimp='mv'
RUN="T001" ; ILDIR=10 ; cext=".L75" ; yyyy=1990 ; idi=1 ; cimp='cp'
#RUN="5c13" ; ILDIR=5 ; cext=".L75" ; yyyy=1990 ; idi=1 ; cimp='cp'


#NEMO_SAVED_FILE="grid_T" "grid_U" "grid_V" "icemod"
NEMO_SAVED_FILE="grid_T icemod"


# Uwe:
ROOT="/nobackup/rossby15/sm_uflad/run/${RUN}/output/nemo"
# Klaus:
#ROOT="/nobackup/rossby15/sm_wyser/run/${RUN}/output/nemo"
#
#Me:
#ROOT="/nobackup/rossby15/x_laubr/run/${RUN}/output/nemo"


if [ ! -d ${ROOT} ]; then
    echo "Mhhhh! Directory ${ROOT} does not exist !!!"
fi


DOUT=/proj/bolinc/users/x_laubr/ORCA1${cext}/ORCA1${cext}-${RUN}-S ; mkdir -p ${DOUT}




dir_job_temp=/scratch/local/${SLURM_JOBID} ; mkdir -p ${dir_job_temp}

cd ${dir_job_temp}/

while [ ${idi} -le ${ILDIR} ]; do
    
    cdi=`printf "%03d" "${idi}"`

    echo " ${cdi} , ${yyyy} "
    
    for gt in ${NEMO_SAVED_FILE}; do

        TT="${yyyy}0101_${yyyy}1231"
        
        fin=${RUN}_1m_${TT}_${gt}.nc


        if [ "${cimp}" = "cp" ]; then
            
            fo=ORCA1${cext}-${RUN}_1m_${TT}_${gt}.nc4
            echo "${NCDIR}/bin/nccopy -k 4 -d 9 ${ROOT}/${cdi}/${fin}  ${DOUT}/${fo}"
            ${NCDIR}/bin/nccopy -k 4 -d 9 ${ROOT}/${cdi}/${fin}  ${DOUT}/${fo}
            
        elif [ "${cimp}" = "mv" ]; then
            
            fo=ORCA1${cext}-${RUN}_1m_${TT}_${gt}.nc
            echo "mv -f  ${ROOT}/${cdi}/${fin}  ${DOUT}/${fo}"
            mv -f  ${ROOT}/${cdi}/${fin}  ${DOUT}/${fo}

        fi
        echo
        
    done
    
    echo
    echo "Done for year $yyyy !"
    echo

    idi=$(($idi+1))
    yyyy=$(($yyyy+1))
    
done

