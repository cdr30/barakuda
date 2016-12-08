#!/bin/bash

fs=RbldNm.tmp

cat > ${fs} <<EOF
#!/bin/bash
#
#######
#SBATCH -A snic2014-10-3
#SBATCH --reservation=dcs
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -J RbldNm
#SBATCH -t 11:00:00
#SBATCH -o out_RbldNm.out
#SBATCH -e err_RbldNm.err
########

list=\`find . -name *_0000.nc | grep -v restart\`

here=\`pwd\`



NCCOPY=/software/apps/netcdf/4.3.2/i1403-impi-5.0.2.044-hdf5-1.8.14/bin/nccopy


cfstm=\`pwd | cut -c 1-12\`
if [ "\${cfstm}" = "/proj/climod" ]; then
    TRASH=/proj/climod/users/\${USER}/tmp/REBUILT
elif [ "\${cfstm}" = "/proj/bolinc" ]; then
    TRASH=/proj/bolinc/users/\${USER}/tmp/poubelle/REBUILT
else
    echo "Where TF are we???"; pwd ; exit
fi
mkdir -p \${TRASH}/



for ff in \${list}; do

    cd \${here}/

    dd=\`dirname \$ff\`
    fr=\`echo \$ff | sed -e "s|_0000.nc||g"\`    ; fr=\`basename \${fr}\`
    fn=\`echo \$ff | sed -e "s|_0000.nc|.nc4|g"\`; fn=\`basename \${fn}\`

    echo ; echo " *** Entering \$dd "
    cd \${dd}/

    if [ -f \${fr}.nc -o -f \${fr}.nc4 ]; then

        echo " File \${fr} already recomposed!"

    else
        echo "  => going to create \${fr}"
        nbp=\`\ls \${fr}_0*.nc | wc -l\`
        echo "rebuild_nemo -t 4 \${fr} \${nbp}"
        echo

        rebuild_nemo -t 4 \${fr} \${nbp}

        echo ; echo

        rm -f nam_rebuild

        mv -f \${fr}_0*.nc \${TRASH}/

        \${NCCOPY} -k 4 -d 9 \${fr}.nc \${fr}.nc4

        rm -f \${fr}.nc

        echo " File \${dd}/\${fr}.nc4 sucessfully created !"
        echo

    fi

done

EOF

chmod +x ${fs}

sbatch ./${fs}

