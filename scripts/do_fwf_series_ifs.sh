#!/bin/bash

# Lagerloef and colleagues note that:
#
# Evaporation from the global ocean is estimated to be ~ 13 Sv, and the
# precipitation sums to ~ 12.2 Sv. The difference of 0.8 Sv compares with the
# estimate of river input of 1.2 Sv. The apparent imbalance of ~ 0.4 Sv excess
# is smaller than the estimated error bars. Sea level rise due to melting
# glaciers is only ~ 0.01 Sv, so cannot account for the imbalance. Global
# groundwater flows are poorly known but generally estimated to be similarly
# small (Cable et al., 1996). Likely sources of error include scant data in the
# southern oceans and a possible underestimate of evaporation in very-high-wind
# conditions. Other surface flux climatologies display similar patterns, but the
# range of the estimates is quite large and unlikely to be significantly
# improved in the near future.
#
# References:
#    Lagerloef, G., Schmitt, R., Schanze, J. Kao, H-Y, 2010, The Ocean and the Global Water Cycle, Oceanography, Vol.23, No.4
#



#RUN=32bG
#Y_INI_EC=1990

#cyear=2038
#NEMO_OUT_D=/nobackup/rossby15/sm_uflad/run/${RUN}/output/nemo
#DIAG_D=/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${RUN}


cmsg="ERROR: $0 => global variable"
if [ -z ${RUN} ]; then echo "${cmsg} RUN is unknown!"; exit ; fi
if [ -z ${Y_INI_EC} ]; then echo "${cmsg} Y_INI_EC is unknown!"; exit ; fi
if [ -z ${cyear} ]; then echo "${cmsg} cyear is unknown!"; exit ; fi
if [ -z ${NEMO_OUT_D} ]; then echo "${cmsg} NEMO_OUT_D is unknown!"; exit ; fi
if [ -z ${DIAG_D} ]; then echo "${cmsg} DIAG_D is unknown!"; exit ; fi


RUN_DIR=`echo ${NEMO_OUT_D} | sed -e "s|/output/nemo||g"`
IFS_OUT_D=`echo ${NEMO_OUT_D} | sed -e "s|/output/nemo|/output/ifs|g"`

echo
echo " RUN_DIR = ${RUN_DIR}"
echo " IFS_OUT_D = ${IFS_OUT_D}"
echo



YDIR=$((${cyear}-${Y_INI_EC}+1))

dir_ece=`printf "%03d" ${YDIR}`
dir_ece=${RUN_DIR}/output/ifs/${dir_ece}

echo ${dir_ece}


F_AREA=${RUN_DIR}/areas.nc
F_MASK=${RUN_DIR}/masks.nc

echo ${F_AREA} ; ls -l ${F_AREA}
echo ${F_MASK} ; ls -l ${F_MASK}


#rm -f *.grb *.nc *.tmp


# Create ifs_area_masked:
echo
CMD="cdo setmisstoc,0 -ifthen -eqc,0 -selvar,A128.msk ${F_MASK} -selvar,A128.srf ${F_AREA} ifs_area_masked.nc"
echo $CMD
${CMD}
echo
ncrename -O -v A128.srf,ifs_area_masked ifs_area_masked.nc -o ifs_area_masked.nc
ncwa -O -a y ifs_area_masked.nc ifs_area_masked.nc # remove y of length !
# Add degenerate time record to masks and appends it to ftreat:
ncecat -O ifs_area_masked.nc -o ifs_area_masked.nc
ncrename -d record,time ifs_area_masked.nc




for VAR in  "e" "lsp" "cp"; do
    
    for cm in "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12"; do

        ftreat=${VAR}_${RUN}_${cyear}${cm}

        grib_copy -w shortName=${VAR} ${dir_ece}/ICMGG${RUN}+${cyear}${cm} ${ftreat}.grb


        BVAR=`echo ${VAR}| tr '[:lower:]' '[:upper:]'`


#  post-processing timestep in seconds
        if [ "${cm}" = "01" ]; then
            pptime=$(cdo showtime -seltimestep,1,2 ${ftreat}.grb | tr -s ' ' ':' | awk -F: '{print ($5-$2)*3600+($6-$3)*60+($7-$4)}' )
            if [ $pptime -le 0 ]; then
                pptime=21600 # default 6-hr output timestep
            fi
            echo " pptime = ${pptime} seconds !"
        fi



    # To netcdf but as a vector not 2D array as when using "-R":
        cdo -t ecmwf -f nc copy ${ftreat}.grb ${ftreat}.nc  ; rm -f ${ftreat}.grb
        ncrename -h -O -v ${BVAR},${VAR} ${ftreat}.nc
        ncrename -h -O -d rgrid,x ${ftreat}.nc

    # To m/s
        ncap2 -h -A -s "${VAR}=${VAR}/${pptime}" ${ftreat}.nc -o ${ftreat}.nc
        ncatted -O -a units,${VAR},o,c,'m of water / s' ${ftreat}.nc

        echo

    # To netcdf monthly:
        ncra -h -O ${ftreat}.nc -O ${ftreat}_m.nc
        mv -f ${ftreat}_m.nc ${ftreat}.nc








#    cdo -t ecmwf -f nc shifttime,-1month -timmean -selvar,${BVAR} ${ftreat}.grb ${ftreat}.nc
#    ncrename -h -O -d rgrid,x ${ftreat}.nc
#ncwa -h -O -a time ${ftreat}.nc -o ${ftreat}.nc
#ncks -h -O -C -x -v time ${ftreat}.nc -o ${ftreat}.nc

# NUmber of time steps for current month:
#    nbrec=`grib_ls ${ftreat}.grb | tail -1 | cut -d ' ' -f1`
#    echo "nbrec = ${nbrec}"

#    rm -f ${ftreat}.grb
#
#    tot_time=`expr ${nbrec} \* ${pptime}`
#
#    echo "Total time in seconds for this month: ${tot_time}" ; echo




# Ok to kg/m2/s == mm/s:
#ncap2 -h -A -s "${VAR}=${BVAR}*1000./${pptime}" ${ftreat}.nc -o ${ftreat}.nc
#ncatted -O -a units,${VAR},o,c,'kg/m^2/s' ${ftreat}.nc
#ncks -h -O -C -x -v ${BVAR} ${ftreat}.nc -o ${ftreat}.nc
#OR:
#    ncrename -h -O -v ${BVAR},${VAR} ${ftreat}.nc


# Convert from a height of water to a flux of water (dependant on the output freq.):
# precip and evap in kg m-2 s-1


#cdo -t $ecearth_table setparam,228.128 -expr,"totp=1000*(lsp+cp)/$pptime" \
#    icmgg_${year}.grb tmp_totp.grb
#cdo -R -t $ecearth_table setreftime,$reftime tmp_totp.grb ${out}_totp.nc
#cdo -R -t $ecearth_table splitvar -setreftime,$reftime \
#    -mulc,1000 -divc,$pptime -selvar,e,lsp,cp icmgg_${year}.grb ${out}_
#cdo -R -t $ecearth_table setparam,80.128 -fldmean \
#    -expr,"totp=1000*(lsp+cp+e)/$pptime" icmgg_${year}.grb tmp_pme.grb
#cdo -t $ecearth_table setreftime,$reftime tmp_pme.grb ${out}_pme.nc


# divide fluxes by PP timestep
#cdo -R -t $ecearth_table splitvar -setreftime,$reftime -divc,$pptime \
#    -selvar,ssr,str,sshf,slhf,tsr,ttr,ewss,nsss,ssrc,strc,tsrc,ttrc,ssrd,strd \
#    icmgg_${year}.grb ${out}_





        echo
        ncks -h -A -v ifs_area_masked ifs_area_masked.nc -o ${ftreat}.nc
        echo


        isign=1
        if [ "${VAR}" = "e" ]; then isign=-1; fi

    # Multiplying ${VAR} and ifs_area_masked:
        ncap2 -h -A -s "${VAR}2d=(${isign}*ifs_area_masked*${VAR})" ${ftreat}.nc -o ${ftreat}.nc


    # Total volume evaporated over ocean during the current month:
        ncap2 -h -A -s "flx_${VAR}_sv=${VAR}2d.total(\$x)*1.E-6" ${ftreat}.nc -o ${ftreat}.nc
        ncatted -O -a units,flx_${VAR}_sv,o,c,'Sv' ${ftreat}.nc

        ncks -h -A -v flx_${VAR}_sv ${ftreat}.nc -o final_${cm}.nc

        rm -f ${ftreat}.nc

#cdo output –div –fldsum –mul ifs_area_masked E –fldsum ifs_area_masked ${ftreat}.nc tamere.nc

    done
    
done

ncrcat -O final_*.nc -o final.nc

ncap2 -h -O -s "time=array(${cyear}.0416667,0.08333333,\$time)" final.nc -o final.nc
ncatted -O -a units,time,o,c,'years' final.nc

ncap2 -h -A -s "flx_p_sv=flx_cp_sv+flx_lsp_sv" final.nc
ncap2 -h -A -s "flx_emp_sv=flx_e_sv-flx_p_sv"  final.nc

rm -f ifs_area_masked.nc final_*.nc

fout=${DIAG_D}/mean_fwf_IFS_32bI_global.nc

if [ ! -f ${fout} ]; then
    mv final.nc ${fout}
else
    ncrcat -h -A ${fout} final.nc -o ${fout}
fi
rm -f final.nc


