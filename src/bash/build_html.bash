#!/usr/bin/env bash

function build_index_html()
{
    echo; echo; echo
    echo "Creating HTML file!"

    ctl='<br><br><br><big><big>' ; ctr='</big></big><br><br>'
    spf='<br><br>'
    img_l='<img style="border: 0px solid" alt="" src="' ; img_r='"> <br><br>'

    cr="${CONFRUN}" ; ff="${FIG_FORM}"

    cd ${DIAG_D}/

    rm -f index.html

    if [ "${JTITLE}" = "" ]; then
        echo "Problem, variable JTITLE is not set!" ; exit
    else
        TITLE="Ocean diagnostics, run ${RUN}, conf ${ORCA}_${JTITLE}"
    fi

    # Starting to configure HTML index file:
    sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${CONFRUN}|g" -e "s|{DATE}|`date`|g" -e "s|{HOST}|${HOST}:`hostname`|g" \
        ${BARAKUDA_ROOT}/src/html/conf_start.html > index.html

    # Climato section
    if ${l_pclim}; then
        cat >> index.html <<EOF
    ${ctl} Diags from climatology (${CLIM_PER}) ${ctr}
    <big> <a href="./temp_sal/index.html"> Temperature and Salinity vs CLIM</a> </big>             ${spf}
    <big> <a href="./ssh/index.html">  Sea Surface Height </a> </big>                              ${spf}
    <big> <a href="./sea_ice/index.html">  Arctic and Antarctic sea-ice extent vs CLIM </a> </big> ${spf}
    <big> <a href="./mld/index.html">  Mixed Layer depth in relevent regions </a> </big>           ${spf}
    <big> <a href="./moc/index.html">  Meridional Overturning Circulation </a> </big>              ${spf}
EOF
        if ${lcomp_to_run}; then
            cat >> index.html <<EOF
    ${ctl} Comparison with run ${RUNREF}, climatology (2004-2007) ${ctr}
    <big> <a href="./temp_sal/index_${RUNREF}.html"> Temperature and Salinity vs ${RUNREF}</a> </big>             ${spf}
    <!--        <big> <a href="./ssh/index_${RUNREF}.html">  Sea Surface Height </a> </big>                       ${spf}
    <big> <a href="./sea_ice/index_${RUNREF}.html">  Arctic and Antarctic sea-ice extent vs ${RUNREF} </a> </big> ${spf}
    <big> <a href="./mld/index_${RUNREF}.html">  Mixed Layer depth in relevent regions </a> </big>                ${spf}
    <big> <a href="./moc/index_${RUNREF}.html">  Meridional Overturning Circulation </a> </big>                   ${spf}
    -->
EOF
        fi
    fi

    if [ ${i_do_movi} -eq 1 ]; then
        cat >> index.html <<EOF
    ${ctl} Evolution of SST and SSS biases (w.r.t observations) ${ctr}
    ${img_l} dsst_${cr}.gif ${img_r}
    ${img_l} dsss_${cr}.gif ${img_r}
EOF
    fi

    # AMOC section:
    cat >> index.html <<EOF
    ${ctl} Atlantic Meridional Overturning Circulation ${ctr}
    ${img_l} amoc_${cr}.${ff} ${img_r}
    ${img_l} amoc_${cr}_comp.${ff} ${img_r}
EOF

    # Temperature section
    cat >> index.html <<EOF
    ${ctl} Temperature time-series ${ctr}
    ${img_l} 3d_thetao_${cr}.${ff} ${img_r}
    ${img_l} mean_tos_${cr}.${ff} ${img_r}
    ${img_l} 3d_thetao_lev_${cr}.${ff} ${img_r}
    ${img_l} 3d_thetao_basins_${cr}.${ff} ${img_r}
    ${img_l} Nino34_${cr}.${ff} ${img_r}
    ${img_l} hov_temperature_${cr}_global.${ff} ${img_r}
    ${img_l} hov_temperature_${cr}_atlantic.${ff} ${img_r}
    ${img_l} hov_temperature_${cr}_pacific.${ff} ${img_r}
    ${img_l} hov_temperature_${cr}_indian.${ff} ${img_r}
EOF

    # Salinity section
    cat >> index.html <<EOF
    ${ctl} Salinity time-series ${ctr}
    ${img_l} 3d_so_${cr}.${ff} ${img_r}
    ${img_l} mean_sos_${cr}.${ff} ${img_r}
    ${img_l} 3d_so_lev_${cr}.${ff} ${img_r}
    ${img_l} 3d_so_basins_${cr}.${ff} ${img_r}
    ${img_l} hov_salinity_${cr}_global.${ff} ${img_r}
    ${img_l} hov_salinity_${cr}_atlantic.${ff} ${img_r}
    ${img_l} hov_salinity_${cr}_pacific.${ff} ${img_r}
    ${img_l} hov_salinity_${cr}_indian.${ff} ${img_r}
EOF


LIST_FW_FIG="zos fwf_fwf fwf_emp fwf_prc fwf_rnf fwf_clv \
fwf_prc_NEMO_IFS_annual fwf_emp_IFS fwf_emp_IFS_annual fwf_evp_IFS \
fwf_evp_IFS_annual fwf_rnf_IFS fwf_rnf_IFS_annual fwf_prc_IFS fwf_emp_ALL_IFS"

    cat >> index.html <<EOF
    ${ctl} Freshwater-flux-related time-series ${ctr}
EOF

    for fd in ${LIST_FW_FIG}; do
        fgn="mean_${fd}_${cr}.${ff}"; fgf="${HTML_DIR}/${fgn}"
        if [ -f ${fgf} ]; then
            echo "${img_l} ${fgn} ${img_r}" >> index.html
    #else
    #    echo "LOLO: ${fgf} not found!!!!"
        fi
    done


    # Sea-ice section
    if [ ${i_do_ice}  -gt 0 ]; then
        if [ ${i_do_movi} -eq 1 ]; then
            cat >> index.html <<EOF
    ${ctl} Evolution of Arctic/Antarctic concentration ${ctr}
    ${img_l} icen_${cr}.gif ${img_r}
    ${img_l} ices_${cr}.gif ${img_r}
EOF
        fi
        cat >> index.html <<EOF
    ${ctl} Arctic/Antarctic sea-ice time-series${ctr}
    ${img_l} seaice_extent_winter_${cr}.${ff} ${img_r}
    ${img_l} seaice_extent_summer_${cr}.${ff} ${img_r}
    ${img_l} seaice_volume_winter_${cr}.${ff} ${img_r}
    ${img_l} seaice_volume_summer_${cr}.${ff} ${img_r}
EOF
    fi

    if [ ${i_do_trsp} -gt 0 ]; then
        # Adding transport section part:
        echo "    ${ctl} Transport through sections${ctr}" >> index.html
        list_section=`cat ${TRANSPORT_SECTION_FILE} | grep '-'`
        for cs in ${list_section}; do
            echo ${cs}
            echo "    ${img_l} transport_vol_${cs}_${cr}.${ff} ${img_r}"  >> index.html
            echo "    ${img_l} transport_heat_${cs}_${cr}.${ff} ${img_r}" >> index.html
            echo "    <br>" >> index.html
        done
    fi

    # Checking if figures with time-series of MLD in specified boxes are here and adding them:
    if [ ${i_do_mean} -eq 1 ]; then
        list_mld_figs=`\ls ${HTML_DIR}/mean_mldr10_1_${cr}*.${ff}`
        if [ ! "${list_mld_figs}" = "" ]; then
            echo "    ${ctl} Horizontally-averaged Mixed-Layer Depth in different regions${ctr}" >> index.html
            for fmld in ${list_mld_figs}; do
                echo "    ${img_l} ${fmld} ${img_r}"  >> index.html
            done
        fi
    fi

    if [ ${i_do_sigt} -eq 1 ]; then
        # Adding transport by sigma class section part:
        echo "${ctl} Transport by sigma class at Nordic sills${ctr}" >> index.html
        list_section=`cat ${DENSITY_SECTION_FILE} | grep '_'`
        for cs in ${list_section}; do
            echo ${cs}
            echo "    ${img_l} transport_sigma_class_${cs}_${cr}.${ff} ${img_r}"  >> index.html
        done
        echo "    ${img_l} tr_sigma_gt278_${cr}.${ff} ${img_r}"  >> index.html
        echo "    ${spf}" >> index.html
    fi

    if [ ${i_do_mht} -eq 1 ]; then
        # Adding meridional heat transport:
        echo "${ctl} Meridional transports${ctr}"  >> index.html
        for coce in "global" "atlantic" "pacific" "indian"; do
            echo "    ${img_l} MHT_${cr}_${coce}.${ff} ${img_r}"     >> index.html
            echo "    ${img_l} MST_${cr}_${coce}.${ff} ${img_r}" >> index.html
        done
        echo "    ${spf}" >> index.html
    fi

    cat ${BARAKUDA_ROOT}/src/html/conf_end.html >> index.html

}




function build_sub_html()
{
    echo; echo; echo
    echo "Creating sub HTML files!"

    ctl='<br><br><br><big><big>' ; ctr='</big></big><br><br>'
    spf='<br><br>'
    img_l='<img style="border: 0px solid" alt="" src="' ; img_r='"> <br><br>'

    cr="${CONFRUN}" ; ff="${FIG_FORM}"

    cd ${DIAG_D}/

       # T, S, SSH and ice HTML page:
    for cdiag in ${DIRS_2_EXP}; do
        cat ${BARAKUDA_ROOT}/src/html/conf_start.html               > index.tmp
        cat ${BARAKUDA_ROOT}/src/html/${cdiag}/index_${cdiag}.html >> index.tmp
        cat ${BARAKUDA_ROOT}/src/html/conf_end.html                >> index.tmp
        sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${cr}|g" \
            -e "s|{DATE}|`date`|g" -e "s|{HOST}|`hostname`|g" -e "s|{COMP2D}|CLIM|g" \
            index.tmp > ${cdiag}/index.html
        rm -f index.tmp
        cd ${cdiag}/ ; ln -sf ../logo.png . ; cd ../
    done

    for var in "sst" "sss" "sections_ts" "ts_100m" "ts_1000m" "ts_3000m"; do
        cat ${BARAKUDA_ROOT}/src/html/conf_start.html               > index.tmp
        cat ${BARAKUDA_ROOT}/src/html/temp_sal/${var}.html         >> index.tmp
        cat ${BARAKUDA_ROOT}/src/html/conf_end.html                >> index.tmp
        sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${cr}|g" \
            -e "s|{DATE}|`date`|g" -e "s|{HOST}|`hostname`|g" -e "s|{COMP2D}|CLIM|g" \
            index.tmp > temp_sal/${var}_CLIM.html
    done

    if ${lcomp_to_run}; then
        for cdiag in ${DIRS_2_EXP_RREF}; do
            cat ${BARAKUDA_ROOT}/src/html/conf_start.html               > index.tmp
            cat ${BARAKUDA_ROOT}/src/html/${cdiag}/index_${cdiag}.html >> index.tmp
            cat ${BARAKUDA_ROOT}/src/html/conf_end.html                >> index.tmp
            sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${cr}|g" \
                -e "s|{DATE}|`date`|g" -e "s|{HOST}|`hostname`|g" -e "s|{COMP2D}|${RUNREF}|g" \
                index.tmp > ${cdiag}/index_${RUNREF}.html
            rm -f index.tmp
            cd ${cdiag}/ ; ln -sf ../logo.png . ; cd ../
        done
        for var in "sst" "sss" "sections_ts" "ts_100m" "ts_1000m" "ts_3000m"; do
            cat ${BARAKUDA_ROOT}/src/html/conf_start.html               > index.tmp
            cat ${BARAKUDA_ROOT}/src/html/temp_sal/${var}.html         >> index.tmp
            cat ${BARAKUDA_ROOT}/src/html/conf_end.html                >> index.tmp
            sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${cr}|g" \
                -e "s|{DATE}|`date`|g" -e "s|{HOST}|`hostname`|g" -e "s|{COMP2D}|${RUNREF}|g" \
                index.tmp > temp_sal/${var}_${RUNREF}.html
        done
    fi

}
