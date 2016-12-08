#!/usr/bin/env bash

#==============================================================
#
#                    B A R A K U D A
#
#    An OCEAN MONITORING python environment for NEMO
#
#             L. Brodeau, 2009-2016
#
#===============================================================

export BARAKUDA_ROOT=`pwd`

# Checking available configs
list_conf=`\ls configs/config_*.sh` ; list_conf=`echo ${list_conf} | sed -e s/'configs\/config_'/''/g -e s/'.sh'/''/g`

# Important bash functions:
. ${BARAKUDA_ROOT}/src/bash/bash_functions.bash

barakuda_init

while getopts C:R:f:y:c:FeEh option ; do
    case $option in
        C) export CONFIG=${OPTARG} ;;
        R) export RUN=${OPTARG} ;;
        f) export IFREQ_SAV_YEARS=${OPTARG} ;;
        y) export YEAR0=${OPTARG} ; export LFORCE_YINI=true ;;
        c) export RUNREF=${OPTARG} ;;
        F) export LFORCEDIAG=true ;;
        e) export ISTAGE=2 ;;
        E) export ISTAGE=2 ; export l_clim_diag=true ;;
        h)  barakuda_usage ;;
        \?) barakuda_usage ;;
    esac
done

barakuda_check

# sourcing configuration file
fconfig=${BARAKUDA_ROOT}/configs/config_${CONFIG}.sh

if [ -f ${fconfig} ]; then
    echo "Sourcing ${fconfig} !"
    . ${fconfig}
else
    echo "PROBLEM: cannot find file ${fconfig} !"; exit
fi

# List of CDFTOOLS executables needed for the diagnostics:
export L_EXEC="cdfmaxmoc.x cdfmoc.x cdfvT.x cdftransportiz.x cdficediags.x cdfmhst.x cdfsigtrp.x"

barakuda_setup

echo
echo " SETTINGS: "
echo "   *** CLIM_DIR = ${CLIM_DIR} "
echo "   *** TMP_DIR  = ${TMP_DIR} "
echo "   *** CONFIG   = ${CONFIG} "
echo "   *** GRID     = ${ORCA} "
echo "   *** RUN      = ${RUN} "
echo "   *** CPREF    = ${CPREF} "
echo "   *** IFREQ_SAV_YEARS = ${IFREQ_SAV_YEARS} "
echo


if [ ${ISTAGE} -eq 1 ]; then
    barakuda_first_last_years ; # look at NEMO files to know what are first and last years available...
    echo ${IFREQ_SAV_YEARS} > ${DIAG_D}/numb_year_per_file.info
    echo ${YEAR_INI}        > ${DIAG_D}/first_year.info
else
    # -> this is stage 2 (plot generation) ISTAGE=2 !
    barakuda_init_plot
fi

cyear_ini=`printf "%04d" ${YEAR_INI}`
cyear_end=`printf "%04d" ${YEAR_END}`

# setup over...
###########################################################################################



jyear=${YEAR_INI}

fcompletion=${DIAG_D}/last_year_done.info
if [ -f ${fcompletion} ]; then jyear=`cat ${fcompletion}`; ((jyear++)); fi

cd ${TMP_DIR}/

barakuda_import_mesh_mask ; # Importing mesh_mask (+basin) files...

if [ ${ISTAGE} -eq 1 ]; then
    # Importing cdftools executables:
    for ex in ${L_EXEC}; do rsync -avP ${BARAKUDA_ROOT}/cdftools_light/bin/${ex} . ; done
fi

sgz=""



# L O O P   A L O N G   Y E A R S
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lcontinue=true

if ${LFORCEDIAG}; then lcontinue=false; fi

while ${lcontinue}; do

    export cyear=`printf "%04d" ${jyear}`

    cpf=""
    if [ ${ISTAGE} -eq 1 ] && [ ${ece_run} -gt 0 ]; then
        iy=$((${jyear}-${YEAR_INI}+1+${YEAR_INI}-${YEAR_INI_F}))
        dir_ece=`printf "%03d" ${iy}`
        echo " *** ${cyear} => dir_ece = ${dir_ece}"
        cpf="${dir_ece}/"
    fi

    TTAG_ann=${cyear}0101_${cyear}1231

    i_get_file=0
    if [ $((${jyear}%${IFREQ_SAV_YEARS})) -eq 0 ]; then
        barakuda_check_year_is_complete  ; # lcontinue might be updated to false!
    fi

    CRTM=${CPREF}${TTAG}
    CRT1=${CPREF}${TTAG_ann}

    if ${lcontinue}; then

        if [ ${ISTAGE} -eq 2 ]; then
            echo; echo "You cannot create figures and HTML pages yet!"
            echo " => finish treating the results first by launching barakuda.sh without the '-e' switch."
            exit
        fi

        echo; echo; echo; echo
        echo "*********************************************************************"
        echo "  Run ${RUN}: Will generate diagnostics and data for year ${cyear}..."
        echo "*********************************************************************"
        echo ; echo

        barakuda_import_files

        # Testing if ALL required files are present now:
        for gt in ${NEMO_SAVED_FILES}; do
            ftt="./${CRT1}_${gt}.nc" ;  check_if_file ${ftt}
        done
        echo; echo "All required files are in `pwd` for year ${cyear} !"; echo

        # Files to work with for current year:
        ft=${CRT1}_grid_T.nc
        fu=${CRT1}_grid_U.nc
        fv=${CRT1}_grid_V.nc
        fj=${CRT1}_${FILE_ICE_SUFFIX}.nc ; # can be icemod or grid_T ....
        fvt=${CRT1}_VT.nc


        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # If coupled EC-Earth simu, attempting to compute ocean-averaged fluxes from IFS too (E, P, E-P)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${ece_run} -eq 2 ] && [ ${NBL} -eq 75 ] && [ ${i_do_ifs_flx} -eq 1 ]; then
            echo; echo; echo "Fluxes of freshwater at the surface from IFS..."
            echo "LAUNCHING: ./src/bash/extract_ifs_surf_fluxes.sh in the background!"
            ${BARAKUDA_ROOT}/src/bash/extract_ifs_surf_fluxes.sh &
            pid_flxl=$! ; echo
        fi

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2D maps of NEMO - OBS for SST and SSS (for movies)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_movi} -eq 1 ]; then
            echo; echo; echo "2D maps of NEMO - OBS for SST and SSS (for movies)"
            echo "CALLING: prepare_movies.py ${ft} ${jyear} sst"
            ${PYTH} ${PYBRKD_EXEC_PATH}/prepare_movies.py ${ft} ${jyear} sst &
            pid_movt=$! ; echo
            echo "CALLING: prepare_movies.py ${ft} ${jyear} sss"
            ${PYTH} ${PYBRKD_EXEC_PATH}/prepare_movies.py ${ft} ${jyear} sss &
            pid_movs=$! ; echo
            echo "CALLING: prepare_movies.py ${fj} ${jyear} ice"
            ${PYTH} ${PYBRKD_EXEC_PATH}/prepare_movies.py ${fj} ${jyear} ice &
            pid_movi=$! ; echo
        fi

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Computing time-series of spatially-averaged variables
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_mean} -eq 1 ]; then
            echo; echo; echo "Global monthly values"
            echo "CALLING: mean.py ${ft} ${jyear}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/mean.py ${ft} ${jyear} &
            pid_mean=$! ; echo
        fi

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Creating VT file if needed
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_trsp} -gt 0 ] || [ ${i_do_mht} -eq 1 ]; then
            if [ ! -f ${fvt} ]; then
                echo; echo; echo " *** doing: ./cdfvT.x ${CPREF}${TTAG_ann} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV}"
                ./cdfvT.x ${CPREF}${TTAG_ann} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} &
                pid_vtvt=$! ; echo
            fi
        fi

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Computing time-series of spatially-averaged variables
        # on boxes (saving the variable on 2D box too...
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_ssx_box} -eq 1 ]; then
            echo; echo; echo "Box monthly values"
            echo "CALLING: ssx_boxes ${ft} ${jyear} ${NN_SST} ${NN_SSS}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/ssx_boxes.py ${ft} ${jyear} ${NN_SST} ${NN_SSS} &
            pid_boxa=$! ; echo
        fi

        
        echo
        echo " Gonna wait for level #1 to be done !"
        wait  ${pid_vtvt} ${pid_boxa}
        echo " .... diag level #1 done...." ; echo
        
        echo


        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Meridional heat and salt transport
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_mht} -eq 1 ]; then
            echo; echo; echo "Meridional transport of heat and salt"
            fo=${DIAG_D}/merid_transport_T_S_${CONFRUN}.nc
            if [ ! -f ${fvt} ]; then
                echo "PROBLEM: file ${fvt} is not here, skipping meridional transports section"
            else
                rm -f merid_heat_trp.dat merid_salt_trp.dat
                echo " *** doing: ./cdfmhst.x ${fvt} ${fo} ${jyear}"
                ./cdfmhst.x ${fvt} ${fo} ${jyear} &
                pid_mhst=$! ; echo
            fi
        fi
        

        #~~~~~~~
        #  MOC
        #~~~~~~~
        if [ ${i_do_amoc} -eq 1 ]; then
            echo; echo; echo "MOC"
            rm -f moc.nc *.tmp
            echo " *** doing: ./cdfmoc.x ${fv} ${NN_V} ${NN_V_EIV}"
            ./cdfmoc.x ${fv} ${NN_V} ${NN_V_EIV} &
            pid_amoc=$! ; echo
        fi

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  Transport by sigma-class
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_sigt} -eq 1 ]; then
            echo; echo
            if [ ! -f ./dens_section.dat ]; then
                if [ -f ${DENSITY_SECTION_FILE} ]; then
                    echo "Copying ${DENSITY_SECTION_FILE} to here: `pwd` !"; cp ${DENSITY_SECTION_FILE} ./dens_section.dat
                else
                    echo; echo "WARNING: Can't do Transport by sigma-class: ${DENSITY_SECTION_FILE} is missing!!!"
                fi
            fi
            echo " *** doing: ./cdfsigtrp.x ${ft} ${fu} ${fv} 24.8 28.6 19 ${jyear} ${DIAG_D}  ${NN_T} ${NN_S} ${NN_U} ${NN_V}"; echo
            ./cdfsigtrp.x ${ft} ${fu} ${fv} 24.8 28.6 19 ${jyear} ${DIAG_D} ${NN_T} ${NN_S} ${NN_U} ${NN_V} &
            pid_sigt=$! ; echo
        fi


        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # VOLUME, HEAT and SALT transports through specified section
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_trsp} -gt 0 ]; then

            echo; echo; echo "Transports of volume, heat and salt through different sections"
            if [ -z ${TRANSPORT_SECTION_FILE} ]; then
                echo "Please specify which TRANSPORT_SECTION_FILE to use into the config file!" ; exit
            fi
            if [ ! -f ./transportiz.dat ]; then
                check_if_file ${TRANSPORT_SECTION_FILE}
                cp ${TRANSPORT_SECTION_FILE} ./transportiz.dat
            fi
            if [ ${i_do_trsp} -eq 1 ]; then z1_trsp="" ; z2_trsp=""; fi
            if [ ! -f ${fvt} ]; then
                echo "PROBLEM: file ${fvt} is not here, skipping transport section!"
            else
                echo " *** doing: ./cdftransportiz.x ${CPREF}${TTAG_ann} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} ${jyear} ${DIAG_D} ${z1_trsp} ${z2_trsp}"
                ./cdftransportiz.x ${CPREF}${TTAG_ann} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} ${jyear} ${DIAG_D} ${z1_trsp} ${z2_trsp} &
                pid_trsp=$! ; echo
            fi
        fi   ; # ${i_do_trsp} -gt 0


        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  Deep Mixed Volume (DMV) on a given box
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ! -z ${i_do_dmv} ] && [ ${i_do_dmv} -gt 0 ]; then

            if [ -z ${FILE_DEF_BOXES} ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi
            echo "CALLING: dmv.py ${ft} ${cyear}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/dmv.py ${ft} ${cyear} &
            pid_dmvl=$! ; echo
        fi


        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Budget and other stuffs on a given rectangular box!
        # It provides time-series depending only on time (not depth)
        # budget_rectangle_box.py
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if [ ${i_do_bb} -gt 0 ]; then

            echo; echo; echo "Budget and other stuffs on rectangular boxes!"

            if [ -z ${FILE_DEF_BOXES} ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi
            echo " *** doing: ${PYTH} ${PYBRKD_EXEC_PATH}/budget_rectangle_box.py ${cyear} 100 uv"
            ${PYTH} ${PYBRKD_EXEC_PATH}/budget_rectangle_box.py ${cyear} 100 uv &
            pid_bbbb=$! ; echo
        fi

        # pid_sigt pid_trsp pid_dmvl pid_bbbb pid_mhst
        echo " Gonna wait for level #2 to be done !"
        wait ${pid_amoc} ; # moc needs to be done to call cdfmaxmoc.x ...
        echo
        echo " .... diag level #2 done...." ; echo ; echo


        #~~~~~~~~~~~
        #  Max AMOC
        #~~~~~~~~~~~
        if [ ${i_do_amoc} -eq 1 ]; then
            if [ -z "${LMOCLAT}" ]; then
                echo "AMOC => specify latitude bands with variable LMOCLAT into the config file!!!"; exit
            fi
            for clat in ${LMOCLAT}; do
                cslat=`echo ${clat} | sed -e s/'-'/' '/g`
                echo "  *** ./cdfmaxmoc.x moc.nc atl ${cslat} 500 1500 ${jyear} ${DIAG_D}"
                ./cdfmaxmoc.x moc.nc atl ${cslat} 500 1500 ${jyear} ${DIAG_D} &
            done
        fi


        #~~~~~~~~~
        # SEA-ICE
        #~~~~~~~~~
        if [ ${i_do_ice} -eq 1 ]; then
            echo; echo; echo "Sea-ice extent and volume..." ; rm -f tmp_ice.nc
            echo "ncks  -A -v ${NN_ICEF} ${fj} -o tmp_ice.nc"
            ncks  -A -v ${NN_ICEF} ${fj} -o tmp_ice.nc
            ncrename -v ${NN_ICEF},ice_frac tmp_ice.nc
            coic=""
            if [ -z ${NN_ICET} ]; then
                coic="oic" ; # means only ice concentration available!
            else
                echo "ncks  -A -v ${NN_ICET} ${fj} -o tmp_ice.nc"
                ncks  -A -v ${NN_ICET} ${fj} -o tmp_ice.nc
                ncrename -v ${NN_ICET},ice_thic tmp_ice.nc
            fi
            echo " *** doing: ./cdficediags.x tmp_ice.nc ${jyear} ${DIAG_D} ${coic}"
            ./cdficediags.x tmp_ice.nc ${jyear} ${DIAG_D} ${coic} &

        fi


        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  Vertical profiles of T,S and density on a given box
        #     => it provides time-series depending on time and depth
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_box_TS_z} -gt 0 ]; then
            if [ -z ${FILE_DEF_BOXES} ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi
            echo "CALLING: prof_TS_z_box.py ${cyear}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/prof_TS_z_box.py ${cyear} &
            echo;echo
        fi

        
        # --- end of stuffs that can be launched in bg ---



        #~~~~~~~~~~~~~~~~~~~~~~~
        #  AMO SST time series
        #~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_amo} -eq 1 ]; then
            diro=${DIAG_D}/amo ; mkdir -p ${diro}
            echo; echo; echo "AMO SST time series "
            echo; echo
            echo "CALLING: amo.py ${cyear}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/amo.py ${cyear}
            echo;echo;
            mv -f  AMO_SST_Atl_${CONFRUN}_${cyear}.nc ${diro}/
            echo
        fi
        # End SST time series

        # Vertical meridional or zonal sections:
        if [ ${i_do_sect} -eq 1 ]; then
            diro=${DIAG_D}/sections ; mkdir -p ${diro}
            if [ -z ${VSECT_NM} ] || [ -z ${VSECT_JI} ] || [ -z ${VSECT_JJ} ]; then
                echo "VSECT_NM, VSECT_JI and VSECT_JJ must be defined in:"
                echo "${fconfig}" ; exit
            fi
            js=0
            for cs in ${VSECT_NM[*]}; do
                cc=`echo ${ft} | sed -e s/".nc"/"_${cs}.nc"/g`
                fo="section_${cc}"
                ncks -h -O -v nav_lat,nav_lon,thetao,so \
                    -d x,${VSECT_JI[${js}]} \
                    -d y,${VSECT_JJ[${js}]} \
                    ${ft} -o ${fo}
                mv -f ${fo} ${diro}/
                ((js++))
            done
            echo
        fi


        #==============================
        # BETA STUFF !!
        #==============================

        if [ ${i_do_zcrit} -gt 0 ]; then
            if [ -z ${FILE_DEF_BOXES} ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi
            echo "CALLING: zcrit_conv.py ${cyear}"
            ${PYTH} ${PYBRKD_EXEC_PATH}/zcrit_conv.py ${cyear}
            echo;echo
        fi


        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Solid freshwater transport associated with sea-ice drift
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_icet} -eq 1 ]; then

            echo; echo; echo "Transports of solid freshwater (sea-ice) through different sections"

            if [ ! "${FILE_ICE_SUFFIX}" = "icemod" ]; then
                echo "ERROR: cannot compute ice transport if ice file is set to ${FILE_ICE_SUFFIX} !"; exit
            fi

            check_if_file ${TRANSPORT_ICE_SECTION_FILE}

            cp ${TRANSPORT_ICE_SECTION_FILE} ./transport_ice.dat

            diro=${DIAG_D}/transport_sections ; mkdir -p ${diro}

            echo "CALLING: /home/x_laubr/DEV/CDFTOOLS/bin/cdficeflux ${fj}"
            /home/x_laubr/DEV/CDFTOOLS/bin/cdficeflux ${fj}

            list_ice=`cat transport_ice.dat | grep '-'`

            for sect in ${list_ice}; do
                fo=${diro}/transport_solid_FW_${sect}_${CONFRUN}.dat
                mv -f section_ice-trp_${sect}.dat  ${sect}.tmp
                echo "# Time       VolTrans(Sv)     (${cyear})" >> ${fo}
                cat ${sect}.tmp | grep -v '\#' | awk -v "Y=${jyear}" '{print Y+($1-0.5)*1./12.," ",$2}'>> ${fo}
                cat ${sect}.tmp | grep -v '\#' | awk '{print $2}'> ${sect}_1.tmp
                # Mean val for the current year:
                fo=${diro}/transport_solid_FW_${sect}_${CONFRUN}_annual.dat
                mean_val1=`cat ${sect}_1.tmp | awk '{ SUM += $1} END { printf("%.15g\n", SUM/12) }'`
                ymid=`echo ${jyear} | awk  '{print $1+0.5}'`
                if [ ${jyear} -eq ${YEAR_INI} ]; then
                    echo "# Annually-averaged transports"  >> ${fo}
                    echo "# time       VolTrans(Sv)" >> ${fo}
                fi
                echo "${ymid}   ${mean_val1}" >> ${fo}
            done
            echo; echo; echo
        fi


        echo
        echo " Waiting for backround jobs for current year (${jyear}) !"
        wait ${pid_movt} ${pid_movs} ${pid_movi} ${pid_mean} ${pid_flxl}
        wait
        echo "  Done waiting for year ${cyear} !"
        if [ ${i_do_movi} -eq 1 ]; then rsync -avP movies ${DIAG_D}/ ; fi
        rm -f *.tmp broken_line_* tmp_ice.nc
        rm -f ${CRT1}_*.nc ; #debug
        echo ; echo

        # DIAGS ARE DONE !!!
        echo "${cyear}" > ${fcompletion}


    fi ; # if ${lcontinue}; then


    #if ${LFORCE_YEND}; then
    #    if [ ${jyear} -eq ${YEARN} ]; then lcontinue=false; fi
    #fi

    wait

    ((jyear++))


# end loop years...
done ; # while ${lcontinue}; do






l_pclim=false


# PREPARING HTML PAGE
# ~~~~~~~~~~~~~~~~~~~

if [ ${ISTAGE} -eq 2 ]; then

    rm -rf ${DIAG_D}/${RUN}
    
    y1=`cat ${DIAG_D}/first_year.info`
    y2=`cat ${DIAG_D}/last_year_done.info`
    nby=$((${y2}-${y1}+1)) ; #lulu

    if [ ${IFREQ_SAV_YEARS} -gt 1 ]; then
        fnamelist=namelist.${cy1m}-${cy2m}
    else
        fnamelist=namelist.${cy2}
    fi


    # Agreement between last year from output files and 'fcompletion' file:
    ydum=`cat ${fcompletion}`
    if [ -z ${YEARN} ]; then
        # (if YEARN is not set...)
        if [ ! ${ydum} -eq ${YEAR_END} ]; then
            echo;
            echo "###################################################################"
            echo "PROBLEM: in ${fcompletion} last_year = ${ydum}"
            echo "         and from stored files files last_year = ${YEAR_END} !"
            echo "###################################################################"
            echo
            exit
        fi
    fi


    

    echo; echo; echo "RUN ${RUN}: creating plots"; echo

    #cd ${BARAKUDA_ROOT}/

    cd ${DIAG_D}/


    if [ ${i_do_movi} -eq 1 ]; then
        echo
        idelay=$((120-${nby}*8))
        if [ ${idelay} -lt 10 ]; then idelay=10; fi
        echo "Creating GIF movies out of the 2D maps of NEMO - OBS for SST and SSS (delay=${idelay})"
        rm -f *_${CONFRUN}.gif
        convert -delay ${idelay} -loop 0 movies/dsst_*.png dsst_${CONFRUN}.gif &
        convert -delay ${idelay} -loop 0 movies/dsss_*.png dsss_${CONFRUN}.gif &
        convert -delay ${idelay} -loop 0 movies/icen_*.png icen_${CONFRUN}.gif &
        convert -delay ${idelay} -loop 0 movies/ices_*.png ices_${CONFRUN}.gif &
    fi

    # 1D plots to perform
    # ~~~~~~~~~~~~~~~~~~~

    DIAG_1D_LIST=""

    if [ ${i_do_mean} -eq 1 ]; then
        DIAG_1D_LIST="${DIAG_1D_LIST} 3d_so mean_sos 3d_thetao mean_tos mean_zos mean_fwf mean_htf mean_mldr10_1"
    fi
    if [ ${i_do_amoc} -eq 1 ]; then DIAG_1D_LIST="${DIAG_1D_LIST} amoc";        fi
    if [ ${i_do_trsp} -gt 0 ]; then DIAG_1D_LIST="${DIAG_1D_LIST} transport_sections" ; fi
    if [ ${i_do_ice}  -eq 1 ]; then DIAG_1D_LIST="${DIAG_1D_LIST} seaice";       fi

    dy=$((${YEAR_END}-${YEAR_INI}+1)) ; export YF2=$((${YEAR_END}+1))


    # Doing 1D plots
    # ~~~~~~~~~~~~~~    

    echo ; echo; echo "Going to perform the following 1D plots:"
    echo "    => ${DIAG_1D_LIST}"; echo

    for fd in ${DIAG_1D_LIST}; do
        echo "CALLING: plot_time_series.py ${fd}"
        ${PYTH} ${PYBRKD_EXEC_PATH}/plot_time_series.py ${fd} ; echo
        echo
    done
    echo ; echo ; echo

    if [ ${i_do_mean} -eq 1 ]; then

         # 5-month-running mean SST anomaly on Nino region 3.4 graph:
        echo "CALLING: plot_enso.py Nino34_${CONFRUN}.nc"
        ${PYTH} ${PYBRKD_EXEC_PATH}/plot_enso.py Nino34_${CONFRUN}.nc
        echo; echo

        # Hovmuller of temperature and salinity
        echo "CALLING: plot_hovm_tz.py"
        ${PYTH} ${PYBRKD_EXEC_PATH}/plot_hovm_tz.py
        echo; echo
        #
    fi


    if [ ${i_do_sigt} -eq 1 ]; then
        # Transport by sigma-class
        echo "CALLING: plot_trsp_sigma.py"
        ${PYTH} ${PYBRKD_EXEC_PATH}/plot_trsp_sigma.py
        echo; echo; echo
    fi


    if [ ${i_do_mht} -eq 1 ]; then
        #
        # Hovmullers of advective meridional heat/salt transport
        echo; echo
        echo "CALLING: plot_hovm_merid_trsp.py"
        ${PYTH} ${PYBRKD_EXEC_PATH}/plot_hovm_merid_trsp.py
        echo; echo; echo
        #
    fi



    echo; echo; echo




    if ${l_clim_diag} ; then

        ###########################################################################
        # Climatology over X years (12-month file average of X consecutive years)
        #   => has to be built with the 'build_clim.sh' script
        ###########################################################################

        echo; echo; echo "Checking for presence of ${DIAG_D}/clim/last_clim..."
        if [ -f ${DIAG_D}/clim/last_clim ]; then
            cat ${DIAG_D}/clim/last_clim
            export CLIM_PER=`cat ${DIAG_D}/clim/last_clim`
            ftcli=${DIAG_D}/clim/mclim_${CONFRUN}_${CLIM_PER}_grid_T.nc4
            ficli=${DIAG_D}/clim/mclim_${CONFRUN}_${CLIM_PER}_${FILE_ICE_SUFFIX}.nc4
            fclvt=${DIAG_D}/clim/aclim_${CONFRUN}_${CLIM_PER}_VT.nc4
            fcmoc=${DIAG_D}/clim/aclim_${CONFRUN}_${CLIM_PER}_MOC.nc4
            fcpsi=${DIAG_D}/clim/mclim_${CONFRUN}_${CLIM_PER}_PSI.nc4
            iclyear=`echo ${CLIM_PER} | sed -e s/'-'/' '/g`
        else
            echo; echo "PROBLEM! => you set l_clim_diag to true but no file 'last_clim' was found in:"
            echo "            ${DIAG_D}/clim/"; echo
            exit
        fi

        echo; echo; echo "Checking for presence of ${ftcli}..."
        if [ -f ${ftcli} ]; then
            echo; echo;
            echo "~~~~~~~~~~~~~~~~~~~~~"
            echo "*   CLIMATO FOUND !!!"
            echo "~~~~~~~~~~~~~~~~~~~~~"
            echo
            echo "  => for years ${CLIM_PER}" ; echo "  => using ${ftcli}"

            list_comp_2d="CLIM"
            l_pclim=true
            lcomp_to_run=false

            if [ ! -z ${RUNREF} ]; then
                lcomp_to_run=true
                list_comp_2d="CLIM ${RUNREF}"
                # Must check if climatology for run ${RUNREF} is there:
                fclim_ref=`echo "${ftcli}" | sed -e "s|${RUN}|${RUNREF}|g"`
                check_if_file ${fclim_ref}
                echo "Going to compare also against run ${fclim_ref}!"
                echo
            fi


            echo; echo
            check_if_file ${F_T_CLIM_3D_12} "name:F_T_CLIM_3D_12"
            check_if_file ${F_S_CLIM_3D_12} "name:F_S_CLIM_3D_12"
            check_if_file ${SST_CLIM_12}    "name:SST_CLIM_12"
            if [ ${i_do_ice}  -gt 0 ]; then check_if_file ${ICE_CLIM_12}    "name:ICE_CLIM_12" ; fi
            echo; echo


        #######################################
        # Diags that don't imply a comparison #
        #######################################

            export COMP2D="CLIM"



            # Lat-Depth AMOC
            # ~~~~~~~~~~~~~~
            if [ -f ${fcmoc} ]; then
                echo; echo
                echo " Ploting lat-depth MOC !"
                cd ${DIAG_D}/
                export DIRS_2_EXP="${DIRS_2_EXP} moc"
                rm -rf moc; mkdir moc; cd moc/
                echo; echo; echo "CALLING: moc.py ${iclyear}"
                ${PYTH} ${PYBRKD_EXEC_PATH}/moc.py ${iclyear}
                cd ../
                echo
            fi


            # March Mixed layer depth in Nordic Seas
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if [ `ipresent_var_in_ncf ${ftcli} ${NN_MLD}` -eq 1 ]; then
                echo; echo
                echo " Performing 2D mapping of March Mixed layer depth in Nordic Seas"
                cd ${DIAG_D}/
                export DIRS_2_EXP="${DIRS_2_EXP} mld"
                rm -rf mld; mkdir mld; cd mld/
                echo; echo; echo "CALLING: mld.py ${iclyear}"; echo
                ${PYTH} ${PYBRKD_EXEC_PATH}/mld.py ${iclyear}
                cd ../
                echo
            fi

            # Sea-ice extent stereographic polar projection South and North
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if [ ${i_do_ice}  -gt 0 ] && [ `ipresent_var_in_ncf ${ficli} ${NN_ICEF}` -eq 1 ]; then
                echo; echo
                echo " Performing 2D Sea-ice extent stereographic polar projection South and North"
                cd ${DIAG_D}/
                export DIRS_2_EXP="${DIRS_2_EXP} sea_ice"
                rm -rf sea_ice; mkdir sea_ice; cd sea_ice/
                echo; echo; echo "CALLING: ice.py ${iclyear}"; echo
                ${PYTH} ${PYBRKD_EXEC_PATH}/ice.py ${iclyear}
                cd ../
                echo
            fi



        # Sea-surface height
        # ~~~~~~~~~~~~~~~~~~
            if [ `ipresent_var_in_ncf ${ftcli} ${NN_SSH}` -eq 1 ]; then
                echo; echo
                echo " SSH map"
                export DIRS_2_EXP="${DIRS_2_EXP} ssh"
                cd ${DIAG_D}/
                rm -rf ssh; mkdir ssh; cd ssh/
                echo; echo; echo "CALLING: ssh.py ${iclyear}"; echo
                ${PYTH} ${PYBRKD_EXEC_PATH}/ssh.py ${iclyear}
                cd ../
                echo
            else
                echo; echo "WARNING: did not find ${NN_SSH} into ${ftcli} !!!!"; echo
            fi




        ##################################################
        # Diags that imply a comparison against "COMP2D" #
        ##################################################

            cd ${DIAG_D}/
            rm -rf temp_sal surf_fluxes

            for COMP2D in ${list_comp_2d}; do

                export COMP2D=${COMP2D}

                echo; echo; echo "Clim. comparisons against ${COMP2D}"

                if [ "${COMP2D}" = "${RUNREF}" ]; then
                    export F_T_CLIM_3D_12=${fclim_ref}; check_if_file ${F_T_CLIM_3D_12} "name:F_T_CLIM_3D_12"
                    export F_S_CLIM_3D_12=${fclim_ref}; check_if_file ${F_S_CLIM_3D_12} "name:F_S_CLIM_3D_12"
                    export SST_CLIM_12=${fclim_ref}   ; check_if_file ${SST_CLIM_12}    "name:SST_CLIM_12"
                    if [ ${i_do_ice}  -gt 0 ]; then export ICE_CLIM_12=${fclim_ref}   ; check_if_file ${ICE_CLIM_12}    "name:ICE_CLIM_12"; fi
                fi


            # Temperature and Salinity
            # ~~~~~~~~~~~~~~~~~~~~~~~~
                echo; echo
                echo " Performing 2D Temperature and Salinity maps and sections"
                cd ${DIAG_D}/
                export DIRS_2_EXP="${DIRS_2_EXP} temp_sal"
                DIRS_2_EXP_RREF="${DIRS_2_EXP_RREF} temp_sal"
                mkdir -p temp_sal; cd temp_sal/
                echo; echo; echo "CALLING: temp_sal.py ${iclyear}"; echo
                ${PYTH} ${PYBRKD_EXEC_PATH}/temp_sal.py ${iclyear}
                cd ../
                echo


            done ; # for COMP2D in ${list_comp_2d}; do

        else
            echo; echo
            echo " NO CLIMATO FOUND ...";
            echo "   => you can use 'build_clim.sh' to build a climato of your run"
            echo; echo
        fi

    fi ; # if ${l_clim_diag}


    
    wait

    # Time for HTML stuff!

    export HTML_DIR=${DIAG_D}/${RUN}
    mkdir -p ${HTML_DIR}
    
    cd ${DIAG_D}/

    # Moving all figures to HTML_DIR:
    for fp in ${FIG_FORM} svg gif; do mv -f *.${fp} ${HTML_DIR}/ >/dev/null 2>/dev/null ; done
    mv -f ./merid_transport/*.${FIG_FORM} ${HTML_DIR}/ >/dev/null 2>/dev/null
    
    . ${BARAKUDA_ROOT}/src/bash/build_html.bash
    
    wait
    # Building main index.html 
    build_index_html
    
    # If climatology built, sub 2D html pages
    if ${l_pclim}; then
        build_sub_html
    fi
    
    #==================================================================
    

    wait ; # likely waiting for the creation of the GIFs....

    echo; echo

    cp ${BARAKUDA_ROOT}/src/html/conf_*.html ${HTML_DIR}/
    cp ${BARAKUDA_ROOT}/src/html/logo.png    ${HTML_DIR}/

    mv -f index.html ${HTML_DIR}/

    cp -r ${DIRS_2_EXP} ${HTML_DIR}/ >/dev/null 2>/dev/null

    echo; echo; echo

    if [ ${ihttp} -eq 1 ]; then
        echo "Preparing to export to remote host!"; echo
        send_dir=`basename ${HTML_DIR}`
        tar cvf ${send_dir}.tar ${send_dir}
        ssh ${RUSER}@${RHOST} "mkdir -p ${RWWWD}"
        scp ${send_dir}.tar ${RUSER}@${RHOST}:${RWWWD}/
        ssh ${RUSER}@${RHOST} "cd ${RWWWD}/; rm -rf ${send_dir}; tar xf ${send_dir}.tar 2>/dev/null; rm -f ${send_dir}.tar; \
            chmod -R a+r ${send_dir}; cd ${send_dir}/; source-highlight -i ${fnamelist} -s fortran -o namelist.html"
        echo; echo
        echo "Diagnostic page installed on  http://${RHOST}${RWWWD}/${send_dir}/ !"
        echo "( Also browsable on local host in ${DIAG_D}/${send_dir} )"
        rm -f ${send_dir}.tar

    else
        if [ ${ihttp} -eq 0 ]; then
            echo "Diagnostic page installed in ${HTML_DIR}/"
            echo " => you can view it with a web browser..."
        else
            echo "Error: \"ihttp\" must be either 0 or 1 !"
        fi
    fi

    echo; echo

    rm -rf *.eps


else
    echo; echo "Diagnostics built! Run \"${0} -e\" to create figure and HTML page..."; echo


fi


#debug:
rm -rf ${TMP_DIR} 2>/dev/null

echo

