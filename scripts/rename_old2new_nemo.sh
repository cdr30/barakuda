#!/bin/sh

if [ "${1}" = "" ]; then
    echo "USAGE: ${0} CONFRUN_DATE"; exit
fi

fold_t=${1}_grid_T_old.nc
fnew_t=${1}_grid_T.nc
mv -f ${fnew_t} ${fold_t}

ncrename -O -v votemper,thetao -v vosaline,so -v sosstsst,tos -v sosaline,sos -v sossheig,zos -v sowaflup,wfo \
    -v soshfldo,rsntds -v sohefldo,tohfls -v somxl010,mldr10_1 -v somixhgt,mldkz5 -v soicecov,ice_cover \
    -v sowafldp,erp \
    ${fold_t} -o ${fnew_t}



fold_u=${1}_grid_U_old.nc
fnew_u=${1}_grid_U.nc
mv -f ${fnew_u} ${fold_u}

ncrename -O -v sozotaux,tauuo -v vozocrtx,uo  ${fold_u} -o ${fnew_u}


fold_v=${1}_grid_V_old.nc
fnew_v=${1}_grid_V.nc
mv -f ${fnew_v} ${fold_v}

ncrename -O -v sometauy,tauvo -v vomecrty,vo  ${fold_v} -o ${fnew_v}

fold_i=${1}_icemod_old.nc
fnew_i=${1}_icemod.nc
mv -f ${fnew_i} ${fold_i}

ncrename -O -v iicethic,sit -v iicevelu,uice_ipa -v iicevelv,vice_ipa -v iicesurt,ist_ipa \
    ${fold_i} -o ${fnew_i}


rm -f ${1}_*_old.nc









#    -v soshfldo,rsntds -v sohefldo,tohfls -v somxl010,somxl010 -v somixhgt,mldkz5 -v soicecov,ice_cover \
#    -v sowafldp,erp \


#exit


#		nav_lon:long_name = "Longitude" ;
#		nav_lat:long_name = "Latitude" ;
#		deptht:long_name = "Vertical T levels" ;
#		time_counter:long_name = "Time axis" ;
#		votemper:long_name = "Temperature" ;
#		vosaline:long_name = "Salinity" ;
#		sosstsst:long_name = "Sea surface temperature" ;
#		sosaline:long_name = "Sea surface salinity" ;
#		sossheig:long_name = "Sea surface height" ;
#		sowaflup:long_name = "Net Upward Water Flux" ;
#		soshfldo:long_name = "Shortwave Radiation" ;
#		sowaflcd:long_name = "Concentration/dilution water flux" ;
#		sohefldo:long_name = "Net Downward Heat Flux" ;
#		somxl010:long_name = "Mixed Layer Depth 0.01 ref.10m" ;
#		somixhgt:long_name = "Mixing layer depth (Turbocline)" ;
#		soicecov:long_name = "Ice fraction" ;
#		sowindsp:long_name = "Wind speed module at 10 m" ;
#		sohefldp:long_name = "Surface Heat Flux: Damping" ;
#		sowafldp:long_name = "Surface Water Flux: Damping" ;
#		sothedep:long_name = "Thermocline Depth (max dT/dz)" ;
#		so20chgt:long_name = "Depth of 20C isotherm" ;
#		so28chgt:long_name = "Depth of 28C isotherm" ;
#		sohtc300:long_name = "Heat content 300 m" ;
#		soicetem:long_name = "Ice surface temperature (ice presence average)" ;
#		soicealb:long_name = "Ice albedo (cell average)" ;



#		deptht:long_name = "Vertical T levels" ;
#		erp:long_name = "Surface Water Flux: Damping" ;
#		ice_cover:long_name = "Ice fraction" ;
#		mldkz5:long_name = "mixing layer depth (Turbocline)" ;
#		mldr10_1:long_name = "Mixed Layer Depth 0.01 ref.10m" ;
#		nav_lat:long_name = "Latitude" ;
#		nav_lon:long_name = "Longitude" ;
#		qla_oce:long_name = "Latent Downward Heat Flux over open ocean" ;
#		qlw_oce:long_name = "Longwave Downward Heat Flux over open ocean" ;
#		qns:long_name = "non solar Downward Heat Flux" ;
#		qsb_oce:long_name = "Sensible Downward Heat Flux over open ocean" ;
#		rsntds:long_name = "surface_net_downward_shortwave_flux" ;
#		runoffs:long_name = "River Runoffs" ;
#		so:long_name = "sea_water_salinity" ;
#		sos:long_name = "sea_surface_salinity" ;
#		taum:long_name = "wind stress module" ;
#		thetao:long_name = "sea_water_potential_temperature" ;
#		time_average_1mo:long_name = "Time axis" ;
#		tohfls:long_name = "surface_net_downward_total_heat_flux" ;
#		tos:long_name = "sea_surface_temperature" ;
#		wfo:long_name = "water_flux_into_sea_water" ;
#		zos:long_name = "sea_surface_height_above_geoid" ;
#
#
#