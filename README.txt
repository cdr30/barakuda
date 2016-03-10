~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Quick getting-started guide to BaraKuda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


I / What do you need to be able to use BaraKuda ?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- A Fortan 90 compiler

- netcdf library with support for the appropriate F90 compiler

- NCO 

- For time-series and 2D plots, the following up-to-date packages:
  => python-netcdf4 (from netCDF4 import Dataset) and Matplotlib
  => for map projections you'll also need the Basemap package
  
  A good idea is to install a shiny python distribution, something like Canopy:
  => https://www.enthought.com/products/canopy/

- NEMO output data! => A directory containing the monthly NEMO output to analyze
               (grid_T, grid_U, grid_V and icemod files) as "*.nc", "*.nc.gz" or ".nc4"
               (it might also be a good idea to save the SBC files too...)

- a NEMO mesh_mask file and the the corresponding basin_mask (ocean basins).
  (variables MM_FILE and BM_FILE into your config/conf_<MYCONF>.sh file)
  To create the NEMO mesh_mask.nc just launch the relevant NEMO experiment with the
  namelist parameter nn_msh set to 1 !
  If you use ORCA1 you can use the "orca1_create_basin_mask_from_meshmask.py" in python/exec
  to generate the basin file!
              tmaskatl(y, x) => "Atlantic Basin" ;
              tmaskpac(y, x) => "Pacific Basin" ;
              tmaskind(y, x) => "Indian Basin" ;



II / Compile CDFTOOLS executables 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 * CDFTOOLS is a set of fortran executables intended to perform a multitude of
   diagnostics on NEMO output file and is developed by Jean-Marc Molines at LEGI
   in Grenoble.  However, this is a slightly modified light version here...

- move to the 'barakuda/cdftools_light' directory

- configure your own 'make.macro' for your system (some templates for gfortran and Intel are provided...)
    => just copy or link your own "macro.your_arch" to "make.macro" !
    => F90 compiler and related netcdf library to use

- compile with 'gmake'

- if that was successful the 'barakuda/bin' directory should contain the following executables
    => cdfcurl.x cdfmaxmoc.x cdfmhst.x cdfmoc.x cdfpsi.x  cdftransportiz.x  cdfzonalmean.x
       cdfhflx.x cdfmeanvar.x cdfmocatl.x cdfmoy.x cdfrmsssh.x cdfvT.x cdficediags.x  cdfmean.x
       cdfmocsig.x cdfmxl.x cdfsigtrp.x cdfw.x



III / Create and configure your own "configs/config_<MY_CONF>.sh" from an existing one
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NEMO output files must be stored in 
<STORE_DIR>/<ORCA_GRID>/<ORCA_GRID>-<RUN>-S
(ex: /nobackup/vagn2/x_laubr/ORCA1/ORCA1-SLB3500-S)

Files there should be monthly averages and of the following form:
==> <RUN NAME>_MM_<YEAR>0101_<YEAR>1231_<GRID_TYPE>.nc(.gz)   (GRID_TYPE=grid_T/grid_U/grid_V/icemod) 

 Gzipped or not! And only these files, nothing else!!!

If you want to perform the "climatology" plots (see section IV) you will need
some   2D and 3D climatologies of T and S interpolated on the relevant ORCA
grid. Usually you should already have them since they are needed for your
simulation. These are the following files in your BaraKuda config file:
F_T_CLIM_3D_12, F_S_CLIM_3D_12, SST_CLIM_12



IV) Create diagnostics
~~~~~~~~~~~~~~~~~~~~~~

Launch "barakuda.sh" 
./barakuda.sh -C <MY_CONF> -R <RUN>
(ex: ./barakuda.sh -C ORCA1_L75_v36_triolith -R SL36C00)

Use the -h switch to see available options



V) Create figures and browsable HTML page
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A/ Once the previous job has finished to run, launch:

   ./barakuda.sh -C <MY_CONF> -R <RUN> -e
   (ex: ./barakuda.sh -C ORCA1_L75_v36_triolith -R SL36C00 -e)

B/ If you want to perform the "climatology" plots (maps, sections, etc, based on a
   monthly climatology of a few years) you will have to:

   i) create the climatology with the "build_clim.sh"
      => EX: $ ./build_clim.sh -C ORCA1_L75_v36_triolith -R SL36C00 -f 10 -i 0090 -e 0099 -k 4
         (check ./build_clim.sh -h to see the options)

   ii) set "l_clim_diag=true" in your config file "configs/config_<MY_CONF>.sh"

C/ If you want to create time-series comparing 2 runs (each already diagnosed, at least stage III):
   
   ./compare_time-series.sh -C <MY_CONF> -R <RUN1>,<RUN2>,...,<RUNn>
   (ex: ./compare_time-series.sh -C ORCA1_L75_v36_triolith -R SL36C00,SL36EIE )
