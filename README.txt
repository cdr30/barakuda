~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Quick getting-started guide to BaraKuda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


I / What do you need to be able to use BaraKuda ?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- A FORTRAN 90 compiler

- netcdf library with support for the appropriate F90 compiler

- NCO

- 'convert' from ImageMagick if you want to generate GIF movies (i_do_movi=1) 

- For time-series and 2D plots, the following up-to-date packages:
  => python-netcdf4 (from netCDF4 import Dataset) and Matplotlib
  => for map projections you'll also need the Basemap package
  
  A good idea is to install a shiny python distribution, something like Canopy:
  => https://www.enthought.com/products/canopy/

  In any case, specify the appropriate "PYTHON_HOME" environment variable in
  your config/conf_<MYCONF>.sh file

- NEMO output data! => A directory containing the MONTHLY-AVERAGED, global
                       (rebuilt), NEMO output to analyze
  (grid_T, grid_U, grid_V and icemod files) as "*.nc", "*.nc.gz" or ".nc4"

- a NEMO mesh_mask file and the the corresponding basin_mask (ocean basins).
  (variables MM_FILE and BM_FILE into your config/conf_<MYCONF>.sh file)
  To create the NEMO mesh_mask.nc just launch the relevant NEMO experiment with the
  namelist parameter nn_msh set to 1 !
  If you use ORCA1 or ORCA025 you can use the
  "<ORCA>_create_basin_mask_from_meshmask.py" in python/exec to generate the
  basin file!
              tmaskatl(y, x) => "Atlantic Basin" ;
              tmaskpac(y, x) => "Pacific Basin" ;
              tmaskind(y, x) => "Indian Basin" ;



II / Compile CDFTOOLS executables 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 * CDFTOOLS is a set of FORTRAN executables intended to perform a multitude of
   diagnostics on NEMO output file and is developed by Jean-Marc Molines at LGGE
   in Grenoble.  However, this is a slightly modified light version here...  SO
   DO NOT USE AN OFFICIAL CDFTOOLS DISTRIBUTION, stick to the one that comes
   with BaraKuda!

- move to the 'barakuda/cdftools_light' directory

- configure your own 'make.macro' for your system (some templates for gfortran
  and Intel are provided...)
    => just copy or link your own "macro.your_arch" to "make.macro" !
    => F90 compiler and related netcdf library to use

- compile with 'gmake'

- if that was successful the 'barakuda/bin' directory should contain the 8
  following executables:
 * cdficediags.x
 * cdfmaxmoc.x
 * cdfmhst.x
 * cdfmoc.x
 * cdfpsi.x
 * cdfsigtrp.x
 * cdftransportiz.x
 * cdfvT.x



III / Create and configure your own "configs/config_<MY_CONF>.sh"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPORTANT: Always refer to the 'configs/config_TEMPLATE.sh' config file as a
reference when creating/re-adjusting yours. It is a symbolic link pointing to
the last officially supported and most up-to-date config file.  It should be
sufficiently well commented for you to be able to adjust your own config file.

NEMO output files must be monthly averages and of the following form:
==> <RUN NAME>_1m_<YEAR>0101_<YEAR>1231_<GRID_TYPE>.nc(.gz)
           (GRID_TYPE=grid_T/grid_U/grid_V/icemod) 

Gzipped or not!

All files for all years must all be saved in the same directory (see
NEMO_OUT_STRCT in the config file). Better if this directory only contains NEMO
output files and nothing else!

Alternatively NEMO files can be saved/organized in sub-directories a la
EC-Earth: (ex: year 1995 of run started in 1990 is the 6th year so files for
1995 are saved into sub-directory (of NEMO_OUT_STRCT) "006" (set 'ece_run' to 1
or 2 then).

If you want to perform the "climatology" plots (see section IV) you will need
some 2D and 3D climatologies of T and S interpolated on the relevant ORCA
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

A/ Once the previous job has finished to run, launch

   * To only generate time-series plots use the "-e" switch:

   ./barakuda.sh -C <MY_CONF> -R <RUN> -e
   (ex: ./barakuda.sh -C ORCA1_L75_v36_triolith -R SL36C00 -e)

   * To generate time-series + 2D climatological plots use the "-E" switch,
     provided you have built a climatology out of your run with the
     "build_clim.sh" script (see point V/B):

   ./barakuda.sh -C <MY_CONF> -R <RUN> -E


B/ To be able to create the "climatology" plots (maps, sections, etc, based on a
   monthly climatology of a few years) you will have to

  i) create the climatology with the "build_clim.sh" script:
   
   ./build_clim.sh -C <MY_CONF> -R <RUN> -i <first_year> -e <last_year>
         (check ./build_clim.sh -h to see the other important options...)

    (ex: ./build_clim.sh -C ORCA1_L75_v36_triolith -R SL36C00 -f 10 -i 0090 -e 0099 -k 4)
      

  ii) the you can tell "barakuda.sh" to create climatology-related plots by
       using the "-E" switch instead of "-e" (see point V/A)


C/ To compare time-series between at least 2 (already diagnosed) runs:
   
   ./compare_time-series.sh -C <MY_CONF> -R <RUN1>,<RUN2>,...,<RUNn>
   (ex: ./compare_time-series.sh -C ORCA1_L75_v36_triolith -R SL36C00,SL36EIE )

