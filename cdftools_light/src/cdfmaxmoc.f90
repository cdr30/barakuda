PROGRAM cdfmaxmoc
  !!---------------------------------------------------------------------------------------------------
  !!              ***  PROGRAM cdfmaxmoc  ***
  !!
  !!   ** Purpose : Compute the maximum of the overturning fonction from a file calculated by cdfmoc
  !!
  !!   ** Method : maxovt 'ovtfile' latmin latmax depmin depmax
  !!               return ovtmaximum and ovt minimum in the defined range.
  !!               Also give location of those extrema
  !!               works for Atlantic and Global MOC
  !!
  !!  * history:
  !!              July 2005 : original : J.M. Molines
  !!              November :  modified and adapted to cdf output R. Dussin.
  !!---------------------------------------------------------------------------------------------------
  !!  $Rev: 319 $
  !!  $Date: 2010-05-19 12:19:07 +0200 (Wed, 19 May 2010) $
  !!  $Id: cdfmaxmoc.f90 319 2010-05-19 10:19:07Z dussin $
  !!--------------------------------------------------------------

  USE netcdf

  USE cdfio
  USE io_ezcdf

  IMPLICIT NONE
  !!
  INTEGER :: jj, jk, jt,    jt_pos             ! dummy loop index  !LB
  INTEGER :: id_moc, id_lat_max, id_depthw_max, idd_t, idv_time, idf_out
  LOGICAL :: lfncout
  INTEGER :: npjglo, npk, nt           ! size of the overturning !LB
  INTEGER :: narg, iargc, istatus       ! line command stuff
  INTEGER :: jmin, jmax, kmin, kmax    ! (latitude, depth) window where to look at extrema
  INTEGER :: imaxloc(3)    ! temporary array to use with minloc/maxloc
  ! added to write in netcdf
  INTEGER :: kx=1, ky=1, kz=1          ! dims of netcdf output file
  INTEGER :: nboutput=6                ! number of values to write in cdf output
  INTEGER :: ncout, ierr               ! for netcdf output
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout
  !
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE ::  zomoc   ! zonal MOC (1,npjglo,jpk)
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::    rlat    ! latitude (1, npjglo)
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE   ::   gdepw    ! depth read in the header

  REAL(KIND=4), DIMENSION(:), ALLOCATABLE ::   ovtmax
  INTEGER,      DIMENSION(:), ALLOCATABLE ::   jlatmax, kdmax


  CHARACTER(len=4)                          ::   clatmin, clatmax
  REAL(KIND=4)                              ::   rlatmin, rlatmax, depmin , depmax
  ! added to write in netcdf
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  dumlon, dumlat
  REAL(KIND=4), DIMENSION (1)               ::  tim ! time counter
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar  ! structure of output
  !
  CHARACTER(LEN=256) :: cdum, cfile, comment, cbasin, cvar, cdep, ctim, cv_dum
  CHARACTER(LEN=512) :: cd_out='.', cf_out


  REAL(4) :: ryear, rmon, rdt

  INTEGER :: idf_moc=0 , idv_moc=0


  ! * main program
  narg=iargc()
  IF (narg == 8) THEN
     CALL getarg(1,cfile)
     CALL getarg(2,cbasin)
     CALL getarg(3,cdum)
     READ(cdum,*) rlatmin
     CALL getarg(4,cdum)
     READ(cdum,*) rlatmax
     CALL getarg(5,cdum)
     READ(cdum,*) depmin
     CALL getarg(6,cdum)
     READ(cdum,*) depmax
     CALL getarg(7,cdum)
     READ(cdum,*) ryear
     CALL getarg(8,cd_out)
  ELSE
     PRINT *,' USAGE: cdfmaxmoc ''ovt_file.nc'' cbasin latmin latmax depmin depmax <year> <DIROUT>'
     PRINT *,'        cbasin is one of atl glo inp ind or pac '
     PRINT *,' Output on standard output by default and maxmoc.nc '
     STOP
  ENDIF



  npjglo = getdim(cfile,'y')
  cdep='none'
  npk    = getdim(cfile,'depth',cdtrue=cdep,kstatus=istatus) !LB
  ctim='none'
  nt     = getdim(cfile,'time', cdtrue=ctim,kstatus=istatus) !LB


  ALLOCATE ( zomoc (1,npjglo,npk) ,gdepw(npk), rlat(1,npjglo), ovtmax(nt), jlatmax(nt), kdmax(nt) )
  gdepw(:)  = -getvar1d(cfile,'depthw',npk)
  rlat(:,:) = getvar(cfile,'nav_lat',1,1,npjglo)

  SELECT CASE (cbasin)
  CASE ('atl')
     cvar='zomsfatl'
  CASE ('glo')
     cvar='zomsfglo'
  CASE ('pac')
     cvar='zomsfpac'
  CASE ('inp')
     cvar='zomsfinp'
  CASE ('ind')
     cvar='zomsfind'
  CASE DEFAULT
     STOP 'basin not found'
  END SELECT


  ! look for jmin-jmax :
  DO jj=1, npjglo
     IF ( rlat(1,jj) <= rlatmin )  jmin = jj
     IF ( rlat(1,jj) <= rlatmax )  jmax = jj
  END DO

  ! look for kmin kmax
  DO jk=1,npk
     IF ( gdepw(jk) <= depmin ) kmin = jk
     IF ( gdepw(jk) <= depmax ) kmax = jk
  END DO


  DO jt = 1, nt !LB

     CALL GETVAR_3D(idf_moc, idv_moc,  cfile, cvar, nt, jt, zomoc)


     !DO jk=1,npk
     !   zomoc (:,:,jk) = getvar(cfile, cvar, jk, 1, npjglo, ktime=jt)
     !END DO

     ! look for max/min overturning
     ovtmax(jt) = MAXVAL(zomoc(1,jmin:jmax,kmin:kmax))

     ! find location of min/max
     !iminloc =MINLOC(zomoc(:,jmin:jmax,kmin:kmax))
     imaxloc =MAXLOC(zomoc(:,jmin:jmax,kmin:kmax))

     ! results from minloc/maxloc is relative to the sub -array given as arguments
     jlatmax(jt) = imaxloc(2)+jmin -1
     kdmax(jt)   = imaxloc(3)+kmin -1
     !jlatmin= iminloc(2)+jmin -1 ;
     !kdmin  = iminloc(3)+kmin -1 ;

     !LB:
     !IF ( jt == 1) PRINT *, '#      Time    Max MOC (Sv)     Lat         Depth'
     !PRINT *, jt, ovtmax, rlat(1,jlatmax), gdepw(kdmax)
     !LB.

  END DO !LB nt



  WRITE( clatmin, '(i2.2)') INT(ABS(rlatmin))
  WRITE( clatmax, '(i2.2)') INT(ABS(rlatmax))

  IF ( rlatmin >= 0. ) THEN
     clatmin = '+'//trim(clatmin)//'N'
  ELSE
     clatmin = '-'//trim(clatmin)//'N'
  END IF
  IF ( rlatmax >= 0. ) THEN
     clatmax = '+'//trim(clatmax)//'N'
  ELSE
     clatmax = '-'//trim(clatmax)//'N'
  END IF



  WRITE( cf_out , '(a,"/max_moc_",a,"_",2a,".nc")' ) trim(cd_out), trim(cbasin), clatmin, clatmax



  PRINT *, trim(cf_out)

  id_moc = 0 ; id_lat_max = 0 ; id_depthw_max = 0

  INQUIRE( FILE=cf_out, EXIST=lfncout )


  IF ( .NOT. lfncout ) THEN

     !! Creating file
     PRINT *, ' Creating file '//trim(cf_out)//' !!!'
     ierr = NF90_CREATE(cf_out, NF90_CLOBBER, idf_out)
     ierr = NF90_DEF_DIM(idf_out, 'time', NF90_UNLIMITED, idd_t)
     ierr = NF90_DEF_VAR(idf_out, 'time', NF90_DOUBLE,    idd_t, idv_time)


     ierr = NF90_DEF_VAR(idf_out, 'moc_'//trim(cbasin), NF90_FLOAT, (/idd_t/), id_moc)
     ierr = NF90_PUT_ATT(idf_out, id_moc, 'long_name', 'Meridional Overturning Circulation')
     ierr = NF90_PUT_ATT(idf_out, id_moc, 'units', 'Sv')

     ierr = NF90_DEF_VAR(idf_out, 'lat_max_'//trim(cbasin), NF90_FLOAT, (/idd_t/), id_lat_max)
     ierr = NF90_PUT_ATT(idf_out, id_lat_max, 'long_name', 'Latitude for maximum of MOC')
     ierr = NF90_PUT_ATT(idf_out, id_lat_max, 'units', 'deg.N')

     ierr = NF90_DEF_VAR(idf_out, 'depthw_max_'//trim(cbasin), NF90_FLOAT, (/idd_t/), id_depthw_max)
     ierr = NF90_PUT_ATT(idf_out, id_depthw_max, 'long_name', 'W-depth for maximum of MOC')
     ierr = NF90_PUT_ATT(idf_out, id_depthw_max, 'units', 'm')


     ierr = NF90_PUT_ATT(idf_out, NF90_GLOBAL, 'About', 'Created by BaraKuda (cdfmaxmoc.f90), contact: brodeau@gmail.com')

     ierr = NF90_ENDDEF(idf_out)
     jt_pos = 0

  ELSE

     !! Opening already existing file
     ierr = NF90_OPEN  (cf_out, NF90_WRITE,   idf_out)

     !! Need IDs of variables to append... NF90_INQ_VARID
     ierr = NF90_INQ_VARID(idf_out, 'moc_'//trim(cbasin), id_moc)

     !! Need IDs of variables to append... NF90_INQ_VARID
     ierr = NF90_INQ_VARID(idf_out, 'lat_max_'//trim(cbasin), id_lat_max)

     !! Need IDs of variables to append... NF90_INQ_VARID
     ierr = NF90_INQ_VARID(idf_out, 'depthw_max_'//trim(cbasin), id_depthw_max)

     ! Get ID of unlimited dimension
     ierr = NF90_INQUIRE(idf_out, unlimitedDimId = idv_time)

     ! Need to know jt_pos, record number of the last time record writen in the file
     ierr = NF90_INQUIRE_DIMENSION(idf_out, idv_time, name=cv_dum, len=jt_pos)

  END IF

  WRITE(*,'("Going to write record ",i4.4," to ",i4.4," into ",a)') jt_pos+1, jt_pos+1+nt, trim(cf_out)

  DO jt = 1, nt

     ! Writing record jt for time vector and 1d fields:
     ierr = NF90_PUT_VAR( idf_out, idv_time, (/ryear+1./12.*(REAL(jt)-1.+0.5)/), start=(/jt_pos+jt/), count=(/1/) )

     !! Default variable is the only one present (index = 1) :
     ierr = NF90_PUT_VAR(idf_out, id_moc, (/ ovtmax(jt) /), start=(/jt_pos+jt/), count=(/1/))
     ierr = NF90_PUT_VAR(idf_out, id_lat_max, (/ rlat(1,jlatmax(jt)) /), start=(/jt_pos+jt/), count=(/1/))
     ierr = NF90_PUT_VAR(idf_out, id_depthw_max, (/ gdepw(kdmax(jt)) /), start=(/jt_pos+jt/), count=(/1/))

  END DO

  ierr = NF90_CLOSE(idf_out)

  PRINT *, ''

  !! Writing into ASCCI file:

  !OPEN(15, FILE='max_moc.dat', FORM='FORMATTED', RECL=512, STATUS='unknown')
  !WRITE(15,*) '#    Time    Max MOC (Sv)    Lat       Depth' !     ('//trim(cbasin)//')'
  !DO jt = 1, nt
  !   rdt = 1./REAL(nt)
  !   rmon = ryear + rdt/2. + (jt-1)*rdt
  !   WRITE(15,'(f10.4,"  ",f10.4,"  ",f10.4,"  ",f10.4)') rmon, ovtmax(jt), rlat(1,jlatmax(jt)), gdepw(kdmax(jt))
  !END DO
  !CLOSE(15)


END PROGRAM cdfmaxmoc
