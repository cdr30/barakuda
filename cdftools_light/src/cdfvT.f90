PROGRAM cdfvT


  !!-------------------------------------------------------------------
  !!                ***  PROGRAM cdfvT  ***
  !!
  !!  **  Purpose:
  !!
  !!  **  Method: Try to avoid 3 d arrays
  !!              Assume that all input files have the same number of time frames
  !!
  !! history :
  !!   L. Brodeau , 2014 for BaraKuda !!
  !!
  !!   Original : J.M. Molines (Nov 2004 ) for ORCA025
  !!              J.M. Molines (apr 2005 ) : use of modules
  !!              J.M. Molines (Feb. 2010 ): handle multiframes input files.
  !!-------------------------------------------------------------------
  !!  $Rev: 317 $
  !!  $Date: 2010-05-17 14:47:12 +0200 (Mon, 17 May 2010) $
  !!  $Id: cdfvT.f90 317 2010-05-17 12:47:12Z molines $
  !!--------------------------------------------------------------



  USE cdfio
  USE io_ezcdf


  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk,jt,jkk                !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc
  INTEGER   :: npiglo,npjglo, npk, nt              !: size of the domain

  INTEGER, DIMENSION(4) ::  ipk, id_varout


  REAL(KIND=4) , DIMENSION (:),   ALLOCATABLE :: vtime ! lolo
  REAL(KIND=4) , DIMENSION (:),   ALLOCATABLE :: vdepth ! lolo
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: xlon, xlat  ! lolo
  REAL(KIND=4) , DIMENSION (:,:,:), ALLOCATABLE :: T_3D, S_3D, U_3D, V_3D, UEIV_3D, VEIV_3D, X_3D_u, X_3D_v   !: lolo

  REAL(KIND=8) , DIMENSION (:,:,:), ALLOCATABLE :: zcumulut, zcumulus  !: Arrays for cumulated values
  REAL(KIND=8) , DIMENSION (:,:,:), ALLOCATABLE :: zcumulvt, zcumulvs  !: Arrays for cumulated values


  CHARACTER(LEN=256) :: cf_t,cf_u,cf_v , cf_out, conf_tag , ctim !:
  TYPE (variable), DIMENSION(4)     :: typvar     !: structure for attributes
  LOGICAL :: lexist                               !: to inquire existence of files

  INTEGER    :: istatus

  CHARACTER(LEN=64) :: cv_t, cv_s, cv_u, cv_v, cv_ueiv, cv_veiv

  INTEGER :: idf_t=0, idv_t=0, idf_s=0, idv_s=0, &
       &     idf_u=0, idv_u=0, idf_v=0, idv_v=0, &
       &     idf_ueiv=0, idv_ueiv=0, idf_veiv=0, idv_veiv=0

  INTEGER :: idf_vt=0, idv_vt=0, idv_vs=0, idv_ut=0, idv_us=0, idf_0=0, idv_0=0



  CHARACTER(LEN=100) :: cv_depth = 'deptht'

  LOGICAL :: leiv = .FALSE.


  !!  Read command line
  narg= iargc()
  IF ( (narg < 5).OR.(narg > 7).OR.(narg == 6) ) THEN
     PRINT *,' Usage : cdftransportiz <CONF_TAG> <name T> <name S> <name U> <name V> (<name Ueiv> <name Veiv>)'
     PRINT *,'    => files are: <CONF_TAG>_grid_T.nc <CONF_TAG>_grid_U.nc <CONF_TAG>_grid_V.nc'     
     PRINT *,' Files mesh_mask.nc must be in te current directory'
     STOP
  ENDIF

  CALL getarg (1, conf_tag)
  CALL getarg (2, cv_t)
  CALL getarg (3, cv_s)
  CALL getarg (4, cv_u)
  CALL getarg (5, cv_v)

  PRINT *, ' Will compute VT using '//trim(cv_u)//' and '//trim(cv_v)

  IF (narg == 7) THEN
     leiv = .TRUE.     
     CALL getarg (6, cv_ueiv)
     CALL getarg (7, cv_veiv)
     IF ( (trim(cv_ueiv) == '0').AND.(trim(cv_veiv) == '0') )   leiv = .FALSE.
  END IF
  
  
  IF ( leiv) &
       & PRINT *, ' and taking eddy-induced velocity into account using '//trim(cv_ueiv)//' and '//trim(cv_veiv)

  !! Initialisation from 1st file (all file are assume to have the same geometry)


  WRITE(cf_out,'(a,"_VT.nc")') trim(conf_tag)
  
  WRITE(cf_t,'(a,"_grid_T.nc")') trim(conf_tag)
  INQUIRE(FILE=cf_t,EXIST=lexist)
  IF ( .NOT. lexist ) THEN
     WRITE(cf_t,'(a,"_grid_T.nc4")') trim(conf_tag)
     INQUIRE(FILE=cf_t,EXIST=lexist)
     IF ( .NOT. lexist ) THEN
        PRINT *,' ERROR : missing grid_T or even gridT file '
        STOP
     ENDIF
  ENDIF

  PRINT *,TRIM(cf_t)
  npiglo= getdim (cf_t,'x')
  npjglo= getdim (cf_t,'y')
  npk   = getdim (cf_t,'depth')

  ctim = 'none'
  nt    = getdim (cf_t,'time',cdtrue=ctim,kstatus=istatus) !LB


  !LB:
  IF (nt == 0) THEN
     PRINT *, 'nt=0, assume 1' ; nt = 1
  END IF
  !LB.


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk


  !! LOLO BIG:
  ALLOCATE( vtime(nt), vdepth(npk), xlon(npiglo,npjglo), xlat(npiglo,npjglo) )

  ALLOCATE( T_3D(npiglo,npjglo,npk), S_3D(npiglo,npjglo,npk), U_3D(npiglo,npjglo,npk), V_3D(npiglo,npjglo,npk) )
  ALLOCATE( X_3D_u(npiglo,npjglo,npk), X_3D_v(npiglo,npjglo,npk) )

  IF ( leiv ) ALLOCATE( UEIV_3D(npiglo,npjglo,npk), VEIV_3D(npiglo,npjglo,npk) )

  ALLOCATE( zcumulut(npiglo,npjglo,npk), zcumulus(npiglo,npjglo,npk) )
  ALLOCATE( zcumulvt(npiglo,npjglo,npk), zcumulvs(npiglo,npjglo,npk) )



  
  vtime  = getvar1d(cf_t, trim(ctim), nt) !LB
  vdepth = getvar1d(cf_t, trim(cv_depth), npk) !LB
  

  CALL GETVAR_2D(idf_0, idv_0, cf_t, 'nav_lon', 0, 0, 0, xlon)
  idf_0=0 ;  idv_0=0
  CALL GETVAR_2D(idf_0, idv_0, cf_t, 'nav_lat', 0, 0, 0, xlat)




  WRITE(cf_t,'(a,"_grid_T.nc")') trim(conf_tag)
  INQUIRE(FILE=cf_t,EXIST=lexist)
  IF ( .NOT. lexist ) THEN
     WRITE(cf_t,'(a,"_grid_T.nc4")') trim(conf_tag)
     INQUIRE(FILE=cf_t,EXIST=lexist)
     IF ( .NOT. lexist ) THEN
        PRINT *,' ERROR : missing gridT or even grid_T file '
        STOP
     ENDIF
  ENDIF

  ! assume U and V file have same time span ...
  WRITE(cf_u,'(a,"_grid_U.nc")') trim(conf_tag)
  INQUIRE(FILE=cf_u,EXIST=lexist)
  IF ( .NOT. lexist ) THEN
     WRITE(cf_u,'(a,"_grid_U.nc4")') trim(conf_tag)
     INQUIRE(FILE=cf_u,EXIST=lexist)
     IF ( .NOT. lexist ) THEN
        PRINT *,' ERROR : missing grid_U or even gridU file '
        STOP
     ENDIF
  ENDIF

  WRITE(cf_v,'(a,"_grid_V.nc")') trim(conf_tag)
  INQUIRE(FILE=cf_v,EXIST=lexist)
  IF ( .NOT. lexist ) THEN
     WRITE(cf_v,'(a,"_grid_V.nc4")') trim(conf_tag)
     INQUIRE(FILE=cf_v,EXIST=lexist)
     IF ( .NOT. lexist ) THEN
        PRINT *,' ERROR : missing grid_V or even gridV file '
        STOP
     ENDIF
  ENDIF






  DO jt=1,nt

     PRINT *, 'jt =', jt

     CALL GETVAR_3D(idf_t, idv_t, cf_t, cv_t, nt, jt, T_3D)
     CALL GETVAR_3D(idf_s, idv_s, cf_t, cv_s, nt, jt, S_3D)
     CALL GETVAR_3D(idf_u, idv_u, cf_u, cv_u, nt, jt, U_3D)
     CALL GETVAR_3D(idf_v, idv_v, cf_v, cv_v, nt, jt, V_3D)
     
     IF ( leiv ) THEN
        CALL GETVAR_3D(idf_ueiv, idv_ueiv, cf_u, cv_ueiv, nt, jt, UEIV_3D)
        CALL GETVAR_3D(idf_veiv, idv_veiv, cf_v, cv_veiv, nt, jt, VEIV_3D)
        U_3D = U_3D + UEIV_3D
        V_3D = V_3D + VEIV_3D
     END IF


     zcumulut(:,:,:) = 0.d0 ;  zcumulvt(:,:,:) = 0.d0 !; total_time = 0.
     zcumulus(:,:,:) = 0.d0 ;  zcumulvs(:,:,:) = 0.d0


     ! temperature
     X_3D_u(1:npiglo-1,:, :) = 0.5 * ( T_3D(1:npiglo-1,:, :) + T_3D(2:npiglo,:, :) )  ! temper at Upoint
     X_3D_v(:,1:npjglo-1, :) = 0.5 * ( T_3D(:,1:npjglo-1, :) + T_3D(:,2:npjglo, :) )  ! temper at Vpoint

     zcumulut(:,:,:) = X_3D_u(:,:,:) * U_3D(:,:,:)
     zcumulvt(:,:,:) = X_3D_v(:,:,:) * V_3D(:,:,:)

     ! salinity
     X_3D_u(1:npiglo-1,:, :) = 0.5 * ( S_3D(1:npiglo-1,:, :) + S_3D(2:npiglo,:, :) )  ! salinity at Upoint
     X_3D_v(:,1:npjglo-1, :) = 0.5 * ( S_3D(:,1:npjglo-1, :) + S_3D(:,2:npjglo, :) )  ! salinity at Vpoint

     zcumulus(:,:,:) = X_3D_u(:,:,:) * U_3D(:,:,:)
     zcumulvs(:,:,:) = X_3D_v(:,:,:) * V_3D(:,:,:)


     !! Printing record jt:
     CALL P3D_T_4v(idf_vt, idv_vt, idv_vs, idv_ut, idv_us, nt, jt, xlon, xlat, vdepth, REAL(vtime,8), &
          &        REAL(zcumulvt,4), REAL(zcumulvs,4), REAL(zcumulut,4), REAL(zcumulus,4),  &
          &        cf_out, 'nav_lon', 'nav_lat', trim(cv_depth), 'time', &
          &        'vomevt', 'vomevs', 'vozout', 'vozous', &
          &        0., 'time', 'm')
     

  END DO ! jt

  PRINT *, '  ***   => '//trim(cf_out)//' written!'; PRINT *, ''

END PROGRAM cdfvT
