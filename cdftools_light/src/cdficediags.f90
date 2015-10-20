PROGRAM cdficediag
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdficediag  ***
  !!
  !!  **  Purpose  :  Compute the Ice volume, area and extend for each hemisphere
  !!
  !!  **  Method   :   Read the ice output and integrates (2D)
  !!                   determine the hemisphere by the sign of ff (coriolis)
  !!
  !! history ;
  !!  Original :  J.M. Molines (Jan. 2006)
  !!              R. Dussin (Jul. 2009) : Add netcdf output
  !!-------------------------------------------------------------------
  !!  $Rev: 319 $
  !!  $Date: 2010-05-19 12:19:07 +0200 (Wed, 19 May 2010) $
  !!  $Id: cdficediags.f90 319 2010-05-19 10:19:07Z dussin $
  !!--------------------------------------------------------------
  !! * Modules used
  USE netcdf

  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, jj, jt !LB
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line
  INTEGER   :: npiglo,npjglo, nt !LB               !: size of the domain
  INTEGER   :: nvpk                                !: vertical levels in working variable
  INTEGER   :: nperio = 4                          !: boundary condition ( periodic, north fold)

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  e1, e2            !:  metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask ,ff         !:   npiglo x npjglo
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  ricethick, riceldfra !: thickness, leadfrac (concentration)

  REAL(KIND=8)      :: zvols,  zareas, zextends,zextends2        !: volume, area extend South hemisphere
  REAL(KIND=8)      :: zvoln,  zarean, zextendn,zextendn2        !: volume, area extend North hemisphere

  REAL(KIND=8), DIMENSION (:), ALLOCATABLE  :: vvolu_n, varea_n, vvolu_s, varea_s



  CHARACTER(LEN=3)   :: conly_c='000' !lolo
  CHARACTER(LEN=256) :: cf_ice , cdum, ctim
  CHARACTER(LEN=256) :: cf_mm='mesh_mask.nc'

  INTEGER    :: istatus

  !! LOLO:
  LOGICAL :: lfncout = .false.
  CHARACTER(LEN=256) :: cd_out = '.', cf_out
  CHARACTER(LEN=64)  :: cv_u, cv_v, cv_ueiv, cv_veiv, cv_dum
  REAL :: ryear
  INTEGER :: jt_pos, idf_out, idd_t, idv_time
  INTEGER :: id_volu_n, id_area_n, id_volu_s, id_area_s
  !! LOLO.



  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg >= 3) THEN
     CALL getarg (1, cf_ice)
     npiglo = getdim (cf_ice,'x')
     npjglo = getdim (cf_ice,'y')
     ctim = 'none'
     nt    = getdim (cf_ice,'time',cdtrue=ctim,kstatus=istatus) !LB
     CALL getarg (2, cdum)  ; READ(cdum,*) ryear
     CALL getarg (3, cd_out)
     IF ( narg == 4) CALL getarg (4, conly_c)
  ELSE
     PRINT *,' Usage : cdficediag.x <ice_file> <year> <DIROUT>'
     PRINT *,' File mesh_mask.nc must be in the current directory'
     STOP
  ENDIF

  ALLOCATE ( zmask(npiglo,npjglo) ,ff(npiglo,npjglo) )
  ALLOCATE ( ricethick(npiglo,npjglo) )
  ALLOCATE ( riceldfra(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo),e2(npiglo,npjglo) )
  !ALLOCATE ( tim(nt) ) !LB
  !tim = getvar1d(cf_ice,trim(ctim),nt) !LB


  e1(:,:) = getvar(cf_mm, 'e1t', 1,npiglo,npjglo)
  e2(:,:) = getvar(cf_mm, 'e2t', 1,npiglo,npjglo)
  ! only the sign of ff is important
  ff(:,:) = getvar(cf_mm, 'gphit' , 1,npiglo,npjglo)


  ! modify the mask for periodic and north fold condition (T pivot, F Pivot ...)
  ! in fact should be nice to use jperio as in the code ...

  zmask(:,:)=getvar(cf_mm,'tmask',1,npiglo,npjglo)
  SELECT CASE (nperio)
  CASE (0) ! closed boundaries
     ! nothing to do
  CASE (4) ! ORCA025 type boundary
     zmask(1:2,:)=0.
     zmask(:,npjglo)=0.
     zmask(npiglo/2+1:npiglo,npjglo-1)= 0.
  CASE (6)
     zmask(1:2,:)=0.
     zmask(:,npjglo)=0.
  CASE DEFAULT
     PRINT *,' Nperio=', nperio,' not yet coded'
     STOP
  END SELECT




  ALLOCATE ( vvolu_n(nt), varea_n(nt), vvolu_s(nt), varea_s(nt) )


  DO jt = 1, nt

     PRINT *, ' jt = ', jt
     
     IF ( conly_c == '000' ) THEN
        ricethick(:,:)= getvar(cf_ice, 'ice_thic',  1 ,npiglo,npjglo, ktime=jt)
     END IF
     riceldfra(:,:)= getvar(cf_ice, 'ice_frac',  1 ,npiglo,npjglo, ktime=jt)
     
     !LB:
     WHERE ( zmask == 0. )
        ricethick = 0.
        riceldfra = 0.
     END WHERE
     !LB.

     ! North : ff > 0
     zvoln=SUM( ricethick (:,:)* e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:) , (ff > 0 ) )

     zarean=SUM( e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:) ,( ff > 0 ) )
     zextendn=SUM( e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:), (riceldfra > 0.15 .AND. ff > 0 ) )
     ! JMM added 22/01/2007 : to compute same extent than the NSIDC
     zextendn2=SUM( e1(:,:) * e2(:,:) * zmask (:,:), (riceldfra > 0.15 .AND. ff > 0 ) )

     ! South : ff < 0
     zvols=SUM( ricethick (:,:)* e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:) ,(ff < 0 ) )

     zareas=SUM( e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:), ( ff < 0 ) )
     zextends=SUM( e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:), (riceldfra > 0.15 .AND. ff < 0  ) )
     zextends2=SUM( e1(:,:) * e2(:,:)* zmask (:,:), (riceldfra > 0.15 .AND. ff < 0  ) )

     vvolu_n(jt) =  zvoln/1.d12
     varea_n(jt) =  zarean/1.d12
     vvolu_s(jt) =  zvols/1.d12
     varea_s(jt) =  zareas/1.d12


     !LB:
     !IF ( jt == 1 ) THEN
     !   OPEN(UNIT=12, FILE='ice.dat', FORM='FORMATTED', RECL=256, STATUS='unknown')
     !   WRITE(12,*) ' #    month  N.Volume (10^3km^3)      N.Area (10^6km^2)     S.Volume (10^3km^3)    S.Area (10^6km^2)'
     !END IF
     !
     !WRITE(12,*) jt, zvoln /1.d12, zarean /1.d12, zvols /1.d12, zareas /1.d12


     !PRINT *,' Northern Hemisphere '
     !PRINT *,'          NVolume (10^9 m3)  ', zvoln /1.d9
     !PRINT *,'          NArea (10^9 m2)    ', zarean /1.d9
     !PRINT *,'          NExtend (10^9 m2)  ', zextendn /1.d9
     !PRINT *,'          NExnsidc (10^9 m2)  ', zextendn2 /1.d9
     !PRINT *
     !PRINT *,' Southern Hemisphere '
     !PRINT *,'          SVolume (10^9 m3)  ', zvols /1.d9
     !PRINT *,'          SArea (10^9 m2)    ', zareas /1.d9
     !PRINT *,'          SExtend (10^9 m2)  ', zextends /1.d9
     !PRINT *,'          SExnsidc (10^9 m2)  ', zextends2 /1.d9


  END DO ! DO jt = 1, nt




  !! Time to create or append X_trsp into the netcdf file for current section
  !! -----------------------------------------------------

  WRITE(cf_out, '(a,"/seaice_diags.nc")') trim(cd_out)

  id_volu_n = 0 ; id_area_n = 0 ; id_volu_s = 0
  
  !! LOLO netcdf
  INQUIRE( FILE=cf_out, EXIST=lfncout )
  
  IF ( .NOT. lfncout ) THEN
     
     !! Creating file
     PRINT *, ' Creating file '//trim(cf_out)//' !!!'
     ierr = NF90_CREATE(cf_out, NF90_CLOBBER, idf_out)
     ierr = NF90_DEF_DIM(idf_out, 'time', NF90_UNLIMITED, idd_t)
     ierr = NF90_DEF_VAR(idf_out, 'time', NF90_DOUBLE,    idd_t, idv_time)
     
     ierr = NF90_DEF_VAR(idf_out, 'volu_ne', NF90_FLOAT, (/idd_t/), id_volu_n)
     ierr = NF90_DEF_VAR(idf_out, 'area_ne', NF90_FLOAT, (/idd_t/), id_area_n)
     ierr = NF90_DEF_VAR(idf_out, 'volu_se', NF90_FLOAT, (/idd_t/), id_volu_s)
     ierr = NF90_DEF_VAR(idf_out, 'area_se', NF90_FLOAT, (/idd_t/), id_area_s)
     
     ierr = NF90_PUT_ATT(idf_out, id_volu_n, 'long_name', 'Total volume of sea-ice in Northern Hemisphere')
     ierr = NF90_PUT_ATT(idf_out, id_area_n, 'long_name', 'Total area of sea-ice in Northern Hemisphere')
     ierr = NF90_PUT_ATT(idf_out, id_volu_s, 'long_name', 'Total volume of sea-ice in Southern Hemisphere')
     ierr = NF90_PUT_ATT(idf_out, id_area_s, 'long_name', 'Total area of sea-ice in Southern Hemisphere')
     
     ierr = NF90_PUT_ATT(idf_out, id_volu_n, 'units', '10^3km^3')
     ierr = NF90_PUT_ATT(idf_out, id_area_n, 'units', '10^6km^2')
     ierr = NF90_PUT_ATT(idf_out, id_volu_s, 'units', '10^3km^3')
     ierr = NF90_PUT_ATT(idf_out, id_area_s, 'units', '10^6km^2')
     
     !ierr = NF90_PUT_ATT(idf_out, NF90_GLOBAL, 'Info', 'Reference temperature for heat transport is '//trim(crt)//' deg.C')
     ierr = NF90_PUT_ATT(idf_out, NF90_GLOBAL, 'About', 'Created by BaraKuda (cdficediags.f90), contact: brodeau@gmail.com')
     
     ierr = NF90_ENDDEF(idf_out)
     jt_pos = 0
     
  ELSE
     
     !! Opening already existing file
     ierr = NF90_OPEN  (cf_out, NF90_WRITE,   idf_out)

     !! Need IDs of variables to append... NF90_INQ_VARID
     ierr = NF90_INQ_VARID(idf_out, 'volu_ne', id_volu_n)
     ierr = NF90_INQ_VARID(idf_out, 'area_ne', id_area_n)
     ierr = NF90_INQ_VARID(idf_out, 'volu_se', id_volu_s)
     ierr = NF90_INQ_VARID(idf_out, 'area_se', id_area_s)
     
     ! Get ID of unlimited dimension
     ierr = NF90_INQUIRE(idf_out, unlimitedDimId = idv_time)
     
     ! Need to know jt_pos, record number of the last time record writen in the file
     ierr = NF90_INQUIRE_DIMENSION(idf_out, idv_time, name=cv_dum, len=jt_pos)
     
  END IF
  
  WRITE(*,'("Going to write record ",i4.4," to ",i4.4," into ",a)') jt_pos+1, jt_pos+1+nt, trim(cf_out)

  DO jt = 1, nt

     ! Writing record jt for time vector and 1d fields:
     ierr = NF90_PUT_VAR( idf_out, idv_time, (/ryear+1./12.*(REAL(jt)-1.+0.5)/), start=(/jt_pos+jt/), count=(/1/) )
     
     ierr = NF90_PUT_VAR(idf_out, id_volu_n, (/ vvolu_n(jt) /), start=(/jt_pos+jt/), count=(/1/))
     ierr = NF90_PUT_VAR(idf_out, id_area_n, (/ varea_n(jt) /), start=(/jt_pos+jt/), count=(/1/))
     ierr = NF90_PUT_VAR(idf_out, id_volu_s, (/ vvolu_s(jt) /), start=(/jt_pos+jt/), count=(/1/))
     ierr = NF90_PUT_VAR(idf_out, id_area_s, (/ varea_s(jt) /), start=(/jt_pos+jt/), count=(/1/))
     
  END DO
  
  ierr = NF90_CLOSE(idf_out)
  
END PROGRAM cdficediag
