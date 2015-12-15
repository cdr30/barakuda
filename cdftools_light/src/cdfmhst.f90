PROGRAM cdfmhst

  !!--------------------------------------------------------------------
  !!               ***  PROGRAM  cdfmhst  ***
  !!
  !!  **  Purpose  : Compute Meridional Heat Salt  Transport.
  !!
  !!  **  Method   : Starts from the mean VT, VS fields computed by cdfvT
  !!                 The program looks for the file "new_maskglo.nc". If it does not exist,
  !!                 only the calculation over all the domain is performed (this is adequate
  !!                 for a basin configuration like NATL4).
  !!
  !!
  !! history :
  !!      Original : J.M. Molines (jan. 2005)
  !!                 J.M. Molines apr. 2005 : use modules
  !!                 A.M. Treguier (april 2006) adaptation to NATL4 case
  !!                 J.M. Molines ( April 2007) : add netcdf output
  !!--------------------------------------------------------------------
  !!  $Rev: 298 $
  !!  $Date: 2010-04-14 18:09:06 +0200 (Wed, 14 Apr 2010) $
  !!  $Id: cdfmhst.f90 298 2010-04-14 16:09:06Z dussin9r $
  !!--------------------------------------------------------------

  USE netcdf

  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE

  LOGICAL, PARAMETER :: l_indo_pac = .FALSE.

  INTEGER   :: jj,jk,jt                         !: dummy loop index
  INTEGER   :: narg, iargc                      !: command line
  INTEGER   :: npiglo, npjglo, npk, nt          !: size of the domain
  INTEGER   :: numout = 10
  INTEGER, DIMENSION(2)          ::  iloc
  LOGICAL    :: llglo = .false.                !: indicator for presence of new_maskglo.nc file

  LOGICAL :: lfncout = .false.

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e3v ,gphiv, zvt, zvs !: mask, metrics

  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE ::  zmask

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat              !: latitude for i = north pole

  REAL(KIND=4), DIMENSION (:),       ALLOCATABLE ::  tim       !LB

  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE :: zwkt, zwks
  REAL(KIND=8), DIMENSION (:,:,:), ALLOCATABLE :: ztrpt, ztrps !LB

  REAL(KIND=8) ,DIMENSION(:,:) , ALLOCATABLE ::  vmerid_heat, vmerid_salt
  REAL(KIND=8) ,DIMENSION(:)   , ALLOCATABLE ::  vback_1d

  REAL(KIND=8) ,DIMENSION(:,:) , ALLOCATABLE ::  merid_heat_glo, merid_heat_atl, merid_heat_pac,&
       &                                       merid_heat_ind, merid_heat_inp
  REAL(KIND=8) ,DIMENSION(:,:) , ALLOCATABLE ::  merid_salt_glo, merid_salt_atl, merid_salt_pac,&
       &                                       merid_salt_ind, merid_salt_inp, zmtrp

  CHARACTER(LEN=256) :: cfileVT

  CHARACTER(LEN=256), PARAMETER :: cfileout='merid_heat_trp.dat', cfileouts='merid_salt_trp.dat'

  ! to be put in namelist eventually
  CHARACTER(LEN=256) :: cf_mm='mesh_mask.nc', cf_bm='new_maskglo.nc'

  ! NC output
  INTEGER            :: jbasin, js, jvar   !: dummy loop index
  INTEGER            :: nbasins, ierr

  REAL(KIND=4), PARAMETER :: rpspval=9999.99
  REAL(KIND=4), DIMENSION(1) :: gdep

  CHARACTER(LEN=256) :: cdum, ctim
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE   :: vcv_name_t, vcv_name_s   !: array of var name for input
  CHARACTER(LEN=4), DIMENSION(5) :: cbasin= (/'_glo','_atl','_pac','_ind','_inp'/)

  ! constants
  REAL(KIND=4),PARAMETER   ::  rau0=1000.,   rcp=4000.


  CHARACTER(LEN=64) :: cv_out = 'mht', cv_dum
  CHARACTER(LEN=1024) :: cf_out = 'mht_1d.nc'
  INTEGER :: idf_out, idd_y, idv_lat, idd_t, idv_time, jt_pos

  INTEGER, DIMENSION(5) :: vid_heat, vid_salt
  REAL :: ryear
  REAL :: minus_one = -1.0


  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg /= 3 ) THEN
     PRINT *,' Usage : cdfmhst  <VTfile> <file_output.nc> <year>'
     PRINT *,' Files mesh_mask.nc and new_maskglo.nc must be in te current directory'
     PRINT *,' NetCDF Output in <file_output.nc> with variables :'
     PRINT *,'                       zomht_glo, zomht_atl, zomht_inp, zomht_pac'
     PRINT *,'                       zomst_glo, zomst_atl, zomst_inp, zomst_pac'
     STOP
  ENDIF

  CALL getarg(1, cfileVT)
  CALL getarg(2, cf_out)
  CALL getarg(3,cdum)    ; READ(cdum,*) ryear




  npiglo= getdim (cfileVT,'x')
  npjglo= getdim (cfileVT,'y')
  npk   = getdim (cfileVT,'depth')

  ctim = 'none'
  nt    = getdim (cfileVT,'time', cdtrue=ctim) !LB

  !LB:
  IF (nt == 0) THEN
     PRINT *, 'nt=0, assume 1' ; nt = 1
  END IF
  !LB.



  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'nt    =', nt

  ! Allocate arrays
  ALLOCATE ( zwkt(npiglo,npjglo), zmask(npiglo,npjglo,5), zvt(npiglo,npjglo) )
  ALLOCATE ( zwks(npiglo,npjglo) ,zvs(npiglo,npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo), gphiv(npiglo,npjglo) )
  ALLOCATE ( ztrpt(npiglo,npjglo,nt) , ztrps(npiglo,npjglo,nt) )
  ALLOCATE ( vmerid_heat(npjglo,5), vmerid_salt(npjglo,5), vback_1d(npjglo) )
  ALLOCATE ( merid_heat_glo(npjglo,nt), merid_heat_atl(npjglo,nt), merid_heat_pac(npjglo,nt) )
  ALLOCATE ( merid_heat_ind(npjglo,nt), merid_heat_inp(npjglo,nt) )
  ALLOCATE ( merid_salt_glo(npjglo,nt), merid_salt_atl(npjglo,nt), merid_salt_pac(npjglo,nt) )
  ALLOCATE ( merid_salt_ind(npjglo,nt), merid_salt_inp(npjglo,nt) )
  ALLOCATE ( zmtrp(npjglo,nt) )
  ALLOCATE ( dumlon(1,npjglo) , dumlat(1,npjglo))
  ALLOCATE ( tim(nt) )




  ! create output fileset
  e1v(:,:)   = getvar(cf_mm, 'e1v', 1,npiglo,npjglo)
  gphiv(:,:) = getvar(cf_mm, 'gphiv', 1,npiglo,npjglo)
  gdep(:)    = getvare3(cf_mm, 'nav_lev' ,1)

  iloc=maxloc(gphiv)
  dumlat(1,:) = gphiv(iloc(1),:)
  dumlon(:,:) = 0.   ! set the dummy longitude to 0



  ztrpt(:,:,:)= 0
  ztrps(:,:,:)= 0


  DO jt = 1, nt !LB

     DO jk = 1,npk

        ! Get temperature and salinity at jk
        zvt(:,:) = getvar(cfileVT, 'vomevt', jk, npiglo,npjglo, ktime=jt)  !LB
        zvs(:,:) = getvar(cfileVT, 'vomevs', jk, npiglo,npjglo, ktime=jt)  !LB

        ! get e3v at level jk
        e3v(:,:)  = getvar(cf_mm, 'e3v_0', jk,npiglo,npjglo, ldiom=.true.)
        zwkt(:,:) = zvt(:,:)*e1v(:,:)*e3v(:,:)
        zwks(:,:) = zvs(:,:)*e1v(:,:)*e3v(:,:)

        ! integrates vertically
        ztrpt(:,:,jt) = ztrpt(:,:,jt) + zwkt(:,:)*rau0*rcp
        ztrps(:,:,jt) = ztrps(:,:,jt) + zwks(:,:)

     END DO  ! jk


     ! ~~~~~~
     ! global
     ! ~~~~~~~
     IF ( jt == 1 ) zmask(:,:,1) = getvar(cf_mm, 'vmask', 1, npiglo, npjglo) !LB

     DO jj=1,npjglo
        merid_heat_glo(jj,jt) = SUM( ztrpt(2:npiglo-1,jj,jt)*zmask(2:npiglo-1,jj,1) )
        merid_salt_glo(jj,jt) = SUM( ztrps(2:npiglo-1,jj,jt)*zmask(2:npiglo-1,jj,1) )
     END DO


  END DO ! jt




  !  Detects newmaskglo file
  INQUIRE( FILE=cf_bm, EXIST=llglo )


  nbasins=1
  IF ( llglo) THEN ! 5 basins
     nbasins=5
  ENDIF


  ! Allocate output variables
  ALLOCATE(vcv_name_t(nbasins), vcv_name_s(nbasins))



  PRINT *, ''

  IF ( llglo ) THEN

     DO jt = 1, nt !LB

        ! Merid mean with mask

        ! Atlantic
        IF ( jt == 1 ) zmask(:,:,2) = getvar(cf_bm,'tmaskatl',1,npiglo,npjglo)

        DO jj=1,npjglo
           merid_heat_atl(jj,jt) = SUM( ztrpt(:,jj,jt)*zmask(:,jj,2) )
           merid_salt_atl(jj,jt) = SUM( ztrps(:,jj,jt)*zmask(:,jj,2) )
        END DO


        ! Pacific
        IF ( jt == 1 ) zmask(:,:,3) = getvar(cf_bm,'tmaskpac',1,npiglo,npjglo)

        DO jj=1,npjglo
           merid_heat_pac(jj,jt)= SUM( ztrpt(:,jj,jt)*zmask(:,jj,3) )
           merid_salt_pac(jj,jt)= SUM( ztrps(:,jj,jt)*zmask(:,jj,3) )
        END DO


        ! Indian
        IF ( jt == 1 ) zmask(:,:,4) = getvar(cf_bm,'tmaskind',1,npiglo,npjglo)

        DO jj=1,npjglo
           merid_heat_ind(jj,jt)= SUM( ztrpt(:,jj,jt)*zmask(:,jj,4) )
           merid_salt_ind(jj,jt)= SUM( ztrps(:,jj,jt)*zmask(:,jj,4) )
        END DO


        IF ( l_indo_pac ) THEN
           ! Indo-Pacific

           IF ( jt == 1 ) THEN
              !PRINT *, ''; PRINT *, 'cdfmhst.f90: opening "tmaskinp" into ', trim(cf_bm)
              zmask(:,:,5) = getvar(cf_bm,'tmaskinp',1,npiglo,npjglo)
           END IF

           DO jj=1,npjglo
              merid_heat_inp(jj,jt)= SUM( ztrpt(:,jj,jt)*zmask(:,jj,5) )
              merid_salt_inp(jj,jt)= SUM( ztrps(:,jj,jt)*zmask(:,jj,5) )
           END DO
           !PRINT *, ''

           !zmask_inp(:,:) = 0.
           !zmask(:,:) = getvar(cf_bm,'tmaskind',1,npiglo,npjglo)
           !WHERE( zmask == 1. ) zmask_inp = 1.
           !zmask(:,:) = getvar(cf_bm,'tmaskpac',1,npiglo,npjglo)
           !WHERE( zmask == 1. ) zmask_inp = 1.
           !DO jj=1,npjglo
           !   merid_heat_inp(jj,jt)= SUM( ztrpt(:,jj,jt)*zmask_inp(:,jj))
           !   merid_salt_inp(jj,jt)= SUM( ztrps(:,jj,jt)*zmask_inp(:,jj))
           !END DO

        END IF


     END DO

  END IF



  tim  = getvar1d(cfileVT, trim(ctim), nt) !LB: nt

  PRINT *, ''

!  DO jt = 1, nt !LB
!
!     PRINT *, ' *** jt =', jt
!
!     ! MHT !
!
!     js=1
!
!     zmtrp(:,jt)=merid_heat_glo(:,jt)/1.e15                        ! GLO
!     WHERE ( zmtrp(:,jt) == 0 ) zmtrp(:,jt) = rpspval
!     js=js+1
!
!     IF ( nbasins == 5 ) THEN
!
!        zmtrp(:,jt)=merid_heat_atl(:,jt)/1.e15                      ! ATL
!        WHERE ( zmtrp(:,jt) == 0 ) zmtrp(:,jt)=rpspval
!        js=js+1
!
!        zmtrp(:,jt)=merid_heat_pac(:,jt)/1.e15                      ! PAC
!        WHERE ( zmtrp(:,jt) == 0 ) zmtrp(:,jt)=rpspval
!        js=js+1
!
!        zmtrp(:,jt)=merid_heat_ind(:,jt)/1.e15                      ! IND
!        WHERE ( zmtrp(:,jt) == 0 ) zmtrp(:,jt)=rpspval
!        js=js+1
!
!        IF ( l_indo_pac ) THEN
!           zmtrp(:,jt)=merid_heat_inp(:,jt)/1.e15                      ! INP
!           WHERE ( zmtrp(:,jt) == 0 ) zmtrp(:,jt)=rpspval
!        END IF
!        js=js+1
!
!     ENDIF
!
!
!!!!!!!!
!     ! MST !
!!!!!!!!
!
!     zmtrp(:,jt)=merid_salt_glo(:,jt)/1.e6                        ! GLO
!     WHERE ( zmtrp(:,jt) == 0 ) zmtrp(:,jt)=rpspval
!     js = js + 1
!
!     IF ( nbasins == 5 ) THEN
!
!        zmtrp(:,jt)=merid_salt_atl(:,jt)/1.e6                      ! ATL
!        WHERE ( zmtrp(:,jt) == 0 ) zmtrp(:,jt)=rpspval
!        js = js + 1
!
!        zmtrp(:,jt)=merid_salt_pac(:,jt)/1.e6                      ! PAC
!        WHERE ( zmtrp(:,jt) == 0 ) zmtrp(:,jt)=rpspval
!        js = js + 1
!
!        zmtrp(:,jt)=merid_salt_ind(:,jt)/1.e6                      ! IND
!        WHERE ( zmtrp(:,jt) == 0 ) zmtrp(:,jt)=rpspval
!        js = js + 1
!
!        zmtrp(:,jt)=merid_salt_inp(:,jt)/1.e6                      ! INP
!        WHERE ( zmtrp(:,jt) == 0 ) zmtrp(:,jt)=rpspval
!        js = js + 1
!
!     ENDIF
!
!  END DO !jt




  !! Plotting annually-averaged profiles:  !LB
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  DO jj = 1, npjglo

     vmerid_heat(jj,1) = SUM(merid_heat_glo(jj,:))/nt/1e15
     vmerid_heat(jj,2) = SUM(merid_heat_atl(jj,:))/nt/1e15
     vmerid_heat(jj,3) = SUM(merid_heat_pac(jj,:))/nt/1e15
     vmerid_heat(jj,4) = SUM(merid_heat_ind(jj,:))/nt/1e15
     IF ( l_indo_pac ) vmerid_heat(jj,5) = SUM(merid_heat_inp(jj,:))/nt/1e15

     vmerid_salt(jj,1) = SUM(merid_salt_glo(jj,:))/nt/1e6
     vmerid_salt(jj,2) = SUM(merid_salt_atl(jj,:))/nt/1e6
     vmerid_salt(jj,3) = SUM(merid_salt_pac(jj,:))/nt/1e6
     vmerid_salt(jj,4) = SUM(merid_salt_ind(jj,:))/nt/1e6
     IF ( l_indo_pac ) vmerid_salt(jj,5) = SUM(merid_salt_inp(jj,:))/nt/1e6

     IF ( (vmerid_salt(jj,2) == 0.0).AND.(jj > 20) ) THEN
        vmerid_salt(jj,4) = 0.0
        vmerid_salt(jj,5) = 0.0
     END IF
     IF ( (vmerid_heat(jj,2) == 0.0).AND.(jj > 20) ) THEN
        vmerid_heat(jj,4) = 0.0
        vmerid_heat(jj,5) = 0.0
     END IF

  END DO


  ! Smoothing Pacific, Indian and Indo-Pacific:
  !DO jbasin = 3,5
  !   vback_1d(:) = vmerid_heat(:,jbasin)
  !   DO jj = 3, npjglo-2
  !      vmerid_heat(jj,jbasin) = (vback_1d(jj) + vback_1d(jj-1) + vback_1d(jj+1) + vback_1d(jj-2) + vback_1d(jj+2))/5.
  !   END DO
  !   vback_1d(:) = vmerid_salt(:,jbasin)
  !   DO jj = 3, npjglo-2
  !      vmerid_salt(jj,jbasin) = (vback_1d(jj) + vback_1d(jj-1) + vback_1d(jj+1) + vback_1d(jj-2) + vback_1d(jj+2))/5.
  !   END DO
  !END DO

  WHERE ( vmerid_heat == 0.0 ) vmerid_heat = sqrt(minus_one)
  WHERE ( vmerid_salt == 0.0 ) vmerid_salt = sqrt(minus_one)






  DO jbasin = 1, nbasins
     vcv_name_t(jbasin) = 'zomht'//TRIM(cbasin(jbasin))
     vcv_name_s(jbasin) = 'zomst'//TRIM(cbasin(jbasin))
     !long_name='Meridional Heat Transport '//TRIM(cbasin(jbasin))
  END DO


  !


  !! LOLO netcdf
  INQUIRE( FILE=cf_out, EXIST=lfncout )


  IF ( .NOT. lfncout ) THEN

     !! Creating file
     PRINT *, ' Creating file '//trim(cf_out)//' !!!'
     ierr = NF90_CREATE(cf_out, NF90_CLOBBER, idf_out)
     ierr = NF90_DEF_DIM(idf_out, 'y',   npjglo, idd_y) !
     ierr = NF90_DEF_VAR(idf_out, 'lat', NF90_FLOAT, idd_y, idv_lat)
     ierr = NF90_DEF_DIM(idf_out, 'time', NF90_UNLIMITED, idd_t)
     ierr = NF90_DEF_VAR(idf_out, 'time', NF90_DOUBLE,    idd_t, idv_time)

     DO jbasin = 1, nbasins
        ierr = NF90_DEF_VAR(idf_out, trim(vcv_name_t(jbasin)), NF90_FLOAT, (/idd_y,idd_t/), vid_heat(jbasin))
        ierr = NF90_DEF_VAR(idf_out, trim(vcv_name_s(jbasin)), NF90_FLOAT, (/idd_y,idd_t/), vid_salt(jbasin))
     END DO

     ierr = NF90_ENDDEF(idf_out)
     ierr = NF90_PUT_VAR(idf_out, idv_lat,  dumlat(1,:))
     jt_pos = 0
  ELSE

     !! Opening already existing file
     ierr = NF90_OPEN  (cf_out, NF90_WRITE,   idf_out)

     !! Need IDs of variables to append... NF90_INQ_VARID
     DO jbasin = 1, nbasins
        ierr = NF90_INQ_VARID(idf_out, trim(vcv_name_t(jbasin)), vid_heat(jbasin))
        ierr = NF90_INQ_VARID(idf_out, trim(vcv_name_s(jbasin)), vid_salt(jbasin))
     END DO


     ! Get ID of unlimited dimension
     ierr = NF90_INQUIRE(idf_out, unlimitedDimId = idv_time)

     ! Need to know jt_pos, record number of the last time record writen in the file
     ierr = NF90_INQUIRE_DIMENSION(idf_out, idv_time, name=cv_dum, len=jt_pos)
     !PRINT *, 'cv_dum, jt_pos = ', trim(cv_dum), jt_pos
  END IF

  WRITE(*,'("Going to write record #",i4.4," into ",a)') jt_pos+1, trim(cf_out)

  !DO jt = 1, nt
  jt = 1 ! saving 1 record per year!

  ! Writing record jt for time vector and 1d fields:
  ierr = NF90_PUT_VAR( idf_out, idv_time, (/ryear+0.5/), start=(/jt_pos+jt/), count=(/1/) )

  DO jbasin = 1, nbasins
     ierr = NF90_PUT_VAR(idf_out, vid_heat(jbasin), vmerid_heat(:,jbasin), start=(/1,jt_pos+jt/), count=(/npjglo,1/))
     ierr = NF90_PUT_VAR(idf_out, vid_salt(jbasin), vmerid_salt(:,jbasin), start=(/1,jt_pos+jt/), count=(/npjglo,1/))
  END DO
  !END DO

  ierr = NF90_PUT_ATT(idf_out, NF90_GLOBAL, 'About', 'Created by BaraKuda (cdfmhst.f90), contact: brodeau@gmail.com')
  ierr = NF90_CLOSE(idf_out)

END PROGRAM cdfmhst
