PROGRAM cdfmean
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfmean  ***
  !!
  !!  **  Purpose  :  Compute the Mean Value over the ocean
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  compute the sum ( V * e1 *e2 * e3 *mask )/ sum( e1 * e2 * e3 *mask )
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (Oct. 2005) 
  !!              R. Dussin (Jul 2009) : add cdf output
  !!-------------------------------------------------------------------
  !!  $Rev: 319 $
  !!  $Date: 2010-05-19 12:19:07 +0200 (Wed, 19 May 2010) $
  !!  $Id: cdfmean.f90 319 2010-05-19 10:19:07Z dussin $
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, ik, jt, jj
  INTEGER   :: imin=0, imax=0, jmin=0, jmax=0      !: domain limitation for computation
  INTEGER   :: kmin=0, kmax=0                      !: domain limitation for computation
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk, nt              !: size of the domain
  INTEGER   :: nvpk                                !: vertical levels in working variable
  INTEGER   :: numout=10                           !: logical unit for output file
  ! added to write in netcdf
  INTEGER :: kx=1, ky=1                ! dims of netcdf output file
  INTEGER :: nvars=2                ! number of values to write in cdf output
  INTEGER :: ncout, ierr               ! for netcdf output
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout
  !
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  e1, e2, e3,  zv   !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask             !:   npiglo x npjglo
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE ::  gdep              !:  depth 
  ! added to write in netcdf
  REAL(KIND=4) :: threedmeanout, pmissing_value
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  dumlon, dumlat, dummymean
  REAL(KIND=4), DIMENSION (1)               ::  tim ! time counter
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: meanout
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: year_mean
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar  ! structure of output
  !
  REAL(KIND=8)      :: zvol, zsum, zvol2d, zsum2d, zsurf
  CHARACTER(LEN=256) :: cfilev , cdum, cf_out, cf_asout
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc',cmask='mask.nc'
  CHARACTER(LEN=256) :: cvar, cvartype, cdep, ctim
  CHARACTER(LEN=20) :: ce1, ce2, ce3, cvmask, cvtype
  CHARACTER(LEN=256) :: cfilout='out.txt'
  ! added to write in netcdf
  CHARACTER(LEN=256) :: cfileoutnc
  CHARACTER(LEN=256) :: cdunits, cdlong_name, cdshort_name 
  ! added to write in netcdf
  LOGICAL :: lwrtcdf=.TRUE., l_pvp
  !LOGICAL :: lwrtcdf=.TRUE., l_treat_ssh, l_pvp


  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmean  ncfile cdfvar T| U | V | F | W  file_ascii_out  [imin imax jmin jmax kmin kmax] '
     PRINT *,' Computes the mean value of the field (3D, weighted) '
     PRINT *,' imin imax jmin jmax kmin kmax can be given in option '
     PRINT *,'    if imin = 0 then ALL i are taken'
     PRINT *,'    if jmin = 0 then ALL j are taken'
     PRINT *,'    if kmin = 0 then ALL k are taken'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on standard output'
     STOP
  ENDIF
  ! Open standard output with recl=256 to avoid wrapping of long lines (ifort)
  !LB:
  !OPEN(6,FORM='FORMATTED',RECL=256)
  !LB.
  
  CALL getarg (1, cfilev)
  CALL getarg (2, cvar)
  CALL getarg (3, cvartype)
  CALL getarg (4, cf_asout) !LB

  IF (narg > 4 ) THEN
     IF ( narg /= 10 ) THEN
        PRINT *, ' ERROR : You must give 6 optional values (imin imax jmin jmax kmin kmax)'
        STOP
     ELSE
        ! input optional imin imax jmin jmax
        CALL getarg ( 5,cdum) ; READ(cdum,*) imin
        CALL getarg ( 6,cdum) ; READ(cdum,*) imax
        CALL getarg ( 7,cdum) ; READ(cdum,*) jmin
        CALL getarg ( 8,cdum) ; READ(cdum,*) jmax
        CALL getarg ( 9,cdum) ; READ(cdum,*) kmin
        CALL getarg ( 10,cdum) ; READ(cdum,*) kmax
     ENDIF
  ENDIF

  cdep='none'
  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth',cdtrue=cdep,kstatus=istatus)
  

  !LB: fix bug SSH with ORCA025.L75
  !l_treat_ssh = .FALSE.
  !IF ( (trim(cvar) == 'sossheig').AND.(npiglo == 1442).AND.(npjglo == 1021) ) THEN
  !   PRINT *, '' ; PRINT *, 'BUG !!! This is sossheig on ORCA025 !!!'; PRINT *, ''
  !   l_treat_ssh = .TRUE.
  !END IF



  
  !LB:
  l_pvp = .FALSE.
  IF ( (kmin == 0).AND.(kmax == 0) ) l_pvp = .TRUE.
  !PRINT *, 'npk, kmin, kmax, l_pvp =', npk, kmin, kmax, l_pvp; PRINT *, ''
  !LB.
  

  IF (istatus /= 0 ) THEN
     npk   = getdim (cfilev,'z',cdtrue=cdep,kstatus=istatus)
     IF (istatus /= 0 ) THEN
       npk   = getdim (cfilev,'sigma',cdtrue=cdep,kstatus=istatus)
        IF ( istatus /= 0 ) THEN
          npk = getdim (cfilev,'nav_lev',cdtrue=cdep,kstatus=istatus)
            IF ( istatus /= 0 ) THEN
              PRINT *,' assume file with no depth'
              npk=0
            ENDIF
        ENDIF
     ENDIF
  ENDIF
  
  ctim = 'none'
  nt    = getdim (cfilev,'time',cdtrue=ctim,kstatus=istatus) !LB

  nvpk  = getvdim(cfilev,cvar)

  IF (npk == 0  ) THEN ; npk = 1              ; ENDIF  ! no depth dimension ==> 1 level
  IF (imin /= 0 ) THEN ; npiglo=imax -imin + 1;  ELSE ; imin=1 ; ENDIF
  IF (jmin /= 0 ) THEN ; npjglo=jmax -jmin + 1;  ELSE ; jmin=1 ; ENDIF
  IF (kmin /= 0 ) THEN ; npk   =kmax -kmin + 1;  ELSE ; kmin=1 ; ENDIF

  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk
  



  !LB:
  OPEN(16, FILE = cf_asout, FORM='FORMATTED', RECL=256)


  WRITE(16, *) 'npiglo=', npiglo
  WRITE(16, *) 'npjglo=', npjglo
  WRITE (16,*) 'npk   =', npk
  WRITE (16,*) 'nt    =', nt
  WRITE (16,*) 'nvpk  =', nvpk
  WRITE (16,*) 'depth dim name is ', TRIM(cdep)



  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo),e2(npiglo,npjglo), e3(npiglo,npjglo) )
  ALLOCATE ( gdep (npk) )
  SELECT CASE (TRIM(cvartype))
     CASE ( 'T' )
        ce1='e1t'
        ce2='e2t'
        ce3='e3t_0'
        cvmask='tmask'
        cdep='gdept'
     CASE ( 'U' )
        ce1='e1u'
        ce2='e2u'
        ce3='e3t_0'
        cvmask='umask'
        cdep='gdept'
     CASE ( 'V' )
        ce1='e1v'
        ce2='e2v'
        ce3='e3t_0'
        cvmask='vmask'
        cdep='gdept'
     CASE ( 'F' )
        ce1='e1f'
        ce2='e2f'
        ce3='e3t_0'
        cvmask='fmask'
        cdep='gdept'
     CASE ( 'W' )
        ce1='e1t'
        ce2='e2t'
        ce3='e3w_0'
        cvmask='tmask'
        cdep='gdepw'
     CASE DEFAULT
        PRINT *, 'this type of variable is not known :', TRIM(cvartype)
        STOP
  END SELECT


  e1(:,:) = getvar(coordhgr, ce1, 1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  e2(:,:) = getvar(coordhgr, ce2, 1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  gdep(:) = getvare3(coordzgr,cdep,npk)


  IF(lwrtcdf) THEN
      ALLOCATE ( typvar(nvars), ipk(nvars), id_varout(nvars) )
      ALLOCATE (dumlon(kx,ky) , dumlat(kx,ky), dummymean(kx,ky) )
      ALLOCATE ( meanout(npk) , year_mean(npk) )

      dumlon(:,:)=0.
      dumlat(:,:)=0.

      ipk(1)=npk ! mean for each level
      ipk(2)=1   ! 3D mean


      ierr=getvaratt (cfilev,cvar,cdunits, &
                pmissing_value, cdlong_name, cdshort_name)

      ! define new variables for output 
      typvar(1)%name='mean_'//TRIM(cvar)
      typvar%units=TRIM(cdunits)
      typvar%missing_value=99999.
      typvar%valid_min= -1000.
      typvar%valid_max= 1000.
      typvar%scale_factor= 1.
      typvar%add_offset= 0.
      typvar%savelog10= 0.
      typvar(1)%long_name='mean_'//TRIM(cdlong_name)
      typvar(1)%short_name='mean_'//TRIM(cdshort_name)
      typvar%online_operation='N/A'
      typvar%axis='ZT'

      typvar(2)%name='mean_3D'//TRIM(cvar)
      typvar(2)%long_name='mean_3D'//TRIM(cdlong_name)
      typvar(2)%short_name='mean_3D'//TRIM(cdshort_name)
      typvar%online_operation='N/A'
      typvar%axis='T'
   ENDIF




   !LB:
   OPEN(numout,FILE=cfilout)
   WRITE(numout,*) '# Time, '//trim(cvar)//',   Depth,    level'
   !LB.


   IF ( nt < 1 ) THEN
      PRINT *, 'ERROR: the file contains no time records! nt =', nt
      STOP
   END IF


   DO jt=1,nt

      zvol=0.d0
      zsum=0.d0

      DO jk = 1,nvpk
         ik = jk+kmin-1
         ! Get velocities v at ik
!         zv(:,:)= getvar(cfilev, cvar,  ik ,npiglo,npjglo,kimin=imin,kjmin=jmin)
          zv(:,:)= getvar(cfilev, cvar,  ik ,npiglo,npjglo,kimin=imin,kjmin=jmin,ktime=jt)
          !!
          !!LB:
          !IF ( l_treat_ssh ) zv(:,npjglo) = zv(:,npjglo-1)
          !!LB.
          !!

         zmask(:,:)=getvar(cmask,cvmask,ik,npiglo,npjglo,kimin=imin,kjmin=jmin)
         !    zmask(:,npjglo)=0.  
         ! get e3 at level ik ( ps...)
         e3(:,:) = getvar(coordzgr, ce3, ik,npiglo,npjglo,kimin=imin,kjmin=jmin, ldiom=.TRUE.)
         !
         !LB:
         WHERE ( zmask == 0. ) zv = 0. ! avoid problem with earth processors that have masked values
         !LB.
         ! 
         !!
         zsurf=SUM(e1 * e2 * zmask)
         zvol2d=SUM(e1 * e2 * e3 * zmask)
         zvol=zvol+zvol2d
         zsum2d=SUM(zv*e1*e2*e3*zmask)
         zsum=zsum+zsum2d
         IF (zvol2d /= 0 )THEN
            WRITE(16,*)' Mean value at level ',ik,'(',gdep(ik),' m) ',zsum2d/zvol2d, 'surface = ',zsurf/1.e6,' km^2'
            !WRITE(numout,9004) gdep(ik), ik, zsum2d/zvol2d
            WRITE(numout,*) jt, zsum2d/zvol2d, gdep(ik), ik
            !IF (lwrtcdf) meanout(jk)=zsum2d/zvol2d
            year_mean(jk) = zsum2d/zvol2d !LB
            !!
         ELSE
            WRITE(16,*) ' No points in the water at level ',ik,'(',gdep(ik),' m) '
            year_mean(jk) = 0.
            IF (lwrtcdf) meanout(jk)=99999.
         ENDIF
      END DO
      WRITE(16,*) ' Mean value over the ocean: ', zsum/zvol, jt
      threedmeanout=zsum/zvol
   END DO
   CLOSE(numout)
9004 FORMAT(f9.2,' ',i2,' ',f9.2)

   CLOSE(16)
   !LB.

   
   !LB:
   IF ( l_pvp ) THEN
      cf_out = trim(cvar)//'_mean_vert_profile.dat'
      OPEN(14, FILE=cf_out)
      DO jk = 1,nvpk
         ik = jk+kmin-1
         WRITE(14,*) gdep(ik), year_mean(ik)
      END DO
      CLOSE(14)
   END IF
   !LB.

   IF(lwrtcdf) THEN
      
      !LB:
      cfileoutnc=trim(cvar)//'_cdfmean.nc'
      !LB.
      
      ! create output fileset
      ncout =create(cfileoutnc,'none',kx,ky,npk,cdep=cdep)
      ierr= createvar(ncout,typvar,nvars,ipk,id_varout )
      ierr= putheadervar(ncout, cfilev ,kx, &
           ky,npk,pnavlon=dumlon,pnavlat=dumlat,pdep=gdep,cdep=cdep)
      tim=getvar1d(cfilev,'time_counter',1)
      ierr=putvar1d(ncout,tim,1,'T')

      ! netcdf output 
      DO jk=1, nvpk
         dummymean(1,1)=meanout(jk)
         ierr = putvar(ncout, id_varout(1), dummymean, jk, kx, ky )
      END DO

      ierr=putvar0d(ncout,id_varout(2), threedmeanout )

      ierr = closeout(ncout)

   ENDIF

 END PROGRAM cdfmean
