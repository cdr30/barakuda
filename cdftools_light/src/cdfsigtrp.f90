PROGRAM cdfsigtrp

  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdfsigtrp  ***
  !!
  !!  **  Purpose: Compute density class Mass Transports  across a section
  !!               PARTIAL STEPS version
  !!
  !!  **  Method:
  !!        -The beginner and end point of the section are given in term of f-points index.
  !!        -The progracdfm works for zonal or meridional sections.
  !!        -The section definitions are given in an ASCII FILE dens_section.dat
  !!            foreach sections, 2 lines : (i) : section name (String, no blank)
  !!                                       (ii) : imin imax jmin jmax for the section
  !!        -Only vertical slices corrsponding to the sections are read in the files.
  !!            read metrics, depth, etc
  !!            read normal velocity (either uo or vo )
  !!            read 2 rows of T and S ( i i+1  or j j+1 )
  !!                compute the mean value at velocity point
  !!                compute sigma0 (can be easily modified for sigmai )
  !!            compute the depths of isopyncal surfaces
  !!            compute the transport from surface to the isopycn
  !!            compute the transport in each class of density
  !!            compute the total transport (for information)
  !!
  !! history :
  !!   Original :  J.M. Molines March 2006
  !!            :  R. Dussin (Jul. 2009) add cdf output
  !!---------------------------------------------------------------------
  !!  $Rev: 322 $
  !!  $Date: 2010-05-20 07:14:33 +0200 (Thu, 20 May 2010) $
  !!  $Id: cdfsigtrp.f90 322 2010-05-20 05:14:33Z molines $
  !!--------------------------------------------------------------
  !! * Modules used
  USE netcdf
  USE cdfio
  USE io_ezcdf
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: nbins                                  !: number of density classes
  INTEGER   :: ji, jk, jt, jclass, jsec,jiso , jbin,jarg  !: dummy loop index
  INTEGER   :: ipos                                   !: working variable
  INTEGER   :: narg, iargc                            !: command line
  INTEGER   :: npiglo,npjglo, npk, nk, nt                                !: vertical size, number of wet layers in the section
  INTEGER   :: numout=11                              !: ascii output

  INTEGER                            :: nsection                    !: number of sections (overall)
  INTEGER ,DIMENSION(:), ALLOCATABLE :: imina, imaxa, jmina, jmaxa  !: sections limits
  INTEGER                            :: imin, imax, jmin, jmax      !: working section limits
  INTEGER                            :: npts                        !: working section number of h-points
  ! added to write in netcdf
  INTEGER :: kx=1, ky=1                ! dims of netcdf output file
  INTEGER :: nboutput=2                ! number of values to write in cdf output
  INTEGER :: ncout, ierr, istatus      ! for netcdf output
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout

  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: gdept, gdepw !: depth of T and W points
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zs, zt       !: salinity and temperature from file
  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: XOUT, tmpm, tmpz   !: temporary arrays

  ! double precision for cumulative variables and densities
  REAL(KIND=8), DIMENSION (:),     ALLOCATABLE ::  eu                 !: either e1v or e2u
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  zu, e3 , zmask     !: velocities e3 and umask
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  zsig ,gdepu        !: density, depth of vel points
  REAL(KIND=8)                                 :: sigma_min, sigma_max,dsigma   !: Min and Max for sigma bining
  REAL(KIND=8)                                 :: sigma,zalfa                   !: current working sigma
  REAL(KIND=8), DIMENSION (:),   ALLOCATABLE   :: sigma_lev, sigma_lev_center   !: built array with sigma levels
  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE   :: hiso                          !: depth of isopycns

  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE   :: zwtrp, zwtrpbin       !: transport arrays
  ! added to write in netcdf
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  dumlon, dumlat
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE   ::  pdep

  REAL(KIND=4) ,DIMENSION(:), ALLOCATABLE  ::  tim !LB

  REAL(KIND=4), DIMENSION (1)               ::  dummy1, dummy2
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar  ! structure of output
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvarin  !: structure for recovering input informations such as iwght
  CHARACTER(LEN=256), DIMENSION(:),ALLOCATABLE :: cvarname !: names of input variables
  INTEGER                                   :: nvarin    !: number of variables in input file
  INTEGER                                   :: iweight   !: weight of input file for further averaging

  CHARACTER(LEN=256) :: cd_out = '.', cf_out
  REAL :: ryear

  CHARACTER(LEN=256), DIMENSION (:), ALLOCATABLE :: csection                     !: section name
  CHARACTER(LEN=256) :: cfilet, cfileu, cfilev, cfilesec='dens_section.dat'      !: files name
  CHARACTER(LEN=256) :: cf_mm='mesh_mask.nc'                                     !: coordinates files
  CHARACTER(LEN=256) :: cfilout='trpsig.txt'                                     !: output file
  CHARACTER(LEN=256) :: cdum, cf_dep_iso, cf_trp_srf, cf_trp_bin !LB             !: dummy string
  ! added to write in netcdf
  CHARACTER(LEN=256) :: cfileoutnc , ctim
  CHARACTER(LEN=256) :: cdunits, cdlong_name, cdshort_name, cdep

  LOGICAL    :: l_merid                     !: flag is true for meridional working section
  LOGICAL    :: lfncout = .false.

  CHARACTER(LEN=80) :: cfor9000, cfor9001, cfor9002, cfor9003, cfor9004


  REAL(4), DIMENSION(:,:),   ALLOCATABLE :: E2U_2D, E1V_2D
  REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: T_3D, S_3D, U_3D, V_3D, E3U_3D, E3V_3D, E3W_3D

  INTEGER, DIMENSION(:), ALLOCATABLE :: id_var_out
  
  CHARACTER(LEN=64) :: cv_t, cv_s, cv_u, cv_v
 
  INTEGER :: jt_pos, idf_out, idd_t, idd_sbins, idv_time, idv_sbins

  INTEGER :: idf_0=0, idv_0=0, &
       &     idf_t=0, idv_t=0, idf_s=0, idv_s=0, idf_u=0, idv_u=0, idf_v=0, idv_v=0










  !!  * Initialisations

  !  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 12  ) THEN
     PRINT *, ''
     PRINT *,' Usage: cdfsigtrp.x <Tfile> <Ufile> <Vfile> s_min s_max nbins <year> <DIROUT> <name T> <name S> <name U> <name V>'
     PRINT *,'         *   s_min, s_max : limit for density bining '
     PRINT *,'         *   nbins        : number of bins to use '
     PRINT *, ''
     PRINT *,' File mesh_mask.nc must be in the current directory'
     PRINT *,' File dens_section.dat must also be in the current directory '
     STOP
  ENDIF

  !! Read arguments
  CALL getarg (1, cfilet)
  CALL getarg (2, cfileu)
  CALL getarg (3, cfilev)
  CALL getarg (4,cdum)           ; READ(cdum,*) sigma_min
  CALL getarg (5,cdum)           ; READ(cdum,*) sigma_max
  CALL getarg (6,cdum)           ; READ(cdum,*) nbins
  CALL getarg (7, cdum)          ; READ(cdum,*) ryear
  CALL getarg (8, cd_out)
  CALL getarg (9,  cv_t)
  CALL getarg (10, cv_s)
  CALL getarg (11, cv_u)
  CALL getarg (12, cv_v)


  ! start for netcdf:
  nvarin=getnvar(cfileu)  ! smaller than cfilet
  ALLOCATE(typvarin(nvarin), cvarname(nvarin)  )
  cvarname(:)=getvarname(cfileu,nvarin,typvarin)

  DO  jarg=1,nvarin
     IF ( TRIM(cvarname(jarg))  == trim(cv_u) ) THEN
        iweight=typvarin(jarg)%iwght
        EXIT  ! loop
     ENDIF
  END DO


  ALLOCATE ( typvar(nboutput), ipk(nboutput), id_varout(nboutput) )
  ALLOCATE (dumlon(kx,ky) , dumlat(kx,ky) )

  dumlon(:,:)=0.
  dumlat(:,:)=0.

  ipk(1)=nbins ! sigma for each level
  ipk(2)=nbins ! transport for each level

  ! define new variables for output
  typvar(1)%name='sigma_class'
  typvar%units='[]'
  typvar%missing_value=99999.
  typvar%valid_min= 0.
  typvar%valid_max= 100.
  typvar%scale_factor= 1.
  typvar%add_offset= 0.
  typvar%savelog10= 0.
  typvar%iwght=iweight
  typvar(1)%long_name='class of potential density'
  typvar(1)%short_name='sigma_class'
  typvar%online_operation='N/A'
  typvar%axis='ZT'

  typvar(2)%name='sigtrp'
  typvar(2)%units='Sv'
  typvar(2)%valid_min= -1000.
  typvar(2)%valid_max= 1000.
  typvar(2)%long_name='transport in sigma class'
  typvar(2)%short_name='sigtrp'

  !end of netcdf...



  INQUIRE( FILE=cfilesec, EXIST=lfncout )
  IF ( .NOT. lfncout ) THEN
     PRINT *, 'ERROR: file '//trim(cfilesec)//' not found!' ; STOP
  END IF




  ! Initialise sections from file
  ! first call to get nsection and allocate arrays
  nsection = 0 ; CALL section_init(cfilesec, csection,imina,imaxa,jmina,jmaxa, nsection)
  ALLOCATE ( csection(nsection), imina(nsection), imaxa(nsection), jmina(nsection),jmaxa(nsection) )
  CALL section_init(cfilesec, csection,imina,imaxa,jmina,jmaxa, nsection)


  ! Allocate and build sigma levels and section array
  ALLOCATE ( sigma_lev(nbins+1), sigma_lev_center(nbins) )

  sigma_lev(1)=sigma_min
  dsigma=( sigma_max - sigma_min) / nbins
  DO jclass =2, nbins+1
     sigma_lev(jclass)= sigma_lev(1) + (jclass-1) * dsigma
  END DO

  DO jclass =1, nbins
     sigma_lev_center(jclass) = 0.5*(sigma_lev(jclass) + sigma_lev(jclass+1))
  END DO

  ctim = 'none'
  nt    = getdim (cfilet,'time',cdtrue=ctim,kstatus=istatus) !LB


  ALLOCATE ( tim(nt), XOUT(nsection, nbins, nt), id_var_out(nsection) ) !LB



  ! Look for vertical size of the domain
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

  ALLOCATE ( gdept(npk), gdepw(npk) )


  ALLOCATE ( T_3D(npiglo,npjglo,npk), S_3D(npiglo,npjglo,npk), &
       &     U_3D(npiglo,npjglo,npk), V_3D(npiglo,npjglo,npk), &
       &     E2U_2D(npiglo,npjglo), E1V_2D(npiglo,npjglo),     &
       &     E3U_3D(npiglo,npjglo,npk), E3V_3D(npiglo,npjglo,npk), E3W_3D(npiglo,npjglo,npk) )



  ! read gdept, gdepw : it is OK even in partial cells, as we never use the bottom gdep
  gdept(:) = getvare3(cf_mm,'gdept', npk)
  gdepw(:) = getvare3(cf_mm,'gdepw', npk)





  !lolo
  ! get e3u, e3v  at all levels
  CALL GETVAR_2D(idf_0, idv_0, cf_mm, 'e2u',   0, 0, 0, E2U_2D) ; idf_0 = 0. ; idv_0 = 0.
  CALL GETVAR_2D(idf_0, idv_0, cf_mm, 'e1v',   0, 0, 0, E1V_2D) ; idf_0 = 0. ; idv_0 = 0.
  CALL GETVAR_3D(idf_0, idv_0, cf_mm, 'e3u_0', 0, 0, E3U_3D) ; idf_0 = 0. ; idv_0 = 0.
  CALL GETVAR_3D(idf_0, idv_0, cf_mm, 'e3v_0', 0, 0, E3V_3D) ; idf_0 = 0. ; idv_0 = 0.
  CALL GETVAR_3D(idf_0, idv_0, cf_mm, 'e3w_0', 0, 0, E3W_3D)




  !! *  Main loop on sections

  DO jt = 1, nt !LB

     PRINT *, 'Time record = ', jt

     CALL GETVAR_3D(idf_t, idv_t,  cfilet, trim(cv_t), nt, jt, T_3D)
     CALL GETVAR_3D(idf_s, idv_s,  cfilet, trim(cv_s),     nt, jt, S_3D)
     CALL GETVAR_3D(idf_u, idv_u,  cfileu, trim(cv_u),     nt, jt, U_3D)
     CALL GETVAR_3D(idf_v, idv_v,  cfilev, trim(cv_v),     nt, jt, V_3D)




     DO jsec=1,nsection

        PRINT *, 'Treating section '//TRIM(csection(jsec))

        l_merid=.FALSE.
        imin=imina(jsec) ; imax=imaxa(jsec) ; jmin=jmina(jsec) ; jmax=jmaxa(jsec)
        IF (imin == imax ) THEN        ! meridional section
           l_merid=.TRUE.
           npts=jmax-jmin

        ELSE IF ( jmin == jmax ) THEN  ! zonal section
           npts=imax-imin

        ELSE
           PRINT *,' Section ',TRIM(csection(jsec)),' is neither zonal nor meridional :('
           PRINT *,' We skip this section .'
           CYCLE
        ENDIF

        !PRINT *, 'Allocating for jt, jsec =', jt, jsec
        ALLOCATE ( zu(npts, npk), zt(npts,npk), zs(npts,npk) ,zsig(npts,0:npk) )
        ALLOCATE ( eu(npts), e3(npts,npk), gdepu(npts, 0:npk), zmask(npts,npk) )
        ALLOCATE ( tmpm(1,npts,2), tmpz(npts,1,2) )
        ALLOCATE ( zwtrp(npts, nbins+1) , hiso(npts,nbins+1), zwtrpbin(npts,nbins) )

        !END IF

        !PRINT *, ' filling with 0...'
        zt = 0. ; zs = 0. ; zu = 0. ; gdepu= 0. ; zmask = 0.  ; zsig=0.d0

        IF (l_merid ) THEN   ! meridional section at i=imin=imax

           !tmpm(:,:,1)=getvar(cf_mm, 'e2u', 1,1,npts, kimin=imin, kjmin=jmin+1)
           !PRINT *, 'LOLO 1 =>', tmpm(:,:,1)
           tmpm(1,:,1) = E2U_2D(imin,jmin+1:jmin+1+npts)
           !PRINT *, 'LOLO 2 =>', tmpm(:,:,1)
           !STOP

           eu(:)=tmpm(1,:,1)  ! metrics varies only horizontally
           DO jk=1,npk
              ! initiliaze gdepu to gdept()
              gdepu(:,jk) = gdept(jk)

              ! vertical metrics (PS case)
              !tmpm(:,:,1)=getvar(cf_mm,'e3u_ps',jk,1,npts, kimin=imin, kjmin=jmin+1, ldiom=.TRUE.)
              !PRINT *, 'LOLO1 =>', tmpm(:,:,1)
              tmpm(1,:,1) = E3U_3D(imin,jmin+1:jmin+1+npts,jk)
              !PRINT *, 'LOLO2 =>', tmpm(:,:,1)

              e3(:,jk)=tmpm(1,:,1)

              !tmpm(:,:,1)=getvar(cf_mm,'e3w_ps',jk,1,npts, kimin=imin, kjmin=jmin+1, ldiom=.TRUE.)
              !tmpm(:,:,2)=getvar(cf_mm,'e3w_ps',jk,1,npts, kimin=imin+1, kjmin=jmin+1, ldiom=.TRUE.)
              tmpm(1,:,1) = E3W_3D(imin,  jmin+1:jmin+1+npts,jk)
              tmpm(1,:,2) = E3W_3D(imin+1,jmin+1:jmin+1+npts,jk)


              IF (jk >= 2 ) THEN
                 DO ji=1,npts
                    gdepu(ji,jk)= gdepu(ji,jk-1) + MIN(tmpm(1,ji,1), tmpm(1,ji,2))
                 END DO
              ENDIF

              ! Normal velocity
              tmpm(1,:,1) = U_3D(imin,jmin+1:jmin+1+npts,jk)

              zu(:,jk)=tmpm(1,:,1)

              ! salinity and deduce umask for the section
              tmpm(1,:,1) = S_3D(imin,  jmin+1:jmin+1+npts,jk)
              tmpm(1,:,2) = S_3D(imin+1,jmin+1:jmin+1+npts,jk)

              zmask(:,jk)=tmpm(1,:,1)*tmpm(1,:,2)
              WHERE ( zmask(:,jk) /= 0 ) zmask(:,jk)=1
              ! do not take special care for land value, as the corresponding velocity point is masked
              zs(:,jk) = 0.5 * ( tmpm(1,:,1) + tmpm(1,:,2) )

              ! limitation to 'wet' points
              IF ( SUM(zs(:,jk))  == 0 ) THEN
                 nk=jk ! first vertical point of the section full on land
                 EXIT  ! as soon as all the points are on land
              ENDIF

              ! temperature
              tmpm(1,:,1) = T_3D(imin,  jmin+1:jmin+1+npts,jk)
              tmpm(1,:,2) = T_3D(imin+1,jmin+1:jmin+1+npts,jk)

              zt(:,jk) = 0.5 * ( tmpm(1,:,1) + tmpm(1,:,2) )

           END DO



        ELSE                   ! zonal section at j=jmin=jmax
           !tmpz(npts,1,2)
           !tmpz(:,:,1)=getvar(cf_mm, 'e1v', 1,npts,1,kimin=imin, kjmin=jmin)
           !PRINT *, 'LOLO1 tmpz =>', tmpz(:,:,1)
           tmpz(:,1,1) = E1V_2D(imin:imin+npts,jmin)
           !PRINT *, 'LOLO2 tmpz =>', tmpz(:,:,1)

           eu=tmpz(:,1,1)
           DO jk=1,npk
              ! initiliaze gdepu to gdept()
              gdepu(:,jk) = gdept(jk)

              ! vertical metrics (PS case)
              !tmpz(:,:,1)=getvar(cf_mm,'e3v_ps',jk, npts, 1, kimin=imin+1, kjmin=jmin, ldiom=.TRUE.)
              tmpz(:,1,1) = E3V_3D(imin+1:imin+1+npts,jmin,jk)

              e3(:,jk)=tmpz(:,1,1)
              !tmpz(:,:,1)=getvar(cf_mm,'e3w_ps',jk,npts,1, kimin=imin+1, kjmin=jmin, ldiom=.TRUE.)
              !tmpz(:,:,2)=getvar(cf_mm,'e3w_ps',jk,npts,1, kimin=imin+1, kjmin=jmin+1, ldiom=.TRUE.)
              tmpz(:,1,1) = E3W_3D(imin+1:imin+1+npts,jmin,  jk)
              tmpz(:,1,2) = E3W_3D(imin+1:imin+1+npts,jmin+1,jk)

              IF (jk >= 2 ) THEN
                 DO ji=1,npts
                    gdepu(ji,jk)= gdepu(ji,jk-1) + MIN(tmpz(ji,1,1), tmpz(ji,1,2))
                 END DO
              ENDIF

              ! Normal velocity
              tmpz(:,1,1) = V_3D(imin+1:imin+1+npts,jmin,  jk)
              zu(:,jk)=tmpz(:,1,1)

              ! salinity and mask
              tmpz(:,1,1) = S_3D(imin+1:imin+1+npts,jmin,  jk)
              tmpz(:,1,2) = S_3D(imin+1:imin+1+npts,jmin+1,jk)

              zmask(:,jk)=tmpz(:,1,1)*tmpz(:,1,2)
              WHERE ( zmask(:,jk) /= 0 ) zmask(:,jk)=1
              ! do not take special care for land value, as the corresponding velocity point is masked
              zs(:,jk) = 0.5 * ( tmpz(:,1,1) + tmpz(:,1,2) )

              ! limitation to 'wet' points
              IF ( SUM(zs(:,jk))  == 0 ) THEN
                 nk=jk ! first vertical point of the section full on land
                 EXIT  ! as soon as all the points are on land
              ENDIF

              ! temperature
              tmpz(:,1,1) = T_3D(imin+1:imin+1+npts,jmin,  jk)
              tmpz(:,1,2) = T_3D(imin+1:imin+1+npts,jmin+1,jk)

              zt(:,jk) = 0.5 * ( tmpz(:,1,1) + tmpz(:,1,2) )
           END DO

        ENDIF  !(l_merid)

        ! compute density only for wet points
        zsig(:,1:nk)=sigma0( zt, zs, npts, nk)*zmask(:,:)
        zsig(:,0)=zsig(:,1)-1.e-4   ! dummy layer for easy interpolation



        WRITE(cfor9000,'(a,i3,a)') '(i7,',npts,'f8.3)'
        WRITE(cfor9001,'(a,i3,a)') '(i7,',npts,'f8.0)'
        WRITE(cfor9002,'(a,i3,a)') '(f7.3,',npts,'f8.0)'
        WRITE(cfor9003,'(a,i3,a)') '(f7.3,',npts,'f8.3)'
        WRITE(cfor9004,'(a,i3,a)') '(f7.3,',npts+1,'f8.3)'





        DO  jiso =1, nbins+1
           sigma=sigma_lev(jiso)
           DO ji=1,npts
              hiso(ji,jiso) = gdept(npk)
              DO jk=1,nk
                 IF ( zsig(ji,jk) < sigma ) THEN
                 ELSE
                    ! interpolate between jk-1 and jk
                    zalfa=(sigma - zsig(ji,jk-1)) / ( zsig(ji,jk) -zsig(ji,jk-1) )
                    IF (ABS(zalfa) > 1.1 .OR. zalfa < 0 ) THEN   ! case zsig(0) = zsig(1)-1.e-4
                       hiso(ji,jiso)= 0.
                    ELSE
                       hiso(ji,jiso)= gdepu(ji,jk)*zalfa + (1.-zalfa)* gdepu(ji,jk-1)
                    ENDIF
                    EXIT
                 ENDIF
              END DO
           END DO
        END DO


        DO jiso = 1, nbins + 1
           sigma=sigma_lev(jiso)
           DO ji=1,npts
              zwtrp(ji,jiso) = 0.d0
              DO jk=1, nk-1
                 IF ( gdepw(jk+1) < hiso(ji,jiso) ) THEN
                    zwtrp(ji,jiso)= zwtrp(ji,jiso) + eu(ji)*e3(ji,jk)*zu(ji,jk)
                 ELSE  ! last box ( fraction)
                    zwtrp(ji,jiso)= zwtrp(ji,jiso) + eu(ji)*(hiso(ji,jiso)-gdepw(jk))*zu(ji,jk)
                    EXIT  ! jk loop
                 ENDIF
              END DO
           END DO
        END DO

        !DO jbin=1, nbins
        !   sigma=sigma_lev(jbin)
        !   DO ji=1, npts
        !      zwtrpbin(ji,jbin) = zwtrp(ji,jbin+1) -  zwtrp(ji,jbin)
        !   END DO
        !   trpbin(jsec,jbin)=SUM(zwtrpbin(:,jbin) )
        !END DO


        DO jbin=1, nbins
           sigma=sigma_lev(jbin)
           DO ji=1, npts
              zwtrpbin(ji,jbin) = zwtrp(ji,jbin+1) -  zwtrp(ji,jbin)
           END DO
           XOUT(jsec,jbin,jt)=SUM(zwtrpbin(:,jbin))/1.e6
        END DO


        ! free memory for the next section
        DEALLOCATE ( zu,zt, zs ,zsig ,gdepu, hiso, zwtrp, zwtrpbin )
        DEALLOCATE ( eu, e3 ,tmpm, tmpz,zmask )

        CLOSE(15)


     END DO  !DO jsec=1,nsection

  END DO   !DO jt = 1, nt



  WRITE(cf_out, '(a,"/transport_by_sigma_class.nc")') trim(cd_out)

  id_var_out(:) = 0


  INQUIRE( FILE=cf_out, EXIST=lfncout )


  IF ( .NOT. lfncout ) THEN

     !! Creating file
     PRINT *, ' Creating file '//trim(cf_out)//' !!!'
     ierr = NF90_CREATE(cf_out, NF90_CLOBBER, idf_out)

     ierr = NF90_DEF_DIM(idf_out, 'time',      NF90_UNLIMITED, idd_t)
     ierr = NF90_DEF_DIM(idf_out, 'sigma_bins', nbins, idd_sbins)

     ierr = NF90_DEF_VAR(idf_out, 'time',       NF90_DOUBLE, idd_t,     idv_time)

     ierr = NF90_DEF_VAR(idf_out, 'sigma_bins', NF90_DOUBLE, idd_sbins, idv_sbins)
     WRITE(cdum,'(f8.4)') dsigma
     ierr = NF90_PUT_ATT(idf_out, idv_sbins, 'long_name', 'Center of sigma bin (dsigma = '//TRIM(cdum)//')')


     DO jsec=1, nsection
        ierr = NF90_DEF_VAR(idf_out, 'sigtrsp_'//TRIM(csection(jsec)), &
             & NF90_FLOAT, (/idd_sbins,idd_t/), id_var_out(jsec))
        ierr = NF90_PUT_ATT(idf_out, id_var_out(jsec), 'long_name', 'Transport by sigma class in section '//TRIM(csection(jsec)))
        ierr = NF90_PUT_ATT(idf_out, id_var_out(jsec), 'units', 'Sv')
     END DO

     ierr = NF90_PUT_ATT(idf_out, NF90_GLOBAL, 'About', 'Created by BaraKuda (cdfsigtrp.f90), contact: brodeau@gmail.com')
     ierr = NF90_ENDDEF(idf_out)
     jt_pos = 0

     ierr = NF90_PUT_VAR( idf_out, idv_sbins, sigma_lev_center(:))

  ELSE

     !! Opening already existing file
     ierr = NF90_OPEN  (cf_out, NF90_WRITE, idf_out)

     !! Need IDs of variables to append... NF90_INQ_VARID
     DO jsec=1, nsection
        ierr = NF90_INQ_VARID(idf_out, 'sigtrsp_'//TRIM(csection(jsec)), id_var_out(jsec))
     END DO


     ! Get ID of unlimited dimension
     ierr = NF90_INQUIRE(idf_out, unlimitedDimId = idv_time)

     ! Need to know jt_pos, record number of the last time record writen in the file
     ierr = NF90_INQUIRE_DIMENSION(idf_out, idv_time, name=cdum, len=jt_pos)

  END IF


  WRITE(*,'("Going to write record ",i4.4," to ",i4.4," into ",a)') jt_pos+1, jt_pos+1+nt, trim(cf_out)

  DO jt = 1, nt

     ! Writing record jt for time vector and 1d fields:
     ierr = NF90_PUT_VAR( idf_out, idv_time, (/ryear+1./12.*(REAL(jt)-1.+0.5)/), start=(/jt_pos+jt/), count=(/1/) )

     DO jsec=1, nsection
        ierr = NF90_PUT_VAR(idf_out, id_var_out(jsec), XOUT(jsec,:,jt), start=(/1,jt_pos+jt/), count=(/nbins,1/))
     END DO

  END DO


  ierr = NF90_CLOSE(idf_out)


CONTAINS



  SUBROUTINE section_init(cdfile,cdsection,kimin,kimax,kjmin,kjmax,knumber)
    IMPLICIT NONE
    ! Arguments
    !   INTEGER, DIMENSION(:),ALLOCATABLE :: kimin,kimax, kjmin,kjmax
    INTEGER, INTENT(INOUT) :: knumber
    INTEGER, DIMENSION(knumber) :: kimin,kimax, kjmin,kjmax
    !   CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cdsection
    CHARACTER(LEN=256), DIMENSION(knumber) :: cdsection
    CHARACTER(LEN=*), INTENT(IN) :: cdfile

    ! Local variables
    INTEGER :: ii, numit=10, jsec
    CHARACTER(LEN=256) :: cline
    LOGICAL :: lfirst

    INTEGER :: idf_0=0, idv_0=0, &
         idf_t=0, idv_t=0, idf_s=0, idv_s=0, idf_u=0, idv_u=0, idf_v=0, idv_v=0





    lfirst=.FALSE.
    IF ( knumber == 0 ) lfirst=.TRUE.

    OPEN(numit, FILE=cdfile, STATUS='old')
    REWIND(numit)
    ii=0

    DO
       READ(numit,'(a)') cline
       IF (INDEX(cline,'EOF') == 0 ) THEN
          READ(numit,*)    ! skip one line
          ii = ii + 1
       ELSE
          EXIT
       ENDIF
    END DO


    knumber=ii
    IF ( lfirst ) RETURN
    !   ALLOCATE( cdsection(knumber) )
    !   ALLOCATE( kimin(knumber), kimax(knumber), kjmin(knumber), kjmax(knumber) )
    REWIND(numit)
    DO jsec=1,knumber
       READ(numit,'(a)') cdsection(jsec)
       READ(numit,*) kimin(jsec), kimax(jsec), kjmin(jsec), kjmax(jsec)
    END DO

    CLOSE(numit)

  END SUBROUTINE section_init


END PROGRAM cdfsigtrp
