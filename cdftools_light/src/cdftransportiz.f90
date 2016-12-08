PROGRAM cdftransportiz

  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdftransportiz  ***
  !!
  !!  **  Purpose: Compute Transports across a section
  !!               PARTIAL STEPS version
  !!
  !!  **  Method: Try to avoid 3 d arrays.
  !!             The begining and end point of the section are given in term of f-points index.
  !!             This program computes the transport across this section for
  !!               (1) Mass transport ( Sv)
  !!               (2) Heat Transport (PW)
  !!               (3) Salt Transport (kT/sec)
  !!             The transport is > 0 left handside of the line
  !!             This program use a zig-zag line going through U and V-points.
  !!             It takes as input : VT files, gridU, gridV files.
  !!             The mesh_mask.nc, mesh_hzr.nc are required.
  !!             It is convenient to use an ASCII file as the standard input to give
  !!             the name and the imin imax jmin jmax for each section required
  !!             The last name of this ASCII file must be EOF
  !!
  !!
  !! history :
  !!   Original :  J.M. Molines (jan. 2005)
  !!               J.M. Molines Apr 2005 : use modules
  !!               J.M. Molines Apr 2007 : merge with Julien Jouanno version (std + file output)
  !!               R. Dussin (Jul. 2009) : add cdf output
  !!---------------------------------------------------------------------
  !!  $Rev: 257 $
  !!  $Date: 2009-07-27 18:25:04 +0200 (Mon, 27 Jul 2009) $
  !!  $Id: cdftransportiz.f90 257 2009-07-27 16:25:04Z forge $
  !!--------------------------------------------------------------
  !! * Modules used

  USE netcdf

  USE cdfio
  USE io_ezcdf


  !! * Local variables
  IMPLICIT NONE

  LOGICAL, PARAMETER :: l_save_broken_lines = .TRUE. !lolo

  INTEGER :: nclass   !: number of depth class
  INTEGER ,DIMENSION (:),ALLOCATABLE ::  imeter  !: limit beetween depth level, in m (nclass -1)
  CHARACTER(len=64), DIMENSION (:),ALLOCATABLE ::  cdepths
  INTEGER ,DIMENSION (:),ALLOCATABLE :: ilev0,ilev1 !: limit in levels  ! nclass
  INTEGER   :: jk, jc, jj, jt                  !: dummy loop index !LB
  INTEGER   :: narg, iargc, istatus                         !: command line
  INTEGER   :: npiglo,npjglo, npk, nt               !: size of the domain !LB
  INTEGER   :: imin, imax, jmin, jmax, ik
  INTEGER   :: numin  = 16

  ! broken line stuff
  INTEGER, PARAMETER :: jpseg=10000
  INTEGER :: i0,j0,i1,j1, i, j
  INTEGER :: n,nn,k, jseg
  INTEGER :: norm_u, norm_v, ist, jst

  REAL(4) ::  rxi0,ryj0, rxi1, ryj1

  REAL(8) :: ref_temp=0

  REAL(4) ::   ai,bi, aj,bj,d
  REAL(4) ::    rxx(jpseg),ryy(jpseg)
  REAL(4), DIMENSION(jpseg) :: gla, gphi

  REAL(8), DIMENSION(jpseg) :: voltrp, heatrp, saltrp
  REAL(8)                   :: voltrpsum, heatrpsum, saltrpsum
  COMPLEX yypt(jpseg), yypti

  REAL(4), DIMENSION (:,:),   ALLOCATABLE ::         e1v, e3v ,gphiv, zv, zvt, zvs !: mask, metrics
  REAL(4), DIMENSION (:,:),   ALLOCATABLE ::         e2u, e3u ,gphiu, zu, zut, zus !: mask, metrics
  REAL(4), DIMENSION (:,:),   ALLOCATABLE ::         glamu, glamv
  REAL(4), DIMENSION (:),     ALLOCATABLE ::         gdepw, gdept
  REAL(4)                                 ::   rd1, rd2
  REAL(4)                                 ::  udum, vdum

  REAL(8),   DIMENSION (:,:), ALLOCATABLE :: zwku,zwkv,    zwkut,zwkvt,   zwkus,zwkvs
  REAL(8),   DIMENSION (:,:,:), ALLOCATABLE :: ztrpu, ztrpv, ztrput,ztrpvt, ztrpus,ztrpvs

  CHARACTER(len=8) :: ctest, crt

  CHARACTER(LEN=256) :: conf_tag, cf_vt , cf_u, cf_v, csection, cf_broken_line,  &
       &                cfilvtrp='vtrp.txt', cfilhtrp='htrp.txt', cfilstrp='strp.txt'

  CHARACTER(LEN=256) :: cf_mm='mesh_mask.nc', cdum
  CHARACTER(LEN=256) ,DIMENSION(4)   :: cvarname   !: array of var name for output

  INTEGER :: nxtarg

  LOGICAL :: ltest=.FALSE., &
       &     leiv  = .TRUE., &    !: weather to use Eddy Induced velocity from GM90
       &     lcontinue = .TRUE.
  
  CHARACTER(LEN=256) :: ctim

  ! constants
  !lolo: # The density of seawater is about 1025 kg/m^3 and the specific heat is about 4000 J/(kg C)
  REAL(4)   ::  rau0=1025.,  rcp=3990.


  REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: U_3D, V_3D, UEIV_3D, VEIV_3D, UT_3D, VT_3D, US_3D, VS_3D, E3V_3D, E3U_3D

  INTEGER :: idf_u=0, idv_u=0, idf_v=0, idv_v=0, idf_ueiv=0, idv_ueiv=0, idf_veiv=0, idv_veiv=0, &
       &     idf_ut=0, idv_ut=0, idf_vt=0, idv_vt=0, &
       &     idf_us=0, idv_us=0, idf_vs=0, idv_vs=0, idf_0=0, idv_0=0


  !! LOLO:
  LOGICAL :: lfncout = .false.
  CHARACTER(LEN=256) :: cd_out = '.', cf_out
  CHARACTER(LEN=64)  :: cv_u, cv_v, cv_ueiv, cv_veiv, cv_dum
  REAL :: ryear
  INTEGER :: ierr, jt_pos, idf_out, idd_t, idv_time
  INTEGER, DIMENSION(:)    , ALLOCATABLE :: id_volu, id_heat, id_salt
  REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: X_trsp   ! lolo
  !! LOLO.



  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 7  ) THEN
     PRINT *,'Usage : cdftransportiz <CONFTAG> <nameU> <nameV> <nameUeiv> & 
          & <nameVeiv> <year> <DIROUT> "limit of level" (<ref temp.degree>)'
     PRINT *,'    => files are: <CONFTAG>_VT.nc <CONFTAG>_grid_U.nc <CONFTAG>_grid_V.nc'
     PRINT *, '  If eddy-induced velocity is not relevant, specify "0" "0" for <nameUeiv> <nameVeiv>'
     PRINT *,' Files mesh_mask.nc must be in te current directory'
     STOP
  ENDIF

  CALL getarg (1, conf_tag)
  CALL getarg (2, cv_u)
  CALL getarg (3, cv_v)
  CALL getarg (4, cv_ueiv)
  CALL getarg (5, cv_veiv)
  CALL getarg (6, cdum)    ; READ(cdum,*) ryear
  CALL getarg (7, cd_out)

  nxtarg=7
  nclass = narg -nxtarg + 1

  leiv = .TRUE.
  IF ( (trim(cv_ueiv) == '0').AND.(trim(cv_veiv) == '0') ) THEN
     leiv = .FALSE.
     PRINT *, ' Not taking eddy-induced velocity into account!'
  ELSE
     PRINT *, ' Taking eddy-induced velocity into account using '//trim(cv_ueiv)//' and '//trim(cv_veiv)
  END IF


  WRITE(cf_vt, '(a,"_VT.nc")')     trim(conf_tag)
  WRITE(cf_u,  '(a,"_grid_U.nc")') trim(conf_tag)
  WRITE(cf_v,  '(a,"_grid_V.nc")') trim(conf_tag)




  ALLOCATE (  cdepths(0:nclass), imeter(nclass -1), ilev0(nclass), ilev1(nclass) )

  cdepths(0) = '0'
  cdepths(nclass) = 'botto'
  DO jk=1, nclass -1
     CALL getarg(nxtarg+jk,cdepths(jk))
     READ(cdepths(jk),*) imeter(jk)
  END DO

  npiglo= getdim (cf_vt,'x')
  npjglo= getdim (cf_vt,'y')
  npk   = getdim (cf_vt,'depth')

  ctim = 'none'
  nt    = getdim (cf_vt,'time',cdtrue=ctim,kstatus=istatus) !LB

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'nt    =', nt !LB
  PRINT *, ''

  ALLOCATE( U_3D(npiglo,npjglo,npk) , V_3D(npiglo,npjglo,npk) , UT_3D(npiglo,npjglo,npk) , VT_3D(npiglo,npjglo,npk) , &
       &    US_3D(npiglo,npjglo,npk) , VS_3D(npiglo,npjglo,npk), E3V_3D(npiglo,npjglo,npk) , E3U_3D(npiglo,npjglo,npk) )

  IF ( leiv ) ALLOCATE( UEIV_3D(npiglo,npjglo,npk) , VEIV_3D(npiglo,npjglo,npk) )


  !! Lolo: allocating array to contain transports:
  ALLOCATE( X_trsp(3,nt,nclass) ) ; ! 3 transports, nt values, nclass levels
  ALLOCATE( id_volu(0:nclass), id_heat(0:nclass), id_salt(0:nclass) )

  ! Allocate arrays
  ALLOCATE( zu (npiglo,npjglo), zut(npiglo,npjglo), zus(npiglo,npjglo) )
  ALLOCATE( zv (npiglo,npjglo), zvt(npiglo,npjglo), zvs(npiglo,npjglo) )
  !
  ALLOCATE ( zwku (npiglo,npjglo), zwkut(npiglo,npjglo), zwkus(npiglo,npjglo) )
  ALLOCATE ( zwkv (npiglo,npjglo), zwkvt(npiglo,npjglo), zwkvs(npiglo,npjglo) )
  !
  ALLOCATE ( ztrpu (npiglo,npjglo,nclass), ztrpv (npiglo,npjglo,nclass))
  ALLOCATE ( ztrput(npiglo,npjglo,nclass), ztrpvt(npiglo,npjglo,nclass))
  ALLOCATE ( ztrpus(npiglo,npjglo,nclass), ztrpvs(npiglo,npjglo,nclass))
  !
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo))
  ALLOCATE ( e2u(npiglo,npjglo),e3u(npiglo,npjglo))
  !
  ALLOCATE ( gphiu(npiglo,npjglo),  gphiv(npiglo,npjglo) )
  ALLOCATE ( glamu(npiglo,npjglo),  glamv(npiglo,npjglo) )
  ALLOCATE ( gdepw(npk), gdept(npk) )
  !

  e1v(:,:) = getvar(cf_mm, 'e1v', 1,npiglo,npjglo)
  e2u(:,:) = getvar(cf_mm, 'e2u', 1,npiglo,npjglo)

  glamv(:,:) =  getvar(cf_mm, 'glamv', 1,npiglo,npjglo)
  glamu(:,:) =  getvar(cf_mm, 'glamu', 1,npiglo,npjglo)

  gphiv(:,:) = getvar(cf_mm, 'gphiv', 1,npiglo,npjglo)
  gphiu(:,:) = getvar(cf_mm, 'gphiu', 1,npiglo,npjglo)

  gdepw(:) = getvare3(cf_mm, 'gdepw_1d',npk)
  gdept(:) = getvare3(cf_mm, 'gdept_1d',npk)




  DO WHILE ( lcontinue )

     OPEN(numin, FILE='transportiz.dat', status='old')
     READ(numin,'(a)') csection
     IF (TRIM(csection) == 'EOF' ) THEN
        CLOSE(numin)
        lcontinue = .FALSE.
        EXIT
     END IF
     
     IF ( .NOT. lcontinue ) EXIT
     
     READ(numin,*) imin, imax, jmin, jmax
     PRINT*, 'Section = ', trim(csection); PRINT*, ' =>', imin, imax, jmin, jmax
     
     READ(numin,'(a)') ctest
     IF (TRIM(ctest) == 'ref_temp' ) THEN
        !READ(numin,'(f)') ref_temp
        READ(numin,*) ref_temp
        PRINT *, '  => reference temperature for heat transport is ', ref_temp
     ELSE
        BACKSPACE(numin)
     END IF
     
     WRITE(crt,'(f7.4)') ref_temp


     ! get e3u, e3v  at all levels
     CALL GETVAR_3D(idf_0, idv_0, cf_mm, 'e3v_0', 0, 0, E3V_3D)
     idf_0 = 0. ; idv_0 = 0.
     CALL GETVAR_3D(idf_0, idv_0, cf_mm, 'e3u_0', 0, 0, E3U_3D)







     DO jt = 1, nt !lolo
        !! -------------

        PRINT *, 'jt =', jt

        !! Reading 3D fields at time jt...
        CALL GETVAR_3D(idf_u, idv_u,  cf_u, trim(cv_u), nt, jt, U_3D)
        CALL GETVAR_3D(idf_v, idv_v,  cf_v, trim(cv_v), nt, jt, V_3D)
        
        IF ( leiv ) THEN
           CALL GETVAR_3D(idf_ueiv, idv_ueiv,  cf_u, trim(cv_ueiv), nt, jt, UEIV_3D)
           CALL GETVAR_3D(idf_veiv, idv_veiv,  cf_v, trim(cv_veiv), nt, jt, VEIV_3D)
        END IF

        

        CALL GETVAR_3D(idf_ut, idv_ut, cf_vt, 'vozout', nt, jt, UT_3D)
        CALL GETVAR_3D(idf_vt, idv_vt, cf_vt, 'vomevt', nt, jt, VT_3D)
        CALL GETVAR_3D(idf_us, idv_us, cf_vt, 'vozous', nt, jt, US_3D)
        CALL GETVAR_3D(idf_vs, idv_vs, cf_vt, 'vomevs', nt, jt, VS_3D)



        ! look for nearest level to imeter
        ik = 1

        ilev0(1)      = 1
        ilev1(nclass) = npk-1

        DO jk = 1, nclass -1
           DO WHILE ( gdepw(ik)  < imeter(jk) )
              ik = ik +1
           END DO

           rd1= ABS(gdepw(ik-1) - imeter(jk) )
           rd2= ABS(gdepw(ik) - imeter(jk) )
           IF ( rd2 < rd1 ) THEN
              ilev1(jk) = ik -1  ! t-levels
              ilev0(jk+1) = ik
           ELSE
              ilev1(jk) = ik -2  ! t-levels
              ilev0(jk+1) = ik -1
           END IF
        END DO



        !! compute the transport
        ztrpu (:,:,:)= 0
        ztrpv (:,:,:)= 0

        ztrput(:,:,:)= 0
        ztrpvt(:,:,:)= 0

        ztrpus(:,:,:)= 0
        ztrpvs(:,:,:)= 0

        DO jc = 1, nclass
           
           !PRINT *, '' ; PRINT *, 'LOLO: ilev0,ilev1 => ', ilev0(jc),ilev1(jc)
           !!PRINT *, '   W-depths  => ', gdepw(ilev0(jc)+1), '=>', gdepw(ilev1(jc)+1)
           !PRINT *, '   T-depths  => ', gdept(ilev0(jc)), '=>', gdept(ilev1(jc))
           
           DO jk = ilev0(jc),ilev1(jc)

              ! Get velocities, temperature and salinity fluxes at jk
              IF ( ltest ) THEN
                 zu (:,:)= udum
                 zv (:,:)= vdum
                 zut(:,:)= udum
                 zvt(:,:)= vdum
                 zus(:,:)= udum
                 zvs(:,:)= vdum
              ELSE
                 zu (:,:) =  U_3D(:,:,jk)
                 zv (:,:) =  V_3D(:,:,jk)
                 IF ( leiv ) THEN
                    zu(:,:) = zu(:,:) + UEIV_3D(:,:,jk)
                    zv(:,:) = zv(:,:) + VEIV_3D(:,:,jk)
                 END IF

                 zut(:,:) = UT_3D(:,:,jk)
                 zvt(:,:) = VT_3D(:,:,jk)
                 zus(:,:) = US_3D(:,:,jk)
                 zvs(:,:) = VS_3D(:,:,jk)

              ENDIF

              e3v(:,:) = E3V_3D(:,:,jk)
              e3u(:,:) = E3U_3D(:,:,jk)

              zwku (:,:) = zu (:,:)*e2u(:,:)*e3u(:,:)
              zwkv (:,:) = zv (:,:)*e1v(:,:)*e3v(:,:)
              zwkut(:,:) = zut(:,:)*e2u(:,:)*e3u(:,:)
              zwkvt(:,:) = zvt(:,:)*e1v(:,:)*e3v(:,:)
              zwkus(:,:) = zus(:,:)*e2u(:,:)*e3u(:,:)
              zwkvs(:,:) = zvs(:,:)*e1v(:,:)*e3v(:,:)

              ! integrates vertically
              ztrpu (:,:,jc) = ztrpu (:,:,jc) + zwku (:,:)
              ztrpv (:,:,jc) = ztrpv (:,:,jc) + zwkv (:,:)
              ztrput(:,:,jc) = ztrput(:,:,jc) + zwkut(:,:) * rau0*rcp
              ztrpvt(:,:,jc) = ztrpvt(:,:,jc) + zwkvt(:,:) * rau0*rcp

              ! At the end instead
              !ztrput(:,:,jc) = ztrput(:,:,jc) + rau0*rcp*(zwkut(:,:) - ref_temp*zu (:,:)*e2u(:,:)*e3u(:,:))
              !ztrpvt(:,:,jc) = ztrpvt(:,:,jc) + rau0*rcp*(zwkvt(:,:) - ref_temp*zv (:,:)*e1v(:,:)*e3v(:,:))

              ztrpus(:,:,jc) = ztrpus(:,:,jc) + zwkus(:,:)
              ztrpvs(:,:,jc) = ztrpvs(:,:,jc) + zwkvs(:,:)

           END DO  ! loop to next level
        END DO    ! next class



        IF (jt == 1) THEN

           IF ( nclass > 3 ) THEN
              PRINT *, 'LOLO=> cdftransportiz.f90 only supports 2 depths currently!!!'
           END IF


           !! Find the broken line between P1 (imin,jmin) and P2 (imax, jmax)
           !! ---------------------------------------------------------------
           ! ... Initialization
           i0=imin; j0=jmin; i1=imax;  j1=jmax
           rxi1=i1;  ryj1=j1; rxi0=i0; ryj0=j0

           ! .. Compute equation:  ryj = aj rxi + bj
           IF ( (rxi1 -rxi0) /=  0 ) THEN
              aj = (ryj1 - ryj0 ) / (rxi1 -rxi0)
              bj = ryj0 - aj * rxi0
           ELSE
              aj=10000.
              bj=0.
           END IF

           ! .. Compute equation:  rxi = ai ryj + bi
           IF ( (ryj1 -ryj0) /=  0 ) THEN
              ai = (rxi1 - rxi0 ) / ( ryj1 -ryj0 )
              bi = rxi0 - ai * ryj0
           ELSE
              ai=10000.
              bi=0.
           END IF

           ! ..  Compute the integer pathway:
           n=0
           ! .. Chose the strait line with the smallest slope
           IF (ABS(aj) <=  1 ) THEN
              ! ... Here, the best line is y(x)
              ! ... If i1 < i0 swap points and remember it has been swapped
              IF (i1 <  i0 ) THEN
                 i  = i0 ; j  = j0
                 i0 = i1 ; j0 = j1
                 i1 = i  ; j1 = j
              END IF

              IF ( j1 >= j0 ) THEN
                 ist = 1     ; jst = 1
                 norm_u =  1 ;  norm_v = -1
              ELSE
                 ist = 1     ; jst = 0
                 norm_u = -1 ; norm_v = -1
              END IF

              ! ... compute the nearest j point on the line crossing at i
              DO i=i0,i1
                 n=n+1
                 IF (n > jpseg) STOP 'n > jpseg !'
                 j=NINT(aj*i + bj )
                 yypt(n) = CMPLX(i,j)
              END DO
           ELSE
              ! ... Here, the best line is x(y)
              ! ... If j1 < j0 swap points and remember it has been swapped
              IF (j1 <  j0 ) THEN
                 i  = i0 ; j  = j0
                 i0 = i1 ; j0 = j1
                 i1 = i  ; j1 = j
              END IF
              IF ( i1 >=  i0 ) THEN
                 ist = 1    ;  jst = 1
                 norm_u = 1 ;  norm_v = -1
              ELSE
                 ist = 0
                 jst = 1
                 norm_u = 1
                 norm_v = 1
              END IF

              ! ... compute the nearest i point on the line crossing at j
              DO j=j0,j1
                 n=n+1
                 IF (n > jpseg) STOP 'n>jpseg !'
                 i=NINT(ai*j + bi)
                 yypt(n) = CMPLX(i,j)
              END DO
           END IF

           !!
           !! Look for intermediate points to be added.
           !  ..  The final positions are saved in rxx,ryy
           rxx(1)=REAL(yypt(1))
           ryy(1)=IMAG(yypt(1))
           nn=1

           DO k=2,n
              ! .. distance between 2 neighbour points
              d=ABS(yypt(k)-yypt(k-1))
              ! .. intermediate points required if d > 1
              IF ( d > 1 ) THEN
                 CALL interm_pt(yypt,k,ai,bi,aj,bj,yypti)
                 nn=nn+1
                 IF (nn > jpseg) STOP 'nn>jpseg !'
                 rxx(nn)=REAL(yypti)
                 ryy(nn)=IMAG(yypti)
              END IF
              nn=nn+1
              IF (nn > jpseg) STOP 'nn>jpseg !'
              rxx(nn)=REAL(yypt(k))
              ryy(nn)=IMAG(yypt(k))
           END DO


           ! Now extract the transport through a section
           ! ... Check whether we need a u velocity or a v velocity
           !   Think that the points are f-points and delimit either a U segment
           !   or a V segment (ist and jst are set in order to look for the correct
           !   velocity point on the C-grid
           !PRINT *, TRIM(csection)
           !PRINT *, 'IMIN IMAX JMIN JMAX', imin, imax, jmin, jmax
           !!


           IF ( l_save_broken_lines ) THEN
              !! LOLO: We want to save the Broken line in an ascci file for further use:
              WRITE(cf_broken_line,'("broken_line_",a,".dat")') trim(csection)
              OPEN(UNIT=19, FILE=trim(cf_broken_line), FORM='FORMATTED', RECL=256, STATUS='unknown')
              PRINT *, 'Saving broken line into file ', trim(cf_broken_line)
              WRITE(19,*)'#       ji          jj'
              DO jseg = 1, nn
                 i0=rxx(jseg)
                 j0=ryy(jseg)
                 WRITE(19,*) i0, j0
              END DO
              CLOSE(19)
              PRINT *, ''
           END IF

        END IF !* IF ( jt == 1 )





        DO jc=1,nclass

           voltrpsum = 0.
           heatrpsum = 0.
           saltrpsum = 0.

           DO jseg = 1, nn-1
              i0=rxx(jseg)
              j0=ryy(jseg)
              IF ( rxx(jseg) ==  rxx(jseg+1) ) THEN
                 gla(jseg)=glamu(i0,j0+jst)   ; gphi(jseg)=gphiu(i0,j0+jst)
                 voltrp(jseg)= ztrpu (i0,j0+jst,jc)*norm_u
                 heatrp(jseg)= ztrput(i0,j0+jst,jc)*norm_u
                 saltrp(jseg)= ztrpus(i0,j0+jst,jc)*norm_u
              ELSE IF ( ryy(jseg) == ryy(jseg+1) ) THEN
                 gla(jseg)=glamv(i0+ist,j0)  ;  gphi(jseg)=gphiv(i0+ist,j0)
                 voltrp(jseg)=ztrpv (i0+ist,j0,jc)*norm_v
                 heatrp(jseg)=ztrpvt(i0+ist,j0,jc)*norm_v
                 saltrp(jseg)=ztrpvs(i0+ist,j0,jc)*norm_v
              ELSE
                 PRINT *,' ERROR :',  rxx(jseg),ryy(jseg),rxx(jseg+1),ryy(jseg+1)
              END IF
              voltrpsum = voltrpsum+voltrp(jseg)
              heatrpsum = heatrpsum+heatrp(jseg)
              saltrpsum = saltrpsum+saltrp(jseg)
           END DO   ! next segment

           !IF ( jt == 1 )  THEN
           !   WRITE(numout+jc,*)  '# LONmin LATmin LONmax LATmax = ', &
           !        &                 gla(1),gphi(1), gla(nn-1), gphi(nn-1)
           !   WRITE(numout+jc,*)  '# Top(m)  Bottom(m) =', gdepw(ilev0(jc)), gdepw(ilev1(jc)+1)
           !   WRITE(numout+jc,*)  '#'
           !   WRITE(numout+jc,*)  '# Time rec.  MassTrans(Sv) HeatTrans(PW) SaltTrans(kt/s)'
           !END IF
           !!lolo
           !WRITE(numout+jc,9002) jt, voltrpsum/1.e6, heatrpsum/1.e15, saltrpsum/1.e6


           heatrpsum = heatrpsum - rau0*rcp*ref_temp*voltrpsum

           X_trsp(:,jt,jc) = (/ REAL(voltrpsum/1.e6,4), REAL(heatrpsum/1.e15,4), REAL(saltrpsum/1.e6,4) /) ! lolo


        END DO ! next class

     END DO ! loop on jt





     !! Time to create or append X_trsp into the netcdf file for current section
     !! -----------------------------------------------------
     
     WRITE(cf_out, '(a,"/transport_sect_",a,".nc")') trim(cd_out), trim(csection)
     IF ( leiv ) WRITE(cf_out, '(a,"/transport_sect_",a,"_eiv.nc")') trim(cd_out), trim(csection)
     
     id_volu = 0 ; id_heat = 0 ; id_salt = 0

     !! LOLO netcdf
     INQUIRE( FILE=cf_out, EXIST=lfncout )


     IF ( .NOT. lfncout ) THEN

        !! Creating file
        PRINT *, ' Creating file '//trim(cf_out)//' !!!'
        ierr = NF90_CREATE(cf_out, NF90_CLOBBER, idf_out)
        ierr = NF90_DEF_DIM(idf_out, 'time', NF90_UNLIMITED, idd_t)
        ierr = NF90_DEF_VAR(idf_out, 'time', NF90_DOUBLE,    idd_t, idv_time)


        ierr = NF90_DEF_VAR(idf_out, 'trsp_volu', NF90_FLOAT, (/idd_t/), id_volu(0))
        ierr = NF90_DEF_VAR(idf_out, 'trsp_heat', NF90_FLOAT, (/idd_t/), id_heat(0))
        ierr = NF90_DEF_VAR(idf_out, 'trsp_salt', NF90_FLOAT, (/idd_t/), id_salt(0))

        ierr = NF90_PUT_ATT(idf_out, id_volu(0), 'long_name', 'TOTAL: Transport of volume')
        ierr = NF90_PUT_ATT(idf_out, id_heat(0), 'long_name', 'TOTAL: Transport of heat')
        ierr = NF90_PUT_ATT(idf_out, id_salt(0), 'long_name', 'TOTAL: Transport of salt')

        ierr = NF90_PUT_ATT(idf_out, id_volu(0), 'units', 'Sv')
        ierr = NF90_PUT_ATT(idf_out, id_heat(0), 'units', 'PW')
        ierr = NF90_PUT_ATT(idf_out, id_salt(0), 'units', 'kt/s')

        IF ( nclass > 1 ) THEN
           DO jc = 1, nclass
              WRITE(cdum,'("_",i2.2)') jc ! suffix for variable_name
              ierr = NF90_DEF_VAR(idf_out, 'trsp_volu'//trim(cdum), NF90_FLOAT, (/idd_t/), id_volu(jc))
              ierr = NF90_DEF_VAR(idf_out, 'trsp_heat'//trim(cdum), NF90_FLOAT, (/idd_t/), id_heat(jc))
              ierr = NF90_DEF_VAR(idf_out, 'trsp_salt'//trim(cdum), NF90_FLOAT, (/idd_t/), id_salt(jc))
              
              WRITE(cdum,'(f7.2,"-",f7.2,"m  (t-points)")') gdept(ilev0(jc)), gdept(ilev1(jc))

              ierr = NF90_PUT_ATT(idf_out, id_volu(jc), 'long_name', 'Transport of volume, '//trim(cdum))
              ierr = NF90_PUT_ATT(idf_out, id_heat(jc), 'long_name', 'Transport of heat, '//trim(cdum))
              ierr = NF90_PUT_ATT(idf_out, id_salt(jc), 'long_name', 'Transport of salt, '//trim(cdum))

              ierr = NF90_PUT_ATT(idf_out, id_volu(jc), 'units', 'Sv')
              ierr = NF90_PUT_ATT(idf_out, id_heat(jc), 'units', 'PW')
              ierr = NF90_PUT_ATT(idf_out, id_salt(jc), 'units', 'kt/s')
           END DO
        END IF

        ierr = NF90_PUT_ATT(idf_out, NF90_GLOBAL, 'Info', 'Reference temperature for heat transport is '//trim(crt)//' deg.C')
        ierr = NF90_PUT_ATT(idf_out, NF90_GLOBAL, 'About', 'Created by BaraKuda (cdftransportiz.f90), contact: brodeau@gmail.com')

        ierr = NF90_ENDDEF(idf_out)
        jt_pos = 0

     ELSE

        !! Opening already existing file
        ierr = NF90_OPEN  (cf_out, NF90_WRITE,   idf_out)

        !! Need IDs of variables to append... NF90_INQ_VARID
        ierr = NF90_INQ_VARID(idf_out, 'trsp_volu', id_volu(0))
        ierr = NF90_INQ_VARID(idf_out, 'trsp_heat', id_heat(0))
        ierr = NF90_INQ_VARID(idf_out, 'trsp_salt', id_salt(0))

        IF ( nclass > 1 ) THEN
           DO jc = 1, nclass
              WRITE(cdum,'("_",i2.2)') jc ! suffix for variable_name
              ierr = NF90_INQ_VARID(idf_out, 'trsp_volu'//trim(cdum), id_volu(jc))
              ierr = NF90_INQ_VARID(idf_out, 'trsp_heat'//trim(cdum), id_heat(jc))
              ierr = NF90_INQ_VARID(idf_out, 'trsp_salt'//trim(cdum), id_salt(jc))
           END DO
        END IF

        ! Get ID of unlimited dimension
        ierr = NF90_INQUIRE(idf_out, unlimitedDimId = idv_time)

        ! Need to know jt_pos, record number of the last time record writen in the file
        ierr = NF90_INQUIRE_DIMENSION(idf_out, idv_time, name=cv_dum, len=jt_pos)

     END IF

     WRITE(*,'("Going to write record ",i4.4," to ",i4.4," into ",a)') jt_pos+1, jt_pos+1+nt, trim(cf_out)

     DO jt = 1, nt

        ! Writing record jt for time vector and 1d fields:
        ierr = NF90_PUT_VAR( idf_out, idv_time, (/ryear+1./12.*(REAL(jt)-1.+0.5)/), start=(/jt_pos+jt/), count=(/1/) )



        IF ( nclass == 1 ) THEN
           !! Default variable is the only one present (index = 1) :
           ierr = NF90_PUT_VAR(idf_out, id_volu(0), (/ X_trsp(1,jt,1) /), start=(/jt_pos+jt/), count=(/1/))
           ierr = NF90_PUT_VAR(idf_out, id_heat(0), (/ X_trsp(2,jt,1) /), start=(/jt_pos+jt/), count=(/1/))
           ierr = NF90_PUT_VAR(idf_out, id_salt(0), (/ X_trsp(3,jt,1) /), start=(/jt_pos+jt/), count=(/1/))
           
        ELSE

           !! Default variable is the sum:
           ierr = NF90_PUT_VAR(idf_out, id_volu(0), (/ SUM(X_trsp(1,jt,:)) /), start=(/jt_pos+jt/), count=(/1/))
           ierr = NF90_PUT_VAR(idf_out, id_heat(0), (/ SUM(X_trsp(2,jt,:)) /), start=(/jt_pos+jt/), count=(/1/))
           ierr = NF90_PUT_VAR(idf_out, id_salt(0), (/ SUM(X_trsp(3,jt,:)) /), start=(/jt_pos+jt/), count=(/1/))

           DO jc = 1, nclass
              ierr = NF90_PUT_VAR(idf_out, id_volu(jc), (/ X_trsp(1,jt,jc) /), start=(/jt_pos+jt/), count=(/1/))
              ierr = NF90_PUT_VAR(idf_out, id_heat(jc), (/ X_trsp(2,jt,jc) /), start=(/jt_pos+jt/), count=(/1/))
              ierr = NF90_PUT_VAR(idf_out, id_salt(jc), (/ X_trsp(3,jt,jc) /), start=(/jt_pos+jt/), count=(/1/))
           END DO
        END IF

     END DO

     ierr = NF90_CLOSE(idf_out)

     PRINT *, ''





  END DO ! loop on sections...






  DEALLOCATE( X_trsp ) !lolo


CONTAINS


  SUBROUTINE interm_pt (ydpt,k,pai,pbi,paj,pbj,ydpti)
    !! -----------------------------------------------------
    !!           SUBROUTINE INTERM_PT
    !!           ********************
    !!
    !!   PURPOSE:
    !!   --------
    !!     Find the best intermediate points on a pathway.
    !!
    !!    ARGUMENTS:
    !!    ----------
    !!      ydpt : complex vector of the positions of the nearest points
    !!         k : current working index
    !!       pai ,pbi : slope and original ordinate of x(y)
    !!       paj ,pbj : slope and original ordinate of y(x)
    !!      ydpti : Complex holding the position of intermediate point
    !!
    !!    AUTHOR:
    !!    -------
    !!      19/07/1999 : Jean-Marc MOLINES
    !!      14/01/2005 : J M M in F90
    !!
    !!--------------------------------------------------------------
    !!
    !! 0. Declarations:
    !! ----------------
    IMPLICIT NONE
    COMPLEX, INTENT(in) :: ydpt(*)
    COMPLEX, INTENT(out) :: ydpti
    REAL(4), INTENT(IN) ::  pai,pbi,paj,pbj
    INTEGER ,INTENT(in) :: k
    ! ... local
    COMPLEX :: ylptmp1, ylptmp2
    REAL(4) ::  za0,zb0,za1,zb1,zd1,zd2
    REAL(4) ::  zxm,zym
    REAL(4) ::  zxp,zyp
    !!
    !! 1. Compute intermediate points
    !! ------------------------------
    !
    ! ... Determines whether we use y(x) or x(y):
    IF (ABS(paj) <=  1) THEN
       ! ..... y(x)
       ! ... possible intermediate points:
       ylptmp1=ydpt(k-1)+(1.,0.)
       ylptmp2=ydpt(k-1)+CMPLX(0.,SIGN(1.,paj))
       !
       ! ... M is the candidate point:
       zxm=REAL(ylptmp1)
       zym=IMAG(ylptmp1)
       za0=paj
       zb0=pbj
       !
       za1=-1./za0
       zb1=zym - za1*zxm
       ! ... P is the projection of M on the strait line
       zxp=-(zb1-zb0)/(za1-za0)
       zyp=za0*zxp + zb0
       ! ... zd1 is the distance MP
       zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       !
       ! ... M is the candidate point:
       zxm=REAL(ylptmp2)
       zym=IMAG(ylptmp2)
       za1=-1./za0
       zb1=zym - za1*zxm
       ! ... P is the projection of M on the strait line
       zxp=-(zb1-zb0)/(za1-za0)
       zyp=za0*zxp + zb0
       ! ... zd2 is the distance MP
       zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       ! ... chose the smallest (zd1,zd2)
       IF (zd2 <=  zd1) THEN
          ydpti=ylptmp2
       ELSE
          ydpti=ylptmp1
       END IF
       !
    ELSE
       !
       ! ... x(y)
       ylptmp1=ydpt(k-1)+CMPLX(SIGN(1.,pai),0.)
       ylptmp2=ydpt(k-1)+(0.,1.)
       zxm=REAL(ylptmp1)
       zym=IMAG(ylptmp1)
       za0=pai
       zb0=pbi
       !
       za1=-1./za0
       zb1=zxm - za1*zym
       zyp=-(zb1-zb0)/(za1-za0)
       zxp=za0*zyp + zb0
       zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       !
       zxm=REAL(ylptmp2)
       zym=IMAG(ylptmp2)
       za1=-1./za0
       zb1=zxm - za1*zym
       zyp=-(zb1-zb0)/(za1-za0)
       zxp=za0*zyp + zb0
       zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       IF (zd2 <=  zd1) THEN
          ydpti=ylptmp2
       ELSE
          ydpti=ylptmp1
       END IF
    END IF
  END SUBROUTINE interm_pt

END PROGRAM cdftransportiz
