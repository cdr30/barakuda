PROGRAM cdfmoc
   !!-------------------------------------------------------------------
   !!               ***  PROGRAM cdfmoc  ***
   !!
   !!  **  Purpose  :  Compute the Meridional Overturning Cell (MOC)
   !!                  PARTIAL STEPS
   !!
   !!  **  Method   :  The MOC is computed from the V velocity field, integrated
   !!                  from the bottom to the surface, then zonally averaged with
   !!                  eventual masking for oceanic basins.
   !!                  The program looks for the file "new_maskglo.nc". If it does not exist,
   !!                  only the calculation over all the domain is performed (this is adequate
   !!                  for a basin configuration like NATL4).
   !!                  In new_maskglo.nc the masking corresponds to the global
   !!                  configuration. MOC for Global, Atlantic, Indo-Pacific, Indian,Pacific ocean
   !!                  Results are saved on moc.nc file with variables name respectively
   !!                  zomsfglo, zomsfatl, zomsfinp, zomsfind, zomsfpac
   !!
   !!
   !! history ;
   !!  Original :  J.M. Molines  (jul.  2005)
   !!              A.M. Treguier (april 2006) adaptation to NATL4 case
   !!-------------------------------------------------------------------
   !!  $Rev: 256 $
   !!  $Date: 2009-07-21 17:49:27 +0200 (Tue, 21 Jul 2009) $
   !!  $Id: cdfmoc.f90 256 2009-07-21 15:49:27Z molines $
   !!--------------------------------------------------------------
   !! * Modules used

   USE cdfio
   USE io_ezcdf

   !! * Local variables
   IMPLICIT NONE
   INTEGER   :: jpbasins
   INTEGER   :: jbasin, jj, jk ,ji, jt              !: dummy loop index
   INTEGER   :: ierr                                !: working integer
   INTEGER   :: narg, iargc                         !: command line
   INTEGER   :: npiglo,npjglo, npk, nt              !: size of the domain
   INTEGER   :: ncout, np
   INTEGER   :: numout=10
   INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout         !
   INTEGER, DIMENSION(2)              ::  iloc

   !REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  e1v, e3v, gphiv, zv !:  metrics, velocity
   REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  e1v, gphiv, zv !:  metrics, velocity
   REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  e3v, V_3D, Veiv_3D
   REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
   REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat              !: latitude for i = north pole
   REAL(KIND=4), DIMENSION (:),       ALLOCATABLE ::  gdepw               !: deptw
   REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  zmask               !:  jpbasins x npiglo x npjglo
   REAL(KIND=4), DIMENSION (:),       ALLOCATABLE ::  tim       !LB

   REAL(KIND=8) ,DIMENSION(:,:,:) ,   ALLOCATABLE ::  zomsf                 !: jpbasins x npjglo x npk

   CHARACTER(LEN=256) :: cf_v , cfileoutnc='moc.nc', ctim, cv_v, cv_veiv

   CHARACTER(LEN=256) :: cf_mm='mesh_mask.nc', cf_bm='new_maskglo.nc'
   TYPE(variable)    ,DIMENSION(:), ALLOCATABLE   :: typvar                   !: structure for attribute

   LOGICAL    :: llglo = .false., & !: indicator for presence of new_maskglo.nc file
      &        leiv  = .FALSE.    !: weather to use Eddy Induced velocity from GM90

   INTEGER    :: istatus

   INTEGER :: idf_0, idv_0, idf_v, idv_v, idf_veiv, idv_veiv

   ! constants

   !!  Read command line and output usage message if not compliant.
   narg= iargc()
   IF ( (narg < 2).OR.(narg > 3) ) THEN
      PRINT *,' Usage : cdfmoc  <V file> <V variable> (<Veiv variable>)'
      PRINT *,' Computes the MOC for oceanic basins as described in new_maskglo.nc'
      PRINT *,' PARTIAL CELLS VERSION'
      PRINT *,' Files mesh_mask.nc and new_maskglo.nc '
      PRINT *,'  must be in the current directory'
      PRINT *,' Output on moc.nc: '
      PRINT *,'      variables zomsfglo  : Global ocean '
      PRINT *,'      variables zomsfatl  : Atlantic Ocean '
      PRINT *,'      variables zomsfinp  : Indo Pacific '
      PRINT *,'      variables zomsfind  : Indian Ocean alone'
      PRINT *,'      variables zomsfpac  : Pacific Ocean alone'
      STOP
   ENDIF

   CALL getarg (1, cf_v)
   CALL getarg (2, cv_v)

   PRINT *, ' Will compute MOC using '//trim(cv_v)//' !!!'

   leiv = .FALSE.
   IF (narg == 3) THEN
      CALL getarg (3, cv_veiv)
      IF ( trim(cv_veiv) /= '0' ) THEN
         leiv = .TRUE.
         PRINT *, ' and taking eddy-induced velocity into account using '//trim(cv_veiv)//' !!!'
      END IF
   END IF


   npiglo= getdim(cf_v,'x',ldexact=.TRUE.)
   npjglo= getdim(cf_v,'y',ldexact=.TRUE.)
   npk   = getdim(cf_v,'depthv')

   ctim = 'none'
   nt    = getdim (cf_v,'time',cdtrue=ctim,kstatus=istatus) !LB

   !LB:
   PRINT *, ''
   PRINT *, 'npiglo =',npiglo
   PRINT *, 'npjglo =',npjglo
   PRINT *, 'npk    =',npk
   PRINT *, 'nt     =',nt  !PM
   PRINT *, ''
   !LB.

   !LB:
   IF (nt == 0) THEN
      PRINT *, 'nt=0, assume 1' ; nt = 1
   END IF
   !LB.




   !  Detects newmaskglo file
   INQUIRE( FILE=cf_bm, EXIST=llglo )
   IF (llglo) THEN
      jpbasins = 5
      PRINT *, 'Basin file "new_maskglo.nc" found!'; PRINT *, ''
   ELSE
      jpbasins = 1
   ENDIF

   ALLOCATE ( typvar(jpbasins), ipk(jpbasins), id_varout(jpbasins) )

   ! define new variables for output
   typvar(1)%name= 'zomsfglo'
   typvar%units='Sverdrup'
   typvar%missing_value=99999.
   typvar%valid_min= -1000.
   typvar%valid_max= 1000.
   typvar%scale_factor= 1.
   typvar%add_offset= 0.
   typvar%savelog10= 0.
   typvar(1)%long_name='Meridional_Overt.Cell_Global'
   typvar(1)%short_name='zomsfglo'
   typvar%online_operation='N/A'
   typvar%axis='TZY'

   ipk(1) = npk  !  2D

   IF (llglo) THEN
      typvar(2)%name= 'zomsfatl'
      typvar(2)%long_name='Meridional_Overt.Cell_Atlantic'
      typvar(2)%short_name='zomsfatl'

      typvar(3)%name= 'zomsfinp'
      typvar(3)%long_name='Meridional_Overt.Cell_IndoPacif'
      typvar(3)%short_name='zomsfinp'

      typvar(4)%name= 'zomsfind'
      typvar(4)%long_name='Meridional_Overt.Cell_Indian'
      typvar(4)%short_name='zomsfind'

      typvar(5)%name= 'zomsfpac'
      typvar(5)%long_name='Meridional_Overt.Cell_pacif'
      typvar(5)%short_name='zomspac'

      ipk(2) = npk  !  2D
      ipk(3) = npk  !  2D
      ipk(4) = npk  !  2D
      ipk(5) = npk  !  2D
   ENDIF


   ! Allocate arrays
   ALLOCATE ( zmask(jpbasins,npiglo,npjglo) )
   ALLOCATE ( zv(npiglo,npjglo) )
   ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo,npk), V_3D(npiglo,npjglo,npk), gphiv(npiglo,npjglo) ,gdepw(npk) )
   ALLOCATE ( zomsf(jpbasins, npjglo, npk) )
   ALLOCATE ( dumlon(1,npjglo) , dumlat(1,npjglo))
   ALLOCATE ( tim(nt) )

   IF ( leiv ) ALLOCATE ( Veiv_3D(npiglo,npjglo,npk) )


   e1v(:,:)   = getvar  (cf_mm, 'e1v', 1,npiglo,npjglo)
   gphiv(:,:) = getvar  (cf_mm, 'gphiv', 1,npiglo,npjglo)
   !gdepw(:)  = getvare3(cf_mm, 'gdepw',npk)
   gdepw(:)   = getvar1d(cf_mm, 'gdepw_1d',npk)
   gdepw(:)   = -1.*  gdepw(:)



   iloc=maxloc(gphiv)
   dumlat(1,:) = gphiv(iloc(1),:)
   dumlon(:,:) = 0.   ! set the dummy longitude to 0

   ! create output fileset
   ncout =create(cfileoutnc, 'none', 1,npjglo,npk,cdep='depthw')
   ierr= createvar(ncout ,typvar,jpbasins, ipk,id_varout )
   ierr= putheadervar(ncout, cf_v,1, npjglo,npk,pnavlon=dumlon,pnavlat=dumlat,pdep=gdepw)

   tim=getvar1d(cf_v,trim(ctim),nt) !LB: nt
   ierr=putvar1d(ncout,tim,nt,'T') !LB: nt


   ! reading the masks
   ! 1 : global ; 2 : Atlantic ; 3 : Indo-Pacif ; 4 : Indian ; 5 : Pacif
   zmask(1,:,:)=getvar(cf_mm,'vmask',1,npiglo,npjglo)
   IF ( llglo ) THEN
      zmask(2,:,:)=getvar(cf_bm,'tmaskatl',1,npiglo,npjglo)
      zmask(4,:,:)=getvar(cf_bm,'tmaskind',1,npiglo,npjglo)
      zmask(5,:,:)=getvar(cf_bm,'tmaskpac',1,npiglo,npjglo)
      zmask(3,:,:)=zmask(5,:,:)+zmask(4,:,:)
      ! ensure that there are no overlapping on the masks
      WHERE(zmask(3,:,:) > 0 ) zmask(3,:,:) = 1
      ! change global mask for GLOBAL periodic condition
      zmask(1,1,:) = 0.
      zmask(1,npiglo,:) = 0.
   ENDIF




   CALL GETVAR_3D(idf_0, idv_0, cf_mm, 'e3v_0', 0, 0, e3v) ; idf_0 = 0. ; idv_0 = 0.





   DO jt = 1, nt !LB

      PRINT *, '  * jt =', jt

      CALL GETVAR_3D(idf_v, idv_v,  cf_v, trim(cv_v), nt, jt, V_3D)

      IF ( leiv ) CALL GETVAR_3D(idf_veiv, idv_veiv,  cf_v, trim(cv_veiv), nt, jt, Veiv_3D)


      ! initialize moc to 0
      zomsf(:,:,:) = 0.

      DO jk = 1,npk-1
         ! Get velocities v at jk
         !zv(:,:) =  getvar(cf_v, trim(cv_v), jk, npiglo, npjglo, ktime=jt)  !LB
         zv(:,:) = V_3D(:,:,jk)

         IF ( leiv ) zv(:,:) = zv(:,:) + Veiv_3D(:,:,jk)

         ! integrates 'zonally' (along i-coordinate)
         DO ji=1,npiglo
            ! For all basins
            DO jbasin = 1, jpbasins
               DO jj=1,npjglo
                  zomsf(jbasin,jj,jk) = zomsf(jbasin,jj,jk) &
                     &              - e1v(ji,jj)*e3v(ji,jj,jk)*zmask(jbasin,ji,jj)*zv(ji,jj)
               ENDDO
            END DO
         END DO
      END DO

      ! integrates vertically   from bottom to surface
      DO jk=npk-1 , 1 , -1
         zomsf(:,:,jk) = zomsf(:,:,jk+1) + zomsf(:,:,jk)/1.e6
      END DO  ! loop to next level

      ! netcdf output
      DO jbasin= 1, jpbasins
         DO jk =1, npk
            ierr = putvar(ncout, id_varout(jbasin),REAL(zomsf(jbasin,:,jk)), jk,1,npjglo, jt) !LB
            !ierr = putvar(ncout, id_varout(1)     ,rotn,                     1 ,npiglo, npjglo, jt)
         END DO
      END DO


   END DO

   ierr = closeout(ncout)

END PROGRAM cdfmoc
