PROGRAM cdfpsi
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfpsi  ***
  !!
  !!  **  Purpose  :  Compute Barotropic Stream Function
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  Compute the 2D fields ztrpu, ztrpv 
  !!                  as the integral on the vertical of u, v on their
  !!                  respective points. 
  !!                  Then integrate from south to north : ==> psiu
  !!                  Then integrate from West to East   : ==> psiv
  !!                  (should be almost the same (if no error ))
  !!   Default (appropriate for global model): output psiu;
  !!                    normalizes the values setting psi (jpi,jpj) = 0
  !!   If option "V" is given as last argument, output psiv, 
  !!                    normalizes values setting psi(jpi,1) = 0.
  !!                    This is appropriate for North Atlantic 
  !!
  !! history ;
  !!  Original :  J.M. Molines (May 2005 )
  !!-------------------------------------------------------------------
  !!  $Rev: 256 $
  !!  $Date: 2009-07-21 17:49:27 +0200 (Tue, 21 Jul 2009) $
  !!  $Id: cdfpsi.f90 256 2009-07-21 15:49:27Z molines $
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE io_ezcdf

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk,jt                         !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk, nt              !: size of the domain
  INTEGER   :: ncout
  INTEGER, DIMENSION(1) ::  ipk, id_varout         !

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask, e1v, zv !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e2u, zu !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         glamf, gphif
  REAL(KIND=4) ,DIMENSION(:),      ALLOCATABLE ::  tim

  REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: E3U_3D, E3V_3D


  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: ztrpu, ztrpv, psiu, psiv

  CHARACTER(LEN=256) :: cfileu ,cfilev, cfileoutnc='psi.nc', ctim
  CHARACTER(LEN=256) :: cf_mm='mesh_mask.nc'
  CHARACTER(LEN=1)  :: coption

  TYPE(variable), DIMENSION(1)  :: typvar         !: structure for attributes

  INTEGER    :: istatus

  INTEGER :: idf_0=0, idv_0=0

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfpsi  Ufile Vfile <V> (optional argument)'
     PRINT *,' Computes the barotropic stream function as the integral of the transport'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc must be in te current directory'
     PRINT *,' Output on psi.nc, variables sobarstf on f-points'
     PRINT *,' Default works well for a global ORCA grid. use V 3rdargument for North Atlantic'
     STOP
  ENDIF

  CALL getarg (1, cfileu  )
  CALL getarg (2, cfilev  )
  CALL getarg (3, coption )

  npiglo= getdim (cfileu,'x')
  npjglo= getdim (cfileu,'y')
  npk   = getdim (cfileu,'depth')
  
  ctim = 'none'
  nt    = getdim (cfileu,'time',cdtrue=ctim,kstatus=istatus) !LB

 ! define new variables for output ( must update att.txt)
  typvar(1)%name= 'sobarstf'
  typvar(1)%units='m3/s'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -300.e6
  typvar(1)%valid_max= 300.e6
  typvar(1)%long_name='Barotropic_Stream_Function'
  typvar(1)%short_name='sobarstf'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'
  ipk(1) = 1  !  2D ( X, Y , T )

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'nt   =', nt !LB

  
  IF ( coption == 'V') PRINT *, ' Use psiv (ex. North Atlantic case)'

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo), e2u(npiglo,npjglo), E3U_3D(npiglo,npjglo,npk), E3V_3D(npiglo,npjglo,npk))
  ALLOCATE ( zu(npiglo,npjglo),ztrpu(npiglo,npjglo), psiu(npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo),ztrpv(npiglo,npjglo), psiv(npiglo,npjglo))
  ALLOCATE ( glamf(npiglo,npjglo), gphif(npiglo,npjglo))
  ALLOCATE ( tim(nt) ) !LB

  glamf(:,:) = getvar(cf_mm, 'glamf',1,npiglo,npjglo)
  gphif(:,:) = getvar(cf_mm, 'gphif',1,npiglo,npjglo)

  ! create output fileset
  ncout =create(cfileoutnc, cfileu, npiglo,npjglo,1)
  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfileu,npiglo, npjglo,1,glamf, gphif)
  tim=getvar1d(cfileu,trim(ctim),nt) !LB
  ierr=putvar1d(ncout,tim,nt,'T') !LB


  e1v(:,:)   = getvar(cf_mm, 'e1v', 1,npiglo,npjglo)
  e2u(:,:)   = getvar(cf_mm, 'e2u', 1,npiglo,npjglo)
  zmask(:,:) = getvar(cf_mm, 'fmask', 1,npiglo,npjglo)



  !lolo
  ! get e3u, e3v  at all levels
  CALL GETVAR_3D(idf_0, idv_0, cf_mm, 'e3u_0', 0, 0, E3U_3D) ; idf_0 = 0. ; idv_0 = 0.
  CALL GETVAR_3D(idf_0, idv_0, cf_mm, 'e3v_0', 0, 0, E3V_3D) ; idf_0 = 0. ; idv_0 = 0.

  ! get rid of the free-slip/no-slip condition
  WHERE ( zmask >= 2 ) zmask = 1



  DO jt=1,nt
     
     PRINT *, 'jt =', jt

     ztrpu(:,:)= 0.d0
     ztrpv(:,:)= 0.d0
     
     DO jk = 1,npk
        PRINT *,'level ',jk
        IF ( coption == 'V' ) THEN
           zv(:,:)= getvar(cfilev, 'vomecrty',  jk ,npiglo,npjglo,  ktime=jt) !LB
           !loloe3v(:,:) = getvar(cf_mm, 'e3v_ps', jk,npiglo,npjglo, ldiom=.true.)
           ztrpv(:,:) = ztrpv(:,:) + zv(:,:)*e1v(:,:)*E3V_3D(:,:,jk)  ! meridional transport of each grid cell
        ELSE
           ! Get zonal velocity  at jk
           zu(:,:)= getvar(cfileu, 'vozocrtx',  jk ,npiglo,npjglo,  ktime=jt) !LB
           ! get e3v at level jk
           !loloe3u(:,:) = getvar(cf_mm, 'e3u_ps', jk,npiglo,npjglo, ldiom=.true.)
           ! integrates vertically 
           ztrpu(:,:) = ztrpu(:,:) + zu(:,:)*e2u(:,:)*E3U_3D(:,:,jk)  ! zonal transport of each grid cell
        ENDIF
     END DO  ! loop to next level
     
     IF (coption == 'V' ) THEN
        ! integrate zonally from east to west 
        psiv(npiglo,:)= 0.0
        DO ji=npiglo-1,1,-1
           psiv(ji,:) = psiv(ji+1,:) - ztrpv(ji,:)  ! psi at f point
        END DO
        psiv(:,:) = psiv(:,:) *zmask(:,:)
        ierr = putvar(ncout, id_varout(1) ,SNGL(psiv),   1, npiglo, npjglo, jt) !LB
        
     ELSE
        ! integrate from the south to the north with zonal transport
        psiu(:,:) = 0.d0
        
        DO jj = 2, npjglo
           psiu(:,jj) = psiu(:,jj-1) - ztrpu(:,jj)   ! psi at f point
        END DO
        psiu(:,:) = (psiu(:,:) -psiu(npiglo,npjglo) ) * zmask(:,:)
        ierr = putvar(ncout, id_varout(1) ,SNGL(psiu),   1, npiglo, npjglo, jt) !LB
     ENDIF
     
  END DO

   istatus = closeout (ncout)

   END PROGRAM cdfpsi
   
