PROGRAM cdfw
  !!---------------------------------------------------------------------------
  !!         ***  PROGRAM  cdfw  ***
  !!
  !!  **  Purpose: Compute the 3D w for given gridU gridV files and variables
  !! 
  !!  **  Method : Use the equation on continuity: Integrate the horizontal 
  !!               divergence from bottom to the top.
  !!               ( Use the same routines than in the code )
  !!                PARTIAL STEPS 
  !!
  !! history :
  !!   Original :  J.M. Molines (June 2005)
  !!---------------------------------------------------------------------
  !!  $Rev: 256 $
  !!  $Date: 2009-07-21 17:49:27 +0200 (Tue, 21 Jul 2009) $
  !!  $Id: cdfw.f90 256 2009-07-21 15:49:27Z molines $
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: ji,jj,jk,jt                              !: dummy loop index
  INTEGER :: npiglo, npjglo, npk, nt                  !: size of the domain
  INTEGER :: narg, iargc, ncout, ierr, istatus              !: 
  INTEGER, DIMENSION(1) ::  ipk, id_varout         ! 
  INTEGER               :: itop = 1 , ibot = 2 , itmp

  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: e1v, e2u, e1t, e2t  !: metrics
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: e3u,e3v,e3t         !:   ""
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: glamt, gphit        !: longitude latitude
  REAL(kind=4), DIMENSION(:)  , ALLOCATABLE  :: gdepw               !: 
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: un, vn
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  ::  hdivn
  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: wn  !: vertical velocity on the top
  !                                              !: and bottom of a cell.
  !                                              !: wn(top) is computed
  REAL(KIND=4) ,DIMENSION(:), ALLOCATABLE  ::  tim !LB

  CHARACTER(LEN=256) :: cfilu, cfilv, ctim
  CHARACTER(LEN=256) :: chgr='mesh_hgr.nc', czgr='mesh_zgr.nc', cfileout='w.nc'
  CHARACTER(LEN=256) :: cvaru='vozocrtx', cvarv='vomecrty', cvarw='vovecrtz'

  TYPE(variable), DIMENSION(1)      :: typvar     !: structure for attributes

  !!
  narg = iargc()
  IF ( narg < 2 ) THEN
     PRINT *,' USAGE : cdfw fileU fileV [varU varV] '
     PRINT *,'        version PARTIAL STEPS '
     PRINT *,'        Produce a cdf file w.nc with vovecrtz  variable'
     PRINT *,'        Need mesh_hgr.nc mesh_zgr.nc'
     PRINT *,'        If no varU and varV variables given, assume vozocrtx, vomecrty'
     STOP
  ENDIF

  CALL getarg(1, cfilu)
  CALL getarg(2, cfilv)
  IF ( narg >2 ) THEN
     CALL getarg(3, cvaru)
     CALL getarg(4, cvarv)
  ENDIF

  npiglo = getdim(cfilu,'x')
  npjglo = getdim(cfilu,'y')
  npk    = getdim(cfilu,'depth')

  ctim = 'none'
  nt    = getdim (cfilu,'time',cdtrue=ctim,kstatus=istatus) !LB

  PRINT *, 'nk     =',npk  !LB
  PRINT *, 'nt     =',nt  !LB

  !LB:
  IF (nt == 0) THEN
     PRINT *, 'nt=0, assume 1' ; nt = 1
  END IF
  !LB.
  


  ! define new variables for output ( must update att.txt)
  typvar(1)%name= TRIM(cvarw)
  typvar(1)%units='m/s'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -1.
  typvar(1)%valid_max= 1.
  typvar(1)%long_name='Vertical_Velocity'
  typvar(1)%short_name=TRIM(cvarw)
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'

  ipk(1) = npk             !  3D

  ! Allocate the memory
  ALLOCATE ( e1v(npiglo,npjglo) ,e2u(npiglo,npjglo) )
  ALLOCATE ( e1t(npiglo,npjglo) ,e2t(npiglo,npjglo) )
  ALLOCATE ( e3u(npiglo,npjglo) ,e3v(npiglo,npjglo) ,e3t(npiglo,npjglo) )
  ALLOCATE ( glamt(npiglo,npjglo), gphit(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo) ,hdivn(npiglo,npjglo) )
  ALLOCATE ( wn(npiglo,npjglo,2) , gdepw(npk) )
  ALLOCATE ( tim(nt) ) !LB

  ! Read the metrics from the mesh_hgr file
  e2u=  getvar(chgr, 'e2u', 1,npiglo,npjglo)
  e1v=  getvar(chgr, 'e1v', 1,npiglo,npjglo)
  e1t=  getvar(chgr, 'e1t', 1,npiglo,npjglo)
  e2t=  getvar(chgr, 'e2t', 1,npiglo,npjglo)

  ! and the coordinates   from the mesh_hgr file
  glamt = getvar(chgr, 'glamt', 1,npiglo,npjglo)
  gphit = getvar(chgr, 'gphit', 1,npiglo,npjglo)

  ! Read the depth of the w points (in the file, it is not a vector but a 1x1xnpk array)
  gdepw(:) = getvare3(czgr,'gdepw',npk)

  ! create output fileset
  ncout =create(cfileout, cfilu, npiglo,npjglo,npk,cdep='depthw')
  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, 'dummy', npiglo, npjglo, npk, glamt, gphit, gdepw )
  
  tim=getvar1d(cfilu,trim(ctim),nt) !LB
  ierr=putvar1d(ncout,tim,nt,'T') !LB

  !tim=getvar1d(cfilu,'time_counter',1)
  !ierr=putvar1d(ncout,tim,1,'T')



  DO jt = 1, nt !LB
     
     PRINT *, 'jt=',jt
     
     wn(:,:,:) = 0.
     
     ! Main level loop from bottom to top
     DO jk = npk-1, 1, -1

        !print *,' jk = ', jk
        
        ! veloccities at level jk
        un(:,:) =  getvar(cfilu, cvaru, jk ,npiglo,npjglo, ktime=jt)
        vn(:,:) =  getvar(cfilv, cvarv, jk ,npiglo,npjglo, ktime=jt)
        
        ! e3 metrics at level jk ( Partial steps)
        e3u(:,:) = getvar(czgr,'e3u_0',jk ,npiglo,npjglo, ldiom=.true.) 
        e3v(:,:) = getvar(czgr,'e3v_0',jk ,npiglo,npjglo, ldiom=.true.) 
        e3t(:,:) = getvar(czgr,'e3t_0',jk ,npiglo,npjglo, ldiom=.true.) 

        
        ! Compute divergence :
        DO jj = 2, npjglo -1
           DO ji = 2, npiglo -1
              hdivn(ji,jj) =   &
                   (  e2u(ji,jj)*e3u(ji,jj) * un(ji,jj) - e2u(ji-1,jj  )*e3u(ji-1,jj)  * un(ji-1,jj ) &       
                   + e1v(ji,jj)*e3v(ji,jj) * vn(ji,jj) - e1v(ji  ,jj-1)*e3v(ji  ,jj-1)  * vn(ji  ,jj-1)  ) &
                   / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj) )
           END DO
        END DO
        
        ! Computation from the bottom
        wn(:,:,itop) = wn(:,:,ibot) - e3t(:,:) * hdivn(:,:)
        
        ! write wn  on file at level jk (This coculd be epensive at it writes from the bottom ...
        ierr = putvar(ncout, id_varout(1) ,wn(:,:,itop), jk ,npiglo, npjglo, jt) !LB
        
        ! swap top and bottom index
        itmp=itop ; itop = ibot ; ibot = itmp 
        
     ENDDO  ! loop to next level

  END DO
  ierr = closeout(ncout)
  
  
END PROGRAM cdfw

