!=================================================================
!  subroutine InitMPhys() 
!=================================================================
!
  SUBROUTINE InitMPhys(nC,Time)

     USE mo_reac
     USE mo_met
     USE mo_mphys
     USE mo_control

     IMPLICIT NONE

     INTEGER ::  nC
     REAL(8) ::  Time

!---  internal variable
     INTEGER ::  ic,it, iFrac, idra, idre
     REAL(8) ::  expo
     REAL(8) ::  radcon(ntFrac), LWClass(ntFrac)

!======================================================
!===  Set Pseudo-Microphysics
!======================================================
!
     IF (LwcAll > 0.e0)  THEN
        CALL DropDis(ntFrac,radcon,MeanRad,LWClass,LwcAll)
     ELSE IF (LwcAll == 0.e0) THEN
        IF (ntFrac <= 1)  THEN
           radcon(1) = MeanRad
        ELSE IF (ntFrac == 2)  THEN
           radcon(1) = 0.1e0 * MeanRad
           radcon(2) = MeanRad
        ELSE IF (ntFrac == 3)  THEN
           radcon(1) = 0.1e0 * MeanRad
           radcon(2) = MeanRad
           radcon(3) = 4.0e0 * MeanRad
        ELSE IF (ntFrac >= 4)  THEN
           expo = 0.e0
           MeanRad = 0.1e0 * MeanRad
           DO iFrac=1,ntFrac
              radcon(iFrac) = MeanRad * 2**expo
              expo = expo + 6.e0/(ntFrac-1)
           END DO
        END IF
        WRITE(*,*) 'InitChem ...  ntFrac =',ntFrac
        WRITE(*,*) '       Droplet Sizes =',   &
&                          (radcon(iFrac),iFrac=1,ntFrac)
     ELSE IF (LwcAll < 0.e0) THEN
        radcon(:) = MeanRad
        WRITE(*,*) 'InitChem ...  ntFrac =',ntFrac
        WRITE(*,*) '       Droplet Sizes =',   &
&                          (radcon(iFrac),iFrac=1,ntFrac)
     END IF
!
!--------------------------------------------------
!--  Set up liquid water content of fractions
!--  from meteorology (only for ntFrac=1!!)
!--------------------------------------------------
!
     IF (LwcAll <= zero)  THEN
        IF (ch2o(1) >= LwcLevel) THEN
           LWClass(:) = ch2o(1)
           met_new%ModAqua = 1
        ELSE
           LWClass(:) = zero
           met_new%ModAqua = 0
        END IF
     END IF
!
!======================================================
!===  Define Microphysics
!======================================================
!
!---  array allocation
     ALLOCATE(mphys(3))
     DO it=1,3
        ALLOCATE(mphys(it)%mpStat(nC,ntFrac))
        ALLOCATE(mphys(it)%mpAnf(nC,2))
        ALLOCATE(mphys(it)%mpEnd(nC,2))
        ALLOCATE(mphys(it)%LwcSum(nC))
     !
        ALLOCATE(mphys(it)%cradius(nC,ntFrac))
        ALLOCATE(mphys(it)%LWC(nC,ntFrac))
        ALLOCATE(mphys(it)%LWCFlux(nC,ntFrac,ntFrac))
        ALLOCATE(mphys(it)%ActFlux(nC,ntFrac,ntAqua))
     END DO
     mp_cur => mphys(1)

!---  set pointer
     idra = 1
     idre = ntFrac

     it = mp_index_old

     mp_new => mphys(it)
     mpAnf  => mp_new%mpAnf (:,:)
     mpEnd  => mp_new%mpEnd (:,:)
     mpStat => mp_new%mpStat(:,:)

     lwc     => mp_new%lwc(:,:)
     cradius => mp_new%cradius(:,:)

     mp_time_old = Time
     mp_time_new = Time
     mp_new%time = mp_time_new

     DO ic=1,nC
        LWC(ic,idra:idre)    = LWClass(1:ntFrac)   
        cradius(ic,1:ntFrac) = radcon (1:ntFrac)
     END DO

     mpAnf(:,:) = idra
     mpEnd(:,:) = idre

     mpStat(:,:) = 0
     mpStat(:,idra:idre) = 1

!--------------------------------------------------
!---  Swich Pointers
!--------------------------------------------------
!
     mp_index_old = mp_index_new
     mp_index_new = it
     mp_time_next = met_time_next

     mp_cur => mphys(1)
     mp_new => mphys(mp_index_new)
     mp_old => mphys(mp_index_old)

!--  initialization is ready
     mpFlag = 1         

!--  check liquid water switching
     FlagAqua = 0         ! the same status
!
!--------------------------------------------------
  END SUBROUTINE InitMPhys
