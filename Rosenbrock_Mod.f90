!=========================================================================!
!                                                                         !         
!                                                                         !
!          Rosenbrock Modul for integrating one time step                 ! 
!                                                                         ! 
!                                                                         !
!=========================================================================!

MODULE Rosenbrock_Mod
  !
  USE Kind_Mod
  USE Sparse_Mod
  USE ErrorROW_Mod
  USE Chemsys_Mod
  USE Factorisation_Mod
  USE Rates_Mod
  USE mo_control
  USE mo_reac
  IMPLICIT NONE
  !
  !
  ! Rosenbrock-Parameter
  TYPE RosenbrockMethod_T
    INTEGER :: Order                                ! Classical approximation order of the method
    INTEGER :: nStage                               ! Number of stages
    INTEGER :: pinterp                              ! Interpolation order
    REAL(RealKind) :: ga                            ! Diagonalentry gamma
    REAL(RealKind) :: pow                           ! needed for sitepsize control pow=1/nstage
    REAL(RealKind), ALLOCATABLE :: Asum(:)          ! Row sum of A
    REAL(RealKind), ALLOCATABLE :: Alpha(:,:)       ! Propagation table, strictly lower triangular
    REAL(RealKind), ALLOCATABLE :: a(:,:)           ! Propagation table, strictly lower triangular (converted Alpha)
    REAL(RealKind), ALLOCATABLE :: Gamma(:,:)       ! Stage table, lower triangular with nonzero diagonal
    REAL(RealKind), ALLOCATABLE :: iGamma(:,:)      ! inverse Stage table
    REAL(RealKind), ALLOCATABLE :: C(:,:)           ! Stage table, lower triangular with nonzero diagonal (converted Gamma)
    REAL(RealKind), ALLOCATABLE :: B(:)             ! Step completion table
    REAL(RealKind), ALLOCATABLE :: m(:)             ! Step completion table(converted B)
    REAL(RealKind), ALLOCATABLE :: Be(:)            ! Step completion table for embedded method of order one less
    REAL(RealKind), ALLOCATABLE :: me(:)            ! Step completion table for embedded method of order one less (converted Be)
    REAL(RealKind), ALLOCATABLE :: binterpt(:,:)    ! Dense output formula
  END TYPE RosenbrockMethod_T
  
  
  TYPE IntArgs
    INTEGER :: nep                                  ! length y vector
    REAL(RealKind) :: Tend                          ! end of integration intervall
    INTEGER :: Tdir                                 ! direction (1 if t is ascending, -1 if t is descending)
    REAL(RealKind), ALLOCATABLE :: f0(:)            ! first rhs eval
    REAL(RealKind), ALLOCATABLE :: threshold(:)     ! ATol/RTol
    REAL(RealKind) :: hmax                          ! max step size
    REAL(RealKind) :: hTspan
    REAL(RealKind), ALLOCATABLE :: ATol(:)          ! absolute tolerances AtolGas
    REAL(RealKind) :: RTol                          ! relative tolerance RtolROW
    REAL(RealKind) :: RTolpow
  END TYPE IntArgs

  TYPE Out
    REAL(RealKind), ALLOCATABLE :: y(:)    ! y-vector at Tend
    INTEGER :: nsteps=0                    ! # succ. steps
    INTEGER :: nfailed=0                   ! # failed steps
    INTEGER :: nRateEvals=0                ! # Rate evaluation
    INTEGER :: npds=0                      ! # Jacobian evaluation
    INTEGER :: ndecomps=0                  ! # LU factorisation
    INTEGER :: nsolves=0                   ! # solved lin algebra
    REAL(RealKind) :: Ttimestep=0.0d0  ! mean Time for one ROW step
  END TYPE Out
  TYPE(Out) :: Output

  TYPE(CSR_Matrix_T) :: Miter                       ! 
  TYPE(CSR_Matrix_T) :: LU_Miter

  TYPE(SpRowIndColInd_T) :: MiterFact
  !TYPE(SpRowColD_T) :: MiterMarko

  TYPE(IntArgs), PUBLIC :: Args                     ! Initial arguments  
  
  REAL(RealKind), PRIVATE :: timerEnd,timerStart
  REAL(RealKind), ALLOCATABLE :: LUvalsFix(:)

  INTEGER, ALLOCATABLE :: LU_Perm(:)

  CONTAINS
  !
  !
  !=======================================================
  !       chose one of the methods in ~/METHODS/
  !=======================================================
  SUBROUTINE SetRosenbrockMethod(RCo,method)
    !-------------------------------------------------------------
    ! Input: method ... string with Rosenbrock method path
    TYPE(RosenbrockMethod_T), INTENT(OUT) :: RCo
    CHARACTER(*) :: method
    !-------------------------------------------------------------
    ! Output:
    !   - RosenbrockMethod_T ... coefficients of the chosen one 
    !-------------------------------------------------------------
    ! Temporary variables: 
    REAL(RealKind), ALLOCATABLE :: ID(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:)
    INTEGER :: INFO
    !
    INTEGER :: i
    !
    method=method(INDEX(method,'/')+1:INDEX(method,'.')-1)
    ! dynamisches inlcude möglich???
    ! oder zeile für zeile READ ?
    SELECT CASE (method)
      CASE ('Ros2AMF')
        INCLUDE 'METHODS/Ros2AMF.fort'
      CASE ('Ros3w')
        INCLUDE 'METHODS/Ros3w.fort'
      CASE ('Ros3Pw')
        INCLUDE 'METHODS/Ros3Pw.fort'
      CASE ('Ros34PW1a')
        INCLUDE 'METHODS/Ros34PW1a.fort'
      CASE ('Ros34PW2')
        INCLUDE 'METHODS/Ros34PW2.fort'
      CASE ('Ros34PW3')
        INCLUDE 'METHODS/Ros34PW3.fort'
      CASE ('TSRosW2P')
        INCLUDE 'METHODS/TSRosW2P.fort'
      CASE ('TSRosW2M')
        INCLUDE 'METHODS/TSRosW2M.fort'
      CASE ('TSRosWRA3PW')
        INCLUDE 'METHODS/TSRosWRA3PW.fort'
      CASE ('TSRosWRA34PW2')
        INCLUDE 'METHODS/TSRosWRA34PW2.fort'
      CASE ('TSRosWRodas3')
        INCLUDE 'METHODS/TSRosWRodas3.fort'
      CASE ('TSRosWSandu3')
        INCLUDE 'METHODS/TSRosWSandu3.fort'
      CASE DEFAULT
        IF (MPI_ID==0) WRITE(*,*) '    Unknown Method:  ',method
        IF (MPI_ID==0) WRITE(*,*) '    Use TSRosWRodas3 instead.'
        INCLUDE 'METHODS/TSRosWRodas3.fort'
    END SELECT
    !INCLUDE RosenbrockMethod
    !
    ! converting the butcher tableau 
    ! automatic transformation to avoid mat*vec in ROW methode
    RCo%pow=1.0d0/(RCo%Order+1.0d0)
    !RCo%pow=1.0d0/RCo%nStage
    ALLOCATE(RCo%iGamma(RCo%nStage,RCo%nStage))
    RCo%iGamma=0.0d0
    ALLOCATE(ID(RCo%nStage,RCo%nStage))
    ID=0.0d0
    DO i=1,RCo%nStage
      RCo%iGamma(i,i)=1.0d0
      ID(i,i)=1.0d0
    END DO
    !  
    ! calculate the inverse matrix of gamma
    !
    ! CAUTION RCo%Gamma (IN) =/= RCo%Gamma (OUT) !!!
    !
    ALLOCATE(IPIV(RCo%nStage))
    CALL dgesv(  RCo%nStage,     &        ! # linear eqations
    &            RCo%nStage,     &        ! # RHS (coloums)
    &            RCo%Gamma,      &        ! Matrix A of A*A^(-1)=ID
    &            RCo%nStage,     &        ! leading dimension of A (nStage)
    &            IPIV,           &        ! pivot indices of dimension nStage
    &            RCo%iGamma,     &        ! Matrix ID of A*A^(-1)=ID
    &            RCo%nStage,     &        ! leading dimension of RHS
    &            INFO)                    ! INFO (integer) if INFO=0 succsessfull	
    !
    IF (INFO/=0.AND.MPI_ID==0) WRITE(*,*) 'Error while calc row-method parameter'
    !       
    ALLOCATE(RCo%a(RCo%nStage,RCo%nStage))
    RCo%a=0.0d0
    RCo%a=RCo%ga*MATMUL(RCo%Alpha, RCo%iGamma)
    !  
    ALLOCATE(RCo%C(RCo%nStage,RCo%nStage))
    RCo%C=0.0d0
    RCo%C=ID-RCo%ga*RCo%iGamma
    FORALL (i=1:RCo%nStage) RCo%C(i,i)=0.0d0
    !  
    ALLOCATE(RCo%m(RCo%nStage))
    RCo%B=MATMUL(RCo%B, RCo%iGamma)
    RCo%m=RCo%ga*RCo%B(:)
    !
    ALLOCATE(RCo%me(RCo%nStage))
    RCo%Be=MATMUL(RCo%Be,RCo%iGamma)
    RCo%me=RCo%ga*RCo%Be(:)
    !
    DEALLOCATE(ID)
    DEALLOCATE(IPIV)
  END SUBROUTINE SetRosenbrockMethod
  !
  !
  !=================================================================
  !         Check some values and set few parameters
  !=================================================================
  SUBROUTINE CheckInputParameters(Tspan,aTol,rTolROW)
    !--------------------------------------------------------------------
    ! Input:
    !   - Tspan
    !   - aTol, rTolROW
    REAL(RealKind) :: Tspan(2)               
    REAL(RealKind) :: aTol(2)
    REAL(RealKind) :: rTolROW
    !
    ALLOCATE(ThresholdStepSizeControl(nDIM))
    ThresholdStepSizeControl(:ntGas)=AtolGas/RTolROW
    ThresholdStepSizeControl(ntGas+1:)=AtolAqua/RTolROW
    ALLOCATE(ATolAll(nDIM))
    ATolAll(:ntGas)=ATolGas
    ATolAll(ntGas+1:)=ATolAqua
    IF ( combustion ) THEN
      ThresholdStepSizeControl(nDIM)=ATolTemp/RTolROW
      ATolAll(nDIM)=ATolTemp
    END IF
    ! Test that Tspan is internally consistent.
    IF ( Tspan(1)>=Tspan(2) ) THEN
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) 'start time >= end time'
      CALL FinishMPI()
      STOP '  STOP  '
    END IF
    !
    ! Test that Rosenbrock tolerance > 0
    IF ( rTolROW<=ZERO ) THEN
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) 'RTol must be positiv scalar'
      CALL FinishMPI()
      STOP '  STOP  '
    END IF
    !
    ! Test that absolute tolerance for gas and aqua species is > 0
    IF ( .NOT.(ALL(aTol>=ZERO)) ) THEN
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) 'ATols must be positive'
      CALL FinishMPI()
      STOP '  STOP  '
    END IF
    !
    ! Test if maximum stepsize is not to small/big
    IF ( (maxStp<=ZERO).OR.(maxStp>Tspan(2)-Tspan(1)) ) THEN
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) 'Maximum stepsize = ',maxStp, ' to low or to high'
      CALL FinishMPI()
      STOP 'STOP'
    END IF
  END SUBROUTINE CheckInputParameters
  !
  !
  !==========================================================
  !   Calculates an initial stepsize based on 2. deriv.
  !==========================================================
  SUBROUTINE InitialStepSize(h,hmin,absh,Jac,Rate,t,y,pow)
    !------------------------------------------------- 
    ! Input:
    !        - public variables
    !        - Tspan 
    !        - y0  ( initial vector )
    !        - Jacobian matrix
    REAL(RealKind), INTENT(IN) :: t, pow
    REAL(RealKind), INTENT(IN) :: y(nspc)
    TYPE(CSR_Matrix_T), INTENT(IN) :: Jac
    REAL(RealKind) :: Rate(neq)
    REAL(RealKind) :: DRatedT(neq)     ! part. derv. rate over temperatur vector
    !-------------------------------------------------
    ! Output:
    !        - initial step size
    REAL(RealKind)  , INTENT(OUT) :: absh
    REAL(RealKind) :: h, hmin 
    REAL(RealKind) :: f0(nspc)
    !-------------------------------------------------
    !
    ! Temp vars:
    REAL(RealKind) :: tdel, rh
    REAL(RealKind), DIMENSION(nspc) ::  wt, DfDt, Tmp, f1, zeros
    REAL(RealKind) :: sqrteps=SQRT(eps)
    !
    zeros=ZERO    
    !
    ! hmin is a small number such that t + hmin is clearly different from t in
    ! the working precision, but with this definition, it is 0 if t = 0.
    hmin=minStp
    !
    !print*,'debug:: in ', SUM(ABS(y))
    !---- Compute an initial step size h using yp=y'(t) 
    CALL MatVecMult(BAT,Rate,y_e,f0)
    wt=MAX(ABS(y),ThresholdStepSizeControl(1:nspc))
    rh=(1.25D0*MAXVAL(ABS(f0(:)/wt(:))))/(RTolRow**pow)
    !
    !print*,'debug:: pw,nspcm', RTolRow**pow, nspc
    !print*,
    !print*, 'debug:: SUM(wt),rh,sum(f0)', sum(ABS(wt)) , rh, SUM(ABS(f0))
    absh=MIN(maxStp,Tspan(2)-Tspan(1))
    IF (absh*rh>1.0D0) THEN
      absh=1.0D0/rh
    END IF
    !
    !---- Compute y''(t) and a better initial step size
    h=absh
    tdel=(t+MIN(sqrteps*MAX(ABS(t),ABS(t+h)),absh))-t
    !print*, 'debug:: tdel=     ', tdel, t+tdel
    !
    CALL Rates((t+tdel),y,Rate,DRatedT)
    Output%nRateEvals=Output%nRateEvals+1
    !
    !print*, 'debug:: sumzzz(bat,rate,yem,f1)=     ', SUM(BAT%Val),SUM(Rate),SUM(y_e)
    CALL MatVecMult(BAT,Rate,y_e,f1)
    !print*, 'debug:: sum(f1)=     ', SUM(f1), SIZE(BAT%VAL),SIZE(rate),SIZE(y_e)
    !
    !stop 'stop'
    DfDt=(f1-f0)/tdel
    CALL MatVecMult(Jac,f0,zeros,Tmp)
    DfDt=DfDt+Tmp
    !
    rh=1.25D0*SQRT(0.5d0*MAXVAL(ABS(DfDt/wt)))/(RTolRow**pow)
    !
    absh=MIN(maxStp,Tspan(2)-Tspan(1))
    IF (absh*rh>1.0d0) THEN
      absh=1.0D0/rh
    END IF
    absh=MAX(absh,hmin)
    !print*, 'debug:: h, absh', h, absh
    !stop
  END SUBROUTINE InitialStepSize
  !
  !
  !=========================================================================
  !    Subroutine Rosenbrock-Method universal for classic and extended case
  !=========================================================================
  SUBROUTINE Rosenbrock(y0,t,h,RCo,err,errind,yNew)
    !--------------------------------------------------------
    ! Input:
    !   - y0............. actual concentrations y 
    !   - t.............. time
    !   - h.............. step size
    !   - RCo............ Rosenbrock method
    !   - Temp........... actual Temperatur (optional for combustion)
    !
    REAL(RealKind), INTENT(IN) :: y0(:)
    REAL(RealKind), INTENT(IN) :: t, h
    TYPE(RosenbrockMethod_T)   :: RCo
    !--------------------------------------------------------
    ! Output:
    !   - ynew........... new concentratinos 
    !   - err............ error calc with embedded formula.
    !   - TempNew........ new temperature (optional for combustion)
    !
    REAL(RealKind), INTENT(INOUT) :: yNew(:)
    REAL(RealKind), INTENT(OUT) :: err
    INTEGER       , INTENT(OUT) :: errind(1,1)
    !-------------------------------------------------------
    ! Temporary variables:
    !
    REAL(RealKind), DIMENSION(nDIM,RCo%nStage) :: k     
    REAL(RealKind), DIMENSION(nDIM)    :: y, yhat, fRhs
    REAL(RealKind), DIMENSION(nspc)    :: hy
    REAL(RealKind), DIMENSION(nspc)    :: Umol, DUmoldT, DcDt
    REAL(RealKind), DIMENSION(nDIMex)    :: bb
    !
    REAL(RealKind) :: Tarr(7)
    REAL(RealKind) :: Rate(neq)
    REAL(RealKind) :: DRatedT(neq)        
    REAL(RealKind) :: invRate(neq)
    REAL(RealKind) :: tt
    REAL(RealKind) :: c_v ! mass average mixture specific heat at constant volume
    REAL(RealKind) :: X
    !
    !
    INTEGER :: iStage, jStage, i, rPtr          ! increments
    !
    ! Initial settings
    k(:,:)=ZERO
    fRhs(:)=ZERO
    Rate(:)=ZERO
    !
    y(:nspc)=MAX(ABS(y0(:nspc)),eps)*SIGN(ONE,y0(:nspc))                   ! concentrations =/= 0
    IF ( combustion ) y(nDIM)=y0(nDIM)
    yNew(:)=y0(:)
    yHat(:)=y0(:)
    !
    !IF (PRESENT(Temp)) yNew(Temp_ind)=y0(Temp_ind)
    !
    !
    !********************************************************************************
    !    _   _             _         _          __  __         _          _       
    !   | | | | _ __    __| |  __ _ | |_  ___  |  \/  |  __ _ | |_  _ __ (_)__ __
    !   | | | || '_ \  / _` | / _` || __|/ _ \ | |\/| | / _` || __|| '__|| |\ \/ /
    !   | |_| || |_) || (_| || (_| || |_|  __/ | |  | || (_| || |_ | |   | | >  < 
    !    \___/ | .__/  \__,_| \__,_| \__|\___| |_|  |_| \__,_| \__||_|   |_|/_/\_\
    !          |_|                                                                
    !
    !********************************************************************************
    !
    ! --- Nessesary for combustion systems

    
    CALL Rates(t,y,Rate,DRatedT)
    Rate(:)=MAX(ABS(Rate(:)),eps)*SIGN(ONE,Rate(:))         ! reaction rates =/= 0
    hy=y(:nspc)/h

    !print*, 'debug:: ', h , t ,Rate(1:3)
    !print*, 'debug::2',y(1:3)
    !
    IF ( combustion ) THEN
      CALL UpdateTempArray(Tarr,y0(nDIM))
      CALL SpcInternalEnergy(Umol,Tarr)
      CALL DiffSpcInternalEnergy(DUmoldT,Tarr)
      CALL MassAveMixSpecHeat(c_v,Tarr,y0(:nspc),DUmoldT)
      CALL DiffConcDt(BAT,Rate,DcDt)
      DUmoldT=DUmoldT*RCo%ga
      c_v=c_v/h
      Umol(:)=-Umol(:)*hy(:)
      X=c_v+SUM(DUmoldT(:)*DcDt(:))
      print*, 'debug     Temparr  :: ',Tarr
      print*, 'debug     Umol     :: ',Umol
      print*, 'debug     DUmoldT  :: ',DUmoldT
      print*, 'debug     c_v      :: ',c_v
      print*, 'debug     DcDt     :: ',DcDt
      print*, 'debug     X        :: ',X
      stop
    END IF
    ! 
    ! --- Update matrix procedure
    IF ( solveLA=='cl') THEN
      !
      ! classic case needs to calculate the Jacobian first
      TimeJacobianA=MPI_WTIME()
      CALL Miter_Classic(BAT,A,Rate,y,h,RCo%ga,Miter)
      TimeJac=TimeJac+(MPI_WTIME()-TimeJacobianA)
      Output%npds=Output%npds+1
      !
      IF (OrderingStrategie==8) THEN
        CALL SetLUvaluesCL(LU_Miter,Miter,LU_Perm)
      END IF
    ELSE !IF ( solveLA=='ex') THEN
      !
      invRate(:)=ONE/Rate(:)
      IF (OrderingStrategie==8) THEN
        CALL SetLUvaluesEX(LU_Miter,nspc,neq,invRate,hy,DRatedT,Umol,X,LUvalsFix)
      ELSE
        CALL Miter_Extended(Miter,nspc,neq,invRate,hy)
      END IF
    END IF
    !
    WRITE(*,*) '----------------------------'
    WRITE(*,*) 'debug h, t, Temp  :: ', h , t, y(nDIM)
    WRITE(*,*) '      rate(1:3)   :: ', rate(1:3), SUM(rate)
    IF(combustion) WRITE(*,*) '   DRatedT(1:3)   :: ', DRatedT(1:3), SUM(DRatedT)
    WRITE(*,*) '      conc(1:3)   :: ', y(1:3), SUM(y)
    WRITE(*,*) '      sum(Miter)  :: ', SUM(Miter%val)
    WRITE(*,*) '      sum(LU)vor  :: ', SUM(LU_Miter%val)

    !
    !****************************************************************************************
    !   _____             _                _             __  __         _ _       
    !  |  ___|__ _   ___ | |_  ___   _ __ (_) ____ ___  |  \/  |  __ _ | |_  _ __ (_)__  __
    !  | |_  / _` | / __|| __|/ _ \ | '__|| ||_  // _ \ | |\/| | / _` || __|| '__|| |\ \/ /
    !  |  _|| (_| || (__ | |_| (_) || |   | | / /|  __/ | |  | || (_| || |_ | |   | | >  < 
    !  |_|   \__,_| \___| \__|\___/ |_|   |_|/___|\___| |_|  |_| \__,_| \__||_|   |_|/_/\_\
    !
    !****************************************************************************************
    !
    ! --- LU - Decomposition ---
    IF (OrderingStrategie==8) THEN
      !
      timerStart=MPI_WTIME()
      CALL SparseLU(LU_Miter)
      WRITE(*,*) '      sum(LU)nach :: ', SUM(LU_Miter%val)
      timerEnd=MPI_WTIME()
      TimeFac=TimeFac+(timerEnd-timerStart)
    ELSE
      CALL FactorizeCoefMat(Miter%Val)
    END IF
    !
    !stop
    !call printsparse(LU_miter,'*')
    !WRITE(*,*) '      sum(LU)2    :: ', SUM(LU_Miter%val)
    !WRITE(*,*) 
    !
    !****************************************************************************************
    !   ____    ___ __        __          _____  _                    ____   _               
    !  |  _ \  / _ \\ \      / /         |_   _|(_) _ __ ___    ___  / ___| | |_  ___  _ __  
    !  | |_) || | | |\ \ /\ / /   _____    | |  | || '_ ` _ \  / _ \ \___ \ | __|/ _ \| '_ \ 
    !  |  _ < | |_| | \ V  V /   |_____|   | |  | || | | | | ||  __/  ___) || |_|  __/| |_) |
    !  |_| \_\ \___/   \_/\_/              |_|  |_||_| |_| |_| \___| |____/  \__|\___|| .__/ 
    !                                                                                |_|    
    !****************************************************************************************
    !
    DO iStage=1,RCo%nStage
      IF ( iStage==1 ) THEN
        IF ( solveLA=='ex' ) THEN
          bb(:neq)=mONE
          bb(neq+1:NSactNR)=y_e(:)
        END IF
        IF (combustion) bb(nDIMex)=ZERO
      ELSE
        tt=t+RCo%Asum(iStage)*h
        y=y0
        DO jStage=1,iStage
          y=y+RCo%a(iStage,jStage)*k(:,jStage)
        END DO
        !
        ! Update Rates at  (t + SumA*h) , and  (y + A*)k
        CALL Rates(tt,y,Rate,DRatedT)
      END IF
      !print*, 'debug :: stage, rates=',iStage,Rate(1:5)
      !
      !
      IF ( solveLA=='cl') THEN
        CALL MatVecMult(BAT,Rate,y_e,fRhs)           
        fRhs=h*fRhs
        DO jStage=1,iStage-1
          fRhs=fRhs+RCo%C(iStage,jStage)*k(:,jStage)
        END DO
      ELSE 
        IF (iStage/=1) THEN
          bb(:neq)=-invRate(:)*Rate(:)
          fRhs=ZERO
          bb(nDIMex)=ZERO
          DO jStage=1,iStage-1
            fRhs(1:nspc)=fRhs(1:nspc)+RCo%C(iStage,jStage)*k(1:nspc,jStage)
            IF (combustion) THEN
              bb(nDIMex) = bb(nDIMex)  +  RCo%C(iStage,jStage)                  &
              &          * ( c_v*k(nDIM,jStage) - SUM(Umol(:)*k(1:nspc,jStage)) )
            END IF
          END DO
          bb(neq+1:NSactNR)=fRhs(1:nspc)/h+y_e(:)
          IF (combustion) bb(nDIMex)=bb(nDIMex)/h
        END IF
      END IF
      !
      print*, 'debug:: ', 'istage=',iStage,SUM(ABS(fRhs(:)))
      !
      ! ---  Solve linear System  ---
      !
      IF (OrderingStrategie==8) THEN  
        timerStart=MPI_WTIME()
        IF ( solveLA=='cl') THEN
          CALL SolveSparse(LU_Miter,fRhs)
          k(:,iStage)=fRhs(:)
        ELSE
          CALL SolveSparse(LU_Miter,bb)
          k(:,iStage)=y0(:)*bb(neq+1:)
        END IF
        TimeSolve=TimeSolve+(MPI_WTIME()-timerStart)
        !print*, 'debug:: ', 'istage=',iStage,SUM(ABS(fRhs(:)))
      ELSE
        timerStart=MPI_WTIME()
        IF ( solveLA=='cl') THEN
          CALL SolveLinAlg(fRhs)
          k(:,iStage)=Mumps_Par%RHS(:)    
        ELSE
          CALL SolveLinAlg(bb)
          k(:nspc,iStage)=y0(:nspc)*Mumps_Par%RHS(neq+1:NSactNR)
          IF (combustion) k(nDIM,iStage)=Mumps_Par%RHS(nDIM)
          print*, 'debug:: ', 'istage=',iStage,k(:,iStage)
        END IF
        TimeSolve=TimeSolve+(MPI_WTIME()-timerStart)
      END IF    
    END DO
      !call dropout
    !  
    DO jStage=1,RCo%nStage
      ynew(:)=ynew(:)+RCo%m(jStage)*k(:,jStage)! new y vector
      yhat(:)=yhat(:)+RCo%me(jStage)*k(:,jStage)! embedded formula for err calc ord-1
    END DO
    !
    !***********************************************************************************************
    !   _____                           _____       _    _                    __               
    !  | ____| _ __  _ __  ___   _ __  | ____| ___ | |_ (_) _ __ ___    __ _ | |_  (_)  ___   _ __  
    !  |  _|  | '__|| '__|/ _ \ | '__| |  _|  / __|| __|| || '_ ` _ \  / _` ||  __|| | / _ \ | '_ \ 
    !  | |___ | |   | |  | (_) || |    | |___ \__ \| |_ | || | | | | || (_| || |_  | || (_) || | | |
    !  |_____||_|   |_|   \___/ |_|    |_____||___/ \__||_||_| |_| |_| \__,_| \__| |_| \___/ |_| |_|
    !                                                                                              
    !***********************************************************************************************
    CALL ERROR(err,errind,ynew,yhat,y0,ATolAll,RTolROW,t)
    !
  END SUBROUTINE Rosenbrock
END MODULE Rosenbrock_Mod
