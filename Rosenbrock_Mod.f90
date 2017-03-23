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
  USE Chemsys_Mod
  USE Factorisation_Mod
  USE Rates_Mod
  USE mo_control
  USE mo_reac
  USE mo_ckinput
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
    INTEGER :: nsteps     = 0              ! # succ. steps
    INTEGER :: nfailed    = 0              ! # failed steps
    INTEGER :: nRateEvals = 0              ! # Rate evaluation
    INTEGER :: npds       = 0              ! # Jacobian evaluation
    INTEGER :: ndecomps   = 0              ! # LU factorisation
    INTEGER :: nsolves    = 0              ! # solved lin algebra
    REAL(RealKind) :: Ttimestep = 0.0d0  ! mean Time for one ROW step
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
  SUBROUTINE InitialStepSize(h,hmin,absh,Jac,Rate,t,Y,pow)
    !------------------------------------------------- 
    ! Input:
    !        - public variables
    !        - Tspan 
    !        - Y0  ( initial vector )
    !        - Jacobian matrix
    REAL(RealKind), INTENT(IN) :: t, pow
    REAL(RealKind), INTENT(IN) :: Y(nspc)
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
    ! DEBUG
    INTEGER :: i  

    zeros = ZERO    

    ! hmin is a small number such that t + hmin is clearlY different from t in
    ! the working precision, but with this definition, it is 0 if t = 0.
    hmin  = minStp

    !---- Compute an initial step size h using Yp=Y'(t) 
    CALL MatVecMult( f0 , BAT , Rate , Y_e )
    !print*, 'debug:: sum(bat),rate =', SUM(BAT%val),SUM(rate)
    !print*, 'debug:: sum(f0)=', SUM(f0)
    !stop
    wt    = MAX( ABS(Y) , ThresholdStepSizeControl(1:nspc) )
    rh    = ( 1.25D0 * MAXVAL( ABS(f0(:)/wt(:)) ) )/(RTolRow**pow)
    absh  = MIN( maxStp , Tspan(2)-Tspan(1) )
    IF ( absh * rh > ONE )  absh = ONE / rh
    !print*, 'Debug:: f0 = ', f0, SUM(Y_e), SUM(Rate)
    !do i=1,nspc
    !  print*, Rate(i)
    !end do 
    !stop

    !---- Compute Y''(t) and a better initial step size
    h     = absh
    !tdel  = ( t + MIN( sqrteps * MAX( ABS(t) , ABS(t+h) ) , absh ) ) - t
    tdel  = t + MIN( sqrteps * MAX( ABS(t) , ABS(t+h) ) , absh )

    CALL Rates( t+tdel , Y , Rate , DRatedT )
    Output%nRateEvals = Output%nRateEvals + 1

    CALL MatVecMult( f1 , BAT , Rate , Y_e )
 
    DfDt  = ( f1 - f0 ) / tdel
    CALL MatVecMult( Tmp , Jac , f0 , zeros )
    DfDt  = DfDt  + Tmp
  
    rh    = 1.25D0  * SQRT( rTWO * MAXVAL( ABS(DfDt/wt) ) ) / RTolRow**pow
   
    absh  = MIN( maxStp , Tspan(2)-Tspan(1) )
    IF ( absh * rh > ONE )  absh = ONE / rh
    absh  = MAX( absh , hmin )
    
  END SUBROUTINE InitialStepSize
  !
  !
  !=========================================================================
  !    Subroutine Rosenbrock-Method universal for classic and extended case
  !=========================================================================
  SUBROUTINE Rosenbrock(Y0,t,h,RCo,err,errind,YNew)
    !--------------------------------------------------------
    ! Input:
    !   - Y0............. actual concentrations Y 
    !   - t.............. time
    !   - h.............. step size
    !   - RCo............ Rosenbrock method
    !   - Temp........... actual Temperatur (optional for combustion)
    !
    REAL(RealKind), INTENT(IN) :: Y0(nDIM)
    REAL(RealKind), INTENT(IN) :: t, h
    TYPE(RosenbrockMethod_T)   :: RCo
    !--------------------------------------------------------
    ! Output:
    !   - Ynew........... new concentratinos 
    !   - err............ error calc with embedded formula.
    !   - TempNew........ new temperature (optional for combustion)
    !
    REAL(RealKind), INTENT(OUT) :: YNew(nDIM)
    REAL(RealKind), INTENT(OUT)   :: err
    INTEGER       , INTENT(OUT)   :: errind(1,1)
    !-------------------------------------------------------
    ! TemporarY variables:
    !
    REAL(RealKind), DIMENSION(nDIM)            :: Y, Yrh, Yhat, fRhs
    REAL(RealKind), DIMENSION(nspc)            :: U, UMat, dUdT, dCdt, d2UdT2 
    REAL(RealKind), DIMENSION(nDIMex)          :: bb
    !
    REAL(RealKind) :: k( nDIM , RCo%nStage )
    !
    REAL(RealKind) :: Tarr(8)
    REAL(RealKind) :: Rate(neq), rRate(neq)
    REAL(RealKind) :: DRatedT(neq)        
    REAL(RealKind) :: C(nspc)       ! molar heat capacities at constant pressure
    REAL(RealKind) :: H_e(nspc)       ! the standardstate molar enthalpY
    REAL(RealKind) :: S(nspc)       ! standard-state entropY at 298 K
    !
    REAL(RealKind) :: dHdT(nspc)    ! EnthaplY derivative in dT [J/mol/K^2]
    REAL(RealKind) :: dGdT(nspc)    ! Gibbs potential derivative in dT [J/mol/K^2]
    REAL(RealKind) :: dTdt
      
    REAL(RealKind) :: tt
    REAL(RealKind) :: cv ! mass average mixture specific heat at constant volume
    REAL(RealKind) :: dcvdT ! mass average mixture specific heat at constant volume
    REAL(RealKind) :: X
    !
    ! fuer verlgeich mit speedchem, andere spc reihenfolge
    !
    INTEGER :: iStg, jStg, i, rPtr          ! increments
    LOGICAL :: dprint=.false.
    INTEGER :: iprnt

    !print*, nDIM
    !print*, Y0
    !stop
    !
    iprnt = MinVAL((/ 5 , neq , nspc /))

    ! Initial settings
    k(:,:)    = ZERO
    fRhs(:)   = ZERO
    Rate(:)   = ZERO
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
   
    ! HIER UNBEDINGT RATE MIT Y0 ALS INPUT SONST SUCHST DU DICH WIEDER DUMM
    Y(1:nspc)   = MAX( ABS(Y0(1:nspc))   , eps ) * SIGN( ONE , Y0(1:nspc) )  ! concentrations =/= 0

    ! HIER MUSS DIE RATE MIT Y UND NICHT MIT Y0 BERECHNET WERDEN
    CALL Rates( t, Y, Rate, DRatedT )
    Rate(1:neq) = MAX( ABS(Rate(1:neq)) , eps ) * SIGN( ONE , Rate(1:neq) )    ! reaction rates =/= 0

    !do i=1,neq
    !  write(876,*) i,rate(i)
    !end do
    !stop
    
    ! UND HIER BLOß NICHT AUF DIE IDEE KOMMEN VON Y(:nspc) AUF Y0(:nspc) ZU ÄNDERN
    Yrh(1:nspc)  = Y(1:nspc) / h

    !
    !cv = ONE  ! for tropospheric mechanism

    IF ( combustion ) THEN
      !                             OUT:      IN:
      CALL UpdateTempArray        ( Tarr    , Y0(nDIM) )       
      
      ! Compute species internal energies in moles [J/mol/K]
      CALL InternalEnergy         ( U       , Tarr)             

      ! Specific heat [J/mol/K]
      CALL DiffSpcInternalEnergy  ( dUdT    , Tarr)              

      ! Derivatives of specific heats at constant volume in [J/mol/K2]
      CALL Diff2SpcInternalEnergy ( d2UdT2  , Tarr)
      
      ! Compute system pressure [Pa]
      CALL pressureHOT            ( Press   , Y0(:nspc) , Tarr(1)  )  
      
      ! Average mixture properties, constant volume specific heats [J/kg/K]
      CALL MassAveMixSpecHeat     ( cv      , Y0(:nspc) , dUdT )
      ! Div. Average mixture properties, constant volume specific heats [J/kg/K2]
      CALL MassAveMixSpecHeat     ( dcvdT   , Y0(:nspc) , d2UdT2 )

      ! Computing molar concentration rate of change imposing mass consv. [mol/cm3/s]
      CALL MatVecMult             ( dCdt    , BAT   , Rate , Y_e )

      dTdt     = - SUM( U * dCdt)/cv
      UMat     = - U * Yrh
      dRatedT  = dRatedT * RCo%ga
      X = cv/h + RCo%ga * ( dTdt * dcvdT  + SUM( dUdT * dCdt ) )
      !
      !
      IF (dprint) THEN
        print*, 'debug     Temparr=  ',Tarr
        print*, 'debug        time=  ',t
        print*, 'debug          cv=  ',cv
        print*, 'debug         SCP=  ',Press
        print*, 'debug     SUM(Y0)=  ',SUM(Y0(1:iprnt))
        print*, 'debug   SUM(Umol)=  ',SUM(U)
        print*, 'debug   SUM(dCdt)=  ',SUM(dCdt)
        print*, 'debug          X =  ',X
        print*, 'debug    U(b-a)r =  ',-SUM(U*dCdt)/cv
      END IF
      !
      !stop
      !
    END IF
   
    ! --- Update matrix procedure
    IF ( CLASSIC ) THEN
    
      ! classic case needs to calculate the Jacobian first
      TimeJacobianA = MPI_WTIME()
      CALL Miter_Classic( BAT , A , Rate , Y , h , RCo%ga , Miter )
      TimeJac       = TimeJac     + (MPI_WTIME()-TimeJacobianA)
      Output%npds   = Output%npds + 1

      IF ( OrderingStrategie == 8 ) THEN
        CALL SetLUvaluesCL( LU_Miter , Miter , LU_Perm , U , DRatedT )
      END IF

    ELSE !IF ( EXTENDED )
     
      rRate  = ONE / Rate
      IF ( OrderingStrategie == 8 ) THEN
        CALL SetLUvaluesEX  ( LU_Miter  , nspc , neq , rRate , Yrh, &
                            & DRatedT   , U    , X   , LUvalsFix    )
        !call writesparsematrix(LU_Miter,'LU_Miter_vor')
        !do i=1,size(LU_Miter%val)
        !  write(888,*) lu_Miter%val(i)
        !end do
        !stop
      ELSE
        CALL Miter_Extended ( Miter     , nspc , neq , rRate , Yrh )
        !call printsparsematrix(Miter,'Miter')
      END IF

    END IF
    !
    !   _____             _                _             __  __         _ _       
    !  |  ___|__ _   ___ | |_  ___   _ __ (_) ____ ___  |  \/  |  __ _ | |_  _ __ (_)__  __
    !  | |_  / _` | / __|| __|/ _ \ | '__|| ||_  // _ \ | |\/| | / _` || __|| '__|| |\ \/ /
    !  |  _|| (_| || (__ | |_| (_) || |   | | / /|  __/ | |  | || (_| || |_ | |   | | >  < 
    !  |_|   \__,_| \___| \__|\___/ |_|   |_|/___|\___| |_|  |_| \__,_| \__||_|   |_|/_/\_\
    !
    ! --- LU - Decomposition ---
    timerStart  = MPI_WTIME()
    IF (OrderingStrategie==8) THEN
      CALL SparseLU( LU_Miter )
      !call printsparse(LU_Miter,'LU_Miter_nach_LU')
        !call writesparsematrix(LU_Miter,'LU_Miter_nach')
    ELSE
      CALL MumpsLU( Miter%Val )
    END IF
    TimeFac     = TimeFac + (MPI_WTIME()-timerStart)


    IF (dprint) THEN
      WRITE(*,*) '----------------------------'
      WRITE(*,*) 'debug h, t, Temp  :: ', h , t, Y(nDIM)
      WRITE(*,*) '    rate(1:3),sum :: ', rate(1:iprnt), SUM(rate)
      WRITE(*,*) '   conc0(1:3),sum :: ', y0(1:iprnt),SUM(Y0)
      !IF(combustion) WRITE(*,*) '   DRatedT(1:3)   :: ', DRatedT(1:3), SUM(DRatedT)
      WRITE(*,*) '      sum(Miter)  :: ', SUM(Miter%val)
      do i=1,neq
        write(789,*), ' Miter%val = ', Miter%val(i)
      end do
      WRITE(*,*) '           Miter  :: ', Miter%val(1:iprnt)
      IF(OrderingStrategie==8) THEN
        WRITE(*,*) '           LUvor  :: ', LU_Miter%val(1:iprnt)
        WRITE(*,*) '      sum(LU)vor  :: ', SUM(LU_Miter%val)
      END IF
      WRITE(*,*)
    END IF

    !
    !****************************************************************************************
      !stop

    !****************************************************************************************
    !   ____    ___ __        __          _____  _                    ____   _               
    !  |  _ \  / _ \\ \      / /         |_   _|(_) _ __ ___    ___  / ___| | |_  ___  _ __  
    !  | |_) || | | |\ \ /\ / /   _____    | |  | || '_ ` _ \  / _ \ \___ \ | __|/ _ \| '_ \ 
    !  |  _ < | |_| | \ V  V /   |_____|   | |  | || | | | | ||  __/  ___) || |_|  __/| |_) |
    !  |_| \_\ \___/   \_/\_/              |_|  |_||_| |_| |_| \___| |____/  \__|\___|| .__/ 
    !                                                                                 |_|    
    !****************************************************************************************
    !
    LOOP_n_STAGES:  DO iStg = 1 , RCo%nStage

      IF ( iStg==1 ) THEN

        IF ( EXTENDED ) THEN
          bb( 1      : neq )          = mONE      ! = -1.0d0
          bb( neq+1  : nsr )          = Y_e       ! emission
          IF ( combustion ) bb(nDIMex)  = ZERO    ! =  0.0d0
          IF (dprint) print*, 'debug:: ', 'istage, vor solve   bb = ',iStg,bb(1:iprnt)
        END IF

      ELSE ! iStage > 1 ==> Update time and concentration

        tt  = t + RCo%Asum(iStg) * h
        Y   = Y0
        DO jStg = 1 , iStg
          Y = Y + RCo%a(iStg,jStg) * k(:,jStg)
        END DO
        IF (dprint) WRITE(*,*) '      conc(:)   :: ',istg, Y(1:iprnt)
        
        ! Update Rates at  (t + SumA*h) , and  (Y + A*)k
        CALL Rates( tt , Y , Rate , DRatedT )

        IF (dprint) print*, 'DEBUG:: sum(rate)=',SUM(Rate)
      END IF
      
      !--- Calculate the right hand side of the linear System
      IF ( CLASSIC ) THEN

        CALL MatVecMult( fRhs , BAT , Rate , Y_e )           
        fRhs = h * fRhs

        DO jStg = 1 , iStg-1
          fRhs  = fRhs + RCo%C(iStg,jStg) * k(:,jStg)
        END DO

        IF (dprint) print*, 'debug:: ', 'istage,  frhs = ',iStg,frhs(1:iprnt)

      ELSE !IF ( EXTENDED ) 

        IF ( iStg/=1 ) THEN

          fRhs(:)     = ZERO
          bb(nDIMex)  = ZERO

          DO jStg = 1 , iStg-1
            fRhs(:nspc)   = fRhs(:nspc) + RCo%C(iStg,jStg)*k(:nspc,jStg)
            IF (combustion) &                                            
              bb(nDIMex)  = bb(nDIMex)  + RCo%C(iStg,jStg)*( cv*k(nDIM,jStg) &
                                                & - SUM( U * k(1:nspc,jStg)) )
          END DO

          bb( 1      : neq ) = -rRate * Rate
          bb( neq+1  : nsr ) = Y_e + fRhs(:nspc)/h
          IF (combustion) bb(nDIMex) = - TWO*SUM(U * Y_e) + bb(nDIMex)/h

          IF (dprint) print*, 'debug:: ', 'istage, vor solve   bb = ',iStg,bb(1:iprnt)
        END IF

      END IF
      !
      !print*, 'debug:: ', 'istage=',iStg,SUM(ABS(fRhs(:)))
      !
      ! ---  solve LGS  ---
      !
      timerStart  = MPI_WTIME()
      SOLVE_LINEAR_SYSTEM: IF ( OrderingStrategie==8 ) THEN  

        IF ( CLASSIC ) THEN

          CALL SolveSparse( LU_Miter , fRhs )
          k( : , iStg ) = fRhs(:)
          IF (dprint) print*, 'debug:: ', 'istage, nach solve  k(:iprint,iStage) = ',iStg,frhs(1:iprnt)

        ELSE !IF ( EXTENDED )

          CALL SolveSparse( LU_Miter , bb)
          k( 1:nspc , iStg ) = Y0(:nspc) * bb(neq+1:nsr)
          IF ( combustion ) k(nDIM,iStg) = bb(nDIM)
          
          IF (dprint) print*, 'debug:: ', 'istage, nach solve  bb = ',iStg,bb(1:iprnt)
          
        END IF

      ELSE

        IF ( CLASSIC ) THEN

          CALL MumpsSolve( fRhs )
          k( 1:nDIM , iStg ) = Mumps_Par%RHS(:)    
          IF (dprint) print*, 'debug:: ', 'istage, nach solve  k(:iprint,iStage) = ',iStg,Mumps_Par%RHS(1:iprnt)

        ELSE !IF ( EXTENDED )

          CALL MumpsSolve( bb )
          k( 1:nspc , iStg ) = Y0(:nspc) * Mumps_Par%RHS(neq+1:nsr)
          IF ( combustion ) k(nDIM,iStg) = Mumps_Par%RHS(nDIMex)

          IF (dprint) print*, 'debug:: ', 'istage, nach solve  bb = ',iStg,Mumps_Par%RHS(1:iprnt)

        END IF

      END IF SOLVE_LINEAR_SYSTEM
      TimeSolve   = TimeSolve + (MPI_WTIME()-timerStart)

    END DO  LOOP_n_STAGES

    
    !--- Update Concentrations (Temperatur)

    YNew(:)   = Y0(:)
    YHat(:)   = Y0(:)

    DO jStg = 1 , RCo%nStage
      Ynew(:) = Ynew(:) + RCo%m(jStg) * k(:,jStg)! new Y vector
      Yhat(:) = Yhat(:) + RCo%me(jStg) * k(:,jStg)! embedded formula for err calc ord-1
    END DO
    !
    IF (dprint) THEN
      print*, ''
      !print*, 'debug:: k(ndim)=',k(nDIM,:)
      print*, 'debug:: coefs  =',RCo%m(:)
      do istg=1,nspc
        print*,'debug::   y0  yn    ', Y0(istg),Ynew(iStg)
      end do
    END IF
    !***********************************************************************************************
    !   _____                           _____       _    _                    __               
    !  | ____| _ __  _ __  ___   _ __  | ____| ___ | |_ (_) _ __ ___    __ _ | |_  (_)  ___   _ __  
    !  |  _|  | '__|| '__|/ _ \ | '__| |  _|  / __|| __|| || '_ ` _ \  / _` ||  __|| | / _ \ | '_ \ 
    !  | |___ | |   | |  | (_) || |    | |___ \__ \| |_ | || | | | | || (_| || |_  | || (_) || | | |
    !  |_____||_|   |_|   \___/ |_|    |_____||___/ \__||_||_| |_| |_| \__,_| \__| |_| \___/ |_| |_|
    !                                                                                              
    !***********************************************************************************************
    CALL ERROR( err , errind , Ynew , Yhat , ATolAll , RTolROW , t )
    !
    CALL MPI_BARRIER(MPI_COMM_WORLD,MPIErr)
    IF (dprint) print*,'debug::   error      ', err
    !print*,'debug::   error      ', err, SUM(Ynew)
    IF (dprint) stop 'rosenbrockmod'
  END SUBROUTINE Rosenbrock


  SUBROUTINE ERROR(err,En_Index,ynew,yhat,ATol,RTol,t)
    !
    REAL(RealKind) :: err
    REAL(RealKind), DIMENSION(:) :: ynew, yhat, ATol
    REAL(RealKind) :: RTol, t
    !
    REAL(RealKind) :: scalTol(nDIM), En_Values(nDIM)
    INTEGER :: En_Index(1,1)
    !
    scalTol       = ATol + MAX( ABS(yhat) , ABS(ynew) ) * RTol ! scaling strategie
    En_Values     = ABS( ynew - yhat ) / scalTol               ! local error est.
    En_Index(1,1) = MAXLOC( En_Values , 1 )                    ! max error component
    !
    IF ( Error_Est == 2 ) THEN
      err = SUM( En_Values*En_Values ) / nspc   ! euclikd norm
    ELSE
      err = MAXVAL( ABS(En_Values) )     ! maximum norm
    END IF
    !
  END SUBROUTINE ERROR

END MODULE Rosenbrock_Mod
