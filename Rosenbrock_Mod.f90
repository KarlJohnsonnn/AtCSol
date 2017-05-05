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

  TYPE(IntArgs), PUBLIC :: Args                     ! Initial arguments  
  
  REAL(RealKind), PRIVATE :: timerStart
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
    !method=method(INDEX(method,'/')+1:INDEX(method,'.')-1)
    ! dynamisches inlcude möglich???
    ! oder zeile für zeile READ ?

    SELECT CASE (TRIM(method(9:)))
      CASE ('Ros2AMF.fort')
        INCLUDE 'METHODS/Ros2AMF.fort'
      CASE ('Ros3w.fort')
        INCLUDE 'METHODS/Ros3w.fort'
      CASE ('Ros3Pw.fort')
        INCLUDE 'METHODS/Ros3Pw.fort'
      CASE ('Ros34PW1a.fort')
        INCLUDE 'METHODS/Ros34PW1a.fort'
      CASE ('Ros34PW2.fort')
        INCLUDE 'METHODS/Ros34PW2.fort'
      CASE ('Ros34PW3.fort')
        INCLUDE 'METHODS/Ros34PW3.fort'
      CASE ('TSRosW2P.fort')
        INCLUDE 'METHODS/TSRosW2P.fort'
      CASE ('TSRosW2M.fort')
        INCLUDE 'METHODS/TSRosW2M.fort'
      CASE ('TSRosWRA3PW.fort')
        INCLUDE 'METHODS/TSRosWRA3PW.fort'
      CASE ('TSRosWRA34PW2.fort')
        INCLUDE 'METHODS/TSRosWRA34PW2.fort'
      CASE ('TSRosWRodas3.fort')
        INCLUDE 'METHODS/TSRosWRodas3.fort'
      CASE ('TSRosWSandu3.fort')
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
    IF ( TempEq ) THEN
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
    CALL SparseID( ID_1  , nDim )
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
    CALL DAXPY_sparse( f0 , BAT , Rate , Y_e )
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

    CALL DAXPY_sparse( f1 , BAT , Rate , Y_e )
 
    DfDt  = ( f1 - f0 ) / tdel
    CALL DAXPY_sparse( Tmp , Jac , f0 , zeros )
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
    !   - Temp........... actual Temperatur (optional for TempEq)
    !
    REAL(RealKind), INTENT(IN) :: Y0(nDIM)
    REAL(RealKind), INTENT(IN) :: t, h
    TYPE(RosenbrockMethod_T)   :: RCo
    !--------------------------------------------------------
    ! Output:
    !   - Ynew........... new concentratinos 
    !   - err............ error calc with embedded formula.
    !   - TempNew........ new temperature (optional for TempEq)
    !
    REAL(RealKind), INTENT(OUT) :: YNew(nDIM)
    REAL(RealKind), INTENT(OUT)   :: err
    INTEGER       , INTENT(OUT)   :: errind(1,1)
    !-------------------------------------------------------
    ! TemporarY variables:
    !
    REAL(RealKind), DIMENSION(nDIM)            :: Y,  Yhat, fRhs
    REAL(RealKind), DIMENSION(nspc)            :: Yrh, U, UMat, dUdT, dCdt, d2UdT2 
    REAL(RealKind), DIMENSION(nspc)            :: Jac_CT, Jac_TC
    REAL(RealKind), DIMENSION(nDIMex)          :: bb
    !
    REAL(RealKind) :: k( nDIM , RCo%nStage )
    !
    REAL(RealKind) :: Tarr(8)
    REAL(RealKind) :: Rate(neq), rRate(neq)
    REAL(RealKind) :: DRatedT(neq)        
    !REAL(RealKind) :: C(nspc)       ! molar heat capacities at constant pressure
    !REAL(RealKind) :: H_e(nspc)       ! the standardstate molar enthalpY
    !REAL(RealKind) :: S(nspc)       ! standard-state entropY at 298 K
    !
    !REAL(RealKind) :: dHdT(nspc)    ! EnthaplY derivative in dT [J/mol/K^2]
    !REAL(RealKind) :: dGdT(nspc)    ! Gibbs potential derivative in dT [J/mol/K^2]
    REAL(RealKind) :: dTdt, Jac_TT
      
    REAL(RealKind) :: tt
    REAL(RealKind) :: cv ! mass average mixture specific heat at constant volume
    REAL(RealKind) :: dcvdT ! mass average mixture specific heat at constant volume
    REAL(RealKind) :: X
    REAL(RealKind) :: Press
    !
    ! fuer verlgeich mit speedchem, andere spc reihenfolge
    !
    INTEGER :: iStg, jStg, i          ! increments
    LOGICAL :: dprint=.false.
    INTEGER :: iprnt

    dprint = DebugPrint   !init run
    !print*, nDIM
    !print*, Y0
    !stop
    !
    iprnt = MinVAL((/ 5 , neq , nspc /))

    ! Initial settings
    k(:,:)  = ZERO
    fRhs(:) = ZERO
    Rate(:) = ZERO
    rRate(:)= ZERO
    bb(:)   = ZERO
    Y(:)    = Y0
    !

    IF (dprint) THEN
      print*, ' '
      print*, '******************************************************************************'
    END IF

    !********************************************************************************
    !    _   _             _         _          __  __         _          _       
    !   | | | | _ __    __| |  __ _ | |_  ___  |  \/  |  __ _ | |_  _ __ (_)__ __
    !   | | | || '_ \  / _` | / _` || __|/ _ \ | |\/| | / _` || __|| '__|| |\ \/ /
    !   | |_| || |_) || (_| || (_| || |_|  __/ | |  | || (_| || |_ | |   | | >  < 
    !    \___/ | .__/  \__,_| \__,_| \__|\___| |_|  |_| \__,_| \__||_|   |_|/_/\_\
    !          |_|                                                                
    !
    !********************************************************************************
   
    ! HIER UNBEDINGT RATE MIT Y0 ALS INPUT
    IF ( .NOT.TempEq ) THEN
      Y   = MAX( ABS(Y0)  , eps ) * SIGN( ONE , Y0 )  ! concentrations =/= 0
      Yrh = Y(1:nspc) / h
      IF (CLASSIC)  CALL Rates( t, Y, Rate, DRatedT )     
      IF (EXTENDED) CALL Rates( t, Y0, Rate, DRatedT )     
      Rate  = MAX( ABS(Rate) , eps ) * SIGN( ONE , Rate )    ! reaction rates =/= 0
      rRate = ONE / Rate
    ELSE
      Yrh = Y0 / h
      CALL Rates( t, Y0, Rate, DRatedT )     
    END IF
    
    ! HIER  NICHT VON Y(:nspc) AUF Y0(:nspc) ÄNDERN
    !Yrh  = Y(1:nspc) / h

    IF ( TempEq ) THEN
      !                          OUT:      IN:
      CALL UpdateTempArray     ( Tarr    , Y0(nDIM) )       
      CALL InternalEnergy      ( U       , Tarr)             
      CALL DiffInternalEnergy  ( dUdT    , Tarr)              
      CALL Diff2InternalEnergy ( d2UdT2  , Tarr)
      CALL MassAveMixSpecHeat  ( cv      , dUdT    , MoleConc=Y0(1:nspc) , rho=rho)
      CALL MassAveMixSpecHeat  ( dcvdT   , d2UdT2  , MoleConc=Y0(1:nspc) , rho=rho)
      CALL DAXPY_sparse        ( dCdt    , BAT   , Rate , Y_e )

      dTdt    = - SUM( U * dCdt) * rRho/cv
      UMat    = RCo%ga*dcvdT*dTdt*Y0(1:nspc) + U*Yrh
      dRatedT = RCo%ga*dRatedT
      X       = cv/(h*rRho) + RCo%ga/cv*dcvdT*dTdt + RCo%ga*SUM(dUdT*dCdt) 
      !
      !
      IF (dprint) THEN
        WRITE(*,*) '------------------------------------------------------------------------------'
        WRITE(*,*) '|    Combustion                                                              |'
        WRITE(*,*) '------------------------------------------------------------------------------'
        WRITE(*,'(A,E23.16,A)') 'debug     Temperature =  ',Tarr(1)   ,'   [K]'
        WRITE(*,'(A,E23.16,A)') 'debug     c_v         =  ',cv        ,'   [J/kg/K]'
        WRITE(*,'(A,E23.16,A,E23.16,A)') 'debug     Density     =  ',rho       ,'   [kg/m3]',rRho       ,'   [kg/m3]'
        WRITE(*,'(A,E23.16,A)') 'debug     SUM(U)      =  ',SUM(U)    ,'   [J/mol/K]'
        WRITE(*,'(A,E23.16,A)') 'debug     X in Matrix =  ',X         ,'   [???]'
        WRITE(*,'(A,E23.16,A)') 'debug     SUM(dCdt)   =  ',SUM(dCdt) ,'   [mol/cm3/sec]'
        WRITE(*,'(A,E23.16,A)') 'debug     dTdt        =  ',dTdt      ,'   [K/sec]'
        WRITE(*,*) '------------------------------------------------------------------------------'
        WRITE(*,*) ''
        do i=1,nspc
          WRITE(*,'(A,I5,A,E22.14,A3,A)') 'debug     dCdt(',i,') = ',dCdt(SCperm(i)),'   ',y_name(SCperm(i))
        end do
        WRITE(*,'(A,E22.14,A3,A)') 'debug      dTdt(last) = ',dTdt,'   ','Temperature'
        WRITE(*,*) ''
        !stop 'debug ros'
      END IF
      !stop
      !
    END IF
   
    ! --- Update matrix procedure
    IF ( CLASSIC ) THEN
    
      ! classic case needs to calculate the Jacobian first
      TimeJacobianA = MPI_WTIME()

      ! d(dcdt)/dc
      CALL Jacobian_CC(  Jac_CC , BAT  , A , Rate , Y )

      IF ( TempEq ) THEN
        CALL Jacobian_CT( Jac_CT , BAT , Rate , DRatedT )
        CALL Jacobian_TC( Jac_TC , Jac_CC , cv , dUdT , dTdt , U , rRho)
        CALL Jacobian_TT( Jac_TT , Jac_CT , cv , dcvdT , dTdt , dUdT , dCdt , U , rRho)
        !
        CALL Miter_Classic( Miter , h , RCo%ga , Jac_CC , Jac_TC , Jac_CT , Jac_TT )
      ELSE
        CALL Miter_Classic( Miter , h , RCo%ga , Jac_CC )
      END IF
      Output%npds = Output%npds + 1

      IF ( useSparseLU ) THEN
        CALL SetLUvaluesCL( LU_Miter , Miter , LU_Perm )
      END IF
      TimeJac = TimeJac + (MPI_WTIME()-TimeJacobianA)

    ELSE !IF ( EXTENDED )
      ! 
      ! g = gamma (Rosenbrock method)
      !          _                                         _
      !         |              |              |             |
      !         |     rRate    |      g*A     |  g*dRatedT  |
      !         |              |              |             |
      !         |--------------+--------------+-------------|
      !         |              |              |             |
      ! miter = |      BAT     |      Yrh     |      0      |
      !         |              |              |             |
      !         !--------------+--------------+-------------|
      !         |              |              |             |
      !         |       0      |     UMat     |      X      |
      !         |_             |              |            _|
      !
      !
      rRate  = ONE / Rate
      IF ( useSparseLU ) THEN
        CALL SetLUvaluesEX ( LU_Miter, rRate , Yrh, DRatedT , UMat , X , LUvalsFix)
      ELSE
        CALL SetLUvaluesEX ( Miter, rRate , Yrh, DRatedT , UMat , X )
      END IF
    END IF
    !
    !                                  _      
    !  _ __  _   _ _ __ ___   ___ _ __(_) ___ 
    ! | '_ \| | | | '_ ` _ \ / _ \ '__| |/ __|
    ! | | | | |_| | | | | | |  __/ |  | | (__ 
    ! |_| |_|\__,_|_| |_| |_|\___|_|  |_|\___|
    !         _                                          _ _   _             
    !      __| | ___  ___ ___  _ __ ___  _ __   ___  ___(_) |_(_) ___  _ __  
    !     / _` |/ _ \/ __/ _ \| '_ ` _ \| '_ \ / _ \/ __| | __| |/ _ \| '_ \ 
    !    | (_| |  __/ (_| (_) | | | | | | |_) | (_) \__ \ | |_| | (_) | | | |
    !     \__,_|\___|\___\___/|_| |_| |_| .__/ \___/|___/_|\__|_|\___/|_| |_|
    !                                   |_|                                   
    ! --- LU - Decomposition ---
    timerStart  = MPI_WTIME()
    IF ( useSparseLU ) THEN
      CALL SparseLU( LU_Miter )
    ELSE
      CALL MumpsLU( Miter%Val )
    END IF
    TimeFac     = TimeFac + (MPI_WTIME()-timerStart)

    IF (dprint) THEN
      print*, '------------------------------------------------------------------------------'
      print*, '|    Rosenbrock Input                                                        |'
      print*, '------------------------------------------------------------------------------'
      print*, 'debug     Stepsize       =  ',h
      print*, 'debug     Time           =  ',t
      print*, 'debug     SUM(Y0)        =  ',SUM(Y0)
      print*, 'debug     SUM(Miter%val) =  ',SUM(Miter%val)
      IF( useSparseLU ) THEN
        print*, 'debug    SUM(LU%val) vor = ', SUM(LU_Miter%val)
      END IF
      print*, '------------------------------------------------------------------------------'
      print*, ''
      print*, ''
      print*, '------------------------------------------------------------------------------'
      print*, '| Before solving Ax=b:  iStage                   b                           |'
      print*, '------------------------------------------------------------------------------'
      print*, ''
    END IF

    !****************************************************************************************
    !   ____    ___ __        __          _____  _                    ____   _               
    !  |  _ \  / _ \\ \      / /         |_   _|(_) _ __ ___    ___  / ___| | |_  ___  _ __  
    !  | |_) || | | |\ \ /\ / /   _____    | |  | || '_ ` _ \  / _ \ \___ \ | __|/ _ \| '_ \ 
    !  |  _ < | |_| | \ V  V /   |_____|   | |  | || | | | | ||  __/  ___) || |_|  __/| |_) |
    !  |_| \_\ \___/   \_/\_/              |_|  |_||_| |_| |_| \___| |____/  \__|\___|| .__/ 
    !                                                                                 |_|    
    !****************************************************************************************
    
    LOOP_n_STAGES:  DO iStg = 1 , RCo%nStage

      IF ( iStg==1 ) THEN

        IF ( EXTENDED ) THEN
          bb( 1     : neq ) = mONE 
          bb( neq+1 : nsr ) = Y_e 
          IF ( TempEq ) bb(nDIMex)  = ZERO
          IF (dprint) WRITE(*,'(A25,I4,A3,*(E15.8,2X))') ' ',iStg,'   ',bb(1:iprnt)
        END IF

      ELSE ! iStage > 1 ==> Update time and concentration

        tt  = t + RCo%Asum(iStg) * h
        Y   = Y0

        DO jStg = 1 , iStg
          Y = Y + RCo%a(iStg,jStg) * k(:,jStg)
        END DO
        
        ! Update Rates at  (t + SumA*h) , and  (Y + A*)k
        CALL Rates( tt , Y , Rate , DRatedT )

        IF (dprint) THEN
          print*, ''
          print*, 'debug::         SUM(concentration) at (Y + a*k)  = ',SUM(Y)
          print*, 'debug::  SUM(Rates) at (t + SumA*h),  (Y + a*k)  = ',SUM(Rate)
        END IF
      END IF
      
      !--- Calculate the right hand side of the linear System
      IF ( CLASSIC ) THEN

        CALL DAXPY_sparse( dCdt , BAT , Rate , Y_e )           
        fRhs(1:nspc) =  h * dCdt

        IF (TempEq) &
        fRhs( nDIM ) = - h * SUM(U*dCdt) * rRho / cv

        DO jStg = 1 , iStg-1
          fRhs  = fRhs + RCo%C(iStg,jStg) * k(:,jStg)
        END DO

        IF (dprint) WRITE(*,'(A25,I4,A3,*(E15.8,2X))') ' ',iStg,'   ',fRhs( 1:iprnt)

      ELSE !IF ( EXTENDED ) THEN

        IF ( iStg/=1 ) THEN

          fRhs = ZERO

          DO jStg = 1 , iStg-1
            fRhs(1:nspc) = fRhs(1:nspc) + RCo%C(iStg,jStg)*k(1:nspc,jStg)
            IF (TempEq)                                                &
            fRhs(nDIM)   = fRhs(nDIM) + RCo%C(iStg,jStg)*                  &
            &               ( cv/rRho*k(nDIM,jStg) + SUM(U*k(1:nspc,jStg)) )
          END DO

          ! right hand side of the extended linear system

          bb( 1      : neq )     = -rRate * Rate
          bb( neq+1  : nsr )     = Y_e + fRhs(1:nspc)/h
          IF (TempEq) bb(nDIMex) = fRhs(nDIM)/h

          IF (dprint) WRITE(*,'(A25,I4,A3,*(E15.8,2X))') ' ',iStg,'   ',bb(1:iprnt)
        END IF

      END IF
      
      timerStart  = MPI_WTIME()
      SOLVE_LINEAR_SYSTEM: IF ( useSparseLU ) THEN  

        IF ( CLASSIC ) THEN

          CALL SolveSparse( LU_Miter , fRhs )
          k( 1:nDIM , iStg ) = fRhs

        ELSE !IF ( EXTENDED ) THEN

          CALL SolveSparse( LU_Miter , bb)
          k( 1:nspc , iStg ) = Y0(1:nspc) * bb(neq+1:nsr)
          IF ( TempEq ) &
          k(  nDIM  , iStg ) = bb(nDIMex)
          
        END IF

      ELSE

        IF ( CLASSIC ) THEN

          CALL MumpsSolve( fRhs )
          k( 1:nDIM , iStg ) = Mumps_Par%RHS(1:nDIM)    

        ELSE !IF ( EXTENDED ) THEN

          CALL MumpsSolve( bb )
          k( 1:nspc , iStg ) = Y0(1:nspc) * Mumps_Par%RHS(neq+1:nsr)
          IF ( TempEq ) &
          k(  nDIM  , iStg ) = Mumps_Par%RHS(nDIMex)

        END IF

      END IF SOLVE_LINEAR_SYSTEM
      TimeSolve   = TimeSolve + (MPI_WTIME()-timerStart)

    END DO  LOOP_n_STAGES
    IF (dprint) THEN
      print*, ''
      print*, '------------------------------------------------------------------------------'
      print*, '| After solving Ax=b:  Species         k( iSpc , : )                         |'
      print*, '------------------------------------------------------------------------------'
      print*, ''
      do istg=1,nspc
        WRITE(*,'(A9,I5,A25,A3,*(E15.8,2X))') 'debug::  ',istg, TRIM(y_name(istg)),'   ',k( istg , : )
      end do
      IF (TempEq) WRITE(*,'(A9,I5,A25,A3,*(E15.8,2X))') 'debug::  ',nDIM,'Temperature','   ',k( nDIM , : )
      print*, '------------------------------------------------------------------------------'
    END IF

    
    !--- Update Concentrations (+Temperatur)

    YNew = Y0
    YHat = Y0

    DO jStg = 1 , RCo%nStage
      YNew = YNew +  RCo%m(jStg) * k(:,jStg)! new Y vector
      YHat = YHat + RCo%me(jStg) * k(:,jStg)! embedded formula for err calc ord-1
    END DO
    !
    IF (dprint) THEN
      print*, ''
      print*, '------------------------------------------------------------------------------'
      print*, '| After Ros step:   Y_Old                     Y_New              Species name|'
      print*, '------------------------------------------------------------------------------'
      do istg=1,nspc
        print*,'debug::  ', Y0(istg),Ynew(iStg),'   ', TRIM(y_name(istg))
      end do
      IF (TempEq) print*,'debug::  ', Y0(nDIM),Ynew(nDIM), '   Temperature'
      print*, '------------------------------------------------------------------------------'
      print*, ''
    END IF

    !***********************************************************************************************
    !   _____                           _____       _    _                    __               
    !  | ____| _ __  _ __  ___   _ __  | ____| ___ | |_ (_) _ __ ___    __ _ | |_  (_)  ___   _ __  
    !  |  _|  | '__|| '__|/ _ \ | '__| |  _|  / __|| __|| || '_ ` _ \  / _` ||  __|| | / _ \ | '_ \ 
    !  | |___ | |   | |  | (_) || |    | |___ \__ \| |_ | || | | | | || (_| || |_  | || (_) || | | |
    !  |_____||_|   |_|   \___/ |_|    |_____||___/ \__||_||_| |_| |_| \__,_| \__| |_| \___/ |_| |_|
    !                                                                                              
    !***********************************************************************************************
    CALL ERROR( err , errind , YNew , YHat , ATolAll , RTolROW , t )
   
    !print*, 't,h, kvecs = ', t,h,SUM(k(:,1))
    !stop
    
    IF (dprint) THEN
      print*,'debug::     Error     =  ', err, '  Error index  =  ', errind
      print*, '------------------------------------------------------------------------------'
      print*, ''
      print*, ' Press ENTER to calculate next step '
      read(*,*) 
    END IF
   
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
      err = MAXVAL( En_Values )     ! maximum norm
    END IF
    !
  END SUBROUTINE ERROR

END MODULE Rosenbrock_Mod
