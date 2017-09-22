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
  USE Mumps_Mod
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
    REAL(dp) :: ga                            ! Diagonalentry gamma
    REAL(dp) :: pow                           ! needed for sitepsize control pow=1/nstage
    REAL(dp), ALLOCATABLE :: Asum(:)          ! Row sum of A
    REAL(dp), ALLOCATABLE :: Alpha(:,:)       ! Propagation table, strictly lower triangular
    REAL(dp), ALLOCATABLE :: a(:,:)           ! Propagation table, strictly lower triangular (converted Alpha)
    REAL(dp), ALLOCATABLE :: Gamma(:,:)       ! Stage table, lower triangular with nonzero diagonal
    REAL(dp), ALLOCATABLE :: iGamma(:,:)      ! inverse Stage table
    REAL(dp), ALLOCATABLE :: C(:,:)           ! Stage table, lower triangular with nonzero diagonal (converted Gamma)
    REAL(dp), ALLOCATABLE :: B(:)             ! Step completion table
    REAL(dp), ALLOCATABLE :: m(:)             ! Step completion table(converted B)
    REAL(dp), ALLOCATABLE :: Be(:)            ! Step completion table for embedded method of order one less
    REAL(dp), ALLOCATABLE :: me(:)            ! Step completion table for embedded method of order one less (converted Be)
    REAL(dp), ALLOCATABLE :: binterpt(:,:)    ! Dense output formula
  END TYPE RosenbrockMethod_T

!  TYPE Out
!    REAL(dp), ALLOCATABLE :: y(:)    ! y-vector at Tend
!    INTEGER :: nsteps     = 0              ! # succ. steps
!    INTEGER :: nfailed    = 0              ! # failed steps
!    INTEGER :: nRateEvals = 0              ! # Rate evaluation
!    INTEGER :: npds       = 0              ! # Jacobian evaluation
!    INTEGER :: ndecomps   = 0              ! # LU factorisation
!    INTEGER :: nsolves    = 0              ! # solved lin algebra
!    REAL(dp) :: Ttimestep = 0.0d0  ! mean Time for one ROW step
!  END TYPE Out
!
!  TYPE(Out) :: Output
  TYPE(RosenbrockMethod_T) :: ROS

  REAL(dp), PRIVATE :: timerStart
  REAL(dp), ALLOCATABLE :: LUvalsFix(:)

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
    REAL(dp), ALLOCATABLE :: ID(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:)
    INTEGER :: INFO

    CHARACTER(20) :: tmethod
    !
    INTEGER :: i
    !
    IF ( TRIM(method(1:7)) == 'bwEuler') THEN
      tmethod = ADJUSTL(method)
    ELSE
      tmethod = TRIM(method(INDEX(method,'/')+1:INDEX(method,'.')-1))
    END IF

    SELECT CASE (tmethod)
      CASE ('bwEuler')       
        INCLUDE 'METHODS/bwEuler.ros'
      CASE ('Ros2AMF')       
        INCLUDE 'METHODS/Ros2AMF.ros'
      CASE ('Ros3w')         
        INCLUDE 'METHODS/Ros3w.ros'
      CASE ('Ros3Pw')        
        INCLUDE 'METHODS/Ros3Pw.ros'
      CASE ('Ros34PW1a')     
        INCLUDE 'METHODS/Ros34PW1a.ros'
      CASE ('Ros34PW2')      
        INCLUDE 'METHODS/Ros34PW2.ros'
      CASE ('Ros34PW3')      
        INCLUDE 'METHODS/Ros34PW3.ros'
      CASE ('TSRosW2P')      
        INCLUDE 'METHODS/TSRosW2P.ros'
      CASE ('TSRosW2M')      
        INCLUDE 'METHODS/TSRosW2M.ros'
      CASE ('TSRosWRA3PW')   
        INCLUDE 'METHODS/TSRosWRA3PW.ros'
      CASE ('TSRosWRA34PW2') 
        INCLUDE 'METHODS/TSRosWRA34PW2.ros'
      CASE ('TSRosWRodas3')  
        INCLUDE 'METHODS/TSRosWRodas3.ros'
      CASE ('TSRosWSandu3')  
        INCLUDE 'METHODS/TSRosWSandu3.ros'
      CASE DEFAULT
        IF (MPI_ID==0) WRITE(*,*) '    Unknown Method:  ',method
        IF (MPI_ID==0) WRITE(*,*) '    Use TSRosWRodas3 instead.'
        INCLUDE 'METHODS/TSRosWRodas3.ros'
    END SELECT
    !INCLUDE RosenbrockMethod
    !
    ! converting the butcher tableau 
    ! automatic transformation to avoid mat*vec in ROW methode
    ROS%pow=ONE/(ROS%Order+ONE)
    !ROS%pow=1.0d0/ROS%nStage
    ALLOCATE(ROS%iGamma(ROS%nStage,ROS%nStage))
    ROS%iGamma=ZERO
    ALLOCATE(ID(ROS%nStage,ROS%nStage))
    ID=ZERO
    DO i=1,ROS%nStage
      ROS%iGamma(i,i)=ONE
      ID(i,i)=ONE
    END DO
    !  
    ! calculate the inverse matrix of gamma
    !
    ! CAUTION ROS%Gamma (IN) =/= ROS%Gamma (OUT) !!!
    !
    ALLOCATE(IPIV(ROS%nStage))
    CALL dgesv(  ROS%nStage,     &        ! # linear eqations
    &            ROS%nStage,     &        ! # RHS (coloums)
    &            ROS%Gamma,      &        ! Matrix A of A*A^(-1)=ID
    &            ROS%nStage,     &        ! leading dimension of A (nStage)
    &            IPIV,           &        ! pivot indices of dimension nStage
    &            ROS%iGamma,     &        ! Matrix ID of A*A^(-1)=ID
    &            ROS%nStage,     &        ! leading dimension of RHS
    &            INFO)                    ! INFO (integer) if INFO=0 succsessfull	
    !
    IF (INFO/=0.AND.MPI_ID==0) WRITE(*,*) 'Error while calc row-method parameter'
    !       
    ALLOCATE(ROS%a(ROS%nStage,ROS%nStage))
    ROS%a=ZERO
    ROS%a=ROS%ga*MATMUL(ROS%Alpha, ROS%iGamma)
    !  
    ALLOCATE(ROS%C(ROS%nStage,ROS%nStage))
    ROS%C=ZERO
    ROS%C=ID-ROS%ga*ROS%iGamma
    FORALL (i=1:ROS%nStage) ROS%C(i,i)=ZERO
    !  
    ALLOCATE(ROS%m(ROS%nStage))
    ROS%B=MATMUL(ROS%B, ROS%iGamma)
    ROS%m=ROS%ga*ROS%B(:)
    !
    IF (.NOT.ROS%nStage==1) THEN
      ALLOCATE(ROS%me(ROS%nStage))
      ROS%Be=MATMUL(ROS%Be,ROS%iGamma)
      ROS%me=ROS%ga*ROS%Be(:)
    END IF
    !
    DEALLOCATE(ID)
    DEALLOCATE(IPIV)
  END SUBROUTINE SetRosenbrockMethod
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
    REAL(dp), INTENT(IN) :: t, pow
    REAL(dp), INTENT(IN) :: Y(nspc)
    TYPE(CSR_Matrix_T), INTENT(IN) :: Jac
    REAL(dp) :: Rate(neq)
    REAL(dp) :: DRatedT(neq)     ! part. derv. rate over temperatur vector
    !-------------------------------------------------
    ! Output:
    !        - initial step size
    REAL(dp)  , INTENT(OUT) :: absh
    REAL(dp) :: h, hmin 
    REAL(dp) :: f0(nspc)
    !-------------------------------------------------
    !
    ! Temp vars:
    REAL(dp) :: tdel, rh
    REAL(dp), DIMENSION(nspc) ::  wt, DfDt, Tmp, f1, zeros
    REAL(dp) :: sqrteps=SQRT(eps)
    ! DEBUG
    INTEGER :: i  

    zeros = ZERO    

    ! hmin is a small number such that t + hmin is clearlY different from t in
    ! the working precision, but with this definition, it is 0 if t = 0.
    hmin  = minStp

    !---- Compute an initial step size h using Yp=Y'(t) 
    f0 = DAXPY_sparse( BAT , Rate , Y_e )
    !print*, 'debug:: sum(bat),rate =', SUM(BAT%val),SUM(rate)
    !print*, 'debug:: sum(f0)=', SUM(f0)
    !stop
    wt    = MAX( ABS(Y) , ThresholdStepSizeControl(1:nspc) )
    rh    = ( 1.25_dp * MAXVAL( ABS(f0(:)/wt(:)) ) )/(RTolRow**pow)
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

    IF (Teq) THEN
      CALL ReactionRatesAndDerivative_ChemKin( t+tdel , Y , Rate , DRatedT )
    ELSE
      CALL ReactionRates_Tropos( t+tdel , Y , Rate )
    END IF
    Out%nRateEvals = Out%nRateEvals + 1

    f1 = DAXPY_sparse( BAT , Rate , Y_e )
    DfDt  = ( f1 - f0 ) / tdel
    
    tmp = DAXPY_sparse( Jac , f0 , zeros )
    DfDt  = DfDt  + Tmp
  
    rh    = 1.25_dp * SQRT( rTWO * MAXVAL( ABS(DfDt/wt) ) ) / RTolRow**pow
   
    absh  = MIN( maxStp , Tspan(2)-Tspan(1) )
    IF ( absh * rh > ONE )  absh = ONE / rh
    absh  = MAX( absh , hmin )
    
  END SUBROUTINE InitialStepSize
  !
  !
  !=========================================================================
  !    Subroutine Rosenbrock-Method universal for classic and extended case
  !=========================================================================
  SUBROUTINE Rosenbrock(YNew,err,ierr,Y0,t,h,Euler)
    USE issa
    USE mo_IO
    !--------------------------------------------------------
    ! Input:
    !   - Y0............. concentrations at Time = t
    !   - t.............. time
    !   - h.............. step size
    !   - RCo............ Rosenbrock method
    !   - Temp........... Temperatur at Time = t (optional for Teq)
    !   - Euler.......... if .true. --> use backward Euler method
    !
    REAL(dp),          INTENT(IN) :: Y0(nDIM)
    REAL(dp),          INTENT(IN) :: t, h
    !TYPE(RosenbrockMethod_T)      :: RCo
    LOGICAL, OPTIONAL, INTENT(IN) :: Euler
    !--------------------------------------------------------
    ! Output:
    !   - Ynew........... new concentratinos 
    !   - err............ error calc with embedded formula.
    !   - TempNew........ new temperature (optional for Teq)
    !
    REAL(dp), INTENT(OUT) :: YNew(nDIM)
    REAL(dp), INTENT(OUT) :: err
    INTEGER , INTENT(OUT) :: ierr(1,1)
    !-------------------------------------------------------
    ! TemporarY variables:
    !
    REAL(dp), DIMENSION(nDIM)   :: Y,  Yhat, fRhs
    REAL(dp), DIMENSION(nspc)   :: Yrh, U, UMat, dUdT, dCdt, d2UdT2 
    REAL(dp), DIMENSION(nspc)   :: Jac_CT, Jac_TC
    REAL(dp), DIMENSION(nDIMex) :: bb
    !
    REAL(dp) :: k( nDIM , ROS%nStage )
    !
    REAL(dp) :: Tarr(10)
    REAL(dp) :: Rate(neq), rRate(neq), Rate_t(neq)
    REAL(dp) :: DRatedT(neq)        
    REAL(dp) :: dTdt, Jac_TT
      
    REAL(dp) :: tt
    REAL(dp) :: cv ! mass average mixture specific heat at constant volume
    REAL(dp) :: dcvdT ! mass average mixture specific heat at constant volume
    REAL(dp) :: X
    REAL(dp) :: Press

    REAL(dp) :: TimeRhsCalc0
    !
    ! fuer verlgeich mit speedchem, andere spc reihenfolge
    !
    INTEGER :: iStg, jStg, i, j          ! increments
    LOGICAL :: dprint=.false.
    INTEGER :: iprnt

    REAL(dp) :: CM(nDIM)

    dprint = DebugPrint   !init run
    !
    !iprnt = MinVAL((/ 5 , neq , nspc /))

    ! Initial settings
    k       = ZERO
    fRhs    = ZERO
    Rate    = ZERO
    rRate   = ZERO
    bb      = ZERO
    Y       = Y0
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
    IF ( .NOT.Teq ) THEN
      Y   = MAX( ABS(Y0)  , eps ) * SIGN( ONE , Y0 )  ! concentrations =/= 0
      Yrh = Y(1:nspc) / h
      CALL ReactionRates_Tropos( t, Y, Rate )
      Rate  = MAX( ABS(Rate) , eps ) * SIGN( ONE , Rate )    ! |r| >= eps
      rRate = ONE / Rate
    ELSE
      Yrh = Y0(1:nspc) / h
      CALL ReactionRatesAndDerivative_ChemKin( t, Y0, Rate, DRatedT )     
    END IF
    Rate_t = Rate

    IF ( Teq ) THEN
      Tarr = UpdateTempArray   ( Y0(nDIM) )       

      !                          OUT:      IN:
      CALL InternalEnergy      ( U       , Tarr)             
      CALL DiffInternalEnergy  ( dUdT    , Tarr)              
      CALL Diff2InternalEnergy ( d2UdT2  , Tarr)
      CALL MassAveMixSpecHeat  ( cv      , dUdT    , MoleConc=Y0(1:nspc) , rho=rho)
      CALL MassAveMixSpecHeat  ( dcvdT   , d2UdT2  , MoleConc=Y0(1:nspc) , rho=rho)

      dCdt = DAXPY_sparse( BAT   , Rate , Y_e )

      dTdt    = - SUM(U * dCdt) * rRho/cv
      UMat    = ROS%ga*dcvdT*dTdt*Y0(1:nspc) + U*Yrh
      dRatedT = ROS%ga*dRatedT
      X       = cv/(h*rRho) + ROS%ga/cv*dcvdT*dTdt + ROS%ga*SUM(dUdT*dCdt) 
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

      IF ( Teq ) THEN
        CALL Jacobian_CT( Jac_CT , BAT , Rate , DRatedT )
        CALL Jacobian_TC( Jac_TC , Jac_CC , cv , dUdT , dTdt , U , rRho)
        CALL Jacobian_TT( Jac_TT , Jac_CT , cv , dcvdT , dTdt , dUdT , dCdt , U , rRho)
        !
        CALL Miter_Classic( Miter , h , ROS%ga , Jac_CC , Jac_TC , Jac_CT , Jac_TT )
      ELSE
        CALL Miter_Classic( Miter , h , ROS%ga , Jac_CC )
      END IF
      Out%npds = Out%npds + 1

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
    timerStart = MPI_WTIME()
    IF ( useSparseLU ) THEN
      CALL SparseLU( LU_Miter )
    ELSE
      CALL MumpsLU( Miter%Val )
    END IF
    TimeFac = TimeFac + (MPI_WTIME()-timerStart)

    IF (dprint) THEN
      print*, '------------------------------------------------------------------------------'
      print*, '|    Rosenbrock Input                                                        |'
      print*, '------------------------------------------------------------------------------'
      print*, 'debug     Stepsize       =  ',h
      print*, 'debug     Time           =  ',t
      print*, 'debug     SUM(Y0)        =  ',SUM(Y0)
      print*, 'debug     SUM(Rate)      =  ',SUM(Rate)
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
    
    LOOP_n_STAGES:  DO iStg = 1 , ROS%nStage

      IF ( iStg==1 ) THEN

        IF ( EXTENDED ) THEN
          bb( 1     : neq ) = mONE 
          bb( neq+1 : nsr ) = Y_e 
          IF ( Teq ) bb(nDIMex)  = ZERO
          IF (dprint) WRITE(*,'(A25,I4,A3,*(E15.8,2X))') ' ',iStg,'   ',bb(1:iprnt)
        END IF

      ELSE ! iStage > 1 ==> Update time and concentration

        tt  = t + ROS%Asum(iStg) * h
        Y   = Y0

        DO jStg = 1 , iStg
          Y = Y + ROS%a(iStg,jStg) * k(:,jStg)
        END DO
        
        ! Update Rates at  (t + SumA*h) , and  (Y + A*)k
        IF (Teq) THEN
          CALL ReactionRatesAndDerivative_ChemKin( tt , Y , Rate , DRatedT )
        ELSE
          CALL ReactionRates_Tropos( tt , Y , Rate )
        END IF

        IF (dprint) THEN
          print*, ''
          print*, 'debug::         SUM(concentration) at (Y + a*k)  = ',SUM(Y)
          print*, 'debug::  SUM(Rates) at (t + SumA*h),  (Y + a*k)  = ',SUM(Rate)
        END IF
      END IF
      
      TimeRhsCalc0 = MPI_WTIME()
      !--- Calculate the right hand side of the linear System
      IF ( CLASSIC ) THEN

        dCdt = DAXPY_sparse( BAT , Rate , Y_e )           
        fRhs(1:nspc) =  h * dCdt
        IF (Teq) fRhs( nDIM ) = - h * SUM(U*dCdt) * rRho / cv

        DO jStg=1,iStg-1; fRhs = fRhs + ROS%C(iStg,jStg)*k(:,jStg); END DO

        IF (dprint) WRITE(*,'(A25,I4,A3,*(E15.8,2X))') ' ',iStg,'   ',fRhs( 1:iprnt)

      ELSE !IF ( EXTENDED ) THEN

        IF ( iStg/=1 ) THEN

          fRhs = ZERO

          DO jStg = 1 , iStg-1
            fRhs(1:nspc) = fRhs(1:nspc) + ROS%C(iStg,jStg)*k(1:nspc,jStg)
            IF (Teq)                                                       &
            fRhs(nDIM)   = fRhs(nDIM) + ROS%C(iStg,jStg)*                  &
            &               ( cv/rRho*k(nDIM,jStg) + SUM(U*k(1:nspc,jStg)) )
          END DO

          ! right hand side of the extended linear system

          bb( 1      : neq )  = -rRate * Rate
          bb( neq+1  : nsr )  = Y_e + fRhs(1:nspc)/h
          IF (Teq) bb(nDIMex) = fRhs(nDIM)/h

          IF (dprint) WRITE(*,'(A25,I4,A3,*(E15.8,2X))') ' ',iStg,'   ',bb(1:iprnt)
        END IF

      END IF
      TimeRhsCalc = TimeRhsCalc + MPI_WTIME() - TimeRhsCalc0
      
      timerStart  = MPI_WTIME()
      SOLVE_LINEAR_SYSTEM: IF ( useSparseLU ) THEN  

        IF ( CLASSIC ) THEN
          CALL SolveSparse( LU_Miter , fRhs )
          k( 1:nDIM , iStg ) = fRhs
        ELSE !IF ( EXTENDED ) THEN
          CALL SolveSparse( LU_Miter , bb)
          k( 1:nspc , iStg ) = Y0(1:nspc) * bb(neq+1:nsr)
          IF ( Teq ) &
          k(  nDIM  , iStg ) = bb(nDIMex)
        END IF

      ELSE

        IF ( CLASSIC ) THEN
          CALL MumpsSolve( fRhs )
          k( 1:nDIM , iStg ) = Mumps_Par%RHS(1:nDIM)    
        ELSE !IF ( EXTENDED ) THEN
          CALL MumpsSolve( bb )
          k( 1:nspc , iStg ) = Y0(1:nspc) * Mumps_Par%RHS(neq+1:nsr)
          IF ( Teq ) &
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
      IF (Teq) WRITE(*,'(A9,I5,A25,A3,*(E15.8,2X))') 'debug::  ',nDIM,'Temperature','   ',k( nDIM , : )
      print*, '------------------------------------------------------------------------------'
    END IF

    
    !--- Update Concentrations (+Temperatur)
    YNew = Y0
    DO jStg=1,ROS%nStage; YNew = YNew + ROS%m(jStg)*k(:,jStg); END DO
    
    IF (dprint) THEN
      print*, ''
      print*, '------------------------------------------------------------------------------'
      print*, '| After Ros step:   Y_Old                     Y_New              Species name|'
      print*, '------------------------------------------------------------------------------'
      do istg=1,nspc
        print*,'debug::  ', Y0(istg),Ynew(iStg),'   ', TRIM(y_name(istg))
      end do
      IF (Teq) print*,'debug::  ', Y0(nDIM),Ynew(nDIM), '   Temperature'
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

    IF (.NOT.EULER) THEN

      timerStart = MPI_WTIME()
      ! embedded formula for err calc ord-1
      YHat = Y0
      DO jStg=1,ROS%nStage; YHat = YHat + ROS%me(jStg)*k(:,jStg); END DO

      CALL ERROR( err , ierr , YNew , YHat , ATolAll , RTolROW , t )
      TimeErrCalc = TimeErrCalc + MPI_WTIME() - timerStart

      ! for analysis and reduction
      IF ( MPI_ID == 0 .AND. err < ONE ) THEN
        IF ( Lehmann ) integrated_rates = integrated_rates + Rate*h
        !CALL ISSA_iter( Rate_t , YNew , t , h )

        IF ( FluxAna ) THEN
          IF ( t - Tspan(1) >= StpFlux*REAL(iStpFlux,dp) ) THEN
            timerStart = MPI_WTIME()
            CALL StreamWriteFluxes(Rate_t,t,h)
            !CALL SequentialWriteFluxes(Rate_t,t,h)
            TimeFluxWrite = TimeFluxWrite + MPI_WTIME() - timerStart
          END IF
        END IF

      END IF

    END IF


   
    IF (dprint) THEN
      print*,'debug::     Error     =  ', err, '  Error index  =  ', ierr
      print*, '------------------------------------------------------------------------------'
      print*, ''
      print*, ' Press ENTER to calculate next step '
      read(*,*) 
    END IF

   
  END SUBROUTINE Rosenbrock


  SUBROUTINE ERROR(err,ierr,ynew,yhat,ATol,RTol,t)
    !
    REAL(dp) :: err
    REAL(dp), DIMENSION(:) :: ynew, yhat, ATol
    REAL(dp) :: RTol, t
    !
    REAL(dp) :: scalTol(nDIM), e_n(nDIM)
    INTEGER :: ierr(1,1)
    !
    scalTol   = ONE / (ATol + MAX( ABS(yhat) , ABS(ynew) ) * RTol )  ! scaling strategie
    e_n       = ABS( ynew - yhat ) * scalTol      ! local error est.
    ierr(1,1) = MAXLOC( e_n , 1 )           ! max error component
    !
    IF ( Error_Est == 2 ) THEN
      err = SUM( e_n*e_n ) * rNspc   ! euclikd norm
    ELSE
      err = MAXVAL( e_n )     ! maximum norm
    END IF
    !
  END SUBROUTINE ERROR


END MODULE Rosenbrock_Mod
