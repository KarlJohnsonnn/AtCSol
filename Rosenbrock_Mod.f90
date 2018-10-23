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
  USE Rates_Mod
  USE Control_Mod
  USE Reac_Mod
  USE CombustionInput_Mod
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
    tmethod = ADJUSTL(method(INDEX(method,'/')+1:INDEX(method,'.')-1))
    SELECT CASE (tmethod)
      CASE ('bwEuler')
        INCLUDE 'METHODS/bwEuler.ros'
      CASE ('Ros2AMF')       
        INCLUDE 'METHODS/Ros2AMF.ros'
      CASE ('Ros3w')         
        INCLUDE 'METHODS/Ros3w.ros'
      CASE ('Ros3Dw')        
        INCLUDE 'METHODS/Ros3Dw.ros'
      CASE ('Ros3Pw')        
        INCLUDE 'METHODS/Ros3Pw.ros'
      CASE ('Ros34PW1a')     
        INCLUDE 'METHODS/Ros34PW1a.ros'
      CASE ('Ros34PW2')      
        INCLUDE 'METHODS/Ros34PW2.ros'
      CASE ('Ros34PW3')      
        INCLUDE 'METHODS/Ros34PW3.ros'
      CASE ('Rodas3')  
        INCLUDE 'METHODS/Rodas3.ros'
      CASE ('TSRosW2P')      
        INCLUDE 'METHODS/TSRosW2P.ros'
      CASE ('TSRosW2M')      
        INCLUDE 'METHODS/TSRosW2M.ros'
      CASE ('TSRosWRA34PW2') 
        INCLUDE 'METHODS/TSRosWRA34PW2.ros'
      CASE ('TSRosWSandu3')  
        INCLUDE 'METHODS/TSRosWSandu3.ros'
      CASE DEFAULT
        WRITE(*,*) '    Unknown Method:  ',method
        WRITE(*,*) '    Use Rodas3 instead.'
        INCLUDE 'METHODS/Rodas3.ros'
    END SELECT
    !INCLUDE RosenbrockMethod
    !
    ! converting the butcher tableau 
    ! automatic transformation to avoid mat*vec in ROW methode
    ALLOCATE(ID(ROS%nStage,ROS%nStage) , ROS%iGamma(ROS%nStage,ROS%nStage))
    ROS%pow    = ONE / (ROS%Order+ONE)
    ROS%iGamma = ZERO
    ID         = ZERO
    DO i=1,ROS%nStage
      ROS%iGamma(i,i) = ONE
      ID(i,i) = ONE
    END DO

    ALLOCATE(RCo%Asum(RCo%nStage))
    DO i=1,RCo%nStage
      RCo%Asum(i) = SUM(RCo%Alpha(i,:))
    END DO
    
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
    IF ( INFO/= 0 ) WRITE(*,*) 'Error while calc row-method parameter'
    !       
    ALLOCATE(ROS%a(ROS%nStage,ROS%nStage))
    ROS%a = ZERO
    ROS%a = ROS%ga*MATMUL(ROS%Alpha, ROS%iGamma)
    !  
    ALLOCATE(ROS%C(ROS%nStage,ROS%nStage))
    ROS%C = ZERO
    ROS%C = ID - ROS%ga * ROS%iGamma
    FORALL (i=1:ROS%nStage) ROS%C(i,i)=ZERO
    !  
    ALLOCATE(ROS%m(ROS%nStage))
    ROS%B = MATMUL(ROS%B, ROS%iGamma)
    ROS%m = ROS%ga * ROS%B
    !
    IF (.NOT.ROS%nStage==1) THEN
      ALLOCATE(ROS%me(ROS%nStage))
      ROS%Be = MATMUL(ROS%Be,ROS%iGamma)
      ROS%me = ROS%ga * ROS%Be(:)
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
  FUNCTION InitialStepSize(Jac,Rate,t,Y_in,pow) RESULT(h)
    !------------------------------------------------- 
    ! Input:
    !        - public variables
    !        - Tspan 
    !        - Y0  ( initial vector )
    !        - Jacobian matrix
    TYPE(CSR_Matrix_T), INTENT(IN) :: Jac
    REAL(dp), INTENT(IN) :: t, pow
    REAL(dp), INTENT(IN) :: Y_in(:)
    REAL(dP), INTENT(INOUT) :: Rate(:)
    !REAL(dp) :: DRatedT(nr)     ! part. derv. rate over temperatur vector
    !-------------------------------------------------
    ! Output:
    !        - initial step size
    REAL(dp) :: h 
    !-------------------------------------------------
    !
    ! Temp vars:
    REAL(dp) :: tdel, rh, absh
    REAL(dp), DIMENSION(nDIM) :: Y
    REAL(dp), DIMENSION(nspc) :: wt, DfDt, Tmp, f0, f1, zeros, Emiss
    REAL(dp), DIMENSION(nr)   :: dRdT
    REAL(dp), ALLOCATABLE     :: thresh(:)
    ! DEBUG
    INTEGER :: i  

    zeros = ZERO

    Y = Y_in 

    ! hmin is a small number such that t + hmin is clearlY different from t in
    ! the working precision, but with this definition, it is 0 if t = 0.

    !---- Compute an initial step size h using Yp=Y'(t) 
    CALL UpdateEmission(Emiss,t)
    f0 = BAT * Rate + Emiss

    ALLOCATE( thresh(nDIM) )
    IF (  Combustion ) THEN
      thresh(iGs)  = AtolGas / RTolROW
      thresh(nDIM) = AtolTemp / RTolROW
    ELSE
      IF ( hasGasSpc  ) thresh(iGs) = AtolGas / RTolROW
      IF ( hasAquaSpc ) thresh(iAs) = AtolAqua / RTolROW
    END IF

    wt   = MAX( ABS(Y(1:nspc)) , thresh(1:nspc) )
    rh   = ( 1.25_dp * MAXVAL(ABS(f0/wt)) )/(RTolRow**pow)
    absh = MIN( maxStp , Tspan(2)-Tspan(1) )
    IF ( absh * rh > ONE )  absh = ONE / rh

    !---- Compute Y''(t) and a better initial step size
    h = absh
    tdel  = t + MIN( SQRT(eps) * MAX( ABS(t) , ABS(t+h) ) , absh )

    IF ( Combustion) THEN
      CALL ReactionRates( t+tdel , Y , Rate , dRdT )
    ELSE
      CALL ReactionRates( t+tdel , Y , Rate)
    END IF
    Out%nRateEvals = Out%nRateEvals + 1

    CALL UpdateEmission(Emiss,tdel)
    f1   = BAT * Rate + Emiss
    DfDt = ( f1 - f0 ) / tdel
    
    tmp  = Jac * f0 
    DfDt = DfDt + Tmp

    rh   = 1.25_dp * SQRT( rTWO * MAXVAL( ABS(DfDt/wt) ) ) / RTolRow**pow
  
    absh = MIN( maxStp , Tspan(2)-Tspan(1) )
    IF ( absh*rh > ONE )  absh = ONE/rh
    h = MAX( absh , minStp )
    
  END FUNCTION InitialStepSize
  !
  !
  !=========================================================================
  !    Subroutine Rosenbrock-Method universal for classic and extended case
  !=========================================================================
  SUBROUTINE Rosenbrock(YNew,err,ierr,Y0,t,h,Euler)
    USE ISSA_mod
    USE IO_Mod
    !--------------------------------------------------------
    ! Input:
    !   - Y0............. concentrations at Time = t
    !   - t.............. time
    !   - h.............. step size
    !   - RCo............ Rosenbrock method
    !   - Temp........... Temperatur at Time = t (optional for  Combustion)
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
    !   - TempNew........ new temperature (optional for  Combustion)
    !
    REAL(dp), INTENT(OUT) :: YNew(nDIM)
    REAL(dp), INTENT(OUT) :: err
    INTEGER , INTENT(OUT) :: ierr(1,1)
    !-------------------------------------------------------
    ! TemporarY variables:
    !
    REAL(dp), DIMENSION(nDIM)   :: Y,  Yhat, fRhs
    REAL(dp), DIMENSION(nspc)   :: Yrh, U, UMat, dUdT, dCdt, d2UdT2, Emiss
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


    ! Initial settings
    k       = ZERO
    fRhs    = ZERO
    Rate    = ZERO
    rRate   = ZERO
    bb      = ZERO
    Y       = Y0
    Emiss   = ZERO
    


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
    IF (  Combustion ) THEN
      Yrh = Y0(1:nspc) / h
      CALL ReactionRates( t, Y0, Rate, DRatedT )     
    ELSE
      Y   = MAX( ABS(Y0)  , eps ) * SIGN( ONE , Y0 )  ! concentrations =/= 0
      Yrh = Y(1:nspc) / h
      CALL ReactionRates( t, Y, Rate )
      Rate  = MAX( ABS(Rate) , eps ) * SIGN( ONE , Rate )    ! |r| >= eps
      rRate = ONE / Rate
    END IF
    Rate_t = Rate

    CALL UpdateEmission(Emiss,t)

    IF (  Combustion ) THEN
      Tarr = UpdateTempArray   ( Y0(nDIM) )       
      CALL InternalEnergy      ( U       , Tarr)             
      CALL DiffInternalEnergy  ( dUdT    , Tarr)              
      CALL Diff2InternalEnergy ( d2UdT2  , Tarr)
      CALL MassAveMixSpecHeat  ( cv      , dUdT    , MoleConc=Y0(1:nspc) , rho=rho)
      CALL MassAveMixSpecHeat  ( dcvdT   , d2UdT2  , MoleConc=Y0(1:nspc) , rho=rho)

      dCdT = BAT * Rate + Emiss

      dTdt    = - SUM(U * dCdt) * rRho/cv
      UMat    = ROS%ga*dcvdT*dTdt*Y0(1:nspc) + U*Yrh
      dRatedT = ROS%ga*dRatedT
      X       = cv/(h*rRho) + ROS%ga/cv*dcvdT*dTdt + ROS%ga*SUM(dUdT*dCdt) 
      !
      
    END IF
   
    ! --- Update matrix procedure
    IF ( CLASSIC ) THEN
    
      ! classic case needs to calculate the Jacobian first
      CALL CPU_Time(TimeJacobianA)

      ! d(dcdt)/dc
      CALL Jacobian_CC(  Jac_CC , BAT  , A , Rate , Y )

      IF (  Combustion ) THEN
        CALL Jacobian_CT( Jac_CT , BAT , Rate , DRatedT )
        CALL Jacobian_TC( Jac_TC , Jac_CC , cv , dUdT , dTdt , U , rRho)
        CALL Jacobian_TT( Jac_TT , Jac_CT , cv , dcvdT , dTdt , dUdT , dCdt , U , rRho)
        !
        CALL Miter_Classic( Miter , h , ROS%ga , Jac_CC , Jac_TC , Jac_CT , Jac_TT )
      ELSE
        CALL Miter_Classic( Miter , h , ROS%ga , Jac_CC )
      END IF
      Out%npds = Out%npds + 1

      CALL SetLUvaluesCL( LU_Miter , Miter , LU_Perm )
      CALL CPU_TIME(TimeJacobianE)
      TimeJac = TimeJac + (TimeJacobianE-TimeJacobianA)

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
      CALL SetLUvaluesEX ( LU_Miter, rRate , Yrh, DRatedT , UMat , X , LUvalsFix)
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
    CALL CPU_TIME(TimeFacA)
    CALL SparseLU(LU_Miter)
    CALL CPU_TIME(TimeFacE)
    TimeFac = TimeFac + (TimeFacE-TimeFacA)

    

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
          bb( neq+1 : nsr ) = Emiss 
          IF (  Combustion ) bb(nDIMex)  = ZERO
          
        END IF

      ELSE ! iStage > 1 ==> Update time and concentration

        tt  = t + ROS%Asum(iStg) * h
        Y   = Y0

        DO jStg = 1 , iStg
          Y = Y + ROS%a(iStg,jStg) * k(:,jStg)
        END DO
        
        ! Update Rates at  (t + SumA*h) , and  (Y + A*)k
        IF ( Combustion) THEN
          CALL ReactionRates( tt , Y , Rate , DRatedT )
        ELSE
          CALL ReactionRates( tt , Y , Rate )
        END IF

        CALL UpdateEmission(Emiss,tt)

      END IF

      !--- Calculate the right hand side of the linear System
      CALL CPU_TIME(TimeRhsA)
      IF ( CLASSIC ) THEN

        dCdt = BAT * Rate + Emiss
        fRhs(1:nspc) =  h * dCdt
        IF (Combustion) THEN
          Tarr = UpdateTempArray ( Y(nDIM) )       
          CALL InternalEnergy    ( U       , Tarr)  
          CALL DiffInternalEnergy( dUdT    , Tarr)              
          CALL MassAveMixSpecHeat( cv      , dUdT    , MoleConc=Y(1:nspc) , rho=rho)
          fRhs( nDIM ) = - h * SUM(U*dCdt) * rRho / cv
        END IF

        DO jStg=1,iStg-1; fRhs = fRhs + ROS%C(iStg,jStg)*k(:,jStg); END DO

      ELSE !IF ( EXTENDED ) THEN

        IF ( iStg/=1 ) THEN

          fRhs = ZERO

          DO jStg = 1 , iStg-1
            fRhs(1:nspc) = fRhs(1:nspc) + ROS%C(iStg,jStg)*k(1:nspc,jStg)
            
            IF ( Combustion) THEN
              Tarr = UpdateTempArray( Y(nDIM) )       
              CALL InternalEnergy( U , Tarr)    
              CALL DiffInternalEnergy( dUdT , Tarr)              
              CALL MassAveMixSpecHeat( cv   , dUdT , MoleConc=Y(1:nspc) , rho=rho)
              fRhs(nDIM)   = fRhs(nDIM) + ROS%C(iStg,jStg)*( cv/rRho*k(nDIM,jStg) + SUM(U*k(1:nspc,jStg)) )
            END IF
          END DO

          ! right hand side of the extended linear system
          bb( 1      : neq )  = -rRate * Rate
          bb( neq+1  : nsr )  = Emiss + fRhs(1:nspc)/h
          IF (Combustion) bb(nDIMex) = fRhs(nDIM)/h
        END IF

      END IF
      CALL CPU_TIME(TimeRhsE)
      TimeRhsCalc = TimeRhsCalc + TimeRhsE - TimeRhsA
      
      CALL CPU_TIME(TimeSolA)
      IF ( CLASSIC ) THEN
        CALL SolveSparse( LU_Miter , fRhs )
        k( 1:nDIM , iStg ) = fRhs
      ELSE !IF ( EXTENDED ) THEN
        CALL SolveSparse( LU_Miter , bb)
        k( 1:nspc , iStg ) = Y0(1:nspc) * bb(neq+1:nsr)
        IF (  Combustion ) &
        k(  nDIM  , iStg ) = bb(nDIMex)
      END IF
      CALL CPU_TIME(TimeSolE)
      TimeSolve   = TimeSolve + (TimeSolE-TimeSolA)

    END DO  LOOP_n_STAGES

    
    !--- Update Concentrations (+Temperatur)
    YNew = Y0
    DO jStg=1,ROS%nStage 
      YNew = YNew + ROS%m(jStg)*k(:,jStg)
    END DO

    !***********************************************************************************************
    !   _____                           _____       _    _                    __               
    !  | ____| _ __  _ __  ___   _ __  | ____| ___ | |_ (_) _ __ ___    __ _ | |_  (_)  ___   _ __  
    !  |  _|  | '__|| '__|/ _ \ | '__| |  _|  / __|| __|| || '_ ` _ \  / _` ||  __|| | / _ \ | '_ \ 
    !  | |___ | |   | |  | (_) || |    | |___ \__ \| |_ | || | | | | || (_| || |_  | || (_) || | | |
    !  |_____||_|   |_|   \___/ |_|    |_____||___/ \__||_||_| |_| |_| \__,_| \__| |_| \___/ |_| |_|
    !                                                                                              
    !***********************************************************************************************

    IF (.NOT.EULER) THEN

      CALL CPU_TIME(TimeErrCalcA)
      ! embedded formula for err calc ord-1
      YHat = Y0
      DO jStg=1,ROS%nStage
        YHat = YHat + ROS%me(jStg)*k(:,jStg) 
      END DO

      CALL ERROR( err , ierr , YNew , YHat , ATolAll , RTolROW , t )
      CALL CPU_TIME(TimeErrCalcE)
      TimeErrCalc = TimeErrCalc + TimeErrCalcE - TimeErrCalcA

    END IF
   
  END SUBROUTINE Rosenbrock


  SUBROUTINE ERROR(err,ierr,ynew,yhat,ATol,RTol,t)
    !
    REAL(dp) :: err
    REAL(dp), DIMENSION(:) :: ynew, yhat, ATol
    REAL(dp) :: RTol, t
    !
    REAL(dp) :: scalTol(nDIM), e_n(nDIM), ymax(nDIM)
    INTEGER :: ierr(1,1)
    !
    ymax      = MAX(ABS(yhat),ABS(ynew))
    scalTol   = ONE / ( ATol + ymax*RTol )  ! scaling strategie
    e_n       = ABS( ynew - yhat ) * scalTol      ! local error est.
    ierr(1,1) = MAXLOC( e_n , 1 )           ! max error component
    !
    IF ( Error_Est == 2 ) THEN
      !err = SUM( e_n*e_n ) 
      !err = SQRT(err) / rNspc

      err = SUM( e_n*e_n ) * rNspc   ! euclikd norm
    ELSE
      err = MAXVAL( e_n )     ! maximum norm
    END IF

    !err = MAX(err,1.0e-10_dp)
    !
  END SUBROUTINE ERROR


END MODULE Rosenbrock_Mod
