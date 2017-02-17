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
  REAL(RealKind), ALLOCATABLE :: ValCopy(:)

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
  SUBROUTINE SetRosenbrockArgs(Rate,y0,Tspan,Atol,RtolROW,pow)
    !--------------------------------------------------------------------
    ! Input:
    !   - y0
    !   - Tspan
    !   - AtolGas, RtolROW
    !   - pow
    REAL(RealKind), INTENT(IN) :: y0(:)       
    REAL(RealKind) :: Tspan(2)               
    REAL(RealKind) :: Atol(2)
    Real(RealKind) :: RtolROW
    Real(RealKind) :: pow
    !--------------------------------------------------------------------
    ! Output:
    !   -  public TYPE(IntArgs)
    REAL(RealKind) :: Rate(neq)                    ! rate vector
    !--------------------------------------------------------------------
    ! Temorary variables:
    !
    REAL(RealKind) :: Tstart,Tend                 ! start-, end-time
    REAL(RealKind) :: checkdir
    REAL(RealKind) :: one=1.0d0                  ! for Tdir check
    !
    Args%hTspan=ABS(Tspan(2)-Tspan(1))
    Args%Tend=Tspan(2)
    Args%nep=nspc
    ALLOCATE(Args%ATol(nspc))
    Args%ATol(1:ntgas)=Atol(1)
    Args%ATol(ntGas+1:nspc)=Atol(2)
    Args%RTol=RtolROW
    Args%RTolpow=RtolRow**pow
    !
    Tstart=Tspan(1)
    Tend=Tspan(2)
    !
    ! Test that Tspan is internally consistent.
    IF (Tstart>=Tend) THEN
      WRITE(*,*) 'start time >= end time'
      CALL FinishMPI()
      STOP 'STOP'
    END IF
    checkdir=Tend-Tstart
    Args%Tdir=SIGN(one,checkdir)
    !
    ! calc first rate 
    CALL Rates(Tstart,y0,Rate)
    !print*, 'debug:: sum(R)=',sum(Rate),sum(y0)
    ALLOCATE(Args%f0(nspc))
    Args%f0=0.0d0
    !
    print*, 'debug :: roargs   ', SUM(Rate), SUM(y0)
    stop
    ! =====================
    ! Neu
    !============
    !CALL FcnRhs(Args%f0,BAT,Rate,y0,ReactionSystem)
    CALL MatVecMult(BAT,Rate,y_emi,Args%f0)
    IF (Args%RTol<=0.0D0) THEN
      WRITE(*,*) 'RTol must be positiv scalar'
      CALL FinishMPI()
      STOP 'STOP'
    END IF
    !
    IF (.NOT.(ALL(Args%ATol>=0.0D0))) THEN
      WRITE(*,*) 'ATols must be positive'
      CALL FinishMPI()
      STOP 'STOP'
    END IF
    ALLOCATE(Args%threshold(nspc))
    Args%threshold=0.0d0
    Args%threshold=Args%ATol/Args%RTol
    ! By default, hmax is 10% of the interval.
    ! less may be better?
    !Args%hmax=0.1D0*Args%hTspan
    Args%hmax=maxStp
    IF (Args%hmax<=0.0D0) THEN
      WRITE(*,*) 'Max step size <= 0'
      CALL FinishMPI()
      STOP 'STOP'
    END IF
  END SUBROUTINE SetRosenbrockArgs
  !
  !
  !==========================================================
  !   Calculates an initial stepsize based on 2. deriv.
  !==========================================================
  SUBROUTINE InitialStepSize(h,hmin,absh,Jac,t,y)
    !------------------------------------------------- 
    ! Input:
    !        - public variables
    !        - Tspan 
    !        - y0  ( initial vector )
    !        - Jacobian matrix
    REAL(RealKind), INTENT(IN) :: t
    REAL(RealKind), INTENT(IN) :: y(:)
    TYPE(CSR_Matrix_T), INTENT(IN) :: Jac
    !-------------------------------------------------
    ! Output:
    !        - initial step size
    REAL(RealKind)  , INTENT(OUT) :: absh
    REAL(RealKind) :: h, hmin 
    !-------------------------------------------------
    !
    ! Temp vars:
    REAL(RealKind) :: tdel, rh
    REAL(RealKind), DIMENSION(nspc) ::  wt, DfDt, Tmp, f1
    REAL(RealKind), DIMENSION(neq) :: Rate
    REAL(RealKind) :: sqrteps=SQRT(eps)
    REAL(RealKind) :: zeros(nspc)
    !
    zeros=ZERO    
    !
    ! hmin is a small number such that t + hmin is clearly different from t in
    ! the working precision, but with this definition, it is 0 if t = 0.
    hmin=minStp
    !
    !---- Compute an initial step size h using yp=y'(t) 
    wt=MAX(ABS(y),Args%threshold)
    rh=(1.25D0*MAXVAL(ABS(Args%f0/wt)))/(Args%RTolpow)
    print*, 'debug:: SUM(wt),rh,sum(f0)', sum(wt) , rh, SUM(Args%f0)
    !
    absh=MIN(Args%hmax,Args%hTspan)
    IF (absh*rh>1.0D0) THEN
      absh=1.0D0/rh
    END IF
    !
    !---- Compute y''(t) and a better initial step size
    h=Args%Tdir*absh
    tdel=(t+Args%Tdir*MIN(sqrteps*MAX(ABS(t),ABS(t+h)),absh))-t
    !
    print*, 'debug:: tdel=     ', tdel
    CALL Rates((t+tdel),y,Rate)
    Output%nRateEvals=Output%nRateEvals+1
    !------NEU
    !
    !CALL FcnRhs(f1,BAT,Rate,y,ReactionSystem)
    CALL MatVecMult(BAT,Rate,y_emi,f1)
     print*, 'debug:: sum(f1)=     ', SUM(f1)
    !
    DfDt=(f1-Args%f0)/tdel
    CALL MatVecMult(Jac,Args%f0,zeros,Tmp)
    !print*, 'ID=',MPI_ID, 'sumDfDT=',SUM(ABS(DfDt)),'tdel=',tdel, 'SUMtmp=',SUM(ABS(tmp)), 'sumf0=',SUM(ABS(Args%f0)), 'sumf1=',SUM(ABS(f1))
    !print*, 'ID=', MPI_ID, 'SumJacvals=',SUM(Jac%Val)
    !print*, 'ID=',MPI_ID,  'SUMtmp=',SUM(ABS(tmp))
    DfDt=DfDt+Tmp
    !
    rh=1.25D0*SQRT(0.5d0*MAXVAL(ABS(DfDt/wt)))/(Args%RTolpow)
    !print*, 'ID=', MPI_ID,'rh=', rh,'rtolpow=',Args%RTolpow, 'sumWT=',SUM(ABS(wt)), 'sumdfdt=',SUM(DfDt)
    !call dropout
    !
    absh=MIN(Args%hmax,Args%hTspan)
    IF (absh*rh>1.0d0) THEN
      absh=1.0D0/rh
    END IF
    absh=MAX(absh,hmin)
    print*, 'debug:: h, absh', h, absh
    stop
  END SUBROUTINE InitialStepSize
  !
  !
  !=======================================================
  !    Subroutine Rosenbrock-Method classic Coef-Matrix
  !=======================================================
  SUBROUTINE ros_classic(ynew,err,errind,y0,t,h,RCo,errVals)
    !--------------------------------------------------------
    ! Input:
    !   - y0............. actual concentrations y 
    !   - t.............. time
    !   - h.............. step size
    !   - RCo............ Rosenbrock method
    REAL(RealKind), INTENT(IN) :: y0(:)
    !REAL(RealKind) :: y0(:)
    REAL(RealKind), INTENT(IN) :: t, h
    TYPE(RosenbrockMethod_T)   :: RCo
    !--------------------------------------------------------
    ! Output:
    !   - ynew........... new concentratinos 
    !   - err............ error calc with embedded formula.
    REAL(RealKind), INTENT(OUT) :: ynew(nspc)
    REAL(RealKind), INTENT(OUT) :: err
    INTEGER,  INTENT(OUT) :: errind(1,1)
    REAL(RealKind), OPTIONAL, INTENT(OUT) :: errVals(nspc)
    !-------------------------------------------------------
    ! Temporary variables:
    REAL(RealKind), DIMENSION(nspc,RCo%nStage) :: k     
    REAL(RealKind), DIMENSION(nspc) :: y,yhat
    REAL(RealKind), DIMENSION(nspc) :: fRhs
    !
    REAL(RealKind) :: Rate(neq)
    REAL(RealKind) :: tt
    !
    REAL(RealKind) :: mpiTEst(SIZE(Miter%Val))
    LOGICAL :: testLOG=.FALSE.
    !
    INTEGER :: iStage, jStage, i, rPtr
    !
    k=ZERO
    fRhs=ZERO
    Rate=ZERO
    !
    !
    y=MAX(ABS(y0(:)),eps)*SIGN(1.0d0,y0(:))
    ynew=y0
    yhat=y0
    !
    ! ---   miter = ( Id - h*gamma*BATransp*Diag(rate)*A*Diag(yVec)^(-1) )
    !
    CALL Rates(t,y,Rate)

    Rate(:)=MAX(ABS(Rate(:)),eps)*SIGN(1.0d0,Rate(:))
 

    TimeJacobianA=MPI_WTIME()
    CALL Miter_Classic(BAT,A,Rate,y,h,RCo%ga,Miter)
    TimeJacobianE=MPI_WTIME()
    TimeJac=TimeJac+(TimeJacobianE-TimeJacobianA)
    !
    ! --- LU - Decomposition ---
    IF (OrderingStrategie==8) THEN
      !
      CALL SetLUvaluesCL(LU_Miter,Miter,LU_Perm)
      timerStart=MPI_WTIME()
      CALL SparseLU(LU_Miter)
      timerEnd=MPI_WTIME()
      TimeFac=TimeFac+(timerEnd-timerStart)
    ELSE
      CALL FactorizeCoefMat(Miter%Val)
    END IF
    !

    !call printsparse(LU_miter,'*')
    WRITE(*,*) '----------------------------'
    WRITE(*,*) 'debug h, t      :: ', h , t
    WRITE(*,*) '      rate(1:3) :: ', rate(1:3), SUM(rate)
    WRITE(*,*) '      conc(1:3) :: ', y(1:3), SUM(y)
    WRITE(*,*) '      sum(Miter):: ', SUM(Miter%val)
    WRITE(*,*) '      sum(LU)   :: ', SUM(LU_Miter%val)
    WRITE(*,*) 
    !
    DO iStage=1,RCo%nStage
      IF (iStage/=1) THEN
        tt=t+RCo%Asum(iStage)*h
        y=y0
        DO jStage=1,iStage
          y=y+RCo%a(iStage,jStage)*k(:,jStage)
        END DO
        !
        ! Update Rates at  (t + SumA*h) , and  (y + A*)k
        CALL Rates(tt,y,Rate)
        !Rate(:)=MAX(ABS(Rate(:)),eps)*SIGN(1.0d0,Rate(:))
      END IF
      !
      ! print*, 'debug :: stage, rates=',iStage,Rate(1:5)
      !
      CALL MatVecMult(BAT,Rate,y_emi,fRhs)           
      fRhs=h*fRhs
      !
      DO jStage=1,iStage-1
        fRhs=fRhs+RCo%C(iStage,jStage)*k(:,jStage)
      END DO
        !print*, 'ID=',MPI_ID,'sumfrhs=',SUM(ABS(frhs))
        !print*, 'debug:: ', 'istage=',iStage,SUM(ABS(fRhs(:)))
      !
      ! ---  Solve linear System  ---
      !
      IF (OrderingStrategie==8) THEN  
        timerStart=MPI_WTIME()
        CALL SolveSparse(LU_Miter,fRhs)
        timerEnd=MPI_WTIME()
        TimeSolve=TimeSolve+(timerEnd-timerStart)
        k(:,iStage)=fRhs(:)
        print*, 'debug:: ', 'istage=',iStage,SUM(ABS(fRhs(:)))
        !print*, 'ID=',MPI_ID, 'istage=',iStage,SUM(ABS(fRhs(:)))
      ELSE
        timerStart=MPI_WTIME()
        !print*, 'vorID =',MPI_ID, 'istage=',iStage,'sum(x)=',SUM(ABS(mumps_par%rhs(:))),'sum(xfrhs=',SUM(ABS(frhs(:))), 'sum matVal=',SUM(ABS(Miter%Val))
        CALL SolveLinAlg(fRhs)
        timerEnd=MPI_WTIME()
        TimeSolve=TimeSolve+(timerEnd-timerStart)
        k(:,iStage)=Mumps_Par%RHS(:)    
        !print*, 'nachID=',MPI_ID, 'istage=',iStage,'sum(x)=',SUM(ABS(mumps_par%rhs(:))),'sum(xfrhs=',SUM(ABS(frhs(:))), 'sum matVal=',SUM(ABS(Miter%Val))
      END IF    
    END DO
      !call dropout
    !  
    DO jStage=1,RCo%nStage
      ynew=ynew+RCo%m(jStage)*k(:,jStage)! new y vector
      yhat=yhat+RCo%me(jStage)*k(:,jStage)! embedded formula for err calc ord-1
    END DO
    !
    ! ---   err calculation  ---
    CALL ERROR(err,errind,ynew,yhat,y0,Args%ATol,Args%RTol,t,errVals)
  END SUBROUTINE ros_classic
  !
  !
  !
  !=======================================================
  !    Subroutine Rosenbrock-Method extended Coef-Matrix
  !=======================================================
  SUBROUTINE ros_extended(ynew,err,errind,y0,t,h,RCo,errVals)
    !--------------------------------------------------------
    ! Input:
    !   - y0............. actual concentrations y 
    !   - t.............. time
    !   - h.............. step size
    !   - RCo............ Rosenbrock method
    REAL(RealKind), INTENT(IN) :: y0(:)
    !REAL(RealKind) :: y0(:)
    REAL(RealKind), INTENT(IN) :: t, h
    TYPE(RosenbrockMethod_T)   :: RCo
    !--------------------------------------------------------
    ! Output:
    !   - YNew........... new concentratinos 
    !   - err............ error calc with embedded formula.
    REAL(RealKind), INTENT(OUT) :: ynew(nspc)
    REAL(RealKind), INTENT(OUT) :: err
    INTEGER,  INTENT(OUT) :: errind(1,1)
    REAL(RealKind), OPTIONAL, INTENT(OUT) :: errVals(nspc)
    !-------------------------------------------------------
    ! Temporary variables:
    REAL(RealKind), DIMENSION(nspc,RCo%nStage) :: k   
    REAL(RealKind), DIMENSION(nspc)   :: y, yhat, hy, kRhs
    REAL(RealKind), DIMENSION(NSactNR) :: bb
    !
    REAL(RealKind) :: Rate(neq)
    REAL(RealKind) :: invRate(neq)
    REAL(RealKind) :: tt
    !
    INTEGER :: iStage, jStage
    INTEGER :: i
    INTEGER :: rPtr
    !
    k=ZERO
    kRhs=ZERO
    Rate=ZERO
    invRate=ZERO
    bb=ZERO
    !
    y=y0
    ynew=y0
    yhat=y0
    !
    ! Update Rates
    !    
    CALL Rates(t,y,Rate)
    !
    IF ( (OrderingStrategie<8) .AND. (MPI_np>1) )THEN
      DO i=1,loc_rateCnt
        rPtr=loc_ratePtr(i)
        invRate(rPtr)=1.0d0/(MAX(ABS(Rate(rPtr)),eps)*SIGN(1.0d0,Rate(rPtr)))
      END DO
      DO i=1,loc_concCnt
        rPtr=loc_concPtr(i)
        hy(rPtr)=MAX(ABS(y(rPtr)),eps)*SIGN(1.0d0,y(rPtr))/h
      END DO
    ELSE
      invRate(:)=1.0d0/(MAX(ABS(Rate(:)),eps)*SIGN(1.0d0,Rate(:)))
      hy(:)=MAX(ABS(y0(:)),eps)*SIGN(1.0d0,y0(:))/h
    END IF

    ! nur für matrix output
    !CALL Miter_Extended(Miter,nspc,neq,invRate,hy)

    !
    !          _                           _ 
    !         |   invRvec    |   g*A_Mat    |
    ! miter = |--------------+--------------|
    !         |_  BAT_Mat    |    Yvec/h   _|
    !
    !
    ! --- LU - Decomposition ---
    IF (OrderingStrategie==8) THEN
      !
      CALL SetLUvaluesEX(LU_Miter,nspc,neq,invRate,hy,ValCopy)  
      timerStart=MPI_WTIME()
      CALL SparseLU(LU_Miter)
      timerEnd=MPI_WTIME()
      TimeFac=TimeFac+(timerEnd-timerStart)
    ELSE
      CALL Miter_Extended(Miter,nspc,neq,invRate,hy)
      CALL FactorizeCoefMat( Miter%Val )
    END IF
    !
    !y=y0   
    !
    DO iStage=1,RCo%nStage
      IF (iStage==1) THEN
        bb(1:neq)=-One
        bb(neq+1:NSactNR)=y_emi(1:nspc)
      ELSE
        tt=t+RCo%Asum(iStage)*h
        y=y0
        DO jStage=1,iStage
          y=y+RCo%a(iStage,jStage)*k(:,jStage)
        END DO
        !
        ! Update Rates at t + a*h
        CALL Rates(tt,y,Rate)
        !Rate(:)=MAX(ABS(Rate(:)),eps)*SIGN(1.0d0,Rate(:))
        kRhs=Zero
        DO jStage=1,iStage-1
          kRhs=kRhs+RCo%C(iStage,jStage)*k(:,jStage)
        END DO
        bb(1:neq)=-invRate(:)*Rate(:)
        bb(neq+1:NSactNR)=kRhs/h+y_emi(1:nspc)
      END IF
      !      
      ! ---  Solve linear System  ---  
      !
      ! Sparse triangular solve LUx=b
      IF (OrderingStrategie==8) THEN  
        timerStart=MPI_WTIME()
        CALL SolveSparse(LU_Miter,bb)
        timerEnd=MPI_WTIME()
        TimeSolve=TimeSolve+(timerEnd-timerStart)
        k(:,iStage)=y0(:)*bb(neq+1:NSactNR)
        !k(:,iStage)=MAX(ABS(y0(:)),eps)*SIGN(1.0d0,y0(:))*bb(neq+1:NSactNR)
        !IF (iStage==1) THEN
        !  DO jstage=1,nspc+neq
        !    IF (t>43100.0d0) WRITE(112,*) jstage, bb(jstage) ! letzte lösung des großes systems speichern
        !  END DO
        !  REWIND(112)
        !END IF
      ELSE
      ! Mumps version of solving LUx=b  
        timerStart=MPI_WTIME()
        CALL SolveLinAlg(bb)
        timerEnd=MPI_WTIME()
        TimeSolve=TimeSolve+(timerEnd-timerStart)
        k(:,iStage)=y0(:)*Mumps_Par%RHS(neq+1:NSactNR)      
        !IF (iStage==1) THEN
        !  DO jstage=1,nspc+neq
        !    IF (t>43100.0d0) WRITE(112,*) jstage, Mumps_Par%rhs(jstage) ! letzte lösung des großes systems speichern
        !  END DO
        !  REWIND(112)
        !END IF
      END IF 
      !print*, '  --  k(1:3,',iStage,')',k(1:3,iStage)
    END DO
    !
    DO jStage=1,RCo%nStage
      ynew=ynew+RCo%m(jStage)*k(:,jStage)   
      yhat=yhat+RCo%me(jStage)*k(:,jStage) 
    END DO
    !
    ! --- error calculation
    !
    !CALL ERROR_MAXNORM(err,errSpc,ynew,yhat,y0,Args%ATol,Args%RTol,t,errVals)
    CALL ERROR(err,errind,ynew,yhat,y0,Args%ATol,Args%RTol,t,errVals)
  END SUBROUTINE ros_extended
  !
  !=======================================================
  !    Subroutine Rosenbrock-Method extended Coef-Matrix
  !=======================================================
  SUBROUTINE ros_temp_extended(ynew,err,errind,y0,t,h,RCo,Temp,errVals)
    !--------------------------------------------------------
    ! Input:
    !   - y0............. actual concentrations y 
    !   - t.............. time
    !   - h.............. step size
    !   - RCo............ Rosenbrock method
    REAL(RealKind), INTENT(IN) :: y0(:)
    REAL(RealKind), INTENT(IN) :: t, h
    REAL(RealKind), INTENT(IN) :: Temp
    TYPE(RosenbrockMethod_T)   :: RCo
    !--------------------------------------------------------
    ! Output:
    !   - YNew........... new concentratinos 
    !   - err............ error calc with embedded formula.
    REAL(RealKind), INTENT(OUT) :: ynew(nspc)
    REAL(RealKind), INTENT(OUT) :: err
    INTEGER,  INTENT(OUT) :: errind(1,1)
    REAL(RealKind), OPTIONAL, INTENT(OUT) :: errVals(nspc)
    !-------------------------------------------------------
    ! Temporary variables:
    REAL(RealKind), DIMENSION(nspc,RCo%nStage) :: k   
    REAL(RealKind), DIMENSION(nspc)   :: y, yhat, hy, kRhs
    REAL(RealKind), DIMENSION(nspc)   :: Umol,DUmoldT,DcDt
    REAL(RealKind), DIMENSION(NSactNR) :: bb
    !
    REAL(RealKind) :: Tarr(7)
    REAL(RealKind) :: Rate(neq)
    REAL(RealKind) :: invRate(neq)
    REAL(RealKind) :: tt
    REAL(RealKind) :: c_v !is the mass average mixture specific heat at constant volume,
    !
    INTEGER :: iStage, jStage
    INTEGER :: i
    !
    k=ZERO
    kRhs=ZERO
    Rate=ZERO
    invRate=ZERO
    bb=ZERO
    !
    y=y0
    ynew=y0
    yhat=y0
    !
    ! Update Rates
    !          _                                  _ 
    !         | invDiagrVec  |   g*A_Mat    | g*~K |
    !         |--------------+--------------+------|
    ! miter = |   BAT_Mat    |  DiagyAec/h  |   0  |
    !         !--------------+--------------+------|
    !         |_     0       |  -U^T*D_c    |   X _|
    !    
    !    X = c_v/h + W^-1 * (g*dU/dt)(dc/dt)
    CALL UpdateTempArray(Tarr,Temp)
   
    CALL Rates(t,y,Rate,Tarr(1))
    invRate(:)=1.0d0/(MAX(ABS(Rate(:)),eps)*SIGN(1.0d0,Rate(:)))
    hy(:)=MAX(ABS(y(:)),eps)*SIGN(1.0d0,y(:))
    !
    CALL SpcInternalEnergy(Umol,Tarr)
    !
    CALL DiffSpcInternalEnergy(DUmoldT,Tarr)
    !
    CALL MassAveMixSpecHeat(c_v,Tarr,y0,DUmoldT)
    !
    CALL DiffConcDt(BAT,Rate,DcDt)
    !
    DUmoldT=DUmoldT*RCo%ga
    c_v=c_v/h
    Umol(:)=-Umol(i)*hy(i)
    hy(:)=hy(:)/h
    !
    ! --- LU - Decomposition ---
    IF (OrderingStrategie==8) THEN
      !
      CALL SetLUvaluesEX(LU_Miter,nspc,neq,invRate,hy,ValCopy)  
      timerStart=MPI_WTIME()
      CALL SparseLU(LU_Miter)
      timerEnd=MPI_WTIME()
      TimeFac=TimeFac+(timerEnd-timerStart)
    ELSE
      CALL Miter_ExtendedTemp(Miter,nspc,neq,invRate,hy,DUmoldT,Umol,c_v,DcDt)
      CALL FactorizeCoefMat(Miter%Val)
    END IF
    !
    y=y0   
    !
    DO iStage=1,RCo%nStage
      IF (iStage==1) THEN
        bb(1:neq)=-One
        bb(neq+1:NSactNR)=y_emi(1:nspc)
      ELSE
        tt=t+RCo%Asum(iStage)*h
        y=y0
        DO jStage=1,iStage
          y=y+RCo%a(iStage,jStage)*k(:,jStage)
        END DO
        !
        ! Update Rates at t + a*h
        CALL Rates(tt,y,Rate)
        kRhs=Zero
        DO jStage=1,iStage-1
          kRhs=kRhs+RCo%C(iStage,jStage)*k(:,jStage)
        END DO
        bb(1:neq)=-invRate(:)*Rate(:)
        bb(neq+1:NSactNR)=kRhs/h+y_emi(1:nspc)
      END IF
      !      
      ! ---  Solve linear System  ---  
      !
      ! Sparse triangular solve LUx=b
      IF (OrderingStrategie==8) THEN  
        timerStart=MPI_WTIME()
        CALL SolveSparse(LU_Miter,bb)
        timerEnd=MPI_WTIME()
        TimeSolve=TimeSolve+(timerEnd-timerStart)
        k(:,iStage)=y0(:)*bb(neq+1:NSactNR)
      ELSE
      ! Mumps version of solving LUx=b  
        timerStart=MPI_WTIME()
        CALL SolveLinAlg(bb)
        timerEnd=MPI_WTIME()
        TimeSolve=TimeSolve+(timerEnd-timerStart)
        k(:,iStage)=y0(:)*Mumps_Par%RHS(neq+1:NSactNR)      
      END IF 
      !print*, '  --  k(1:3,',iStage,')',k(1:3,iStage)
    END DO
    !
    DO jStage=1,RCo%nStage
      ynew=ynew+RCo%m(jStage)*k(:,jStage)   
      yhat=yhat+RCo%me(jStage)*k(:,jStage) 
    END DO
    !
    ! --- error calculation
    !
    CALL ERROR(err,errind,ynew,yhat,y0,Args%ATol,Args%RTol,t,errVals)
  END SUBROUTINE ros_temp_extended
  !
  !
  !=======================================================
  !    Subroutine Rosenbrock-Method classic Coef-Matrix
  !=======================================================
  SUBROUTINE BackwardEuler(ynew,y0,t,h)
    !--------------------------------------------------------
    ! Input:
    !   - y0............. actual concentrations y 
    !   - t.............. time
    !   - h.............. step size
    !   - RCo............ Rosenbrock method
    REAL(RealKind), INTENT(IN) :: y0(:)
    REAL(RealKind), INTENT(IN) :: t, h
    TYPE(RosenbrockMethod_T)   :: RCo
    !--------------------------------------------------------
    ! Output:
    !   - ynew........... new concentratinos 
    !   - err............ error calc with embedded formula.
    REAL(RealKind), INTENT(OUT) :: ynew(nspc)
    !-------------------------------------------------------
    ! Temporary variables:
    REAL(RealKind), DIMENSION(nspc) :: k     
    REAL(RealKind), DIMENSION(nspc) :: y
    REAL(RealKind), DIMENSION(nspc) :: fRhs
    !
    REAL(RealKind) :: Rate(neq)
    !
    !
    k=ZERO
    fRhs=ZERO
    Rate=ZERO
    !
    y=y0
    ynew=y0
    !
    ! ---   miter = ( Id - h*gamma*BATransp*Diag(rate)*A*Diag(yVec)^(-1) )
    !
    CALL Rates(t,y,Rate)
    TimeJacobianA=MPI_WTIME()
    CALL Miter_Classic(BAT,A,Rate,y,h,1.0d0,Miter)
    TimeJacobianE=MPI_WTIME()
    TimeJac=TimeJac+(TimeJacobianE-TimeJacobianA)
    !
    ! --- LU - Decomposition ---
    IF (FactorisationStrategie==1) THEN
      !
      !WRITE(*,*) LU_Perm(1:10)
      !STOP
      CALL SetLUvaluesCL(LU_Miter,Miter,LU_Perm)
      timerStart=MPI_WTIME()
      CALL SparseLU(LU_Miter)
      timerEnd=MPI_WTIME()
      TimeFac=TimeFac+(timerEnd-timerStart)
    ELSE
      !CALL FactorizeCoefMat(Miter%Val)
    END IF
    !
    !
    CALL MatVecMult(BAT,Rate,y_emi,fRhs)           
    fRhs=h*fRhs
    !
    !
    ! ---  Solve linear System  ---
    !
    IF (SolLinSystemStrategie==1) THEN  
      timerStart=MPI_WTIME()
      CALL SolveSparse(LU_Miter,fRhs)
      timerEnd=MPI_WTIME()
      TimeSolve=TimeSolve+(timerEnd-timerStart)
      k(:)=fRhs(:)
    ELSE
      timerStart=MPI_WTIME()
      CALL SolveLinAlg(fRhs)
      timerEnd=MPI_WTIME()
      TimeSolve=TimeSolve+(timerEnd-timerStart)
      k(:)=Mumps_Par%RHS(:)    
    END IF    
    !  
    ynew = ynew + h*k(:)         ! new y vector
    !
  END SUBROUTINE BackwardEuler
END MODULE Rosenbrock_Mod
