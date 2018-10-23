!=========================================================================!
!                                                                         !         
!                   Module for time integration with diagonal             !
!                      implicite Rosenbrock-Wanner-Methods                !         
!                                                                         !         
!=========================================================================!
!
MODULE Integration_Mod
  !
  USE Kind_Mod
  USE Control_Mod
  USE IO_Mod
  USE Reac_Mod
  USE CombustionInput_Mod, ONLY: Density
  USE Sparse_Mod, ONLY: Jac_CC
  USE Rosenbrock_Mod
  USE NetCDF_Mod
  USE Meteo_Mod, ONLY: Temp, cv
  IMPLICIT NONE
  !
  TYPE PI_Param
    REAL(dp) :: Kp
    REAL(dp) :: KI
    REAL(dp) :: ThetaMAX
    REAL(dp) :: rho
  END TYPE PI_Param
 
  INTEGER :: iBar=-1              ! waitbar increment
  TYPE(PI_Param) :: PI_norm, PI_rej
  !
  CONTAINS
  !
  !======================================================================= 
  !===================    Time Integration Routine  ======================
  !======================================================================= 
  SUBROUTINE Integrate(y_iconc, h0, Tspan, ICNTL, RCNTL)
    !--------------------------------------------------------------------
    ! Input:
    !   - y0 ............. Initial vector
    !   - R0 ............. reaction rates at time=Tspan(1)
    !   - Tspan .......... (/ SimulationTimeStart , SimulationTimeEnd /)
  
    REAL(dp), INTENT(INOUT) :: y_iconc(nDIM)
    REAL(dp), INTENT(INOUT) :: h0
    REAL(dp) :: Tspan(2)
    INTEGER  :: ICNTL(:)
    REAL(dp) :: RCNTL(:)
    !-------------------------------------------------------------------
    ! Output:
    !   - Output.......... struct Out (above) contains output statistics
    !   - Unit............ .simul data 
    !-------------------------------------------------------------------
    ! Temporary variables:
    !
    REAL(dp) :: Y0(nDIM), Y(nDIM)       ! old, current y vector
    !
    REAL(dp) :: t             ! current time
    REAL(dp) :: timepart
    REAL(dp) :: time_int = 0.0d0
    REAL(dp) :: Atol(2)

    REAL(dp) :: h, tnew, tmp, hOld
    REAL(dp) :: error, errorOld
    REAL(dp) :: tmp_tb
    REAL(dp) :: actLWC, zen, Temperature
    INTEGER  :: ierr(1,1)
   ! 
    INTEGER :: i, k

    INTEGER, PARAMETER :: maxnsteps = 50000
    CHARACTER(20) :: tmethod=''
    
    LOGICAL :: done=.FALSE.
    LOGICAL :: failed
    !
    !-----------------------------
    ! LSODE Parameter
    INTEGER :: ITOL , ITASK , ISTATE, IOPT, MF, LIW, LRW
    INTEGER,  ALLOCATABLE :: IWORK(:)
    REAL(dp), ALLOCATABLE :: RWORK(:)
    !
    REAL(dp) :: RTOL1 , ATOL1(nDIM)

    
    time_int = MPI_WTIME()

    Y0 = y_iconc
    Y  = y_iconc

    IF ( Combustion ) THEN
      Y0 = MAX( ABS(Y0) , eps ) * SIGN( ONE , Y0 )    ! |y| >= eps
      Y  = MAX( ABS(Y)  , eps ) * SIGN( ONE , Y  )    ! |y| >= eps
    END IF

    t = Tspan(1)
    tnew = Tspan(1)
    h = h0

    ! this is for the waitbar
    tmp_tb = (tEnd-tBegin) * 0.01_dp

    Atol = [RCNTL(1), RCNTL(2)]
    
    IF (INDEX(ODEsolver,'LSODE')>0) THEN
      tmethod = ADJUSTL(ODEsolver)
    ELSE
      tmethod = ADJUSTL(ODEsolver(INDEX(ODEsolver,'/')+1 : INDEX(ODEsolver,'.')-1))
    END IF
    
    SELECT CASE (TRIM(tmethod))

      CASE ('bwEuler')

        h = maxStp

        BACKWARD_EULER: DO
          
          !-- Stretch the step if within 10% of tfinal-t.
          IF ( 1.05_dp*h >= Tspan(2)-t ) THEN
            h    = Tspan(2) - t
            done = .TRUE.
          END IF

          CALL Rosenbrock(  Y             & ! new concentration
          &               , error         & ! error value
          &               , ierr          & ! max error component
          &               , Y0            & ! current concentration 
          &               , t             & ! current time
          &               , h             & ! stepsize
          &               , Euler=.TRUE. )  ! new concentration 

          tnew = t + h
          IF (done) THEN
            tnew = Tspan(2)  ! Hit end point exactly.
            h    = tnew - t  ! Purify h.
          END IF

          Out%nRateEvals = Out%nRateEvals + 1
          Out%nSolves    = Out%nSolves + 1
          Out%nsteps     = Out%nsteps + 1

          IF ( NetCdfPrint ) THEN 
            TimeNetCDFA = MPI_WTIME()
            IF ( done .OR. StpNetCDF < ZERO ) THEN 
              IF ( Combustion ) Temperature = Y(nDIM)
              CALL SetOutputNCDF( NetCDF, tnew , h , Y )
              CALL StepNetCDF( NetCDF )
            END IF
            TimeNetCDF  = TimeNetCDF + (MPI_WTIME() - TimeNetCDFA)
          END IF

          !-- Call progress bar.
          IF ( WaitBar .AND. t-Tspan(1) >= iBar*tmp_tb ) THEN
            iBar = iBar + 1         ! iBar runs from 0 to 100
            CALL Progress(iBar)
          END IF

          !-- Termination condition for the main loop.
          IF ( done ) EXIT
            
          !-- Advance the integration one step.
          t   = tnew
          Y0  = Y

        END DO BACKWARD_EULER
    
      CASE('LSODE')

        time_int=MPI_WTIME()

        LIW    = 20 + nDIM
        LRW    = 22 + 9*nDIM + nDIM**2
        ITOL   = 2              ! 2 if atol is an array
        RTOL1  = RtolRow        ! relative tolerance parameter
        ATOL1(1:nspc) = Atol(1)  ! atol dimension nspc+1
        ATOL1(nDim)   = Atol(2)  ! temperature tolerance

        ITASK  = 1             ! 1 for normal computation of output values of y at t = TOUT
        ISTATE = 1             ! This is the first call for a problem
        IOPT   = 1             ! optional inputs are used (see Long Desciption Part 1.)
        MF     = 22            ! Stiff method, internally generated full Jacobian.
        
        ALLOCATE(IWORK(LIW))
        ALLOCATE(RWORK(LRW))
        IWORK  = 0
        RWORK  = ZERO

        !RWORK(5) = 1.0d-6    ! first step size, if commented -> calculate h0
        RWORK(6) = maxStp
        RWORK(7) = minStp

        t        = Tspan(1)
        timepart = (Tspan(2)-Tspan(1)) / REAL(nOutP - 1, dp)
        tnew     = t + timepart

        MAIN_LOOP_LSODE: DO k = 1 , nOutP-1
          
          IWORK(6) = maxnsteps
           
          CALL DLSODE ( FRhs  , nDIM  , Y     , t     , tnew           &
          &           , ITOL  , RTOL1 , ATOL1 , ITASK , ISTATE , IOPT  &
          &           , RWORK , LRW   , IWORK , LIW   , dummy  , MF    )
  
  
          Out%nsteps     = Out%nsteps     + IWORK(11)
          Out%nRateEvals = Out%nRateEvals + IWORK(12)
          Out%npds       = Out%npds       + IWORK(13)
          Out%ndecomps   = Out%ndecomps   + IWORK(21)

          ! save data
          IF ( NetCdfPrint ) THEN 
            TimeNetCDFA = MPI_WTIME()
            ! Time to save a step?
            IF ( done .OR. StpNetCDF < ZERO ) THEN 
              IF ( Combustion ) Temperature = Y(nDIM)
              CALL SetOutputNCDF( NetCDF, tnew , h , Y )
              CALL StepNetCDF( NetCDF )
            END IF
            TimeNetCDF  = TimeNetCDF + (MPI_WTIME() - TimeNetCDFA)
          END IF


          !-- Call progress bar.
          IF ( WaitBar .AND. t-Tspan(1) >= iBar*tmp_tb ) THEN
            iBar = iBar + 1         ! iBar runs from 0 to 100
            CALL Progress(iBar)
          END IF

          t    = tnew
          tnew = tnew + timepart
        
        END DO MAIN_LOOP_LSODE

        CALL Progress(100) ! last * needs an extra call


      CASE('LSODES')

        time_int=MPI_WTIME()

        LIW    = 30
        LRW    = 20 + 21 * nDIM + 4*(BAT%nnz+A%nnz+2*nDIM)
        ITOL   = 2              ! 2 if atol is an array
        RTOL1  = RtolRow        ! relative tolerance parameter
        IF (Combustion) THEN
          ATOL1(1:nspc) = Atol(1)  ! atol dimension nspc+1
          ATOL1(nDim)   = Atol(2)  ! temperature tolerance
        ELSE
          ATOL1(1:ns_Gas)  = Atol(1)  ! atol dimension nspc+1
          ATOL1(ns_gas+1:) = Atol(2)  ! temperature tolerance
        END IF

        ITASK  = 1             ! 1 for normal computation of output values of y at t = TOUT
        ISTATE = 1             ! This is the first call for a problem
        IOPT   = 1             ! optional inputs are used (see Long Desciption Part 1.)
        MF     = 222            ! Stiff method, internally generated full Jacobian.
        
        ALLOCATE(IWORK(LIW))
        ALLOCATE(RWORK(LRW))
        IWORK  = 0
        RWORK  = ZERO

        !RWORK(5) = 1.0d-6    ! first step size, if commented -> calculate h0
        RWORK(6) = maxStp
        RWORK(7) = minStp

        t        = Tspan(1)
        timepart = (Tspan(2)-Tspan(1)) / REAL(nOutP - 1, dp)
        tnew     = t + timepart

        MAIN_LOOP_LSODES: DO k = 1 , nOutP-1
          
          IWORK(6) = maxnsteps
           
          CALL DLSODES( FRhs  , nDIM  , Y     , t     , tnew           &
          &           , ITOL  , RTOL1 , ATOL1 , ITASK , ISTATE , IOPT  &
          &           , RWORK , LRW   , IWORK , LIW   , dummy  , MF    )
  
  
          Out%nsteps     = Out%nsteps     + IWORK(11)
          Out%nRateEvals = Out%nRateEvals + IWORK(12)
          Out%npds       = Out%npds       + IWORK(13)
          Out%ndecomps   = Out%ndecomps   + IWORK(21)

          ! save data
          IF ( NetCdfPrint ) THEN 
            TimeNetCDFA = MPI_WTIME()
            IF ( done .OR. StpNetCDF < ZERO ) THEN 
              IF ( Combustion ) Temperature = Y(nDIM)
              CALL SetOutputNCDF( NetCDF, tnew , h , Y )
              CALL StepNetCDF( NetCDF )
            END IF
            TimeNetCDF  = TimeNetCDF + (MPI_WTIME() - TimeNetCDFA)
          END IF


          !-- Call progress bar.
          IF ( WaitBar .AND. t-Tspan(1) >= iBar*tmp_tb ) THEN
            iBar = iBar + 1         ! iBar runs from 0 to 100
            CALL Progress(iBar)
          END IF

          t    = tnew
          tnew = tnew + timepart
        
        END DO MAIN_LOOP_LSODES

        CALL Progress(100) ! last * needs an extra call

      

      CASE DEFAULT

        !IF (MPI_master) WRITE(*,'(10X,A)',ADVANCE='NO') 'Start Integration.............      '
       
            
        MAIN_LOOP_ROSENBROCK_METHOD: DO 
          
          h = MIN( maxStp, MAX( minStp , h ) )
          
          !-- Stretch the step if within 5% of tfinal-t.
          IF ( 1.05_dp * h >= Tspan(2) - t ) THEN
            h = ABS(Tspan(2) - t)
            done = .TRUE.
          END IF

          DO    ! Evaluate the formula.
            
            !-- LOOP FOR ADVANCING ONE STEP.
            failed = .FALSE.      ! no failed attempts
            
            ! Rosenbrock Timestep 
            CALL Rosenbrock(  Y             & ! new concentration
            &               , error         & ! error value
            &               , ierr          & ! max error component
            &               , Y0            & ! current concentration 
            &               , t             & ! current time
            &               , h             & ! stepsize
            &               , Euler=.FALSE. ) ! new concentration 

            ! test Willi
            !WHERE (Y<1.0e-4_dp)
            !  Y = 1.0e-20_dp
            !END WHERE
            
            tnew  = t + h
            
            IF (done) THEN
              tnew  = Tspan(2)    ! Hit end point exactly.
              h     = tnew-t      ! Purify h.
            END IF
            Out%ndecomps   = Out%ndecomps   + 1
            Out%nRateEvals = Out%nRateEvals + ROS%nStage
            Out%nSolves    = Out%nSolves    + ROS%nStage
            
            
            failed = (error > ONE)
            
            IF (failed) THEN      ! failed step
              ! Accept the solution only if the weighted error is no more than the
              ! tolerance one.  Estimate an h that will yield an error of rtol on
              ! the next step or the next try at taking this step, as the case may be,
              ! and use 0.8 of this value to avoid failures.
              
              Out%nfailed  = Out%nfailed+1
              IF ( h <= minStp ) THEN
                STOP '....Integration_Mod '
              END IF
              
              h    = MAX( minstp , h * MAX( rTEN , 0.8_dp * error**(-ROS%pow) ) )
              done = .FALSE.
            ELSE                  ! successful step
              EXIT
            END IF
          END DO
          Out%nsteps = Out%nsteps + 1
         
          ! --- save to NetCDF file
          IF ( NetCdfPrint ) THEN
            IF ( done .OR. StpNetCDF < ZERO ) THEN 
              TimeNetCDFA = MPI_WTIME()
              CALL SetOutputNCDF( NetCDF, tnew , h , Y )
              CALL StepNetCDF( NetCDF )
              !WRITE(999,'(*(Es14.6))') tnew, ( Y(Diag_Index(i)), i=1,SIZE(Diag_Index) ) 
              TimeNetCDF  = TimeNetCDF + (MPI_WTIME() - TimeNetCDFA)
            END IF
          END IF 

          !-- Call progress bar.
          IF ( WaitBar .AND. t-tBegin > iBar*tmp_tb ) CALL Progress(iBar)
          
          !-- If there were no failures compute a new h.
          tmp = 1.25_dp * error**ROS%pow
          IF ( TWO * tmp > ONE ) THEN
            h = h / tmp
          ELSE
            h = h * 4
          END IF
          
          !-- Advance the integration one step.
          t  = tnew
          Y0 = Y
          
          !-- for PI stepsize control
          errorOld  = error
          hOld      = h

          !-- Termination condition for the main loop.
          IF ( done ) EXIT  MAIN_LOOP_ROSENBROCK_METHOD

        END DO MAIN_LOOP_ROSENBROCK_METHOD  ! MAIN LOOP
        
        !WRITE(*,*) ' No other methods implemented jet. Use Rosenbrock, LSODE[S] or backward Euler'

    END SELECT
    ! save values for next initial vector
    y_iconc = Y
    h0 = h

    TimeIntegration = TimeIntegration + MPI_WTIME() - time_int


    ! --- stop timer and print output statistics
    Timer_Finish = MPI_WTIME() - Timer_Start + Time_Read

  END SUBROUTINE Integrate
  
 
  ! The Progressbar
  SUBROUTINE Progress(j)
    !
    INTEGER(4)    :: j,k
    CHARACTER(69) :: bar="          Start Integration...........    ???% |                    |"

    iBar = iBar + 1         ! iBar runs from 0 to 100
    !
    IF (MPI_ID==0) THEN
      WRITE(bar(43:45),'(I3)') j
      !
      DO k=1,j/5
        bar(48+k:48+k)="*"
      END DO
      ! print the progress bar.
      WRITE(*,'(A1,A69,$)') char(13), bar
    END IF
  END SUBROUTINE Progress
  !
  !
  SUBROUTINE PI_StepsizeControl(hnew,Tol,er,erOld,h,hOld,PI_Par,RCo)
    REAL(dp) :: hnew
    !
    !
    REAL(dp) :: Tol       ! rel tol
    REAL(dp) :: er        ! lokal error step n 
    REAL(dp) :: erOld     ! lokal error step n-1
    REAL(dp) :: h         ! stepsize n
    REAL(dp) :: hOld      ! stepsize n-1
    TYPE(PI_Param) :: PI_Par    ! PI control parameter
    TYPE(RosenbrockMethod_T)  :: RCo
    !
    REAL(dp) :: htmp
    !
    !
    htmp = ( Tol/er )**(0.7d0*ROS%pow)* ( erOld/Tol )**(0.4d0*ROS%pow) * h
    !htmp = ( Tol/er )**PI_Par%KI * ( erOld/Tol )**PI_Par%Kp * h
    !
    ! limitation term
    IF (htmp>PI_Par%ThetaMAX*h) THEN
      hnew = PI_Par%ThetaMAX*h
    ELSE
      hnew = htmp
    END IF
  END SUBROUTINE PI_StepsizeControl

  SUBROUTINE FRhs(NEQ1, T, Y, YDOT)
    USE Rates_Mod
    INTEGER        :: NEQ1
    REAL(dp) :: T , Y(NEQ1) , YDOT(NEQ1)
    !
    REAL(dp) :: dCdt(nspc), Emiss(nspc)
    REAL(dp) :: Rate(neq) , DRate(neq)
    REAL(dp) :: Tarr(10)
    REAL(dp) :: U(nspc) , dUdT(nspc)
    REAL(dp) :: c_v
    
    IF (Combustion) THEN         ! OUT:   IN:
      ! MASS CONSERVATION
      CALL ReactionRates( T , Y , Rate , DRate)
      CALL UpdateEmission(Emiss,T)
      dCdt = BAT * Rate + Emiss   ! [mol/cm3/s]
    
      YDOT(1:nspc) = dCdt

      ! ENERGY CONSERVATION
      Tarr = UpdateTempArray ( Y(NEQ1) )
      CALL InternalEnergy    (  U   , Tarr )
      CALL DiffInternalEnergy( dUdT , Tarr) 
      CALL MassAveMixSpecHeat( c_v  , dUdT               &
                             &      , MoleConc=Y(1:nspc) &
                             &      , rho=rho          ) 
    
      YDOT(NEQ1) = - SUM( U * dCdt ) * rRho / c_v   ! rRho=mega/rho 
    ELSE
      ! MASS CONSERVATION
      CALL ReactionRates( T , Y , Rate )
      CALL UpdateEmission(Emiss,T)
      dCdt = BAT * Rate + Emiss   ! [mol/cm3/s]
      
      YDOT(1:nspc) = dCdt
    END IF

  END SUBROUTINE FRhs

     
  SUBROUTINE Jacobian_LSODES (NEQ1, T, Y, J, IA, JA, PDJ)
    INTEGER  ::  NEQ1, J 
    REAL(dp) :: T, Y(NEQ1), PDJ(NEQ1)
    INTEGER  :: IA(:), JA(:)


 !1   PDJ(1) = -RK1
     !PDJ(2) = RK1
     !RETURN
 !2   PDJ(2) = -RK3*Y(3) - RK15*Y(12) - RK2
     !PDJ(3) = RK2 - RK3*Y(3)
     !PDJ(4) = RK3*Y(3)
     !PDJ(5) = RK15*Y(12)
     !PDJ(12) = -RK15*Y(12)
     !RETURN
 !3   PDJ(2) = -RK3*Y(2)
     !PDJ(3) = -RK5 - RK3*Y(2) - RK7*Y(10)
     !PDJ(4) = RK3*Y(2)
     !PDJ(6) = RK7*Y(10)
     !PDJ(10) = RK5 - RK7*Y(10)
     !RETURN
 !4   PDJ(2) = RK11*RK14
     !PDJ(3) = RK11*RK14
     !PDJ(4) = -RK11*RK14 - RK4
     !PDJ(9) = RK4
     !RETURN
 !5   PDJ(2) = RK19*RK14
     !PDJ(5) = -RK19*RK14 - RK16
     !PDJ(9) = RK16
     !PDJ(12) = RK19*RK14
     !RETURN
 !6   PDJ(3) = RK12*RK14
     !PDJ(6) = -RK12*RK14 - RK8
     !PDJ(9) = RK8
     !PDJ(10) = RK12*RK14
     !RETURN
 !7   PDJ(7) = -RK20*RK14 - RK18
     !PDJ(9) = RK18
     !PDJ(10) = RK20*RK14
     !PDJ(12) = RK20*RK14
     !RETURN
 !8   PDJ(8) = -RK13*RK14 - RK10
     !PDJ(10) = RK13*RK14
     !PDJ(11) = RK10
 !9   RETURN
 !10  PDJ(3) = -RK7*Y(3)
     !PDJ(6) = RK7*Y(3)
     !PDJ(7) = RK17*Y(12)
     !PDJ(8) = RK9
     !PDJ(10) = -RK7*Y(3) - RK17*Y(12) - RK6 - RK9
     !PDJ(12) = RK6 - RK17*Y(12)
 !11  RETURN
 !12  PDJ(2) = -RK15*Y(2)
     !PDJ(5) = RK15*Y(2)
     !PDJ(7) = RK17*Y(10)
     !PDJ(10) = -RK17*Y(10)
     !PDJ(12) = -RK15*Y(2) - RK17*Y(10)
     !RETURN
  END SUBROUTINE Jacobian_LSODES

  ! dummy routine for ode solver LSODE
  SUBROUTINE dummy()
    IMPLICIT NONE
  END SUBROUTINE dummy

END MODULE Integration_Mod
