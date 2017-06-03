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
  USE mo_MPI
  USE mo_control
  USE mo_IO
  USE mo_reac
  USE mo_ckinput, ONLY: Density
  USE Sparse_Mod, ONLY: Jac_CC
  !USE Sparse2_Mod
  !USE Chemsys_Mod
  USE Rosenbrock_Mod
  !USE Rates_Mod
  USE Factorisation_Mod
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
  !
  TYPE(PI_Param) :: PI_norm, PI_rej
  !
  CONTAINS
  !
  !======================================================================= 
  !===================    Time Integration Routine  ======================
  !======================================================================= 
  SUBROUTINE Integrate(y_iconc, R0, Tspan, Atol, RtolRow, method)
    !--------------------------------------------------------------------
    ! Input:
    !   - y0 ............. Initial vector
    !   - R0 ............. reaction rates at time=Tspan(1)
    !   - Tspan .......... (/ SimulationTimeStart , SimulationTimeEnd /)
    !   - Atol ........... abs. tolerance for gas spc
    !   - RtolRow ........ rel. tolerance for Rosenbrock-Wanner-Method
    !   - method.......... Rosenbrock-Wanner-Method
    !   - PrintSpc ....... print spc PrintSpc(1:3)
    REAL(dp) :: y_iconc(nspc)
    REAL(dp) :: R0(neq)
    REAL(dp) :: Tspan(2)
    REAL(dp) :: Atol(2)
    REAL(dp) :: RtolROW
    CHARACTER(*) :: method
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
    REAL(dp) :: StartTimer

    REAL(dp) :: h, hmin, absh, tnew, tmp, hOld
    REAL(dp) :: error, errorOld
    REAL(dp) :: tmp_tb
    REAL(dp) :: actLWC, zen, wetRad
    INTEGER  :: errind(1,1)
    REAL(dp) :: ErrVals(nspc)
    ! 
    INTEGER :: iBar=-1              ! waitbar increment
    INTEGER :: i, k

    INTEGER, PARAMETER :: maxnsteps = 50000
    
    LOGICAL :: done=.FALSE.
    LOGICAL :: failed
    !
    CHARACTER(10) :: swMethod
    !
    !-----------------------------
    ! LSODE Parameter
    INTEGER :: ITOL , ITASK , ISTATE, IOPT, MF, LIW, LRW
    INTEGER, ALLOCATABLE :: IWORK(:)
    REAL(dp), ALLOCATABLE :: RWORK(:)
    !
    REAL(dp) :: RTOL1 , ATOL1(nDIM)

    
    Y0(1:nspc)  = y_iconc(:)
    Y(1:nspc)   = y_iconc(:)
   
    IF ( TempEq ) THEN
      !--- initial temperature
      Y0(nDIM) = Temperature0     ! = 750 [K] aus speedchem debug
      Y(nDIM)  = Temperature0     ! = 750 [K]
    END IF

    ALLOCATE(loc_RatePtr(neq))
    loc_RatePtr = [(i , i=1,neq)]
    loc_rateCnt = neq
    
    ! this is for the waitbar
    tmp_tb = (Tspan(2)-Tspan(1)) * 0.01_dp
    
    SELECT CASE (TRIM(method(1:7)))

      CASE('METHODS')
    
        !---- Calculate a first stepsize based on 2nd deriv.
        CALL InitialStepSize( h , hmin , absh , Jac_CC , R0   &
        &                   , Tspan(1) , Y(1:nspc) , RCo%pow  )
     
        !
        !===============================================================================
        !=================================THE MAIN LOOP=================================
        !===============================================================================
        TimeIntegrationA = MPI_WTIME()
        IF (MPI_ID==0) WRITE(*,*) '  Start Integration............. '; WRITE(*,*) ' '
       
        t     = Tspan(1)
        tnew  = Tspan(1)
            
        MAIN_LOOP_ROSENBROCK_METHOD: DO 
          !
          absh  = MIN( maxStp, MAX( minStp , absh ) )
          h     = absh
          !
          !-- Stretch the step if within 5% of tfinal-t.
          IF ( 1.05_dp * absh >= Tspan(2) - t ) THEN
            h     = Tspan(2) - t
            absh  = ABS(h)
            done  = .TRUE.
          END IF
          DO                                ! Evaluate the formula.
            !
            !-- LOOP FOR ADVANCING ONE STEP.
            failed  = .FALSE.                   ! no failed attempts
            !
            ! Rosenbrock Timestep 
            CALL Rosenbrock(  Y0            &       ! current concentration
            &               , t             &       ! current time
            &               , h             &       ! stepsize
            &               , RCo           &       ! Rosenbrock parameter
            &               , error         &       ! error value
            &               , errind        &       ! max error component
            &               , Y             &       ! new concentration 
            &               , Euler=.FALSE. )       ! new concentration 
            !
            tnew  = t + h
            !
            IF (done) THEN
              tnew  = Tspan(2)         ! Hit end point exactly.
              h     = tnew-t                        ! Purify h.
            END IF
            Output%ndecomps   = Output%ndecomps   + 1
            Output%nRateEvals = Output%nRateEvals + RCo%nStage
            Output%nSolves    = Output%nSolves    + RCo%nStage
            !
            !
            failed = (error > ONE)
            !
            IF (failed) THEN               !failed step
              ! Accept the solution only if the weighted error is no more than the
              ! tolerance rtol.  Estimate an h that will yield an error of rtol on
              ! the next step or the next try at taking this step, as the case may be,
              ! and use 0.8 of this value to avoid failures.
              !
              Output%nfailed  = Output%nfailed+1
              IF ( absh <= hmin ) THEN
                WRITE(*,*) ' Stepsize to small: ', h
                CALL FinishMPI()
                STOP '....Integration_Mod '
              END IF
              !
              absh  = MAX( hmin , absh * MAX( rTEN , 0.8_dp * error**(-RCo%pow) ) )
              h     = absh
              done  = .FALSE.
            ELSE                            !succ. step
              EXIT
            END IF
          END DO
          Output%nsteps = Output%nsteps + 1
          !
          !
          IF ( NetCdfPrint ) THEN 
            TimeNetCDFA = MPI_WTIME()
            IF ( (t - Tspan(1) >= StpNetCDF*iStpNetCDF) )  THEN
              iStpNetCDF  = iStpNetCDF + 1
              zen = Zenith(t)
              IF ( ntAqua > 0 ) THEN
                actLWC  = pseudoLWC(t)
                wetRad  = (Pi34*actLWC/SPEK(1)%Number)**(rTHREE)*0.1_dp
              END IF
              ! save data in NetCDF File
              CALL SetOutputNCDF(  Y,    yNcdf, t ,  actLWC)
              !
              CALL StepNetCDF   ( t                                &
              &                 , yNcdf                            &
              &                 , itime_NetCDF                     &
              &                 , (/ actLWC                        &
              &                    , h                             &
              &                    , SUM(Y(1:ntGas))               &
              &                    , SUM(Y(ntGas+1:ntGas+ntAqua))  &
              &                    , wetRad                        &
              &                    , zen                        /) &
              &                 , errind                           &
              &                 , SUM(y(iDiag_Schwefel))           &
              &                 , error                            )
              !
              TimeNetCDF  = TimeNetCDF + (MPI_WTIME() - TimeNetCDFA)
            END IF
          END IF
          !
          !-- Call progress bar.
          IF ( Bar .AND. t-Tspan(1) >= iBar*tmp_tb ) THEN
            iBar = iBar + 1         ! iBar runs from 0 to 100
            CALL Progress(iBar)
          END IF
          !
          !-- Termination condition for the main loop.
          IF ( done ) EXIT
          !
          !-- If there were no failures compute a new h.
          tmp = 1.25_dp * error**RCo%pow
          IF ( TWO * tmp > ONE ) THEN
            absh  = absh / tmp
          ELSE
            absh  = absh * TWO
          END IF
          !
          !-- Advance the integration one step.
          t   = tnew
          Y0  = Y
          !
          !-- for PI stepsize control
          errorOld  = error
          hOld      = h

        END DO MAIN_LOOP_ROSENBROCK_METHOD  ! MAIN LOOP
    
      CASE('LSODE')

        TimeIntegrationA=MPI_WTIME()

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


        IF (MPI_ID==0) WRITE(*,*) '  Start Integration............. '; WRITE(*,*) ' '

        t        = Tspan(1)
        timepart = (Tspan(2)-Tspan(1)) / REAL(nOutP - 1, dp)
        tnew     = t + timepart

        MAIN_LOOP_LSODE: DO k = 1 , nOutP-1
          
          IWORK(6) = maxnsteps
           
          CALL DLSODE ( FRhs  , nDIM  , Y     , t        , tnew            &
          &           , ITOL  , RTOL1 , ATOL1 , ITASK    , ISTATE  , IOPT  &
          &           , RWORK , LRW   , IWORK , LIW      , dummy   , MF    )
  
  
          Output%nsteps     = Output%nsteps     + IWORK(11)
          Output%nRateEvals = Output%nRateEvals + IWORK(12)
          Output%npds       = Output%npds       + IWORK(13)
          Output%ndecomps   = Output%ndecomps   + IWORK(21)

          ! save data
          IF ( NetCdfPrint ) THEN 
            TimeNetCDFA = MPI_WTIME()
            iStpNetCDF  = iStpNetCDF + 1
            zen = Zenith(t)
            IF ( ntAqua > 0 ) THEN
              actLWC  = pseudoLWC(t)
              wetRad  = (Pi34*actLWC/SPEK(1)%Number)**(rTHREE)*0.1_dp
            END IF

            ! save data in NetCDF File
            CALL SetOutputNCDF(  Y , yNcdf, t , actLWC)
          
            CALL StepNetCDF (   t                      &
                            & , yNcdf                  &
                            & , itime_NetCDF           &
                            & , (/  actLWC                        &  ! LWC
                            &     , RWORK(11)                     &  ! Stepsize
                            &     , SUM(Y(1:ntGas))               &    ! Sum gas conc
                            &     , SUM(Y(ntGas+1:ntGas+ntAqua))  &    ! Sum aqua conc
                            &     , wetRad                        &    ! wet radius
                            &     , zen               /)          &  ! zenith
                            & , errind                 &  ! max error index = 0
                            & , ZERO                   &  ! Sum sulfuric conc
                            & , error                 )   ! error 
            
            TimeNetCDF  = TimeNetCDF + (MPI_WTIME() - TimeNetCDFA)
          END IF


          !-- Call progress bar.
          IF ( Bar .AND. t-Tspan(1) >= iBar*tmp_tb ) THEN
            iBar = iBar + 1         ! iBar runs from 0 to 100
            CALL Progress(iBar)
          END IF

          t    = tnew
          tnew = tnew + timepart
        
        END DO MAIN_LOOP_LSODE

        CALL Progress(100) ! last * needs an extra call


      CASE('LSODES')

        TimeIntegrationA=MPI_WTIME()

        LIW    = 30
        LRW    = 20 + 21 * nDIM + 4*(BAT%nnz+A%nnz+2*nDIM)
        ITOL   = 2              ! 2 if atol is an array
        RTOL1  = RtolRow        ! relative tolerance parameter
        ATOL1(1:nspc) = Atol(1)  ! atol dimension nspc+1
        ATOL1(nDim)   = Atol(2)  ! temperature tolerance

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


        IF (MPI_ID==0) WRITE(*,*) '  Start Integration............. '; WRITE(*,*) ' '

        t        = Tspan(1)
        timepart = (Tspan(2)-Tspan(1)) / REAL(nOutP - 1, dp)
        tnew     = t + timepart

        MAIN_LOOP_LSODES: DO k = 1 , nOutP-1
          
          IWORK(6) = maxnsteps
           
          CALL DLSODES( FRhs  , nDIM  , Y     , t        , tnew            &
          &           , ITOL  , RTOL1 , ATOL1 , ITASK    , ISTATE  , IOPT  &
          &           , RWORK , LRW   , IWORK , LIW      , dummy   , MF    )
  
  
          Output%nsteps     = Output%nsteps     + IWORK(11)
          Output%nRateEvals = Output%nRateEvals + IWORK(12)
          Output%npds       = Output%npds       + IWORK(13)
          Output%ndecomps   = Output%ndecomps   + IWORK(21)

          ! save data
          IF ( NetCdfPrint ) THEN 
            TimeNetCDFA = MPI_WTIME()
            iStpNetCDF  = iStpNetCDF + 1
            zen = Zenith(t)
            IF ( ntAqua > 0 ) THEN
              actLWC  = pseudoLWC(t)
              wetRad  = (Pi34*actLWC/SPEK(1)%Number)**(rTHREE)*0.1_dp
            END IF

            ! save data in NetCDF File
            CALL SetOutputNCDF(  Y , yNcdf, t , actLWC)
          
            CALL StepNetCDF (   t                      &
                            & , yNcdf                  &
                            & , itime_NetCDF           &
                            & , (/  actLWC                        &  ! LWC
                            &     , RWORK(11)                     &  ! Stepsize
                            &     , SUM(Y(1:ntGas))               &    ! Sum gas conc
                            &     , SUM(Y(ntGas+1:ntGas+ntAqua))  &    ! Sum aqua conc
                            &     , wetRad                        &    ! wet radius
                            &     , zen               /)          &  ! zenith
                            & , errind                 &  ! max error index = 0
                            & , ZERO                   &  ! Sum sulfuric conc
                            & , error                 )   ! error 
            
            TimeNetCDF  = TimeNetCDF + (MPI_WTIME() - TimeNetCDFA)
          END IF


          !-- Call progress bar.
          IF ( Bar .AND. t-Tspan(1) >= iBar*tmp_tb ) THEN
            iBar = iBar + 1         ! iBar runs from 0 to 100
            CALL Progress(iBar)
          END IF

          t    = tnew
          tnew = tnew + timepart
        
        END DO MAIN_LOOP_LSODES

        CALL Progress(100) ! last * needs an extra call


      CASE ('bwEuler')

        h = maxStp
        
        BACKWARD_EULER: DO
          
          !-- Stretch the step if within 10% of tfinal-t.
          IF ( 1.05_dp*h >= Tspan(2)-t ) THEN
            h    = Tspan(2) - t
            done = .TRUE.
          END IF

          CALL Rosenbrock(  Y0            &       ! current concentration
          &               , t             &       ! current time
          &               , h             &       ! stepsize
          &               , RCo           &       ! Rosenbrock parameter
          &               , error         &       ! error value
          &               , errind        &       ! max error component
          &               , Y             &       ! new concentration 
          &               , Euler=.TRUE.  )       ! specifies backward Euler method

          tnew = t + h
          IF (done) THEN
            tnew = Tspan(2)         ! Hit end point exactly.
            h    = tnew - t         ! Purify h.
          END IF

          Output%nRateEvals = Output%nRateEvals + 1
          Output%nSolves    = Output%nSolves + 1
          Output%nsteps     = Output%nsteps + 1

          IF ( NetCdfPrint ) THEN 
            IF ( (t - Tspan(1) >= StpNetCDF*iStpNetCDF) )  THEN
              iStpNetCDF  = iStpNetCDF + 1
              zen = Zenith(t)
              IF ( ntAqua > 0 ) THEN
                actLWC  = pseudoLWC(t)
                wetRad  = (Pi34*actLWC/SPEK(1)%Number)**(rTHREE)*0.1_dp
              END IF
              ! save data in NetCDF File
              TimeNetCDFA = MPI_WTIME()
              CALL SetOutputNCDF(  Y,    yNcdf, t ,  actLWC)
              !
              CALL StepNetCDF   ( t                                &
              &                 , yNcdf                            &
              &                 , itime_NetCDF                     &
              &                 , (/ actLWC                        &
              &                    , h                             &
              &                    , SUM(Y(1:ntGas))               &
              &                    , SUM(Y(ntGas+1:ntGas+ntAqua))  &
              &                    , wetRad                        &
              &                    , zen                        /) &
              &                 , errind                           &
              &                 , SUM(y(iDiag_Schwefel))           &
              &                 , error                            )
              !
              TimeNetCDF  = TimeNetCDF + (MPI_WTIME() - TimeNetCDFA)
            END IF
          END IF

          !-- Call progress bar.
          IF ( Bar .AND. t-Tspan(1) >= iBar*tmp_tb ) THEN
            iBar = iBar + 1         ! iBar runs from 0 to 100
            CALL Progress(iBar)
          END IF

          !-- Termination condition for the main loop.
          IF ( done ) EXIT
            
          !-- Advance the integration one step.
          t   = tnew
          Y0  = Y

        END DO BACKWARD_EULER

      CASE DEFAULT
        
        WRITE(*,*) ' No other methods implemented jet. Use Rosenbrock, LSODE[S] or backward Euler'

    END SELECT
    TimeIntegrationE  = MPI_WTIME() - TimeIntegrationA

    !
    ! DEALLOCATE Mumps instance
    IF ( useMUMPS ) THEN
      mumps_par%JOB = -2
      CALL DMUMPS( mumps_par )
    END IF
  END SUBROUTINE Integrate
  
 
  ! The Progressbar
  SUBROUTINE Progress(j)
    !
    INTEGER(4)    :: j,k
    CHARACTER(29) :: bar="  ???% |                    |"
    !
    IF (MPI_ID==0) THEN
      WRITE(bar(3:5),'(I3)') j
      !
      DO k=1,j/5
        bar(8+k:8+k)="*"
      END DO
      ! print the progress bar.
      WRITE(*,'(A1,A29,$)') char(13), bar
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
    htmp = ( Tol/er )**(0.7d0*RCo%pow)* ( erOld/Tol )**(0.4d0*RCo%pow) * h
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
    REAL(dp) :: dCdt(nspc)
    REAL(dp) :: Rate(neq) , DRate(neq)
    REAL(dp) :: Tarr(8)
    REAL(dp) :: U(nspc) , dUdT(nspc)
    REAL(dp) :: c_v


    ! MASS CONSERVATION
    CALL ReactionRatesAndDerivative( T , Y , Rate , DRate)
    CALL DAXPY_sparse( dCdt , BAT , Rate , Y_e )  ! [mol/cm3/s]
    
    YDOT(1:nspc) =  dCdt
    
    ! ENERGY CONSERVATION
    IF (TempEq) THEN         ! OUT:   IN:
      CALL UpdateTempArray   ( Tarr , Y(NEQ1) )
      CALL InternalEnergy    (  U   , Tarr )
      CALL DiffInternalEnergy( dUdT , Tarr) 
      CALL MassAveMixSpecHeat( c_v  , dUdT               &
                             &      , MoleConc=Y(1:nspc) &
                             &      , rho=rho          ) 
    
      YDOT(NEQ1) = - SUM( U * dCdt ) * rRho / c_v   ! rRho=mega/rho 
    END IF

  END SUBROUTINE FRhs

  ! dummy routine for ode solver LSODE
  SUBROUTINE dummy()
    IMPLICIT NONE
  END SUBROUTINE dummy

END MODULE Integration_Mod
