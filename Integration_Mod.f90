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
  USE Sparse_Mod, ONLY: A,B,BA,BAT
  USE Sparse2_Mod
  USE Chemsys_Mod
  USE Rosenbrock_Mod
  USE Rates_Mod
  USE Factorisation_Mod
  USE NetCDF_Mod
  USE Meteo_Mod, ONLY: Temp
  IMPLICIT NONE
  !
  TYPE PI_Param
    REAL(RealKind) :: Kp
    REAL(RealKind) :: KI
    REAL(RealKind) :: ThetaMAX
    REAL(RealKind) :: rho
  END TYPE PI_Param
  !
  TYPE(PI_Param) :: PI_norm, PI_rej
  !
  CONTAINS
  !
  !======================================================================= 
  !===================    Time Integration Routine  ======================
  !======================================================================= 
  SUBROUTINE Integrate(y_iconc, Tspan, Atol, RtolRow, method)
    !--------------------------------------------------------------------
    ! Input:
    !   - y0 ............. Initial vector
    !   - Tspan .......... (/ SimulationTimeStart , SimulationTimeEnd /)
    !   - Atol ........ abs. tolerance for gas spc
    !   - RtolRow ........ rel. tolerance for Rosenbrock-Wanner-Method
    !   - method...........Rosenbrock-Wanner-Method
    !   - PrintSpc ....... print spc PrintSpc(1:3)
    REAL(RealKind) :: y_iconc(nspc)
    REAL(RealKind) :: Tspan(2)
    REAL(RealKind) :: Atol(2)
    REAL(RealKind) :: RtolROW
    CHARACTER(*) :: method
    !-------------------------------------------------------------------
    ! Output:
    !   - Output.......... struct Out (above) contains output statistics
    !   - Unit............ .simul data 
    !-------------------------------------------------------------------
    ! Temporary variables:
    TYPE(RosenbrockMethod_T) :: RCo  
    !
    TYPE(CSR_Matrix_T)     :: Id , tmpJacCC        ! compressed row
    TYPE(SpRowIndColInd_T) :: MiterFact         ! sparse row-ind col-ind (MUMPS)
    TYPE(SpRowColD_T)      :: temp_LU_Dec        ! sparse-LU matrix format
    !
    INTEGER, ALLOCATABLE :: InvPermu(:)
    INTEGER, ALLOCATABLE :: PivOrder(:)
    ! 
    REAL(RealKind) :: Y0(nDIM)
    REAL(RealKind) :: Y(nDIM)       ! current y vector
    !
    REAL(RealKind) :: t             ! current time
    REAL(RealKind) :: timepart
    REAL(RealKind) :: StartTimer

    REAL(RealKind) :: Rate(neq)
    REAL(RealKind) :: DRatedT(neq)     ! part. derv. rate over temperatur vector
    REAL(RealKind) :: h, hmin, absh, tnew, tmp, hOld
    REAL(RealKind) :: error, errorOld
    REAL(RealKind) :: tmp_tb
    REAL(RealKind) :: actLWC, zen, wetRad
    INTEGER        :: errind(1,1)
    REAL(RealKind) :: ErrVals(nspc)
    ! 
    INTEGER :: iBar=-1              ! waitbar increment
    INTEGER :: i, k
    LOGICAL :: Bar=.FALSE.

    INTEGER, PARAMETER :: maxnsteps = 50000
    !
    ! for NetCDF
    INTEGER :: iStpNetCDF 
    INTEGER :: itime_NetCDF
    REAL(RealKind), ALLOCATABLE :: yNcdf(:)     ! current output vector
    !
    LOGICAL :: done=.FALSE.
    LOGICAL :: failed
    !
    CHARACTER(14) :: fmt0
    CHARACTER(17) :: fmt1
    !
    !-----------------------------
    ! LSODE Parameter
    INTEGER :: ITOL , ITASK , ISTATE, IOPT, MF, LIW, LRW
    INTEGER, ALLOCATABLE :: IWORK(:)
    REAL(RealKind), ALLOCATABLE :: RWORK(:)
    !
    REAL(RealKind) :: RTOL1 , ATOL1(nDIM)

    !
    !
    Y0(1:nspc)  = y_iconc(:)
    Y(1:nspc)   = y_iconc(:)
    !
    IF ( TempEq ) THEN
      !--- initial temperature
      Y0(nDIM)      = Temperature0     ! = 750 [K] aus speedchem debug
      Y(nDIM)       = Temperature0     ! = 750 [K]
      
      !--- malloc gibbs energy, derivates
      ALLOCATE( GFE(nspc)   , DGFEdT(nspc)   )
      ALLOCATE( DelGFE(neq) , DDelGFEdT(neq) )
      GFE    = ZERO;  DGFEdT    = ZERO
      DelGFE = ZERO;  DDelGFEdT = ZERO
    END IF

    !------------------------------------------------------------------------------------------|
    !                                                                                          |
    !                    Print Init parameters, concs, temp, pressure,....                     |
    !                                                                                          |
    !------------------------------------------------------------------------------------------|
    
    IF (MPI_ID==0) THEN
       
      IF ( Ladebalken==1 ) Bar = .TRUE.
      WRITE(*,*) '  Initial values:   '
      WRITE(*,*)
      fmt0 = '  [molec/cm3] '
      IF ( TempEq ) fmt0 = '  [mol/cm3]'
      fmt1 = '(A34,2X,E16.10,A)'
      IF (ntGas >0) WRITE(*,fmt1)  '      sum initval (gaseous)    = ', SUM(Y0(1:ntGas)),      fmt0
      IF (ntAqua>0) WRITE(*,fmt1)  '      sum initval (aqueous)    = ', SUM(Y0(ntGas+1:nspc)), fmt0
      WRITE(*,fmt1)  '      sum emissions (gaseous)  = ', SUM(Y_e),              fmt0

      IF ( TempEq ) THEN
        WRITE(*,fmt1) '          Initial Temperature  = ', Temperature0,'  [K]'
        WRITE(*,fmt1) '          Initial Pressure     = ', Pressure0,'  [Pa]'
        WRITE(*,fmt1) '          Reactor denstiy      = ', rho,'  [kg/cm3]'
      END IF
      WRITE(*,*)

    ELSE
      NetCdfPrint = .FALSE.
    END IF
    !
    ALLOCATE(loc_RatePtr(neq))
    loc_RatePtr = [(i , i=1,neq)]
    loc_rateCnt = neq
  
    CALL CheckInputParameters( Tspan , Atol , RtolROW )


    !------------------------------------------------------------------------------------------|
    !                                                                                          |
    !                                Initialize NetCDF output file                             |
    !                                                                                          |
    !------------------------------------------------------------------------------------------|

    IF ( NetCdfPrint ) THEN 
      StartTimer  = MPI_WTIME()

      itime_NetCDF = 0
      iStpNetCDF   = 1
    
      ! aqua phase in [mol/m3] and [mol/l]
      OutNetcdfANZ = COUNT(OutNetcdfPhase=='g') + 2*COUNT(OutNetcdfPhase=='a')
      OutNetcdfDIM = OutNetcdfANZ
      IF ( TempEq ) OutNetcdfDIM = OutNetcdfANZ + 1
      ALLOCATE(yNcdf(OutNetcdfDIM))          ! output array , +1 for Temperature
      yNcdf = ZERO
          
      !--  Netcdf Output File
      IF ( ErrorLog == 1 ) ErrVals = ZERO
      IF ( ntAqua > 0 ) THEN
        actLWC = pseudoLWC(Tspan(1))
        wetRad = (Pi34*actLWC/SPEK(1)%Number)**(rTHREE)*1.0d-1
      ELSE
        actLWC = ZERO
        wetRad = ZERO
      END IF
     

      errind(1,1) = 0
      CALL InitNetcdf( Tspan(1) , Tspan(2) )

      WRITE(*,*) '  Init NetCDF................... done'
      !--  print initial values to NetCDF file
      CALL SetOutputNCDF( Y0, yNcdf, Tspan(1), actLWC )
      CALL StepNetCDF   ( Tspan(1),                        &
                        & yNcdf(:),                        &
                        & itime_NetCDF,                    &
                        & (/ actLWC,                        &
                        &    h,                             &
                        &    SUM(Y0(1:ntGas)),              &
                        &    SUM(Y0(ntGas+1:ntGas+ntAqua)), &
                        &    wetRad,                        &
                        &    ZERO   /),                    &
                        &  errind,                         &
                        &  SUM(Y0(iDiag_Schwefel)),        &
                        &  ZERO                            )
      TimeNetCDF = MPI_WTIME() - StartTimer
    END IF
    
    ! this is for the waitbar
    tmp_tb = (Tspan(2)-Tspan(1)) * 1.0d-2
    
    !do i=1,SIZE(a%val); write(777,*) A%val(i); end do
    !do i=1,SIZE(b%val); write(778,*) B%val(i); end do
    !call printsparse1(A,'alpha_MCMC')
      
    
    ! build nessesarry matrecies for Rate Calculation and FRhs
    CALL SymbolicAdd( BA , B , A )      ! symbolic addition:    BA = B + A
    CALL SparseAdd  ( BA , B , A, '-' )  ! numeric subtraction:  BA = B - A
    CALL TransposeSparse( BAT , BA )    ! transpose BA:        BAT = Transpose(BA) 
    !do k=1,size(ba%val); print*, 'BAVAL=',BA%Val(k); end do
    !stop

    SELECT CASE (TRIM(method(1:7)))

      CASE('METHODS')
       
        !---- Get Rosenbrock Parameters
        CALL SetRosenbrockMethod( RCo , method )  
    
        !----------------------------------------------------------
        ! ----------- Beginning with symbolic phase --------------
        !----------------------------------------------------------
      
        StartTimer = MPI_WTIME()            ! start timer for symb phase
    
        ! we need to calculate the Jacobian for both versions 'cl' and 'ex' to
        ! calculate an initial stepsize based on 2nd derivative (copy of MATLABs ode23s)
        ! symbolic mult for Jacobian Jac = [ BAT * Ones_nR * A * Ones_nS ], 
        ! add diagonal entries and set them to zero
        CALL SymbolicMult( BAT , A , tmpJacCC )  
        CALL SparseID(Id,nspc)                    
        CALL SymbolicAdd(Jac_CC,Id,tmpJacCC)
        CALL Free_Matrix_CSR(Id)
        CALL Free_Matrix_CSR(tmpJacCC)
    
        ! Set symbolic structure of iteration matrix for Row-Method
        ! also assign constant matrix parts like (beta-alpha)^T, alpha
        IF ( CLASSIC ) THEN
          CALL BuildSymbolicClassicMatrix(  Miter , Jac_CC  , RCo%ga )
        ELSE
          CALL BuildSymbolicExtendedMatrix( Miter , A , BAT , RCo%ga ) 
        END IF
    
        ! Choose ordering/factorisation strategie and do symb LU fact
        IF ( useMUMPS ) THEN 
          ! Use MUMPS to factorise and solve     
          ! Convert compressed row format to row index format for MUMPS
          CALL CompRowToIndRow( Miter , MiterFact )
          CALL InitMumps( MiterFact ) 
    
          IF ( MatrixPrint ) THEN
            CALL CSRToSpRowColD   ( temp_LU_Dec , Miter ) 
            CALL PermuToInvPer    ( InvPermu   , Mumps_Par%SYM_PERM )
            CALL SymbLU_SpRowColD ( temp_LU_Dec , InvPermu )        
            CALL RowColDToCSR     ( LU_Miter   , temp_LU_Dec , nspc , neq ) 
          END IF
    
        ELSE IF ( useSparseLU ) THEN
          ! Permutation given by Markowitz Ordering strategie
          CALL CSRToSpRowColD( temp_LU_Dec , Miter) 

          IF (OrderingStrategie==8) THEN

            CALL SymbLU_SpRowColD_M ( temp_LU_Dec )        

          ELSE 

            ALLOCATE(PivOrder(temp_LU_Dec%n))
            PivOrder = -90
            
            ! Pivots bis neq in reihenfolge 1,2,...,neq
            IF ( EXTENDED ) THEN
              PivOrder(     1 : neq      ) = [(i , i = 1     , neq  )]
              PivOrder( neq+1 : neq+nDIM ) = [(i , i = neq+1 , neq+nDIM )]
            ELSE
              PivOrder(     1 : nDIM     ) = [(i , i = 1     , nDim )]
            END IF
            CALL SymbLU_SpRowColD ( temp_LU_Dec , PivOrder)

          END IF
          CALL RowColDToCSR       ( LU_Miter , temp_LU_Dec , nspc , neq )
    
          ! Get the permutation vector LU_Perm and map values of Miter
          ! to the permuted LU matrix
          CALL Get_LU_Permutaion  ( LU_Perm , LU_Miter , Miter , nspc , neq )
          
          ! For the extended case one can save the values of alpha and 
          ! (beta-alpha)^T and just copy them for each iteration in ROS
          IF ( EXTENDED ) THEN
            ALLOCATE(LUValsFix(LU_Miter%nnz))
            LUvalsFix = LU_Miter%Val
          END IF
        END IF
        CALL Free_SpRowColD( temp_LU_Dec )

        IF (MatrixPrint) THEN
          CALL WriteSparseMatrix(Miter,'MATRICES/Miter_'//TRIM(BSP)//'_'//solveLA,neq,nspc)
          CALL WriteSparseMatrix(LU_Miter,'MATRICES/LU_Miter_'//TRIM(BSP)//'_'//solveLA,neq,nspc)
          WRITE(*,*) '  Writing matrices to file: ','MATRICES/Miter_'//TRIM(BSP)//'_'//solveLA
          WRITE(*,*) '                            ','MATRICES/LU_Miter_'//TRIM(BSP)//'_'//solveLA
          stop 'after writesparsematrix'
        END IF
        
        IF (MPI_ID==0) WRITE(*,*) '  SYMBOLIC PHASE................ done'
        TimeSymbolic = MPI_WTIME() - StartTimer   ! stop timer for symbolic matrix calculations
    
        ! ---- Calculate first reaction rates
        CALL Rates( Tspan(1) , Y0 , Rate , DRatedT )
        Y = MAX( ABS(Y) , eps ) * SIGN( ONE , Y )    ! |y| >= eps
        
        ! ---- Calculate values of Jacobian
        StartTimer  = MPI_WTIME()
        CALL Jacobian_CC(Jac_CC , BAT , A , Rate , Y )
        TimeJac     = TimeJac + MPI_WTIME() - StartTimer
        Output%npds = Output%npds + 1
       
        !---- Calculate a first stepsize based on 2nd deriv.
        CALL InitialStepSize( h , hmin , absh , Jac_CC , Rate  &
        &                   , Tspan(1) , Y(1:nspc) , RCo%pow  )
     
        !
        ! future stuff (PI-stepsize control)
        ! set for normal case
        PI_norm%Kp        = 0.13d0
        PI_norm%KI        = ONE/15.0d0        
        PI_norm%ThetaMAX  = TWO         
        PI_norm%rho       = 1.2d0 
        ! set for rejected case
        PI_rej%Kp         = ZERO
        PI_rej%KI         = rFIVE        
        PI_rej%ThetaMAX   = TWO        
        PI_rej%rho        = 1.2d0 
        !===============================================================================
        !=================================THE MAIN LOOP=================================
        !===============================================================================
        TimeIntegrationA=MPI_WTIME()
        IF (MPI_ID==0) WRITE(*,*) '  Start Integration............. '; WRITE(*,*) ' '
       
        t     = Tspan(1)
        tnew  = Tspan(1)
        
            
        MAIN_LOOP_ROSENBROCK_METHOD: DO 
          !
          absh  = MIN( maxStp, MAX( minStp , absh ) )
          h     = absh
          !
          !-- Stretch the step if within 5% of tfinal-t.
          IF ( 1.05d0 * absh >= Tspan(2) - t ) THEN
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
            &               , Y             )       ! new concentration 
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
            IF ( PI_StepSize .AND. Output%nSteps>1 ) THEN
              failed = (error > h*PI_rej%rho*RtolRow)
            ELSE
              failed = (error > ONE)
            END IF
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
              IF ( PI_StepSize .AND. Output%nSteps>1  ) THEN
                CALL  PI_StepsizeControl(absh,RtolRow,error,errorOld,h,hOld,PI_rej,RCo)
              ELSE
                absh  = MAX( hmin , absh * MAX( rTEN , 0.8d0 * (ONE/error)**RCo%pow) )
              END IF
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
            IF ( (t - Tspan(1) >= StpNetCDF*iStpNetCDF) )  THEN
              iStpNetCDF  = iStpNetCDF + 1
              zen = Zenith(t)
              IF ( ntAqua > 0 ) THEN
                actLWC  = pseudoLWC(t)
                wetRad  = (Pi34*actLWC/SPEK(1)%Number)**(rTHREE)*1.0d-1
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
          IF ( PI_StepSize ) THEN
            CALL PI_StepsizeControl(absh,RtolRow,error,errorOld,h,hOld,PI_norm,RCo)
          ELSE
            tmp = 1.25d0  * error**RCo%pow
            IF ( TWO * tmp > ONE ) THEN
              absh  = absh / tmp
            ELSE
              absh  = absh * TWO
            END IF
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
    
        TimeIntegrationE  = MPI_WTIME() - TimeIntegrationA
        !
        ! DEALLOCATE Mumps instance
        IF ( useMUMPS ) THEN
          mumps_par%JOB = -2
          CALL DMUMPS( mumps_par )
        END IF


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
        timepart = (Tspan(2)-Tspan(1)) / REAL(nOutP - 1, RealKind)
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
            iStpNetCDF  = iStpNetCDF + 1
            zen = Zenith(t)
            IF ( ntAqua > 0 ) THEN
              actLWC  = pseudoLWC(t)
              wetRad  = (Pi34*actLWC/SPEK(1)%Number)**(rTHREE)*1.0d-1
            END IF

            ! save data in NetCDF File
            TimeNetCDFA = MPI_WTIME()
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

        TimeIntegrationE  = MPI_WTIME() - TimeIntegrationA
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
        timepart = (Tspan(2)-Tspan(1)) / REAL(nOutP - 1, RealKind)
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
            iStpNetCDF  = iStpNetCDF + 1
            zen = Zenith(t)
            IF ( ntAqua > 0 ) THEN
              actLWC  = pseudoLWC(t)
              wetRad  = (Pi34*actLWC/SPEK(1)%Number)**(rTHREE)*1.0d-1
            END IF

            ! save data in NetCDF File
            TimeNetCDFA = MPI_WTIME()
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

        TimeIntegrationE  = MPI_WTIME() - TimeIntegrationA
        CALL Progress(100) ! last * needs an extra call



      CASE DEFAULT
        
        WRITE(*,*) ' No other methods implemented jet. Use Rosenbrock or LSODE'

    END SELECT

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
    REAL(RealKind) :: hnew
    !
    !
    REAL(RealKind) :: Tol       ! rel tol
    REAL(RealKind) :: er        ! lokal error step n 
    REAL(RealKind) :: erOld     ! lokal error step n-1
    REAL(RealKind) :: h         ! stepsize n
    REAL(RealKind) :: hOld      ! stepsize n-1
    TYPE(PI_Param) :: PI_Par    ! PI control parameter
    TYPE(RosenbrockMethod_T)  :: RCo
    !
    REAL(RealKind) :: htmp
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
    INTEGER        :: NEQ1
    REAL(RealKind) :: T , Y(NEQ1) , YDOT(NEQ1)
    !
    REAL(RealKind) :: dCdt(nspc)
    REAL(RealKind) :: Rate(neq) , DRate(neq)
    REAL(RealKind) :: Tarr(8)
    REAL(RealKind) :: U(nspc) , dUdT(nspc)
    REAL(RealKind) :: cv


    ! MASS CONSERVATION
    CALL Rates( T , Y , Rate , DRate)
    CALL DAXPY_sparse( dCdt , BAT , Rate , Y_e )  ! [mol/cm3/s]
    
    YDOT(1:nspc) =  dCdt
    
    ! ENERGY CONSERVATION
    IF (TempEq) THEN         ! OUT:   IN:
      CALL UpdateTempArray   ( Tarr , Y(NEQ1) )
      CALL InternalEnergy    (  U   , Tarr )
      CALL DiffInternalEnergy( dUdT , Tarr) 
      CALL MassAveMixSpecHeat( cv   , dUdT               &
                             &      , MoleConc=Y(1:nspc) &
                             &      , rho=rho          ) 
    
      YDOT(NEQ1) = - SUM( U * dCdt ) * rRho / cv    ! rRho=mega/rho 
    END IF

  END SUBROUTINE FRhs

  ! dummy routine for ode solver LSODE
  SUBROUTINE dummy()
    IMPLICIT NONE
  END SUBROUTINE dummy

END MODULE Integration_Mod
