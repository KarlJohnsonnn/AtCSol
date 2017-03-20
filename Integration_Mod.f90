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
  USE mo_ckinput, ONLY: rhoY
  USE Sparse_Mod, ONLY: A,B,BA,BAT
  USE Sparse2_Mod
  USE Chemsys_Mod
  USE Rosenbrock_Mod
  USE Rates_Mod
  USE Factorisation_Mod
  USE ErrorROW_Mod
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
    TYPE(CSR_Matrix_T)     :: Id, Jac_C           ! compressed row
    TYPE(SpRowIndColInd_T) :: MiterFact         ! sparse row-ind col-ind (MUMPS)
    TYPE(SpRowColD_T)      :: MiterMarko        ! sparse-LU matrix format
    !
    INTEGER, ALLOCATABLE :: InvPermu(:)
    ! 
    REAL(RealKind) :: Y0(nDIM)
    REAL(RealKind) :: Y(nDIM)       ! current y vector
    REAL(RealKind) :: Y_mol(nspc)       ! current y vector
    !
    REAL(RealKind) :: t             ! current time
    REAL(RealKind) :: StartTimer
    REAL(RealKind) :: Rate(neq) , dummyrate(neq)    ! rate vector
    REAL(RealKind) :: DRatedT(neq)     ! part. derv. rate over temperatur vector
    REAL(RealKind) :: h, hmin, absh, tnew, tmp, hOld
    REAL(RealKind) :: error, errorOld
    REAL(RealKind) :: tmp_ta,tmp_tb
    REAL(RealKind) :: actLWC, zen, wetRad
    INTEGER        :: errind(1,1)
    REAL(RealKind) :: ErrVals(nspc)
    !
    REAL(RealKind) :: Temp0   ! Temperature old
    REAL(RealKind) :: Temp    ! Temperature new
    ! 
    INTEGER :: iBar=-1              ! waitbar increment
    INTEGER :: i
    LOGICAL :: Bar=.FALSE.
    !
    ! for NetCDF
    INTEGER :: iStpNetCDF 
    INTEGER :: itime_NetCDF
    REAL(RealKind), ALLOCATABLE :: yNcdf(:)     ! current output vector
    !
    LOGICAL :: done=.FALSE.
    LOGICAL :: failed
    !
    INTEGER :: i_error
    CHARACTER(14) :: fmt0
    !
    !
    Y0(1:nspc)  = y_iconc(:)
    Y(1:nspc)   = y_iconc(:)
    !
    IF ( combustion ) THEN
      !--- initial temperature
      Y0(nDIM)      = Temperature0     ! = 750 [K] aus speedchem debug
      Y(nDIM)       = Temperature0     ! = 750 [K]
      
      !--- malloc gibbs energy, derivates
      ALLOCATE( GFE(nspc)   , DGFEdT(nspc)   )
      ALLOCATE( DelGFE(neq) , DDelGFEdT(neq) )
      GFE        =  ZERO
      DGFEdT     =  ZERO
      DelGFE     =  ZERO
      DDelGFEdT  =  ZERO
      CALL rhoY( rho , Y0(1:nspc) )  ! Initialising reactor density
      print*, 'DEBUG::INTEGRmod  scrho,scp,temp0    ', rho, Press, Y0(nDIM)
    END IF
    !
    !
    IF (MPI_ID==0) THEN
      IF ( Ladebalken==1 ) Bar = .TRUE.
      WRITE(*,*) '  Initial values:   '
      WRITE(*,*)
      fmt0 = '  [molec/cm3] '
      IF ( combustion ) fmt0 = '  [mol/cm3]   '
      WRITE(*,'(A34,2X,E23.14,A13)')  '      sum initval (gaseous)    = ', SUM(Y0(1:ntGas)),      fmt0
      WRITE(*,'(A34,2X,E23.14,A13)')  '      sum initval (aqueous)    = ', SUM(Y0(ntGas+1:nspc)), fmt0
      WRITE(*,'(A34,2X,E23.14,A13)')  '      sum emissions (gaseous)  = ', SUM(Y_e),              fmt0
      IF ( combustion ) THEN
        WRITE(*,'(A34,2X,E23.14,A13)') '          Initial Temperature  = ', Temperature0
      END IF
      WRITE(*,*)
    ELSE
      NetCdfPrint = .FALSE.
    END IF
    !
    !IF (MPI_np<2) THEN
      ALLOCATE(loc_RatePtr(neq))
      FORALL( i=1:neq ) loc_RatePtr(i)=i
      loc_rateCnt=neq
    !ELSE
      ! get local rate pointer for each process
    !END IF
    !
    !---- Check input parameters 
    CALL CheckInputParameters( Tspan , Atol , RtolROW )
    !
    !---- Get Rosenbrock Parameters
    CALL SetRosenbrockMethod( RCo , method )  
    !
    !----------------------------------------------------------
    ! ----------- Beginning with symbolic phase --------------
    !----------------------------------------------------------
    !call printsparsematrix(A,'A')
    !call printsparsematrix(B,'B')
    !
    StartTimer = MPI_WTIME()            ! start timer for symb phase
    CALL SymbolicAdd( BA , B , A )      ! symbolic addition:    BA = B + A
    CALL SparseAdd  ( BA , B , A, '-')  ! numeric subtraction:  BA = B - A
    CALL TransposeSparse( BAT , BA )    ! transpose BA:        BAT = Transpose(BA) 
    !  

    !call printsparsematrix(BA,'BA')
    !call printsparsematrix(BAT,'BAT')
    ! we need to calculate the Jacobian for both versions 'cl' and 'ex' to
    ! calculate an initial stepsize based on 2nd derivative (copy of MATLABs ode23s)
    !
    CALL SymbolicMult( BAT , A , Jac_C )   ! symbolic mult for Jacobian Jac = [ BAT * Ones_nR * A * Ones_nS ]
    !
    ! Set symbolic structure of iteration matrix for Row-Method
    IF ( CLASSIC ) THEN
      CALL BuildSymbolicClassicMatrix(  Miter , Jac_C   , RCo%ga )
    ELSE
      CALL BuildSymbolicExtendedMatrix( Miter , A , BAT , RCo%ga ) 
    END IF
    !
    ! Choose ordering/factorisation strategie and do symb LU fact
    !
    IF ( OrderingStrategie<8 .OR. ParOrdering>=0 ) THEN 
      ! Use MUMPS to factorise and solve     
      ! Convert compressed row format to row index format for MUMPS
      CALL CompRowToIndRow( Miter , MiterFact )
      CALL InitMumps( MiterFact ) 

      IF ( MatrixPrint ) THEN
        CALL CSRToSpRowColD   ( MiterMarko , Miter ) 
        CALL PermuToInvPer    ( InvPermu   , Mumps_Par%SYM_PERM )
        CALL SymbLU_SpRowColD ( MiterMarko , InvPermu )        
        CALL RowColDToCSR     ( LU_Miter   , MiterMarko , nspc , neq ) 
      END IF

    ELSE
      ! Permutation given by Markowitz Ordering strategie
      CALL CSRToSpRowColD     ( MiterMarko , Miter) 
      CALL SymbLU_SpRowColD_M ( MiterMarko )        
      CALL RowColDToCSR       ( LU_Miter , MiterMarko , nspc , neq )

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
    CALL Free_SpRowColD( MiterMarko )
    !
    TimeSymbolic = MPI_WTIME() - StartTimer   ! stop timer for symbolic matrix calculations
    IF (MPI_ID==0) WRITE(*,*) '  SYMBOLIC PHASE................ done'
    
   
    CALL Rates( Tspan(1) , Y0 , Rate , DRatedT )      ! Calculate first reaction rates
    
    Rate      = MAX( ABS(Rate) , eps ) * SIGN( ONE , Rate )
    Y(1:nspc) = MAX( ABS(Y(1:nspc)) , eps ) * SIGN( ONE , Y(1:nspc) )
    
    ! ----calc values of Jacobian
    StartTimer  = MPI_WTIME()
    CALL JacobiMatrix( BAT , A , Rate , Y , Jac_C )
    TimeJac     = TimeJac + MPI_WTIME() - StartTimer
    Output%npds = Output%npds + 1
   
    !---- Save matrix structures for matlab spy
    IF ( MatrixPrint ) CALL SaveMatricies(A,B,BAT,Miter,LU_Miter,BSP)
  
    !---- calculate a first stepsize based on 2nd deriv.
    CALL InitialStepSize( h , hmin , absh , Jac_C , Rate  &
    &                   , Tspan(1) , Y(1:nspc) , RCo%pow  )
 
    !-----------------------------------------------------------------------
    !-- Initialize Netcdf Output File
    itime_NetCDF = 0

    ! aqua phase in [mol/m3] and [mol/l]
    OutNetcdfANZ = COUNT(OutNetcdfPhase=='g') + 2*COUNT(OutNetcdfPhase=='a')
    OutNetcdfDIM = OutNetcdfANZ
    IF ( combustion ) OutNetcdfDIM = OutNetcdfANZ + 1
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
    !
    IF ( NetCdfPrint ) THEN 
      iStpNetCDF  = 1
      StartTimer  = MPI_WTIME()
      errind(1,1) = 0
      CALL InitNetcdf( tAnf , tEnd )

      WRITE(*,*) '  Init NetCDF................... done'
      !--  print initial values to NetCDF file
      CALL SetOutputNCDF( Y, yNcdf, Tspan(1), actLWC )
      CALL StepNetCDF   ( Tspan(1),                        &
                        & yNcdf(:),                        &
                        & itime_NetCDF,                    &
                        & (/ actLWC,                       &
                        &    h,                            &
                        &    SUM(y(1:ntGas)),              &
                        &    SUM(y(ntGas+1:ntGas+ntAqua)), &
                        &    wetRad,                       &
                        &    ZERO   /),                    &
                        &  errind,                         &
                        &  SUM(y(iDiag_Schwefel)),         &
                        &  ZERO                            )
      TimeNetCDF = MPI_WTIME() - StartTimer
    END IF
    !
    ! this is for the waitbar
    tmp_tb = (Tspan(2)-Tspan(1)) * 1.0d-2
    !
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
    IF (MPI_ID==0) WRITE(*,*) '  Start Integration............. '
    IF (MPI_ID==0) WRITE(*,*) ' '
    !
    t     = Tspan(1)
    tnew  = Tspan(1)
    !
    DO 
      !
      absh  = MIN( maxStp, MAX( minStp , absh ) )
      h     = absh
      !
      !-- Stretch the step if within 5% of tfinal-t.
      IF ( 1.05 * absh >= Tspan(2) - t ) THEN
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
        !print*, 'DEBUG::INTEGR     temp=y(end)=',y0(nDim)
        CALL Rosenbrock(  Y0            &       ! old concentration
        &               , t             &       ! time
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
    END DO  ! MAIN LOOP
    TimeIntegrationE  = MPI_WTIME()
    !
    ! DEALLOCATE Mumps instance
    IF (  OrderingStrategie<8 ) THEN
      mumps_par%JOB = -2
      CALL DMUMPS( mumps_par )
    END IF
    !
  END SUBROUTINE Integrate
  !
  !
  ! The Progressbar
  SUBROUTINE Progress(j)
    !
    INTEGER(4) :: j,k
    CHARACTER(29) :: bar="  ???% |                    |"
    !
    IF (MPI_ID==0) THEN
      WRITE(unit=bar(3:5),fmt="(i3)") j
      !
      DO k=1,j/5
        bar(8+k:8+k)="*"
      END DO
      ! print the progress bar.
      WRITE(*,fmt="(a1,a29,$)") char(13), bar
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
    print*, 'htemp= ',htmp, h , hold
    !
    ! limitation term
    IF (htmp>PI_Par%ThetaMAX*h) THEN
      hnew = PI_Par%ThetaMAX*h
    ELSE
      hnew = htmp
    END IF
  END SUBROUTINE PI_StepsizeControl
END MODULE Integration_Mod
