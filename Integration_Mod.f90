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
  SUBROUTINE Integrate(y_iconc, Tspan, Atol, RtolRow, vers, method)
    !--------------------------------------------------------------------
    ! Input:
    !   - y0 ............. Initial vector
    !   - Tspan .......... (/ SimulationTimeStart , SimulationTimeEnd /)
    !   - Atol ........ abs. tolerance for gas spc
    !   - RtolRow ........ rel. tolerance for Rosenbrock-Wanner-Method
    !   - vers ........... version of solving the linear system
    !   - method...........Rosenbrock-Wanner-Method
    !   - PrintSpc ....... print spc PrintSpc(1:3)
    REAL(RealKind) :: y_iconc(nspc)
    REAL(RealKind) :: Tspan(2)
    REAL(RealKind) :: Atol(2)
    REAL(RealKind) :: RtolROW
    CHARACTER(2) :: vers
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
    REAL(RealKind) :: y0(nDIM)
    REAL(RealKind) :: y(nDIM)       ! current y vector
    !
    REAL(RealKind) :: t             ! current time
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
    !
    ! for NetCDF
    INTEGER :: iStpNetCDF 
    INTEGER :: tncdf_ind
    REAL(RealKind), ALLOCATABLE :: yNcdf(:)     ! current output vector
    !
    LOGICAL :: done=.FALSE.
    LOGICAL :: failed
    !
    INTEGER :: i_error
    !
    !
    y0(1:nspc)=y_iconc(:)
    y(1:nspc)=y_iconc(:)
    !
    IF ( combustion ) THEN
      y0(nDIM)=750.0d0   ! = 750 [K] aus speedchem debug
      y(nDIM)=750.0d0    ! = 750 [K]
      ALLOCATE(GFE(nspc),DGFEdT(nspc))
      ALLOCATE(DelGFE(neq),DDelGFEdT(neq))
      GFE(:)=ZERO
      DGFEdT(:)=ZERO
      DelGFE(:)=ZERO
      DDelGFEdT(:)=ZERO
    END IF
    !
    !
    IF (MPI_ID==0) THEN
      WRITE(*,*) '  Initial values:   '
      WRITE(*,*)
      WRITE(*,'(A34,2X,E23.14,A13)')  '      sum initval (gaseous)    = ', SUM(y0(1:ntGas)), '  [molec/cm3] '
      WRITE(*,'(A34,2X,E23.14,A13)')  '      sum initval (aqueous)    = ', SUM(y0(ntGas+1:nspc)), '  [molec/cm3] '
      WRITE(*,'(A34,2X,E23.14,A13)')  '      sum emissions (gaseous)  = ', SUM(y_e), '  [molec/cm3] '
      WRITE(*,*)
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
    CALL CheckInputParameters(Tspan,Atol,RtolROW)
    !
    !---- Get Rosenbrock Parameters
    CALL SetRosenbrockMethod(RCo,method)  
    !
    !----------------------------------------------------------
    ! ----------- Beginning with symbolic phase --------------
    !----------------------------------------------------------
    !call printsparse(A,'*')
    !call printsparse(B,'*')

    !
    TimeSymbolic=MPI_WTIME()       ! start timer for symb phase
    CALL SymbolicAdd(BA,B,A)       ! symbolic addition:    BA = B + A
    CALL SparseAdd(BA,B,A,'-')     ! numeric subtraction:  BA = B - A
    CALL TransposeSparse(BAT,BA)   ! transpose BA:        BAT = Transpose(BA) 
    !  
    !call printsparse(BAT,'*')
    !stop
    ! we need to calculate the Jacobian for both versions 'cl' and 'ex' to
    ! calculate an initial stepsize based on 2nd derivative (copy of MATLABs ode23s)
    !
    CALL SymbolicMult(BAT,A,Jac_C)   ! symbolic mult for Jacobian Jac = [ BAT * Ones_nR * A * Ones_nS ]
    !
    ! Set symbolic structure of iteration matrix for Row-Method
    IF ( vers=='cl') THEN
      CALL BuildSymbolicClassicMatrix(Miter,Jac_C,RCo%ga)
    ELSE
      CALL BuildSymbolicExtendedMatrix(Miter,A,BAT,RCo%ga)
    END IF
    !
    ! Choose ordering/factorisation strategie and do symb LU fact
    !
    IF ( OrderingStrategie<8 .OR. ParOrdering>=0 ) THEN 
      ! Use MUMPS to factorise and solve     
      ! Convert compressed row format to row index format for MUMPS
      CALL CompRowToIndRow(Miter,MiterFact)
      CALL InitMumps(MiterFact) 
      IF ( MatrixPrint ) THEN
        CALL CSR_Matrix_To_SpRowColD(MiterMarko,Miter) 
        CALL PermuToInvPer(InvPermu,Mumps_Par%SYM_PERM)
        CALL SymbLU_SpRowColD(MiterMarko,InvPermu)        
        CALL RowColD_To_CRF_Matrix(LU_Miter,MiterMarko,nspc,neq,vers)
      END IF
    ELSE
      ! Permutation given by Markowitz Ordering strategie
      CALL CSR_Matrix_To_SpRowColD(MiterMarko,Miter) 
      CALL SymbLUMarko_SpRowColD(MiterMarko)        
      CALL RowColD_To_CRF_Matrix(LU_Miter,MiterMarko,nspc,neq,vers)
      ! Get the permutation vector LU_Perm and map the
      ! Miter values to the permuted LU factorisation
      CALL Get_LU_Permutaion(LU_Miter,Miter,nspc,neq,LU_Perm)
      !
      IF (vers=='ex') THEN
        ALLOCATE(LUValsFix(LU_Miter%RowPtr(LU_Miter%m+1)-1))
        LUvalsFix(:)=LU_Miter%Val(:)
      END IF
    END IF
    CALL Kill_Matrix_SpRowColD(MiterMarko)
    !
    !

    !
    TimeSymbolic=MPI_WTIME()-TimeSymbolic   ! stop timer for symbolic matrix calculations
    IF (MPI_ID==0) WRITE(*,*) '  SYMBOLIC PHASE................ done'
    !
    !  
    !
    t=Tspan(1)
    CALL Rates(Tspan(1),y0,Rate,DRatedT)      ! Calculate first reaction rates
    dummyrate(:)=Rate(:)

    !print*, 'debug :: roargs srate,sy  ', SUM(Rate), SUM(y0)
    !stop
    
    Rate(:)=MAX(ABS(Rate(:)),eps)*SIGN(ONE,Rate(:))
    y(1:nspc)=MAX(ABS(y(1:nspc)),eps)*SIGN(ONE,y(1:nspc))
    y0(1:nspc)=MAX(ABS(y0(1:nspc)),eps)*SIGN(ONE,y0(1:nspc))
    !
    ! ----calc values of Jacobian
    TimeJacobianA=MPI_WTIME()
    CALL JacobiMatrix(BAT,A,Rate,y,Jac_C)
    !CALL JacobiMatrix(BAT,A,Rate,y,Jac_C)
    TimeJac=TimeJac+MPI_WTIME()-TimeJacobianA
    Output%npds=Output%npds+1
    !
    !---- Save matrix structures for matlab spy
    IF (MPI_ID==0.AND.MatrixPrint) CALL SaveMatricies(A,B,BAT,Miter,LU_Miter,BSP)
    !
    !---- calculate a first stepsize based on 2nd deriv.
    CALL InitialStepSize(h,hmin,absh,Jac_C,dummyrate,t,y(1:nspc),RCo%pow)
    !CALL InitialStepSize(h,hmin,absh,Jac_C,Rate,t,y(1:nspc),RCo%pow)
    ! 
    !call printsparse(Jac_C,'*')
    !print*, 'debug :: h0=', h,absh
    !stop 'integrationmod'
    !-----------------------------------------------------------------------
    !-- Initialize Netcdf Output File
    ! 
    tncdf_ind=0
    OutNetcdfANZ=COUNT(OutNetcdfPhase(:)=='g')+2*COUNT(OutNetcdfPhase(:)=='a')
    OutNetcdfDIM=OutNetcdfANZ
    IF ( combustion ) OutNetcdfDIM=OutNetcdfANZ+1
    ALLOCATE(yNcdf(OutNetcdfDIM))          ! output array , +1 for Temperature
    yNcdf(:)=0.0d0
    !
    !-----------------------------------------------------------------------
    !--  Netcdf Output File
    IF (ErrorLog==1) ErrVals=0.0d0
    IF ( ntAqua > 0 ) THEN
      actLWC=pseudoLWC(Tspan(1))
      wetRad=(3.0d0/4.0d0/PI*actLWC/SPEK(1)%Number)**(1.0d0/3.0d0)*1.0d-1
    ELSE
      actLWC=ZERO
      wetRad=ZERO
    END IF
    !
    IF (MPI_ID==0.AND.NetCdfPrint) THEN 
      iStpNetCDF=1
      TimeNetCDF=MPI_WTIME()
      errind(1,1)=0
      CALL InitNetcdf(tAnf, tEnd)
      WRITE(*,*) '  Init NetCDF................... done'
      !--  print initial values to NetCDF file
      CALL SetOutputNCDF( y, yNcdf, Tspan(1), actLWC )
      CALL StepNetCDF   ( Tspan(1)                          &
      &                 , yNcdf(:)                          &
      &                 , tncdf_ind                         &
      &                 , (/ actLWC                         &
      &                    , h                              &
      &                    , SUM(y(1:ntGas))                &
      &                    , SUM(y(ntGas+1:ntGas+ntAqua))   &
      &                    , wetRad                         &
      &                    , ZERO                        /) &
      &                 , errind                            &
      &                 , SUM(y(iDiag_Schwefel))            &
      &                 , ZERO                              )
      TimeNetCDF=MPI_WTIME()-TimeNetCDF
    END IF
    !
    !--- Print the head of the .simul file 
    !CALL PrintHeadSimul(ChemFile)
    !
    tmp_tb=(Tspan(2)-Tspan(1))*1.0d-2
    !
    !
    ! Set PI Parameter
    ! set for normal case
    PI_norm%Kp=0.13d0
    PI_norm%KI=1.0d0/15.0d0        
    PI_norm%ThetaMAX=2.0d0         
    PI_norm%rho=1.2d0 
    ! set for rejected case
    PI_rej%Kp=0.0d0
    PI_rej%KI=1.0d0/5.0d0        
    PI_rej%ThetaMAX=2.0d0        
    PI_rej%rho=1.2d0 
    !===============================================================================
    !=================================THE MAIN LOOP=================================
    !===============================================================================
    TimeIntegrationA=MPI_WTIME()
    IF (MPI_ID==0) WRITE(*,*) '  Start Integration............. '
    IF (MPI_ID==0) WRITE(*,*) ' '
    !
    tnew=0.0d0
    !
    DO 
      !
      absh=MIN(maxStp,MAX(minStp,absh))
      h=absh
      !
      !-- Stretch the step if within 5% of tfinal-t.
      IF (1.05*absh>=Tspan(2)-t) THEN
        h=Tspan(2)-t
        absh=ABS(h)
        done=.TRUE.
      END IF
      DO                                ! Evaluate the formula.
        !
        !-- LOOP FOR ADVANCING ONE STEP.
        failed=.FALSE.                   ! no failed attempts
        !
        ! Rosenbrock Timestep 
        print*, 'DEBUG::INTEGR     temp=y(end)=',y0(nDim)
        CALL Rosenbrock(  y0            &       ! old concentration
        &               , t             &       ! time
        &               , h             &       ! stepsize
        &               , RCo           &       ! Rosenbrock parameter
        &               , error         &       ! error value
        &               , errind        &       ! max error component
        &               , y             )       ! new concentration 
        !
        tnew=t+h
        !
        IF (done) THEN
          tnew=Tspan(2)         ! Hit end point exactly.
          h=tnew-t                        ! Purify h.
        END IF
        Output%ndecomps=Output%ndecomps+1
        Output%nRateEvals=Output%nRateEvals+RCo%nStage
        Output%nSolves=Output%nSolves+RCo%nStage
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
          Output%nfailed=Output%nfailed+1
          IF (absh<=hmin) THEN
            WRITE(*,*) ' Stepsize to small: ', h
            CALL FinishMPI()
            STOP '....Integration_Mod '
          END IF
          !
          IF (PI_StepSize .AND. Output%nSteps>1) THEN
            CALL PI_StepsizeControl(absh,RtolRow,error,errorOld,h,hOld,PI_rej,RCo)
          ELSE
            absh=MAX(hmin,absh*MAX(0.1d0,0.8d0*(1.0d0/error)**RCo%pow))
            !absh=MIN(Args%hmax,absh)
          END IF
          h=absh
          done=.FALSE.
        ELSE                            !succ. step
          EXIT
        END IF
      END DO
      Output%nsteps=Output%nsteps+1
      !
      tmp_ta=t-Tspan(1)
      !
      !SPEK(1)%wetRadius=(3.0d0/4.0d0/PI*actLWC/SPEK(1)%Number)**(1.0d0/3.0d0)
      !print*, t, y(Positionspeciesall('OH')), SPEK(1)%wetRadius
      !
      IF ( (tmp_ta >= StpNetCDF*iStpNetCDF) )  THEN
        iStpNetCDF=iStpNetCDF+1
        IF ( (MPI_ID==0) .AND. (NetCdfPrint) ) THEN !.AND.t>10.0d0) THEN
          zen=Zenith(t)
          IF ( ntAqua > 0 ) THEN
            actLWC=pseudoLWC(t)
            wetRad=(3.0d0/4.0d0/PI*actLWC/SPEK(1)%Number)**(1.0d0/3.0d0)*1.0d-1
          END IF
          ! save data in NetCDF File
          TimeNetCDFA=MPI_WTIME()
          CALL SetOutputNCDF(  y,    yNcdf, t ,  actLWC)
          !
          CALL StepNetCDF   ( t                                &
          &                 , yNcdf                            &
          &                 , tncdf_ind                        &
          &                 , (/ actLWC                        &
          &                    , h                             &
          &                    , SUM(y(1:ntGas))               &
          &                    , SUM(y(ntGas+1:ntGas+ntAqua))  &
          &                    , wetRad                        &
          &                    , zen                        /) &
          &                 , errind                           &
          &                 , SUM(y(iDiag_Schwefel))           &
          &                 , error                            )
          !
          TimeNetCDFA=MPI_WTIME()-TimeNetCDFA
          TimeNetCDF=TimeNetCDF+TimeNetCDFA
        END IF
      END IF
      !
      !-- Call progress bar.
      IF (MPI_ID==0.AND.Ladebalken==1.AND.tmp_ta>=iBar*tmp_tb) THEN
        iBar=iBar+1
        CALL Progress(iBar)
      END IF
      !
      !-- Termination condition for the main loop.
      IF (done) EXIT
      !
      !-- If there were no failures compute a new h.
      IF (PI_StepSize) THEN
        CALL PI_StepsizeControl(absh,RtolRow,error,errorOld,h,hOld,PI_norm,RCo)
      ELSE
        tmp=1.25d0*(error)**RCo%pow
        IF (2.0d0*tmp>1.0d0) THEN
          absh=absh/tmp
        ELSE
          absh=2.0d0*absh
        END IF
      END IF
      !
      !-- Advance the integration one step.
      t=tnew
      y0=y
      !
      !-- for PI stepsize control
      errorOld=error
      hOld=h
    END DO  ! MAIN LOOP
    TimeIntegrationE=MPI_WTIME()
    !
    IF (OrderingStrategie<8) THEN
      mumps_par%JOB=-2
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
