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
  USE mo_reac, ONLY: y_total
  USE Sparse_Mod
  USE Sparse2_Mod
  USE Chemsys_Mod
  USE Rosenbrock_Mod
  USE Rates_Mod
  USE Factorisation_Mod
  USE ErrorROW_Mod
  USE NetCDF_Mod
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
  SUBROUTINE Integrate(y0,Tspan,Atol,RtolRow,vers,method)
    !--------------------------------------------------------------------
    ! Input:
    !   - y0 ............. Initial vector
    !   - Tspan .......... (/ SimulationTimeStart , SimulationTimeEnd /)
    !   - Atol ........ abs. tolerance for gas spc
    !   - RtolRow ........ rel. tolerance for Rosenbrock-Wanner-Method
    !   - vers ........... version of solving the linear system
    !   - method...........Rosenbrock-Wanner-Method
    !   - PrintSpc ....... print spc PrintSpc(1:3)
    REAL(RealKind) :: y0(:)
    REAL(RealKind) :: Tspan(2)
    REAL(RealKind) :: Atol(2)
    REAL(RealKind) :: RtolROW
    CHARACTER(2) :: vers
    CHARACTER(*) :: method
    !-------------------------------------------------------------------
    ! Output:
    !   - Output.......... struct Out (above) contains output statistics
    !   - Unit............ .simul data 
    TYPE(RosenbrockMethod_T) :: RCo  
    !INTEGER :: Unit=90
    !-------------------------------------------------------------------
    ! Temporary variables:
    TYPE(CSR_Matrix_T) :: Jac!, Id
    REAL(RealKind) :: t            ! current time
    REAL(RealKind) :: y(nspc)        ! current y vector
    REAL(RealKind) :: hy0(nspc)
    REAL(RealKind) :: Rate(neq)     ! rate vector
    REAL(RealKind) :: h, hmin, absh, tnew, temp, hOld
    REAL(RealKind) :: error, errorOld
    REAL(RealKind) :: tmp_ta,tmp_tb
    REAL(RealKind) :: actLWC, zen
    INTEGER        :: errind(1,1)
    REAL(RealKind) :: ErrVals(nspc)
    ! 
    INTEGER :: i=-1
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
    !do i=1,size(y0)
    !  print*, y_name(i),y0(i)
    !end do
    !
    !
    tncdf_ind=0
    !
    TimeSymbolic=MPI_WTIME()       ! start timer for symb phase
    !
    y_total=SUM(y0)
    IF (MPI_ID==0) THEN
      WRITE(*,*) '  Initial values:   '
      WRITE(*,*) 
      WRITE(*,'(A35,2X,E23.14,A13)')  '      sum initval (gaseous)    = ', SUM(y0(1:ntGas)), '  [molec/cm3] '
      WRITE(*,'(A35,2X,E23.14,A13)')  '      sum initval (aqueous)    = ', SUM(y0(ntGas+1:)), '  [molec/cm3] '
      WRITE(*,'(A35,2X,E23.14,A13)')  '      sum emissions (gaseous)  = ', SUM(y_e), '  [molec/cm3] '
      WRITE(*,*) 
    END IF
    !---- transpose (B-A) matrix from chemsys
    CALL SymbolicAdd(B,A,BA)
    CALL SparseAdd(B,A,BA,'-')
    CALL TransposeSparse(BA,BAT) 
    !
    !call printsparse(A,'*')
    !stop
    !  
    method=ADJUSTL(TRIM(method))
    CALL SetRosenbrockMethod(RCo,method)  
    !
    !CALL GetConstReactionRates(ReactionRateConst,y0,Tspan(1),RefTemp)
    !---- Initialize stuff for Rosenbrock
    CALL SetRosenbrockArgs(Rate,y0,Tspan,Atol,RtolROW,RCo%pow)
    t=Tspan(1)
    y=y0
    !
    !---- locate nonzeros in Jacobian matrix
    CALL SymbolicMult(BAT,A,Jac)
    !
    ! ----calc values of Jacobian

    CALL Rates(t,y,Rate)
    Rate(:)=MAX(ABS(Rate(:)),eps)*SIGN(1.0d0,Rate(:))
    y(:)=MAX(ABS(y(:)),eps)*SIGN(1.0d0,y(:))

    !print*, 'nID=', MPI_ID, 'SumJacColInd=',SUM(Jac%ColInd), 'sumf0=',sum(abs(Rate)),&
    !  &           'sumY=',SUM(ABS(y)), 'sumAvls=', SUM(ABS(A%val)), 'sumBATvals=',SUM(ABS(BAT%val))

    TimeJacobianA=MPI_WTIME()
    !CALL JacobiMatrix(BAT,A,Args%f0,y,Jac)
    CALL JacobiMatrix(BAT,A,Rate,y,Jac)

    TimeJacobianE=MPI_WTIME()
    TimeJac=TimeJac+TimeJacobianE-TimeJacobianA
    Output%npds=Output%npds+1
    !
    !---- Set nonzero entries in coef Matrix
    CALL SetMiterStructure(Miter,vers,A,BAT,Jac,RCo%ga)
    !
    !---- Choose ordering/factorisation strategie and do symb LU fact
    CALL SetLU_MiterStructure(LU_Miter,Miter,vers)
    !
    IF (MPI_ID==0) WRITE(*,*) '  Symbolic phase................ done'
    !---- Save matrix structures for matlab spy
    IF (MPI_ID==0.AND.MatrixPrint) CALL SaveMatricies(A,B,BAT,Miter,LU_Miter,BSP)
    !
    CALL InitialStepSize(h,hmin,absh,Jac,t,y)
    hy0=y0/h
    !---- stop timer for symbolic matrix calculations
    TimeSymbolic=MPI_WTIME()-TimeSymbolic     
    ! 
    !print*, 'Start h= ', h
    !
    sqrtNSPC=SQRT(REAL(nspc,KIND=RealKind))
    !-----------------------------------------------------------------------
    !-- Initialize Netcdf Output File
    OutNetcdfANZ=COUNT(OutNetcdfPhase(:)=='g')+2*COUNT(OutNetcdfPhase(:)=='a')
    ! 
    ALLOCATE(yNcdf(OutNetcdfANZ))                        ! output array
    yNcdf(:)=0.0d0
    !
    actLWC=pseudoLWC(Tspan(1))                  ! in [l/m3]
    !-----------------------------------------------------------------------
    !--  Netcdf Output File
    IF (ErrorLog==1) ErrVals=0.0d0
    !
    IF (MPI_ID==0.AND.NetCdfPrint) THEN 
      iStpNetCDF=1
      TimeNetCDF=MPI_WTIME()
      errind(1,1)=0
      CALL InitNetcdf(tAnf, tEnd)
      WRITE(*,*) '  Init NetCDF................... done'
      !--  print initial values to NetCDF file
      CALL SetOutputNCDF( y, yNcdf, Tspan(1), actLWC )
      WRITE(*,*) '  Init Values...................     ',yNcdf(1:3)
      CALL StepNetCDF   ( Tspan(1)                          &
      &                 , yNcdf                             &
      &                 , tncdf_ind                         &
      &                 , (/ actLWC                         &
      &                    , h                              &
      &                    , SUM(y(1:ntGas))                &
      &                    , SUM(y(ntGas+1:ntGas+ntAqua))   &
      &                    , SPEK(1)%wetRadius              &
      &                    , 0.0d0                       /) &
      &                 , errind                            &
      &                 , SUM(y(iDiag_Schwefel))            &
      &                 , 0.0d0                           )
      TimeNetCDF=MPI_WTIME()-TimeNetCDF
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,MPIErr)
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
      absh=MIN(Args%hmax,MAX(hmin,absh))
      h=Args%Tdir*absh
      !tnew=t+h
      !print*, 'h=',h,'t=',t
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
        !print*, 'ID= ', MPI_ID,'h=',h,'t=',t,'tnew=',tnew
        SELECT CASE (vers)
          CASE ('cl')
            CALL ros_classic(y,error,errind,y0,t,h,RCo,errVals)
            Output%npds=Output%npds+1
          CASE ('ex')
            CALL ros_extended(y,error,errind,y0,t,h,RCo,errVals)
        END SELECT
        !
        tnew=t+h
        !print*, '---------'
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
          failed = (error > 1.0d0)
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
          h=Args%Tdir*absh
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
          actLWC = pseudoLWC(t)
          zen    = Zenith(t)
          SPEK(1)%wetRadius=(3.0d0/4.0d0/PI*actLWC/SPEK(1)%Number)**(1.0d0/3.0d0)*1.0d-1
          ! save data in NetCDF File
          TimeNetCDFA=MPI_WTIME()
          CALL SetOutputNCDF(  y,    yNcdf, t ,  actLWC)
          !
          CALL StepNetCDF   ( t                                  &
          &                 , yNcdf                              &
          &                 , tncdf_ind                          &
          &                 , (/ actLWC                          &
          &                    , h                               &
          &                    , SUM(y(1:ntGas))                 &
          &                    , SUM(y(ntGas+1:ntGas+ntAqua))    &
          &                    , SPEK(1)%wetRadius               &
          &                    , zen                             /) &
          &                    , errind                          &
          &                    , SUM(y(iDiag_Schwefel))          &
          &                    , error                         )
          !
          TimeNetCDFA=MPI_WTIME()-TimeNetCDFA
          TimeNetCDF=TimeNetCDF+TimeNetCDFA
        END IF
      END IF
      !CALL MPI_BARRIER(MPI_COMM_WORLD,MPIErr)
      !
      !print*,'ID=',MPI_ID, 'h=', h,'t=',t
      !-- Call progress bar.
      IF (MPI_ID==0.AND.Ladebalken==1.AND.tmp_ta>=i*tmp_tb) THEN
        i=i+1
        CALL Progress(i)
      END IF
      !
      !-- Termination condition for the main loop.
      IF (done) EXIT
      !
      !-- If there were no failures compute a new h.
      !IF (failed) THEN
        IF (PI_StepSize) THEN
          CALL PI_StepsizeControl(absh,RtolRow,error,errorOld,h,hOld,PI_norm,RCo)
        ELSE
          temp=1.25d0*(error)**RCo%pow
          !IF (5.0d0*temp>1.0d0) THEN
          IF (2.0d0*temp>1.0d0) THEN
            absh=absh/temp
          ELSE
            !absh=5.0d0*absh
            absh=2.0d0*absh
          END IF
          !absh=MIN(Args%hmax,absh)
        END IF
      !END IF
      !
      !-- Advance the integration one step.
      !CALL MPI_BARRIER(MPI_COMM_WORLD,MPIErr)
      t=tnew
      y0=y
      !
      !-- for PI stepsize control
      errorOld=error
      hOld=h
    END DO  ! MAIN LOOP
    TimeIntegrationE=MPI_WTIME()
    !
    IF (MPI_ID==0) THEN
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) '  Nonzeros in matrices:'
      WRITE(*,*) ' '
      WRITE(*,*) '             alpha  = ',SIZE(A%Val)
      WRITE(*,*) '     + (beta-alpha) = ',SIZE(BA%Val)
      WRITE(*,*) '     +  L+U factors = ',SIZE(LU_Miter%Val) 
      WRITE(*,*) '     -------------------------------------'
      WRITE(*,*) '     =     sum(NNZ) = ',SIZE(A%Val)+SIZE(BA%Val)+SIZE(LU_Miter%Val)
      WRITE(*,*) ' '
    END IF

    !
    !CALL printsparse(Miter,'Miter_INORGex')
    !stop
    call MPI_Barrier( MPI_COMM_WORLD, i_error)
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
  !
  SUBROUTINE SetMiterStructure(Miter,vers,A,BAT,Jac,g)
    !-----------------------------------------------------------
    ! Input: 
    !   - version (extended or classic)
    CHARACTER(2), INTENT(IN) :: vers
    REAL(RealKind), INTENT(IN) :: g
    ! Input (optional):
    !   - Matrix A (left side of chemical system)
    !   - Matrix BAT = (B - A)^T 
    !   - Jacobian 
    TYPE(CSR_Matrix_T), INTENT(IN) :: A, BAT, Jac
    !-----------------------------------------------------------
    ! Output:
    !   - Miter with nonzero structure 
    TYPE(CSR_Matrix_T), INTENT(OUT) :: Miter
    ! Temporary variables:
    TYPE(CSR_Matrix_T) :: Id
    !
    !
    SELECT CASE (vers)
      CASE ('cl')
        CALL SparseID(Id,nspc)
        CALL SymbolicAdd(Id,Jac,Miter)
        CALL Kill_Matrix_CSR(Id)
        !CALL Kill_Matrix_CSR(Jac)
      CASE ('ex')
        IF (ckTEMP) THEN
          CALL SymbolicExtendedMatrixTemp(A,BAT,Miter)
          CALL Miter0_ExtendedTemp(BAT,A,g,Miter)
        ELSE
          CALL SymbolicExtendedMatrix(A,BAT,Miter)
          CALL Miter0_Extended(BAT,A,g,Miter)
        END IF
      CASE DEFAULT
        STOP 'change solveLA in .run to "ex" or "cl"'
    END SELECT
  END SUBROUTINE SetMiterStructure
  !
  !
  SUBROUTINE SetLU_MiterStructure(LU_Miter,Miter,version)
    TYPE(CSR_Matrix_T), INTENT(OUT) :: LU_Miter
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: Miter
    CHARACTER(2) :: version
    INTEGER, ALLOCATABLE :: InvPermu(:)
    !
    TYPE(SpRowIndColInd_T) :: MiterFact
    TYPE(SpRowColD_T) :: MiterMarko
    INTEGER :: i
    !
    IF (OrderingStrategie<8.OR.ParOrdering>=0) THEN 
      ! Use MUMPS to factorise and solve     
      ! Convert compressed row format to row index format for MUMPS
      CALL CompRowToIndRow(Miter,MiterFact)
      CALL InitMumps(MiterFact) 
      CALL CSR_Matrix_To_SpRowColD(MiterMarko,Miter) 
      CALL PermuToInvPer(InvPermu,Mumps_Par%SYM_PERM)
      CALL SymbLU_SpRowColD(MiterMarko,InvPermu)        
      CALL RowColD_To_CRF_Matrix(LU_Miter,MiterMarko)
      
      ! set values to LU_Miter and get the LU permutation vector
      CALL SetLU_Val_INI(LU_Miter,Miter,LU_Perm)
      CALL Kill_Matrix_SpRowIndColInd(MiterFact)
    ELSE
      ! Permutation given by Markowitz Ordering strategie
      CALL CSR_Matrix_To_SpRowColD(MiterMarko,Miter) 
      CALL SymbLUMarko_SpRowColD(MiterMarko)        
      CALL RowColD_To_CRF_Matrix(LU_Miter,MiterMarko)
      
      !STOP 'Integration_Mod'
      ! set values to LU_Miter and get the LU permutation vector
      CALL SetLU_Val_INI(LU_Miter,Miter,LU_Perm)
        !call printsparse(LU_Miter,'*')
        !print*, 'diagptr=', LU_Miter%DiagPtr
        !stop 'integration'
      ! 
    END IF
    IF (version=='ex') THEN
      !          _                           _ 
      !         |              |   g*A_Mat    |
      !         |--------------+--------------|     ~=~   ValCopy = const.
      !         |_  BAT_Mat    |             _|
      !
      ALLOCATE(ValCopy(LU_Miter%RowPtr(LU_Miter%m+1)-1))
      ValCopy=LU_Miter%Val
      ! 
      CALL Kill_Matrix_SpRowColD(MiterMarko)
    END IF
  END SUBROUTINE SetLU_MiterStructure
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
