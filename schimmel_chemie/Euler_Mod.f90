!=========================================================================!
!                                                                         !         
!                   Module for time integration with diagonal             !
!                      implicite Rosenbrock-Wanner-Methods                !         
!                                                                         !         
!=========================================================================!
!
MODULE Euler_Mod
  !
  USE Kind_Mod
  USE mo_MPI
  USE mo_control
  USE mo_IO
  USE Sparse_Mod
  USE Sparse2_Mod
  USE Chemsys_Mod
  USE Rosenbrock_Mod
  USE Rates_Mod
  !USE Factorisation_Mod
  USE Integration_Mod
  IMPLICIT NONE
  
  !
  CONTAINS
  !
  !======================================================================= 
  !===================    Time Integration Routine  ======================
  !======================================================================= 
  SUBROUTINE IntegrateEuler(y0,Tspan,h)
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
    !-------------------------------------------------------------------
    ! Output:
    !   - Output.......... struct Out (above) contains output statistics
    !   - Unit............ .simul data 
    INTEGER :: Unit=90
    !-------------------------------------------------------------------
    ! Temporary variables:
    TYPE(CSR_Matrix_T) :: Jac!, Id
    REAL(RealKind) :: t,tnew            ! current time, tnew
    REAL(RealKind) :: y(nspc)        ! current y vector
    REAL(RealKind) :: hy0(nspc)
    REAL(RealKind) :: Rate(neq)     ! rate vector
    REAL(RealKind) :: h
    REAL(RealKind) :: tmp_ta,tmp_tb
    ! 
    INTEGER :: i=0,ii=0
    INTEGER :: iOutAqua,iOUt
    REAL(RealKind) :: yOut(Nspc)
    !
    LOGICAL :: done=.FALSE.
    !
    !INTEGER :: mUnit            !Mapping Unit
    !CHARACTER(50) :: mName      !Mapping file name
    !
    !
    !
    INTERFACE
      REAL FUNCTION pseudoLWC(zeit)
        USE Kind_Mod
        USE mo_reac
        USE mo_control
        REAL(RealKind) ::zeit
      END FUNCTION pseudoLWC
    END INTERFACE
    ! 
    !
    TimeSymbolic=MPI_WTIME()                    ! start timer for symb phase
    !---- transpose (B-A) matrix from chemsys
    CALL SymbolicAdd(B,A,BA)
    CALL SparseAdd(B,A,BA,'-')
    CALL TransposeSparse(BA,BAT) 
    !
    !---- Initialize stuff for Rosenbrock
    t=Tspan(1)
    y=y0
    hy0=y0/h
    !
    !
    !---- locate nonzeros in Jacobian matrix
    CALL SymbolicMult(BAT,A,Jac)
    !
    !---- calc rate vector
    CALL Rates(t,y,Rate)
    Output%nRateEvals=Output%nRateEvals+1
    !
    ! ----calc values of Jacobian
    TimeJacobianA=MPI_WTIME()
    CALL JacobiMatrix(BAT,A,Rate,y,Jac)
    TimeJacobianE=MPI_WTIME()
    TimeJac=TimeJac+TimeJacobianE-TimeJacobianA
    Output%npds=Output%npds+1
    !
    !---- Set nonzero entries in coef Matrix
    CALL SetMiterStructure(Miter,'cl',A,BAT,Jac,1.0d0)
    !
    !---- Choose ordering/factorisation strategie and do symb LU fact
    CALL SetLU_MiterStructure(LU_Miter,Miter,'cl')
    !
    !
    TimeSymbolic=MPI_WTIME()-TimeSymbolic       ! end timer for symb phase
    ! 
    !--- Print the head of the .simul file 
    CALL PrintHeadSimul(ChemFile)
    !
    !
    tmp_tb=Tspan(2)-Tspan(1)
    !
    !===============================================================================
    !=================================THE MAIN LOOP=================================
    !===============================================================================
    TimeIntegrationA=MPI_WTIME()
    !
    DO 
      !
      !-- Stretch the step if within 10% of tfinal-t.
      IF (1.1d0*h>=Tspan(2)-t) THEN
        h=Tspan(2)-t
        h=ABS(h)
        done=.TRUE.
      END IF
      !
      !
      CALL BackwardEuler(y,y0,t,h)
      !
      tnew=t+h
      IF (done) THEN
        tnew=Tspan(2)         ! Hit end point exactly.
        h=tnew-t                        ! Purify h.
      END IF
      Output%ndecomps=Output%ndecomps+1
      Output%nRateEvals=Output%nRateEvals+1
      Output%nSolves=Output%nSolves+1
      !
      !
      Output%nsteps=Output%nsteps+1
      !
      tmp_ta=tnew-Tspan(1)
      iOutAqua=SIZE(ListGas2)
      yOut=y
      DO iOut=iOutAqua+1,SIZE(ListAqua2)+iOutAqua
        yOut(iOut)=yOut(iOut)/LWCconst/mol2part
      END DO
      !
      IF (tmp_ta>=ii*0.01*tmp_tb.AND.NetCdfPrint.EQV..TRUE..AND.MPI_ID==0) THEN
        !-- Save data in .simul
        !-- Call progress bar.
        IF (tmp_ta>=i*0.1*tmp_tb.AND.Ladebalken==1) THEN
          CALL Progress(i)
          i=i+1
        END IF
        ii=ii+1
      END IF
      !
      !-- Termination condition for the main loop.
      IF (done) EXIT
      !
      !-- Advance the integration one step.
      t=tnew
      y0=y
      !IF (Output%nSteps==5) STOP 'hlaudhfkauhs'
    END DO  ! MAIN LOOP
    TimeIntegrationE=MPI_WTIME()
    !
    CALL CloseSimulFile(Unit)
    ALLOCATE(Output%y(nspc))
    Output%y=y
  END SUBROUTINE IntegrateEuler
END MODULE Euler_Mod
