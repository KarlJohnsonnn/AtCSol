!=======================================
!=======================================
! vorlaeufiges Hauptprogramm fuer Tests
! box model chemical kinetics simulation
!=======================================
!=======================================
!
!
PROGRAM AtCSol
  !
  USE Kind_Mod,         ONLY: dp
  USE IO_Mod,           ONLY: Logo, Print_Run_Param, Output_Statistics
  USE InitRoutines_Mod, ONLY: Initialize
  USE Integration_Mod,  ONLY: Integrate

  USE Reac_Mod,     ONLY: InitValAct
  USE Control_Mod,  ONLY: Simulation,               &
                          tBegin, tEnd, dt_output,  &
                          Timer_Finish,  &
                          ICNTL, RCNTL
  !

  REAL(dp) :: DT = 0.0
  REAL(dp) :: T0 = 0.0, TNEXT = 0.0
  REAL(dp) :: H0 = 1.e-7_dp
  REAL(dp), ALLOCATABLE :: Y(:)

  REAL(dp) :: Timer0

  !
  !=======================================================================
  !===                            MAIN Programm
  !=======================================================================

  CALL CPU_TIME(Timer0)

  !-----------------------------------------------------------------------
  ! --- Print the program logo  
  CALL Logo
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !--- Read all input files and build all symbolic structures
  CALL Initialize
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !--- print input parameter, method, tols, concs, temp, pressure,....
  CALL Print_Run_Param
  !-----------------------------------------------------------------------

  IF (Simulation) THEN
    
    Y = InitValAct

    T0 = tBegin
    DT = dt_output
    
    IF (DT < 0.0_dp) THEN
      TNEXT = tEnd
    ELSE 
      TNEXT = T0 + DT
    END IF

    H0 = MIN(H0, TNEXT-T0)


    ! --- start integration loop
    INTEGRATION_LOOP: DO
      
      CALL Integrate ( Y            &  ! current state vector
      &              , H0           &  ! stepsize
      &              , [T0, TNEXT]  &  ! integration invervall
      &              , ICNTL        &  ! integer controle units
      &              , RCNTL        )  ! real value controle units
      
      IF (TNEXT == tEnd) EXIT INTEGRATION_LOOP

      ! --- advance one delta t
      T0    = TNEXT
      TNEXT = TNEXT + DT
      
      ! --- hit end point exactly
      IF ( TNEXT > tEnd ) TNEXT = tEnd

    END DO INTEGRATION_LOOP

    !-----------------------------------------------------------------------
    ! --- stop timer and print output statistics
    CALL CPU_TIME(Timer_Finish)
    Timer_Finish = Timer_Finish - Timer0
    CALL Output_Statistics
    !-----------------------------------------------------------------------

  END IF


END PROGRAM AtCSol
