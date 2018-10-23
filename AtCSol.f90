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
  USE Kind_Mod
  USE IO_Mod,           ONLY: Logo, Print_Run_Param, Output_Statistics
  USE InitRoutines_Mod, ONLY: Initialize
  USE Integration_Mod,  ONLY: Integrate

  USE Reac_Mod,     ONLY: InitValAct, Temperature0, Combustion
  USE Control_Mod,  ONLY: Simulation,               &
                          tBegin, tEnd, StpNetCDF,  &
                          Timer_Finish,             &
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
    DT = StpNetCDF
    
    IF (DT < 0.0_dp) THEN
      TNEXT = tEnd
    ELSE 
      TNEXT = StpNetCDF
    END IF

    H0 = MIN(H0, TNEXT-T0)


    !-----------------------------------------------------------------------
    ! --- Start the integration routine 
    !-----------------------------------------------------------------------

    INTEGRATION_LOOP: DO
      
      CALL Integrate ( Y            &  ! initial concentrations activ species
      &              , H0           &  ! stepsize
      &              , [T0, TNEXT]  &  ! integration invervall
      &              , ICNTL        &  ! integer controle units
      &              , RCNTL        )  ! real value controle units
      
      IF (TNEXT == tEnd) EXIT INTEGRATION_LOOP
      T0    = TNEXT
      TNEXT = TNEXT + DT
      
      ! --- Hit end point exactly.
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
