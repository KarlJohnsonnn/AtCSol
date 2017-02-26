!=======================================
!=======================================
! vorlaeufiges Hauptprogramm fuer Tests
!=======================================
!=======================================
!
!
PROGRAM Main_ChemKin
  !
  USE Sparse_Mod
  USE Kind_Mod
  USE Chemsys_Mod
  USE Integration_Mod
  USE Rosenbrock_Mod
  USE mo_control
  USE mo_unirnk
  USE mo_reac, ONLY: nspc,combustion, MW, rMW
  USE mo_MPI
  USE mo_IO
  USE mo_ckinput
  !USE Euler_Mod
  USE NetCDF_Mod
  IMPLICIT NONE
  !
  CHARACTER(80) :: Filename0=''        ! *.run file
  !
  REAL(RealKind) :: Atol(2)
  INTEGER :: i,nSchwefel
  CHARACTER(80) :: constH=''
  CHARACTER(80) :: neuTolR=''
  CHARACTER(80) :: neuTolA=''
  CHARACTER(80) :: neuROW=''
  CHARACTER(16) :: tmpchar0
  REAL(RealKind) :: h 
  REAL(RealKind), PARAMETER :: HR=3600.0d0
  REAL(RealKind) :: tmpMW0
  INTEGER :: i_error, linc
  !
  !
  !================================================================
  !===                     MAIN Programm
  !================================================================
  !
  !----------------------------------------------------------------
  ! ---  Initialize MPI library 
  CALL StartMPI()
  !
  !----------------------------------------------------------------
  ! --- Print the program logo  
  CALL Logo()
  !
  !----------------------------------------------------------------
  ! --- Read run control parameters (which runfile)
  CALL getarg(1,FileName0)             
  IF (FileName0=='')   THEN
    WRITE(*,*) 'Input RUNFilename:'
    READ(*,*)   FileName0
  END IF
  FileName0=TRIM(ADJUSTL(FileName0))
  !
  !
  !
  !================================================================
  !===                     Initialization
  !================================================================
  !
  !----------------------------------------------------------------
  ! --- Read run-file
  CALL InitRun(FileName0)
  IF (MPI_ID==0) WRITE(*,*) '  Initialize run-file .......... done'
  !
  !----------------------------------------------------------------
  ! --- Initialize all reaction types
  CALL InitNReacType()
  !
  Tspan=(/tAnf , tEnd/)
  !
  !----------------------------------------------------------------
  ! --- set cloud intervall
  LWCBounds(1)=0.0d0*HR
  LWCBounds(2)=LWCBounds(1) + 1.d0*HR
  LWCBounds(3)=LWCBounds(2) + 0.25d0*HR
  LWCBounds(4)=LWCBounds(3) + 9.5d0*HR
  LWCBounds(5)=LWCBounds(4) + 0.25d0*HR
  LWCBounds(6)=LWCBounds(5) + 1.d0*HR
  !----------------------------------------------------------------
!  --- read the .sys data, save coefs in sparse matrix
  OPEN(UNIT=89,FILE=ADJUSTL(TRIM(ChemFile))//'.chem',STATUS='UNKNOWN')
  IF (MPI_ID==0) WRITE(*,'(A38)',ADVANCE='NO') '  Read sys-file ................ done'
  IF (ChemFile(1:2)=='ck') THEN
    ! *)
    IF (MPI_ID==0) WRITE(*,*) ' ---->  Solve Gas Energy Equation '
    combustion=.TRUE.
    Time_Read=MPI_WTIME()
    CALL Read_Elements(ChemFile,969)
    CALL Read_Species(ChemFile,969)
    CALL Read_Reaction(ChemFile,969)
    CALL Read_ThermoData(SwitchTemp,DataFile,696,nspc)
    print*, 'djfhalskdjhfaisdh ', nspc,TRIM(ChemFile)//'.mw'
    OPEN(UNIT=998,FILE=TRIM(ChemFile)//'.mw',STATUS='UNKNOWN')
    ALLOCATE(MW(nspc),rMW(nspc))
    MW=ZERO
    DO i=1,nspc
      READ(998,*) tmpChar0, tmpMW0
      MW(PositionSpeciesAll(tmpChar0))=REAL(tmpMW0,RealKind)
    END DO
    rMW(:)=ONE/MW(:)
    DO i=1,nspc
      print*, 'debug:: main    ',i, MW(i),rMW(i)
    END DO
    CLOSE(linc)
    !
    CALL PrintHeadSpecies(ChemFile,89)
    CALL PrintSpecies(ListGas2,89)
    CALL PrintHeadReactions(89)
    CALL PrintReactions(ReactionSystem,89,.TRUE.)       ! .TRUE. in 3rd agument for chemkin input file
    CALL PrintFinalReactions(89)
    !
    ALLOCATE(InitValAct(ntGas),y_e(ntGas))
    ALLOCATE(InitValKat(ntKat))
    CALL Read_GASini(InitFile,InitValAct,InitValKat)
    CALL Read_EMISS(InitFile,y_e)
    CALL GetSpeciesNames(ChemFile,y_name)
  ELSE
    IF (MPI_ID==0) WRITE(*,*) ' ---->  Fix Temperature'
    Time_Read=MPI_WTIME()
    CALL ReadSystem(ChemFile)
    !CALL InputDatFile(DataFile)
    !   
    !----------------------------------------------------------------
    ! ---  build the coeficient matrices and write .chem
    CALL PrintHeadSpecies(ChemFile,89)
    CALL PrintSpecies(ListGas2,89)
    CALL PrintSpecies(ListAqua2,89)
    CALL PrintSpecies(ListSolid2,89)
    CALL PrintSpecies(ListPartic2,89)
    CALL PrintSpecies(ListNonReac2,89)
    CALL PrintHeadReactions(89)
    ! 
    !----------------------------------------------------------------
    ! --- Build the reaction system
    CALL AllListsToArray( ReactionSystem           &
    &                    ,ListRGas                 &
    &                    ,ListRHenry               &
    &                    ,ListRAqua                &
    &                    ,ListRDiss                )
    !
    !CALL CheckConstants(ReactionSystem)
    !----------------------------------------------------------------
    ! --- print reactions and build A, B and (B-A) structure
    CALL PrintReactions(ReactionSystem,89)
    CALL PrintFinalReactions(89)
    !
    !----------------------------------------------------------------
    ! --- Input of initial data and thermodynamic properties
    CALL InputChemicalData(InitFile,DataFile,MetFile)
    ! 
  END IF
  CLOSE(89)
  ! Stop Input timer
  Time_Read=MPI_WTIME()-Time_Read
  IF (MPI_ID==0) WRITE(*,*) '  Read ini-file ................ done'
  IF (MPI_ID==0) WRITE(*,*) '  Print chem-file .............. done'
  !
  !----------------------------------------------------------------
  ! --- Split the rate array in mpi_np parts 
  CALL BuildPartitions(MyParties,neq,MPI_np)
  !
  !
  !----------------------------------------------------------------
  ! ---  print some species (set output spc in InitFile under DIAG)
  CALL Read_Diag(OutNetcdfSpc,OutNetcdfPhase,InitFile)
  !
  !
  !-----------------------------------------------------------------------
  ! --- Start timer and set absolut tolerance for species
  Timer_Start=MPI_WTIME()
  Atol=(/AtolGas,AtolAqua/)
  ! 
  NSactNR=nspc+neq
  IF ( combustion ) THEN
    nDIM=nspc+1
    nDIMcl=nspc+1
    nDIMex=nspc+neq+1
  ELSE
    nDIM=nspc
    nDIMcl=nspc
    nDIMex=nspc+neq
  END IF
    !
  ! 
  CALL getarg(2,neuTolR)             
  CALL getarg(3,neuTolA)             
  IF (neuTolR/='')   THEN
    READ(neuTolR,*) RtolROW
  END IF
  IF (neuTolA/='')   THEN
    READ(neuTolA,*) AtolGas
    AtolAqua=AtolGas
  END IF
  !
  CALL getarg(4,neuROW)
  IF (neuROW/='')   THEN
    READ(neuRow,*) RosenbrockMethod
    RosenbrockMethod=ADJUSTL('METHODS/'//RosenbrockMethod)
  END IF
  !
  ! hole alle schwefelstoffe
  nSchwefel=0
  DO i=1,nspc
    IF (INDEX(y_name(i),'S')>0) THEN
      nSchwefel=nSchwefel+1
    END IF
  END DO
  ALLOCATE(iDiag_Schwefel(nSchwefel))
  nSchwefel=0
  DO i=1,nspc
    IF (INDEX(y_name(i),'S')>0) THEN
      nSchwefel=nSchwefel+1
      iDiag_Schwefel(nSchwefel)=PositionSpeciesAll(y_name(i))
    END IF
  END DO
  !
    !print*, 'debug :: nach schwefel'
    !do i=1,nspc
      !print*, y_name(i),InitValAct(i)
      !WRITE(111,*) y_name(i)
    !end do
    !stop
  !
  CALL Print_Run_Param()
  !-----------------------------------------------------------------------
  ! --- choose between implicit Euler or some Rosenbrock-Wanner methode
  IF (ImpEuler==1) THEN
    !---  Read run control parameters (constant stepsize)
    !  
    !CALL getarg(2,constH)             
    !IF (constH=='')   THEN
    !  WRITE(*,*) 'Input constant stepsize:'
    !  READ(*,*)   constH
    !END IF
    !!
    !READ(ConstH,*) h
    !CALL IntegrateEuler(y_iconc(1:nspc),Tspan,h)
    IF (MPI_ID==0) THEN
      WRITE(*,*) ' Currently not available....'
    END IF
    CALL Dropout()
  ELSE
    CALL Integrate (  InitValAct(1:nspc)   &     ! initial concentrations activ species
    &               , Tspan                &     ! integration invervall
    &               , Atol                 &     ! abs. tolerance gas species
    &               , RtolROW              &     ! rel. tolerance Rosenbrock method
    &               , solveLA              &     ! method solving the linear systems
    &               , RosenbrockMethod     )     ! Rosenbrock methode
  END IF
  !
  !---------------------------------------------------------------
  ! --- stop timer and print output statistics
  Timer_Finish=MPI_WTIME()
  TimeIntegrationE=(TimeIntegrationE-TimeIntegrationA)
  Timer_Finish=(Timer_Finish-Timer_Start)+Time_Read
  !
  ! Print statistics
  CALL Output_Statistics( Time_Read, TimeSymbolic, TimeFac, TimeSolve, TimeRates,         &
  &                       TimeJac, TimeIntegrationE,Timer_Finish,TimeRateSend, TimeNetCDF )
  WRITE(*,*) ''
  !---------------------------------------------------------------
  ! --- Close MPI 
  CALL FinishMPI()
END PROGRAM Main_ChemKin
