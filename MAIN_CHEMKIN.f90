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
  USE mo_reac, ONLY: nspc
  USE mo_MPI
  USE mo_IO
  USE mo_ckinput
  USE Euler_Mod
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
  REAL(RealKind) :: h 
  REAL(RealKind), PARAMETER :: HR=3600.0d0
  INTEGER :: i_error
  !
  !
  !================================================================
  !===    Initialization
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
  !----------------------------------------------------------------
  ! --- Read run-file
  CALL InitRun(FileName0)
  IF (MPI_ID==0) WRITE(*,'(3X,A35)',ADVANCE='NO') 'Initialize run-file .......... done'
  !
  !----------------------------------------------------------------
  ! --- Initialize all reaction types
  CALL InitNReacType()
  !
  Tspan=(/tAnf , tEnd/)
  !
  !----------------------------------------------------------------
  !  --- read the .sys data, save coefs in sparse matrix
  OPEN(UNIT=89,FILE=ADJUSTL(TRIM(ChemFile))//'.chem',STATUS='UNKNOWN')
  IF (ChemFile(1:2)=='ck') THEN
    ! 
    IF (MPI_ID==0) WRITE(*,*) ' ---->  CHEMKIN Files : '
    ckTEMP=.TRUE.
    CALL Read_Elements(ChemFile,969)
    CALL Read_Species(ChemFile,969)
    CALL Read_Reaction(ChemFile,969)
    CALL Read_ThermoData(DataFile,696,nspc)
    !
    CALL PrintHeadSpecies(ChemFile,89)
    CALL PrintSpecies(ListGas2,89)
    CALL PrintHeadReactions(89)
    CALL PrintReactions(ReactionSystem,89,.TRUE.)       ! .TRUE. in 3rd agument for chemkin input file
    CALL PrintFinalReactions(89)
    !
    !CALL PrintReactionSystem(ReactionSystem)
  ELSE
    IF (MPI_ID==0) WRITE(*,*) ''
    Time_Read=MPI_WTIME()
    CALL ReadSystem(ChemFile)
    CALL ReadThermoData('Thermo')
    !CALL InputDatFile(DataFile)
    IF (MPI_ID==0) WRITE(*,*) '  Read sys-file ................ done'
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
    IF (MPI_ID==0) WRITE(*,*) '  Print chem-file .............. done'
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
  END IF
  CLOSE(89)
  !
  call MPI_Barrier( MPI_COMM_WORLD, i_error)
  !----------------------------------------------------------------
  ! --- Split the rate array in mpi_np parts 
  CALL BuildPartitions(MyParties,neq,MPI_np)
  !
  !----------------------------------------------------------------
  !
  !
  LWCBounds(1)=0.0d0*HR
  LWCBounds(2)=LWCBounds(1) + 1.d0*HR
  LWCBounds(3)=LWCBounds(2) + 0.25d0*HR
  LWCBounds(4)=LWCBounds(3) + 9.5d0*HR
  LWCBounds(5)=LWCBounds(4) + 0.25d0*HR
  LWCBounds(6)=LWCBounds(5) + 1.d0*HR
  !LWCBounds(1)=0.0d0
  !DO i=2,6
  !  LWCBounds(i)=Lwcbounds(i-1)+HR
  !END DO
  !
  ! LWCBounds verschieben auf x-Achse?
  !LWCBounds(:)=LWCBounds(:) - 3.0d0*HR
  !
  !----------------------------------------------------------------
  !  read .ini data
  CALL InitChem(1,Tspan(1))
  IF (MPI_ID==0) WRITE(*,*) '  Read ini-file ................ done'
  !CALL InputChemicalData(DataFile,MetFile)    ! alt
  !
  !
  ! Stop Input timer
  Time_Read=MPI_WTIME()-Time_Read
  !
  !----------------------------------------------------------------
  !   print some species (set output spc in InitFile under DIAG)
  !----------------------------------------------------------------
  !
  CALL Read_Diag(OutNetcdfSpc,OutNetcdfPhase,InitFile)
  !
  !
  !-----------------------------------------------------------------------
  ! --- Start timer and set absolut tolerance for species
  Timer_Start=MPI_WTIME()
  Atol=(/AtolGas,AtolAqua/)
  ! 
  ! print name list
  !do i=1,size(OutNetcdfSpc)
  !  print*, i, OutNetcdfSpc(i) , y_name(OutNetcdfSpc(i))
  !END do 
  !stop
  !
  ! für testzwecke!

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
  !
  !
  !
  CALL Print_Run_Param()
  !-----------------------------------------------------------------------
  ! --- choose between implicit Euler or some Rosenbrock-Wanner methode
  IF (ImpEuler==1) THEN
    !---  Read run control parameters (constant stepsize)
    !  
    CALL getarg(2,constH)             
    IF (constH=='')   THEN
      WRITE(*,*) 'Input constant stepsize:'
      READ(*,*)   constH
    END IF
    !
    READ(ConstH,*) h
    CALL IntegrateEuler(y_iconc(1:nspc),Tspan,h)
  ELSE
    CALL Integrate (  y_iconc(1:nspc)   &     ! initial concentrations activ species
    &               , Tspan             &     ! integration invervall
    &               , Atol              &     ! abs. tolerance gas species
    &               , RtolROW           &     ! rel. tolerance Rosenbrock method
    &               , solveLA           &     ! method solving the linear systems
    &               , RosenbrockMethod  )     ! Rosenbrock methode
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
!
!

! LWC = 10^-4 [l/m3]
