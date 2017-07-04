   SUBROUTINE InitRun(RunFile)
!==================================================
!===  Reading and Setting of Run Control Parameters
!==================================================
!
      USE mo_control
!
!-----------------------------------------------------------------
      IMPLICIT NONE
      
      CHARACTER(80) :: RunFile

      
!
!-----------------------------------------------------------------
!--- NAMELISTS
      NAMELIST /SCENARIO/  Bsp,                                            &
&               LwcLevelmin, LwcLevelmax, Ladebalken,  pHSet,              &
&               constLWC, ErrorLog, MatrixPrint, NetCdfPrint, Temperature0,&
&               Pressure0, DebugPrint, ChemKin
!
      NAMELIST /FILES/  MetFile, SysFile,MWeights, InitFile, DataFile,    &
&               NetcdfFileName, MetUnit, ChemUnit, MWUnit, InitUnit,       &
&               DataUnit, Targets
!
      NAMELIST /TIMES/  tAnf, tEnd, idate, rlat, rlon, Dust, StpNetcdf,  &
&               minStp, maxStp, nOutP
!
      NAMELIST /NUMERICS/  RtolROW, AtolGas, AtolAqua, AtolTemp, PI_StepSize,      &
&               solveLA,  ODEsolver, ImpEuler, Error_Est, Vectorized
!
      NAMELIST /ORDERING/  OrderingStrategie, ParOrdering
!      

!
!===================================================================
!===  Set and Read Default Values
!==================================================
!
!--- ROpen run control file
      OPEN(UNIT=15,FILE=RunFile)

!-----------------------------------------------------------------
!---  Scenario
!-----------------------------------------------------------------
!
!--- CHARACTER(20) : Suffix
      Bsp = ''                                     ! Identifier of scenario
!
!-- REAL(8) : Set Levels and Parameters for Processes
      LwcLevelmin = 1.0d-12       ! Lower level for LWC      {l/m3]
      LwcLevelmin = 3.0d-10       ! Lower level for LWC      {l/m3]
      constLWC = .FALSE.          ! true = konstanter LWC wert für ges. Simulation

!--- INTEGER : Control Parameter
      pHSet     = 1          ! Initial pH by charge balance (1=on, 0=off)
      Error_Est = 2          ! default for error estimation is euklid norm
      ErrorLog  = 0 
      MatrixPrint = .FALSE.       ! 0 = Print no matrix, 0 /= Print all matrices -> no simulation
      DebugPrint  = .FALSE.
      NetCdfPrint = .FALSE.
      ChemKin     = .FALSE.
      Temperature0= 750.0d0
      Pressure0   = 2.0d5

!--- Read SCENARIO namelist
      READ(15,SCENARIO)
      
      IF (ChemKin) Teq = .TRUE.

!-----------------------------------------------------------------
!---  Files
!-----------------------------------------------------------------
!
!--- CHARACTER(80) : Files
      MetFile  = 'MET/initial'             ! Meteorology file
      SysFile  = 'CHEM/'//TRIM(RunFile)//'.sys'    ! Chemical mechanism
      DataFile = 'CHEM/'//TRIM(RunFile)//'.dat'    ! Gas and aqueous phase data
      InitFile = 'INI/'//TRIM(RunFile)//'.ini'     ! Initial concentrations
      Targets  = 'RED/'//TRIM(RunFile)//'.ctrl'    ! Controlfile for mechanism reduction

      NetcdfFileName = TRIM(Bsp)//'.nc'            ! Netcdf output file

!
!--- INTEGER : Unit Numbers
      MetUnit  = 15           ! Meteorology file
      ChemUnit = 10           ! Chemical mechanism
      MWUnit   = 19           ! Chemical mechanism
      InitUnit = 12           ! Initial concentrations
      DataUnit = 13           ! Gas and aqueous phase data
      
!--- Read FILES namelist
      READ(15,FILES)
      
!--- Adjust Filenames
      Bsp  = ADJUSTL(Bsp)

      MetFile  = ADJUSTL(MetFile)
      SysFile  = ADJUSTL(SysFile)
      DataFile = ADJUSTL(DataFile)
      MWeights = ADJUSTL(MWeights)
      InitFile = ADJUSTL(InitFile)
      Targets  = ADJUSTL(Targets)

      NetcdfFileName=ADJUSTL(NetcdfFileName)

!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!
!--- REAL(8): Times ( < 0 in seconds, > 0 in hours, = 0  meteorology)
      tAnf    =  0.0_dp        ! Model start time    [in h}
      tEnd    =  0.0_dp        ! model end time      [in h}

!--- REAL(8): Times in seconds.
      StpNetcdf   = 0.0_dp      ! Time step for Netcdf output      [in sec]

      nOutP = 100

!---  Photolysis (Here: FEBUKO chemistry-case I)
      idate = 011027          ! Date: yymmdd  (21.June 2001)
      rlat  = 50.65_dp       ! latitude  [grad] (Schmuecke)
      rlon  = 10.77_dp       ! longitude [grad] (Schmuecke)
      Dust  = 1.0_dp           ! dust factor (damping of photolysis)

!--- Read TIMES namelist
      READ(15,TIMES)

      ! minimum output steps are 2
      IF (nOutP<2) nOutP = 2

!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!

!--- REAL(RealKind) : Tolerances for ROW Scheme
      RtolROW  = 1.0d-5               ! Relative tolerance For ROW
      AtolGas  = 1.0d-7               ! Absolute tolerance for gas phase
      AtolAqua = 1.0d-7               ! Absolute tolerance for liquid phase
      AtolTemp = 1.0d-7               ! Absolute tolerance for temperature
      PI_StepSize = .FALSE.
      solveLA  = 'ex'                 ! method of solving linear algebra
      ODEsolver  = 'ROS34PW3'  ! ROW scheme
      ImpEuler = 0                    ! 1 for implicit euler integration
      Vectorized = .TRUE.
      
!
!--- Read NUMERICS namelist
      READ(15,NUMERICS)

      IF ( solveLA=='cl' ) CLASSIC  = .TRUE.
      IF ( solveLA=='ex' ) EXTENDED = .TRUE.
!
!-----------------------------------------------------------------
!---  Linear Algebra
!-----------------------------------------------------------------
!
!--- Control Parameter
      READ(15,ORDERING)

      IF (OrderingStrategie<8 .OR. ParOrdering>=0 ) useMUMPS = .TRUE.
      IF (OrderingStrategie>=8) useSparseLU = .TRUE.

      !OrderingStrategie=7
      !WRITE(*,*) '    ORDERINGSTRATEGIE FIXED TO MUMPS = 7'

      CLOSE(15)

   END SUBROUTINE InitRun
