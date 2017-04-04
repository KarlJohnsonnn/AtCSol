   SUBROUTINE InitRun(RunFile)
!==================================================
!===  Reading and Setting of Run Control Parameters
!==================================================
!
      !USE mo_bdf
      !USE Kind_Mod,    ONLY: RealKind
      !USE mo_micfys
      USE mo_control,  ONLY:      &
!
!-----------------------------------------------------------------
!---  Scenario
!-----------------------------------------------------------------
!
!--- CHARACTER(RealKind) : Identifier for scenario
                 Bsp               &
!
!--- CHARACTER(80) : Files
&               ,MetFile           & ! Meteorology file
&               ,ChemFile          & ! Chemical mechanism
&               ,InitFile          & ! Initial concentrations
&               ,DataFile          & ! Gas and Aqeous DATA
&               ,NetcdfFileName    & ! Output NetCdfFile
!
!--- INTEGER : Unit Numbers
&               ,MetUnit           & ! Meteorology file
&               ,ChemUnit          & ! Chemical mechanism
&               ,InitUnit          & ! Initial concentrations
&               ,DataUnit          & ! Initial concentrations
!
!-- REAL(RealKind) : Set Levels and Parameters for Processes
&               ,LwcLevelmin       & ! lower Level for LWC      {l/m3]
&               ,LwcLevelmax       & ! upper Level for LWC      {l/m3]
&               ,constLWC          & ! wähle konstanten lwc wert (=0) oder lineare funktion (/=0)

!-- LOGICAL :: Print system matrices
&               ,MatrixPrint       &
&               ,DebugPrint        &
&               ,NetCdfPrint       &


!--- INTEGER : Control Parameter
&               ,pHSet             & ! Initial pH by charge balance (1=on, 0=off)
&               ,Ladebalken        &
&               ,Error_Est         &
&               ,ErrorLog          &

!--- REAL : combustion stuff
&               ,Temperature0      &
&               ,Pressure0         &

!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!
!--- REAL(8): Times ( < 0 in seconds, > 0 in hours, = 0  meteorology)
&               ,tAnf              & ! Model start time
&               ,tEnd              & ! model end time
&               ,StpNetcdf         & ! Time step for Netcdf output 
!
!---  Photolysis
&               ,idate             & ! Date: yymmdd
&               ,rlat              & ! latitude  [rad]
&               ,rlon              & ! longitude [rad]
&               ,dust              & ! dust factor (damping of photolysis)
&               ,minStp            & ! minimal time step size
&               ,maxStp            & ! maximal time step size
&               ,nOutP             & !  Number of intermediate output points per simulation
!
!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!
!--- INTEGER : Control Parameter
&               ,ImpEuler          & ! implicit euler 

!--- REAL(RealKind) : Tolerances for ROW Scheme
&               ,RtolROW           & ! Relative tolerance For ROW
&               ,AtolGas           & ! Absolute tolerance for gas phase
&               ,AtolAqua          & ! Absolute tolerance for liquid phase
&               ,AtolTemp          & ! Absolute tolerance for temperature
&               ,PI_StepSize       & ! logical for pi stepsize control

!--- CHARACTER(2) : Control Parameter
&               ,solveLA           & ! how to solve linear algebra 'cl' or 'ex'
&               ,ODEsolver         & ! choose ode solver (rosenbrock,lsode)  (integer 1,...,15)
!
!-----------------------------------------------------------------
!---  Linear Algebra
!-----------------------------------------------------------------
!
!--- Control Parameter
&               ,OrderingStrategie      & ! MUMPS ordering strategie
&               ,ParOrdering            & ! MUMPS ordering type for parallel symbolic phase
&               ,CLASSIC,EXTENDED 
   
!-----------------------------------------------------------------
      IMPLICIT NONE
      
      CHARACTER(80) :: RunFile

      
!
!-----------------------------------------------------------------
!--- NAMELISTS
      NAMELIST /SCENARIO/  Bsp,                                            &
&               LwcLevelmin, LwcLevelmax, Ladebalken,  pHSet,              &
&               constLWC, ErrorLog, MatrixPrint, NetCdfPrint, Temperature0,&
&               Pressure0, DebugPrint
!
      NAMELIST /FILES/  MetFile, ChemFile, InitFile, DataFile, NetcdfFileName, &
&               MetUnit, ChemUnit, InitUnit, DataUnit
!
      NAMELIST /TIMES/  tAnf, tEnd, idate, rlat, rlon, Dust, StpNetcdf,  &
&               minStp, maxStp, nOutP
!
      NAMELIST /NUMERICS/  RtolROW, AtolGas, AtolAqua, AtolTemp, PI_StepSize,      &
&               solveLA,  ODEsolver, ImpEuler, Error_Est
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
      constLWC = 0             ! 1 = konstanter LWC wert für ges. Simulation

!--- INTEGER : Control Parameter
      pHSet     = 1          ! Initial pH by charge balance (1=on, 0=off)
      Error_Est = 2          ! default for error estimation is euklid norm
      ErrorLog  = 0 
      MatrixPrint = .FALSE.       ! 0 = Print no matrix, 0 /= Print all matrices -> no simulation
      DebugPrint  = .FALSE.
      NetCdfPrint = .FALSE.
      Temperature0= 750.0d0
      Pressure0   = 2.0d5

!--- Read SCENARIO namelist
      READ(15,SCENARIO)

!-----------------------------------------------------------------
!---  Files
!-----------------------------------------------------------------
!
!--- CHARACTER(80) : Files
      MetFile    = 'MET/initial'             ! Meteorology file
      ChemFile   = 'CHEM/'//TRIM(RunFile)//'.sys'    ! Chemical mechanism
      DataFile   = 'CHEM/'//TRIM(RunFile)//'.dat'    ! Gas and aqueous phase data
      InitFile   = 'INI/'//TRIM(RunFile)//'.ini'     ! Initial concentrations

      NetcdfFileName = TRIM(Bsp)//'.nc'            ! Netcdf output file

!
!--- INTEGER : Unit Numbers
      MetUnit  = 15           ! Meteorology file
      ChemUnit = 10           ! Chemical mechanism
      InitUnit = 12           ! Initial concentrations
      DataUnit = 13           ! Gas and aqueous phase data
      
!--- Read FILES namelist
      READ(15,FILES)
      
!--- Adjust Filenames
      Bsp  = ADJUSTL(Bsp)

      MetFile  = ADJUSTL(MetFile)
      ChemFile = ADJUSTL(ChemFile)
      DataFile = ADJUSTL(DataFile)
      InitFile = ADJUSTL(InitFile)

      NetcdfFileName=ADJUSTL(NetcdfFileName)

!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!
!--- REAL(8): Times ( < 0 in seconds, > 0 in hours, = 0  meteorology)
      tAnf    =  0.d0        ! Model start time    [in h}
      tEnd    =  0.d0        ! model end time      [in h}

!--- REAL(8): Times in seconds.
      StpNetcdf   = 0.d0      ! Time step for Netcdf output      [in sec]

      nOutP = 100

!---  Photolysis (Here: FEBUKO chemistry-case I)
      idate = 011027          ! Date: yymmdd  (21.June 2001)
      rlat  = 5.065d+01       ! latitude  [grad] (Schmuecke)
      rlon  = 1.077d+01       ! longitude [grad] (Schmuecke)
      Dust  = 1.0d0           ! dust factor (damping of photolysis)

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

      !OrderingStrategie=7
      !WRITE(*,*) '    ORDERINGSTRATEGIE FIXED TO MUMPS = 7'

      CLOSE(15)

   END SUBROUTINE InitRun
