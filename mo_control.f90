!-------------------------------------------------------------
!--  Modul for the Saving Control Values
!-------------------------------------------------------------
 MODULE mo_control
   USE Kind_Mod
!
!-----------------------------------------------------------------
!---  Scenario
!-----------------------------------------------------------------
!
!--- Identifier for scenario
      CHARACTER(20) :: Bsp           ! Identifier for scenario

!--- Files
      CHARACTER(80) :: MetFile           & ! Meteorology file
&                     ,ChemFile          & ! Chemical mechanism
&                     ,InitFile          & ! Initial concentrations
&                     ,DataFile          & ! Gas and Aqueous DATA
&                     ,NetcdfFileName    & ! NetCDF output file
&                     ,RosenbrockMethod    ! Method for Rosenbrock Integration
!
!--- Unit Numbers
      INTEGER :: MetUnit           & ! Meteorology file
&               ,ChemUnit          & ! Chemical mechanism
&               ,InitUnit          & ! Initial concentrations
&               ,DataUnit            ! Gas and Aqueous DATA
!
!-- Set Levels and Parameters for Processes
      REAL(RealKind) :: &
&                LwcLevelmin       & ! Lower level for LWC
&               ,LwcLevelmax         ! Upper level for LWC


!-- print matrices of the system (alpha, beta, miter, lu_miter, permutvector,..)
      LOGICAL :: MatrixPrint   
      LOGICAL :: NetCdfPrint   

!--- Control Parameter
      INTEGER :: pHSet             & ! Initial pH by charge balance (1=on, 0=off)
                ,constLWC          & ! with cloud constLWC>=1
&               ,Ladebalken        & ! ladebalken im terminal bei simulation (=1, default=0)
&               ,Error_Est         & ! error estimation 1 = inf norm  , 2 = euklid norm
&               ,ErrorLog            ! if = 0 do not print error log 
!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!--- Times:  < 0 in seconds, > 0 in hours, = 0  meteorology
REAL(RealKind) :: tAnf              & ! Model start time
&               , tEnd              & ! Model end time
&               , StpNetcdf           ! Time step for Netcdf output
!

!---  Photolysis
      INTEGER :: iDate               ! Current date
      REAL(RealKind) :: rlat, rlon          ! Latitude, Longitude


!--  dust factor for damping of photolysis rates, measured JNO2
      REAL(RealKind) :: Dust = 1.0d0
      REAL(RealKind) :: minStp = 1.0d-25
      REAL(RealKind) :: maxStp = 100.0d0

!-- initialize timers
      REAL(RealKind) :: Timer_Start=0.0d0, Timer_Finish=0.d0
      REAL(RealKind) :: Tspan(2)
  

      REAL(RealKind) :: TimeFac=0.0d0
      REAL(RealKind) :: TimeSolve=0.0d0
      REAL(RealKind) :: TimeRates=0.0d0
      REAL(RealKind) :: TimeRateSend=0.0d0
      REAL(RealKind) :: TimeJac=0.0d0
      REAL(RealKind) :: Time_Read=0.0d0
      REAL(RealKind) :: TimeSymbolic=0.0d0
      REAL(RealKind) :: TimeNetCDF=0.0d0
  
      REAL(RealKind) :: TimeIntegrationA=0.0d0
      REAL(RealKind) :: TimeIntegrationE=0.0d0
      REAL(RealKind) :: TimeRateA=0.0d0
      REAL(RealKind) :: TimeRateE=0.0d0
      REAL(RealKind) :: TimeRateSendA=0.0d0
      REAL(RealKind) :: TimeRateSendE=0.0d0
      REAL(RealKind) :: TimeJacobianA=0.0d0
      REAL(RealKind) :: TimeJacobianE=0.0d0
      REAL(RealKind) :: TimeNetCDFA=0.0d0
!
!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!
!--- Control Parameter
      INTEGER :: ImpEuler
!
!--- Tolerances for ROW Scheme
      REAL(RealKind) :: RtolROW        & ! Relative tolerance for ROW method
&                     , AtolGas        & ! Absolute tolerance for gas phase
&                     , AtolAqua       & ! Absolute tolerance for liquid phase
&                     , AtolTemp         ! Absolute tolerance for Temperatur
      REAL(RealKind), ALLOCATABLE :: ThresholdStepSizeControl(:)
      REAL(RealKind), ALLOCATABLE :: ATolAll(:)

!--- Logical variable for PI setsize controler
      LOGICAL :: PI_StepSize
      
!--- Control Parameter for solving linear algebra  
      CHARACTER(2) :: solveLA  ! method of solving the linear system
      ! 'cl' for classic jacobian calc or
      ! 'ex' for extended matrix without calc of jac
      
      INTEGER :: ROWcoef       ! chose ROW method
!
!-----------------------------------------------------------------
!---  Linear Algebra
!-----------------------------------------------------------------
!
!--- Control Parameter
      INTEGER :: OrderingStrategie = 7  ! default = automatic choice
      INTEGER :: ParOrdering = -1       ! < 0 serial ordering , 0,1,2 parallel
!
!-----------------------------------------------------------------
!---  Set Constants and Unit Conversion
!-----------------------------------------------------------------
!
!---  constants
    REAL(RealKind), PARAMETER :: mol2part   = 6.02295d17         &
&                              , GasConst_R = 0.082056d0         &   ! [in l*atm/mol/K]
&                              , hour       = 3600.d0            &
&                              , secday     = 4.32d04            &
&                              , ZERO       = 0.d0               &
&                              , ONE        = 1.d0               &
&                              , mONE       = -1.d0              &
&                              , Pi         = 4.0d0*ATAN(1.0d0)  &
&                              , DR         = Pi / 180.d0        &
&                              , PiHalf     = 2.0d0*ATAN(1.0d0)  &
&                              , eps        = EPSILON(1.0d0)     &
&                              , epsY       = 1.0d-7

!--- Unit Conversion
    INTEGER :: GasUnit, AquaUnit, GasRateUnit
    REAL(RealKind) :: GasFac

    REAL(RealKind) :: ConvGas, ConvAir


!-----------------------------------------------------------------
!---  Output control 
!-----------------------------------------------------------------


!-- number of Qt-Output variable
    INTEGER :: nQtGas  = 0      & ! number of gas phase species
&             ,nQtAqua = 0        ! number of aqueous phase species

    INTEGER, PARAMETER :: nphys = 9     & ! microphysical output values
&                        ,nmet  = 12       ! meteorological output values

!--  Qt-Output species
    INTEGER, ALLOCATABLE :: IndQtGas(:), IndQtAqua(:)    ! Flags for Qt-Output of each species
!
!-- control of output for different impactor stages (in measurements)
    INTEGER ::  nOutImp = 0                     ! number of output impactor stages
    INTEGER, ALLOCATABLE ::  IndOutImp(:,:)     ! start and end fraction indices
    REAL(RealKind), ALLOCATABLE ::  OutImp(:,:)        ! boundary radii of impactor stages
!

 END MODULE mo_control

