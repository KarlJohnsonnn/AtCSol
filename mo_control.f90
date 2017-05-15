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
&                     ,MWeights          & ! molecular weights of species 
&                     ,NetcdfFileName    & ! NetCDF output file
&                     ,ODEsolver    ! Method for Rosenbrock Integration
!
!--- Unit Numbers
      INTEGER :: MetUnit           & ! Meteorology file
&               ,ChemUnit          & ! Chemical mechanism
&               ,MWUnit            & ! molecular weights unit
&               ,InitUnit          & ! Initial concentrations
&               ,DataUnit            ! Gas and Aqueous DATA
!
!-- Set Levels and Parameters for Processes
      REAL(dp) :: &
&                LwcLevelmin       & ! Lower level for LWC
&               ,LwcLevelmax         ! Upper level for LWC


!-- print matrices of the system (alpha, beta, miter, lu_miter, permutvector,..)
      LOGICAL :: MatrixPrint   
      LOGICAL :: DebugPrint   
      LOGICAL :: NetCdfPrint   
      LOGICAL :: constLWC   

!--- Control Parameter
      INTEGER :: pHSet             & ! Initial pH by charge balance (1=on, 0=off)
&               ,Ladebalken        & ! ladebalken im terminal bei simulation (=1, default=0)
&               ,Error_Est         & ! error estimation 1 = inf norm  , 2 = euklid norm
&               ,ErrorLog            ! if = 0 do not print error log 
    
!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!--- Times:  < 0 in seconds, > 0 in hours, = 0  meteorology
REAL(dp) :: tAnf              & ! Model start time
&               , tEnd              & ! Model end time
&               , StpNetcdf           ! Time step for Netcdf output
INTEGER        :: nOutP
!
!--- Initial temperature and pressure for combustion mechanisms
      REAL(dp) :: Temperature0
      REAL(dp) :: Pressure0  

!---  Photolysis
      INTEGER :: iDate               ! Current date
      REAL(dp) :: rlat, rlon          ! Latitude, Longitude


!--  dust factor for damping of photolysis rates, measured JNO2
      REAL(dp) :: Dust = 1.0d0
      REAL(dp) :: minStp = 1.0d-25
      REAL(dp) :: maxStp = 100.0d0

!-- initialize timers
      REAL(dp) :: Timer_Start=0.0d0, Timer_Finish=0.d0
      REAL(dp) :: Tspan(2)
  

      REAL(dp) :: TimeFac=0.0d0
      REAL(dp) :: TimeSolve=0.0d0
      REAL(dp) :: TimeRates=0.0d0
      REAL(dp) :: TimeRateSend=0.0d0
      REAL(dp) :: TimeJac=0.0d0
      REAL(dp) :: Time_Read=0.0d0
      REAL(dp) :: TimeSymbolic=0.0d0
      REAL(dp) :: TimeNetCDF=0.0d0
      REAL(dp) :: TimeErrCalc=0.0d0
      REAL(dp) :: TimeRhsCalc=0.0d0
  
      REAL(dp) :: TimeIntegrationA=0.0d0
      REAL(dp) :: TimeIntegrationE=0.0d0
      REAL(dp) :: TimeRateA=0.0d0
      REAL(dp) :: TimeRateE=0.0d0
      REAL(dp) :: TimeRateSendA=0.0d0
      REAL(dp) :: TimeRateSendE=0.0d0
      REAL(dp) :: TimeJacobianA=0.0d0
      REAL(dp) :: TimeJacobianE=0.0d0
      REAL(dp) :: TimeNetCDFA=0.0d0
!
!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!
!--- Control Parameter
      INTEGER :: ImpEuler
!
!--- Tolerances for ROW Scheme
      REAL(dp) :: RtolROW        & ! Relative tolerance for ROW method
&                     , AtolGas        & ! Absolute tolerance for gas phase
&                     , AtolAqua       & ! Absolute tolerance for liquid phase
&                     , AtolTemp         ! Absolute tolerance for Temperatur
      REAL(dp), ALLOCATABLE :: ThresholdStepSizeControl(:)
      REAL(dp), ALLOCATABLE :: ATolAll(:)

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
    REAL(dp), PARAMETER :: HOUR       = 3600.d0            &
&                              , secday     = 4.32d04            &
&                              , hourday    = 24.0d0             &
&                              , mONE       = -1.d0              &
&                              , Pi         = 4.0d0*ATAN(1.0d0)  &
&                              , DR         = Pi / 180.d0        &  
&                              , PiHalf     = 2.0d0*ATAN(1.0d0)  & 
&                              , Pi34       = 3.0d0/4.0d0/Pi     & 
&                              , eps        = EPSILON(1.0d0)     &  ! such that 1+eps>1 with working precision
&                              , small      = TINY(1.0d0)        &  ! smallest pos. real value
&                              , epsY       = 1.0d-7             &
&                              , SI_am      = 1.66053892173d-27  &  ! Atomic mass unit  [kg]  
&                              , SI_na      = 6.0221412927d+23   &  ! Avogadro's number [1/mol]
&                              , SI_kB      = 1.380648813d-23    &  ! Bolzmann constant [J/K]
&                              , SI_Gas     = SI_na * SI_kB         ! Gas constant      [J/mol/K]
!
!--- Real number constants
  REAL(dp), PARAMETER ::   ZERO    =     0.0d0   & 
&                       ,  ONE     =     1.0d0   & 
&                       ,  TWO     =     2.0d0  ,   rTWO     =   ONE/tWO     &
&                       ,  THREE   =     3.0d0  ,   rTHREE   =   ONE/THREE   &
&                       ,  FOUR    =     4.0d0  ,   rFOUR    =   ONE/FOUR    &
&                       ,  FIVE    =     5.0d0  ,   rFIVE    =   ONE/FIVE    &
&                       ,  SIX     =     6.0d0  ,   rSIX     =   ONE/SIX     &
&                       ,  SEVEN   =     7.0d0  ,   rSEVEN   =   ONE/SEVEN   &
&                       ,  EIGHT   =     8.0d0  ,   rEIGHT   =   ONE/EIGHT   &
&                       ,  NINE    =     9.0d0  ,   rNINE    =   ONE/NINE    &
&                       ,  TEN     =    10.0d0  ,   rTEN     =   ONE/TEN     &
&                       ,  ELEVN   =    11.0d0  ,   rELEVN   =   ONE/ELEVN   &
&                       ,  TWELV   =    12.0d0  ,   rTWELV   =   ONE/TWELV   &
&                       , TWENTY   =    20.0d0  ,   rTWENTY  =   ONE/TWENTY  &
&                       , mTHIRTY  =   -30.0d0  &
&                       , rm300    = mONE/300.d0,   r300     =   ONE/300.d0
!
!--- Orders of magnitude
  REAL(dp), PARAMETER ::   nano    =     1.0d-09    &
&                             , micro    =     1.0d-06    &
&                             , milli    =     1.0d-03    &
&                             , kilo     =     1.0d+03    &
&                             , mega     =     1.0d+06    &
&                             , tera     =     1.0d+09
!
!--- Natural logarithms
  REAL(dp), PARAMETER ::    ln10   =     LOG(TEN)    &
&                             ,  rln10   = ONE/LOG(TEN)
!
!--- minimum values if there is no sun
  REAL(dp), PARAMETER ::    EyChiZmin  =  9.357d-14
!
!
!--- Unit Conversion constants
!
!     Pressure
  REAL(dp)      , PARAMETER :: bar_to_dyncm2 = 1.0d06
  REAL(dp)      , PARAMETER :: dyncm2_to_Pa  = 1.0d-01
  REAL(dp)      , PARAMETER :: Pa_to_dyncm2  = 1.0d+01
  REAL(dp)      , PARAMETER :: bar_to_Pa     = 1.0d05
  REAL(dp)      , PARAMETER :: atm_to_Pa     = 101325.0d0
!
!     Energy
  REAL(dp)      , PARAMETER :: cal_to_joule  = 4.184d0
  REAL(dp)      , PARAMETER :: joule_to_cal  = ONE / cal_to_joule
  REAL(dp)      , PARAMETER :: joule_to_kcal = milli * joule_to_cal
  REAL(dp)      , PARAMETER :: kcal_to_joule = kilo * cal_to_joule
  REAL(dp)      , PARAMETER :: joule_to_erg  = TEN * mega
  REAL(dp)      , PARAMETER :: erg_to_joule  = rTEN * micro
!
!--- Physical constants ********************************************

!     Universal gas constant [J / mol / K]
  REAL(dp), PARAMETER :: R         = 8.31446210000000D0     
  REAL(dp), PARAMETER :: rR        = ONE/R

!
!     Universal gas constant, calorie units [cal / mol / K]
  REAL(dp), PARAMETER :: Rcal      = R * joule_to_cal
  REAL(dp), PARAMETER :: rRcal     = ONE/Rcal
!
!     Universal gas constant, CGS units [erg / mol / K]
  REAL(dp), PARAMETER :: Rerg      = R * joule_to_erg
  REAL(dp), PARAMETER :: rRerg      = ONE/Rerg
!
!     Standard pressure [Pa]
  REAL(dp), PARAMETER :: Patm      = atm_to_Pa
  REAL(dp), PARAMETER :: rPatm     = ONE/Patm


!--- Unit Conversion
    INTEGER :: GasUnit, AquaUnit, GasRateUnit
    REAL(dp) :: GasFac

    REAL(dp) :: ConvGas, ConvAir


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
    REAL(dp), ALLOCATABLE ::  OutImp(:,:)        ! boundary radii of impactor stages
!
!-- logicals for classic or extended matrix case
    LOGICAL :: CLASSIC  = .FALSE.
    LOGICAL :: EXTENDED = .FALSE.

    LOGICAL :: useMUMPS = .FALSE.
    LOGICAL :: useSparseLU = .FALSE.

    LOGICAL :: TempEq  = .FALSE.
    LOGICAL :: ChemKin = .FALSE.

    LOGICAL :: Vectorized = .FALSE.
    
 END MODULE mo_control

