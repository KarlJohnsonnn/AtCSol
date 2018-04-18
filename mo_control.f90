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
      CHARACTER(80) :: Bsp        = ''      ! Identifier for scenario

!--- Files
      CHARACTER(80) :: RunFile    = ''    & ! Simulation data file
&                    , MetFile    = ''    & ! Meteorology file
&                    , SysFile    = ''    & ! Chemical mechanism
&                    , ChemFile   = ''    & ! Chemical mechenism (spaccim input)
&                    , InitFile   = ''    & ! Initial concentrations
&                    , DataFile   = ''    & ! Gas and Aqueous DATA
&                    , MWFile     = ''    & ! molecular weights of species 
&                    , NetcdfFile = ''    & ! NetCDF output file
&                    , ODEsolver  = ''    & ! Method for Rosenbrock Integration
&                    , TargetFile = ''      ! file for reductions analysis (target species)

      CHARACTER(7)  :: OutputPath   = 'OUTPUT/'      ! path to output folder
      CHARACTER(19) :: FluxMetaFile = 'OUTPUT/fluxmeta.dat' ! meta data for unformatted flux data
      CHARACTER(17) :: FluxFile     = 'OUTPUT/fluxes.dat'   ! flux data (unformatted)
!
!--- Unit Numbers
      INTEGER, PARAMETER :: RunUnit      = 101  & 
&                         , MetUnit      = 102  & 
&                         , SysUnit      = 103  & 
&                         , ChemUnit     = 104  & 
&                         , MWUnit       = 105  & 
&                         , InitUnit     = 106  & 
&                         , DataUnit     = 107  & 
&                         , FluxMetaUnit = 109  &
&                         , FluxUnit     = 110  & 
&                         , TikZUnit     = 111   
!
!-- Set Levels and Parameters for Processes
      REAL(dp) :: LWCLevelmin    & ! Lower level for LWC
&               , LWCLevelmax    & ! Upper level for LWC
&               , Temperature0   & ! Initial temperature
&               , Pressure0        ! Initial pressure for combustion mechanisms


!-- print matrices of the system (alpha, beta, miter, lu_miter, permutvector,..)
      LOGICAL :: MatrixPrint     & ! print certain matrices to file
&              , DebugPrint      & ! debugging rosenbrock steps
&              , NetCdfPrint     & ! print out concs and other stuff to netcdf file
&              , constLWC        & ! true if lwc value is fixed
&              , Lehmann         & ! prints out pathway analysis file 
&              , Teq             & ! if Teq=.TRUE. simulate combustion mechanism
&              , pHSet           & ! Initial pH by charge balance (1=on, 0=off)
&              , WaitBar         & ! ladebalken im terminal bei simulation (=1, default=0)
&              , FluxAna         & ! writing flux data and analyse after simulaiton -> print new reaction file
&              , Simulation        ! calculation of species concentration 

      INTEGER :: Error_Est         ! error estimation 1 = inf norm  , 2 = euklid norm
    
!--- Output of reaction fluxes
      REAL(dp)                    :: StpFlux
      INTEGER                     :: iStpFlux

      INTEGER,          PARAMETER :: newReac_nr    = 1977
      CHARACTER(LEN=*), PARAMETER :: newReac_name  = 'BigRates.dat'
      INTEGER, ALLOCATABLE        :: newReac_List(:)
!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!--- Times:  < 0 in seconds, > 0 in hours, = 0  meteorology
      REAL(dp) :: tBegin              & ! Model start time
&               , tEnd              & ! Model end time
&               , StpNetcdf           ! Time step for Netcdf output
      INTEGER  :: nOutP

!--- NetCDF globals      
      INTEGER, ALLOCATABLE :: iNcdfGas(:), iNcdfAqua(:), iNcdfSolid(:), iNcdfParti(:)
      INTEGER, ALLOCATABLE :: iNCout_G(:), iNCout_A_l(:), iNCout_A_m3(:), iNCout_S(:), iNCout_P(:)
      INTEGER  :: nNcdfGas=0, nNcdfAqua=0, nNcdfSolid=0, nNcdfParti=0 ! number of output spc for each phase
!

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
  

      REAL(dp) :: Time_Read=0.0d0
      REAL(dp) :: TimeRates=0.0d0
      REAL(dp) :: TimeSymbolic=0.0d0
      REAL(dp) :: TimeFac=0.0d0
      REAL(dp) :: TimeSolve=0.0d0
      REAL(dp) :: TimeJac=0.0d0
      REAL(dp) :: TimeNetCDF=0.0d0
      REAL(dp) :: TimeErrCalc=0.0d0
      REAL(dp) :: TimeFluxWrite=0.0d0
      REAL(dp) :: TimeRhsCalc=0.0d0
      REAL(dp) :: TimeReduction=0.0d0
  
      REAL(dp) :: TimeIntegration=0.0d0
      REAL(dp) :: TimeRateA=0.0d0
      REAL(dp) :: TimeRateE=0.0d0
      REAL(dp) :: TimeJacobianA=0.0d0
      REAL(dp) :: TimeJacobianE=0.0d0
      REAL(dp) :: TimeNetCDFA=0.0d0
!

!--- type for some statistics
      TYPE Output_T
        REAL(dp), ALLOCATABLE :: y(:)    ! y-vector at Tend
        INTEGER :: nsteps     = 0              ! # succ. steps
        INTEGER :: nfailed    = 0              ! # failed steps
        INTEGER :: nRateEvals = 0              ! # Rate evaluation
        INTEGER :: npds       = 0              ! # Jacobian evaluation
        INTEGER :: ndecomps   = 0              ! # LU factorisation
        INTEGER :: nsolves    = 0              ! # solved lin algebra
        REAL(dp) :: Ttimestep = 0.0d0  ! mean Time for one ROW step
      END TYPE Output_T
    
      TYPE(Output_T) :: Out
!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!
!
!--- Tolerances for ROW Scheme
      REAL(dp) :: RtolROW        & ! Relative tolerance for ROW method
&               , AtolGas        & ! Absolute tolerance for gas phase
&               , AtolAqua       & ! Absolute tolerance for liquid phase
&               , AtolTemp         ! Absolute tolerance for Temperatur

      INTEGER  :: Ordering     ! default = 8
      INTEGER  :: ParOrdering  ! < 0 serial ordering , 0,1,2 parallel

      
      REAL(dp), ALLOCATABLE :: ATolAll(:)

!--- Logical variable for PI setsize controler
      LOGICAL :: PI_StepSize
      
!--- Control Parameter for solving linear algebra  
      ! 'cl' for classic jacobian calc or
      ! 'ex' for extended matrix without calc of jac
      CHARACTER(2) :: LinAlg = '??'

!
!-----------------------------------------------------------------
!---  Set Constants and Unit Conversion
!-----------------------------------------------------------------
!
!---  constants
    REAL(dp), PARAMETER :: HOUR       = 3600.0             &
&                        , hourday    = 86400.0            &
&                        , secday     = 4.32d04            &
&                        , mONE       = -1.d0              &
&                        , Pi         = 4.0d0*ATAN(1.0d0)  &
&                        , DR         = Pi / 180.d0        &  
&                        , PiHalf     = 2.0d0*ATAN(1.0d0)  & 
&                        , Pi34       = 3.0d0/4.0d0/Pi     & 
&                        , eps        = EPSILON(1.0d0)     &  ! such that 1+eps>1 with working precision
&                        , small      = TINY(1.0d0)        &  ! smallest pos. real value
&                        , big        = HUGE(1.0d0)        &  ! largest pos. real value
&                        , epsY       = 1.0d-7             &
&                        , SI_am      = 1.66053892173d-27  &  ! Atomic mass unit  [kg]  
&                        , SI_na      = 6.0221412927d+23   &  ! Avogadro's number [1/mol]
&                        , SI_kB      = 1.380648813d-23    &  ! Bolzmann constant [J/K]
&                        , SI_Gas     = SI_na * SI_kB         ! Gas constant      [J/mol/K]
!
!--- Real number constants
  REAL(dp), PARAMETER ::   ZERO    =     0.0d0   & 
&                      ,  ONE     =     1.0d0   & 
&                      ,  TWO     =     2.0d0  ,   rTWO     =   ONE/tWO     &
&                      ,  THREE   =     3.0d0  ,   rTHREE   =   ONE/THREE   &
&                      ,  FOUR    =     4.0d0  ,   rFOUR    =   ONE/FOUR    &
&                      ,  FIVE    =     5.0d0  ,   rFIVE    =   ONE/FIVE    &
&                      ,  SIX     =     6.0d0  ,   rSIX     =   ONE/SIX     &
&                      ,  SEVEN   =     7.0d0  ,   rSEVEN   =   ONE/SEVEN   &
&                      ,  EIGHT   =     8.0d0  ,   rEIGHT   =   ONE/EIGHT   &
&                      ,  NINE    =     9.0d0  ,   rNINE    =   ONE/NINE    &
&                      ,  TEN     =    10.0d0  ,   rTEN     =   ONE/TEN     &
&                      ,  ELEVN   =    11.0d0  ,   rELEVN   =   ONE/ELEVN   &
&                      ,  TWELV   =    12.0d0  ,   rTWELV   =   ONE/TWELV   &
&                      , TWENTY   =    20.0d0  ,   rTWENTY  =   ONE/TWENTY  &
&                      , mTHIRTY  =   -30.0d0  &
&                      , rm300    = mONE/300.d0,   r300     =   ONE/300.d0
!
!--- Orders of magnitude
  REAL(dp), PARAMETER ::  nano    =     1.0d-09    &
&                      , micro    =     1.0d-06    &
&                      , milli    =     1.0d-03    &
&                      , kilo     =     1.0d+03    &
&                      , mega     =     1.0d+06    &
&                      , giga     =     1.0d+09
!
!--- Natural logarithms
  REAL(dp), PARAMETER ::   ln10   =     LOG(TEN)    &
&                      ,  rln10   = ONE/LOG(TEN)
!
!--- minimum values if there is no sun
  REAL(dp), PARAMETER ::    EyChiZmin  =  9.357d-14
!
!
!--- Unit Conversion constants
!
!     Pressure
  REAL(dp), PARAMETER :: bar_to_dyncm2 = 1.0d06
  REAL(dp), PARAMETER :: dyncm2_to_Pa  = 1.0d-01
  REAL(dp), PARAMETER :: Pa_to_dyncm2  = 1.0d+01
  REAL(dp), PARAMETER :: bar_to_Pa     = 1.0d05
  REAL(dp), PARAMETER :: atm_to_Pa     = 101325.0d0
!
!     Energy
  REAL(dp), PARAMETER :: cal_to_joule  = 4.184d0
  REAL(dp), PARAMETER :: joule_to_cal  = ONE / cal_to_joule
  REAL(dp), PARAMETER :: joule_to_kcal = milli * joule_to_cal
  REAL(dp), PARAMETER :: kcal_to_joule = kilo * cal_to_joule
  REAL(dp), PARAMETER :: joule_to_erg  = TEN * mega
  REAL(dp), PARAMETER :: erg_to_joule  = rTEN * micro
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


!-----------------------------------------------------------------
!---  Output control 
!-----------------------------------------------------------------

!
!-- input chemical reaction mechnism
    INTEGER, PARAMETER :: LenLine=400
    INTEGER, PARAMETER :: LenName=100
    INTEGER, PARAMETER :: LenType=20

!-- logicals for classic or extended matrix case
    LOGICAL :: CLASSIC  = .FALSE.
    LOGICAL :: EXTENDED = .FALSE.

    LOGICAL :: useMUMPS = .FALSE.
    LOGICAL :: useSparseLU = .FALSE.

    LOGICAL :: ChemKin = .FALSE.

    REAL(dp), ALLOCATABLE :: integrated_rates(:)
    REAL(dp), ALLOCATABLE :: mixing_ratios_spc(:,:)

!-- Type declaration 
    TYPE List
      INTEGER,   ALLOCATABLE  :: List(:)
      Real(dp),  ALLOCATABLE  :: ListE(:)
      Real(dp),  ALLOCATABLE  :: ListP(:)
      INTEGER                :: len
    END TYPE List

    TYPE Chain
      CHARACTER(100), ALLOCATABLE :: sName(:,:)
      INTEGER,        ALLOCATABLE :: sIdx(:,:)
      Type(List),     ALLOCATABLE :: rIdx(:)
    END TYPE Chain

    TYPE LUMP_SPC
      CHARACTER(LenName)              :: SuperSpc
      CHARACTER(LenName), ALLOCATABLE :: cSingleSpc(:)
      INTEGER,            ALLOCATABLE :: iSingleSpc(:)
    END TYPE LUMP_SPC

    TYPE(LUMP_SPC), ALLOCATABLE :: LUMP(:)

    TYPE Families_T
      CHARACTER(LenName), ALLOCATABLE :: Name(:)
      INTEGER,            ALLOCATABLE :: Index(:)
    END TYPE Families_T
    


    REAL(dp) :: eps_red       ! threshold for reduction procedure

    
    INTEGER, ALLOCATABLE :: maxErrorCounter(:)
   
 END MODULE mo_control

