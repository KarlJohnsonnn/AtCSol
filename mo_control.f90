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
&                     ,AVSFile           & ! AVS   output file
&                     ,NetcdfFileName    & ! NetCDF output file
&                     ,OutFile           & ! ASCII output file
&                     ,DiagFile          & ! Diagnose file
&                     ,StatFile          & ! Statistics file
&                     ,JacFile           & ! Jacobian structure file
&                     ,QtListFile        & ! File with listed species for Qt-Output 
&                     ,ErrFile           & ! Output error file
&                     ,RosenbrockMethod    ! Method for Rosenbrock Integration
!
!--- Unit Numbers
      INTEGER :: MetUnit           & ! Meteorology file
&               ,ChemUnit          & ! Chemical mechanism
&               ,InitUnit          & ! Initial concentrations
&               ,DataUnit          & ! Gas and Aqueous DATA
&               ,AVSUnit           & ! AVS   output file
&               ,OutUnit           & ! ASCII output file
&               ,DiagUnit          & ! Diagnose file
&               ,StatUnit          & ! Statistics file
&               ,JacUnit           & ! Jacobian structure file
&               ,ErrUnit             ! Output error file
!
!-- Set Levels and Parameters for Processes
      REAL(RealKind) :: &
&                LwcLevelmin       & ! Lower level for LWC
&               ,LwcLevelmax       & ! Upper level for LWC
&               ,DryLevel          & ! Lower level for DryMass
&               ,MinScav           & ! Minimum radius for gas uptake
&               ,MinChem           & ! Minimum radius for chemistry
&               ,LwcAll            & ! Liqid Water Content
&               ,MeanRad             ! Mean Droplet Radius


!-- print matrices of the system (alpha, beta, miter, lu_miter, permutvector,..)
      LOGICAL :: MatrixPrint   
      LOGICAL :: NetCdfPrint   

!--- Dimensions and Resolution
      INTEGER :: ntfrac            & ! Number of fractions
&               ,ResFrac           & ! Resolution of chemistry spectrum
&               ,FirstChem         & ! First fraction for chemistry
&               ,LastChem          & ! Last  fraction for chemistry
&               ,MaxFrac           & ! Maximum fraction for chemistry
&               ,nReso             & ! Resolution for smoothing   
&               ,nSmooth           & ! Number of fine fractions for smoothing
&               ,nCell               ! Number of grid cells

!--- Control Parameter
      INTEGER :: nCouple           & ! Coupling flag
&               ,nDep              & ! Flag for deposition gas phase
&               ,nDepAqua          & ! Flag for deposition aqueous phase
&               ,nIO               & ! Restart  flag
&               ,mJacIO            & ! Flag for setting sparse Jacobian structure
&               ,MinLwcSet         & ! Set LWC allways to LWCLevel  (1=on, 0=off)
&               ,pHSet             & ! Initial pH by charge balance (1=on, 0=off)
                ,constLWC          & ! with cloud constLWC>=1
&               ,mAcoeff           & ! model for mixed activity coefficients 
&               ,mPitz             & ! determination Pitzer activity coefficients
&               ,mLR               & ! determination Long Range activity coefficients(AIMOFAC)
&               ,mAMR              & ! determination middle range activity coefficients (MR in AIOMFAC)
&               ,mMR               & ! determination middle range activity coefficients (MR in LIFAC)
&               ,mUni              & ! determination UNIFAC activity coefficients
&               ,IowMode           & ! consideration of Ion-Organics-Water interactions
&               ,rqMode            & ! UNIFAC: Setting of RK and QK for ions
&               ,LRMode            & ! LongRange : Setting to take Input
&               ,mofac             & ! output of activity coefficients for mol fractions
&               ,Ladebalken        & ! ladebalken im terminal bei simulation (=1, default=0)
&               ,Error_Est         & ! error estimation 1 = inf norm  , 2 = euklid norm
&               ,ErrorLog            ! if = 0 do not print error log 
!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!--- Times:  < 0 in seconds, > 0 in hours, = 0  meteorology
REAL(RealKind) :: &
&                DtCoupl           & ! Coupling time step
&               ,tAnf              & ! Model start time
&               ,tEnd              & ! Model end time
&               ,tSimul            & ! Real start time
&               ,dStart            & ! Period of start phase
&               ,dPostRun          & ! Period for postrun
&               ,ErrTime           & ! Restart time for error treating
&               ,StpAVS            & ! Time step for AVS output
&               ,StpNetcdf         & ! Time step for Netcdf output
&               ,StpDiagNcdf       & ! Time step for diagnose Netcdf output
&               ,StpOut            & ! Time step for ASCII output
&               ,StpDiag             ! Time step for diagnose
!--- Checkpoints
     INTEGER :: PtrCp               ! Current checkpoint position 
     INTEGER :: nCheck
     INTEGER, ALLOCATABLE :: iCheck(:)
     REAL(RealKind), ALLOCATABLE :: tCheck(:)

!--- Current day for several days runs 
      INTEGER :: ntag

!---  Photolysis
      INTEGER :: iDate               ! Current date
      REAL(RealKind) :: rlat, rlon          ! Latitude, Longitude


!--  dust factor for damping of photolysis rates, measured JNO2
      REAL(RealKind) :: Dust = 1.0d0
      REAL(RealKind) :: JNO2          
      REAL(RealKind) :: minStp = 1.0d-25
      REAL(RealKind) :: maxStp = 500.0d0

!-----------------------------------------------------------------
!---  Meteorology
!-----------------------------------------------------------------
!
!--- Control Parameter
      INTEGER :: mpMod             & ! Microphysical model
&               ,mCase             & ! Microphysical scenario
&               ,mMicPhys          & ! Microphysic on/off
&               ,ItemEnd           & ! Last trajectory section
&               ,nfRaoult          & ! Flag for feedback in Raoult term
&               ,nfChem            & ! Flag for chemistry for surfacetension 
&               ,nfMass            & ! Flag for aerosol mass feedback
&               ,stModel             ! Flag for the Surface tension parameterization

!--- Microphysical control parameters
      REAL(RealKind) :: StartWind,        & ! Wind speed in Goldlauter
                 VdepDrop,         & ! Deposition velocity of drops
                 FixedEpsi           ! Fixed solutable aerosol part
                                     ! EPSILON 
!
!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!
!--- Control Parameter
      INTEGER :: ITolMod        & ! Setting of BDF tolerances
&               ,MaxOrdROW      & ! Maximum order of ROW scheme
&               ,LinSolv        & ! Flag for linear algebra approach
&               ,ImpEuler
!
!--- Tolerances for ROW Scheme
      REAL(RealKind) :: &
&                RtolROW        & ! Relative tolerance for ROW method
&               ,AtolGas        & ! Absolute tolerance for gas phase
&               ,AtolAqua         ! Absolute tolerance for liquid phase

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
      INTEGER :: FactorisationStrategie = 0  ! 
      INTEGER :: SolLinSystemStrategie = 0  ! 
      
!-----------------------------------------------------------------
!---  Entrainment Control 
!-----------------------------------------------------------------
!
!---   Flags to control entrainment
    INTEGER ::  ExMet  = 0      & ! temperature           (0=off, 1=read)
&              ,ExLwc  = 0      & ! total humitity        (0=off, 1=read)
&              ,ExGas  = 0      & ! gas phase             (0=off, 1=read, 2=activation, 3=initial)
&              ,ExPart = 0      & ! particle distribution (0=off, 1=read, 2=activation, 3=initial)
&              ,ExComp = 0        ! particle composition  (0=off, 1=read, 2=activation, 3=initial)

!---   Entrainment parameter 
    REAL(RealKind) :: &
&              MuFac   = 1.0d0    & ! Factor for modification of all entrainment parameters
&             ,MuMet0  = 1.0d-4   & ! temperature
&             ,MuLwc0  = 1.0d-3   & ! total humitity
&             ,MuGas0  = 1.0d-4   & ! gas phase
&             ,MuPart0 = 1.0d-5     ! particle distribution

    CHARACTER(80) :: ExFile

!-----------------------------------------------------------------
!---  Set Constants and Unit Conversion
!-----------------------------------------------------------------
!
!---  constants
    REAL(RealKind), PARAMETER :: mol2part  = 6.02295d17     &
&                        ,GasConst_R = 0.082056d0    &   ! [in l*atm/mol/K]
&                        ,hour       = 3600.d0       &
&                        ,secday     = 4.32d04      &
&                        ,ZERO       = 0.d0          &
&                        ,ONE        = 1.d0          &
&                        ,Pi         = 4.0d0*ATAN(1.0d0) &
&                        ,DR         = Pi / 180.d0  &
&                        ,PiHalf     = 2.0d0*ATAN(1.0d0) &
&                        ,eps        = EPSILON(1.0d0) &
&                        ,epsY       = 1.0d-7

!--- Unit Conversion
    INTEGER :: GasUnit, AquaUnit, GasRateUnit
    REAL(RealKind) :: GasFac

    REAL(RealKind) :: ConvGas, ConvAir

!--- Activation control
    INTEGER :: nact, nact_Old
    REAL(RealKind) :: ActLevel

!-----------------------------------------------------------------
!---  Output control 
!-----------------------------------------------------------------

!--- Checkpoint control
    REAL(RealKind), PARAMETER :: tEps = 5.0d-04  & ! Time tolerance for output control
&                        ,tBig = 1.0d30     ! Infinity time

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

!---  Flags for diagnose and ASCII-Output
    INTEGER, ALLOCATABLE :: FracDiag(:,:)
    INTEGER, ALLOCATABLE :: FracOut(:,:)

    !
    ! 
    REAL(RealKind) :: sqrtNSPC

 END MODULE mo_control

