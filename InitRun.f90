   SUBROUTINE InitRun(RunFile)
!==================================================
!===  Reading and Setting of Run Control Parameters
!==================================================
!
      !USE mo_bdf
      USE Kind_Mod
      USE mo_micfys
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
&               ,QtListFile        & ! File with listed species for Qt-Output 
&               ,AVSFile           & ! AVS   output file
&               ,OutFile           & ! ASCII output file
&               ,DiagFile          & ! Diagnose file
&               ,StatFile          & ! Statistics file
&               ,JacFile           & ! Sparse Jacobian File
&               ,ErrFile           & ! Output error file
&               ,NetcdfFileName    & ! Output NetCdfFile
!
!--- INTEGER : Unit Numbers
&               ,MetUnit           & ! Meteorology file
&               ,ChemUnit          & ! Chemical mechanism
&               ,InitUnit          & ! Initial concentrations
&               ,DataUnit          & ! Initial concentrations
&               ,AVSUnit           & ! AVS   output file
&               ,OutUnit           & ! ASCII output file
&               ,DiagUnit          & ! Diagnose file
&               ,StatUnit          & ! Statistics file
&               ,JacUnit           & ! File for sparse Jacobian
&               ,ErrUnit           & ! Output error file
!
!-- REAL(RealKind) : Set Levels and Parameters for Processes
&               ,LwcLevelmin       & ! lower Level for LWC      {l/m3]
&               ,LwcLevelmax       & ! upper Level for LWC      {l/m3]
&               ,constLWC          & ! wähle konstanten lwc wert (=0) oder lineare funktion (/=0)
&               ,DryLevel          & ! Level for dry mass [g/m3]
&               ,MinScav           & ! Minimum drop radius for scavening [µm]
&               ,MinChem           & ! Minimum drop radius for chemistry [µm]
&               ,LwcAll            & ! Liqid Water Content
&               ,MeanRad           & ! Mean Droplet Radius

!-- LOGICAL :: Print system matrices
&               ,MatrixPrint       &
&               ,NetCdfPrint       &

!--- INTEGER : Dimensions and Resolution
&               ,ntfrac            & ! Number of fractions
&               ,ResFrac           & ! Resolution of chemistry spectrum
&               ,FirstChem         & ! First fraction for chemistry
&               ,LastChem          & ! Last  fraction for chemistry
&               ,MaxFrac           & ! Maximum fraction for chemistry
&               ,nReso             & ! Resolution for smoothing
&               ,nSmooth           & ! Number of fine fractions for smoothing
&               ,nCell             & ! Number of grid cells

!--- INTEGER : Control Parameter
&               ,nCouple           & ! Coupling flag
&               ,nDep              & ! Flag for deposition gas phase
&               ,nDepAqua          & ! Flag for deposition aqueous phase
&               ,nIO               & ! Restart  flag
&               ,mJacIO            & ! Flag for setting sparse Jacobian
&               ,MinLwcSet         & ! Set all LWC to LWCLevel      (1=on, 0=off)
&               ,pHSet             & ! Initial pH by charge balance (1=on, 0=off)
&               ,mAcoeff           & ! model for mixed activity coefficients
&               ,mPitz             & ! determination Pitzer activity coefficients
&               ,IowMode           & ! consideration of Ion-Organics-Water interactions
&               ,mUni              & ! determination UNIFAC activity coefficients
&               ,mLR               & ! determination Long Range activity coefficients (LR in LIFAC & AIOMFAC)
&               ,mMR               & ! determination Middle Range activity coefficients (MR in LIFAC)
&               ,mAMR              & ! determination middle range activity coefficients (MR in AIOMFAC)
&               ,rqMode            & ! UNIFAC: Setting of RK and QK for ions
&               ,mofac             & ! output of activity coefficients for mol fraction
&               ,Ladebalken        &
&               ,Error_Est         &
&               ,ErrorLog          &

!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!
!--- REAL(8): Times ( < 0 in seconds, > 0 in hours, = 0  meteorology)
&               ,DtCoupl           & ! Coupling time step
&               ,tAnf              & ! Model start time
&               ,tEnd              & ! model end time
&               ,tSimul            & ! Real start time
&               ,dStart            & ! Period of start phase
&               ,dPostRun          & ! Period of postrun
&               ,ErrTime           & ! Restart time for error treating
&               ,StpAVS            & ! Time step for AVS output
&               ,StpOut            & ! Time step for ASCII output
&               ,StpNetcdf         & ! Time step for Netcdf output 
&               ,StpDiagNcdf       & ! Time step for diagnose Netcdf output 
&               ,StpDiag           & ! Time step for diagnose 
&               ,StpDiag           & ! Time step for diagnose 
!
!--- INTEGER : nCheck, iCheck
!--- REAL(RealKind) : tCheck
&               ,nCheck            & ! Number of checkpoints  
&               ,iCheck            & ! Meaning of checkpoints        
&               ,tCheck            & ! Time of checkpoints

!---  Photolysis
&               ,idate             & ! Date: yymmdd
&               ,rlat              & ! latitude  [rad]
&               ,rlon              & ! longitude [rad]
&               ,dust              & ! dust factor (damping of photolysis)
&               ,JNO2              & ! measured NO2 photolysis rate [1/s]
&               ,minStp            & ! minimal time step size
&               ,maxStp            & ! maximal time step size
!
!-----------------------------------------------------------------
!---  Meteorology
!-----------------------------------------------------------------
!
!--- INTEGER : Control Parameter
&               ,mpMod             & ! microphysical model
&               ,mCase             & ! microphysical scenario
&               ,mMicPhys          & ! Microphysic on/off
&               ,nfRaoult          & ! Flag for feedback in Raoult term
&               ,nfMass            & ! Flag for aerosol mass feedback
&               ,stModel           & ! Surface tension model 
&               ,nfChem            & ! Flag for chemistry for surface tension calculations
&               ,ItemEnd           & ! trajectory stop section

!--- REAL(RealKind) : Microphysical control parameters
&               ,StartWind         & ! wind speed in Goldlauter
&               ,VdepDrop          & ! Deposition velocity of drops
&               ,FixedEpsi         & ! Fixed solutable aerosol part
                                     ! EPSILON
!
!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!
!--- INTEGER : Control Parameter
&               ,ITolMod           & ! Setting of BDF tolerances 
&               ,MaxOrdROW         & ! Maximum order of BDF scheme
&               ,LinSolv           & ! Flag for linear algebra approach
&               ,ImpEuler          & ! implicit euler 

!--- REAL(RealKind) : Tolerances for ROW Scheme
&               ,RtolROW           & ! Relative tolerance For ROW
&               ,AtolGas           & ! Absolute tolerance for gas phase
&               ,AtolAqua          & ! Absolute tolerance for liquid phase
&               ,PI_StepSize       & ! logical for pi stepsize control

!--- CHARACTER(2) : Control Parameter
&               ,solveLA           & ! how to solve linear algebra 'cl' or 'ex'
&               ,RosenbrockMethod  & ! choose ROW scheme  (integer 1,...,15)
!
!-----------------------------------------------------------------
!---  Linear Algebra
!-----------------------------------------------------------------
!
!--- Control Parameter
&               ,OrderingStrategie      & ! MUMPS ordering strategie
&               ,ParOrdering            & ! MUMPS ordering type for parallel symbolic phase
&               ,FactorisationStrategie & ! MUMPS factorisation / own factorisation
&               ,SolLinSystemStrategie    &
   
                    
!-----------------------------------------------------------------
!---  Entrainment Control
!-----------------------------------------------------------------
!                
!--- INTEGER :: Flags to control entrainment
&               ,ExMet         & ! temperature           (0=off, 1=read)
&               ,ExLwc         & ! total humitity        (0=off, 1=read)
&               ,ExGas         & ! gas phase             (0=off, 1=read, 2=activation, 3=initial)
&               ,ExPart        & ! particle distribution (0=off, 1=read, 2=activation, 3=initial)
&               ,ExComp        & ! particle composition  (0=off, 1=read, 2=activation, 3=initial)
!
!--- REAL(RealKind) :: Entrainment parameter 
&               ,MuFac         & ! Factor for modification of all entrainment parameters
&               ,MuMet0        & ! temperature
&               ,MuLwc0        & ! total humitity
&               ,MuGas0        & ! gas phase
&               ,MuPart0       & ! particle distribution
!
!--- CHARACTER(80) :: Entrainment file
                ,ExFile        
!
!-----------------------------------------------------------------
      IMPLICIT NONE
      
      CHARACTER(80) :: RunFile

!--- Internal parameters
      INTEGER :: i
      INTEGER :: iHelp(11)
      INTEGER :: Ipoint1, Ipoint2, Ipoint3, Ipoint4, Ipoint5
      
      INTEGER :: maxord
!
      REAL(RealKind) :: tHelp(11)
      REAL(RealKind) :: Cpoint1, Cpoint2, Cpoint3, Cpoint4, Cpoint5
      CHARACTER(2) :: iFrac_String 
!
!-----------------------------------------------------------------
!--- NAMELISTS
      NAMELIST /SCENARIO/  Bsp,                                          &
&               LwcLevelmin, LwcLevelmax, DryLevel, MinScav, MinChem, LwcAll, MeanRad,   &
&               ntfrac, nCouple, ResFrac, nDep, nDepAqua, Ladebalken,    &  
&               FirstChem, MaxFrac, nIO, mJacIO, MinLwcSet, pHSet,       &
&               mAcoeff, mUni, mPitz, nReso, nSmooth, mMR, constLWC,     &
&               IowMode, rqMode, mofac, mAMR, mLR, ErrorLog, MatrixPrint, NetCdfPrint
!
      NAMELIST /FILES/  MetFile, ChemFile, InitFile, DataFile, AVSFile,  &
&               OutFile, DiagFile, StatFile, ErrFile, JacFile,           &
&               QtListFile,NetcdfFileName,                               &
&               MetUnit, ChemUnit, InitUnit, DataUnit, AVSUnit, OutUnit, &
&               DiagUnit,StatUnit, ErrUnit, JacUnit
!
      NAMELIST /METEO/ mpMod, mCase, mMicPhys, nfRaoult, nfMass, nfChem, &
&               stModel, ItemEnd, StartWind, VdepDrop, FixedEpsi

      NAMELIST /TIMES/ DtCoupl, tAnf, tEnd, tSimul, dStart, StpAVS,      &
&               StpOut, StpDiag, ErrTime, idate, rlat, rlon, Dust, JNO2, &
&               StpNetcdf, StpDiagNcdf,                                  &
&               dPostRun, Cpoint1, Cpoint2, Cpoint3, Cpoint4, Cpoint5,   &
&               Ipoint1, Ipoint2, Ipoint3, Ipoint4, Ipoint5, minStp, maxStp
!
      NAMELIST /NUMERICS/  RtolROW, AtolGas, AtolAqua, PI_StepSize, solveLA,          &
&               ITolMod, MaxOrdROW, LinSolv, RosenbrockMethod, ImpEuler,Error_Est
!
      NAMELIST /ORDERING/  OrderingStrategie, ParOrdering,               &
&              FactorisationStrategie, SolLinSystemStrategie
!      
      NAMELIST /ENTRAIN/  ExFile,                                        &
&               MuFac, MuMet0, MuLwc0, MuGas0, MuPart0,                  &
&               ExMet, ExLwc, ExGas, ExPart, ExComp


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
      DryLevel = 1.0d-30       ! Lower level for dry mass [g/m3]
      MinScav  = 0.1d0        ! Minimum drop radius for scavening [µm]
      MinChem  = 0.1d0        ! Minimum drop radius for chemistry [µm]
      LwcAll   = -1.0d0        ! Liqid Water Content (Dummy)
      MeanRad  = 10.0d0        ! Mean Droplet Radius [µm]

!--- INTEGER : Dimensions and Resolution
      ntfrac    = JMAX       ! Number of fractions
      ResFrac   = 1          ! Resolution of chemistry spectrum
      FirstChem = 1          ! First fraction for chemistry
      LastChem  = JMAX       ! Last  fraction for chemistry
      MaxFrac   = JMAX       ! Maximum fraction for chemistry
      nSmooth   = 1          ! Number of fine fractions for smoothing
      nReso     = 2          ! Resolution for smoothing   
      nCell     = 1          ! Number of grid cells

!--- INTEGER : Control Parameter
      nCouple   = 2          ! Coupling flag
      nDep      = 0          ! Flag for deposition gas phase
      nDepAqua  = 1          ! Flag for deposition aqueous phase
      nIO       = 0          ! Restart  flag
      mJacIO    = 0          ! Setting of sparse Jacobian   (0=calculate, 1=store, 2=read)
      MinLwcSet = 0          ! Set all LWC to LWCLevel      (1=on, 0=off)
      pHSet     = 1          ! Initial pH by charge balance (1=on, 0=off)
      mAcoeff   = 0          ! model for mixed activity coefficients
      mPitz     = 0          ! Approach for activity coefficients (0=off,1=Harvie,2=Pitzer)
      mUni      = 0          ! determination UNIFAC activity coefficients
      mLR       = 0          ! determination  Long Range activity coefficients
      mMR       = 0          ! determination LIFAC Middle Range activity coefficients
      mAMR      = 0          ! determination AIOMFAC Middle Range activity coefficients
      IowMode   = 1          ! consideration of Ion-Organics-Water interactions
      rqMode    = 0          ! UNIFAC: Setting of RK and QK for ions
      mofac     = 1          ! output of activity coefficients for mol fraction
      Error_Est = 2          ! default for error estimation is euklid norm
      ErrorLog  = 0 
      MatrixPrint = .FALSE.       ! 0 = Print no matrix, 0 /= Print all matrices -> no simulation
      NetCdfPrint = .FALSE.
!--- Read SCENARIO namelist
      READ(15,SCENARIO)

!--- Reset wrong values
      nCouple = MIN(ntFrac,nCouple)
      ntFrac  = 1 + (MaxFrac-1)/ResFrac

!-----------------------------------------------------------------
!---  Files
!-----------------------------------------------------------------
!
!--- CHARACTER(80) : Files
      MetFile    = 'MET/initial'             ! Meteorology file
      ChemFile   = 'CHEM/'//TRIM(RunFile)//'.sys'    ! Chemical mechanism
      DataFile   = 'CHEM/'//TRIM(RunFile)//'.dat'    ! Gas and aqueous phase data
      JacFile    = 'CHEM'//TRIM(RunFile)//'.str'     ! Sparse Jacobian File
      InitFile   = 'INI/'//TRIM(RunFile)//'.ini'     ! Initial concentrations
      QtListFile = 'WITHOUT'                         ! File with listed species for Qt-Output 

      AVSFile  = 'AVS/'//TRIM(Bsp)                 ! AVS   output file
      NetcdfFileName = TRIM(Bsp)//'.nc'            ! Netcdf output file
      OutFile  = 'OUT_'//TRIM(Bsp)                 ! ASCII output file
      DiagFile = 'Diag_'//TRIM(Bsp)                ! Diagnose file
      StatFile = 'Stat_'//TRIM(Bsp)                ! Statistics file
      ErrFile  = 'INPUT/'//TRIM(Bsp)//'.err'       ! Input/Output error file

      !CALL SYSTEM('mkdir DIAG/'//TRIM(Bsp))
      !CALL SYSTEM('mkdir DIAG/!'//TRIM(Bsp))
!
!--- INTEGER : Unit Numbers
      MetUnit  = 15           ! Meteorology file
      ChemUnit = 10           ! Chemical mechanism
      InitUnit = 12           ! Initial concentrations
      DataUnit = 13           ! Gas and aqueous phase data
      AVSUnit  = 21           ! AVS   output file
      OutUnit  = 39           ! ASCII output file
      DiagUnit = 41           ! Diagnose file
      StatUnit = 38           ! Statistics file
      JacUnit  = 31           ! Sparse Jacobian File
      ErrUnit  = 32           ! Output error file
      
!--- Read FILES namelist
      READ(15,FILES)
      
!--- Adjust Filenames
      Bsp  = ADJUSTL(Bsp)

      MetFile  = ADJUSTL(MetFile)
      ChemFile = ADJUSTL(ChemFile)
      DataFile = ADJUSTL(DataFile)
      InitFile = ADJUSTL(InitFile)
      ErrFile  = ADJUSTL(ErrFile)

      AVSFile  = ADJUSTL(AVSFile)
      NetcdfFileName=ADJUSTL(NetcdfFileName)
      OutFile  = ADJUSTL(OutFile)
      DiagFile = ADJUSTL(DiagFile)
      StatFile = ADJUSTL(StatFile)
      JacFile  = ADJUSTL(JacFile)

!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!
!--- REAL(8): Times ( < 0 in seconds, > 0 in hours, = 0  meteorology)
      tAnf    =  0.d0        ! Model start time    [in h}
      tEnd    =  0.d0        ! model end time      [in h}
      tSimul  =  0.d0        ! Real start time     [in h]  

!--- REAL(8): Times in seconds.
      DtCoupl  = 10.d0        ! Coupling time step            [in sec]
      dStart   =  0.d0        ! Period of start phase         [in sec]
      dPostRun =  0.d0        ! Period of post run            [in sec]
      StpAVS   = 10.d0        ! Time step for AVS output      [in sec]
      StpNetcdf   = 0.d0      ! Time step for Netcdf output      [in sec]
      StpDiagNcdf = 0.d0      ! Time step for diagnose Netcdf output [in sec]
      StpOut   =  0.d0        ! Time step for ASCII output    [in sec]
      StpDiag  =  0.d0        ! Time step for diagnose        [in sec]
      ErrTime  =  1.d30       ! Input time for error treating [in sec]

!--- INTEGER: Checkpoint meaning
      Ipoint1 = 0             ! Control check points
      Ipoint2 = 0         
      Ipoint3 = 0         
      Ipoint4 = 0         
      Ipoint5 = 0         

!--- REAL(8): Checkpoint times
      Cpoint1 = 0.d0          ! Control check points
      Cpoint2 = 0.d0 
      Cpoint3 = 0.d0 
      Cpoint4 = 0.d0 
      Cpoint5 = 0.d0 

!---  Photolysis (Here: FEBUKO chemistry-case I)
      idate = 011027          ! Date: yymmdd  (21.June 2001)
      rlat  = 5.065d+01       ! latitude  [grad] (Schmuecke)
      rlon  = 1.077d+01       ! longitude [grad] (Schmuecke)
      Dust  = 1.0d0           ! dust factor (damping of photolysis)
      JNO2  = -1.0d0           ! measured NO2 photolysis rate [1/s]

!--- Read TIMES namelist
      READ(15,TIMES)

      IF (StpDiagNcdf>0.d0) THEN
         DO i=1,MaxFrac
            WRITE(iFrac_String,'(I2)') i 
            CALL SYSTEM('mkdir DIAG/'//TRIM(Bsp)//'/Frac_'//TRIM(ADJUSTL(iFrac_String)))
         END DO
      END IF
!
!--- Reset wrong walues
      IF (StpAVS>0.0e0) StpAVS  = MAX(StpAVS,DtCoupl)
      IF (StpNetcdf>0.0e0) StpNetcdf = MAX(StpNetcdf,DtCoupl)
      IF (StpOut>0.0e0) StpOut  = MAX(StpOut,DtCoupl)
      IF (StpDiag>0.0e0) StpDiag = MAX(StpDiag,DtCoupl)
      IF (StpDiagNcdf>0.0e0) StpDiagNcdf = MAX(StpDiagNcdf,DtCoupl)

!-----------------------------------------------------------------
!--- Set Checkpoints
      nCheck = 0
      IF (dStart > 0.0d0)  THEN   
         nCheck = nCheck + 1                ! Initialization phase
         iHelp(ncheck) = 1
         tHelp(ncheck) = 0.d0
      END IF
!
      IF (Cpoint1 > 0.0d0)  THEN             ! Control check points
         nCheck = nCheck + 1  
         iHelp(ncheck) = Ipoint1
         tHelp(ncheck) = Cpoint1
      ELSE IF (Ipoint1 > 0)  THEN      
         nCheck = nCheck + 1  
         iHelp(ncheck) = - Ipoint1 
         tHelp(ncheck) = 1.d20
      END IF
      IF (Cpoint2 > 0.0d0)  THEN             ! Control check points
         nCheck = nCheck + 1  
         iHelp(ncheck) = Ipoint2
         tHelp(ncheck) = Cpoint2
      ELSE IF (Ipoint2 > 0)  THEN      
         nCheck = nCheck + 1  
         iHelp(ncheck) = - Ipoint2 
         tHelp(ncheck) = 1.d20
      END IF
      IF (Cpoint3 > 0.0d0)  THEN             ! Control check points
         nCheck = nCheck + 1  
         iHelp(ncheck) = Ipoint3 
         tHelp(ncheck) = Cpoint3
      ELSE IF (Ipoint3 > 0)  THEN      
         nCheck = nCheck + 1  
         iHelp(ncheck) = - Ipoint3 
         tHelp(ncheck) = 1.0d20
      END IF
      IF (Cpoint4 > 0.0d0)  THEN             ! Control check points
         nCheck = nCheck + 1  
         iHelp(ncheck) = Ipoint4 
         tHelp(ncheck) = Cpoint4
      ELSE IF (Ipoint4 > 0)  THEN      
         nCheck = nCheck + 1  
         iHelp(ncheck) = - Ipoint4 
         tHelp(ncheck) = 1.d20
      END IF
      IF (Cpoint5 > 0.0d0)  THEN             ! Control check points
         nCheck = nCheck + 1  
         tHelp(ncheck) = Ipoint5 
         tHelp(ncheck) = Cpoint5
      ELSE IF (Ipoint5 > 0)  THEN      
         nCheck = nCheck + 1  
         iHelp(ncheck) = - Ipoint5 
         tHelp(ncheck) = 1.0d20
      END IF
!
      IF (dPostRun > 0.0d0)  THEN   
         nCheck = nCheck + 1                ! Postrun phase
         iHelp(ncheck) = 10
         tHelp(ncheck) = 1.0d20
      END IF
!
      nCheck = nCheck + 1                   ! for Infinity
      iHelp(nCheck) =  9
      tHelp(nCheck) = 1.0d20
!
!---  Allocation and copying of checkpoint arrays
      ALLOCATE (iCheck(nCheck))
      ALLOCATE (tCheck(nCheck))
!
!---  Store checkpoints
      iCheck(:) = iHelp(1:nCheck)
      tCheck(:) = tHelp(1:nCheck)
      
!-----------------------------------------------------------------
!---  Meteorology
!-----------------------------------------------------------------
!
!--- INTEGER : Control Parameter
      mpMod    = 1            ! 1D, Fixed grid
      mCase    = 0            ! FEBUKO scenario
      nfRaoult = 0            ! Feedback in Raoult term (off)
      nfChem   = 0            ! Chemistry for surfacetension (off)
      nfMass   = 0            ! Feedback for aerosol mass (off)
      ItemEnd  = 0            ! trajectory stop section
      stModel  = 1            ! Surface tension model: Set to water surface tension

!--- REAL(RealKind) : Microphysical control parameters
      FixedEpsi =  0.d0      ! Fixed soluble aerosol part
      VdepDrop  =  0.d0      ! Drop deposition velocity
      StartWind = -1.d0      ! wind speed in Goldlauter
!                             (take measured value)

!--- Read METEO namelist
      !READ(15,METEO)
                            
!--- Set default values for moving run
      IF (mpMod == 2)  THEN
         nCouple = MIN(1,nCouple)
      END IF
!
!---  Check feedback values
      IF (mpMod == 1 .AND. nfMass == 1) THEN
         WRITE(*,*)  'PT1 ', 'Mass Feedback not implemented in Fixed Mode:'  &
&                          , ' mfMass = 0  is used!'
         nfMass = 0
      END IF

!--- Set default trajectory stop section
      IF (ItemEnd == 0)  THEN
         IF (mCase == 0)  THEN
            ItemEnd = 22
         ELSE IF (mCase == 1)  THEN
            ItemEnd = 4
         ELSE 
            ItemEnd = 1.d5
         END IF
      END IF

!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!
!--- INTEGER : Control Parameter
      ITolMod   = 1          ! Setting of BDF tolerances 
      MaxOrdROW = 3          ! Maximum order of ROW scheme
      LinSolv   = 3          ! Flag for linear algebra approach

!--- REAL(RealKind) : Tolerances for ROW Scheme
      RtolROW  = 1.0d-5               ! Relative tolerance For ROW
      AtolGas  = 1.0d-7               ! Absolute tolerance for gas phase
      AtolAqua = 1.0d-7               ! Absolute tolerance for liquid phase
      PI_StepSize = .FALSE.
      solveLA  = 'ex'                 ! method of solving linear algebra
      RosenbrockMethod  = 'ROS34PW3'  ! ROW scheme
      ImpEuler = 0                    ! 1 for implicit euler integration
      
!
!--- Read NUMERICS namelist
      READ(15,NUMERICS)
!
!--- Set NUMERICS values
      maxord = MaxOrdROW
      
!
!-----------------------------------------------------------------
!---  Linear Algebra
!-----------------------------------------------------------------
!
!--- Control Parameter
      READ(15,ORDERING)

      !OrderingStrategie=7
      !WRITE(*,*) '    ORDERINGSTRATEGIE FIXED TO MUMPS = 7'
!-----------------------------------------------------------------
!---  Entrainment Control
!-----------------------------------------------------------------
!                
!--- INTEGER :: Flags to control entrainment
      ExMet   = 0          ! temperature           (0=off, 1=read)
      ExLwc   = 0          ! total humitity        (0=off, 1=read)
      ExGas   = 0          ! gas phase             (0=off, 1=read, 2=activation, 3=initial)
      ExPart  = 0          ! particle distribution (0=off, 1=read, 2=activation, 3=initial)
      ExComp  = 0          ! particle composition  (0=off, 1=read, 2=activation, 3=initial)
!
!--- REAL(RealKind) :: Entrainment parameter 
      MuFac   = 1.0d0       ! Factor for modification of all entrainment parameters
      MuMet0  = 1.0d-4      ! temperature
      MuLwc0  = 1.0d-3      ! total humitity
      MuGas0  = 1.0d-4      ! gas phase
      MuPart0 = 1.0d-5      ! particle distribution
!
!--- CHARACTER(80) :: Entrainment file 
      ExFile = 'MET/entrain.dat'    
!
!--- Read ENTRAIN namelist
      READ(15,ENTRAIN)
!
!--- Set entrainment file 
      ExFile = ADJUSTL(ExFile)
!
      CLOSE(15)

!==================================================
!===  Set Values
!==================================================
!
!---  Transformation of times  (hours ==> seconds)
    !  IF (tAnf < 0.0d0) THEN
    !     tAnf = ABS(tAnf) 
    !  ELSE IF (tAnf > 0.0d0) THEN
    !     tAnf = 3600.0d0 * tAnf 
    !  END IF

    !  IF (tEnd < 0.0d0) THEN
    !     tEnd = ABS(tEnd) 
    !  ELSE IF (tAnf > 0.0d0) THEN
    !     tEnd = 3600.0d0 * tEnd 
    !  END IF

    !  IF (tSimul < 0.0d0) THEN
    !     tSimul = ABS(tSimul) 
    !  ELSE IF (tSimul > 0.0d0) THEN
    !     tSimul = 3600.e0 * tSimul 
    !  END IF
!
!==================================================
!===  Write control values
!==================================================
      !
      !---  File names

      


!----------------------------------------------------------------------------------------
   END SUBROUTINE InitRun
