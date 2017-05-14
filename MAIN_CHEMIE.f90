!=======================================
!=======================================
! vorlaeufiges Hauptprogramm fuer Tests
!=======================================
!=======================================
!
!
PROGRAM MAIN_CHEMIE
  !
  USE Sparse_Mod
  USE Kind_Mod
  USE Chemsys_Mod
  USE Integration_Mod
  USE Rosenbrock_Mod
  USE mo_control
  USE mo_unirnk
  USE mo_reac, ONLY: nspc, MW, rMW
  USE mo_MPI
  USE mo_IO
  USE mo_ckinput
  !USE Euler_Mod
  USE NetCDF_Mod
  IMPLICIT NONE
  !
  CHARACTER(80)   :: Filename0 = ''        ! *.run file
  !
  REAL(RealKind)  :: Atol(2)
  INTEGER         :: i,j,nSchwefel
  INTEGER         :: io_err,  STAT
  !
  ! use other ROS methode or tolerance
  CHARACTER(80)   :: neuTolR  = ''
  CHARACTER(80)   :: neuTolA  = ''
  CHARACTER(80)   :: neuROW   = ''

  ! temp variables for molecular weights
  CHARACTER(16)   :: tmpchar0 = '-'
  REAL(RealKind)  :: tmpMW0
  INTEGER         :: tmpPos, tmpCnt

  ! convertion from mole to mass to conc
  REAL(RealKind), ALLOCATABLE   :: MoleFrac(:), MassFrac(:), MoleConc(:)
  REAL(RealKind)                :: Press_in_dyncm2
  !
  !================================================================
  !===                     MAIN Programm
  !================================================================

  !----------------------------------------------------------------
  ! ---  Initialize MPI library 
  CALL StartMPI

  !----------------------------------------------------------------
  ! --- Print the program logo  
  CALL Logo

  !----------------------------------------------------------------
  ! --- Read run control parameters (which runfile)
  CALL getarg( 1 , FileName0 )             
  IF ( FileName0 == '' ) THEN
    WRITE(*,*) 'Input RUNFilename:'
    READ(*,*)   FileName0
  END IF
  FileName0 = TRIM(ADJUSTL(FileName0))

  !================================================================
  !===                     Initialization
  !================================================================
  !
  !----------------------------------------------------------------
  ! --- Read run-file
  CALL InitRun( FileName0 )
  IF ( MPI_ID == 0 ) WRITE(*,*) '  Initialize run-file .......... done'
  !
  !----------------------------------------------------------------
  ! --- Initialize all reaction types
  !CALL InitNReacType()
  !
  Tspan = (/ tAnf , tEnd /)
 
  !----------------------------------------------------------------
  ! --- set cloud intervall
  LWCBounds(1) = tAnf * HOUR
  LWCBounds(2) = LWCBounds(1) + 1.00d0  * HOUR
  LWCBounds(3) = LWCBounds(2) + 0.25d0  * HOUR
  LWCBounds(4) = LWCBounds(3) + 9.50d0  * HOUR
  LWCBounds(5) = LWCBounds(4) + 0.25d0  * HOUR
  LWCBounds(6) = LWCBounds(5) + 1.00d0  * HOUR

  !----------------------------------------------------------------
  !  --- read the .sys data, save coefs in sparse matrix
  Time_Read = MPI_WTIME()
  OPEN ( UNIT=89 , FILE=ADJUSTL(TRIM(ChemFile))//'.chem' , STATUS='UNKNOWN' )
  IF ( MPI_ID == 0 ) WRITE(*,'(A38)',ADVANCE='NO') '  Reading sys-file .............     '

  IF ( TempEq.AND.ChemKin ) THEN
    
    IF ( MPI_ID == 0 ) WRITE(*,*) ' ---->  Solve Gas Energy Equation '

    CALL Read_Elements    ( ChemFile    , 969 )
    CALL Read_Species     ( ChemFile    , 969 )
    CALL Read_Reaction    ( ChemFile    , 969 )

    CALL PrintHeadSpecies   ( ChemFile    , 89 ) 
    CALL PrintSpecies       ( ListGas2    , 89 )
    CALL PrintHeadReactions ( 89 )
    CALL PrintReactions     ( ReactionSystem , 89 , .TRUE. )
    CALL PrintFinalReactions( 89 )

    CALL GetSpeciesNames( ChemFile , y_name )
    CALL Read_ThermoData( SwitchTemp , DataFile , 696 , nspc )

    !--- richtigen index holen, da TB unsortiert eingelesen wurde
    CALL GatherTBindex

    IF (Vectorized) THEN
      CALL GatherReactionTypeIndex
    END IF
   
    !CALL PrintReactionSystem(ReactionSysteM)
    !--- Read molecular mass
   
    CALL Read_MolecularWeights(MW,MWeights,MWUnit,nspc)

    ALLOCATE( rMW(nspc) )  
    rMW = ZERO
    rMW = ONE / MW
    
    ALLOCATE( MoleFrac(ntGas)   , MassFrac(ntGas) )
    ALLOCATE( InitValAct(ntGas) , y_e(ntGas) )
    ALLOCATE( InitValKat(ntKat) )
    MoleFrac    = ZERO           ! mole fraction 
    MassFrac    = ZERO           ! mole fraction 
    InitValAct  = ZERO           ! mol/cm3 

    CALL Read_GASini    ( InitFile , MoleFrac , InitValKat )
    CALL Read_EMISS     ( InitFile , y_e )

    !Press = Pressure0               ! initial pressure in [Pa]
    Press_in_dyncm2 = Pressure0 * Pa_to_dyncm2

    MassFrac = MoleFr_To_MassFr( MoleFrac ) 
    
    MoleConc = MoleFr_To_MoleConc( MoleFrac,               &
                                 & Press = Press_in_dyncm2,&
                                 & Temp  = Temperature0    )
    
    ! Initialising reactor density
    rho  = Density( MoleConc )
    !rRho = kilo/rho       ! in [cm3/g]
    rRho = mega/rho       ! in [cm3/g]

    !InitValAct = MassFrac
    InitValAct = MoleConc


  ELSE
    IF ( MPI_ID==0 ) WRITE(*,*) ' ---->  Fix Temperature'
    CALL ReadSystem( ChemFile )
  
    !----------------------------------------------------------------
    ! ---  build the coeficient matrices and write .chem
    CALL PrintHeadSpecies ( ChemFile , 89 )

    IF ( ntGas   > 0 ) CALL PrintSpecies( ListGas2     , 89 )
    IF ( ntAqua  > 0 ) CALL PrintSpecies( ListAqua2    , 89 )
    IF ( ntSolid > 0 ) CALL PrintSpecies( ListSolid2   , 89 ) 
    IF ( ntPart  > 0 ) CALL PrintSpecies( ListPartic2  , 89 )
    IF ( ntKat   > 0 ) CALL PrintSpecies( ListNonReac2 , 89 )

    CALL PrintHeadReactions( 89 )
   
    !----------------------------------------------------------------
    ! --- Build the reaction system
    CALL AllListsToArray( ReactionSystem  &
    &                    , ListRGas       &
    &                    , ListRHenry     &
    &                    , ListRAqua      &
    &                    , ListRDiss      )
    !
    !----------------------------------------------------------------
    ! --- print reactions and build A, B and (B-A) structure
    CALL PrintReactions( ReactionSystem , 89 )
    CALL PrintFinalReactions( 89 )

    !----------------------------------------------------------------
    ! --- Input of initial data and thermodynamic properties
    CALL InputChemicalData( InitFile , DataFile , MetFile )
  END IF
  CLOSE(89)

  !----------------------------------------------------------------
  ! --- Read species for diagnose (print species concentrations to NetCDF file)
  CALL Read_Diag( OutNetcdfSpc , OutNetcdfPhase , InitFile )

  Time_Read = MPI_WTIME() - Time_Read  ! Stop Input timer

  IF (MPI_ID==0) WRITE(*,*) '  Read ini-file ................ done'
  IF (MPI_ID==0) WRITE(*,*) '  Print chem-file .............. done'

  !----------------------------------------------------------------
  ! --- Split the rate array in mpi_np parts 
  CALL BuildPartitions( MyParties , neq , MPI_np )
  
  !----------------------------------------------------------------
  ! --- this is for the new mass action product routine 
  CALL GatherSpeciesOrder(A)

  !-----------------------------------------------------------------------
  ! --- Start timer and set absolut tolerance for species
  Timer_Start = MPI_WTIME()
  
  !-----------------------------------------------------------------------
  ! --- Dimension initialisation for the unknowns and matrices

  nsr = nspc + neq

  IF ( TempEq ) THEN
    nDIM    = nspc + 1
    nDIMcl  = nspc + 1
    nDIMex  = nsr  + 1
    Atol    = (/ AtolGas , AtolTemp /)
  ELSE
    nDIM    = nspc
    nDIMcl  = nspc
    nDIMex  = nsr
    Atol    = (/ AtolGas , AtolAqua /)
  END IF

  rNspc = ONE/REAL(nspc,KIND=RealKind)

  !-----------------------------------------------------------------------
  ! --- If more than one argument is passed set new tolerance and ROS methode
  CALL getarg ( 2 , neuTolR )             
  IF ( neuTolR /= '' )  READ( neuTolR , * ) RtolROW

  CALL getarg ( 3 , neuTolA )             
  IF ( neuTolA /= '' )  READ( neuTolA , * ) AtolGas

  CALL getarg ( 4 , neuTolA )             
  IF ( neuTolA /= '' )  READ( neuTolA , * ) AtolAqua
 
  CALL getarg ( 5 , neuROW )
  IF ( neuROW /= '' )   THEN  
    READ( neuRow , * ) ODEsolver
    ODEsolver = ADJUSTL('METHODS/'//ODEsolver)
  END IF
  
  !-----------------------------------------------------------------------
  !--- get all sulphuric species
  nSchwefel = 0
  DO i = 1 , nspc
    IF ( INDEX(y_name(i),'S') > 0 ) nSchwefel = nSchwefel + 1
  END DO

  ALLOCATE(iDiag_Schwefel(nSchwefel))
  nSchwefel = 0
  DO i=1,nspc
    IF ( INDEX(y_name(i),'S') > 0 ) THEN
      nSchwefel = nSchwefel + 1
      iDiag_Schwefel(nSchwefel) = PositionSpeciesAll( y_name(i) ) 
    END IF
  END DO
  
  !-----------------------------------------------------------------------
  !--- print input parameter, method, tols, etc.
  CALL Print_Run_Param()
  !print*, 'NTypes%AquaPhoto=',NTypes%AquaPhoto
  !print*, 'NTypes%Spec1MCM=',NTypes%Spec1MCM
  !print*, 'NTypes%Spec1MCM=',NTypes%GasConst
  !stop
  !-----------------------------------------------------------------------
  CALL Integrate (  InitValAct(1:nspc)   &  ! initial concentrations activ species
  &               , Tspan                &  ! integration invervall
  &               , Atol                 &  ! abs. tolerance gas species
  &               , RtolROW              &  ! rel. tolerance Rosenbrock method
  &               , ODEsolver     )         ! methode for solving the ode
  !---------------------------------------------------------------
  ! --- stop timer and print output statistics
  Timer_Finish = MPI_WTIME()
  Timer_Finish = Timer_Finish - Timer_Start + Time_Read
  !
  ! Print statistics
  CALL Output_Statistics( Time_Read       , TimeSymbolic  , TimeFac       &
  &                     , TimeSolve       , TimeRates     , TimeJac       &
  &                     , TimeIntegrationE, Timer_Finish  , TimeRateSend  &
  &                     , TimeNetCDF      , TimeErrCalc                   )
  WRITE(*,*) ''
  !---------------------------------------------------------------
  ! --- Close MPI 
  CALL FinishMPI()
END PROGRAM MAIN_CHEMIE
