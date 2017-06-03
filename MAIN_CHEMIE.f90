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
  USE NetCDF_Mod
  IMPLICIT NONE
  !
  CHARACTER(80)   :: Filename0 = ''        ! *.run file
  !
  REAL(dp) :: Atol(2)
  INTEGER  :: i,j,nSchwefel
  INTEGER  :: io_err,  STAT

  ! NetCDF stuff
  REAL(dp) :: StartTimer
  REAL(dp) :: actLWC, zen, wetRad
  INTEGER  :: errind(1,1)
  REAL(dp), ALLOCATABLE :: ErrVals(:), Y(:)
  ! reaction rate array + part. derv. rate over temperatur vector
  REAL(dp), ALLOCATABLE :: Rate(:), DRatedT(:)  
  !
  ! use other ROS methode or tolerance
  CHARACTER(80) :: neuTolR  = ''
  CHARACTER(80) :: neuTolA  = ''
  CHARACTER(80) :: neuROW   = ''

  ! temp variables for molecular weights
  CHARACTER(16) :: tmpchar0 = '-'
  REAL(dp)  :: tmpMW0
  INTEGER   :: tmpPos, tmpCnt

  ! convertion from mole to mass to conc
  REAL(dp), ALLOCATABLE   :: MoleFrac(:), MassFrac(:), MoleConc(:)
  REAL(dp)                :: Press_in_dyncm2

  ! matricies for symbolic phase
  TYPE(CSR_Matrix_T)     :: Id , tmpJacCC  ! compressed row
  TYPE(SpRowColD_T)      :: temp_LU_Dec    ! sparse-LU matrix format
  ! permutation vector/ pivot order for LU-decomp
  INTEGER, ALLOCATABLE :: InvPermu(:)
  INTEGER, ALLOCATABLE :: PivOrder(:)

  ! format string
  CHARACTER(14) :: fmt0=''
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


    !CALL PM_int(AtomicMatrix)
    !stop 'MAIN after PM'

    !--- richtigen index holen, da TB unsortiert eingelesen wurde
    CALL GatherTBindex

    IF (Vectorized) CALL GatherReactionTypeIndex
   
    !--- Read initial values
    ALLOCATE( InitValAct(ntGas) , y_e(ntGas) )
    ALLOCATE( InitValKat(ntKat) )

    IF ( MWeights /= '' ) THEN
      CALL Read_MolecularWeights(MW,MWeights,MWUnit,nspc)

      ALLOCATE( rMW(nspc) )  
      rMW = ONE / MW
      
      ALLOCATE( MoleFrac(ntGas)   , MassFrac(ntGas) )
      MoleFrac    = ZERO           ! mole fraction 
      MassFrac    = ZERO           ! mass fraction 

      CALL Read_GASini    ( InitFile , MoleFrac , InitValKat )

      !Press = Pressure0               ! initial pressure in [Pa]
      Press_in_dyncm2 = Pressure0 * Pa_to_dyncm2

      MassFrac = MoleFr_To_MassFr( MoleFrac ) 
      
      MoleConc = MoleFr_To_MoleConc( MoleFrac,               &
                                   & Press = Press_in_dyncm2,&
                                   & Temp  = Temperature0    )
    ELSE
      IF (MPI_ID == 0 ) THEN
        WRITE(*,*) ''
        WRITE(*,*) '  No molecular weights are given.  '
        WRITE(*,*) '       ---> Initial Values in [mole/cm3] '
        WRITE(*,*) ''
      END IF
      CALL Read_GASini    ( InitFile , MoleConc , InitValKat )
    END IF
    CALL Read_EMISS     ( InitFile , y_e )
    
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
   
    !-----------------------------------------------------------------------
    ! --- Build the reaction system
    !-----------------------------------------------------------------------
    CALL AllListsToArray( ReactionSystem  , ListRGas       &
    &                    , ListRHenry     , ListRAqua      &
    &                    , ListRDiss                       )
    !
    !-----------------------------------------------------------------------
    ! --- print reactions and build A, B and (B-A) structure
    !-----------------------------------------------------------------------
    CALL PrintReactions( ReactionSystem , 89 )
    CALL PrintFinalReactions( 89 )

    !-----------------------------------------------------------------------
    ! --- Input of initial data and thermodynamic properties
    CALL InputChemicalData( InitFile , DataFile , MetFile )
    !-----------------------------------------------------------------------
  END IF
  CLOSE(89)

  IF (MPI_ID==0) WRITE(*,*) '  Read ini-file ................ done'
  IF (MPI_ID==0) WRITE(*,*) '  Print chem-file .............. done'

  !-----------------------------------------------------------------------
  ! --- Read species for diagnose (print species concs to NetCDF file)
  CALL Read_Diag( OutNetcdfSpc , OutNetcdfPhase , InitFile )
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! --- Timers
  Time_Read   = MPI_WTIME() - Time_Read  ! Stop Input timer
  Timer_Start = MPI_WTIME()              ! Start timer for integration
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! --- Split the rate array in mpi_np parts 
  CALL BuildPartitions( MyParties , neq , MPI_np )
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  ! --- this is for the new mass action product routine 
  CALL GatherSpeciesOrder(A)
  !-----------------------------------------------------------------------

  
  !-----------------------------------------------------------------------
  ! --- Dimension initialisation for the unknowns and matrices
  !-----------------------------------------------------------------------

  nsr = nspc + neq

  IF ( TempEq ) THEN
    nDIM    = nspc + 1
    nDIMcl  = nspc + 1
    nDIMex  = nsr  + 1
    Atol    = (/ AtolGas , AtolTemp /)

    !--- malloc gibbs energy, derivates
    ALLOCATE( GFE(nspc)   , DGFEdT(nspc)   )
    ALLOCATE( DelGFE(neq) , DDelGFEdT(neq) )
    GFE    = ZERO;  DGFEdT    = ZERO
    DelGFE = ZERO;  DDelGFEdT = ZERO
  ELSE
    nDIM    = nspc
    nDIMcl  = nspc
    nDIMex  = nsr
    Atol    = (/ AtolGas , AtolAqua /)
  END IF

  rNspc = ONE/REAL(nspc,KIND=dp)  ! rNspc for error calculation

  !---------------------------------------------------------------------------
  ! --- If more than one argument is passed set new tolerance and ROS methode
  !---------------------------------------------------------------------------
  CALL getarg ( 2 , neuTolR ); IF ( neuTolR /= '' )  READ( neuTolR , * ) RtolROW
  CALL getarg ( 3 , neuTolA ); IF ( neuTolA /= '' )  READ( neuTolA , * ) AtolGas
  CALL getarg ( 4 , neuTolA ); IF ( neuTolA /= '' )  READ( neuTolA , * ) AtolAqua
  CALL getarg ( 5 , neuROW );  IF ( neuROW  /= '' )  READ( neuRow  , * ) ODEsolver
  

  !-----------------------------------------------------------------------
  ! --- Get all sulphuric species
  !-----------------------------------------------------------------------
  ! count species
  nSchwefel = 0
  DO i=1,nspc
    IF (INDEX(y_name(i),'S')>0) nSchwefel = nSchwefel + 1
  END DO

  ! save species index to iDiag_Schwefel
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
  !
  !--- Print Init parameters, concs, temp, pressure,....
  IF (MPI_ID==0) THEN
    IF ( Ladebalken==1 ) Bar = .TRUE.
    WRITE(*,*) '  Initial values:   '
    WRITE(*,*)
    fmt0 = '  [molec/cm3] '
    IF ( TempEq ) fmt0 = '  [mol/cm3]'
    IF (ntGas>0)  WRITE(*,798) SUM(InitValAct(1:ntGas)) , fmt0
    IF (ntAqua>0) WRITE(*,799) SUM(InitValAct(ntGas+1:nspc)) , fmt0
    WRITE(*,800) SUM(Y_e) , fmt0

    IF ( TempEq ) THEN
      WRITE(*,801) Temperature0
      WRITE(*,802) Pressure0
      WRITE(*,803) rho
    END IF
    WRITE(*,*)
  ELSE
    ! if more than 1 process, output data only on master process
    NetCdfPrint = .FALSE.
  END IF

  !-----------------------------------------------------------------------
  ! ---  check the consistency of input parameters
  !-----------------------------------------------------------------------
  CALL CheckInputParameters( Tspan , Atol , RtolROW )
  ALLOCATE( ErrVals(nspc) , Y(nspc)  )
  ALLOCATE( Rate(neq) , DRatedT(neq) )

  !-----------------------------------------------------------------------
  ! --- Initialize NetCDF output file
  !-----------------------------------------------------------------------

  IF ( NetCdfPrint ) THEN 
    StartTimer  = MPI_WTIME()

    itime_NetCDF = 0
    iStpNetCDF   = 1
  
    ! aqua phase in [mol/m3] and [mol/l]
    OutNetcdfANZ = COUNT(OutNetcdfPhase=='g') + 2*COUNT(OutNetcdfPhase=='a')
    OutNetcdfDIM = OutNetcdfANZ
    IF ( TempEq ) OutNetcdfDIM = OutNetcdfANZ + 1
    ALLOCATE(yNcdf(OutNetcdfDIM))          ! output array , +1 for Temperature
    yNcdf = ZERO
        
    !--  Netcdf Output File
    IF ( ErrorLog == 1 ) ErrVals = ZERO
    IF ( ntAqua > 0 ) THEN
      actLWC = pseudoLWC(Tspan(1))
      wetRad = (Pi34*actLWC/SPEK(1)%Number)**(rTHREE)*1.0d-1
    ELSE
      actLWC = ZERO
      wetRad = ZERO
    END IF

    errind(1,1) = 0
    CALL InitNetcdf( Tspan(1) , Tspan(2) )

    WRITE(*,*) '  Init NetCDF................... done'
    !--  print initial values to NetCDF file
    CALL SetOutputNCDF( InitValAct, yNcdf, Tspan(1), actLWC )
    CALL StepNetCDF   ( Tspan(1) , yNcdf(:) , itime_NetCDF ,  &
                      & (/ actLWC , ZERO ,                        &
                      &    SUM(InitValAct(1:ntGas)),              &
                      &    SUM(InitValAct(ntGas+1:ntGas+ntAqua)), &
                      &    wetRad , ZERO   /),                    &
                      &  errind , SUM(InitValAct(iDiag_Schwefel)),&
                      &  ZERO                            )
    TimeNetCDF = MPI_WTIME() - StartTimer
  END IF

  !-----------------------------------------------------------------------
  ! --- Symbolic phase 
  !    - generating sparse matrecies of stoechiometic coefs (beta-alpha)^T 
  !    - generating sparse matrecies of Jacobian matrix (if nessesarry)
  !    - generating sparse matrecies of LU decomposition of the Jacobian
  !-----------------------------------------------------------------------
  StartTimer = MPI_WTIME()            ! start timer for symb phase
  
  CALL SymbolicAdd( BA , B , A )      ! symbolic addition:    BA = B + A
  CALL SparseAdd  ( BA , B , A, '-' )  ! numeric subtraction:  BA = B - A
  CALL TransposeSparse( BAT , BA )    ! transpose BA:        BAT = Transpose(BA) 
 
  ! if either Rosenbrock method or backward Euler was choosen
  IF ( MAXVAL(INDEX(TRIM(ODEsolver(1:7)),['METHODS','bwEuler']))>0) THEN
    !---- Get Rosenbrock Parameters
    CALL SetRosenbrockMethod( RCo , ODEsolver )  

    ! we need to calculate the Jacobian for both versions 'cl' and 'ex' to
    ! calculate an initial stepsize based on 2nd derivative (copy of MATLABs ode23s)
    ! symbolic mult for Jacobian Jac = [ BAT * Ones_nR * A * Ones_nS ], 
    ! add diagonal entries and set them to zero
    CALL SymbolicMult( BAT , A , tmpJacCC )  
    CALL SparseID( Id , nspc )                    
    CALL SymbolicAdd( Jac_CC , Id , tmpJacCC )
    CALL Free_Matrix_CSR( Id )

    ! Set symbolic structure of iteration matrix for Row-Method
    ! also assign constant matrix parts like (beta-alpha)^T, alpha
    IF ( CLASSIC ) THEN
      CALL BuildSymbolicClassicMatrix(  Miter , Jac_CC  , RCo%ga )
      IF (DebugPrint) CM_1 = Copy_CSR(Jac_CC)
    ELSE
      CALL BuildSymbolicExtendedMatrix( Miter , A , BAT , RCo%ga ) 
    END IF

    ! Choose ordering/factorisation strategie and do symb LU fact
    IF ( useMUMPS ) THEN 
      ! Use MUMPS to factorise and solve     
      ! Convert compressed row format to row index format for MUMPS
      CALL CompRowToIndRow( Miter , MiterFact )
      CALL InitMumps( MiterFact ) 

      IF ( MatrixPrint ) THEN
        CALL CSRToSpRowColD   ( temp_LU_Dec , Miter ) 
        CALL PermuToInvPer    ( InvPermu   , Mumps_Par%SYM_PERM )
        CALL SymbLU_SpRowColD ( temp_LU_Dec , InvPermu )        
        CALL RowColDToCSR     ( LU_Miter   , temp_LU_Dec , nspc , neq ) 
      END IF

    ELSE IF ( useSparseLU ) THEN
      ! Permutation given by Markowitz Ordering strategie
      CALL CSRToSpRowColD( temp_LU_Dec , Miter) 

      IF (OrderingStrategie==8) THEN
        CALL SymbLU_SpRowColD_M ( temp_LU_Dec )        
      ELSE 
        ALLOCATE(PivOrder(temp_LU_Dec%n))
        PivOrder = -90
        
        ! Pivots bis neq in reihenfolge 1,2,...,neq
        IF ( EXTENDED ) THEN
          PivOrder(     1 : neq      ) = [(i , i = 1     , neq  )]
          PivOrder( neq+1 : neq+nDIM ) = [(i , i = neq+1 , neq+nDIM )]
        ELSE
          PivOrder(     1 : nDIM     ) = [(i , i = 1     , nDim )]
        END IF
        CALL SymbLU_SpRowColD ( temp_LU_Dec , PivOrder)
      END IF

      ! converting back to csr format
      CALL RowColDToCSR( LU_Miter , temp_LU_Dec , nspc , neq )

      ! Get the permutation vector LU_Perm and map values of Miter
      ! to the permuted LU matrix
      CALL Get_LU_Permutaion( LU_Perm , LU_Miter , Miter , nspc , neq )
      
      ! For the extended case one can save the values of alpha and 
      ! (beta-alpha)^T and just copy them for each iteration in ROS
      IF ( EXTENDED ) ALLOCATE(LUValsFix(LU_Miter%nnz)); LUvalsFix = LU_Miter%Val
    END IF

    IF (MPI_ID==0) THEN
      IF (MatrixPrint) THEN
        CALL WriteSparseMatrix(tmpJacCC,'MATRICES/JAC_'//TRIM(BSP)//'_'//solveLA     , neq, nspc)
        CALL WriteSparseMatrix(BA,      'MATRICES/BA_'//TRIM(BSP)//'_'//solveLA      , neq, nspc)
        CALL WriteSparseMatrix(Miter,   'MATRICES/Miter_'//TRIM(BSP)//'_'//solveLA   , neq, nspc)
        CALL WriteSparseMatrix(LU_Miter,'MATRICES/LU_Miter_'//TRIM(BSP)//'_'//solveLA, neq, nspc)
        STOP 'after writesparsematrix'
      END IF
      WRITE(*,*) '  Symbolic-phase................ done'
      WRITE(*,*) ' '
    END IF

    CALL Free_Matrix_CSR( tmpJacCC )
    CALL Free_SpRowColD( temp_LU_Dec )

    TimeSymbolic = MPI_WTIME() - StartTimer   ! stop timer for symbolic matrix calculations

    ! ---- Calculate first reaction rates
    CALL ReactionRatesAndDerivative( Tspan(1) , InitValAct , Rate , DRatedT )
    Y = MAX( ABS(InitValAct) , eps ) * SIGN( ONE , InitValAct )    ! |y| >= eps
    
    ! ---- Calculate values of Jacobian
    StartTimer  = MPI_WTIME()
    CALL Jacobian_CC(Jac_CC , BAT , A , Rate , Y )
    TimeJac     = TimeJac + MPI_WTIME() - StartTimer
    Output%npds = Output%npds + 1
  END IF

  !-----------------------------------------------------------------------
  !--- Routines for pathway analysis
  !-----------------------------------------------------------------------
  !IF (MPI_ID==0) THEN
  !  WRITE(*,*)
  !  DO i = 1 , SIZE(OutNetcdfSpc)
  !    CALL SearchReactions(y_name(OutNetcdfSpc(i)))
  !  END DO
  !  WRITE(*,*)
  !END IF
  !stop 'MAIN'

  
  !-----------------------------------------------------------------------
  ! --- Start the integration routine 
  !-----------------------------------------------------------------------
  CALL Integrate (  InitValAct(1:nspc)   &  ! initial concentrations activ species
  &               , Rate                 &  ! reaction rates at time=t0
  &               , Tspan                &  ! integration invervall
  &               , Atol                 &  ! abs. tolerance gas species
  &               , RtolROW              &  ! rel. tolerance Rosenbrock method
  &               , ODEsolver     )         ! methode for solving the ode
  !---------------------------------------------------------------
  ! --- stop timer and print output statistics
  Timer_Finish = MPI_WTIME() - Timer_Start + Time_Read
  !
  ! Print statistics
  CALL Output_Statistics( Time_Read       , TimeSymbolic  , TimeFac       &
  &                     , TimeSolve       , TimeRates     , TimeJac       &
  &                     , TimeIntegrationE, Timer_Finish  , TimeRateSend  &
  &                     , TimeNetCDF      , TimeErrCalc   , TimeRhsCalc   )
  WRITE(*,*) ''
  !---------------------------------------------------------------
  ! --- Close MPI 
  CALL FinishMPI()

  !================================================================
  !==  FORMAT Statements
  !================================================================
  !
  798  FORMAT('      Sum Initval (gaseous)      = ', E16.10, A)
  799  FORMAT('      Sum Initval (aqueous)      = ', E16.10, A)
  800  FORMAT('      Sum Emissions (gaseous)    = ', E16.10, A)
  801  FORMAT('      Initial Temperature        = ', E16.10,'  [K]') 
  802  FORMAT('      Initial Pressure           = ', E16.10,'  [Pa]')
  803  FORMAT('      Reactor denstiy            = ', E16.10,'  [kg/cm3]')
END PROGRAM MAIN_CHEMIE
