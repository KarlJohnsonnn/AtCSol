!=======================================
!=======================================
! vorlaeufiges Hauptprogramm fuer Tests
! box model chemical kinetics simulation
!=======================================
!=======================================
!
!
PROGRAM chemie
  !
  USE Kind_Mod
  USE Sparse_Mod
  USE Chemsys_Mod
  USE Integration_Mod
  USE Rosenbrock_Mod
  USE mo_control
  USE mo_reac
  USE mo_MPI
  USE mo_IO
  USE mo_ckinput
  USE NetCDF_Mod
  USE Cycles_Mod
  USE mo_reduction
  USE fparser
  USE issa
  IMPLICIT NONE
  !
  CHARACTER(80)   :: Filename0 = ''        ! *.run file
  CHARACTER(80)   :: ChemFile  = ''        ! *.chem file (output)
  !
  REAL(dp) :: Atol(2)
  INTEGER  :: i,j,jj,nSchwefel
  INTEGER  :: io_err,  STAT

  ! NetCDF stuff
  REAL(dp) :: StartTimer
  REAL(dp) :: LWC, zen
  INTEGER  :: errind(1,1)
  REAL(dp), ALLOCATABLE :: ErrVals(:), Y(:)
  ! reaction rate array + part. derv. rate over temperatur vector
  REAL(dp), ALLOCATABLE :: Rate(:), DRatedT(:)  , wetRad(:)
  !
  ! use other ROS methode or tolerance
  CHARACTER(80) :: neuTolR  = ''
  CHARACTER(80) :: neuTolA  = ''
  CHARACTER(80) :: neuROW   = ''
  CHARACTER(2)  :: newSolveLA  = ''

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

  ! new testing stuff 
  TYPE(CSR_Matrix_T) :: E_HG, P_HG, S_HG, R_HG, S_HG_transp
  CHARACTER(1) :: inpt=''
  
  ! timer
  REAL(dp) :: t_1,t_2

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
    WRITE(*,*) 'Input RUNFilename: '; READ(*,*)   FileName0
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
  Tspan = (/ tAnf , tEnd /)
 
  !----------------------------------------------------------------
  ! --- set cloud intervall
  LWCBounds(1) = tAnf * HOUR
  LWCBounds(2) = LWCBounds(1) + 1.00_dp  * HOUR
  LWCBounds(3) = LWCBounds(2) + 0.25_dp  * HOUR
  LWCBounds(4) = LWCBounds(3) + 9.50_dp  * HOUR
  LWCBounds(5) = LWCBounds(4) + 0.25_dp  * HOUR
  LWCBounds(6) = LWCBounds(5) + 1.00_dp  * HOUR

  !----------------------------------------------------------------
  !  --- read the .sys data, save coefs in sparse matrix
  Time_Read = MPI_WTIME()

  ChemFile = ADJUSTL(TRIM(SysFile(:INDEX(SysFile,'.sys')-1)))
  OPEN ( UNIT=89 , FILE=ADJUSTL(TRIM(ChemFile))//'.chem' , STATUS='UNKNOWN' )

  IF ( MPI_ID == 0 ) WRITE(*,'(A33)',ADVANCE='NO') '   Reading sys-file .............'

  IF ( Teq.AND.ChemKin ) THEN
    
    IF ( MPI_ID == 0 ) WRITE(*,*) 'done   ---->  Solve Gas Energy Equation '

    CALL Read_Elements    ( SysFile    , 969 )
    CALL Read_Species     ( SysFile    , 969 )
    CALL Read_Reaction    ( SysFile    , 969 )

    CALL PrintHeadSpecies   ( ChemFile    , 89 ) 
    CALL PrintSpecies       ( ListGas2    , 89 )
    CALL PrintHeadReactions ( 89 )
    CALL PrintReactions     ( ReactionSystem , 89 , .TRUE. )
    CALL PrintFinalReactions( 89 )

    CALL GetSpeciesNames( ChemFile , y_name )
    CALL Read_ThermoData( SwitchTemp , DataFile , 696 , nspc )


    !--- richtigen index holen, da TB unsortiert eingelesen wurde
    CALL Setup_ThirdBodyIndex
    CALL Setup_ReactionIndex
   
    !--- Read initial values
    ALLOCATE( InitValAct(ns_GAS) , y_e(ns_GAS) , InitValKat(ns_KAT) )

    IF ( MWeights /= '' ) THEN
      CALL Read_MolecularWeights(MW,MWeights,MWUnit,nspc)

      ALLOCATE( rMW(nspc) )  
      rMW = ONE / MW
      
      ALLOCATE( MoleFrac(ns_GAS) , MassFrac(ns_GAS) )
      MoleFrac = ZERO         ! mole fraction 
      MassFrac = ZERO         ! mass fraction 

      CALL Read_GASini( InitFile , MoleFrac , InitValKat )

      !Press = Pressure0               ! initial pressure in [Pa]
      Press_in_dyncm2 = Pressure0 * Pa_to_dyncm2

      MassFrac = MoleFr_To_MassFr  ( MoleFrac ) 
      MoleConc = MoleFr_To_MoleConc( MoleFrac                &
      &                            , Press = Press_in_dyncm2 &
      &                            , Temp  = Temperature0    )
    ELSE
      IF (MPI_ID == 0 ) THEN
        WRITE(*,*)
        WRITE(*,*) '  No molecular weights are given.  '
        WRITE(*,*) '  Make sure the initial values are given in [mole/cm3] !'
        WRITE(*,*)
      END IF
      CALL Read_GASini( InitFile , MoleConc , InitValKat )
    END IF
    CALL Read_EMISS( InitFile , y_e )
    
    ! Initialising reactor density
    rho  = Density( MoleConc )
    rRho = mega/rho       ! in [cm3/g]
    InitValAct = MoleConc

  ELSE

    IF ( MPI_ID==0 ) WRITE(*,*) 'done  ---->  Fix Temperature'
    CALL ReadSystem( SysFile )
  
    !----------------------------------------------------------------
    ! ---  build the coeficient matrices and write .chem
    CALL PrintHeadSpecies ( ChemFile , 89 )

    IF ( ns_GAS    > 0 ) CALL PrintSpecies( ListGas2     , 89 )
    IF ( ns_AQUA   > 0 ) CALL PrintSpecies( ListAqua2    , 89 )
    IF ( ns_SOLID  > 0 ) CALL PrintSpecies( ListSolid2   , 89 ) 
    IF ( ns_PARTIC > 0 ) CALL PrintSpecies( ListPartic2  , 89 )
    IF ( ns_KAT    > 0 ) CALL PrintSpecies( ListNonReac2 , 89 )

    CALL PrintHeadReactions( 89 )
   
    !-----------------------------------------------------------------------
    ! --- Build the reaction system
    !-----------------------------------------------------------------------
    CALL AllListsToArray( ReactionSystem            &  
    &                   , ListRGas    , ListRHenry  &
    &                   , ListRAqua   , ListRDiss   &
    &                   , ListRSolid  , ListRPartic &
    &                   , ListRMicro  )
    
    !-----------------------------------------------------------------------
    ! --- print reactions and build A, B and (B-A) structure
    !-----------------------------------------------------------------------
    CALL PrintReactions( ReactionSystem , 89 )
    CALL PrintFinalReactions( 89 )

    !-----------------------------------------------------------------------
    ! --- initialize fpraser for reactions with special rate formula
    !-----------------------------------------------------------------------
    IF ( nr_special > 0 ) THEN
     
      ! Initialize function parser for n special functions
      CALL initf( nr_special ) 
      
      ! Parse and bytecompile ith function string 
      DO i = 1,nr_special
        CALL parsef ( i, ReactionSystem(iR%iSPECIAL(i))%Special%Formula    &
        &              , ReactionSystem(iR%iSPECIAL(i))%Special%cVariables )
      END DO
    END IF

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
  ! --- this is for the new mass action product routine 
  CALL Setup_SpeciesOrder(A)
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  ! --- Dimension initialisation for the unknowns and matrices
  !-----------------------------------------------------------------------
  !
  nsr = nspc + neq

  IF ( Teq ) THEN
    nDIM    = nspc + 1
    nDIMcl  = nspc + 1
    nDIMex  = nsr  + 1
    Atol    = [ AtolGas , AtolTemp ]

    !--- malloc gibbs energy, derivates
    ALLOCATE( GFE(nspc)   , DGFEdT(nspc)   )
    ALLOCATE( DelGFE(neq) , DDelGFEdT(neq) )
    GFE    = ZERO;  DGFEdT    = ZERO
    DelGFE = ZERO;  DDelGFEdT = ZERO
  ELSE
    nDIM    = nspc
    nDIMcl  = nspc
    nDIMex  = nsr
    Atol    = [ AtolGas , AtolAqua ]
  END IF

  rNspc = ONE/REAL(nspc,KIND=dp)  ! rNspc for error calculation

  ! true if mechanism contains photolytic reactions
  PHOTO = (nr_G_photo+nr_A_photo) > 0

  !---------------------------------------------------------------------------
  ! --- If more than one argument is passed set new tolerance and ROS methode
  !---------------------------------------------------------------------------

  ! if you want to change the calc of the linear systems without opening the run-file
  CALL getarg ( 2 , newSolveLA )
  IF (newSolveLA /= '') THEN
    SELECT CASE (newSolveLA)
      CASE ('cl'); CLASSIC = .TRUE.;  EXTENDED = .FALSE.; SolveLA = newSolveLA
      CASE ('ex'); CLASSIC = .FALSE.; EXTENDED = .TRUE.;  SolveLA = newSolveLA
      CASE DEFAULT;  WRITE(*,*) '  Linear Algebra either "cl" or "ex" !'
    END SELECT
  END IF

  !-----------------------------------------------------------------------
  ! --- Get all sulphuric species
  !-----------------------------------------------------------------------
  ! count species
  ALLOCATE(iDiag_Schwefel(0));  nSchwefel = 0
  DO i=1,nspc
    IF ( INDEX(y_name(i),'S') > 0 ) THEN
      nSchwefel = nSchwefel + 1
      iDiag_Schwefel = [iDiag_Schwefel , PositionSpeciesAll(y_name(i))]
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
    IF ( Teq ) fmt0 = '  [mol/cm3]'
    IF ( ns_GAS  > 0 ) WRITE(*,798) SUM(InitValAct(1:ns_GAS)) , fmt0
    IF ( ns_AQUA > 0 ) WRITE(*,799) SUM(InitValAct(ns_GAS+1:nspc)) , fmt0
    WRITE(*,800) SUM(Y_e) , fmt0

    IF ( Teq ) THEN
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
  ALLOCATE( wetRad(ntFrac) ) 

  !-----------------------------------------------------------------------
  ! --- Initialize NetCDF output file
  !-----------------------------------------------------------------------

  IF ( NetCdfPrint ) THEN 
    StartTimer  = MPI_WTIME()

    itime_NetCDF = 0
    iStpNetCDF   = 1
  
    ! aqua phase in [mol/m3] and [mol/l]
    OutNetcdfANZ = COUNT(OutNetcdfPhase=='g') + 2*ntFrac*COUNT(OutNetcdfPhase=='a')
    OutNetcdfDIM = OutNetcdfANZ
    IF ( Teq ) OutNetcdfDIM = OutNetcdfANZ + 1
    ALLOCATE(yNcdf(OutNetcdfDIM))          ! output array , +1 for Temperature
    yNcdf = ZERO
        
    !--  Netcdf Output File
    IF ( ErrorLog == 1 ) ErrVals = ZERO
    IF ( ns_AQUA > 0 ) THEN
      LWC = pseudoLWC(Tspan(1))
      wetRad(:) = (Pi34*LWC/Frac%Number(:))**(rTHREE)*0.1_dp
    ELSE
      LWC = ZERO
      wetRad(:) = ZERO
    END IF

    errind(1,1) = 0
    CALL InitNetcdf( Tspan(1) , Tspan(2) )

    WRITE(*,*) '  Init NetCDF................... done'
    !--  print initial values to NetCDF file
    CALL SetOutputNCDF( InitValAct, yNcdf, Tspan(1), LWC )
    CALL StepNetCDF   ( Tspan(1) , yNcdf(:) , itime_NetCDF ,      &
                      & (/ LWC , ZERO ,                           &
                      &    SUM(InitValAct(1:ns_GAS)),              &
                      &    SUM(InitValAct(ns_GAS+1:ns_GAS+ns_AQUA)), &
                      &    wetRad , ZERO   /),                    &
                      &  errind , SUM(InitValAct(iDiag_Schwefel)),&
                      &  ZERO                                     )
    TimeNetCDF = MPI_WTIME() - StartTimer
  END IF

  !-----------------------------------------------------------------------
  ! --- Symbolic phase 
  !    - generating sparse matrecies of stoechiometic coefs (beta-alpha)^T 
  !    - generating sparse matrecies of Jacobian matrix (if nessesarry)
  !    - generating sparse matrecies of LU decomposition of the Jacobian
  !-----------------------------------------------------------------------
  WRITE(*,'(A33)',ADVANCE='NO') '   Symbolic-phase................'
  StartTimer = MPI_WTIME()            ! start timer for symb phase
  
  CALL SymbolicAdd( BA , B , A )      ! symbolic addition:    BA = B + A
  CALL SparseAdd  ( BA , B , A, '-' )  ! numeric subtraction:  BA = B - A
  CALL TransposeSparse( BAT , BA )    ! transpose BA:        BAT = Transpose(BA) 

  IF ( MatrixPrint ) THEN
    ! write educt and product matrix to file
    P_HG = Copy_CSR(B)
    CALL TransposeSparse( E_HG , Copy_CSR(A) )
    P_HG%Val = ONE
    E_HG%Val = ONE

    CALL SymbolicMult( E_HG , P_HG , S_HG )
    CALL SymbolicMult( P_HG , E_HG , R_HG )
    S_HG%Val = ONE;   R_HG%Val = ONE

    CALL TransposeSparse( S_HG_transp , S_HG )

    IF( Targets/='' ) THEN
      CALL Read_Target_Spc(Target_Index,Target_Names,Targets)
      CALL Find_Elem_Circuits(S_HG,Target_Index)

      ! issa routines
      CALL ISSA_nue(BAT,Cyclic_Set,ReactionSystem)


      !WRITE(*,*) '  Continue? [y/n]'
      !READ(*,*) inpt
      !IF (inpt/='y') STOP
      STOP ' main '
    END IF
    
    !STOP ' kreise gefunden '

    !CALL WriteSparseMatrix(P_HG,'MATRICES/P_HG_'//BSP,neq,nspc)
    !CALL WriteSparseMatrix(E_HG,'MATRICES/E_HG_'//BSP,neq,nspc)
    !CALL WriteSparseMatrix(S_HG,'MATRICES/S_HG_'//BSP,neq,nspc)
    !CALL WriteSparseMatrix(R_HG,'MATRICES/R_HG_'//BSP,neq,nspc)
  END IF

  ! if either Rosenbrock method or backward Euler was choosen
  IF ( MAXVAL(INDEX(TRIM(ODEsolver(1:7)),['METHODS','bwEuler']))>0) THEN
    !---- Get Rosenbrock Parameters
    CALL SetRosenbrockMethod( ROS , ODEsolver )  

    ! we need to calculate the Jacobian for both versions 'cl' and 'ex' to
    ! calculate an initial stepsize based on 2nd derivative (copy of MATLABs ode23s)
    ! symbolic mult for Jacobian Jac = [ BAT * Ones_nR * A * Ones_nS ], 
    ! add diagonal entries and set them to zero
    CALL SymbolicMult( BAT , A , tmpJacCC )  
    Id = SparseID( nspc )                    
    CALL SymbolicAdd( Jac_CC , Id , tmpJacCC )
    CALL Free_Matrix_CSR( Id )

    ! gephi testing
    !CALL CSR_to_GephiGraph(tmpJacCC,y_name,'Test1')
    !stop ' after gephi plot'

    ! Set symbolic structure of iteration matrix for Row-Method
    ! also assign constant matrix parts like (beta-alpha)^T, alpha
    IF ( CLASSIC ) THEN
      CALL BuildSymbolicClassicMatrix(  Miter , Jac_CC  , ROS%ga )
      IF (DebugPrint) CM_1 = Copy_CSR(Jac_CC)
    ELSE !IF ( EXTENDED ) THEN
      CALL BuildSymbolicExtendedMatrix( Miter , A , BAT , ROS%ga ) 
    END IF

    ! Choose ordering/factorisation strategie and do symb LU fact
    IF ( useMUMPS ) THEN 
      ! Use MUMPS to factorise and solve     
      ! Convert compressed row format to row index format for MUMPS
      MiterFact = CSR_to_SpRowIndColInd( Miter )
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

      IF ( OrderingStrategie == 8 ) THEN
        CALL SymbLU_SpRowColD_M ( temp_LU_Dec )        
      ELSE 
        ALLOCATE(PivOrder(temp_LU_Dec%n))
        PivOrder = -90
        
        ! Pivots bis neq in reihenfolge 1,2,...,neq
        IF ( CLASSIC ) THEN
          PivOrder(     1 : nDIM     ) = [(i , i = 1     , nDim )]
        ELSE !IF ( EXTENDED ) THEN
          PivOrder(     1 : neq      ) = [(i , i = 1     , neq  )]
          PivOrder( neq+1 : neq+nDIM ) = [(i , i = neq+1 , neq+nDIM )]
        END IF
        CALL SymbLU_SpRowColD(temp_LU_Dec , PivOrder)
      END IF

      ! converting back to csr format
      CALL RowColDToCSR( LU_Miter , temp_LU_Dec , nspc , neq )

      ! Get the permutation vector LU_Perm and map values of Miter
      ! to the permuted LU matrix
      CALL Get_LU_Permutaion( LU_Perm , LU_Miter , Miter , nspc , neq )
      
      ! For the extended case one can save the values of alpha and 
      ! (beta-alpha)^T and just copy them for each iteration in ROS
      IF ( EXTENDED ) LUvalsFix = LU_Miter%Val
    END IF

    IF (MPI_ID==0) THEN
      !IF (MatrixPrint) THEN
      !  CALL WriteSparseMatrix(tmpJacCC,'MATRICES/JAC_'//TRIM(BSP), neq, nspc)
      !  CALL WriteSparseMatrix(BA,      'MATRICES/BA_'//TRIM(BSP), neq, nspc)
      !  CALL WriteSparseMatrix(Miter,   'MATRICES/Miter_'//TRIM(BSP), neq, nspc)
      !  CALL WriteSparseMatrix(LU_Miter,'MATRICES/LU_Miter_'//TRIM(BSP), neq, nspc)
      !  !STOP 'after writesparsematrix'
      !END IF
      WRITE(*,*) 'done'
      WRITE(*,*) ' '
      WRITE(*,*) '  Matrix Statistics: '
      WRITE(*,*) ' '
      CALL Matrix_Statistics(A,B,BA,BAT,S_HG,tmpJacCC,Miter,LU_Miter)
      !STOP 'after writesparsematrix'
    END IF

    CALL Free_Matrix_CSR( tmpJacCC )
    CALL Free_SpRowColD( temp_LU_Dec )

    TimeSymbolic = MPI_WTIME() - StartTimer   ! stop timer for symbolic matrix calculations

    ! ---- Calculate first reaction rates
    IF (Teq) THEN
      CALL ReactionRatesAndDerivative_ChemKin( Tspan(1) , InitValAct , Rate , DRatedT )
    ELSE
      CALL ReactionRates_Tropos( Tspan(1) , InitValAct , Rate )
    END IF
    Y = MAX( ABS(InitValAct) , eps ) * SIGN( ONE , InitValAct )    ! |y| >= eps
    
    ! ---- Calculate values of Jacobian
    StartTimer  = MPI_WTIME()
    CALL Jacobian_CC( Jac_CC , BAT , A , Rate , Y )
    TimeJac     = TimeJac + MPI_WTIME() - StartTimer
    Output%npds = Output%npds + 1
  END IF

  !-----------------------------------------------------------------------
  !--- Routines for pathway analysis
  !-----------------------------------------------------------------------
  ALLOCATE(mixing_ratios_spc(nspc,3),integrated_rates(neq))
  mixing_ratios_spc = ZERO;  integrated_rates  = ZERO
  mixing_ratios_spc(:,1) = InitValAct
  
  !-----------------------------------------------------------------------
  ! --- Start the integration routine 
  !-----------------------------------------------------------------------
  CALL Integrate (  InitValAct     &  ! initial concentrations activ species
  &               , Rate           &  ! reaction rates at time=t0
  &               , Tspan          &  ! integration invervall
  &               , Atol           &  ! abs. tolerance gas species
  &               , RtolROW        &  ! rel. tolerance Rosenbrock method
  &               , ODEsolver      )  ! methode for solving the ode
  !---------------------------------------------------------------
  ! --- stop timer and print output statistics
  Timer_Finish = MPI_WTIME() - Timer_Start + Time_Read
  !
  ! Print statistics
  CALL Output_Statistics( Time_Read       , TimeSymbolic  , TimeFac       &
  &                     , TimeSolve       , TimeRates     , TimeJac       &
  &                     , TimeIntegrationE, Timer_Finish  , TimeRateSend  &
  &                     , TimeNetCDF      , TimeErrCalc   , TimeRhsCalc   )
  !---------------------------------------------------------------


  IF ( Lehmann ) THEN
    ! writing data for pahtway analysis
    mixing_ratios_spc(:,:) = (mixing_ratios_spc / (6.02e20_dp)) * 1.e09_dp ! in ppb
    integrated_rates = (integrated_rates / (6.02e20_dp)) * 1.e09_dp ! in ppb
    CALL WriteAnalysisFile(ReactionSystem,y_name,mixing_ratios_spc,integrated_rates)
  END IF

  ! --- Close MPI 
  CALL FinishMPI()

  !================================================================
  !==  FORMAT Statements
  !================================================================
  !
  798  FORMAT('      Sum Initval (gaseous)      =  ', Es8.2, A)
  799  FORMAT('      Sum Initval (aqueous)      =  ', Es8.2, A)
  800  FORMAT('      Sum Emissions (gaseous)    =  ', Es8.2, A)
  801  FORMAT('      Initial Temperature        =  ', Es8.2,'  [K]') 
  802  FORMAT('      Initial Pressure           =  ', Es8.2,'  [Pa]')
  803  FORMAT('      Reactor denstiy            =  ', Es8.2,'  [kg/cm3]')
END PROGRAM chemie
