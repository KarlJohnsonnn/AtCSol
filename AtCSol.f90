!=======================================
!=======================================
! vorlaeufiges Hauptprogramm fuer Tests
! box model chemical kinetics simulation
!=======================================
!=======================================
!
!
PROGRAM AtCSol
  !
  USE Kind_Mod
  USE Rates_Mod
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
  USE fparser
  USE issa
  IMPLICIT NONE
  !
  CHARACTER(80)   :: Filename0 = ''        ! *.run file
  !
  REAL(dp), ALLOCATABLE :: Atol(:)
  INTEGER  :: i,nSchwefel

  ! NetCDF stuff
  REAL(dp) :: StartTimer
  INTEGER  :: errind(1,1)
  REAL(dp), ALLOCATABLE :: ErrVals(:), Y(:)
  ! reaction rate array + part. derv. rate over temperatur vector
  REAL(dp), ALLOCATABLE :: Rate(:), DRatedT(:)  , wetRad(:)
  !
  ! use other ROS methode or tolerance
  CHARACTER(80) :: neuTolR  = ''
  CHARACTER(80) :: neuTolA  = ''
  CHARACTER(80) :: neuROW   = ''
  CHARACTER(2)  :: newLinAlg  = ''
  CHARACTER(1)  :: simul  = ''


  ! convertion from mole to mass to conc
  REAL(dp), ALLOCATABLE   :: MoleFrac(:), MassFrac(:), MoleConc(:)
  REAL(dp)                :: Press_in_dyncm2

  ! matricies for symbolic phase
  TYPE(CSR_Matrix_T)     :: Id , tmpJacCC  ! compressed row
  TYPE(SpRowColD_T)      :: temp_LU_Dec,temp_LU_Dec2    ! sparse-LU matrix format
  ! permutation vector/ pivot order for LU-decomp
  INTEGER, ALLOCATABLE :: InvPermu(:)
  INTEGER, ALLOCATABLE :: PivOrder(:)

  ! format string
  CHARACTER(14) :: fmt0 = '  [molec/cm3] '
  CHARACTER(8) ::  unit  = ''

  ! new testing stuff 
  TYPE(CSR_Matrix_T) :: A_T, P_HG, S_HG, R_HG
  CHARACTER(1) :: inpt=''
  TYPE(List), ALLOCATABLE :: ReacCyc(:)
  INTEGER,    ALLOCATABLE :: Target_Spc(:) 
  
  ! timer
  REAL(dp) :: t_1,t_2

  LOGICAL :: more_than_one=.FALSE.

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
    WRITE(*,777,ADVANCE='NO') 'Input RUNFilename: '; READ(*,*)   FileName0
  END IF
  RunFile = TRIM(ADJUSTL(FileName0))

  !================================================================
  !===                     Initialization
  !================================================================
  !
  !----------------------------------------------------------------
  ! --- Read run-file
  CALL InitRun
  IF ( MPI_master ) WRITE(*,777) 'Initialize run-file .......... done'
 
  !----------------------------------------------------------------
  ! --- set cloud intervall
  LWCb = Set_pseudoLWCbounds()

  !----------------------------------------------------------------
  !  --- read the .sys data, save coefs in sparse matrix
  Time_Read = MPI_WTIME()

  IF ( MPI_master ) WRITE(*,777,ADVANCE='NO') 'Reading sys-file .............'

  IF ( ChemKin ) THEN

    CALL Read_Elements( SysFile , SysUnit )
    CALL Read_Species ( SysFile , SysUnit )
    CALL Read_Reaction( SysFile , SysUnit )

    IF ( MPI_master ) WRITE(*,*) 'done   ---->  Solve Gas Energy Equation '

   !-----------------------------------------------------------------------
   ! --- print reactions and build A, B and (B-A) structure
   !-----------------------------------------------------------------------
    CALL Print_ChemFile   ( ReactionSystem , ChemFile , ChemUnit , .TRUE. )

    CALL GetSpeciesNames( ChemFile , y_name )
    CALL Read_ThermoData( SwitchTemp , DataFile , DataUnit , nspc )

    ! set boundaries for gaseous species
    bGs(1) = 1
    bGs(2) = ns_GAS
    iGs = [(i, i=bGs(1),bGs(2))]
    hasGasSpc = .TRUE.

    !--- gather correct indices cause ChemKin mechanisms read in unsorted
    CALL Setup_ThirdBodyIndex
    CALL Setup_ReactionIndex
   
    !--- Read initial values
    ALLOCATE( InitValAct(ns_GAS) , InitValKat(ns_KAT) , y_emi(ns_GAS) , y_depos(ns_GAS))
    y_emi   = ZERO
    y_depos = ZERO 

		!--- malloc gibbs energy, derivates
    ALLOCATE( GFE(nspc)   , DGFEdT(nspc)   &
    &       , DelGFE(neq) , DDelGFEdT(neq) )
		GFE      = ZERO; 	DGFEdT    = ZERO
		DelGFE   = ZERO; 	DDelGFEdT = ZERO

    IF ( MWFile /= '' ) THEN
      CALL Read_MolecularWeights(MW,MWFile,MWUnit,nspc)

      rMW = [ONE / MW]
      
      ALLOCATE( MoleFrac(ns_GAS) , MassFrac(ns_GAS) )
      MoleFrac = ZERO;	MassFrac  = ZERO  

      MoleFrac = 1.0e-20_dp
      CALL Read_INI_file( InitFile , MoleFrac, InitValKat , 'GAS' , 'INITIAL' )

      !Press = Pressure0               ! initial pressure in [Pa]
      Press_in_dyncm2 = Pressure0 * Pa_to_dyncm2

      MassFrac = MoleFr_To_MassFr  ( MoleFrac ) 
      MoleConc = MoleFr_To_MoleConc( MoleFrac                &
      &                            , Press = Press_in_dyncm2 &
      &                            , Temp  = Temperature0    )
    ELSE
      IF (MPI_master ) THEN
        WRITE(*,*)
        WRITE(*,777) '    No molecular weights are given.  '
        WRITE(*,777) '    Make sure the initial values are given in [mole/cm3] !'
        WRITE(*,*)
      END IF
      ALLOCATE( MoleConc(ns_GAS) )
      MoleConc = 1.0e-20_dp
      CALL Read_INI_file( InitFile , MoleConc, InitValKat , 'GAS' , 'INITIAL' )
    END IF
    !CALL Read_EMISS( InitFile , y_emi )
    
    ! Initialising reactor density
    rho  = Density( MoleConc )
    rRho = mega/rho       ! in [cm3/g]
    InitValAct = [MoleConc , Temperature0]

    

  ELSE

    CALL ReadSystem( SysFile )
    IF ( MPI_master ) WRITE(*,*) 'done  ---->  Fix Temperature'
  
   !-----------------------------------------------------------------------
   ! --- print reactions and build A, B and (B-A) structure
   !-----------------------------------------------------------------------
    CALL Print_ChemFile( ReactionSystem , ChemFile , ChemUnit , .FALSE. )

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

  END IF

  IF (MPI_ID==0) THEN
    WRITE(*,777) 'Reading ini-file ............. done'
    WRITE(*,777) 'Printing chem-file ........... done'
    WRITE(*,*)
    WRITE(*,'(10X,A,I6)') '    Number of Reactions = ', neq
    WRITE(*,'(10X,A,I6)') '    Number of Species   = ', nspc
  END IF


  !-----------------------------------------------------------------------
  ! --- Read species for diagnose (print species concs to NetCDF file)
  CALL Read_Diag( NetCDF%spc_Pos , NetCDF%spc_Phase , InitFile )
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
    nDIM = nspc+1; nDIMcl = nspc+1; nDIMex = nsr+1
  ELSE
    nDIM = nspc; 	 nDIMcl = nspc;		nDIMex = nsr
  END IF

  rNspc = ONE/REAL(nspc,KIND=dp)  ! rNspc for error calculation

  ! true if mechanism contains photolytic reactions
  PHOTO = (nr_G_photo+nr_A_photo) > 0

  !----------------------------------------------------------------------------
  !--- If one wants to change the calculation of the linear systems without 
  !    opening the run-file, start the simulation with a second parameter
  !----------------------------------------------------------------------------
  CALL getarg ( 2 , newLinAlg )
  IF ( newLinAlg /= '' ) THEN
    SELECT CASE (newLinAlg)
      CASE ('cl'); CLASSIC = .TRUE.;  EXTENDED = .FALSE.; LinAlg = newLinAlg
      CASE ('ex'); CLASSIC = .FALSE.; EXTENDED = .TRUE.;  LinAlg = newLinAlg
      CASE DEFAULT;  WRITE(*,777) '    Linear Algebra either "cl" or "ex" !'
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
  !--- print input parameter, method, tols, concs, temp, pressure,....
  !-----------------------------------------------------------------------
  CALL Print_Run_Param
  
  IF (MPI_ID==0) THEN
    WRITE(*,777) 'Initial values:   '
    WRITE(*,*)
    
    IF ( Teq ) fmt0 = '  [mol/cm3]'
    IF ( hasGasSpc )   WRITE(*,798) 'gaseous', SUM(InitValAct( bGs(1):bGs(2) )) , fmt0
    IF ( hasAquaSpc )  WRITE(*,798) 'aqueous', SUM(InitValAct( bAs(1):bAs(2) )) , fmt0
    IF ( hasSolidSpc ) WRITE(*,798) 'solid  ', SUM(InitValAct( bSs(1):bSs(2) )) , fmt0
    IF ( hasPartiSpc ) WRITE(*,798) 'parti  ', SUM(InitValAct( bPs(1):bPs(2) )) , fmt0
    WRITE(*,800) SUM(y_emi) , '  [molec/cm3/s]'

    IF ( Teq ) THEN
      WRITE(*,801) Temperature0
      WRITE(*,802) Pressure0
      WRITE(*,803) rho
    END IF
    WRITE(*,*)
  ELSE
    ! if more than 1 process, output data only on master process
    WaitBar     = .FALSE.
    NetCdfPrint = .FALSE.
  END IF

  !-----------------------------------------------------------------------
  ! ---  Allocate arrays for some varables
  !-----------------------------------------------------------------------
	ID_1 = SparseID( nDim )
  ALLOCATE( ErrVals(nspc) , Y(nspc)   &
  &       , Rate(neq) , DRatedT(neq)  &
  &       , wetRad(nFrac)             )

	
	! All species get the same abs tolerance
	ALLOCATE( ATolAll(nDIM) )
	IF ( hasGasSpc )  ATolAll(iGs) = ATolGas
  IF ( hasAquaSpc ) ATolAll(iAs) = ATolAqua
  
	IF ( Teq ) ATolAll(nDIM) = ATolTemp
  
	IF ( Teq ) THEN
		Atol = [AtolGas , AtolTemp]		! abs. tolerance for ChemKin scenario
	ELSE
		Atol = [AtolGas]		! abs. tolerance for Tropos scenario
		IF ( hasAquaSpc ) Atol = [AtolGas , AtolAqua]		! abs. tolerance for Tropos scenario
  END IF

  !-----------------------------------------------------------------------
  ! --- Initialize NetCDF output file
  !-----------------------------------------------------------------------
  IF ( NetCdfPrint ) THEN 
    StartTimer = MPI_WTIME()
    errind = 0
    CALL InitNetcdf
    CALL SetOutputNCDF( NetCDF, Tspan(1) , ZERO , errind , ZERO , InitValAct , Temperature0 )
    CALL StepNetCDF( NetCDF )
    TimeNetCDF = MPI_WTIME() - StartTimer
  END IF

  !-----------------------------------------------------------------------
  ! --- Symbolic phase 
  !    - generating sparse matrecies of stoechiometic coefs (beta-alpha)^T 
  !    - generating sparse matrecies of Jacobian matrix (if nessesarry)
  !    - generating sparse matrecies of LU decomposition of the Jacobian
  !-----------------------------------------------------------------------
  WRITE(*,777,ADVANCE='NO') 'Symbolic-phase................'
  StartTimer = MPI_WTIME()              ! start timer for symb phase
  
  CALL SymbolicAdd( BA , B , A )      	! symbolic addition:    BA = B + A
  CALL SparseAdd  ( BA , B , A, '-' )   ! numeric subtraction:  BA = B - A
  CALL TransposeSparse( BAT , BA )    	! transpose BA:        BAT = Transpose(BA) 

  

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

    ALLOCATE(maxErrorCounter(nDIM))
    maxErrorCounter = 0

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
    IF ( useSparseLU ) THEN
      ! Permutation given by Markowitz Ordering strategie
      temp_LU_Dec = CSR_to_SpRowColD(Miter) 


      !temp_LU_Dec2 = Copy_SpRowColD(temp_LU_Dec) 

      IF ( Ordering == 8 ) THEN

        ! test for better permutation 

        CALL SymbLU_SpRowColD_M( temp_LU_Dec )

        !IF (more_than_one) THEN
        !  DO i=1,10
        !    WRITE(*,*) ' nonzeros after fact :: ' , i , temp_LU_Dec2%nnz
        !    
        !    CALL Free_SpRowColD ( temp_LU_Dec2 )
        !
        !    temp_LU_Dec2 = Copy_SpRowColD(temp_LU_Dec) 
        !    
        !    CALL SymbLU_SpRowColD_M ( temp_LU_Dec2 )        
        !  END DO
        !END IF
      
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
      LU_Miter = RowColD_to_CSR( temp_LU_Dec , nspc , neq )

      ! Get the permutation vector LU_Perm and map values of Miter
      ! to the permuted LU matrix
      CALL Get_LU_Permutaion( LU_Perm , LU_Miter , Miter , nspc , neq )
      
      ! For the extended case one can save the values of alpha and 
      ! (beta-alpha)^T and just copy them for each iteration in ROS
      IF ( EXTENDED ) LUvalsFix = LU_Miter%Val
    END IF

    IF (MPI_ID==0) THEN
      
      WRITE(*,*) 'done'     ! Symbolic phase
      WRITE(*,*) 

      IF (MatrixPrint) THEN
        CALL WriteSparseMatrix(A,'MATRICES/alpha_'//TRIM(BSP), neq, nspc)
        CALL WriteSparseMatrix(B,'MATRICES/beta_'//TRIM(BSP), neq, nspc)
        CALL WriteSparseMatrix(tmpJacCC,'MATRICES/JAC_'//TRIM(BSP), neq, nspc)
        CALL WriteSparseMatrix(BA,'MATRICES/BA_'//TRIM(BSP), neq, nspc)
        CALL WriteSparseMatrix(Miter,'MATRICES/Miter_'//TRIM(BSP), neq, nspc)
        CALL WriteSparseMatrix(LU_Miter,'MATRICES/LU_Miter_'//TRIM(BSP), neq, nspc)
        WRITE(*,777,ADVANCE='NO') '  Continue? [y/n]';  READ(*,*) inpt
        IF (inpt/='y') STOP
      END IF

      WRITE(*,777) 'Matrix Statistics: '
      WRITE(*,*) 
      CALL Matrix_Statistics(A,B,BA,BAT,S_HG,tmpJacCC,Miter,LU_Miter)
    END IF

    CALL Free_Matrix_CSR( A_T )
    CALL Free_Matrix_CSR( P_HG )
    CALL Free_Matrix_CSR( S_HG )
    CALL Free_Matrix_CSR( R_HG )
    CALL Free_Matrix_CSR( tmpJacCC )
    CALL Free_SpRowColD( temp_LU_Dec )

    TimeSymbolic = MPI_WTIME() - StartTimer   ! stop timer for symbolic matrix calculations
    ! ---- Calculate first reaction rates
    RateCnt = 0
    IF (Teq) THEN
      CALL ReactionRates( Tspan(1) , [InitValAct,Temperature0] , Rate , DRatedT )
    ELSE
      CALL ReactionRates( Tspan(1) , InitValAct , Rate )
    END IF
    Y = MAX( ABS(InitValAct) , eps ) * SIGN( ONE , InitValAct )    ! |y| >= eps
    
    ! ---- Calculate values of Jacobian
    StartTimer = MPI_WTIME()
    CALL Jacobian_CC( Jac_CC , BAT , A , Rate , Y )
    Out%npds = Out%npds + 1
    TimeJac  = TimeJac + MPI_WTIME() - StartTimer
  END IF

  !WRITE(*,777,ADVANCE='NO') 'Simulation? [y/n]   ';  READ(*,*) simul
  !IF ( simul=='y' ) THEN
    ! open file to save the fluxes 
    IF ( MPI_master .AND. FluxAna ) THEN
      iStpFlux = 0
      CALL OpenFile_wStream(FluxUnit,FluxFile);       CLOSE(FluxUnit)
      CALL OpenFile_wSeq(FluxMetaUnit,FluxMetaFile);  CLOSE(FluxMetaUnit)
    END IF
    
    !-----------------------------------------------------------------------
    ! --- Start the integration routine 
    !-----------------------------------------------------------------------
    CALL Integrate ( InitValAct &  ! initial concentrations activ species
    &              , Rate       &  ! reaction rates at time=t0
    &              , Tspan      &  ! integration invervall
    &              , Atol       &  ! abs. tolerance of species
    &              , RtolROW    &  ! rel. tolerance Rosenbrock method
    &              , ODEsolver  )  ! methode for solving the ode system
    
    ! --- stop timer and print output statistics
    Timer_Finish = MPI_WTIME() - Timer_Start + Time_Read
    
    CALL Output_Statistics

    IF (NetCdfPrint) CALL TikZ_finished

  !END IF


  !************************************************************************************************
  !************************************************************************************************
  IF ( MPI_master .AND. FluxAna ) THEN
    StartTimer = MPI_WTIME()
    CALL Logo2()

    DO i=1,nr           ! Name,iReac,Mech,Class,Type,Param)
      CALL WriteReaction( TRIM(ReactionSystem(i)%Line1)       &
      &                 , i                                   &
      &                 , TRIM(BSP)                           &
      &                 , TRIM(ReactionSystem(i)%Type)        &
      &                 , TRIM(ReactionSystem(i)%Line3)       )
    END DO

    !-----------------------------------------------------------------------
    ! --- Read species groups in order to combine them (Species Lumping)
    IF( TargetFile/='' ) THEN
      CALL ISSA_structure( ReacCyc , Target_Spc , TargetFile )
      CALL ISSA_screening( ReactionSystem, ReacCyc, Target_Spc )
    ELSE
      WRITE(*,777) '    ** NO IMPORTANT SPECIES ARE DECLARED **'
    END IF
    TimeReduction = MPI_WTIME()-StartTimer
    CALL ConvertTime(TimeReduction,unit)
    WRITE(*,'(32X,A,1X,F10.4,A)') 'Time ISSA reduction = ', TimeReduction, unit
    WRITE(*,*); WRITE(*,*)
  END IF

  ! --- Close MPI 
  IF ( MPI_master ) THEN
    WRITE(*,*); WRITE(*,*)
    WRITE(*,777) '************ ********** *********** ********** ************'
    WRITE(*,777) '************ **********     DONE    ********** ************'
    WRITE(*,777) '************ ********** *********** ********** ************'
    WRITE(*,*); WRITE(*,*)
  END IF

  CALL ShowMaxErrorCounter()

  CALL FinishMPI()

  !================================================================
  !==  FORMAT Statements
  !================================================================
  !
  777  FORMAT(10X,A)
  798  FORMAT(10X,'    Sum Initval (',A7,')      =  ', Es8.2, A)
  800  FORMAT(10X,'    Sum Emissions (gaseous)    =  ', Es8.2, A)
  801  FORMAT(10X,'    Temperature                =  ', Es8.2,'  [K]') 
  802  FORMAT(10X,'    Pressure                   =  ', Es8.2,'  [Pa]')
  803  FORMAT(10X,'    Reactor denstiy            =  ', Es8.2,'  [kg/cm3]')
END PROGRAM AtCSol
