MODULE InitRoutines_Mod   
   

  USE Kind_Mod
  USE Rates_Mod
  USE Sparse_Mod
  USE Chemsys_Mod
  USE Rosenbrock_Mod
  USE Control_Mod
  USE Reac_Mod
  USE MPI_Mod
  USE IO_Mod
  USE CombustionInput_Mod
  USE NetCDF_Mod
  USE Cycles_Mod
  USE fparser
  USE ISSA_Mod
  IMPLICIT NONE


  CONTAINS

    SUBROUTINE Initialize()

      CHARACTER(80)   :: Filename0 = ''        ! *.run file
      LOGICAL :: PrintToScreen = .True.


      INTEGER  :: i,nSchwefel

      ! NetCDF stuff
      REAL(dp) :: StartTimer
      REAL(dp), ALLOCATABLE :: Y(:)
      ! reaction rate array + part. derv. rate over temperatur vector
      REAL(dp), ALLOCATABLE :: Rate(:), DRatedT(:)  , wetRad(:)
      !
      CHARACTER(2)  :: newLinAlg  = ''
      !
      REAL(dp), ALLOCATABLE :: Atol(:)


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
      LOGICAL :: done=.FALSE.
      
      ! timer
      REAL(dp) :: t_1,t_2

      INTEGER  :: ICNTL(30)
      REAL(dp) :: RCNTL(30)

      !----------------------------------------------------------------
      ! --- Read run control parameters (which runfile)
      CALL getarg( 1 , FileName0 )             
      IF ( FileName0 == '' ) THEN
        WRITE(*,777,ADVANCE='NO') 'Input RUNFilename: '
        READ(*,*)   FileName0
      END IF
      RunFile = TRIM(ADJUSTL(FileName0))

      !================================================================
      !===                     Initialization
      !================================================================
      !
      !----------------------------------------------------------------
      ! --- Read run-file
      CALL InitRun
      IF (PrintToScreen) WRITE(*,777) 'Initialize run-file .......... done'
    
      !----------------------------------------------------------------
      ! --- set cloud intervall
      LWCb = Set_pseudoLWCbounds()

      !----------------------------------------------------------------
      !  --- read the .sys data, save coefs in sparse matrix
      Time_Read = MPI_WTIME()

      IF (PrintToScreen) WRITE(*,777,ADVANCE='NO') 'Reading sys-file .............'

      IF ( Combustion ) THEN

        CALL Read_Elements( SysFile , SysUnit )
        CALL Read_Species ( SysFile , SysUnit )
        CALL Read_Reaction( SysFile , SysUnit )

        IF (PrintToScreen) WRITE(*,*) 'done   ---->  Solve Gas Energy Equation '

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

        !--- gather correct indices cause Combustion mechanisms read in unsorted
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
          IF ( PrintToScreen ) THEN
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
        WRITE(*,*) 'done  ---->  Fix Temperature'
      
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
      IF (PrintToScreen) THEN
        WRITE(*,777) 'Reading ini-file ............. done'
        WRITE(*,777) 'Printing chem-file ........... done'
        WRITE(*,*)
        WRITE(*,'(10X,A,I6)') '    Number of Reactions = ', neq
        WRITE(*,'(10X,A,I6)') '    Number of Species   = ', nspc
      END IF
      !-----------------------------------------------------------------------
      ! --- Timers
      Time_Read   = MPI_WTIME() - Time_Read  ! Stop Input timer
      Timer_Start = MPI_WTIME()              ! Start timer for integration
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! --- Read species for diagnose (print species concs to NetCDF file)
      CALL Read_Diag( NetCDF%spc_Pos , NetCDF%spc_Phase , InitFile )
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! --- this is for the new mass action product routine 
      CALL Setup_SpeciesOrder(A)
      !-----------------------------------------------------------------------

      IF (KPP_Conversion) CALL SYS_TO_KPP(ReactionSystem)


      !-----------------------------------------------------------------------
      ! --- Dimension initialisation for the unknowns and matrices
      !
      nsr = nspc + neq

      IF ( Combustion ) THEN
        nDIM = nspc+1; nDIMcl = nspc+1; nDIMex = nsr+1
      ELSE
        nDIM = nspc; 	 nDIMcl = nspc;		nDIMex = nsr
      END IF

      rNspc = ONE/REAL(nspc,KIND=dp)  ! rNspc for error calculation
      !-----------------------------------------------------------------------

      ! true if mechanism contains photolytic reactions
      PHOTO = (nr_G_photo+nr_A_photo) > 0


      !----------------------------------------------------------------------------
      !--- If one wants to change the calculation of the linear systems without 
      !    opening the run-file, start the simulation with a second parameter
      CALL getarg ( 2 , newLinAlg )
      IF ( newLinAlg /= '' ) THEN
        SELECT CASE (newLinAlg)
          CASE ('cl'); CLASSIC = .TRUE.;  EXTENDED = .FALSE.; LinAlg = newLinAlg
          CASE ('ex'); CLASSIC = .FALSE.; EXTENDED = .TRUE.;  LinAlg = newLinAlg
          CASE DEFAULT;  WRITE(*,777) '    Linear Algebra either "cl" or "ex" !'
        END SELECT
      END IF
      !----------------------------------------------------------------------------

      
      IF ( PrintToScreen ) THEN
        WRITE(*,777) 'Initial values:   '
        WRITE(*,*)
        
        IF ( Combustion ) fmt0 = '  [mol/cm3]'
        IF ( hasGasSpc )   WRITE(*,798) 'gaseous', SUM(InitValAct( bGs(1):bGs(2) )) , fmt0
        IF ( hasAquaSpc )  WRITE(*,798) 'aqueous', SUM(InitValAct( bAs(1):bAs(2) )) , fmt0
        IF ( hasSolidSpc ) WRITE(*,798) 'solid  ', SUM(InitValAct( bSs(1):bSs(2) )) , fmt0
        IF ( hasPartiSpc ) WRITE(*,798) 'parti  ', SUM(InitValAct( bPs(1):bPs(2) )) , fmt0
        WRITE(*,800) SUM(y_emi) , '  [molec/cm3/s]'

        IF ( Combustion ) THEN
          WRITE(*,801) Temperature0
          WRITE(*,802) Pressure0
          WRITE(*,803) rho
        END IF
        WRITE(*,*)
      END IF


      !-----------------------------------------------------------------------
      ! ---  Allocate arrays for some varables
      ID_1 = SparseID( nDim )
      ALLOCATE( Y(nspc)   &
      &       , Rate(neq) , DRatedT(neq)  &
      &       , wetRad(nFrac)             )
      !-----------------------------------------------------------------------

      
      ! All species get the same abs tolerance
      ALLOCATE( ATolAll(nDIM) )
      IF ( hasGasSpc )  ATolAll(iGs) = ATolGas
      IF ( hasAquaSpc ) ATolAll(iAs) = ATolAqua
      
      IF ( Combustion ) ATolAll(nDIM) = ATolTemp
      
      IF ( Combustion ) THEN
        Atol = [AtolGas , AtolTemp]		! abs. tolerance for Combustion scenario
      ELSE
        Atol = [AtolGas]		! abs. tolerance for Tropos scenario
        IF ( hasAquaSpc ) Atol = [AtolGas , AtolAqua]		! abs. tolerance for Tropos scenario
      END IF

      !-----------------------------------------------------------------------
      ! --- Initialize NetCDF output file
      IF ( NetCdfPrint ) THEN 
        StartTimer = MPI_WTIME()
        CALL InitNetcdf
        CALL SetOutputNCDF( NetCDF, Tspan(1) , ZERO , InitValAct )
        CALL StepNetCDF( NetCDF )
        TimeNetCDF = MPI_WTIME() - StartTimer
      END IF
      !-----------------------------------------------------------------------

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



      !-----------------------------------------------------------------------
      ! if either Rosenbrock method or backward Euler was choosen
      IF ( INDEX(TRIM(ODEsolver(1:7)),'METHODS')>0) THEN
        !---- Get Rosenbrock Parameters
        CALL SetRosenbrockMethod( ROS , ODEsolver )  

        ! we need to calculate the Jacobian for both versions 'cl' and 'ex' to
        ! calculate an initial stepsize based on 2nd derivative (copy of MATLABs ode23s)
        ! symbolic mult for Jacobian Jac = [ BAT * Ones_nR * A * Ones_nS ], 
        ! add diagonal entries and set them to zero

        ! Set symbolic structure of iteration matrix for Row-Method
        ! also assign constant matrix parts like (beta-alpha)^T, alpha
        IF ( CLASSIC ) THEN
          CALL SymbolicMult( BAT , A , tmpJacCC )  
          Id = SparseID( nspc )                    
          CALL SymbolicAdd( Jac_CC , Id , tmpJacCC )
          CALL Free_Matrix_CSR( Id )
          CALL BuildSymbolicClassicMatrix(  Miter , Jac_CC  , ROS%ga )
          IF (DebugPrint) CM_1 = Copy_CSR(Jac_CC)
        ELSE !IF ( EXTENDED ) THEN
          CALL BuildSymbolicExtendedMatrix( Miter , A , BAT , ROS%ga ) 
        END IF

        ! Choose ordering/factorisation strategie and do symb LU fact
        IF ( useSparseLU ) THEN
          ! Permutation given by Markowitz Ordering strategie
          temp_LU_Dec = CSR_to_SpRowColD(Miter) 

          IF ( Ordering == 8 ) THEN
            CALL SymbLU_SpRowColD_M( temp_LU_Dec )
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

        CALL Free_Matrix_CSR( tmpJacCC )
        CALL Free_SpRowColD( temp_LU_Dec )

        TimeSymbolic = MPI_WTIME() - StartTimer   ! stop timer for symbolic matrix calculations
        ! ---- Calculate first reaction rates
        RateCnt = 0
        IF (Combustion) THEN
          CALL ReactionRates( Tspan(1) , [InitValAct,Temperature0] , Rate , DRatedT )
        ELSE
          CALL ReactionRates( Tspan(1) , InitValAct , Rate )
        END IF
        Y = MAX( ABS(InitValAct) , eps ) * SIGN( ONE , InitValAct )    ! |y| >= eps
        
        IF (CLASSIC) THEN
          ! ---- Calculate values of Jacobian
          StartTimer = MPI_WTIME()
          CALL Jacobian_CC( Jac_CC , BAT , A , Rate , Y )
          Out%npds = Out%npds + 1
          TimeJac  = TimeJac + MPI_WTIME() - StartTimer
        END IF
      END IF

      !================================================================
      !==  FORMAT Statements
      !================================================================
      !
      777  FORMAT(10X,A)
      798  FORMAT(10X,'    Sum Initval (',A7,')      =  ', Es16.8, A)
      800  FORMAT(10X,'    Sum Emissions (gaseous)    =  ', Es16.8, A)
      801  FORMAT(10X,'    Temperature                =  ', Es16.8,'  [K]') 
      802  FORMAT(10X,'    Pressure                   =  ', Es16.8,'  [Pa]')
      803  FORMAT(10X,'    Reactor denstiy            =  ', Es16.8,'  [kg/cm3]')
      
    END SUBROUTINE Initialize
   
    SUBROUTINE InitRun()
    !==================================================
    !===  Reading and Setting of Run Control Parameters
    !==================================================


      INTEGER        :: io_stat
      CHARACTER(400) :: io_msg = ''
      
!-----------------------------------------------------------------
!--- NAMELISTS
      NAMELIST /SCENARIO/  Bsp ,     &
      &                    WaitBar , &
      &                    Combustion , &
      &                    Simulation, &
      &                    Reduction, &
      &                    KPP_Conversion

      NAMELIST /FILES/  SysFile ,    &
      &                 DataFile ,   &
      &                 InitFile ,   &
      &                 MetFile ,    &
      &                 TargetFile , &
      &                 MWFile

      NAMELIST /TIMES/  tBegin , tEnd

      NAMELIST /METEO/  LWCLevelmin , &
      &                 LWCLevelmax , &
      &                 pHSet ,       &
      &                 iDate ,       &
      &                 rlat ,        &
      &                 rlon ,        &
      &                 dust ,        &
      &                 Temperature0 ,&
      &                 Pressure0
      
      NAMELIST /NUMERICS/  RtolROW ,     &
      &                    AtolGas ,     &
      &                    AtolAqua ,    & 
      &                    AtolTemp,     &
      &                    PI_StepSize , &
      &                    minStp ,      &
      &                    maxStp ,      &
      &                    LinAlg ,      &  
      &                    ODEsolver ,   &
      &                    Error_Est ,   &
      &                    Ordering ,    &
      &                    ParOrdering

      NAMELIST /OUTPUT/  NetCdfFile , &
      &                  StpNetcdf ,  &
      &                  StpFlux ,    &
      &                  nOutP ,      &
      &                  DebugPrint , &
      &                  MatrixPrint, &
      &                  FluxDataPrint

!
!===================================================================
!===  Set and Read Simulation Values
!===================================================================
!
!--- Open run control file
      OPEN(UNIT=RunUnit,FILE=RunFile,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'opening run-file')

!-----------------------------------------------------------------
!---  Scenario
!-----------------------------------------------------------------
!
!--- Set Default Values for SCENARIO Namelist

      WaitBar  = .TRUE.
      Combustion  = .FALSE.
      Simulation = .TRUE.
      Reduction  = .FALSE.
      KPP_Conversion = .FALSE.

!--- Read SCENARIO namelist
      READ(RunUnit,SCENARIO,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading SCENARIO list')
      
      IF (Combustion) Combustion = .TRUE.

!-----------------------------------------------------------------
!---  Files
!-----------------------------------------------------------------
      READ(RunUnit,FILES,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading FILES list')
      
!--- Adjust Filenames
      CALL FileNameCheck(SysFile,'SysFile')
      CALL FileNameCheck(DataFile,'DataFile')
      CALL FileNameCheck(InitFile,'InitFile')
      ChemFile   = ADJUSTL(SysFile(:INDEX(SysFile,'.sys')-1)//'.chem')
      MWFile     = ADJUSTL(MWFile)
      TargetFile = ADJUSTL(TargetFile)

      IF ( TRIM(BSP) == '' ) THEN
			  Bsp = ADJUSTL(SysFile(INDEX(SysFile,'/')+1:INDEX(SysFile,'.sys')-1))
			ELSE
        Bsp = ADJUSTL(Bsp)
			END IF

!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!
!--- Set Default Values for TIMES Namelist

      tBegin = 0.0_dp
      tEnd   = 0.0_dp

!--- Read TIMES namelist
      READ(RunUnit,TIMES,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading TIMES list')

			IF ( tBegin >= tEnd ) THEN
    	  WRITE(*,*) '  tBegin >= tEnd  '
    	   STOP 
			ELSE
				Tspan = [tBegin , tEnd]
    	END IF
!
!-----------------------------------------------------------------
!---  Meteorology
!-----------------------------------------------------------------
!
!--- Set Default Values for METEO Namelist
      pHSet       = .TRUE.
      LwcLevelmin = 2.0e-08_dp 
      LwcLevelmax = 3.0e-04_dp 
      constLWC    = .FALSE. 

      idate  = 011027
      rlat   = 50.65_dp
      rlon   = 10.77_dp
      Dust   = 1.0_dp

      Temperature0 = 280.0_dp
      Pressure0    = 200000.0_dp

      REWIND(RunUnit)
      READ(RunUnit,METEO,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading METEO list')

      IF ( LWCLevelmin ==  LWCLevelmax ) THEN
        constLWC = .TRUE.
        LWCconst = LWCLevelmin
      END IF

!
!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!
!--- Set Default Values for NUMERICS Namelist
      RtolROW     = 1.0e-5_dp   ! Relative tolerance For ROW
      AtolGas     = 1.0e-7_dp   ! Absolute tolerance for gas phase
      AtolAqua    = 1.0e-7_dp   ! Absolute tolerance for liquid phase
      AtolTemp    = 1.0e-7_dp   ! Absolute tolerance for temperature
      Error_Est   = 2           ! error estimation default 2-norm
      PI_StepSize = .FALSE.
      minStp      = 1.0e-20_dp  ! minimum timestep of ROW method in [sec]
      maxStp      = 250.0_dp     ! maximum timestep of ROW method in [sec]
      LinAlg      = 'cl'        ! method of solving linear algebra (classic)
      ODEsolver   = 'ROS34PW3'  ! ROW scheme
      Ordering    = 8           ! sparse LU, no numerical pivoting
      ParOrdering = -1          ! -1 = serial ordering, 0,1,2 = parallel ordering
      eps_red     = 0.11_dp
      
!--- Read NUMERICS namelist
      READ(RunUnit,NUMERICS,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading NUMERICS list')

      IF      ( LinAlg == 'cl' ) THEN
        CLASSIC  = .TRUE.
      ELSE IF ( LinAlg == 'ex' ) THEN
        EXTENDED = .TRUE.
      ELSE
        WRITE(*,*) '  Check run-file: LinAlg either "cl" or "ex" !'
				  STOP
      END IF
      

			 ! Test that Rosenbrock tolerance > 0
	    IF ( RtolROW <= ZERO ) THEN
	      WRITE(*,*) '  RtolROW must be positiv scalar!'
	        STOP
	    END IF

			! Test that absolute tolerance for gas and aqua species is > 0
			IF ( .NOT.(AtolGas*AtolAqua*AtolTemp) >= ZERO ) THEN
				WRITE(*,*) '  ATols must be positive!'
				 STOP
			END IF

			! Test if maximum stepsize is not to small/big
			IF ( maxStp <= ZERO ) THEN
			  WRITE(*,*) '  Maximum stepsize = ',maxStp, ' to low!'
				 STOP
			ELSE IF ( maxStp > tEnd-tBegin ) THEN
				WRITE(*,*) '  Maximum stepsize = ',maxStp, ' to high!'
				 STOP
			END IF
			IF ( minStp <= 1.e-50_dp ) THEN
			  WRITE(*,*) '  Minimums stepsize = ',minStp, ' to low!'
				 STOP
		  ELSE IF ( minStp > maxStp ) THEN
			  WRITE(*,*) '  Minimums stepsize = ', minStp, ' > ', maxStp
				 STOP
			END IF


      useSparseLU = .TRUE.


!-----------------------------------------------------------------
!---  Output of Data
!-----------------------------------------------------------------
!
!--- Set Default Values for OUTPUT Namelist
      StpNetcdf     = -1.0_dp      ! Time step for Netcdf output      [in sec]
      StpFlux       = -1.0_dp
      nOutP         = 100
      MatrixPrint   = .FALSE.
      DebugPrint    = .FALSE.
      NetCdfPrint   = .TRUE.
      FluxDataPrint = .FALSE.
!
!--- Read OUTPUT namelist
      READ(RunUnit,OUTPUT,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading OUTPUT list')

      NetCDFFile = ADJUSTL(NetCDFFile)
      IF ( TRIM(NetCdfFile) == '' ) NetCdfPrint = .FALSE.   ! no output if no filename is declared
      IF ( nOutP < 2 ) nOutP = 2                          ! minimum output steps are 2

      
   END SUBROUTINE InitRun


  SUBROUTINE InitReduction
    USE Control_Mod

    INTEGER, PARAMETER :: ReductionUnit = 112
    INTEGER        :: io_stat
    CHARACTER(400) :: io_msg = ''


    NAMELIST /SCENARIO/  TargetFile , &
     &                   FluxFile, &
     &                   FluxMetaFile, &
     &                   Red_TStart ,  &
     &                   Red_TEnd ,    &
     &                   eps_red 

    OPEN( FILE='REDUCTION/Reduction.init' , UNIT=ReductionUnit &
    &   , IOSTAT=io_stat , IOMSG=io_msg )

    CALL ErrorCheck(io_stat,io_msg,'opening reduction.init file')

    TargetFile   = ''
    FluxFile     = ''
    FluxMetaFile = ''
    Red_TStart   = 0.0d0 
    Red_TEnd     = 0.0d0 
    eps_red      = 0.11d0 


    READ(ReductionUnit,SCENARIO,IOSTAT=io_stat,IOMSG=io_msg)
    CALL ErrorCheck(io_stat,io_msg,'reading SCENARIO list')

    CALL FileNameCheck('REDUCTION/'//TRIM(TargetFile),'TargetFile')
    CALL FileNameCheck('REDUCTION/'//TRIM(FluxFile),'FluxDataFile')
    CALL FileNameCheck('REDUCTION/'//TRIM(FluxMetaFile),'FluxMetaDataFile')
  
    TargetFile   = 'REDUCTION/'//TRIM(TargetFile)
    FluxFile     = 'REDUCTION/'//TRIM(FluxFile)
    FluxMetaFile = 'REDUCTION/'//TRIM(FluxMetaFile)


    IF ( Red_TStart >= Red_TEnd ) THEN
   	  WRITE(*,*) '  tBegin >= tEnd  '
   	   STOP 
    END IF

    IF ( eps_red <= 0.0_dp ) THEN
      WRITE(*,*) '  reduction parameter eps_red <= 0  ---> increase value'
       STOP 
    ELSE IF ( eps_red > 1.0 ) THEN
      WRITE(*,*) '  reduction parameter  eps_red > 1  ---> decrease value'
       STOP 
    END IF


    CLOSE(ReductionUnit)

  END SUBROUTINE InitReduction



  SUBROUTINE ErrorCheck(io_stat,io_msg,cause)
    INTEGER      :: io_stat
    CHARACTER(*) :: io_msg, cause
    IF ( io_stat>0 ) WRITE(*,*) '   ERROR while '//cause//'  ::  ',io_stat,'  '//TRIM(io_msg)
  END SUBROUTINE ErrorCheck

  SUBROUTINE FileNameCheck(Name,miss)
    CHARACTER(*) :: Name
    CHARACTER(*) :: miss
    LOGICAL      :: ex

    INQUIRE(FILE=TRIM(Name), EXIST=ex)
    
    IF ( TRIM(Name) == '' .OR. .NOT.ex ) THEN
      WRITE(*,*); WRITE(*,*)
      WRITE(*,'(10X,A)') 'ERROR    Missing:  '//TRIM(miss)
      WRITE(*,'(10X,A)') '         FileName: '//TRIM(Name)
      WRITE(*,*); WRITE(*,*)
       STOP
    ELSE
      Name = ADJUSTL(Name)
    END IF
  END SUBROUTINE FileNameCheck

END MODULE InitRoutines_Mod
