MODULE InitRoutines_Mod   
   

   USE MPI_Mod
   IMPLICIT NONE


   CONTAINS
   
   SUBROUTINE InitRun()
!==================================================
!===  Reading and Setting of Run Control Parameters
!==================================================
      USE Control_Mod


      INTEGER        :: io_stat
      CHARACTER(400) :: io_msg = ''
      
!-----------------------------------------------------------------
!--- NAMELISTS
      NAMELIST /SCENARIO/  Bsp ,     &
      &                    WaitBar , &
      &                    ChemKin , &
      &                    Simulation , &
      &                    Reduction , &
      &                    KPP_Conversion , &
      &                    Lumping

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
      ChemKin  = .FALSE.
      Simulation = .TRUE.
      Reduction  = .FALSE.
      KPP_Conversion = .FALSE.
      Lumping = .FALSE.

!--- Read SCENARIO namelist
      READ(RunUnit,SCENARIO,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading SCENARIO list')
      
      IF (ChemKin) Teq = .TRUE.

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
        CALL FinishMPI(); STOP 
      ELSE
        Tspan = [tBegin , tEnd]
        Tspan_tot = [tBegin , tEnd]
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
				CALL FinishMPI();  STOP
      END IF
      

			 ! Test that Rosenbrock tolerance > 0
	    IF ( RtolROW <= ZERO ) THEN
	      WRITE(*,*) '  RtolROW must be positiv scalar!'
	      CALL FinishMPI();  STOP
	    END IF

			! Test that absolute tolerance for gas and aqua species is > 0
			IF ( .NOT.(AtolGas*AtolAqua*AtolTemp) >= ZERO ) THEN
				WRITE(*,*) '  ATols must be positive!'
				CALL FinishMPI(); STOP
			END IF

			! Test if maximum stepsize is not to small/big
			IF ( maxStp <= ZERO ) THEN
			  WRITE(*,*) '  Maximum stepsize = ',maxStp, ' to low!'
				CALL FinishMPI(); STOP
			ELSE IF ( maxStp > tEnd-tBegin ) THEN
				WRITE(*,*) '  Maximum stepsize = ',maxStp, ' to high!'
				CALL FinishMPI(); STOP
			END IF
			IF ( minStp <= 1.e-50_dp ) THEN
			  WRITE(*,*) '  Minimums stepsize = ',minStp, ' to low!'
				CALL FinishMPI(); STOP
		  ELSE IF ( minStp > maxStp ) THEN
			  WRITE(*,*) '  Minimums stepsize = ', minStp, ' > ', maxStp
				CALL FinishMPI(); STOP
			END IF


      ! setting up the factorization strategy
      !IF ( Ordering <  8 .OR. ParOrdering >= 0 ) useMUMPS = .TRUE.
      !IF ( Ordering >= 8 ) useSparseLU = .TRUE.
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
     !&                   FluxFile, &
     !&                   FluxMetaFile, &
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

    FluxFile     = 'flux_'//TRIM(BSP)//'.dat'
    FluxMetaFile = 'fluxmeta_'//TRIM(BSP)//'.dat'   

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
   	  CALL FinishMPI(); STOP 
    END IF

    IF ( eps_red <= 0.0_dp ) THEN
      WRITE(*,*) '  reduction parameter eps_red <= 0  ---> increase value'
      CALL FinishMPI(); STOP 
    ELSE IF ( eps_red > 1.0 ) THEN
      WRITE(*,*) '  reduction parameter  eps_red > 1  ---> decrease value'
      CALL FinishMPI(); STOP 
    END IF


    CLOSE(ReductionUnit)

  END SUBROUTINE InitReduction

  SUBROUTINE InitLumping
    USE Control_Mod

    INTEGER, PARAMETER :: LumpingUnit = 126
    INTEGER        :: io_stat
    CHARACTER(400) :: io_msg = ''


    NAMELIST /SCENARIO/  PreserveFile, &
     &                   eps_k , &
     &                   eps_tau 


    OPEN( FILE='LUMPING/Lumping.init' , UNIT=LumpingUnit &
    &   , IOSTAT=io_stat , IOMSG=io_msg )

    CALL ErrorCheck(io_stat,io_msg,'opening Lumping.init file')

    ! default values
    PreserveFile = ''
    eps_k        = rTWO
    eps_tau      = ONE

    READ(LumpingUnit,SCENARIO,IOSTAT=io_stat,IOMSG=io_msg)
    CALL ErrorCheck(io_stat,io_msg,'reading SCENARIO list')

    CALL FileNameCheck('LUMPING/'//TRIM(PreserveFile),'PreserveFile')

    PreserveFile = 'LUMPING/'//TRIM(PreserveFile)

    !IF ( eps_red <= 0.0_dp ) THEN
    !  WRITE(*,*) '  reduction parameter eps_red <= 0  ---> increase value'
    !  CALL FinishMPI(); STOP 
    !ELSE IF ( eps_red > 1.0 ) THEN
    !  WRITE(*,*) '  reduction parameter  eps_red > 1  ---> decrease value'
    !  CALL FinishMPI(); STOP 
    !END IF

    CLOSE(LumpingUnit)

  END SUBROUTINE InitLumping

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
      CALL FinishMPI(); STOP
    ELSE
      Name = ADJUSTL(Name)
    END IF
  END SUBROUTINE FileNameCheck

END MODULE InitRoutines_Mod
