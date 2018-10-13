MODULE IO_Mod
  IMPLICIT NONE
  !
  CONTAINS
  SUBROUTINE Logo()
    USE MPI_Mod
      WRITE(*,777) ! in Georgi16
      WRITE(*,777) '************ ********** ************ ********** ************'
      WRITE(*,777) "*                                                          *"
      WRITE(*,777) "*         _                ____      ____            ___   *"
      WRITE(*,777) "*        dM.              6MMMMb/   6MMMMb\          `MM   *"
      WRITE(*,777) "*       ,MMb      /      8P    YM  6M'    `           MM   *"
      WRITE(*,777) "*       d'YM.    /M     6M      Y  MM         _____   MM   *"
      WRITE(*,777) "*      ,P `Mb   /MMMMM  MM         YM.       6MMMMMb  MM   *"
      WRITE(*,777) "*      d'  YM.   MM     MM          YMMMMb  6M'   `Mb MM   *"
      WRITE(*,777) "*     ,P   `Mb   MM     MM              `Mb MM     MM MM   *"
      WRITE(*,777) "*     d'    YM.  MM     MM               MM MM     MM MM   *"
      WRITE(*,777) "*    ,MMMMMMMMb  MM     YM      6        MM MM     MM MM   *"
      WRITE(*,777) "*    d'      YM. YM.  ,  8b    d9  L    ,M9 YM.   ,M9 MM   *"
      WRITE(*,777) "*  _dM_     _dMM_ YMMM9   YMMMM9   MYMMMM9   YMMMMM9 _MM_  *"
      WRITE(*,777) "*                                                          *"
      WRITE(*,777) '************ ********** ************ ********** ************'
      WRITE(*,*) 
      WRITE(*,777) "               Atmospheric Chemistry Solver                 " 
      WRITE(*,*) 
      WRITE(*,777) "      An experimental Fortran program for the analysis      "
      WRITE(*,777) "         of complex chemical multiphase mechanisms          " 
      WRITE(*,*) ; WRITE(*,*)
    777 FORMAT(10X,A)
  END SUBROUTINE Logo

  SUBROUTINE Logo2()
    USE MPI_Mod
      WRITE(*,*)
      WRITE(*,777) '******** ********** ************ ********** ********'
      WRITE(*,777) "*                                                  *"
      WRITE(*,777) "*       `MM'  6MMMMb\   6MMMMb\       dM.          *"
      WRITE(*,777) "*        MM  6M'    `  6M'    `      ,MMb          *"
      WRITE(*,777) "*        MM  MM        MM            d'YM.         *"
      WRITE(*,777) "*        MM  YM.       YM.          ,P `Mb         *"
      WRITE(*,777) "*        MM   YMMMMb    YMMMMb      d'  YM.        *"
      WRITE(*,777) "*        MM       `Mb       `Mb    ,P   `Mb        *"
      WRITE(*,777) "*        MM        MM        MM    d'    YM.       *"
      WRITE(*,777) "*        MM        MM        MM   ,MMMMMMMMb       *"
      WRITE(*,777) "*        MM  L    ,M9  L    ,M9   d'      YM.      *"
      WRITE(*,777) "*       _MM_ MYMMMM9   MYMMMM9  _dM_     _dMM_     *"
      WRITE(*,777) "*                                                  *"
      WRITE(*,777) '******** ********** ************ ********** ********'
      WRITE(*,*)
    777 FORMAT(10X,A)
  END SUBROUTINE Logo2
  !
  !
  SUBROUTINE Print_Run_Param()
    USE MPI_Mod
    USE Control_Mod
    USE Reac_Mod

    IF ( INDEX(SysFile,'.sys')==0)  SysFile = TRIM(SysFile)//'.sys'

    IF (Simulation) THEN
      WRITE(*,*)
      WRITE(*,777)   'Run - Paramter:'
      WRITE(*,*)
      WRITE(*,777)   '    Mechanism:             '//TRIM(SysFile)
      IF (NetCdfFile /= '') THEN
        WRITE(*,777)   '    NetCDF-File:           '//TRIM(NetCdfFile)
      ELSE
        WRITE(*,777)   '    NetCDF-File:           *** no NetCDF output ***'
      END IF
      WRITE(*,777)   '    Initials:              '//TRIM(InitFile)
      WRITE(*,777)   '    ODE solver:            '//TRIM(ODEsolver)
      IF (ODEsolver/='LSODE') THEN
        IF ( CLASSIC ) THEN 
          WRITE(*,777)   '    Linear Algebra:        Classic'
        ELSE
          WRITE(*,777)   '    Linear Algebra:        Extended'
        END IF
        IF (Error_Est==2) THEN
          WRITE(*,777)   '    Error Estimation:      Euklid Norm'
        ELSE
          WRITE(*,777)   '    Error Estimation:      Maximum Norm'
        END IF
        WRITE(*,777)   '    Solve Linear Systems:  Sparse LU, Markowitz Ordering Algorithm'
      END IF
      WRITE(*,777)
      WRITE(*,777)   'Tolerance:   '
      WRITE(*,777)
      WRITE(*,'(10X,A,2X,Es8.2)')   '    Relative Rosenbrock        = ',RtolROW
      WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute (gaseous species) = ',AtolGas
      IF (ns_AQUA>0) WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute (aqueous species) = ',AtolAqua
      IF ( Combustion ) THEN
        WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute Temperature       = ',AtolTemp
      END IF
      WRITE(*,*) 
    END IF
    777 FORMAT(10X,A)
  END SUBROUTINE Print_Run_Param
  !
  !
  SUBROUTINE Output_Statistics
    USE Kind_Mod
    USE MPI_Mod
    USE Control_Mod
    !
    REAL(dp) :: maxTRead,maxTSymb,maxTFac,maxTSolve,maxTRates,maxTJac &
    &               , maxTInte,maxTAll,maxTSend,maxtNcdf,maxTErr,maxTRhs,maxTFlux
    CHARACTER(8) :: unit(12)
    !
    maxTRead  = Time_Read
    maxTRates  = TimeRates
    maxTSymb  = TimeSymbolic
    maxTInte  = TimeIntegration
    maxTFac  = TimeFac
    maxTSolve  = TimeSolve
    maxTJac  = TimeJac
    maxTAll  = Timer_Finish
    maxTNcdf  = TimeNetCDF
    maxTErr  = TimeErrCalc
    maxTRhs  = TimeRhsCalc
    maxTFlux  = TimeFluxWrite
    maxTInte = maxTInte - maxTNcdf - maxTFlux

    CALL ConvertTime(maxTRead,unit(1)(:))
    CALL ConvertTime(maxTSymb,unit(2)(:))
    CALL ConvertTime(maxTNcdf,unit(3)(:))
    CALL ConvertTime(maxTFlux,unit(4)(:))
    CALL ConvertTime(maxTFac,unit(5)(:))
    CALL ConvertTime(maxTRhs,unit(6)(:))
    CALL ConvertTime(maxTSolve,unit(7)(:))
    CALL ConvertTime(maxTRates,unit(8)(:))
    CALL ConvertTime(maxTJac,unit(9)(:))
    CALL ConvertTime(maxTErr,unit(10)(:))
    CALL ConvertTime(maxTInte,unit(11)(:))
    CALL ConvertTime(maxTAll,unit(12)(:))
    !
      ! print the statistics
      299 format(10X,A,3X,F13.6,A)
      298 format(10X,A,3X,I10)
      777 FORMAT(10X,A)
      WRITE(*,*);  WRITE(*,*);    WRITE(*,*)
      WRITE(*,777) 'Statistics (Numbers):'; 
      WRITE(*,*)
      WRITE(*,298) '    successful time steps   =', Out%nsteps
      WRITE(*,298) '    failed time step        =', Out%nfailed
      WRITE(*,298) '    rate evaluations        =', Out%nRateEvals
      WRITE(*,298) '    Jacobian calculations   =', Out%npds
      WRITE(*,298) '    LU factorisations       =', Out%ndecomps
      WRITE(*,298) '    solved linear systems   =', Out%nsolves
      WRITE(*,*);  WRITE(*,*)
      WRITE(*,777)   'Statistics (Time):'
      WRITE(*,*)
      WRITE(*,299) '    reading mechanism       =', maxTRead,unit(1)
      WRITE(*,299) '    symbolic phase          =', maxTSymb,unit(2)
      WRITE(*,299) '    writing NetCDF-File     =', maxTNcdf,unit(3)
      WRITE(*,299) '    writing flux-dataset    =', maxTFlux,unit(4) ; WRITE(*,*)
      !WRITE(*,777) '    ------------------------+----------------------------'
      WRITE(*,299) '            factorisation   =', maxTFac  ,unit(5)
      WRITE(*,299) '          + right-hand side =', maxTRhs  ,unit(6)
      WRITE(*,299) '          + linear systems  =', maxTSolve,unit(7)
      WRITE(*,299) '          + reaction rates  =', maxTRates,unit(8)
      WRITE(*,299) '          + Jacobian        =', maxTJac  ,unit(9)
      WRITE(*,299) '          + error calc      =', maxTErr  ,unit(10)
      WRITE(*,777) '    ------------------------=----------------------'
      WRITE(*,299) '    integration             =', maxTInte,unit(11); WRITE(*,*)
      WRITE(*,299) '    total runtime           =', maxTAll,unit(12)
      WRITE(*,*);  WRITE(*,*);  WRITE(*,*)
  END SUBROUTINE
  
 
  SUBROUTINE SaveMatricies(aMat,bMat,cMat,dMat,eMat,fName)
    USE MPI_Mod
    USE Sparse_Mod

    TYPE(CSR_Matrix_T)   :: aMat,bMat,cMat,dMat,eMat
    CHARACTER(*)         :: fName
    !
    ! local stuff
    INTEGER, ALLOCATABLE :: InvPermu(:)
    INTEGER              :: mUnit, i
    CHARACTER(50)        :: mName
    !
    !
      ! only if MatrixPrint=True
      CALL WriteSparseMatrix(aMat,TRIM('matrixOut/alpha'//fName))
      CALL WriteSparseMatrix(bMat,TRIM('matrixOut/beta'//fName))
      CALL WriteSparseMatrix(cMat,TRIM('matrixOut/_beta-alpha_T'//fName))
      CALL WriteSparseMatrix(dMat,TRIM('matrixOut/Miter0'//fName))
      
      ! ordering>=8 --> Markowitz count (early minimum degree)
      ! Print structure of LU matrix and Permutation vector
      !
      CALL WriteSparseMatrix(eMat,TRIM('matrixOut/LUmiterStructure'//fName))
      CALL PrintPerm(eMat%Permu,eMat%InvPer,TRIM('matrixOut/Permu'//fName))
      !
      STOP 'MatrixPrint=True --> stop nach print in Integration_Mod'
  END SUBROUTINE SaveMatricies
  !
  !
  SUBROUTINE DebugPrint1(yvec,rvec,stepsize,time)
    USE Kind_Mod
    REAL(dp) :: yvec(:),rvec(:)
    REAL(dp) :: stepsize, time
    !
    WRITE(*,*) '----------------------------'
    WRITE(*,*) 'debug h, t      :: ', stepsize , time
    WRITE(*,*) '      rate(1:3) :: ', rvec(1:3)
    WRITE(*,*) '      conc(1:3) :: ', yvec(1:3)
    WRITE(*,*) 
  END SUBROUTINE
  !
  !
  SUBROUTINE CSR_to_GephiGraph(Matrix,Vnames,FileName)
    USE csv_file
    USE Sparse_Mod, ONLY: CSR_Matrix_T

    TYPE(CSR_Matrix_T) :: Matrix
    CHARACTER(*)       :: FileName
    CHARACTER(*)       :: Vnames(:)

    CHARACTER(80), ALLOCATABLE :: idx(:,:)
    INTEGER :: i, jj, cnt
    
    ALLOCATE(idx(Matrix%nnz,2))

    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName))//'.csv',STATUS='UNKNOWN')
    WRITE(99,'(A13)') 'Source,Target'
    cnt = 0
    DO i = 1 , Matrix%m
      DO jj = Matrix%RowPtr(i) , Matrix%RowPtr(i+1)-1
        cnt = cnt + 1

        ! indizes
        !idx(cnt,1) = i
        !idx(cnt,2) = Matrix%ColInd(jj)
        !CALL csv_write(99,idx(cnt,:),.TRUE.)

        ! strings
        idx(cnt,1) = TRIM(ADJUSTL(Vnames(i)))
        idx(cnt,2) = TRIM(ADJUSTL(Vnames(Matrix%ColInd(jj))))
        CALL csv_write(99,idx(cnt,:),.TRUE.)
       
      END DO
    END DO
    CLOSE(99)
  END SUBROUTINE CSR_to_GephiGraph

  SUBROUTINE  Matrix_Statistics(A,B,BA,BAT,S_HG,Jac,M,LUM)
    
    USE Sparse_Mod, ONLY: CSR_Matrix_T
    TYPE(CSR_Matrix_T) :: A,B,BA,BAT,S_HG,Jac,M,LUM
    297 format(10X,A)
    298 format(10X,A18,3(I12,A2))

    WRITE(*,297) '                 |     rows    |    colums   |      nnz    |'
    WRITE(*,297) ' ----------------+-------------+-------------+-------------+-'
    WRITE(*,298) '           alpha |', A%m,   ' |',A%n,     ' |',A%nnz,   ' |'
    WRITE(*,298) '            beta |', B%m,   ' |',B%n,     ' |',B%nnz,   ' |'
    WRITE(*,298) '  (beta-alpha)^T |', BAT%m, ' |',BAT%n,   ' |',BAT%nnz, ' |'
    !WRITE(*,298) '   Species Graph |', S_HG%m,' |',S_HG%n,  ' |',S_HG%nnz,' |'
    WRITE(*,298) '  Jacobian (= J) |', Jac%m, ' |',Jac%n,   ' |',Jac%nnz, ' |'
    WRITE(*,298) '       I - h*g*J |', M%m,   ' |',M%n,     ' |',M%nnz,   ' |'
    WRITE(*,298) '   LU(I - h*g*J) |', LUM%m, ' |',LUM%n,   ' |',LUM%nnz, ' |'
    WRITE(*,*)
    WRITE(*,*)

  END SUBROUTINE Matrix_Statistics


  !
  SUBROUTINE WriteAnalysisFile(RS,species_names,mixing_ratios,IntRate)
    USE Kind_Mod
    USE Control_Mod
    USE Reac_Mod
    USE ChemSys_Mod, ONLY: ReactionStruct_T

    TYPE(ReactionStruct_T), INTENT(IN) :: RS(:)
    CHARACTER(*),           INTENT(IN) :: species_names(:)
    REAL(dp),               INTENT(IN) :: mixing_ratios(:,:)
    REAL(dp),               INTENT(IN) :: IntRate(:)

    INTEGER  :: i, ieq
    INTEGER  :: lat, alt, day, month
    REAL(dp) :: TimeStp
    
    lat   = 45
    alt   = 20
    day   = 21
    month = 6
    TimeStp = Tspan(2)-Tspan(1)

    OPEN(UNIT=99,FILE='pathway_analysis.txt',STATUS='UNKNOWN')

    ! Headline (Denotation or conditions of model run)
    WRITE( 99 , '(2(I3,A5),2(I0.2,A1))' ) lat , ' N;  ' , alt , ' km;  ' , day , '.' , month , '.'

    ! Timestep = Length of the time interval of interest
    IF (TimeStp<3600.d0) THEN
      WRITE( 99 , '(A11,F5.3,A2)' ) 'Timestep:  ', TimeStp,' s'
    ELSE
      WRITE( 99 , '(A11,F0.3,A2)' ) 'Timestep:  ', TimeStp/3600.0d0,' h'
    END IF
    
    ! Number of chemical species
    WRITE( 99 , '(1X,I0,A)' ) nspc, ' species:  Mixing rations: c_0, c_mean, dc [ppb]:'

    ! The three data columns corresponding to chemical species contain:
    !   - mixing ratio at the beginning of the time interval [ppb]
    !   - mean mixing ratio during the time interval [ppb]
    !   - change of the mixing ratio during the time interval [ppb]
    DO i = 1 , nspc
      WRITE( 99 , '(3(Es14.7e2,3X),4X,A)' ) mixing_ratios(i,:) , TRIM(species_names(i))
    END DO

    ! Number of chemical reactions
    WRITE( 99 , '(I0,A)' ) neq, ' reactions:  Integrated reaction rates [ppb]:'
    
    ! The data column corresponding to reactions contains:
    !   - reaction rate integrated over the time interval of interest [ppb]
    DO i = 1 , neq
      IF ( Combustion) THEN
        ieq = INDEX(TRIM(RS(i)%Line1), ' => ')
        WRITE( 99 , '(Es11.4e2,4X,A,I0,A,A)' ) IntRate(i), TRIM(RS(i)%TypeConstant)//'_',i,':  ', &
        &       TRIM(TRIM(RS(i)%Line1(:ieq))//' -> '//TRIM(ADJUSTL(RS(i)%Line1(ieq+4:))))
      ELSE
        ieq = INDEX(TRIM(RS(i)%Line1), ' = ')
        IF ( INDEX(TRIM(RS(i)%TypeConstant), 'PHOT') > 0 ) THEN
          WRITE( 99 , '(Es11.4e2,4X,A,I0,A,A)' ) IntRate(i), TRIM(RS(i)%TypeConstant)//'_',i,':  ', &
          &       TRIM(TRIM(RS(i)%Line1(:ieq))//'  +  hv  -> '//TRIM(ADJUSTL(RS(i)%Line1(ieq+3:))))
        ELSE
          WRITE( 99 , '(Es11.4e2,4X,A,I0,A,A)' ) IntRate(i), TRIM(RS(i)%TypeConstant)//'_',i,':  ', &
          &       TRIM(TRIM(RS(i)%Line1(:ieq))//' -> '//TRIM(ADJUSTL(RS(i)%Line1(ieq+3:))))
        END IF
      END IF
    END DO
    CLOSE(99)

    WRITE(*,*)
    WRITE(*,*) '  Analysis File written: pathway_analysis.txt' 
    WRITE(*,*)
  END SUBROUTINE WriteAnalysisFile

  SUBROUTINE file_err(filename,io_stat,io_msg)
    CHARACTER(Len=*), INTENT(in) :: filename
    INTEGER         , INTENT(in) :: io_stat
    CHARACTER(Len=*), INTENT(in), OPTIONAL :: io_msg
    IF (io_stat /= 0) THEN
      WRITE(*,"(79('!'))")
      WRITE(*,'(A,I0)')    'ERROR operating on file:  '//TRIM(filename)//'  with io status:  ',io_stat 
      IF (PRESENT(io_msg)) WRITE(*,'(A)')       'Message:  '//TRIM(io_msg)
      WRITE(*,"(79('!'))")
      WRITE(*,*)'Exit ...'
      STOP
    END IF
  END SUBROUTINE file_err

  SUBROUTINE ShowMaxErrorCounter()
    USE Control_Mod, ONLY: maxErrorCounter, BSP, LinAlg, RtolROW, AtolGas, AtolAqua, AtolTemp
    USE Reac_Mod,    ONLY: nspc, nr, y_name
    USE ChemSys_Mod,ONLY: ReactionSystem
    INTEGER :: i

    OPEN(UNIT=99,FILE='OUTPUT/LocalErrorMaxima.log',STATUS='UNKNOWN')
    WRITE(99,*) '         Mechanism: ', TRIM(BSP)
    WRITE(99,*) '    Linear Algebra: ', LinAlg
    WRITE(99,*) '  Tolerance   rel.: ', RtolROW
    WRITE(99,*) '      (gas)   abs.: ', AtolGas
    WRITE(99,*) '      (aqua)  abs.: ', AtolAqua
    WRITE(99,*) '      (temp)  abs.: ', AtolTemp
    WRITE(99,*)
    DO i=1,nspc      
      IF ( maxErrorCounter(i) > 0 ) THEN
        WRITE(99,'(A,2X,I8,5X,A)') ' number of local error maximum = ' , maxErrorCounter(i), TRIM(y_name(i))
      END IF
    END DO
    CLOSE(99)
  END SUBROUTINE ShowMaxErrorCounter

  SUBROUTINE SequentialReadNewReactionList(UnitNr)
    USE Control_Mod
    INTEGER,      INTENT(IN) :: UnitNr
    INTEGER :: i, io_err, nLines
    CHARACTER(LenLine) :: Line

    nLines = 0
    DO 
      READ(UnitNr,'(A)',IOSTAT=io_err) Line
      nLines = nLines + 1
      IF (io_err>0.OR.io_err<0) EXIT
    END DO
    REWIND(UnitNr)
    nLines = nLines - 1
    
    !write(*,*) ' Number of Reactions = ', nLines
    ALLOCATE(newReac_List(nLines))
    DO i=1,nLines
      READ(UnitNr,*,IOSTAT=io_err) newReac_List(i)
      IF (io_err>0.OR.io_err<0) WRITE(*,*)  ' ERROR: ',io_err
    END DO
    CLOSE(UnitNr)
  END SUBROUTINE SequentialReadNewReactionList

  SUBROUTINE SYS_TO_KPP(RS)
    USE Kind_Mod
    USE Control_Mod, ONLY: BSP, Tspan, Temperature0, StpNetcdf
    USE Reac_Mod,    ONLY: y_name, RO2, InitValAct, InitValKat, Diag_Name, iNcdfGas
    USE ChemSys_Mod, ONLY: ReactionStruct_T, Duct_T, PositionSpeciesAll
    USE Meteo_Mod,   ONLY: N2, O2, H2O

    TYPE Dupe_T
      INTEGER :: iReaction, nDuplicates
      INTEGER, ALLOCATABLE :: iDuplicates(:)
    END TYPE Dupe_T

    TYPE(ReactionStruct_T), INTENT(IN) :: RS(:)

    CHARACTER(400) :: Reaction, LeftSide, RightSide, renameSpc, renameSpc2(4)
    CHARACTER(400) :: Reactions1, Reactions2, RateExpression
    
    REAL(dp), ALLOCATABLE :: Param(:)
    INTEGER,  ALLOCATABLE :: Duplicates(:,:)  !, iIndex(:), jIndex(:)
    TYPE(Dupe_T), ALLOCATABLE :: DupeList(:)
    
    INTEGER, PARAMETER :: File_Unit = 112
    INTEGER, PARAMETER :: File_Unit2 = 113
    INTEGER :: iReac, jReac, iSpc, iEnd, iDupe, nSpc, nReac, totDupes=0, cntDupes=1
    INTEGER :: iBrac1, iBrac2, iCol1, iCol2, iSCol1, iSCol
    INTEGER :: io_stat
    INTEGER :: RO2_4,RO2_rest
    CHARACTER(200) :: io_msg

    
    nReac = SIZE(RS)
    ! saves up to five duplicate reactions per reaction
    ! if more are found, increas the value 5 in the allocation
    ALLOCATE(Duplicates(nReac,5))
    Duplicates = -999

    ! --------------------------------------
    ! --- Check for duplicate reactions in .sys file
    DO iReac = 1 , nReac
      
      iDupe = 0

      INNER_LOOP: DO jReac = iReac+1 , nReac

        IF (      RS(iReac)%nActEd   /=      RS(jReac)%nActEd  ) CYCLE INNER_LOOP
        IF (      RS(iReac)%nActPro  /=      RS(jReac)%nActPro ) CYCLE INNER_LOOP
        IF ( TRIM(RS(iReac)%Type)    /= TRIM(RS(jReac)%Type)   ) CYCLE INNER_LOOP

        DO iSpc = 1 , RS(iReac)%nActEd ! does not matter if iReac or jReac because same nActEd
          IF ( RS(iReac)%Educt(iSpc)%iSpecies /= RS(jReac)%Educt(iSpc)%iSpecies ) CYCLE INNER_LOOP
          IF ( RS(iReac)%Educt(iSpc)%Koeff    /= RS(jReac)%Educt(iSpc)%Koeff    ) CYCLE INNER_LOOP
        END DO
        DO iSpc = 1 , RS(iReac)%nActPro ! does not matter if iReac or jReac because same nActPro
          IF ( RS(iReac)%Product(iSpc)%iSpecies /= RS(jReac)%Product(iSpc)%iSpecies ) CYCLE INNER_LOOP
          IF ( RS(iReac)%Product(iSpc)%Koeff    /= RS(jReac)%Product(iSpc)%Koeff    ) CYCLE INNER_LOOP
        END DO

        iDupe = iDupe + 1
        Duplicates(iReac,iDupe) = jReac

      END DO INNER_LOOP
    END DO

    ! turn into list, save only reactions where duplicates were found
    totDupes = COUNT(Duplicates(:,1)/=-999 )

    IF (totDupes == 0) WRITE(*,*)  '         No duplicate reactions found'

    ALLOCATE( DupeList(totDupes) )
    iDupe = 1
    DO iReac = 1 , nReac
      IF ( Duplicates(iReac,1) /= -999 ) THEN
        DupeList(iDupe)%iReaction   = iReac
        DupeList(iDupe)%nDuplicates = COUNT(Duplicates(iReac,:) /= -999)
        DupeList(iDupe)%iDuplicates = [ Duplicates( iReac , 1:DupeList(iDupe)%nDuplicates ) ]
        iDupe = iDupe + 1
      END IF
    END DO
    DEALLOCATE(Duplicates)


    ! --------------------------------------
    ! --- Writing the KPP equation file
    !
    OPEN(UNIT=File_Unit, FILE='OUTPUT/'//TRIM(BSP)//'.eqn', ACTION='write')
    WRITE(File_Unit,'(A)') '#EQUATIONS'

    jReac = 0
    DO iReac=1,nReac
      
      ! Renaming all species names to SPC1,...SPCnspc and construct:
      !  LeftSide  = a_iReac1 * SPC1 + a_iReac2 * SPC2 + ...   
      ! and 
      !  RightSide = b_iReac1 * SPC1 + b_iReac2 * SPC2 + ....

      LeftSide  = BuildReactionString( RS(iReac)%Educt   )
      RightSide = BuildReactionString( RS(iReac)%Product )

      ! concatinate LeftSide and RightSide
      Reaction = TRIM(ADJUSTL(LeftSide))//' = '//TRIM(ADJUSTL(RightSide))

      ! construct reaction rate expression (functions, parameters)
      RateExpression = BuildRateExpression( iReac )

      IF ( iReac_is_not_a_duplicate(iReac) ) THEN
        jReac = jReac + 1
        WRITE(File_Unit,'(A,I0,A)') '{',jReac,'.}  '//TRIM(Reaction)//' :   '//TRIM(RateExpression)//'  ;'      
      END IF

    END DO
    CLOSE(File_Unit)
    
    WRITE(*,*);  WRITE(*,'(10X,A)') 'OUTPUT/'//TRIM(BSP)//'.eqn file written.'



    ! ------------------------------------------------------------
    ! Writing the KPP species file
    !
    OPEN(UNIT=File_Unit, FILE='OUTPUT/'//TRIM(BSP)//'.spc', ACTION='write')

    WRITE(File_Unit,'(A)') '#include atoms'
    WRITE(File_Unit,'(A)')
    WRITE(File_Unit,'(A)') '#DEFVAR'

    nSpc = SIZE(InitValAct)
    DO iSpc=1,nSpc
      WRITE(renameSpc,'(A,I0)') 'SPC',iSpc
      WRITE(File_Unit,'(A)') TRIM(renameSpc)//' = IGNORE ;    {'//TRIM(y_name(iSpc))//'}'
      renameSpc = ''
    END DO
    DO iSpc=1,SIZE(InitValKat)
      WRITE(renameSpc,'(A,I0)') 'SPC',iSpc+nSpc
      WRITE(File_Unit,'(A)') TRIM(renameSpc)//' = IGNORE ;    {'//TRIM(y_name(iSpc+nSpc))//'}'
      renameSpc = ''
    END DO
    
    CLOSE(File_Unit)
    WRITE(*,'(10X,A)') 'OUTPUT/'//TRIM(BSP)//'.spc file written.'


    ! ------------------------------------------------------------
    ! Writing the KPP species file
    !
    OPEN(UNIT=646, FILE='OUTPUT/'//TRIM(BSP)//'_AtCSol_to_KPP.spc', ACTION='write')

    WRITE(646,'(A)') '# This file contains the names of Species in KPP format and AtCSol SMILES notation.'
    WRITE(646,'(A)') '# Serves as input file for the Matlab Program PrintConcentrations.m'
    WRITE(646,'(A)')

    nSpc = SIZE(InitValAct)
    DO iSpc=1,nSpc
      WRITE(renameSpc,'(A,I0)') 'SPC',iSpc
      WRITE(646,'(A)') TRIM(renameSpc)//' = '//TRIM(y_name(iSpc))
      renameSpc = ''
    END DO
    DO iSpc=1,SIZE(InitValKat)
      WRITE(renameSpc,'(A,I0)') 'SPC',iSpc+nSpc
      WRITE(646,'(A)') TRIM(renameSpc)//' = '//TRIM(y_name(iSpc+nSpc))
      renameSpc = ''
    END DO
    
    CLOSE(646)
    WRITE(*,'(10X,A)') 'OUTPUT/'//TRIM(BSP)//'_AtCSol_to_KPP.spc file written.'


    ! ------------------------------------------------------------
    ! Writing the KPP definiton (.def) file
    ! containing initial values, monitoring species, rate functions, etc.
    !
    OPEN(UNIT=File_Unit, FILE='OUTPUT/'//TRIM(BSP)//'.def', ACTION='write')

    ! define the species you want to analyse
    ! be sure to alter the MAX_MONITOR varialbe 
    ! in kpp.2.2.3/src/gen.c if more than 8 species want to be diagnosed
    !
    WRITE(File_Unit,'(A)') '#MONITOR            {Screen Output}'
    WRITE(File_Unit,'(A)') 

    DO iSPc=1,SIZE(iNcdfGas)
      WRITE(renameSpc,'(A,I0)') 'SPC',PositionSpeciesAll(y_name(iNcdfGas(iSpc)))
      WRITE(File_Unit,'(A)') TRIM(renameSpc)//';    {'//TRIM(y_name(iNcdfGas(iSpc)))//'}' 
      renameSpc = ''
    END DO

    WRITE(File_Unit,'(A)')
    WRITE(File_Unit,'(A)') '#INITVALUES             {Initial Values}'   
    WRITE(File_Unit,'(A)') '  CFACTOR = 1.0    ;              {Conversion Factor} '
    WRITE(File_Unit,'(A)')
    !DO iSpc=1,SIZE(InitValAct)
    !  WRITE(renameSpc,'(A,I0)') 'SPC',iSpc
    !  ! write initialvalues for KPP (copy fort.401 into mcm_32_Initialize.f90)
    !  IF (InitValAct(iSpc)>1.0e-8_dp) THEN
    !    WRITE(File_Unit,'(2X,A,D16.8,A)') TRIM(renameSpc)//' = ', InitValAct(iSpc), ';'
    !  END IF
    !  renameSpc = ''
    !END DO
    !WRITE(File_Unit,'(A)')
    !WRITE(File_Unit,'(A)')

    CALL Include_RateFunctions()


    ! define RCONST calculation (reaction rate constant)
    !
    WRITE(File_Unit,'(A)')    '#INLINE F90_RCONST'

    WRITE(File_Unit,'(A)')    '  REAL(dp) :: M0, M, N2, O2, RO2, H2O'
    WRITE(File_Unit,'(A)')    '  ! variables'
    WRITE(File_Unit,'(A)')    '  REAL(dp), PARAMETER     :: PiHalf = 2.0_dp*ATAN(1.0_dp)'
    WRITE(File_Unit,'(A)')    '  REAL(dp), PARAMETER :: Pres = 850.d0 ! hPa'
    WRITE(File_Unit,'(A)')    '  REAL(dp), PARAMETER :: p0   = 1013.25d0 ! hPa Normaldruck'
    WRITE(File_Unit,'(A)')    '  REAL(dp), PARAMETER :: RefTemp = 298.15D0'
  
    WRITE(File_Unit,'(A)')    '  REAL(dp)                ::  chi'
    WRITE(File_Unit,'(A)') 
    WRITE(File_Unit,'(A)') 
    WRITE(File_Unit,'(A,D16.8)')    '  N2 = ', N2
    WRITE(File_Unit,'(A,D16.8)')    '  O2 = ', O2
    WRITE(File_Unit,'(A,D16.8)')    '  M0 = ', N2 + O2
    WRITE(File_Unit,'(A)')    '  M   = M0 * RefTemp / TEMP * Pres / p0'
    WRITE(File_Unit,'(A,D16.8)')    '  H2O = ', H2O
    WRITE(File_Unit,'(A)') 
    WRITE(File_Unit,'(A)') 
    WRITE(File_Unit,'(A)')    '  ! --- Update photo reactions "J(.)" ---'
    WRITE(File_Unit,'(A)')    '  chi = Zenith(TIME)'

    ! ------------------------------------------------
    ! --- Writing the operation for RO2 factor to file
    ! 
    RO2_rest = MODULO(SIZE(RO2),4)
    RO2_4    = SIZE(RO2) - RO2_rest
    
    WRITE(File_Unit,'(2X,A)') 'RO2 = &'

    ! four species next to each other
    DO iSpc = 1 , RO2_4 , 4
      WRITE(File_Unit,140) 'SPC',RO2(iSpc),  'SPC',RO2(iSpc+1),&
                           'SPC',RO2(iSpc+2),'SPC',RO2(iSpc+3)
    END DO
    
    ! one species per line
    DO iSpc = RO2_4+1 , RO2_4+RO2_rest-1
      WRITE(renameSpc,'(A,I0)') 'SPC',RO2(iSpc)
      WRITE(File_Unit,141, IOSTAT=io_stat, IOMSG=io_msg) TRIM(renameSpc)
    END DO

    ! last species without '&'
    WRITE(renameSpc,'(A,I0)') 'SPC',RO2(SIZE(RO2))
    WRITE(File_Unit,142, IOSTAT=io_stat, IOMSG=io_msg) TRIM(renameSpc)
    WRITE(File_Unit,'(A)') '#ENDINLINE'

    !--- RO2 factors
    140 FORMAT(4X,3('C(ind_',A3,I0,') + '),'C(ind_',A3,I0,') + &')
    141 FORMAT(4X,'C(ind_',A,') + & ')
    142 FORMAT(4X,'C(ind_',A,')')


    ! ------------------------------------------------
    ! --- Define star time, end time and timestep (times at which a timstep is saved)
    !
    WRITE(File_Unit,'(A)')
    WRITE(File_Unit,'(A)')
    WRITE(File_Unit,'(A)') '#INLINE F90_INIT'
    WRITE(File_Unit,'(A,D16.8)') '  TSTART = ', Tspan(1)
    WRITE(File_Unit,'(A,D16.8)') '  TEND   = ', Tspan(2)
    WRITE(File_Unit,'(A,D16.8)') '  DT     = ', StpNetcdf
    WRITE(File_Unit,'(A,D16.8)') '  TEMP   = ', Temperature0

    ! --- define initial values here because of better readability
    DO iSpc=1,SIZE(InitValAct)
      WRITE(renameSpc,'(A,I0)') 'SPC',iSpc
      IF (InitValAct(iSpc)>1.0e-8_dp) THEN
        WRITE(File_Unit,'(A,D16.8,A)') "  IF ( TRIM(SPC_NAMES(i)) == '"//TRIM(renameSpc)//"' ) "//&
        &            "  VAR(i) = (",InitValAct(iSpc),"_dp)*CFACTOR      ! "//TRIM(y_name(iSpc))
      END IF
      renameSpc = ''
    END DO

    WRITE(File_Unit,'(A)') '#ENDINLINE'


    CLOSE(File_Unit)
    WRITE(*,'(10X,A)') 'OUTPUT/'//TRIM(BSP)//'.def file written.'


    ! ------------------------------------------------------------
    ! Writing the KPP definiton (.def) file
    ! containing initial values, monitoring species, rate functions, etc.
    !
    OPEN(UNIT=File_Unit, FILE='OUTPUT/'//TRIM(BSP)//'.kpp', ACTION='write')

    WRITE(File_Unit,'(A)')    ' #INTEGRATOR rosenbrock'
    WRITE(File_Unit,'(A)')    ' #LANGUAGE   Fortran90'
    WRITE(File_Unit,'(A)')    ' #DOUBLE     on'
    WRITE(File_Unit,'(A)')    ' #DRIVER     none'
    WRITE(File_Unit,'(A)')    ' #HESSIAN    off'
    WRITE(File_Unit,'(A)')    ' #MEX        off'
    WRITE(File_Unit,'(A)')    ' #STOICMAT   off'
    WRITE(File_Unit,'(A)')    ' #EQNTAGS    on'
    WRITE(File_Unit,'(A)')    ' #INCLUDE '//TRIM(BSP)//'.spc'
    WRITE(File_Unit,'(A)')    ' #INCLUDE '//TRIM(BSP)//'.def'
    WRITE(File_Unit,'(A)')    ' #INCLUDE '//TRIM(BSP)//'.eqn'

    CLOSE(File_Unit)
    WRITE(*,'(10X,A)') 'OUTPUT/'//TRIM(BSP)//'.kpp file written.'
    WRITE(*,*)

    contains

      FUNCTION BuildReactionString(Ducts) RESULT(ReacString)
        TYPE(Duct_T), INTENT(IN) :: Ducts(:)
        CHARACTER(400) :: ReacString
        INTEGER :: nSpc, iSpc, iEnd

        ReacString  = ''

        IF ( Ducts(1)%Koeff /= 1.0_dp ) THEN
          WRITE(ReacString,'(F5.3,A,I0)') Ducts(1)%Koeff,' SPC',Ducts(1)%iSpecies
        ELSE
          WRITE(ReacString,'(A,I0)') 'SPC',Ducts(1)%iSpecies
        END IF
        nSpc = SIZE(Ducts)
        IF ( nSpc >1 ) THEN
          DO iSpc = 2 , nSpc
            iEnd = LEN_TRIM(ReacString)
            IF ( Ducts(iSpc)%iSpecies == 0 ) THEN
              WRITE(ReacString(iEnd+1:),'(" + ",A)') 'SPCdummy'
            ELSE
              IF ( Ducts(iSpc)%Koeff /= 1.0_dp ) THEN
                WRITE(ReacString(iEnd+1:),'(" + ",F5.3,A,I0)') Ducts(iSpc)%Koeff,' SPC',Ducts(iSpc)%iSpecies
              ELSE
                WRITE(ReacString(iEnd+1:),'(" + ",A,I0)') 'SPC',Ducts(iSpc)%iSpecies
              END IF
            END IF
          END DO
        END IF
      END FUNCTION BuildReactionString


      FUNCTION BuildRateExpression(iReacIN) RESULT(RateExpr)
        INTEGER :: iReacIN

        CHARACTER(10)  :: Factor
        CHARACTER(400) :: tmpRateExpr, RateExpr

        INTEGER :: nDupes, i, iR
        INTEGER, ALLOCATABLE :: iDupes(:)

        IF ( DupeList(cntDupes)%iReaction == iReacIN ) THEN
          iDupes = DupeList(cntDupes)%iDuplicates
          nDupes = DupeList(cntDupes)%nDuplicates
          IF ( cntDupes < totDupes ) cntDupes = cntDupes + 1
        ELSE
          nDupes = 0
        END IF

        iR = iReacIN
        RateExpr = ''
        DO i = 0,nDupes

          IF ( i > 0 ) iR = iDupes(i)

          ! check if there are facotrs involved, if so, multiply the rate expression with the factor
          SELECT CASE ( TRIM(RS(iR)%Factor) )
            CASE('$H2');    Factor = '*H2'
            CASE('$O2N2');  Factor = '*O2*N2'
            CASE('$M');     Factor = '*M'
            CASE('$O2');    Factor = '*O2'
            CASE('$N2');    Factor = '*N2'
            CASE('$H2O');   Factor = '*H2O'
            CASE('$O2O2');  Factor = '*O2*O2'  
            CASE('$RO2');   Factor = '*RO2'
            CASE('$RO2aq'); Factor = '*RO2aq'
            CASE('$aH2O');  Factor = '*aH2O'
            CASE DEFAULT;   Factor = ''
          END SELECT

          ! depending on the type of rate constant, write the rate expression string to tmpRateExpr
          SELECT CASE ( TRIM(RS(iR)%TypeConstant) )
            CASE('CONST');    WRITE(tmpRateExpr,100)  RS(iR)%Constants,TRIM(Factor)
            CASE('TEMP0');    WRITE(tmpRateExpr,101)  RS(iR)%Constants,TRIM(Factor)
            CASE('TEMP1');    WRITE(tmpRateExpr,102)  RS(iR)%Constants,TRIM(Factor)
            CASE('TEMP2');    WRITE(tmpRateExpr,103)  RS(iR)%Constants,TRIM(Factor)
            CASE('TEMP3');    WRITE(tmpRateExpr,104)  RS(iR)%Constants,TRIM(Factor)
            CASE('TEMP4');    WRITE(tmpRateExpr,105)  RS(iR)%Constants,TRIM(Factor)
            CASE('SPEC3');    WRITE(tmpRateExpr,130)  RS(iR)%Constants,TRIM(Factor)
            CASE('SPEC2MCM'); WRITE(tmpRateExpr,131)  RS(iR)%Constants,TRIM(Factor)
            CASE('SPEC3MCM'); WRITE(tmpRateExpr,132)  RS(iR)%Constants,TRIM(Factor)
            CASE('SPEC4MCM'); WRITE(tmpRateExpr,133)  RS(iR)%Constants,TRIM(Factor)
            CASE('SPEC5MCM'); WRITE(tmpRateExpr,134)  RS(iR)%Constants,TRIM(Factor)
            CASE('SPEC6MCM'); WRITE(tmpRateExpr,135)  RS(iR)%Constants,TRIM(Factor)
            CASE('SPEC7MCM'); WRITE(tmpRateExpr,136)  RS(iR)%Constants,TRIM(Factor)
            CASE('TROEMCM');  WRITE(tmpRateExpr,110)  RS(iR)%Constants,TRIM(Factor)
            CASE('PHOTMCM');  WRITE(tmpRateExpr,120)  RS(iR)%Constants,TRIM(Factor)
          END SELECT

          ! if there where duplicate reactions, add the rate expressions together
          IF ( LEN_TRIM(RateExpr) > 0 ) THEN
            RateExpr = TRIM(ADJUSTL(RateExpr))//' + '//TRIM(ADJUSTL(tmpRateExpr))
          ELSE
            RateExpr = TRIM(ADJUSTL(tmpRateExpr))
          END IF

        END DO

        !   Definiton of the rate expression format
        !
        !---  constant
        100 FORMAT(D16.8,A)    ! CONST
        !--- temperature dependent arrhenius
        101 FORMAT(D16.8,'*TEMP**',D16.8,'*EXP(-(',F9.2,'/TEMP)',A)       ! TEMP0
        102 FORMAT(D16.8,'*EXP(-1.0D0*(',F9.2,'/TEMP))',A)                ! TEMP1
        103 FORMAT(D16.8,'*TEMP*TEMP*EXP(-1.0D0*(',F9.2,'/TEMP))',A)      ! TEMP2
        104 FORMAT(D16.8,'*EXP(',D16.8,'*(1.0D0/TEMP-1.0D0/298.15D0))',A) ! TEMP3
        105 FORMAT(D16.8,'*TEMP*EXP(-1.0D0*(',F9.2,'/TEMP))',A)           ! TEMP4
        !--- troe pressure dependent reactions
        110 FORMAT('k_TROEMCM(',10(D16.8,','),'TEMP,M)',A)  ! TROEMCM
        !--- mcm version photolysis
        120 FORMAT('k_PHOTOMCM(',3(D16.8,','),'chi)',A)      ! PHOTMCM
        !--- Special Types  (Gas Phase: Density-Dependent)
        130 FORMAT('k_SPEC3(',6(D16.8,','),'TEMP,M)',A)   ! SPEC3
        131 FORMAT('k_SPEC2MCM(',3(D16.8,','),'TEMP)',A)     ! SPEC2MCM
        132 FORMAT('k_SPEC3MCM(',2(D16.8,','),'TEMP,M)',A)   ! SPEC3MCM
        133 FORMAT('k_SPEC4MCM(',4(D16.8,','),'H2O,TEMP)',A) ! SPEC4MCM
        134 FORMAT('k_SPEC5MCM(',4(D16.8,','),'TEMP,M)',A)   ! SPEC5MCM
        135 FORMAT('k_SPEC6MCM(',4(D16.8,','),'TEMP)',A)     ! SPEC6MCM
        136 FORMAT('k_SPEC7MCM(',6(D16.8,','),'TEMP)',A)     ! SPEC7MCM

      END FUNCTION BuildRateExpression


      FUNCTION iReac_is_not_a_duplicate(iR) RESULT(isDupe)
        INTEGER, INTENT(IN) :: iR
        LOGICAL :: isDupe
        INTEGER :: i, j

        isDupe = .TRUE.
        DO i = 1 , totDupes
          DO j = 1 , DupeList(i)%nDuplicates
            IF ( iR == DupeList(i)%iDuplicates(j) ) isDupe = .FALSE.
          END DO
        END DO
      END FUNCTION iReac_is_not_a_duplicate

      SUBROUTINE Include_RateFunctions()
      
        INTEGER :: i, nFcn, io_stat, NstationFiles
        CHARACTER(LEN=400) :: io_msg, Line
        CHARACTER(LEN=100) :: file
        CHARACTER(LEN=100), ALLOCATABLE :: stationFileNames(:)
        CHARACTER(LEN=100), ALLOCATABLE :: RateFcns(:)
      
        ! get the files
        CALL SYSTEM('ls ./KPP_RATEFCN > fileContents.txt')
        OPEN(UNIT=31, FILE='fileContents.txt',ACTION="read")
        ! count how many .fcn file are available
        nFcn = 0
        DO
         READ(31,'(A)',IOSTAT=io_stat) file
         IF (io_stat/=0) EXIT
         file = ADJUSTR(file)
         IF ( file(97:100) == '.fcn' ) nFcn = nFcn + 1
        END DO

        ! save them into a list
        ALLOCATE(RateFcns(nFcn))
        REWIND(31)
        i = 1
        DO 
         READ(31,'(A)',IOSTAT=io_stat) file
         IF (io_stat/=0) EXIT
         file = ADJUSTR(file)
         IF ( file(97:100) == '.fcn' ) THEN
           WRITE(RateFcns(i), '(A)') ADJUSTL(file) 
           i= i + 1
         END IF
        END DO


        WRITE(File_Unit,'(A)') '#INLINE F90_RATES'
        WRITE(File_Unit,'(A)')

        ! copy content of rate functions (KPP_RATEFCN/*.fcn) into .def file

        DO i = 1 , nFcn
          OPEN(UNIT=32, FILE='KPP_RATEFCN/'//TRIM(RateFcns(i)), IOSTAT=io_stat, IOMSG=io_msg)
          IF (io_stat/=0) EXIT
          DO 
            READ(32,'(A)', IOSTAT=io_stat, IOMSG=io_msg ) Line
            IF (io_stat/=0) EXIT
            WRITE(File_Unit,'(A)') TRIM(Line)
          END DO
          WRITE(File_Unit,*)
          CLOSE(32)
        END DO

        WRITE(File_Unit,'(A)') '#ENDINLINE'
        WRITE(File_unit,'(A)')

        
      END SUBROUTINE Include_RateFunctions

  END SUBROUTINE SYS_TO_KPP


  SUBROUTINE OpenFile_wStream(UnitNr,FileName)
    INTEGER,      INTENT(IN) :: UnitNr
    CHARACTER(*), INTENT(IN) :: FileName
    INTEGER :: io_stat
    OPEN(unit=UnitNr, file=FileName, status='replace', action='write', access='stream', iostat=io_stat)
    CALL file_err(FileName,io_stat)
  END SUBROUTINE OpenFile_wStream

  SUBROUTINE OpenFile_wSeq(UnitNr,FileName)
    INTEGER,      INTENT(IN) :: UnitNr
    CHARACTER(*), INTENT(IN) :: FileName
    INTEGER :: io_stat
    OPEN(unit=UnitNr, file=FileName, status='replace', action='write', access='sequential', iostat=io_stat)
    CALL file_err(FileName,io_stat)
  END SUBROUTINE OpenFile_wSeq

  SUBROUTINE OpenFile_rSeq(UnitNr,FileName)
    INTEGER,      INTENT(IN) :: UnitNr
    CHARACTER(*), INTENT(IN) :: FileName
    INTEGER :: io_stat
    OPEN(unit=UnitNr, file=FileName, status='old', action='read', access='sequential', iostat=io_stat)
    CALL file_err(FileName,io_stat)
  END SUBROUTINE OpenFile_rSeq

  SUBROUTINE OpenFile_rStream(UnitNr,FileName)
    INTEGER,      INTENT(IN) :: UnitNr
    CHARACTER(*), INTENT(IN) :: FileName
    INTEGER :: io_stat
    OPEN(unit=UnitNr, file=FileName, status='old', action='read', access='stream', iostat=io_stat)
    CALL file_err(FileName,io_stat)
  END SUBROUTINE OpenFile_rStream

  SUBROUTINE ConvertTime(t,fmt)
    USE Kind_Mod
    REAL(dp),   INTENT(INOUT) :: t    ! in seconds
    CHARACTER(8), INTENT(OUT) :: fmt

    IF ( t > 60.0_dp) THEN
      t   = t/60.0_dp
      fmt = ' [min]'
      IF ( t > 60.0_dp) THEN
        t   = t/60.0_dp
        fmt = ' [hours]'
        IF ( t > 24.0_dp) THEN
          t   = t/24.0_dp
          fmt = ' [days]'
        END IF
      END IF
    ELSE
      fmt = ' [sec]'
    END IF
    
  END SUBROUTINE ConvertTime
END MODULE IO_Mod

