MODULE IO_Mod
  IMPLICIT NONE
  !
  CONTAINS
  SUBROUTINE Logo()
    USE MPI_Mod
    IF (MPI_master) THEN
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
    END IF
    777 FORMAT(10X,A)
  END SUBROUTINE Logo

  SUBROUTINE Logo2()
    USE MPI_Mod
    IF (MPI_master) THEN
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
    END IF
    777 FORMAT(10X,A)
  END SUBROUTINE Logo2
  !
  !
  SUBROUTINE Print_Run_Param()
    USE MPI_Mod
    USE Control_Mod
    USE Reac_Mod

    IF ( INDEX(SysFile,'.sys')==0)  SysFile = TRIM(SysFile)//'.sys'

    IF (MPI_master.AND.Simulation) THEN
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
      IF ( Teq ) THEN
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
    CALL GetMaxTimes(maxTRead,Time_Read)
    CALL GetMaxTimes(maxTRates,TimeRates)
    CALL GetMaxTimes(maxTSymb,TimeSymbolic)
    CALL GetMaxTimes(maxTInte,TimeIntegration)
    CALL GetMaxTimes(maxTFac,TimeFac)
    CALL GetMaxTimes(maxTSolve,TimeSolve)
    CALL GetMaxTimes(maxTJac,TimeJac)
    CALL GetMaxTimes(maxTAll,Timer_Finish)
    CALL GetMaxTimes(maxTNcdf,TimeNetCDF)
    CALL GetMaxTimes(maxTErr,TimeErrCalc)
    CALL GetMaxTimes(maxTRhs,TimeRhsCalc)
    CALL GetMaxTimes(maxTFlux,TimeFluxWrite)
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
    IF (MPI_master) THEN
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
    END IF
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
    IF (MPI_master) THEN
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
      CALL FinishMPI()
      STOP 'MatrixPrint=True --> stop nach print in Integration_Mod'
    END IF
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
      IF ( ChemKin) THEN
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
    USE Control_Mod, ONLY: BSP
    USE Reac_Mod,    ONLY: y_name
    USE ChemSys_Mod, ONLY: ReactionStruct_T  

    TYPE(ReactionStruct_T), INTENT(IN) :: RS(:)

    CHARACTER(400) :: Reaction, LeftSide, RightSide, renameSpc
    CHARACTER(400) :: Reactions1, Reactions2
    CHARACTER(10)  :: Type
    REAL(dp), ALLOCATABLE :: Param(:)
    
    INTEGER, PARAMETER :: File_Unit = 112
    INTEGER, PARAMETER :: File_Unit2 = 113
    INTEGER :: iReac, iSpc, iEnd, nSpc
    INTEGER :: iBrac1, iBrac2, iCol1, iCol2, iSCol1, iSCol2
    INTEGER :: io_stat
    CHARACTER(200) :: io_msg

    !OPEN(UNIT=File_Unit, FILE='OUTPUT/'//TRIM(BSP)//'.eqn', ACTION='write')
    OPEN(UNIT=File_Unit, FILE='OUTPUT/mcm_32.eqn_pre', ACTION='write')
    WRITE(File_Unit,'(A)') '#EQUATIONS'


    DO iReac=1,SIZE(RS)
      
      ! Renaming all species names to SPC1,...SPCnspc
      
      LeftSide  = ''
      IF ( RS(iReac)%Educt(1)%Koeff /= 1.0_dp ) THEN
        WRITE(LeftSide,'(F5.3,A,I0)') RS(iReac)%Educt(1)%Koeff,' SPC',RS(iReac)%Educt(1)%iSpecies
      ELSE
        WRITE(LeftSide,'(A,I0)') 'SPC',RS(iReac)%Educt(1)%iSpecies
      END IF
      nSpc = RS(iReac)%nActEd
      IF ( nSpc >1 ) THEN
        DO iSpc = 2 , nSpc
          iEnd = LEN_TRIM(LeftSide)
          IF ( RS(iReac)%Educt(iSpc)%iSpecies == 0 ) THEN
            WRITE(LeftSide(iEnd+1:),'(" + ",A)') 'SPCdummy'
          ELSE
            IF ( RS(iReac)%Educt(iSpc)%Koeff /= 1.0_dp ) THEN
              WRITE(LeftSide(iEnd+1:),'(" + ",F5.3,A,I0)') RS(iReac)%Educt(iSpc)%Koeff,' SPC',RS(iReac)%Educt(iSpc)%iSpecies
            ELSE
              WRITE(LeftSide(iEnd+1:),'(" + ",A,I0)') 'SPC',RS(iReac)%Educt(iSpc)%iSpecies
            END IF
          END IF
        END DO
      END IF

      RightSide = ''
      IF ( RS(iReac)%Product(1)%Koeff /= 1.0_dp ) THEN
        WRITE(RightSide,'(F5.3,A,I0)') RS(iReac)%Product(1)%Koeff,' SPC',RS(iReac)%Product(1)%iSpecies
      ELSE
        WRITE(RightSide,'(A,I0)') 'SPC',RS(iReac)%Product(1)%iSpecies
      END IF
      nSpc = RS(iReac)%nActPro
      IF ( nSpc >1 ) THEN
        DO iSpc = 2 , nSpc
          iEnd = LEN_TRIM(RightSide)
          IF ( RS(iReac)%Product(iSpc)%iSpecies == 0 ) THEN
            WRITE(RightSide(iEnd+1:),'(" + ",A)') 'SPCdummy'
          ELSE
            IF ( RS(iReac)%Product(iSpc)%Koeff /= 1.0_dp ) THEN
              WRITE(RightSide(iEnd+1:),'(" + ",F5.3,A,I0)') RS(iReac)%Product(iSpc)%Koeff,' SPC',RS(iReac)%Product(iSpc)%iSpecies
            ELSE
              WRITE(RightSide(iEnd+1:),'(" + ",A,I0)') 'SPC',RS(iReac)%Product(iSpc)%iSpecies
            END IF
          END IF
        END DO
      END IF

      Reaction = TRIM(ADJUSTL(LeftSide))//' = '//TRIM(ADJUSTL(RightSide))
      Type     = ADJUSTL(RS(iReac)%TypeConstant)
      Param    = [RS(iReac)%Constants]

      SELECT CASE( TRIM(Type) )
        CASE('CONST'); WRITE(File_Unit,100) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('TEMP0'); WRITE(File_Unit,101) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('TEMP1'); WRITE(File_Unit,102) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('TEMP2'); WRITE(File_Unit,103) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('TEMP3'); WRITE(File_Unit,104) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('TEMP4'); WRITE(File_Unit,105) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('SPEC3'); WRITE(File_Unit,130) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('SPEC2MCM'); WRITE(File_Unit,131) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('SPEC3MCM'); WRITE(File_Unit,132) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('SPEC4MCM'); WRITE(File_Unit,133) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('SPEC5MCM'); WRITE(File_Unit,134) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('SPEC6MCM'); WRITE(File_Unit,135) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('SPEC7MCM'); WRITE(File_Unit,136) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('TROEMCM'); WRITE(File_Unit,110) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
        CASE('PHOTMCM'); WRITE(File_Unit,120) '{',iReac,'.}  '//TRIM(Reaction)//' :   ',Param
      END SELECT

      DEALLOCATE(Param)

    END DO
    CLOSE(File_Unit)


    ! check for duplicate reactions
    OPEN(UNIT=File_Unit, FILE='OUTPUT/mcm_32.eqn_pre', STATUS='OLD', ACTION='READ')
    REWIND(File_Unit)
    READ(File_Unit,'(A400)',IOSTAT=io_stat,IOMSG=io_msg) Reaction     ! dummy (#EQUATIONS line)

    ! open new file (final form)
    OPEN(UNIT=File_Unit2, FILE='OUTPUT/mcm_32.eqn', ACTION='WRITE')
    WRITE(File_Unit2,'(A)') '#EQUATIONS'



    iReac = 0

    DO 
      READ(File_Unit,'(A400)',IOSTAT=io_stat,IOMSG=io_msg) Reactions1
      READ(File_Unit,'(A400)',IOSTAT=io_stat,IOMSG=io_msg) Reactions2

      iBrac1 = INDEX(Reactions1,'}'); iCol1 = INDEX(Reactions1,':')
      iBrac2 = INDEX(Reactions2,'}'); iCol2 = INDEX(Reactions2,':')


      iEnd = INDEX( ADJUSTL(TRIM(Reactions1(iBrac1+1:iCol1-1))) &
      &           , ADJUSTL(TRIM(Reactions2(iBrac2+1:iCol2-1))) ) 

      ! same reaction
      IF ( iEnd > 0 ) THEN
        iSCol1 = INDEX(Reactions1,';')
        iSCol2 = INDEX(Reactions2,';')

        Reaction = TRIM(ADJUSTL(Reactions1(iBrac1+1:iCol1-1)))//' :   '// &
        &          TRIM(ADJUSTL(Reactions1(iCol1+1:iSCol1-1)))//' + '//&
        &          TRIM(ADJUSTL(Reactions2(iCol2+1:)))

        iReac = iReac + 1
        WRITE(File_Unit2,'("{",I0,".}  ",A)') iReac, TRIM(ADJUSTL(Reaction))

      ELSE
        iReac = iReac + 1
        WRITE(File_Unit2,'("{",I0,".}  ",A)') iReac, TRIM(ADJUSTL(Reactions1(iBrac1+1:)))
        BACKSPACE(File_Unit)
      END IF

      IF ( io_stat < 0 ) EXIT

    END DO

    CLOSE(File_Unit)
    CLOSE(File_Unit2)

    ! check for duplicate reactions
    OPEN(UNIT=File_Unit, FILE='OUTPUT/mcm_32.eqn', STATUS='OLD', ACTION='READ')
    REWIND(File_Unit)
    READ(File_Unit,'(A400)',IOSTAT=io_stat,IOMSG=io_msg) Reaction     ! dummy (#EQUATIONS line)

    ! open new file (final form)
    OPEN(UNIT=File_Unit2, FILE='OUTPUT/mcm_32_final.eqn', ACTION='WRITE')
    WRITE(File_Unit2,'(A)') '#EQUATIONS'

    iReac = 0

    DO 
      READ(File_Unit,'(A400)',IOSTAT=io_stat,IOMSG=io_msg) Reactions1
      READ(File_Unit,'(A400)',IOSTAT=io_stat,IOMSG=io_msg) Reactions2

      iBrac1 = INDEX(Reactions1,'}'); iCol1 = INDEX(Reactions1,':')
      iBrac2 = INDEX(Reactions2,'}'); iCol2 = INDEX(Reactions2,':')


      iEnd = INDEX( ADJUSTL(TRIM(Reactions1(iBrac1+1:iCol1-1))) &
      &           , ADJUSTL(TRIM(Reactions2(iBrac2+1:iCol2-1))) ) 

      ! same reaction
      IF ( iEnd > 0 ) THEN
        iSCol1 = INDEX(Reactions1,';')
        iSCol2 = INDEX(Reactions2,';')

        Reaction = TRIM(ADJUSTL(Reactions1(iBrac1+1:iCol1-1)))//' :   '// &
        &          TRIM(ADJUSTL(Reactions1(iCol1+1:iSCol1-1)))//' + '//&
        &          TRIM(ADJUSTL(Reactions2(iCol2+1:)))

        iReac = iReac + 1
        WRITE(File_Unit2,'("{",I0,".}  ",A)') iReac, TRIM(ADJUSTL(Reaction))

      ELSE
        iReac = iReac + 1
        WRITE(File_Unit2,'("{",I0,".}  ",A)') iReac, TRIM(ADJUSTL(Reactions1(iBrac1+1:)))
        BACKSPACE(File_Unit)
      END IF


      IF ( io_stat < 0 ) EXIT

    END DO

    CLOSE(File_Unit)
    CLOSE(File_Unit2)
    
    
    WRITE(*,*)
    WRITE(*,'(10X,A)') 'OUTPUT/mcm_32.eqn file written.'


    !OPEN(UNIT=File_Unit, FILE='OUTPUT/'//TRIM(BSP)//'.spc', ACTION='write')

    OPEN(UNIT=File_Unit, FILE='OUTPUT/mcm_32.spc', ACTION='write')

    WRITE(File_Unit,'(A)') '#include atoms'
    WRITE(File_Unit,'(A)')
    WRITE(File_Unit,'(A)') '#DEFVAR'
    
    DO iSpc=1,SIZE(y_name)
      WRITE(renameSpc,'(A,I0)') 'SPC',iSpc
      WRITE(File_Unit,'(A)') TRIM(renameSpc)//' = IGNORE ;    {'//TRIM(y_name(iSpc))//'}'
      renameSpc = ''
    END DO
    
    CLOSE(File_Unit)

    WRITE(*,'(10X,A)') 'OUTPUT/mcm_32.spc file written.'
    WRITE(*,*)



    !---  constant
    100 FORMAT(A,I0,A,Es9.2,' ;')    ! CONST
    !--- temperature dependent arrhenius
    101 FORMAT(A,I0,A,Es8.2,'*TEMP**',Es9.2,'*EXP(-(',F7.1,'/TEMP)   ;')       ! TEMP0
    102 FORMAT(A,I0,A,Es8.2,'*EXP(-1.0D0*(',F7.1,'/TEMP))   ;')                ! TEMP1
    103 FORMAT(A,I0,A,Es8.2,'*TEMP*TEMP*EXP(-1.0D0*(',F7.1,'/TEMP))   ;')      ! TEMP2
    104 FORMAT(A,I0,A,Es8.2,'*EXP(',Es9.2,'*(1.0D0/TEMP-1.0D0/298.15D0))   ;') ! TEMP3
    105 FORMAT(A,I0,A,Es8.2,'*TEMP*EXP(-1.0D0*(',F7.1,'/TEMP))   ;')           ! TEMP4
    !--- troe pressure dependent reactions
    110 FORMAT(A,I0,A,'k_TROEMCM(',10(Es9.2,','),'TEMP,M)   ;')  ! TROEMCM
    !--- mcm version photolysis
    120 FORMAT(A,I0,A,'k_PHOTOMCM(',3(Es9.2,','),'chi)   ;')      ! PHOTMCM
    !--- Special Types  (Gas Phase: Density-Dependent)
    130 FORMAT(A,I0,A,   'k_SPEC3(',6(Es9.2,','),'TEMP,M)   ;')   ! SPEC3
    131 FORMAT(A,I0,A,'k_SPEC2MCM(',3(Es9.2,','),'TEMP)   ;')     ! SPEC2MCM
    132 FORMAT(A,I0,A,'k_SPEC3MCM(',2(Es9.2,','),'TEMP,M)   ;')   ! SPEC2MCM
    133 FORMAT(A,I0,A,'k_SPEC4MCM(',4(Es9.2,','),'H2O,TEMP)   ;') ! SPEC2MCM
    134 FORMAT(A,I0,A,'k_SPEC5MCM(',6(Es9.2,','),'TEMP,M)   ;')   ! SPEC2MCM
    135 FORMAT(A,I0,A,'k_SPEC6MCM(',6(Es9.2,','),'TEMP)   ;')     ! SPEC2MCM
    136 FORMAT(A,I0,A,'k_SPEC7MCM(',6(Es9.2,','),'TEMP)   ;')     ! SPEC2MCM

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

