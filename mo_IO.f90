MODULE mo_IO
  IMPLICIT NONE
  !
  CONTAINS
  SUBROUTINE Logo()
    USE mo_MPI
    IF (MPI_master) THEN
      WRITE(*,777) 
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
  !
  !
  SUBROUTINE Print_Run_Param()
    USE mo_MPI
    USE mo_control
    USE mo_reac

    IF ( INDEX(SysFile,'.sys')==0)  SysFile = TRIM(SysFile)//'.sys'

    IF (MPI_master) THEN
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
    USE mo_MPI
    USE mo_control
    !
    REAL(dp) :: maxTRead,maxTSymb,maxTFac,maxTSolve,maxTRates,maxTJac &
    &               , maxTInte,maxTAll,maxTSend,maxtNcdf,maxTErr,maxTRhs
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
      WRITE(*,299) '    reading mechanism       =', maxTRead,' [sec]'
      WRITE(*,299) '    symbolic phase          =', maxTSymb,' [sec]'
      WRITE(*,299) '    writing NetCDF-File     =', maxTNcdf,' [sec]'; WRITE(*,*)
      !WRITE(*,777) '    ------------------------+----------------------------'
      WRITE(*,299) '            factorisation   =', maxTFac  ,' [sec]'
      WRITE(*,299) '          + right-hand side =', maxTRhs  ,' [sec]'
      WRITE(*,299) '          + linear systems  =', maxTSolve,' [sec]'
      WRITE(*,299) '          + reaction rates  =', maxTRates,' [sec]'
      WRITE(*,299) '          + Jacobian        =', maxTJac  ,' [sec]'
      WRITE(*,299) '          + error calc      =', maxTErr  ,' [sec]'
      WRITE(*,777) '    ------------------------=----------------------'
      WRITE(*,299) '    integration             =', maxTInte-maxTNcdf,' [sec]'; WRITE(*,*)
      WRITE(*,299) '    total runtime           =', maxTAll,' [sec]'
      WRITE(*,*);  WRITE(*,*);  WRITE(*,*)
    END IF
  END SUBROUTINE
  
 
  SUBROUTINE SaveMatricies(aMat,bMat,cMat,dMat,eMat,fName)
    USE mo_MPI
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
    USE Kind_MOd
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
  !
  SUBROUTINE WriteAnalysisFile(RS,species_names,mixing_ratios,IntRate)
    USE Kind_Mod
    USE mo_control
    USE mo_reac
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

!  SUBROUTINE StreamWriteFluxes(Rate,t,h)
!    USE Kind_Mod
!    USE mo_control, ONLY: FluxUnit, FluxFile, FluxMetaUnit, FluxMetaFile, iStpFlux
!    REAL(dp) :: Rate(:)
!    REAL(dp) :: t , h
!
!    INTEGER :: io_stat, io_pos
!    CHARACTER(100) :: io_msg
!
!    OPEN(unit=FluxUnit,      file=FluxFile,  status='old',   action='write', &
!    &    position='append', access='stream', iostat=io_stat, iomsg=io_msg    )
!    CALL file_err(FluxFile,io_stat,io_msg)
!    INQUIRE(FluxUnit, POS=io_pos)
!    WRITE(FluxUnit) Rate,t,h
!    CLOSE(FluxUnit)
!
!    iStpFlux   = iStpFlux + 1
!    OPEN(unit=FluxMetaUnit, file=FluxMetaFile, status='old', action='write', position='append')
!    WRITE(FluxMetaUnit,*) iStpFlux, io_pos 
!    CLOSE(FluxMetaUnit)
!  END SUBROUTINE StreamWriteFluxes


  SUBROUTINE SequentialWriteFluxes(Rate,t,h)
    USE Kind_Mod
    USE mo_control, ONLY: FluxUnit, FluxFile, FluxMetaUnit, FluxMetaFile, iStpFlux
    REAL(dp) :: Rate(:)
    REAL(dp) :: t , h

    INTEGER :: io_stat, io_pos
    CHARACTER(100) :: io_msg

    OPEN(unit=FluxUnit,      file=FluxFile,  status='old',   action='write', &
    &    position='append', access='sequential', iostat=io_stat, iomsg=io_msg )
    CALL file_err(FluxFile,io_stat,io_msg)
    WRITE(FluxUnit,'(*(1X,E16.10))') Rate, t, h
    CLOSE(FluxUnit)

  END SUBROUTINE SequentialWriteFluxes


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

  SUBROUTINE SequentialReadNewReactionList(UnitNr)
    USE mo_control
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
END MODULE mo_IO

