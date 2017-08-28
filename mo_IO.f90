MODULE mo_IO
  USE Kind_Mod
  USE mo_MPI
  USE mo_control
  USE Rosenbrock_Mod
  USE Factorisation_Mod
  !
  !
  !
  CONTAINS
  SUBROUTINE Logo()
    IF (MPI_ID==0) THEN
      WRITE(*,*) ''
      WRITE(*,*) '===================================================================================='
      WRITE(*,*) '|      ___           ___           ___           ___                       ___     |'
      WRITE(*,*) '|     /\  \         /\__\         /\  \         /\__\          ___        /\  \    |'
      WRITE(*,*) '|    /::\  \       /:/  /        /::\  \       /::|  |        /\  \      /::\  \   |'
      WRITE(*,*) '|   /:/\:\  \     /:/__/        /:/\:\  \     /:|:|  |        \:\  \    /:/\:\  \  |'
      WRITE(*,*) '|  /:/  \:\  \   /::\  \ ___   /::\~\:\  \   /:/|:|__|__      /::\__\  /::\~\:\  \ |'
      WRITE(*,*) '| /:/__/ \:\__\ /:/\:\  /\__\ /:/\:\ \:\__\ /:/ |::::\__\  __/:/\/__/ /:/\:\ \:\__\|'
      WRITE(*,*) '| \:\  \  \/__/ \/__\:\/:/  / \:\~\:\ \/__/ \/__/~~/:/  / /\/:/  /    \:\~\:\ \/__/|'
      WRITE(*,*) '|  \:\  \            \::/  /   \:\ \:\__\         /:/  /  \::/__/      \:\ \:\__\  |'
      WRITE(*,*) '|   \:\  \           /:/  /     \:\ \/__/        /:/  /    \:\__\       \:\ \/__/  |'
      WRITE(*,*) '|    \:\__\         /:/  /       \:\__\         /:/  /      \/__/        \:\__\    |'
      WRITE(*,*) '|     \/__/         \/__/         \/__/         \/__/                     \/__/    |'
      WRITE(*,*) '===================================================================================='
      WRITE(*,*) ''
    END IF
  END SUBROUTINE Logo
  !
  !
  SUBROUTINE Print_Run_Param()
    IF (MPI_ID==0) THEN
      WRITE(*,*)   ''
      WRITE(*,*)   '  Run - Paramter:'
      WRITE(*,*)   ''
      WRITE(*,*)   '     Mechanism:             ', TRIM(SysFile)
      IF (NetCdfPrint) THEN
        WRITE(*,*)   '     NetCDF-File:           ', 'NetCDF/'//TRIM(BSP)//'.nc'
      ELSE
        WRITE(*,*)   '     NetCDF-File:           ', '*** no NetCDF output ***'
      END IF
      WRITE(*,*)   '     Initials:              ', InitFile
      WRITE(*,*)   '     ODE solver:            ', ODEsolver
      IF (ODEsolver/='LSODE') THEN
        IF (solveLA=='cl') THEN 
          WRITE(*,*)   '     Linear Algebra:        Classic'
        ELSE
          WRITE(*,*)   '     Linear Algebra:        Extended'
        END IF
        IF (Error_Est==2) THEN
          WRITE(*,*)   '     Error Estimation:      Euklid Norm'
        ELSE
          WRITE(*,*)   '     Error Estimation:      Maximum Norm'
        END IF
        IF (OrderingStrategie==8) THEN
          WRITE(*,*)   '     Solve Linear Systems:  Sparse LU, Markowitz Ordering Algorithm'
        ELSE 
          WRITE(*,*)   '     Solve Linear Systems:  MUMPS, Ordering Stragegie:  ',OrderingStrategie 
        END IF
      END IF
      WRITE(*,*)   ''
      IF(ImpEuler/=1) THEN
        WRITE(*,*)   '  Tolerance:   '
        WRITE(*,*)   ''
        WRITE(*,'(A34,2X,Es8.2)')   '      Relative Rosenbrock        = ',RtolROW
        WRITE(*,'(A34,2X,Es8.2)')   '      Absolute (gaseous species) = ',AtolGas
        IF (ns_AQUA>0) WRITE(*,'(A34,2X,Es8.2)')   '      Absolute (aqueous species) = ',AtolAqua
        IF ( Teq ) THEN
          WRITE(*,'(A34,2X,Es8.2)')   '      Absolute Temperature       = ',AtolTemp
        END IF
      END IF
      WRITE(*,*)   ''
    END IF
  END SUBROUTINE Print_Run_Param
  !
  !
  SUBROUTINE Output_Statistics(TRead,TSymb,TFac,TSolve,TRates,TJac,TInte,TAll,TSend,TNcdf, TErr, TRhs)
    !
    REAL(dp) :: TRead,TSymb,TFac,TSolve,TRates,TJac,TInte,TAll,TSend, TNcdf, TErr, TRhs
    REAL(dp) :: maxTRead,maxTSymb,maxTFac,maxTSolve,maxTRates,maxTJac &
    &               , maxTInte,maxTAll,maxTSend,maxtNcdf,maxTErr,maxTRhs
    !
    CALL GetMaxTimes(maxTRead,TRead)
    CALL GetMaxTimes(maxTRates,TRates)
    CALL GetMaxTimes(maxTSymb,TSymb)
    CALL GetMaxTimes(maxTInte,TInte)
    CALL GetMaxTimes(maxTFac,TFac)
    CALL GetMaxTimes(maxTSolve,TSolve)
    CALL GetMaxTimes(maxTSend,TSend)
    CALL GetMaxTimes(maxTJac,TJac)
    CALL GetMaxTimes(maxTAll,TAll)
    CALL GetMaxTimes(maxTNcdf,TNcdf)
    CALL GetMaxTimes(maxTErr,TErr)
    CALL GetMaxTimes(maxTRhs,TRhs)
    !
    IF (MPI_ID==0) THEN
      ! print the statistics
      299 format(A40,3X,F13.6,A6)
      298 format(A40,3X,I10)
      WRITE(*,*)   ''
      WRITE(*,*)   '                                    DONE' 
      WRITE(*,*)   ''
      WRITE(*,*)   ''
      WRITE(*,*)   ''
      WRITE(*,*)   '==========================================================================='
      WRITE(*,*)   '================================ Statistics ==============================='
      WRITE(*,*)   '==========================================================================='
      WRITE(*,*)   ''
      WRITE(*,298) ' Number of successful time steps   =', Output%nsteps
      WRITE(*,298) ' Number of failed time step        =', Output%nfailed
      WRITE(*,298) ' Number of rate evaluations        =', Output%nRateEvals
      WRITE(*,298) ' Number of Jacobian calculations   =', Output%npds
      WRITE(*,298) ' Number of LU factorisations       =', Output%ndecomps
      WRITE(*,298) ' Number of solved linear systems   =', Output%nsolves
      WRITE(*,*)   ''
      WRITE(*,*)   ''
      WRITE(*,*)   ' ========================================================================='
      WRITE(*,*)   ' ======================= max times of all processes ====================== '
      WRITE(*,*)   ' ========================================================================= '
      WRITE(*,*)   ''
      WRITE(*,299) ' Time to read the input system     =', maxTRead,' [sec]'
      WRITE(*,299) ' Time for symbolic phase           =', maxTSymb,' [sec]'
      WRITE(*,299) ' Time writing NetCDF-File          =', maxTNcdf,' [sec]'
      WRITE(*,299) ' Time for integration              =', maxTInte-maxTNcdf,' [sec]'
      WRITE(*,*)   ' -------------------------------------+-----------------------------------'
      WRITE(*,299) '          - Factorisation          =', maxTFac  ,' [sec]'
      WRITE(*,299) '          - Right hand side calc   =', maxTRhs  ,' [sec]'
      WRITE(*,299) '          - Solve linear Systems   =', maxTSolve,' [sec]'
      WRITE(*,299) '          - Rates                  =', maxTRates,' [sec]'
      IF(MPI_np>1) WRITE(*,299) '          - Ratessend              =', maxTSend ,' [sec]'
      WRITE(*,299) '          - Jacobian               =', maxTJac  ,' [sec]'
      WRITE(*,299) '          - Error calculation      =', maxTErr  ,' [sec]'
      WRITE(*,*)   ' = = = = = = = = = = = = = = = = = ======= = = = = = = = = = = = = = = = ='
      WRITE(*,299) ' Total runtime                     =', maxTAll,' [sec]'
      WRITE(*,*)  '===========================================================================' 
      WRITE(*,*)  ''
      WRITE(*,*)  ''
    END IF
  END SUBROUTINE
  !
  !
  SUBROUTINE PrintHeadSimul(ChemFile) 
    CHARACTER(*) ::ChemFile
    INTEGER :: Unit=90
    !---- save data:
    !   - evaluated times => t Vector
    !   - evaluated concentrations => y vector
    !   - calculated step sizes => h vector
    !   - more?
    !
    IF (MPI_Id==0) THEN
      !
      OPEN(UNIT=90,FILE=ADJUSTL(TRIM(ChemFile))//'.simul',STATUS='UNKNOWN')
      WRITE(Unit,*) ' ============================================================'
      WRITE(Unit,*) ' ====== Simulation of chemical systems / TESTVERSION ========'
      WRITE(Unit,*) ' ========     Output -  Chemical Reaction Data      ========='
      WRITE(Unit,*) ' ============================================================'
      WRITE(Unit,*) '  '
      WRITE(Unit,*) '  step        time          stepsize          concentrations   '
      WRITE(Unit,*) '=====================================================================>'
      WRITE(Unit,*) ' '
    END IF
  END SUBROUTINE PrintHeadSimul
  !
  !
  SUBROUTINE SaveTimeStp(Unit,nsteps,tnew,h,y)
    INTEGER :: Unit
    INTEGER :: nsteps
    REAL(dp) :: tnew, h
    REAL(dp) :: y(:)
    IF (MPI_ID==0) THEN
      ! Output Label
      !199 format(I6,3X,E18.12,4X,E18.12,4X,2(E16.6,1X) )
      WRITE(Unit,*) nsteps,tnew,h,y(:)
    END IF
  END SUBROUTINE SaveTimeStp
  !
  !
  SUBROUTINE SaveMatricies(aMat,bMat,cMat,dMat,eMat,fName)
    TYPE(CSR_Matrix_T)   :: aMat,bMat,cMat,dMat,eMat
    CHARACTER(*)         :: fName
    !
    ! local stuff
    INTEGER, ALLOCATABLE :: InvPermu(:)
    INTEGER              :: mUnit, i
    CHARACTER(50)        :: mName
    !
    !
    IF (MPI_ID==0) THEN
      ! only if MatrixPrint=True
      CALL WriteSparseMatrix(aMat,TRIM('matrixOut/alpha'//fName))
      CALL WriteSparseMatrix(bMat,TRIM('matrixOut/beta'//fName))
      CALL WriteSparseMatrix(cMat,TRIM('matrixOut/_beta-alpha_T'//fName))
      CALL WriteSparseMatrix(dMat,TRIM('matrixOut/Miter0'//fName))
      !
      IF (OrderingStrategie<8) THEN
        !
        ! ordering<8 --> MUMPS auto choice ordering
        ! Print structure of LU matrix and Permutation vector
        !
        CALL WriteSparseMatrix(eMat,TRIM('matrixOut/LUmiterStructure'//fName))
        ALLOCATE(InvPermu(eMat%m))
        CALL PermuToInvPer(InvPermu,Mumps_Par%SYM_PERM)
        CALL PrintPerm(Mumps_Par%SYM_PERM,InvPermu,TRIM('matrixOut/Permu'//fName))
        !
        ! Print a Mapping supported by MUMPS
        !
        IF (Mumps_Par%ICNTL(18)==1) THEN
          mUnit=55+MPI_ID
          !
          WRITE(mName,'(I1)') mUnit-55  ! convert integer to char to get a uniq filename
          mName=ADJUSTL('matrixOut/'//fName//'Mapping_'//mName//'.SparseMat')
          OPEN(UNIT=mUnit,FILE=TRIM(mName),STATUS='UNKNOWN')
          DO i=1,Mumps_Par%NZ_loc
            WRITE(mUnit,'(1X,I5,1X,I5,10X,E23.14)')     Mumps_Par%IRN_loc(i)     &
            &                                         , Mumps_Par%JCN_loc(i)     &
            &                                         , 1.0d0
          END DO
          CLOSE(mUnit)
        END IF
      ELSE
        !
        ! ordering>=8 --> Markowitz count (early minimum degree)
        ! Print structure of LU matrix and Permutation vector
        !
        CALL WriteSparseMatrix(eMat,TRIM('matrixOut/LUmiterStructure'//fName))
        CALL PrintPerm(eMat%Permu,eMat%InvPer,TRIM('matrixOut/Permu'//fName))
      END IF
      !
      CALL FinishMPI()
      STOP 'MatrixPrint=True --> stop nach print in Integration_Mod'
    END IF
  END SUBROUTINE SaveMatricies
  !
  !
  SUBROUTINE CloseSimulFile(Unit)
    INTEGER :: Unit
    !
    CLOSE(90)
  END SUBROUTINE CloseSimulFile
  !
  !
  SUBROUTINE DebugPrint1(yvec,rvec,stepsize,time)
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
    
    TYPE(CSR_Matrix_T) :: A,B,BA,BAT,S_HG,Jac,M,LUM

    298 format(A18,3(I12,A2))

    WRITE(*,*)   '                |     rows    |    colums   |      nnz    |'
    WRITE(*,*)   '  --------------+-------------+-------------+-------------+-'
    WRITE(*,298) '           alpha |', A%m,   ' |',A%n,     ' |',A%nnz,   ' |'
    WRITE(*,298) '            beta |', B%m,   ' |',B%n,     ' |',B%nnz,   ' |'
    WRITE(*,298) '  (beta-alpha)^T |', BAT%m, ' |',BAT%n,   ' |',BAT%nnz, ' |'
    WRITE(*,298) '   Species Graph |', S_HG%m,' |',S_HG%n,  ' |',S_HG%nnz,' |'
    WRITE(*,298) '  Jacobian (= J) |', Jac%m, ' |',Jac%n,   ' |',Jac%nnz, ' |'
    WRITE(*,298) '       I - h*g*J |', M%m,   ' |',M%n,     ' |',M%nnz,   ' |'
    WRITE(*,298) '   LU(I - h*g*J) |', LUM%m, ' |',LUM%n,   ' |',LUM%nnz, ' |'
    WRITE(*,*)

  END SUBROUTINE Matrix_Statistics
  !
  !
  SUBROUTINE WriteAnalysisFile(RS,species_names,mixing_ratios,IntRate)
    USE mo_reac, ONLY: neq, nspc

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
END MODULE mo_IO

