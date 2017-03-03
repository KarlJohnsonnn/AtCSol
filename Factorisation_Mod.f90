MODULE Factorisation_Mod

! ZEITEN FÜR SOLVE UND FACT MITTELN UND AUSGEBEN LASSEN 
  USE Sparse_Mod
  USE Kind_Mod
  USE mo_control
  USE mo_MPI
  USE mo_reac
  IMPLICIT NONE

  TYPE(DMUMPS_STRUC) :: Mumps_Par
  REAL(RealKind), PRIVATE :: timerEnd,timerStart
  INTEGER, ALLOCATABLE :: MAPPING(:)
  INTEGER, ALLOCATABLE :: loc_diagPtr(:)
  INTEGER, ALLOCATABLE :: glob_diagPtr(:)
  CONTAINS
  
SUBROUTINE InitMumps(A,givenPermutaion)

  TYPE(SpRowIndColInd_T), INTENT(IN) :: A
  INTEGER, INTENT(IN), OPTIONAL :: givenPermutaion(:)
  INTEGER :: i, loc_NZcnt , loc_rateCnt
  
  timerStart=MPI_WTIME()
  ! Init phase
  Mumps_Par%COMM=MPI_COMM_SELF
  Mumps_Par%SYM=0  ! unsymmetric
  !Mumps_Par%SYM=1  ! assumed to be symmetric positive definite
  Mumps_Par%PAR=1   ! host is also involved in the parallel steps of the factorization and solve phases

  !
  ! initializes an instance of the package
  Mumps_Par%JOB=-1  
  Mumps_Par%ICNTL(1)=4  !4
  Mumps_Par%ICNTL(2)=4  !4 
  CALL DMUMPS(Mumps_Par)
  
  !
  ! set parameters for analysis phase
  Mumps_Par%N=A%m
  Mumps_Par%NZ=SIZE(A%RowInd)
  Mumps_Par%IRN=>A%RowInd
  Mumps_Par%JCN=>A%ColInd
  Mumps_Par%A=>A%Val
  ALLOCATE(Mumps_Par%Rhs(Mumps_Par%N))
  Mumps_Par%Rhs=0.0d0
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  !----->  ICNTL(.) <------ 
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  !   
  Mumps_Par%ICNTL(1)=0  !4
  Mumps_Par%ICNTL(2)=0  !4 
  !Mumps_Par%ICNTL(3)=0  !4
  !Mumps_Par%ICNTL(4)=2  !4
  Mumps_Par%ICNTL(5)=0                          ! assembled matrix format
  Mumps_Par%ICNTL(6)=0                           
  Mumps_Par%ICNTL(8)=0                         ! auto scaling strategie
  !Mumps_Par%ICNTL(11)=1                         ! calc statistic related to an error analysis of the linear system solved
  Mumps_Par%ICNTL(14)=50
  !
  ! Analysis phase of MUMPS
  Mumps_Par%JOB=1
  Mumps_Par%ICNTL(7)=OrderingStrategie          ! in RUN/*.run
  !CALL DMUMPS(Mumps_Par)
  !
  IF (MPI_np>1) THEN
    ! SYMBOLISCHE PHASE SERIELL + MAPPING
    !Mumps_Par%ICNTL(18)=1  !  MUMPS returns a mapping 
    Mumps_Par%ICNTL(28)=1   !  serial computation of the ordering is performed
    !
    IF (ParOrdering>=0) THEN
      ! SYMBOLISCHE PHASE PARALLEL + MAPPING
      Mumps_Par%ICNTL(28)=2  !  parallel ordering is performed, in this case ICNTL(7) is meaningless
      Mumps_Par%ICNTL(29)=ParOrdering ! 1 = PT-Scotch , 2 = ParMetis
      IF (ParOrdering>2) Mumps_Par%ICNTL(29)=0  !  auto select PT-SCOTCH(=1) or ParMetis(=2) to reorder the input matrix
    END IF
  ELSE
    ! SYMBOLISCHE PHASE SERIELL
    Mumps_Par%ICNTL(18)=0   !  MUMPS returns no mapping 
    Mumps_Par%ICNTL(28)=1   !  serial computation of the ordering is performed
  END IF
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  ! ----->  CNTL(.) <------ 
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  ! 
  !Mumps_Par%CNTL(1)=0.0d0      ! relative threshold for numerical pivoting
  Mumps_Par%CNTL(1)=0.1d0      ! relative threshold for numerical pivoting
  Mumps_Par%CNTL(4)=1.0d-4      ! determines the threshold for static pivoting
  Mumps_Par%CNTL(5)=1.0d20      ! defines the fixation for null pivots
  CALL DMUMPS(Mumps_Par)
  !
  Mumps_Par%ICNTL(1)=0  !4
  Mumps_Par%ICNTL(2)=0  !4 
  Mumps_Par%ICNTL(3)=0  !4
  Mumps_Par%ICNTL(4)=0  !4
  ! Allocate local parts of the matrix and distribute the entries
  !              ( only used if MPI_np > 1 )
  !IF (Mumps_Par%ICNTL(18)==1) THEN
  !  ALLOCATE(MAPPING(Mumps_Par%NZ))
  !  ALLOCATE(glob_diagPtr(Mumps_Par%N))
  !  MAPPING(:)=-1
  !  glob_diagPtr(:)=-1
  !  !
  !  ! broadcast MAPPING vector from proc0 to proc1,...,proc(n-1)
  !  IF (MPI_ID==0) MAPPING(:)=Mumps_Par%MAPPING(:)
  !  CALL BroadcastIVector(MAPPING)
  !  !
  !  ! allocate local matricies
  !  Mumps_Par%NZ_loc=COUNT(MAPPING==MPI_ID)
  !  ALLOCATE(Mumps_Par%IRN_loc(Mumps_Par%NZ_loc))
  !  ALLOCATE(Mumps_Par%JCN_loc(Mumps_Par%NZ_loc))
  !  ALLOCATE(Mumps_Par%A_loc(Mumps_Par%NZ_loc))
  !  !
  !  loc_NZcnt=0
  !  loc_diagCnt=0
  !  DO i=1,Mumps_Par%NZ
  !    !
  !    ! == LOCAL MEMORY DECLARATIONS HERE ==
  !    IF (MAPPING(i)==MPI_ID) THEN
  !      loc_NZcnt=loc_NZcnt+1
  !      Mumps_Par%IRN_loc(loc_NZcnt)=A%RowInd(i)
  !      Mumps_Par%JCN_loc(loc_NZcnt)=A%ColInd(i)
  !      Mumps_Par%A_loc(loc_NZcnt)=A%Val(i)
  !      !
  !      ! search the diagonal entries and count them
  !      IF (A%RowInd(i)==A%ColInd(i)) THEN
  !        loc_diagCnt=loc_diagCnt+1
  !        IF (A%ColInd(i)<=neq) THEN
  !          loc_rateCnt=loc_rateCnt+1
  !        ELSE
  !          loc_concCnt=loc_concCnt+1
  !        END IF
  !      END IF
  !    END IF
  !  END DO
  !  !
  !  ALLOCATE(loc_diagPtr(loc_diagCnt))
  !  ALLOCATE(loc_ratePtr(loc_rateCnt))
  !  ALLOCATE(loc_concPtr(loc_concCnt))
  !  loc_diagPtr(:)=-1
  !  loc_ratePtr(:)=-1
  !  loc_concPtr(:)=-1
  !  loc_diagCnt=0
  !  loc_rateCnt=0
  !  loc_concCnt=0
  !  loc_NZcnt=0
  !  DO i=1,Mumps_Par%NZ
  !    !
  !    ! == LOCAL MEMORY DECLARATIONS HERE ==
  !    IF (MAPPING(i)==MPI_ID) THEN
  !      loc_NZcnt=loc_NZcnt+1
  !      !
  !      ! search the diagonal entries and store their position
  !      IF (A%RowInd(i)==A%ColInd(i)) THEN
  !        loc_diagCnt=loc_diagCnt+1
  !        loc_diagPtr(loc_diagCnt)=loc_NZcnt
  !        glob_diagPtr(loc_diagCnt)=i
  !        IF (A%ColInd(i)<=neq) THEN
  !          loc_rateCnt=loc_rateCnt+1
  !          loc_ratePtr(loc_rateCnt)=A%ColInd(i)
  !        ELSE
  !          loc_concCnt=loc_concCnt+1
  !          loc_concPtr(loc_concCnt)=A%ColInd(i)
  !        END IF
  !      END IF
  !    END IF
  !  END DO
  !  !
  !  ! deallocate global global row,col index arrays
  !  !DEALLOCATE(Mumps_Par%IRN)
  !  !DEALLOCATE(Mumps_Par%JCN)
  !ELSE
  !  !
  !  ! Allocation for serial use
  !  ALLOCATE(glob_diagPtr(Mumps_Par%N))
  !  glob_diagPtr(:)=-1
  !  !
  !  loc_diagCnt=0
  !  DO i=1,Mumps_Par%NZ
  !    IF(A%RowInd(i)==A%ColInd(i)) THEN
  !      loc_diagCnt=loc_diagCnt+1
  !      glob_diagPtr(loc_diagCnt)=i
  !    END IF
  !  END DO
  !  LowerRateLim=1
  !  UpperRateLim=neq
  !END IF
  !
  ! print time for analyse phase
  timerEnd=MPI_WTime()
  IF (Mumps_Par%MyId==0) THEN
    WRITE(*,'(A25,3X,F12.6)') 'Time DMUMPS Analyse',timerEnd-timerStart
  END IF
  !CALL dropout   
END SUBROUTINE InitMumps

SUBROUTINE FactorizeCoefMat(DiagRate,DiagYConc)
  REAL(RealKind), INTENT(IN) :: DiagRate(:)
  REAL(RealKind), INTENT(IN), OPTIONAL :: DiagYConc(:)
  !
  INTEGER :: i, idx_R, idx_D, idx_C
  !
  Mumps_Par%JOB=2
  !
  ! extended matrix version
  !IF (solveLA=='ex') THEN
  !  IF (Mumps_Par%ICNTL(18)==1) THEN
  !    DO i=1,loc_rateCnt
  !     ! print*,MPI_ID, 'i= ',i,'glob_diagPtr(i)',glob_diagPtr(i),'loc_ratePtr(i)',loc_ratePtr(i),'Val=',DiagRate(loc_ratePtr(i))
  !      idx_D=loc_diagPtr(i)
  !      idx_R=loc_ratePtr(i)
  !      Mumps_Par%A_loc(idx_D)=DiagRate(idx_R)
  !    END DO
  !    DO i=1,loc_concCnt
  !      !print*,MPI_ID, 'i= ',i,'glob_diagPtr(i)',glob_diagPtr(i+neq),'loc_concPtr(i)',loc_concPtr(i),'Val=',DiagYConc(loc_concPtr(i)-neq)
  !      idx_D=loc_diagPtr(i+loc_rateCnt)
  !      idx_C=loc_concPtr(i)-neq
  !      Mumps_Par%A_loc(idx_D)=DiagYConc(idx_C)
  !    END DO
  !  ELSE
  !    DO i=1,neq
  !      Mumps_Par%A(glob_diagPtr(i))=DiagRate(i)
  !    !print*, 'i= ',i,'glob_diagPtr(i)',glob_diagPtr(i),'Val=',DiagRate(i)
  !    END DO
  !    DO i=1,nspc
  !      !print*, 'i= ',i,'glob_diagPtr(i)',glob_diagPtr(i+neq),'Val=',DiagYConc(i)
  !      Mumps_Par%A(glob_diagPtr(i+neq))=DiagYConc(i)
  !    END DO
  !  END IF
  !ELSE
  !
  ! classic matrix version
    Mumps_Par%A=DiagRate(:)     ! ganze matrix wird hier übergeben
  !END IF
  timerStart=MPI_WTime()
  NumFac=NumFac+1
  CALL DMUMPS(Mumps_Par)
  timerEnd=MPI_WTime()
  TimeFac=TimeFac+(timerEnd-timerStart)
  !IF (Mumps_Par%MyId==0) THEN
  !  WRITE(*,'(A25,3X,F12.6)') 'Time DMUMPS Factorize',timefact
  !:END IF
END SUBROUTINE FactorizeCoefMat

SUBROUTINE SolveLinAlg(Rhs)
  REAL(RealKind) :: Rhs(:)

  Mumps_Par%JOB=3
  Mumps_Par%Rhs(:)=Rhs(:)
  Mumps_Par%ICNTL(14)=20
  timerStart=MPI_WTime()
  CALL DMUMPS(Mumps_Par)
  timerEnd=MPI_WTime()
  TimeSolve=TimeSolve+(timerEnd-timerStart)
  !IF (Mumps_Par%MyId==0) THEN
  !  WRITE(*,*) 'Time DMUMPS Solve',tEnd-tStart
  !END IF
END SUBROUTINE SolveLinAlg
END MODULE Factorisation_Mod