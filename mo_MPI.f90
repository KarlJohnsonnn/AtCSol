MODULE mo_MPI
  !
  USE MPI
  USE Kind_Mod
  USE mo_control
  USE mo_reac, ONLY: LowerRateLim, UpperRateLim
  IMPLICIT NONE
  INTEGER :: MPI_dp
  !
  INTEGER :: MPI_ID
  INTEGER :: MPI_np
  INTEGER :: MPIErr
  LOGICAL :: MPI_master
  !
  INTEGER, ALLOCATABLE :: MyParties(:,:)
  !
  CONTAINS
  !
  !
  ! MPI Initialization
  SUBROUTINE StartMPI()
    CALL MPI_INIT( MPIErr )
    CALL CheckMPIErr() 
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, MPI_ID, MPIErr )
    CALL CheckMPIErr() 
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, MPI_np, MPIErr )
    CALL CheckMPIErr() 
    MPI_dp=MPI_Double_Precision
    IF ( MPI_id == 0 ) THEN
      MPI_master = .TRUE.
    ELSE
      MPI_master = .FALSE.
    END IF
  END SUBROUTINE StartMPI
  !
  !
  ! MPI Finalize
  SUBROUTINE FinishMPI()
    CALL MPI_FINALIZE( MPIErr )
    CALL CheckMPIErr() 
  END SUBROUTINE FinishMPI
  !
  !
  SUBROUTINE CheckMPIErr()
    IF (MPIErr/=MPI_SUCCESS) then
      WRITE(*,*) 'MPI error, ID = ',MPI_ID
      WRITE(*,*) 'MPIerr = ', MPIErr
      CALL MPI_FINALIZE(MPIErr)
      STOP 'MPI  error'
    END IF
  END SUBROUTINE CheckMPIErr
  !
  !
  SUBROUTINE dropout
    CALL FinishMPI()
    STOP 'Stoppp!"'
  END SUBROUTINE dropout
  !
  !
  ! This routine will build a matrix (#Process x 3) to
  ! store the integer values:
  !      Split(i,1) ... part of Reaction System left bound  [idxl]
  !      Split(i,2) ... part of Reaction System right bound [idxr]
  !      Split(i,3) ... number of elements to calculate in the i-th part
  SUBROUTINE BuildPartitions(Parties,nreac,nproc)
    INTEGER :: nreac, nproc
    INTEGER, ALLOCATABLE :: Parties(:,:)
    !
    INTEGER :: spltbnd
    INTEGER :: i
    !
    !
    ALLOCATE(Parties(nproc,3))
    Parties=0
    !
    LowerRateLim=1
    IF (nproc>1) THEN
      spltbnd=nreac/nproc
      Parties(1,2)=spltbnd
      Parties(1,3)=spltbnd
      LowerRateLim=MPI_ID*spltbnd+1
      UpperRateLim=(MPI_ID+1)*spltbnd
      !
      DO i=2,nproc-1
        Parties(i,1)=Parties(i-1,2)
        Parties(i,2)=Parties(i,1)+spltbnd
        Parties(i,3)=spltbnd
      END DO
      Parties(nproc,1)=Parties(nproc-1,2)
    END IF
    Parties(nproc,2)=nreac
    Parties(nproc,3)=nreac-Parties(nproc,1)
    IF (MPI_ID==nproc-1) UpperRateLim=nreac
    !WRITE(*,*) 'lowerLimit=',lowerratelim,'UpperLimit=',upperratelim
  END SUBROUTINE BuildPartitions
  !
  !
  SUBROUTINE GatherAllPartitions(Vec,Partitions)
    REAL(dp) :: Vec(:)
    INTEGER :: Partitions(:,:)
    !
    ! collect the parts of Rate on all processes
    ! -- IN_PLACE --> no need for sendcount and sendtype
    ! -- ":" = proc0,..,procN-1
    CALL MPI_AllGatherV( MPI_IN_PLACE        &  ! sendbuffer 
    &                  , 0                   &  ! sendcount 
    &                  , MPI_DATATYPE_NULL   &  ! sendtype
    &                  , Vec                 &  ! recvbuffer
    &                  , Partitions(:,3)     &  ! recvcounts    
    &                  , Partitions(:,1)     &  ! displacements  
    &                  , MPI_DOUBLE          &  ! recvtype
    &                  , MPI_COMM_WORLD      &  ! communicator
    &                  , MPIErr              )  ! error code-\
    CALL CheckMPIErr() ! <------------------check------------/
  END SUBROUTINE GatherAllPartitions
  !
  !
  ! This routine will gather the maximum time value of all processes
  SUBROUTINE GetMaxTimes(maxTime,inTime)
    REAL(dp) :: inTime
    REAL(dp) :: maxTime
    !
    !
    CALL MPI_Reduce( inTime          &    ! sendbuffer
    &              , maxTime         &    ! recvbuffer
    &              , 1               &    ! # elements in sendbuffer
    &              , MPI_DOUBLE      &    ! datatype
    &              , MPI_MAX         &    ! operation
    &              , 0               &    ! root process
    &              , MPI_COMM_WORLD  &    ! communicator
    &              , MPIErr          )    ! error code-\
    CALL CheckMPIErr() ! <------------------check------/
  END SUBROUTINE GetMaxTimes
  !
  !
  ! This routine will broadcaste a vector from root process to all other processes
  SUBROUTINE BroadcastIVector(Vec)
    INTEGER :: Vec(:)
    ! 
    INTEGER :: lenVec
    !
    lenVec=SIZE(Vec)
    CALL MPI_Bcast( Vec               &   ! broadcasting vec to all other processes
    &             , lenVec            &   ! number of elements in vec
    &             , MPI_INT           &   ! datatype
    &             , 0                 &   ! send from root
    &             , MPI_COMM_WORLD    &   ! communicator
    &             , MPIErr            )   ! error code-\
    CALL CheckMPIErr() ! <------------------check------/
  END SUBROUTINE BroadcastIVector
  !
  !
  ! Send array of real values from i  to process 0
  SUBROUTINE SendReal(SendV,ResV,i)
    REAL(dp) :: SendV(:)
    REAL(dp) :: ResV(SIZE(SendV))
    INTEGER :: i,nSendV
     
    nSendV=SIZE(SendV)
    IF (MPI_ID==i) THEN
      CALL MPI_SEND(SendV, nSendV, MPI_DOUBLE_PRECISION, 0, 999, MPI_COMM_WORLD, MPIErr)
    ELSE IF (MPI_ID==0) THEN
      CALL MPI_RECV(ResV, nSendV, MPI_DOUBLE_PRECISION, i, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE, MPIErr)
    END IF
    CALL CheckMPIErr()
  END SUBROUTINE SendReal

  FUNCTION GatherAquaFractions(sbuf) RESULT(rbuf)
    USE mo_reac,    ONLY: nFrac
    USE mo_control, ONLY: nNcdfAqua
    REAL(dp), INTENT(IN) :: sbuf(nNcdfAqua)    
    REAL(dp)             :: rbuf(nFrac*nNcdfAqua)

    CALL MPI_GATHER( sbuf            &  ! sendbuffer 
    &              , nNcdfAqua       &  ! sendcount 
    &              , MPI_DOUBLE      &  ! sendtype
    &              , rbuf            &  ! recvbuffer
    &              , nNcdfAqua       &  ! receivecount 
    &              , MPI_DOUBLE      &  ! recvtype
    &              , 0               &  ! send to root
    &              , MPI_COMM_WORLD  &  ! communicator
    &              , MPIErr          )  ! error code-
    CALL CheckMPIErr

  END FUNCTION GatherAquaFractions

END MODULE mo_MPI
