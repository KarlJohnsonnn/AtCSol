MODULE issa
  
  USE Kind_Mod
  USE Sparse_Mod, ONLY: CSR_Matrix_T

  IMPLICIT NONE



  CONTAINS  


!*********************************************************************************************
!
! 8888888888 8888888 888      8888888888                  8888888        d88P  .d88888b.  
! 888          888   888      888                           888         d88P  d88P" "Y88b 
! 888          888   888      888                           888        d88P   888     888 
! 8888888      888   888      8888888                       888       d88P    888     888 
! 888          888   888      888                           888      d88P     888     888 
! 888          888   888      888             888888        888     d88P      888     888 
! 888          888   888      888                           888    d88P       Y88b. .d88P 
! 888        8888888 88888888 8888888888                  8888888 d88P         "Y88888P" 
!
!*********************************************************************************************

  FUNCTION Read_Target_Spc(FileName) RESULT(Idx)
    USE ChemSys_Mod, ONLY: PositionSpeciesAll
    USE mo_control,  ONLY: LenLine
    USE mo_reac,     ONLY: y_name
    ! This routine will read the target species for ISSA reduction algorithm.
    ! OUT:
    INTEGER, ALLOCATABLE :: Idx(:)
    ! IN:
    CHARACTER(*)  :: FileName
    ! TEMP:
    INTEGER       :: iPos, j
    CHARACTER(LenLine) :: Line
     
    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName)),STATUS='UNKNOWN')
    REWIND(99)

    ALLOCATE(Idx(0))
    iPos = FindSection(99,'TARGET_SPC')
    IF (iPos==-1) RETURN ! no target species declared
    DO
      READ(99,'(A)') Line;  Line = ADJUSTL(Line)
      IF (TRIM(Line) == 'END_TARGET_SPC') EXIT
      IF (TRIM(Line) == '' ) CYCLE
      iPos = PositionSpeciesAll(Line)
      IF ( iPos > 0 ) Idx = [ Idx , iPos ]
    END DO
    CLOSE(99)
    
    WRITE(*,777) REPEAT('*',39)
    WRITE(*,777) '********** Important Species **********'
    WRITE(*,777) REPEAT('*',39)
    WRITE(*,*)
    DO j=1,SIZE(Idx)
      WRITE(*,'(10X,A,I2,A,I6,5X,A)') '    - S_imp(',j,') = ',Idx(j),TRIM(y_name(Idx(j)))
    END DO
    WRITE(*,*); WRITE(*,*)

    
    777 FORMAT(10X,A)
  END FUNCTION Read_Target_Spc


  FUNCTION Read_Spc_Families(FileName) RESULT(Fam)
    USE ChemSys_Mod, ONLY: PositionSpeciesAll
    USE mo_control,  ONLY: Families_T, LenLine
    ! This routine will read the species families for ISSA reduction algorithm.
    ! OUT: globl type(Families_T) declared at start of this file
    TYPE(Families_T), ALLOCATABLE :: Fam(:)
    ! IN:
    CHARACTER(*)       :: FileName
    ! TEMP:
    CHARACTER(LenLine) :: Line, locLine
    INTEGER            :: i, j, locSp, nFam
     
    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName)),STATUS='UNKNOWN')
    REWIND(99)

    ! count number of superspecies
    i = FindSection(99,'SPC_FAMILIES')
    IF (i==-1) RETURN ! no species families declared
    nFam = FindSection(99,'END_SPC_FAMILIES')
    REWIND(99)

    ALLOCATE(Fam(nFam))
    i = FindSection(99,'SPC_FAMILIES')
    i = 0

    DO 
      READ(99,'(A)') Line;  Line = ADJUSTL(Line)
      IF (TRIM(Line) == 'END_SPC_FAMILIES') EXIT
      IF (TRIM(Line) == '' ) CYCLE
      i = i + 1

      ! count single species for i-th family
      nFam = 0
      locLine = Line
      DO 
        IF ( TRIM(locLine) == '' ) EXIT
        nFam = nFam + 1
        locLine = ADJUSTL(locLine(INDEX(locLine,' ')+1:))
      END DO
      ALLOCATE( Fam(i)%Name(nFam), Fam(i)%Index(nFam) )
      Fam(i)%Index = -11
      
      ! save all single species for i-th family
      locLine = ADJUSTL(Line)
      DO j=1,nFam
        locSp = INDEX(locLine,' ')
        Fam(i)%Name(j)  = TRIM(locLine(:locSp-1))
        Fam(i)%Index(j) = PositionSpeciesAll(Fam(i)%Name(j))
        locLine = ADJUSTL(locLine(locSp+1:))
      END DO
    END DO

    WRITE(*,777) REPEAT('*',39)
    WRITE(*,777) '*********** Species Familes ***********'
    WRITE(*,777) REPEAT('*',39)
    WRITE(*,*)
    DO j=1,SIZE(Fam)
      WRITE(*,'(10X,A,I2,A,*(A))')      '    - Family(',j,')  -->  ',    &
      &    ( TRIM(Fam(j)%Name(i))//' , ' , i=1,SIZE(Fam(j)%Name)-1 ),&
      &      TRIM(Fam(j)%Name(SIZE(Fam(j)%Name)))
    END DO
    WRITE(*,*); WRITE(*,*)
    777 FORMAT(10X,A)
  END FUNCTION Read_Spc_Families

  FUNCTION FindSection(iUnit,SectionName) RESULT(iLine)
    USE mo_control, ONLY: lenLine
    ! OUT:
    INTEGER      :: iLine
    ! IN: 
    INTEGER      :: iUnit
    CHARACTER(*) :: SectionName
    ! TEMP: 
    CHARACTER(LenLine) :: Line
    INTEGER      :: io_err
    
    iLine = 0
    DO
      READ(iUnit,'(A)',IOSTAT=io_err) Line
      IF (io_err==0) THEN
        IF (TRIM(ADJUSTL(Line))==SectionName)  EXIT
        iLine = iLine + 1
      ELSEIF (io_err>0.OR.io_err<0) THEN
        WRITE(*,777) '    ****  No '//TRIM(SectionName)//' found.  ****'
        iLine = -1
        EXIT
      END IF
    END DO
     777 FORMAT(10X,A)
  END FUNCTION FindSection


  ! writing reaction rates , time and stepsize h to file via stream access
  SUBROUTINE StreamWriteFluxes(Rate,t,h)
    USE Kind_Mod
    USE mo_control, ONLY: FluxUnit, FluxFile, FluxMetaUnit, FluxMetaFile, iStpFlux
    REAL(dp) :: Rate(:)
    REAL(dp) :: t , h

    INTEGER :: io_stat, io_pos
    CHARACTER(100) :: io_msg

    OPEN(unit=FluxUnit,      file=FluxFile,  status='old',   action='write', &
    &    position='append', access='stream', iostat=io_stat, iomsg=io_msg    )
    CALL file_err(FluxFile,io_stat,io_msg)
    INQUIRE(FluxUnit, POS=io_pos)
    WRITE(FluxUnit) Rate,t,h
    CLOSE(FluxUnit)

    iStpFlux   = iStpFlux + 1
    OPEN(unit=FluxMetaUnit, file=FluxMetaFile, status='old', action='write', position='append')
    WRITE(FluxMetaUnit,*) iStpFlux, io_pos 
    CLOSE(FluxMetaUnit)

    CONTAINS

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
  END SUBROUTINE StreamWriteFluxes

  SUBROUTINE WriteReaction(Name,iReac,Mech,Class,Param)
    CHARACTER(*) :: Name, Mech, Class, Param
    INTEGER      :: iReac
    INTEGER, PARAMETER :: Reac_Unit = 112

    IF ( iReac == 1 ) THEN
      OPEN(unit=Reac_Unit, file='Reactions_'//Mech//'.txt', action='write')
    ELSE
      OPEN(unit=Reac_Unit, file='Reactions_'//Mech//'.txt', status='old', action='write', position='append')
    END IF

    WRITE(Reac_Unit,'(A)') Class // ' ::: ' // Name // ' ::: ' // Param
    CLOSE(Reac_Unit)

  END SUBROUTINE WriteReaction


!******************************************************************************************************************
!
!  .d8888b. 88888888888 8888888b.  888     888  .d8888b. 88888888888 888     888 8888888b.  8888888888 
! d88P  Y88b    888     888   Y88b 888     888 d88P  Y88b    888     888     888 888   Y88b 888        
! Y88b.         888     888    888 888     888 888    888    888     888     888 888    888 888        
!  "Y888b.      888     888   d88P 888     888 888           888     888     888 888   d88P 8888888    
!     "Y88b.    888     8888888P"  888     888 888           888     888     888 8888888P"  888        
!       "888    888     888 T88b   888     888 888    888    888     888     888 888 T88b   888        
! Y88b  d88P    888     888  T88b  Y88b. .d88P Y88b  d88P    888     Y88b. .d88P 888  T88b  888        
!  "Y8888P"     888     888   T88b  "Y88888P"   "Y8888P"     888      "Y88888P"  888   T88b 8888888888  
!
!******************************************************************************************************************

  SUBROUTINE Get_Rk(R_k,A_T)
    USE mo_control, ONLY: List, nCycles_red, Cyclic_Set_red
    USE Sparse_Mod, ONLY: B, WriteSparseMatrix
    USE mo_unirnk

    TYPE(List), ALLOCATABLE :: R_k(:)
    TYPE(CSR_Matrix_T) :: A_T

    INTEGER, ALLOCATABLE :: ReacSet(:), Perm(:)
    INTEGER :: iSpc
    INTEGER :: newLen
    INTEGER :: i, jj, j , iC, iR, iS, iSpcE, iSpcP, iR_a, iS_b

    ALLOCATE(R_k(nCycles_red+1))
 
    DO iC = 1,nCycles_red

      ALLOCATE(R_k(iC)%List(0))

      DO iS = 1 , Cyclic_Set_red(iC)%len-1
        iSpcE = Cyclic_Set_red(iC)%List(iS)
        iSpcP = Cyclic_Set_red(iC)%List(iS+1)

        ALLOCATE(ReacSet(0))
        DO jj = A_T%RowPtr(iSpcE),A_T%RowPtr(iSpcE+1)-1
          iR_a = A_T%ColInd(jj)
          DO iS_b = B%RowPtr(iR_a),B%RowPtr(iR_a+1)-1
            IF ( B%ColInd(iS_b) == iSpcP ) THEN
              ReacSet = [ ReacSet, iR_a]
            END IF
          END DO
        END DO
        
        R_k(iC)%List = [R_k(iC)%List , ReacSet]
        DEALLOCATE(ReacSet)

      END DO

      ALLOCATE(Perm(SIZE(R_k(iC)%List)))
      CALL unirnk( R_k(iC)%List , Perm , newLen )
      R_k(iC)%len   = newLen
      R_k(iC)%List = [ R_k(iC)%List(Perm(1:newLen)) ]
      DEALLOCATE(Perm)

      !WRITE(*,'(A,I0,A,*(I0,3X))') ' R_k(',iC,')%List = ', R_k(iC)%List
    END DO
   !stop 'issa_MOD'


  END SUBROUTINE Get_Rk


  SUBROUTINE ISSA_structure(R_k,nue,A_T,cycles,RS)
    USE mo_reac,    ONLY: nr, nspc, y_name
    USE mo_control, ONLY: List, nCycles_red
    USE ChemSys_Mod,ONLY: ReactionStruct_T
    USE Sparse_Mod

    TYPE(List), ALLOCATABLE, INTENT(OUT) :: R_k(:)
    TYPE(CSR_Matrix_T)     , INTENT(IN) :: nue, A_T
    TYPE(List)             :: cycles(:)
    TYPE(ReactionStruct_T) :: RS(:)

    INTEGER,    ALLOCATABLE :: tnue_pos(:),tnue_neg(:)
    INTEGER,    ALLOCATABLE :: tnue_pos0(:),tnue_neg0(:), Perm(:)
    INTEGER,    ALLOCATABLE :: NonCyclicRemainder(:)
    INTEGER                 :: nNonCycRem=0


    INTEGER :: i, jj, j , k, kk, iR, iS, iSpc
    INTEGER :: cnt_pos, cnt_neg, nS, N_pos0, N_neg0

    nCycles_red = SIZE(cycles)
    
    CALL Get_Rk(R_k,A_T)

    ! this routine generates the nue_plus and nue_minus matricies 
    ! decribed by Mauersberger [2005] in section 2.2 

    ALLOCATE(tnue_pos0(0) , tnue_neg0(0) )

    DO k = 1, nCycles_red

      ALLOCATE( tnue_pos(0) , tnue_neg(0) )
      DO iS = 1, cycles(k)%len-1
        iSpc = cycles(k)%List(iS)
        cnt_pos = 0
        cnt_neg = 0
        DO jj = nue%RowPtr(iSpc), nue%RowPtr(iSpc+1) - 1
          j = nue%ColInd(jj)
          DO kk = 1,R_k(k)%len
            iR = R_k(k)%List(kk)
            IF ( j == iR ) THEN
              IF ( nue%Val(jj) > 0.0_dp ) THEN
                cnt_pos  = cnt_pos + 1
                tnue_pos = [tnue_pos , jj]
              ELSE
                cnt_neg  = cnt_neg + 1
                tnue_neg = [tnue_neg , jj]
              END IF
            END IF
          END DO
        END DO
      END DO

      tnue_pos0 = [tnue_pos0,nue%ColInd(tnue_pos)]
      tnue_neg0 = [tnue_neg0,nue%ColInd(tnue_neg)]

      DEALLOCATE( tnue_pos, tnue_neg )

    END DO

    ! last matricies nue_pos(0) and nue_neg(0) (NON-CYCLIC REMAINDER)
    ALLOCATE(Perm(SIZE(tnue_pos0)))
    CALL unirnk( tnue_pos0 , Perm , N_pos0 )
    tnue_pos0 = [ tnue_pos0(Perm(1:N_pos0)) ] 
    DEALLOCATE(Perm)

    ALLOCATE(Perm(SIZE(tnue_neg0)))
    CALL unirnk( tnue_neg0 , Perm , N_neg0 )
    tnue_neg0 = [ tnue_neg0(1:N_neg0) ] 
    DEALLOCATE(Perm)

    ALLOCATE(tnue_pos(0),tnue_neg(0))
    cnt_pos = 1
    cnt_neg = 1
    DO iR = 1, nr
      IF ( iR == tnue_pos0(cnt_pos) ) THEN
        IF ( cnt_pos < N_pos0 ) cnt_pos = cnt_pos + 1
      ELSE
        tnue_pos = [tnue_pos , iR]
      END IF
      IF ( iR == tnue_neg0(cnt_neg) ) THEN
        IF ( cnt_neg < N_neg0 ) cnt_neg = cnt_neg + 1
      ELSE
        tnue_neg = [tnue_neg , iR]
      END IF
    END DO
    DEALLOCATE(tnue_pos0, tnue_neg0)
    cnt_pos = SIZE(tnue_pos)
    cnt_neg = SIZE(tnue_neg)

    ! get species indices out of residual reactions set
    ALLOCATE(tnue_pos0(0))
    nS = 0
    DO i = 1, cnt_pos
      iR = tnue_pos(i)
      DO jj = BA%RowPtr(iR),BA%RowPtr(iR+1)-1
        IF ( BA%Val(jj) > 0.0_dp ) THEN
          nS  = nS + 1
          tnue_pos0 = [tnue_pos0 , BA%ColInd(jj)]
        END IF
      END DO
    END DO

    ALLOCATE(Perm(SIZE(tnue_pos0)))
    CALL unirnk( tnue_pos0 , Perm , nS )
    tnue_pos0 = [ tnue_pos0(Perm(1:nS)) ] 
    DEALLOCATE(Perm)

    NonCyclicRemainder = [tnue_pos0]

    ALLOCATE(tnue_neg0(0))
    nS = 0
    DO i = 1, cnt_neg
      iR = tnue_neg(i)
      DO jj = BA%RowPtr(iR),BA%RowPtr(iR+1)-1
        IF ( BA%Val(jj) < 0.0_dp ) THEN
          nS  = nS + 1
          tnue_neg0 = [tnue_neg0 , BA%ColInd(jj)]
        END IF
      END DO
    END DO
    ALLOCATE(Perm(SIZE(tnue_neg0)))
    CALL unirnk( tnue_neg0 , Perm , nS )
    tnue_neg0 = [ tnue_neg0(Perm(1:nS)) ] 
    DEALLOCATE(Perm)
    NonCyclicRemainder = [NonCyclicRemainder,tnue_neg0]


    ALLOCATE(Perm(SIZE(NonCyclicRemainder)))
    CALL unirnk( NonCyclicRemainder , Perm , cnt_pos )
    NonCyclicRemainder = [ NonCyclicRemainder(Perm(1:cnt_pos)) ] 
    DEALLOCATE(Perm)
    nNonCycRem = cnt_pos

    R_k(nCycles_red+1)%len  = nNonCycRem
    R_k(nCycles_red+1)%List = [NonCyclicRemainder]

    !write(*,*) ' noncyclic remainder = ', NonCyclicRemainder
    !stop


  END SUBROUTINE ISSA_structure

  SUBROUTINE SortVecAsc_R(v,q,m)
    USE mo_control, ONLY: big

    REAL(dp), INTENT(INOUT) :: v(:)
    INTEGER,  ALLOCATABLE, OPTIONAL :: q(:)
    INTEGER,               OPTIONAL :: m
    !
    INTEGER :: i,n,im(1)
    REAL(dp), ALLOCATABLE :: tv(:)
   
    IF (PRESENT(m)) THEN
      n = m
    ELSE
      n = SIZE(v)
    END IF
    ALLOCATE(tv(n));  tv = 0
   
    IF (PRESENT(q).AND..NOT.ALLOCATED(q)) ALLOCATE(q(n))
    
    DO i = 1 , n
      im       = MINLOC(v)
      tv(i)    = v(im(1))
      v(im(1)) = big
      IF (PRESENT(q)) q(i) = im(1)
    END DO
    v = tv
  END SUBROUTINE SortVecAsc_R

  SUBROUTINE SortVecAsc_I(v,q)

    INTEGER, INTENT(INOUT) :: v(:)
    INTEGER, ALLOCATABLE, OPTIONAL :: q(:)
    !
    INTEGER :: i,n,im(1),maxv
    INTEGER, ALLOCATABLE :: tv(:)
   
    maxv = MAXVAL(v) + 1
    n = SIZE(v)
    ALLOCATE(tv(n));  tv = 0
   
    IF (PRESENT(q).AND..NOT.ALLOCATED(q)) ALLOCATE(q(n))
    
    DO i = 1 , n
      im       = MINLOC(v)
      tv(i)    = v(im(1))
      v(im(1)) = maxv
      IF (PRESENT(q)) q(i) = im(1)
    END DO
    v = tv
  END SUBROUTINE SortVecAsc_I

  SUBROUTINE ShowChain(C1)
    USE mo_control, ONLY: Chain
    TYPE(Chain) :: C1(:)
    INTEGER :: k, ii, iR, j
   
    write(*,*)
    write(*,*)
    write(*,*) REPEAT('-',30)//' Print the Reaction-Cycles '//REPEAT('-',30)
    write(*,*)

    DO k=1,SIZE(C1)
      write(*,'(A,2X,I0,A,2X,I0)') '  Cycle-Number: ',k, ' with Length: ',SIZE(C1(k)%sIdx,1)
      DO ii=1,SIZE(C1(k)%sName(:,:),1)
        write( * , '(A6,A)' ) '     [' , TRIM(C1(k)%sName(ii,1))//' , '//TRIM(C1(k)%sName(ii,2))//']  involved Reactions: ' 
        DO j=1,SIZE(C1(k)%rIdx(ii)%List)
          write(*,'(10X)',ADVANCE='NO') 
          write(*,'(I0)' ,ADVANCE='NO')  C1(k)%rIdx(ii)%List(j) 
          write(*,'(5X,A21,2F4.1,A2)') 'with Coefficients: [ ', C1(k)%rIdx(ii)%ListE(j), C1(k)%rIdx(ii)%ListP(j), ' ]'
        END DO
      END DO
      write(*,*)
    END DO
     777 FORMAT(10X,A)
  END SUBROUTINE ShowChain

  SUBROUTINE Progress3(j,k)
    INTEGER(4)  :: j,k
    ! print the progress bar.
    WRITE(*,'(A1,A,I0,A,I0,A,$)') char(13),'    Cycle :: (',j,'/',k,')  processed.'
  END SUBROUTINE Progress3

  

!******************************************************************************************************************
!
!                                                                                                                                                                     
!  .d8888b.   .d8888b.  8888888b.  8888888888 8888888888 888b    888 8888888 888b    888  .d8888b.  
! d88P  Y88b d88P  Y88b 888   Y88b 888        888        8888b   888   888   8888b   888 d88P  Y88b 
! Y88b.      888    888 888    888 888        888        88888b  888   888   88888b  888 888    888 
!  "Y888b.   888        888   d88P 8888888    8888888    888Y88b 888   888   888Y88b 888 888        
!     "Y88b. 888        8888888P"  888        888        888 Y88b888   888   888 Y88b888 888  88888 
!       "888 888    888 888 T88b   888        888        888  Y88888   888   888  Y88888 888    888 
! Y88b  d88P Y88b  d88P 888  T88b  888        888        888   Y8888   888   888   Y8888 Y88b  d88P 
!  "Y8888P"   "Y8888P"  888   T88b 8888888888 8888888888 888    Y888 8888888 888    Y888  "Y8888P88
!
!******************************************************************************************************************


  SUBROUTINE ISSA_screening(RS,R_k,S_imp_ini)
    USE mo_control
    USE mo_reac
    USE mo_IO
    USE mo_unirnk
    USE ChemSys_Mod
    USE Sparse_Mod

    TYPE(ReactionStruct_T), INTENT(IN) :: RS(:)         ! struct of reaction system 
    TYPE(List),             INTENT(IN) :: R_k(:)        ! cyclic sets of reactions
    INTEGER,                INTENT(IN) :: S_imp_ini(:)  ! initial set of importent Species

    ! TEMP:
    INTEGER        :: Positions(iStpFlux)
    REAL(dp)       :: Rates(nr,iStpFlux)
    REAL(dp)       :: Rate(nr)
    REAL(dp)       :: time(iStpFlux), dt(iStpFlux)
    INTEGER        :: i, j, jj, k, kk, l, dummy, cnt
    
    INTEGER        :: io_stat = 0
    CHARACTER(200) :: io_msg  = ''

    INTEGER, ALLOCATABLE  :: S_imp(:), R_imp(:)    ! Sets of importent Species/Reactions
    INTEGER, ALLOCATABLE  :: S_imp_new(:)
    INTEGER               :: nS_imp=0, nR_k=0, nIter=0
    INTEGER               :: iS=0, iSpc=0, iReac=0

    TYPE(CSR_Matrix_T)    :: pos_BAT, neg_BAT, neg_BA

    INTEGER               :: rp_iSpc(2), rp_iSpcP1(2) ! rowpointer: 1 for source terms, 2 for sink terms
    INTEGER,  ALLOCATABLE :: source_reacs(:), sink_reacs(:)   ! stoech. coefficientes of cycle k
    LOGICAL               :: any_sources, any_sinks

    INTEGER               :: inR_f, inR_g

    REAL(dp), ALLOCATABLE :: f_iJk(:), g_iJk(:)   ! valuation coefficientes (noncyclic reactions)
    REAL(dp)              :: sum_f_ijk, sum_g_ijk
    INTEGER,  ALLOCATABLE :: CF_ik(:), CG_ik(:)  ! max. member groups of redundant reactions (index sets)
    
    INTEGER,  ALLOCATABLE :: f_p(:),  g_p(:)    ! permuation vector of reaction set
    INTEGER,  ALLOCATABLE :: f_iR(:), g_iR(:)    ! permuation vector of reaction set

    REAL(dp), PARAMETER   :: eps_red = 0.11
    INTEGER,  PARAMETER   :: maxIter = 42


    INTEGER, ALLOCATABLE :: NewReactionSet(:), Perm(:)
    INTEGER              :: newLen

          
    WRITE(*,777) REPEAT('*',39)
    WRITE(*,777) '**** Iterative Screening Procedure ****'
    WRITE(*,777) REPEAT('*',39)
    WRITE(*,*)

    ! reading meta data (positions of record)
    CALL OpenFile_rSeq(FluxMetaUnit,FluxMetaFile)
    DO i = 1,iStpFlux
      READ(FluxMetaUnit,*,IOSTAT=io_stat,IOMSG=io_msg) dummy , Positions(i)
      IF ( io_stat>0 ) WRITE(*,'(10X,A,I0,A)') '   ERROR :: ',io_stat,'  '//TRIM(io_msg)
      IF ( io_stat<0 ) EXIT
    END DO
    CLOSE(FluxMetaUnit)

    ! reading unformatted binary file
    CALL OpenFile_rStream(FluxUnit,FluxFile)
    DO i = 1,iStpFlux
      READ(FluxUnit,POS=Positions(i),IOSTAT=io_stat,IOMSG=io_msg) Rates(:,i) , time(i) , dt(i)
      IF ( io_stat>0 ) WRITE(*,'(10X,A,I0,A)') '   ERROR :: ',io_stat,'  '//TRIM(io_msg)
      IF ( io_stat<0 ) EXIT
    END DO
    CLOSE(FluxUnit)

   

    ! calculate time-averaged rate of reactions
    Rate = 0_dp
    DO i = 1,iStpFlux
      Rate = Rate + ABS(Rates(:,i))*dt(i)
    END DO
    Rate = Rate/REAL(iStpFlux)


    OPEN(unit=199, file='Rates_'//BSP//'.txt', action='write')
    DO j=1,SIZE(Rate) 
      WRITE(199,'(A,I0,A,Es16.8,A)') ' integr. rate(',j,') = ',rate(j), '    '//TRIM(RS(j)%Type)//'    '//TRIM(RS(j)%Line1)
    END DO
    CLOSE(199)


    ! vorbereitend: positive und negative elemente separieren
    CALL SeperatePosNegValues( pos_BAT , neg_BAT , BAT)
    neg_BAT%Val = ABS(neg_BAT%Val)
    CALL TransposeSparse(neg_BA,neg_BAT)

    S_imp  = [S_imp_ini]
    nS_imp = SIZE(S_imp_ini)
    ALLOCATE( R_imp(0) )

    nR_k = SIZE(R_k)

    ! ITERATIVE LOOP
    nIter = 1

    ITERATIVE_PROCEDURE: DO !WHILE (nIter <= maxIter)
      
      CALL Progress_Step(nIter,nS_imp)

      !******************************************************************
      ! (a) For the actual group of important species (index set S_imp) 
      !     the valuation coefficients f_ijk, g_ijk are calculated. At
      !     the start S_imp contains the target species only.
      !****************************************************************** 

      DO iS = 1 , nS_imp

        iSpc       = S_imp(iS)
        rp_iSpc(1) = pos_BAT%RowPtr(iSpc);  rp_iSpcP1(1) = pos_BAT%RowPtr(iSpc+1)-1
        rp_iSpc(2) = neg_BAT%RowPtr(iSpc);  rp_iSpcP1(2) = neg_BAT%RowPtr(iSpc+1)-1
      
        DO k = 1, nR_k

          ! --- SOURCE terms
          ALLOCATE(source_reacs(0)); any_sources = .FALSE.
          DO jj = rp_iSpc(1) , rp_iSpcP1(1)
            IF ( ANY( pos_BAT%ColInd(jj) == R_k(k)%List ) ) THEN
              source_reacs = [source_reacs , jj]
              any_sources = .TRUE.
            END IF
          END DO
          
          IF (any_sources) THEN
            f_iR  = [ pos_BAT%ColInd(source_reacs) ]
            inR_f = SIZE(f_iR)
            f_iJk = [ pos_BAT%Val(source_reacs) * Rate(f_iR) ]
            f_iJk = f_iJk / SUM(f_iJk)

            CALL SortVecAsc_R( f_iJk , f_p , inR_f )
            f_iR = f_iR(f_p)
            DEALLOCATE(f_p)
            !DO j=1,inR_f;  WRITE(*,'(A,Es16.8)') ' f_iJk = ',f_iJk(j);  END DO
          ELSE
            inR_f = 0
          END IF


          ! --- SINK terms
          ALLOCATE(sink_reacs(0)); any_sinks = .FALSE.
          DO jj = rp_iSpc(2) , rp_iSpcP1(2)
            IF ( ANY( neg_BAT%ColInd(jj) == R_k(k)%List ) ) THEN
              sink_reacs = [sink_reacs , jj]
              any_sinks = .TRUE.
            END IF
          END DO
          
          IF (any_sinks) THEN
            g_iR = [ neg_BAT%ColInd(sink_reacs) ]
            inR_g = SIZE(g_iR)
            g_iJk = [ pos_BAT%Val(sink_reacs) * Rate(g_iR) ]
            g_iJk = g_iJk / SUM(g_iJk)

            CALL SortVecAsc_R( g_iJk , g_p , inR_g )
            g_iR = g_iR(g_p)
            DEALLOCATE(g_p)
            !DO j=1,inR_g;  WRITE(*,'(A,Es16.8)') ' g_iJk = ',g_iJk(j);  END DO
          ELSE
            inR_g = 0
          END IF

          
          !*******************************************************************
          ! (b) The maximum member groups of redundant reactions (index sets 
          !     F_ik, G_ik) with the property SUM(f_ijk) < eps_red, and 
          !     SUM(f_ijk) < eps_red are determined. Especially reactions with 
          !     f_ijk=0, g_ijk=0 are always part of F_ik, G_ik, respectivley.
          !     The threshold value with 0 <= eps_red <= 1 controls the 
          !     reduction intensity.
          !*******************************************************************

          ALLOCATE(CF_ik(0),CG_ik(0))
          sum_f_ijk = 0.0_dp
          DO j = 1,inR_f
            IF ( sum_f_iJk < eps_red ) THEN
              sum_f_ijk = sum_f_ijk + f_iJk(j)
            ELSE
              CF_ik = [(f_iR(l) , l=j-1,inR_f)]
              EXIT
            END IF
          END DO

          sum_g_ijk = 0.0_dp
          DO j = 1,inR_g
            IF ( sum_g_ijk < eps_red ) THEN
              sum_g_ijk = sum_g_ijk + g_iJk(j)
            ELSE
              CG_ik = [(g_iR(l) , l=j-1,inR_g)]
              EXIT
            END IF
          END DO

          !*******************************************************************
          ! (c) The important reactions (index set R_imp) of important species 
          !     in S_imp are calculated by R_imp = set_add(CF_ik , CG_ik)
          !     where i in S_imp.
          !*******************************************************************

          R_imp = [ R_imp , CF_ik , CG_ik ]

          DEALLOCATE( CF_ik, CG_ik, source_reacs, sink_reacs )

        END DO   !  k = 1, nCycles
        !READ(*,*)
      END DO     !  i = 1 , nS_imp

      !*******************************************************************************
      ! (d) The reactants of R_imp (species i with nue_ij<0 for j in R_imp) form
      !     the new set of important species S_imp_new with S_imp_new >= S_imp. If
      !     S_imp_new > S_imp then the iteration goes on with S_imp = S_imp_new
      !     in step (a). In the other case the iteration is finished; S_imp and R_imp
      !     contain the important species and reactions of the reduced mechanism.
      !*******************************************************************************
      ALLOCATE(S_imp_new(0))
      DO j = 1,SIZE(R_imp)
        jj = R_imp(j)
        S_imp_new = [S_imp_new ,                                          &
        &           neg_BA%ColInd(neg_BA%RowPtr(jj):neg_BA%RowPtr(jj+1)-1)]
      END DO
      

      ! --- set addition -> sort -> remove duplicates
      S_imp_new = [S_imp , S_imp_new]
      CALL Sort_And_Remove_Duplicates(S_imp_new,nS_imp)

      IF ( nS_imp == SIZE(S_imp) ) THEN
        EXIT ITERATIVE_PROCEDURE
      ELSE
        S_imp = [S_imp_new]
        DEALLOCATE(S_imp_new)
        nIter = nIter + 1
      END IF

    END DO ITERATIVE_PROCEDURE

    CALL Sort_And_Remove_Duplicates(R_imp,newLen)

    !DO j=1,SIZE(R_imp); WRITE(*,'(A,I0,A,I0)') ' R_imp(',j,') = ',R_imp(j); END DO
    !DO j=1,SIZE(S_imp); WRITE(*,'(A,I0,A,I0)') ' S_imp(',j,') = ',S_imp(j); END DO
    
    !-----------------------------------------------------------------------
    ! --- Writing a new sys-File with less reactions
    CALL Print_SysFile(ReactionSystem,R_imp,'CHEM/'//TRIM(BSP)//'_red.sys')
    write(*,*); write(*,*) 
    write(*,777) '    Printing reduced system  ::  CHEM/'//TRIM(BSP)//'_red.sys'
    WRITE(*,*)
    write(*,'(10X,A,I6)') '    Reduced system :: Number of Reactions = ',newLen
    write(*,'(10X,A,I6)') '                      Number of Species   = ',nS_imp
    !stop
    777 FORMAT(10X,A)


    !STOP ' UNDER CONSTRUCTION '

    CONTAINS

      SUBROUTINE Progress_Step(j,k)
        INTEGER :: j,k
        ! print the progress bar.
        WRITE(*,'(14X,2(A,I0,4X))') 'Iteration ::  ',j,'nSpecies = ',k
        !WRITE(*,'(A1,14X,2(A,I0,4X),$)') char(13),'Iteration ::  ',j,'nSpecies = ',k
      END SUBROUTINE Progress_Step
      
      FUNCTION Is_ImpSpc_In_Cycle(iSpc,List,n) RESULT(iRow)
        INTEGER             :: iRow
        INTEGER, INTENT(IN) :: List(:)
        INTEGER, INTENT(IN) :: iSpc, n
        INTEGER :: i

        iRow=0
        DO i=1,n
          IF ( List(i)==iSpc ) THEN
            iRow=i;  RETURN
          END IF
        END DO
      END FUNCTION Is_ImpSpc_In_Cycle

      SUBROUTINE Sort_And_Remove_Duplicates(Array,n)
        INTEGER, ALLOCATABLE, INTENT(INOUT) :: Array(:)
        INTEGER, ALLOCATABLE                :: Perm(:)
        INTEGER                             :: n

        ALLOCATE(Perm(SIZE(Array)))
        CALL unirnk( Array , Perm , n )
        Array = [Array(Perm(1:n))]
        DEALLOCATE(Perm)
      END SUBROUTINE Sort_And_Remove_Duplicates

  END SUBROUTINE ISSA_screening

  
END MODULE issa
