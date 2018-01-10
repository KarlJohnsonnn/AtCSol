MODULE issa
  
  USE Kind_Mod
  USE Sparse_Mod, ONLY: CSR_Matrix_T

  IMPLICIT NONE

  TYPE(CSR_Matrix_T) :: nue_pos_0, nue_neg_0

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
    
    WRITE(*,'(7X,A)') REPEAT('*',39)
    WRITE(*,'(7X,A)') '********** Important Species **********'
    WRITE(*,'(7X,A)') REPEAT('*',39)
    WRITE(*,*)
    DO j=1,SIZE(Idx)
      WRITE(*,'(10X,A,I2,A,I6,5X,A)') '- S_imp(',j,') = ',Idx(j),TRIM(y_name(Idx(j)))
    END DO
    WRITE(*,*); WRITE(*,*)
    
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

    WRITE(*,'(7X,A)') REPEAT('*',39)
    WRITE(*,'(7X,A)') '*********** Species Familes ***********'
    WRITE(*,'(7X,A)') REPEAT('*',39)
    WRITE(*,*)
    DO j=1,SIZE(Fam)
      WRITE(*,'(10X,A,I2,A,*(A))')      '- Family(',j,')  -->  ',    &
      &    ( TRIM(Fam(j)%Name(i))//' , ' , i=1,SIZE(Fam(j)%Name)-1 ),&
      &      TRIM(Fam(j)%Name(SIZE(Fam(j)%Name)))
    END DO
    WRITE(*,*); WRITE(*,*)

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
        WRITE(*,*) '  No '//TRIM(SectionName)//' found.'
        iLine = -1
        EXIT
      END IF
    END DO
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


!******************************************************************************************************************
!
!8888888 .d8888b.   .d8888b.        d8888       .d8888b. 88888888888 8888888b.  888     888  .d8888b. 88888888888 
!  888  d88P  Y88b d88P  Y88b      d88888      d88P  Y88b    888     888   Y88b 888     888 d88P  Y88b    888     
!  888  Y88b.      Y88b.          d88P888      Y88b.         888     888    888 888     888 888    888    888     
!  888   "Y888b.    "Y888b.      d88P 888       "Y888b.      888     888   d88P 888     888 888           888     
!  888      "Y88b.     "Y88b.   d88P  888          "Y88b.    888     8888888P"  888     888 888           888     
!  888        "888       "888  d88P   888            "888    888     888 T88b   888     888 888    888    888     
!  888  Y88b  d88P Y88b  d88P d8888888888      Y88b  d88P    888     888  T88b  Y88b. .d88P Y88b  d88P    888     
!8888888 "Y8888P"   "Y8888P" d88P     888       "Y8888P"     888     888   T88b  "Y88888P"   "Y8888P"     888    
!
!******************************************************************************************************************


!  SUBROUTINE ISSA_nue(nue,cyc,RS)
!    USE mo_reac,     ONLY: neq, nspc, y_name
!    USE mo_control,  ONLY: List, Chain, ZERO
!    USE Sparse_Mod,  ONLY: CSR_Matrix_T, New_CSR, CSR_to_SpRowIndColInd,  &
!    &                   SpRowIndColInd_to_CSR, PrintSparse2, PrintSparse3
!    USE ChemSys_Mod, ONLY: ReactionStruct_T
!    ! IN:
!    TYPE(CSR_Matrix_T)      :: nue
!    TYPE(List), ALLOCATABLE :: cyc(:)
!    TYPE(ReactionStruct_T), ALLOCATABLE :: RS(:)
!    ! TEMP:
!    INTEGER :: i, ii, iR, iE, iP, j, jj,  k, n_c, len_c, nE, nP, nnz_k
!    INTEGER,     ALLOCATABLE :: tsi(:), perm(:)
!    INTEGER,     ALLOCATABLE :: ri_p(:), ci_p(:), ri_m(:), ci_m(:)
!    REAL(dp),    ALLOCATABLE :: v_p(:), v_m(:)
!    TYPE(Chain), ALLOCATABLE :: R_k(:)
!    TYPE(SpRowIndColInd_T)   :: sp_nue
!
!    ! Debugging Formats
!
!    n_c = SIZE(cyc)
!
!    write(*,*) '_____________________________________________'
!    write(*,*)
!    write(*,*) '  Testing ISSA stuff: building nue+ and nue- '
!    write(*,*) '_____________________________________________'
!    write(*,*)
!
!    ! build nue_p and nue_m (see atmospheric env. 39 (2005) page: 4343)
!    ALLOCATE( nue_p(0:n_c) , nue_m(0:n_c) , R_k(0:n_c) )
!
!    sp_nue = CSR_to_SpRowIndColInd(nue)
!
!    LOOP_OVER_ALL_CYCLES: DO k = 1,n_c ! Kette Nummer: k
!      
!      ! First we need to detect the reactions of each species cycle
!      len_c = cyc(k)%len
!      tsi   = cyc(k)%List
!
!      ALLOCATE( ri_p(0), ci_p(0), v_p(0), ri_m(0), ci_m(0), v_m(0) )
!
!      !ALLOCATE( R_k(k)%sName(len_c-1,2), R_k(k)%sIdx(len_c-1,2), R_k(k)%rIdx(len_c-1) )
!
!      !write(*,*) REPEAT('-',32)//'  Browse Cycle: ',k,'  with length: ',len_c
!      !write(*,'(5X,*(A))')  ( TRIM(y_name(cyc(k)%List(j)))//' -> ' , j=1,len_c-1 ) ,&
!      !&                       TRIM(y_name(cyc(k)%List(len_c)))
!      !write(*,*) 
!
!      DO ii = 1,len_c-1
!        !ALLOCATE( R_k(k)%rIdx(ii)%List(0), R_k(k)%rIdx(ii)%ListE(0), R_k(k)%rIdx(ii)%ListP(0) ) 
!        DO iR = 1,neq
!
!          nE = SIZE(RS(iR)%Educt)
!          nP = SIZE(RS(iR)%Product)
!
!          DO iE = 1,nE
!            IF ( tsi(ii) == RS(iR)%Educt(iE)%iSpecies ) THEN
!              !R_k(k)%sName(ii,1) = RS(iR)%Educt(iE)%Species
!              !R_k(k)%sIdx(ii,1)  = RS(iR)%Educt(iE)%iSpecies
!            
!              DO iP = 1,nP
!                IF ( tsi(ii+1) == RS(iR)%Product(iP)%iSpecies ) THEN
!                  !R_k(k)%sName(ii,2)   = TRIM(RS(iR)%Product(iP)%Species)
!                  !R_k(k)%sIdx(ii,2)    = RS(iR)%Product(iP)%iSpecies
!                  !R_k(k)%rIdx(ii)%List = [R_k(k)%rIdx(ii)%List , iR]
!                  !R_k(k)%rIdx(ii)%ListE = [-RS(iR)%Educt(iE)%Koeff,  R_k(k)%rIdx(ii)%ListE]
!                  !R_k(k)%rIdx(ii)%ListP = [RS(iR)%Product(iP)%Koeff, R_k(k)%rIdx(ii)%ListP]
!                  !write(*,'(A14,I0,A)',    ADVANCE='NO') '  Reaction :: ',iR,':  '//TRIM(RS(iR)%Line1)
!                  !write(*,'(A9,F4.1,1X,A)',ADVANCE='NO') '  Educt: ', RS(iR)%Educt(iP)%Koeff,TRIM(RS(iR)%Educt(iE)%Species)
!                  !write(*,'(A15,F4.1,1X,A)')             '  and Product: ', RS(iR)%Product(iP)%Koeff,TRIM(RS(iR)%Product(iP)%Species)
!
!                  !write(*,'(A,I0,A3,2X,*(I0,5X))') '  spc in (beta-alpha)_(iR,j) = ',iR,' | ',nue%ColInd(nue%RowPtr(iR):nue%RowPtr(iR+1)-1)
!                  !write(*,'(A,I0,A3,*(F4.1,2X))')  '  val in (beta-alpha)_(iR,j) = ',iR,' | ',nue%Val(nue%RowPtr(iR):nue%RowPtr(iR+1)-1)
!
!                  DO jj = nue%RowPtr(iR),nue%RowPtr(iR+1)-1
!                    IF (nue%ColInd(jj)==RS(iR)%Educt(iE)%iSpecies .OR. nue%ColInd(jj)==RS(iR)%Product(iP)%iSpecies) THEN
!                      !write(*,'(A,I0,A1,I0,A4,F5.1)') '  spc in (beta-alpha)_(',iR,',',nue%ColInd(jj),') = ',nue%Val(jj)
!                      IF ( nue%Val(jj) > ZERO ) THEN
!                        ci_p = [ci_p, iR]
!                        ri_p = [ri_p, nue%ColInd(jj)]
!                        v_p  = [v_p, ABS(nue%Val(jj))]
!                      ELSE
!                        ci_m = [ci_m, iR]
!                        ri_m = [ri_m, nue%ColInd(jj)]
!                        v_m  = [v_m, ABS(nue%Val(jj))]
!                      END IF
!                    END IF
!                  END DO
!
!                END IF 
!              END DO
!
!            END IF
!          END DO
!
!        END DO
!      END DO
!      DEALLOCATE(tsi)
!
!      CALL Progress3(k,n_c)
!
!      ! build the sparse matrices nue_plus_ijk and nue_minus_ijk
!      CALL SortVecAsc_I(ri_p,perm)
!      ci_p = ci_p(perm)
!      v_p  = v_p(perm)
!      nnz_k = SIZE(ri_p)
!      nue_p(k) = New_CSR(nue%n,nue%m,nnz_k,ri_p,ci_p,v_p)
!      DEALLOCATE(ri_p,ci_p,v_p,perm)
!
!      CALL SortVecAsc_I(ri_m,perm)
!      ci_m = ci_m(perm)
!      v_m  = v_m(perm)
!      nnz_k = SIZE(ri_m)
!      nue_m(k) = New_CSR(nue%n,nue%m,nnz_k,ri_m,ci_m,v_m)
!      DEALLOCATE(ri_m,ci_m,v_m,perm)
!
!      !CALL PrintSparse3(nue_p(k)%m, nue_p(k)%n, nue_p(k)%RowPtr, nue_p(k)%ColInd , nue_p(k)%Val , 'nue_p')
!      !CALL PrintSparse3(nue_m(k)%m, nue_m(k)%n, nue_m(k)%RowPtr, nue_m(k)%ColInd , nue_m(k)%Val , 'nue_m')
!
!      !write(*,*) '  Cycle: ',k, '  with nnz_p = ',nue_p(k)%nnz, '  and nnz_m = ',nue_m(k)%nnz
!
!    END DO LOOP_OVER_ALL_CYCLES
!    WRITE(*,*)
!
!    !CALL ShowChain(R_k)
!
!  END SUBROUTINE ISSA_nue

  SUBROUTINE ISSA_nue(nue)
    USE Sparse_Mod, ONLY: New_CSR, WriteSparseMatrix, Copy_CSR
    TYPE(CSR_Matrix_T) :: nue

    INTEGER, ALLOCATABLE :: tnue_pos(:)
    INTEGER, ALLOCATABLE :: tnue_neg(:)
    INTEGER :: cnt_pos, cnt_neg
    INTEGER :: i, jj, j 


    ! this routine generates the nue_plus and nue_minus matricies 
    ! decribed by Mauersberger [2005] in section 2.2 
    nue_pos_0 = New_CSR(nue%m,nue%n)
    nue_neg_0 = New_CSR(nue%m,nue%n)

    ALLOCATE( tnue_pos(0), tnue_neg(0) )

    DO i = 1 , nue%m
      cnt_pos = 0
      cnt_neg = 0
      DO jj = nue%RowPtr(i) , nue%RowPtr(i+1) - 1
        j = nue%ColInd(jj)

        IF ( nue%Val(jj) > 0.0_dp ) THEN
          cnt_pos = cnt_pos + 1
          tnue_pos = [tnue_pos , jj]
        ELSE
          cnt_neg = cnt_neg + 1
          tnue_neg = [tnue_neg , jj]
        END IF

      END DO
      nue_pos_0%RowPtr(i+1) = nue_pos_0%RowPtr(i) + cnt_pos
      nue_neg_0%RowPtr(i+1) = nue_neg_0%RowPtr(i) + cnt_neg
    END DO

    nue_pos_0%nnz = nue_pos_0%RowPtr(nue_pos_0%m+1) - 1
    nue_neg_0%nnz = nue_neg_0%RowPtr(nue_neg_0%m+1) - 1

    nue_pos_0%ColInd = [ nue%ColInd(tnue_pos) ]
    nue_neg_0%ColInd = [ nue%ColInd(tnue_neg) ]

    nue_pos_0%Val = [ nue%Val(tnue_pos) ]
    nue_neg_0%Val = [ nue%Val(tnue_neg) ]

  END SUBROUTINE ISSA_nue

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
  END SUBROUTINE ShowChain

  SUBROUTINE Progress3(j,k)
    INTEGER(4)  :: j,k
    ! print the progress bar.
    WRITE(*,'(A1,A,I0,A,I0,A,$)') char(13),'    Cycle :: (',j,'/',k,')  processed.'
  END SUBROUTINE Progress3

  

!******************************************************************************************************************
!
! 8888888888 888     888     888 Y88b   d88P                         d8888 888b    888        d8888 
! 888        888     888     888  Y88b d88P                         d88888 8888b   888       d88888 
! 888        888     888     888   Y88o88P                         d88P888 88888b  888      d88P888 
! 8888888    888     888     888    Y888P                         d88P 888 888Y88b 888     d88P 888 
! 888        888     888     888    d888b                        d88P  888 888 Y88b888    d88P  888 
! 888        888     888     888   d88888b        888888        d88P   888 888  Y88888   d88P   888 
! 888        888     Y88b. .d88P  d88P Y88b                    d8888888888 888   Y8888  d8888888888 
! 888        88888888 "Y88888P"  d88P   Y88b                  d88P     888 888    Y888 d88P     888
!
!******************************************************************************************************************


  SUBROUTINE FluxAnalysis(S_imp)
    USE mo_control
    USE mo_reac
    USE mo_IO
    USE mo_unirnk
    USE ChemSys_Mod
    USE Sparse_Mod, ONLY: A     ! educt matrix (nr x nspc)

    ! TEMP:
    INTEGER        :: Positions(iStpFlux)
    REAL(dp)       :: Rates(nr,iStpFlux)
    REAL(dp)       :: Rate(nr)
    REAL(dp)       :: time(iStpFlux), dt(iStpFlux)
    INTEGER        :: i, j, jj, k, l, dummy
    
    INTEGER        :: io_stat = 0
    CHARACTER(200) :: io_msg  = ''

    INTEGER, ALLOCATABLE  :: S_imp(:), R_imp(:)    ! Sets of importent Species/Reactions
    INTEGER, ALLOCATABLE  :: S_imp_new(:)
    INTEGER               :: nS_imp=0, nR_imp=0, nIter=0
    INTEGER               :: iS=0

    INTEGER               :: rp_iS, rp_iS1, innz_f, innz_g
    REAL(dp), ALLOCATABLE :: f_ij0(:), g_ij0(:)   ! valuation coefficientes (noncyclic reactions)
    REAL(dp)              :: sum_f_ij0, sum_g_ij0 
    INTEGER,  ALLOCATABLE :: CF_i0(:), CG_i0(:)  ! max. member groups of redundant reactions (index sets)
    
    INTEGER,  ALLOCATABLE :: f_p(:),    g_p(:)    ! permuation vector of reaction set
    INTEGER,  ALLOCATABLE :: f_iR(:), g_iR(:)    ! permuation vector of reaction set

    REAL(dp)              :: threshold
    REAL(dp), PARAMETER   :: eps_red = 0.11
    INTEGER,  PARAMETER   :: versatz = 44


    INTEGER, ALLOCATABLE :: NewReactionSet(:), Perm(:)
    INTEGER              :: newLen

          
    WRITE(*,'(7X,A)') REPEAT('*',39)
    WRITE(*,'(7X,A)') '**** Iterative Screening Procedure ****'
    WRITE(*,'(7X,A)') REPEAT('*',39)
    WRITE(*,*)

    ! reading meta data (positions of record)
    CALL OpenFile_rSeq(FluxMetaUnit,FluxMetaFile)
    DO i = 1,iStpFlux
      READ(FluxMetaUnit,*,IOSTAT=io_stat,IOMSG=io_msg) dummy , Positions(i)
      IF ( io_stat>0 ) WRITE(*,*) '   ERROR :: ',io_stat,'  '//TRIM(io_msg)
      IF ( io_stat<0 ) EXIT
    END DO
    CLOSE(FluxMetaUnit)

    ! reading unformatted binary file
    CALL OpenFile_rStream(FluxUnit,FluxFile)
    DO i = 1,iStpFlux
      READ(FluxUnit,POS=Positions(i),IOSTAT=io_stat,IOMSG=io_msg) Rates(:,i) , time(i) , dt(i)
      IF ( io_stat>0 ) WRITE(*,*) '   ERROR :: ',io_stat,'  '//TRIM(io_msg)
      IF ( io_stat<0 ) EXIT
    END DO
    CLOSE(FluxUnit)


    nS_imp = SIZE(S_imp)
    ALLOCATE( R_imp(0) )

    ! calculate time-averaged rate of reactions
    Rate = 0_dp
    DO i = 1,iStpFlux-versatz
      Rate = Rate + Rates(:,i)*dt(i)
    END DO
    Rate = Rate/REAL(iStpFlux)

    DO j=1,SIZE(Rate); WRITE(*,'(A,I0,A,Es16.8)') ' rate(',j,') = ',rate(j); END DO

    WRITE(*,'(10X,A)') 'Starting iterative ISSA procedure'


    ! ITERATIVE LOOP
    ITERATIVE_PROCEDURE: DO

    !******************************************************************
    ! (a) For the actual group of important species (index set S_imp) 
    !     the valuation coefficients f_ijk, g_ijk are calculated. At
    !     the start S_imp contains the target species only.
    !****************************************************************** 
      DO i = 1 , nS_imp

        iS = S_imp(i)

        ! --- calculate valuation coefficient f_ijk
        rp_iS  = nue_pos_0%RowPtr(iS)
        rp_iS1 = nue_pos_0%RowPtr(iS+1)-1
        innz_f = rp_iS1 - rp_iS + 1
        f_iR   = [ nue_pos_0%ColInd(rp_iS:rp_iS1) ]

        f_ij0 = [ nue_pos_0%Val(rp_iS:rp_iS1) * Rate(nue_pos_0%ColInd(rp_iS:rp_iS1)) ]
        write(*,*) ' sum fij0 = ',SUM(f_ij0)
        f_ij0 = f_ij0 / SUM(f_ij0)

        CALL SortVecAsc_R( f_ij0 , f_p , innz_f )
        f_iR = f_iR(f_p)
        DEALLOCATE(f_p)
        !write(6,*)  SIZE(f_ij0,1)
        !DO jj=1,SIZE(f_ij0);  write(*,*)  f_ij0(jj); END DO

        ! --- calculate valuation coefficient g_ijk
        rp_iS  = nue_neg_0%RowPtr(iS)
        rp_iS1 = nue_neg_0%RowPtr(iS+1)-1
        innz_g = rp_iS1 - rp_iS + 1
        g_iR   = [ nue_neg_0%ColInd(rp_iS:rp_iS1) ]

        g_ij0 = [ nue_neg_0%Val(rp_iS:rp_iS1) * Rate(nue_neg_0%ColInd(rp_iS:rp_iS1)) ]
        write(*,*) ' sum gij0 = ',SUM(g_ij0)
        g_ij0 = g_ij0 / SUM(g_ij0)

        CALL SortVecAsc_R( g_ij0 , g_p , innz_g )
        g_iR = g_iR(g_p)
        DEALLOCATE(g_p)
        !write(*,*)  SIZE(g_ij0,1)
        !DO jj=1,SIZE(g_ij0);  write(*,*)  g_ij0(jj); END DO

        !*******************************************************************
        ! (b) The maximum member groups of redundant reactions (index sets 
        !     F_ik, G_ik) with the property SUM(f_ijk) < eps_red, and 
        !     SUM(f_ijk) < eps_red are determined. Especially reactions with 
        !     f_ijk=0, g_ijk=0 are always part of F_ik, G_ik, respectivley.
        !     The threshold value with 0 <= eps_red <= 1 controls the 
        !     reduction intensity.
        !*******************************************************************

        ALLOCATE(CF_i0(0),CG_i0(0))
        sum_f_ij0 = 0.0_dp
        DO j = 1,innz_f
          IF ( sum_f_ij0 < eps_red ) THEN
            sum_f_ij0 = sum_f_ij0 + f_ij0(j)
          ELSE
            CF_i0 = [(f_iR(l) , l=j-1,innz_f)]
            WRITE(*,*) ' NEW REACIONS :: CF_i0 :: ',CF_i0
            EXIT
          END IF
        END DO

        sum_g_ij0 = 0.0_dp
        DO j = 1,innz_g
          IF ( sum_g_ij0 < eps_red ) THEN
            sum_g_ij0 = sum_g_ij0 + g_ij0(j)
          ELSE
            CG_i0 = [(g_iR(l) , l=j-1,innz_g)]
            WRITE(*,*) ' NEW REACIONS :: CG_i0 :: ',CG_i0
            EXIT
          END IF
        END DO

        !*******************************************************************
        ! (c) The important reactions (index set R_imp) of important species 
        !     in S_imp are calculated by R_imp = set_add(CF_ik , CG_ik)
        !     where i in S_imp.
        !*******************************************************************

        R_imp = [ R_imp , CF_i0 , CG_i0 ]
        DEALLOCATE(CF_i0, CG_i0)
      END DO 


      !*******************************************************************
      ! (d) The reactants of R_imp (species i with nue_ij<0 for j in R_imp)
      !     form the new set of important species S*_imp with S*imp >= S_imp
      !     If S*_imp > S_imp then the iteration goes on with S_imp = S*_imp
      !     in step (a). In the other case the iteration is finished; S_imp
      !     and R_imp contain the important species and reactions of the 
      !     reduced mechanism.
      !*******************************************************************
      ALLOCATE(S_imp_new(0))
      DO j = 1,SIZE(R_imp)
        jj = R_imp(j)
        rp_iS   = A%RowPtr(jj)
        rp_iS1  = A%RowPtr(jj+1)-1
        S_imp_new = [S_imp_new , A%ColInd(rp_iS:rp_iS1)]
      END DO
      

      ! --- set addition -> sort -> remove duplicates
      S_imp_new = [S_imp , S_imp_new]
      ALLOCATE(Perm(SIZE(S_imp_new)))

      CALL unirnk( S_imp_new , Perm , nS_imp )

      !WRITE(*,'(4X,A,I0)') '      nIter = ', nIter
      !WRITE(*,'(I4,A,*(I0,3X))') SIZE(S_imp_new),'      S_imp = ', (S_imp_new(l), l=1,SIZE(S_imp_new))
      !WRITE(*,'(I4,A,*(I0,3X))') SIZE(R_imp),    '      R_imp = ', (R_imp(l), l=1,SIZE(R_imp))


      S_imp_new = [S_imp_new(Perm(1:nS_imp))]
      DEALLOCATE(Perm)

      IF ( nS_imp == SIZE(S_imp)) THEN
        EXIT ITERATIVE_PROCEDURE
      ELSE
        nIter = nIter + 1
        S_imp = [S_imp_new]
        DEALLOCATE(S_imp_new)
        !CALL Progress_Step(nIter,nS_imp,SIZE(R_imp))
      END IF

    END DO ITERATIVE_PROCEDURE
    
    WRITE(*,*); WRITE(*,*)
    WRITE(*,'(10X,A)') 'Iterative procedure done '
    !WRITE(6,'(3X,I0)') S_imp


    !DO j=1,SIZE(f_ij0); WRITE(*,'(A,I0,A,Es16.8)') ' f_ij0_sorted(',j,') = ',f_ij0(j); END DO
    !WRITE(*,*) ' size CF_i0 = ', SIZE(CF_i0)
    !WRITE(*,'(A,I0,*(3X,I0))') ' CF_i0 = ',CF_i0

    !DO j=1,SIZE(g_ij0); WRITE(*,'(A,I0,A,Es16.8)') ' g_ij0_sorted(',j,') = ',g_ij0(j); END DO
    !WRITE(*,*) ' size CG_i0 = ', SIZE(CG_i0)
    !WRITE(*,'(A,I0,*(3X,I0))') ' CG_i0 = ',CG_i0
    !stop

  
    !DO j=1,SIZE(f_ij0); WRITE(*,'(A,I0,A,Es16.8)') ' f_ij0_sorted(',j,') = ',f_ij0(j); END DO



!    ! this loop is for all the reaction cycles
!    DO k = 1 , nCycles
!
!      DO i = 1 , nS_imp
!
!        iS = S_imp(i)
!        
!        !sum_nuep_r = SUM()
!
!        DO j = 1 , nr
!
!        END DO
!
!      END DO
!
!    END DO




    !ALLOCATE(NewReactionSet(0))
    !DO i = 1,iStpFlux
    !
    !  CALL Update_ReactionSet(NewReactionSet,Rates(:,i)*dt(i))
    !
    !  CALL Progress4(i,iStpFlux)
    !
    !END DO

    DO j=1,SIZE(R_imp); WRITE(*,'(A,I0,A,I0)') ' R_imp(',j,') = ',R_imp(j); END DO

    ALLOCATE(Perm(SIZE(R_imp)))
    CALL unirnk(R_imp,Perm,newLen)
    R_imp = [ R_imp(Perm(1:newLen)) ]
    DEALLOCATE(Perm)


    !-----------------------------------------------------------------------
    ! --- Writing a new sys-File with less reactions
    CALL Print_SysFile(ReactionSystem,R_imp,'CHEM/'//TRIM(BSP)//'_red.sys')
    write(*,*); write(*,*) 
    write(*,'(A)') '  Printing New Reaction System done  ::  CHEM/'//TRIM(BSP)//'_red.sys'
    !stop


    CONTAINS
      SUBROUTINE Update_ReactionSet(NewSet,Rates)
        USE mo_reac,  ONLY: neq
        REAL(dp) :: Rates(neq)
        INTEGER,  ALLOCATABLE  :: perm(:)
        INTEGER,  ALLOCATABLE  :: NewSet(:), NewReacs(:)
        INTEGER :: i 
        REAL(dp) :: sum_Rates

        CALL SortVecAsc_R(Rates,perm)
        i = 0
        sum_Rates = 0.0_dp

        DO
          IF ( sum_Rates < eps_red ) THEN
            i = i + 1
            !sum_Rates = sum_Rates + Rates(i)     
            sum_Rates = sum_Rates + ABS(Rates(i))/SUM(ABS(Rates))
            !write(*,*) Rates(i),Rates(i)/SUM(Rates)
            perm(i)   = -1
          ELSE
            EXIT
          END IF
        END DO

        i = 1
        ALLOCATE(NewReacs(0))
        DO 
          IF ( perm(i) /= -1 ) THEN
            NewReacs = [NewReacs , perm(i)]
          END IF
          IF (i == neq ) EXIT
          i = i + 1
        END DO
        NewSet = [NewSet , NewReacs]

      END SUBROUTINE Update_ReactionSet

      SUBROUTINE Progress_Step(j,k,l)
        INTEGER :: j,k,l
        ! print the progress bar.
        WRITE(*,'(A1,10X,3(A,I0,4X),$)') char(13),'Iteration ::  ',j,'nReactions = ',k,'nSpecies = ',l
        !WRITE(*,'(A1,A,I0,A,I0,A,$)') char(13),'       Step  ',j,'  of  ',k,'  processed.'
      END SUBROUTINE Progress_Step

  END SUBROUTINE FluxAnalysis
  
  

END MODULE issa
