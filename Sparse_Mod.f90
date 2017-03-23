!================================================================================!
!                                                                                !
!                         This module includes a collection                      !
!                    of sparse matrix calculations for chemical                  !
!                       reaction systems, main format CSR                        !
!                                                                                !
!================================================================================!
!                                                                                
MODULE Sparse_Mod
  ! Contains:
  !   - SYMBOLIC MATRIX * MATRIC
  !   - SYMBOLIC MATRIX + MATRIX
  !   - TRANSPOSE MATRIX
  !   - SPARSE MATRIX * MATRIC
  !   - SPARSE MATRIX +- MATRIC
  !   - SORT ALGORITHM FOR SYMBOLIC MATRIX*MATRIX CALC
  !   - SPARSE JACOBIAN MATRIX CALC
  !   - SPARSE MITER MATRIX CALC
  !   - CONVERT COMPRESSED ROW FORMAT TO ROW-INDEX, COL-INDEX
  !   - PRINT SPARSE MATRIX (compressed Row format)
  !   - PRINT SPARSE MATRIX (ROW-INDEX COLUMN-INDEX)
  !   - SPARSE MATRIX*VECTOR+VECTOR (SPARSE-MATRIX, DENSE VECTORS)
  ! 
  USE mo_unirnk
  USE mo_MPI
  USE mo_control
  USE Kind_Mod
  USE mo_reac,          ONLY: combustion
  !
  IMPLICIT NONE 
  ! 
  INTEGER, PARAMETER, PRIVATE :: inilen=400
  !
  !
  TYPE CSR_Matrix_T            !Compressed Rowindex, standart columnindex
    INTEGER               :: m=0, n=0, nnz=0
    INTEGER, ALLOCATABLE  :: RowPtr(:)
    INTEGER, ALLOCATABLE  :: ColInd(:)
    INTEGER, ALLOCATABLE  :: DiagPtr(:)         ! position of diag entries
    INTEGER, ALLOCATABLE  :: DiagPtr_P(:)       ! permuted pos of diag entries
    INTEGER, ALLOCATABLE  :: DiagPtr_R(:)       ! permuted pos of rate entries
    INTEGER, ALLOCATABLE  :: DiagPtr_C(:)       ! permuted pos of conenctations
    INTEGER, ALLOCATABLE  :: RowVectorPtr(:)    ! for combustion matrix
    INTEGER, ALLOCATABLE  :: ColVectorPtr(:)    ! for combustion matrix
    INTEGER, ALLOCATABLE  :: Permu(:)           ! permutation vector (markowitz)
    INTEGER, ALLOCATABLE  :: InvPer(:)          ! invers permutation (inv markowitz)
    REAL(RealKind), ALLOCATABLE :: Val(:)
  END TYPE CSR_Matrix_T
  !
  TYPE SpRowIndColInd_T     !standart Rowindex, standart columnindex
    INTEGER :: m=0,n=0
    INTEGER, POINTER :: RowInd(:)
    INTEGER, POINTER :: ColInd(:)
    REAL(RealKind), POINTER :: Val(:)
  END TYPE SpRowIndColInd_T
  !
  TYPE SpRowColD_T
    INTEGER :: m=0,n=0
    INTEGER, POINTER :: RowPtr(:,:)=>NULL()
    INTEGER, POINTER :: ColInd(:)=>NULL()
    INTEGER, POINTER :: Permu(:)=>NULL()
    INTEGER, POINTER :: InvPer(:)=>NULL()
    INTEGER :: ep=1
    INTEGER :: last=0
    INTEGER :: len=0
    INTEGER :: nnz=0
  END TYPE SpRowColD_T
  !
  ! (/mincounts , anzahl gleicher mincounts/)
  INTEGER, ALLOCATABLE :: MarkowitzCounts(:,:)   
  !  
  ! global matrices containing chemical reaktion data (stoech.coefs)
  TYPE(CSR_Matrix_T) ::  A                & ! coef matrix of educts
  &                    , B                & ! coef matrix of products
  &                    , BA               & ! B-A
  &                    , BAT                ! Transpose(B-A)
  !
  !
  CONTAINS
  !
  !
  !
  SUBROUTINE New_CSR(newA,m,n,nnz)
    INTEGER, INTENT(IN) :: m, n
    INTEGER, INTENT(IN), OPTIONAL :: nnz
    !
    TYPE(CSR_Matrix_T) :: newA
    !
    newA%m=m
    newA%n=n
    !
    ALLOCATE(newA%RowPtr(m+1))
    newA%RowPtr=0
    newA%RowPtr(1)=1
    !
    IF (PRESENT(nnz)) THEN
      ALLOCATE(newA%ColInd(nnz))
      newA%ColInd=-1
      ALLOCATE(newA%Val(nnz))
      newA%Val=ZERO
      newA%nnz=nnz
    END IF
  END SUBROUTINE New_CSR
  !
  !
  SUBROUTINE Kill_Matrix_CSR(A)
    TYPE(CSR_Matrix_T) :: A
    !
    IF (ALLOCATED(A%RowPtr))      DEALLOCATE(A%RowPtr)
    IF (ALLOCATED(A%ColInd))      DEALLOCATE(A%ColInd)
    IF (ALLOCATED(A%DiagPtr))    DEALLOCATE(A%DiagPtr)
    IF (ALLOCATED(A%DiagPtr_R))  DEALLOCATE(A%DiagPtr_R)
    IF (ALLOCATED(A%DiagPtr_C))  DEALLOCATE(A%DiagPtr_C)
    IF (ALLOCATED(A%Permu))      DEALLOCATE(A%Permu)
    IF (ALLOCATED(A%InvPer))     DEALLOCATE(A%InvPer)  
    IF (ALLOCATED(A%Val))         DEALLOCATE(A%Val)
  END SUBROUTINE Kill_Matrix_CSR
  !
  !
  SUBROUTINE Free_SpRowColD(A)
    TYPE(SpRowColD_T) :: A
    !
    IF (ASSOCIATED(A%RowPtr))   NULLIFY(A%RowPtr)
    IF (ASSOCIATED(A%ColInd))   NULLIFY(A%ColInd)
    IF (ASSOCIATED(A%Permu ))   NULLIFY(A%Permu)
    IF (ASSOCIATED(A%InvPer))   NULLIFY(A%InvPer)
  END SUBROUTINE Free_SpRowColD
  !
  !
  SUBROUTINE  Kill_Matrix_SpRowIndColInd(A)
    TYPE(SpRowIndColInd_T) :: A
    !
    IF (ASSOCIATED(A%RowInd))   NULLIFY(A%RowInd)
    IF (ASSOCIATED(A%ColInd))   NULLIFY(A%ColInd)
    IF (ASSOCIATED(A%Val))      NULLIFY(A%Val)
  END SUBROUTINE Kill_Matrix_SpRowIndColInd   
  !
  !
  SUBROUTINE SparseID(Mat,dim)
    TYPE(CSR_Matrix_T), INTENT(OUT) :: Mat
    INTEGER, INTENT(IN) :: dim
    INTEGER :: i
    !
    CALL New_CSR(Mat,dim,dim,dim)
    DO i=1,dim
      Mat%RowPtr(i+1)=Mat%RowPtr(i)+1
      Mat%ColInd(i)=i
    END DO
  END SUBROUTINE SparseID
  !
  !
  SUBROUTINE RowColDToCSR(CSR,SpRow,m,n)
    TYPE(CSR_Matrix_T), INTENT(OUT) :: CSR
    TYPE(SpRowColD_T), INTENT(IN) :: SpRow
    INTEGER :: m, n
    !
    INTEGER :: i,j,jj,nzrA
    ! 
    CALL New_CSR(CSR,SpRow%n,SpRow%m,SpRow%nnz)
    
   
    IF ( EXTENDED ) THEN
      IF ( combustion ) THEN
        ALLOCATE(CSR%DiagPtr(n+m+1))
        ALLOCATE(CSR%DiagPtr_P(n+m+1))
      ELSE
        ALLOCATE(CSR%DiagPtr(n+m))
        ALLOCATE(CSR%DiagPtr_P(n+m))
      END IF
      ALLOCATE(CSR%DiagPtr_R(n))
      ALLOCATE(CSR%DiagPtr_C(m))
      CSR%DiagPtr_R = -11
      CSR%DiagPtr_C = -11
    ELSE
      ALLOCATE(CSR%DiagPtr(m))
      ALLOCATE(CSR%DiagPtr_P(m))
    END IF
    CSR%DiagPtr   = -11
    CSR%DiagPtr_P = -11
    
    NzrA=0
    DO i=1,CSR%m
      CSR%RowPtr(i+1)=CSR%RowPtr(i)
      DO jj=SpRow%RowPtr(1,i),SpRow%RowPtr(2,i)
        j=SpRow%ColInd(jj) 
        NzrA=NzrA+1
        CSR%ColInd(NzrA)=j
        IF (i==j) THEN
          CSR%DiagPtr(i)=CSR%RowPtr(i+1)
        END IF
        CSR%RowPtr(i+1)=CSR%RowPtr(i+1)+1
      END DO
    END DO
    
    IF (ASSOCIATED(SpRow%Permu)) THEN
      IF (.NOT.ALLOCATED(CSR%Permu)) THEN
        ALLOCATE(CSR%Permu(CSR%n))
        CSR%Permu(:)=-16
      END IF
      CSR%Permu(:)=SpRow%Permu(:)
    END IF
    IF (ASSOCIATED(SpRow%InvPer)) THEN
      IF (.NOT.ALLOCATED(CSR%InvPer)) THEN
        ALLOCATE(CSR%InvPer(CSR%n))
        CSR%InvPer(:)=-16
      END IF
      CSR%InvPer(:)=SpRow%InvPer(:)
    END IF

    ! permutate pointer to diagonal, rate and conc entries
    CSR%DiagPtr_P = CSR%DiagPtr( CSR%Permu )
    CSR%DiagPtr_R = CSR%DiagPtr_P(   1:n   )
    CSR%DiagPtr_C = CSR%DiagPtr_P( n+1:n+m )
    CSR%nnz = SpRow%nnz
  END SUBROUTINE RowColDToCSR
  !
  !
  SUBROUTINE CSRToSpRowColD(SpRowCol,CSR)
    TYPE(SpRowColD_T) :: SpRowCol
    TYPE(CSR_Matrix_T) :: CSR
    !
    INTEGER :: nzrCSR
    INTEGER :: AddLen
    INTEGER :: Start,End
    INTEGER :: i,j,jj
    !
    AddLen=10
    SpRowCol%m=CSR%m
    SpRowCol%n=CSR%n
    ALLOCATE(SpRowCol%RowPtr(2,SpRowCol%n))
    nzrCSR=SIZE(CSR%ColInd)
    SpRowCol%len=10*nzrCSR+AddLen*SpRowCol%n
    ALLOCATE(SpRowCol%ColInd(SpRowCol%len))
    SpRowCol%ColInd=0
    ALLOCATE(SpRowCol%Permu(SpRowCol%n))
    SpRowCol%Permu=0
    ALLOCATE(SpRowCol%InvPer(SpRowCol%n))
    SpRowCol%InvPer=0
    !
    Start=1
    DO i=1,CSR%n
      SpRowCol%RowPtr(1,i)=Start
      End=Start+(CSR%RowPtr(i+1)-CSR%RowPtr(i)-1)
      SpRowCol%RowPtr(2,i)=End
      DO jj=CSR%RowPtr(i),CSR%RowPtr(i+1)-1
        j=CSR%ColInd(jj)
        SpRowCol%ColInd(Start)=j
        Start=Start+1
      END DO  
      Start=Start+AddLen-1
    END DO  
    SpRowCol%ep=SpRowCol%RowPtr(2,SpRowCol%n)+1
    SpRowCol%last=SpRowCol%n  
  END SUBROUTINE CSRToSpRowColD
  !
  !
  SUBROUTINE Sort_SpRowColD(A)
    TYPE(SpRowColD_T) :: A
    INTEGER :: i,jj
    DO i=1,A%m
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        A%ColInd(jj)=A%Permu(A%ColInd(jj))
      END DO
      CALL SortVec(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
    END DO
  END SUBROUTINE Sort_SpRowColD
  !
  !
  SUBROUTINE SymbLU_SpRowColD(A,Permu)
    TYPE(SpRowColD_T) :: A
    INTEGER :: Permu(:)
    !
    INTEGER :: RowPiv(A%n)
    INTEGER :: i,j,l,jj,ip,iPiv
    LOGICAL :: ins
    !
    DO i=1,A%n
      A%InvPer(i)=i
      A%Permu(i)=i
    END DO
    !
    DO i=1,A%n
      ip=A%Permu(Permu(i))
      CALL Swap(A%InvPer(i),A%InvPer(ip))
      CALL Swap(A%RowPtr(1,i),A%RowPtr(1,ip))
      CALL Swap(A%RowPtr(2,i),A%RowPtr(2,ip))
      A%Permu(A%InvPer(i))=i
      A%Permu(A%InvPer(ip))=ip
      IF (A%last==i) THEN
        A%last=ip
      ELSE IF (A%last==ip) THEN
        A%last=i
      END IF
      !
      !   Update
      iPiv=0
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        IF (A%Permu(A%ColInd(jj))>i) THEN
          iPiv=iPiv+1 
          RowPiv(iPiv)=A%ColInd(jj)
        END IF
      END DO
      IF (iPiv>0) THEN
        DO j=i+1,A%n
          DO jj=A%RowPtr(1,j),A%RowPtr(2,j)
            IF (A%Permu(A%ColInd(jj))==i) THEN
              DO l=1,iPiv
                CALL Insert_SpRowColD(A,j,RowPiv(l),ins)
              END DO
              EXIT
            END IF
          END DO
        END DO
      END IF
    END DO
    DO i=1,A%n
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        A%ColInd(jj)=A%Permu(A%ColInd(jj))
      END DO
      CALL SortVec(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
    END DO
  END SUBROUTINE SymbLU_SpRowColD
  !
  !
  SUBROUTINE SymbLU_SpRowColD_M(A)
    TYPE(SpRowColD_T) :: A
    !
    INTEGER :: r(A%n),c(A%n),RowPiv(A%n)
    INTEGER :: i,j,l,jj,ip,ip1(1),iPiv
    REAL(RealKind) :: md 
    LOGICAL :: ins
    !
    ! für testzwecke
    ALLOCATE(MarkowitzCounts(A%n,2))
    MarkowitzCounts=0
    !
    c=0
    DO i=1,A%n
      A%InvPer(i)=i
      A%Permu(i)=i
      ! Compute initial Markowitz count
      r(i)=A%RowPtr(2,i)-A%RowPtr(1,i)+1
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        c(A%ColInd(jj))=c(A%ColInd(jj))+1
      END DO
    END DO
    !
    ! MAIN LOOP
    DO i=1,A%n
      ip=0
      md=1.d99
      DO j=i,A%n
        IF (((r(j)-1)*(c(j)-1)<=md)) THEN
          md=(r(j)-1)*(c(j)-1)
          ip=j
          !
          ! neu test
          MarkowitzCounts(i,1)=md
          IF ( (r(j)-1)*(c(j)-1) == md )  MarkowitzCounts(i,2)=MarkowitzCounts(i,2)+1
          !
        END IF
      END DO
      !print*, 'markocnts :::   ',i, MarkowitzCounts(i,1), MarkowitzCounts(i,2)
      ! 
      ! wird nie erreicht?
      IF (ip==0) THEN
        ip1(:)=MINLOC((r(i:A%n)-1)*(c(i:A%n)-1))+(i-1)
        ip=ip1(1)
        MarkowitzCounts(i,1)=ip1(1)
        MarkowitzCounts(i,2)=SIZE(ip1)
        print*, 'ip1(1)=',ip1(1), SIZE(ip1)
      END IF
      !
      CALL Swap(A%InvPer(i),A%InvPer(ip))
      CALL Swap(r(i),r(ip))
      CALL Swap(c(i),c(ip))
      CALL Swap(A%RowPtr(1,i),A%RowPtr(1,ip))
      CALL Swap(A%RowPtr(2,i),A%RowPtr(2,ip))
      A%Permu(A%InvPer(i))=i
      A%Permu(A%InvPer(ip))=ip
      IF (A%last==i) THEN
        A%last=ip
      ELSE IF (A%last==ip) THEN
        A%last=i
      END IF
      !
      ! Update
      iPiv=0
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        IF (A%Permu(A%ColInd(jj))>i) THEN
          iPiv=iPiv+1 
          RowPiv(iPiv)=A%ColInd(jj)
          c(A%Permu(A%ColInd(jj)))=c(A%Permu(A%ColInd(jj)))-1 
        END IF
      END DO
      IF (iPiv>0) THEN
        DO j=i+1,A%n
          DO jj=A%RowPtr(1,j),A%RowPtr(2,j)
            IF (A%Permu(A%ColInd(jj))==i) THEN
              r(j)=r(j)-1 
              c(i)=c(i)-1
              DO l=1,iPiv
                CALL Insert_SpRowColD(A,j,RowPiv(l),ins)
                IF (ins) THEN
                  c(A%Permu(RowPiv(l)))=c(A%Permu(RowPiv(l)))+1
                  r(j)=r(j)+1
                END IF
              END DO
              EXIT
            END IF
          END DO
          IF (c(i)==1) EXIT
        END DO
      END IF
    END DO
    DO i=1,A%n
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        A%ColInd(jj)=A%Permu(A%ColInd(jj))
      END DO
      CALL SortVec(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
      A%nnz=A%nnz+SIZE(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
    END DO
  END SUBROUTINE SymbLU_SpRowColD_M
  !
  !
  SUBROUTINE Swap(i,j)
    INTEGER :: i,j
    !
    INTEGER :: iTemp
    !  
    iTemp=i
    i=j
    j=iTemp
  END SUBROUTINE Swap
  !
  !
  SUBROUTINE SortVec(vec)
    INTEGER :: vec(:)
    !
    INTEGER :: i,itemp,j,n
    n=SIZE(Vec)
    DO i=1,n
      DO j=1,n-i
        IF (vec(j)>vec(j+1)) THEN
          itemp=vec(j)
          vec(j)=vec(j+1)
          vec(j+1)=itemp
        END IF
      END DO
    END DO
  END SUBROUTINE SortVec
  !
  !
  SUBROUTINE Insert_SpRowColD(A,iA,jA,ins)
    TYPE(SpRowColD_T) :: A
    INTEGER :: iA,jA
    LOGICAL :: ins
    !
    INTEGER :: itemp,j,l
    INTEGER, ALLOCATABLE :: iWork(:)
    !
    ! Test ob Element (ia,ja) bereits enthalten
    IF (iA==0.OR.jA==0) THEN
      WRITE(*,*) 'iA',iA
      WRITE(*,*) 'jA',jA
      CALL FinishMPI()
      STOP 'STOP'
    END IF  
    ins=.TRUE.
    DO j=A%RowPtr(1,iA),A%RowPtr(2,iA)
      IF (jA==A%ColInd(j)) THEN
        ins=.FALSE.
      END IF
    END DO
    !
    ! Test auf freien Speicherplatz in der ia-ten
    ! Zeile von a
    ! 
    IF (ins) THEN
      IF (A%ColInd(A%RowPtr(2,iA)+1)/=0) THEN
        ! ja-te Zeile von a wird nach hinten
        itemp=A%ep
        DO l=A%RowPtr(1,iA),A%RowPtr(2,iA)
          A%ColInd(A%ep)=A%ColInd(l)
          A%ColInd(l)=0
          A%ep=A%ep+1
        END DO
        A%RowPtr(2,iA)=A%ep-1
        A%RowPtr(1,iA)=itemp
        ! A%ep=A%ep+1
        A%last=iA
      ENDIF
      A%RowPtr(2,iA)=A%RowPtr(2,iA)+1
      A%ColInd(A%RowPtr(2,iA))=jA
      IF (iA==A%last) A%ep=A%ep+1
    END IF
    !
    IF (A%ep>=A%len-A%m) THEN
      CALL gcmat_SpRowColD(A)
    END IF
    !
    IF (A%ep>=A%len-A%m) THEN
      !   Speicherplatz von A nicht ausreichend
      ALLOCATE(iWork(A%ep))
      iWork(1:A%ep)=A%ColInd(1:A%ep)
      DEALLOCATE(A%ColInd)
      A%len=2*A%len
      ALLOCATE(A%ColInd(A%len))
      A%ColInd(1:A%ep)=iWork(1:A%ep)
      DEALLOCATE(iWork)  
    END IF
  END SUBROUTINE Insert_SpRowColD
  !
  ! 
  SUBROUTINE gcmat_SpRowColD(A)
   !   Externe Variable
   TYPE (SpRowColD_T) :: A
   ! 
   !   gcmat komprimiert eine zeilenorientierte, dynamische
   !   Speicherstruktur einer schwachbesetzten Matrix a
   !   der Ordnung n. Die Spaltenindizes der Nichtnull-
   !   elemente der i-ten Zeile von a sind in A%ColInd(A%RowPtr(i,1)),
   !   A%ColInd(A%RowPtr(1,i)+1)...,A%ColInd(A%RowPtr(2,i)) ent-
   !   halten.
   ! 
   !    Beschreibung der Parameter
   ! 
   ! 
   !    A%n      (i/o) integer
   !                 Dimension of matrix a
   !
   !    A%RowPtr (i/o) integer(2,n)
   !                 Pointerfeld zur Beschreibung von a. A%RowPtr muss
   !                 durch das rufende Programm belegt werden.
   ! 
   !    A%ColInd (i/o) integer(len)
   !   Externe Variable
   !   TYPE (SpRowColD_T) :: A
   ! 
   !   gcmat komprimiert eine zeilenorientierte, dynamische
   !   Speicherstruktur einer schwachbesetzten Matrix a
   !   der Ordnung n. Die Spaltenindizes der Nichtnull-
   !   elemente der i-ten Zeile von a sind in A%ColInd(A%RowPtr(i,1)),
   !   A%ColInd(A%RowPtr(1,i)+1)...,A%ColInd(A%RowPtr(2,i)) ent-
   !   halten.
   ! 
   !    Beschreibung der Parameter
   ! 
   ! 
   !    A%n      (i/o) integer
   !                 Dimension of matrix a
   !
   !    A%RowPtr (i/o) integer(2,n)
   !                 Pointerfeld zur Beschreibung von a. A%RowPtr muss
   !                 durch das rufende Programm belegt werden.
   ! 
   !    A%ColInd (i/o) integer(len)
   !                 Feld zur dynamischen Verwaltung der Spaltenindizes
   !                 der Nichtnullelemente von a. 
   ! 
   !    A%ep     (i/o) integer
   !                 Pointer der auf das erste freie Feld in a verweist,
   !                 d.h A%ColInd(ep),...,A%ColInd(len) sind frei verfuegbar.
   ! 
   !    A%len    (i)   integer
   !                 Gesamtlaenge des Feldes A%ColInd. 
   ! 
   !  Interne Variable
   !
   INTEGER i,iz,j,l,pointr,rowlen,ep,len,m
   !
   m=A%m
   ep=A%ep
   len=A%len
   pointr=1
   i=1
   !
   DO 
     IF (i>=ep) EXIT
       IF (A%ColInd(i).ne.0) THEN
         !
         ! Ermittlung der aktuellen Zeile sowie deren Laenge
         DO l=1,m
           IF (A%RowPtr(1,l).le.i.and.i.le.A%RowPtr(2,l)) THEN
             iz=l
           END IF
         END DO
         rowlen=A%RowPtr(2,iz)-A%RowPtr(1,iz)
         ! 
         ! Setzen der neuen Anfangs- und Endadresse der
         ! aktuellen Zeile
         A%RowPtr(1,iz)=pointr
         A%RowPtr(2,iz)=pointr+rowlen
         DO j=pointr,pointr+rowlen
           A%ColInd(j)=A%ColInd(i)
           i=i+1
         END DO
         i=i-1
         pointr=A%RowPtr(2,iz)+1
       ENDIF
      i=i+1
    END DO
    !
    !   Belegung des freien Teils von A%ColInd mit 0
    !
    ep=pointr
    DO i=1,m
      IF (A%RowPtr(1,i).gt.ep) THEN
        A%RowPtr(1,i)=ep
        A%RowPtr(2,i)=A%RowPtr(1,i)-1
        ep=ep+inilen
      END IF
    END DO
    !        
    DO i=pointr,len
      A%ColInd(i)=0
    END DO
    A%ep=ep
  END SUBROUTINE gcmat_SpRowColD
  !
  !
  SUBROUTINE SpDeallocate_SpRowColD(A)
    TYPE(SpRowColD_T) :: A
    !
    A%m=0
    A%n=0
    IF (ASSOCIATED(A%RowPtr)) DEALLOCATE(A%RowPtr)
    IF (ASSOCIATED(A%ColInd)) DEALLOCATE(A%ColInd)
    IF (ASSOCIATED(A%Permu))  DEALLOCATE(A%Permu)
    IF (ASSOCIATED(A%InvPer)) DEALLOCATE(A%InvPer)
  END SUBROUTINE SpDeallocate_SpRowColD
  !
  !
  SUBROUTINE SymbTransposeSparse(A,B)
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    TYPE(CSR_Matrix_T), INTENT(OUT) :: B
    !
    INTEGER :: i,j                                 ! zählvariabl schleIFen
    INTEGER :: indx
    !
    B%m=A%n
    B%n=A%m
    !
    ALLOCATE(B%RowPtr(B%m+1))
    B%RowPtr=0
    B%RowPtr(1)=1
    !
    DO i=1,A%m
      DO j=A%RowPtr(i),A%RowPtr(i+1)-1
        B%RowPtr(A%ColInd(j)+1)=B%RowPtr(A%ColInd(j)+1)+1
      END DO
    END DO
    !
    DO i=1,B%m
      B%RowPtr(i+1)=B%RowPtr(i)+B%RowPtr(i+1)
    END DO
    !
    ALLOCATE(B%ColInd(SIZE(A%ColInd)))
    B%ColInd=0
    DO i=1,A%m
      DO j=A%RowPtr(i),A%RowPtr(i+1)-1
        indx=A%ColInd(j)
        B%ColInd(B%RowPtr(indx))=i
        B%RowPtr(indx)=B%RowPtr(indx)+1
      END DO
    END DO
    DO i=B%m,1,-1
      B%RowPtr(i+1)=B%RowPtr(i)
    END DO
    B%RowPtr(1)=1
    B%nnz=B%RowPtr(B%m+1)-1
  END SUBROUTINE SymbTransposeSparse
  !
  !
  ! Transpose Matrix
  SUBROUTINE TransposeSparse(MatAT,MatA)
    TYPE(CSR_Matrix_T), INTENT(OUT) :: MatAT
    TYPE(CSR_Matrix_T), INTENT(IN)  :: MatA
    !
    INTEGER :: i,j                                 ! zählvariabl schleIFen
    INTEGER :: indx
    !
    !
    CALL New_CSR(MatAT,MatA%n,MatA%m)
    !
    !
    DO i=1,MatA%m
      DO j=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        MatAT%RowPtr(MatA%ColInd(j)+1)=MatAT%RowPtr(MatA%ColInd(j)+1)+1
      END DO
    END DO
    !
    DO i=1,MatAT%m
      !print*, 'RowPtr=',MatAT%RowPtr(i)
      MatAT%RowPtr(i+1)=MatAT%RowPtr(i)+MatAT%RowPtr(i+1)
    END DO
      !print*, 'RowPtr=',MatAT%RowPtr(MatAT%m+1)
    !
    ALLOCATE(MatAT%ColInd(SIZE(MatA%ColInd)))
    ALLOCATE(MatAT%Val(SIZE(MatA%Val)))
    MatAT%ColInd=0
    MatAT%Val=ZERO
    DO i=1,MatA%m
      DO j=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        indx=MatA%ColInd(j)
        MatAT%ColInd(MatAT%RowPtr(indx))=i
        MatAT%Val(MatAT%RowPtr(indx))=MatA%Val(j)
        MatAT%RowPtr(indx)=MatAT%RowPtr(indx)+1
      END DO
    END DO
    DO i=MatAT%m,1,-1
      MatAT%RowPtr(i+1)=MatAT%RowPtr(i)
    END DO
    MatAT%RowPtr(1)=1
    MatAT%nnz=MatA%RowPtr(MatA%m+1)-1
  END SUBROUTINE TransposeSparse
  !
  !
  ! SYMBOLIC MATRIX * MATRIX
  SUBROUTINE SymbolicMult(A,B,C)
    ! A*B=C
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    TYPE(CSR_Matrix_T), INTENT(IN) :: B
    TYPE(CSR_Matrix_T), INTENT(OUT) :: C
    !
    INTEGER ::  indx(MAX(A%n,A%m,B%n))
    INTEGER :: i, j, jj, k
    INTEGER :: istart, length, iTemp
    !
      !symbolic matrix multiply c=a*b
    CALL New_CSR(C,A%m,B%n)
    !
    !main loop
    indx=0
    DO i=1,A%m
      iStart=-1
      Length=0
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        DO k=B%RowPtr(j),B%RowPtr(j+1)-1
          IF (indx(B%ColInd(k))==0) THEN
            indx(B%ColInd(k))=istart
            istart=B%ColInd(k)
            length=length+1
          END IF
        END DO
      END DO
      C%RowPtr(i+1)=C%RowPtr(i)+Length
      !
      DO j=C%RowPtr(i),C%RowPtr(i+1)-1
        iTemp=iStart
        istart=indx(istart)
        indx(iTemp)=0
      END DO
      indx(i) = 0
    END DO
    !==========================================
    ALLOCATE(C%ColInd(C%RowPtr(C%m+1)-1))
    C%ColInd=0
    !
    DO i=1,A%m
      iStart=-1
      Length=0
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        DO k=B%RowPtr(j),B%RowPtr(j+1)-1
          IF (indx(B%ColInd(k))==0) THEN
            indx(B%ColInd(k))=istart
            istart=B%ColInd(k)
            length=length+1
          END IF
        END DO
      END DO
      C%RowPtr(i+1)=C%RowPtr(i)+length
      DO j=C%RowPtr(i),C%RowPtr(i+1)-1
       C%ColInd(j)=istart
       istart=indx(istart)
       indx(C%ColInd(j))=0
      END DO
      indx(i) = 0
    END DO
    !
    DO i=1,C%m
      CALL Sort(C%Colind(C%RowPtr(i):C%RowPtr(i+1)-1))
    END DO
    ALLOCATE(C%Val(C%RowPtr(C%m+1)-1))
    C%Val=ZERO
    C%nnz=C%RowPtr(C%m+1)-1
  END SUBROUTINE SymbolicMult
  !
  !
  ! SPARSE MATRIX*MATRIX CALC
  SUBROUTINE SparseMult(A,B,C)
    ! A * B = C
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: A
    TYPE(CSR_Matrix_T), INTENT(IN) :: B
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: C
    !
    INTEGER :: i,j,jj,k,kk
    REAL(RealKind) :: ajj
    REAL(RealKind) :: temp(MAX(A%m,A%n,B%n))
    temp=ZERO
    !
    DO i=1,A%m
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        ajj=A%Val(jj)
        DO kk=B%RowPtr(j),B%RowPtr(j+1)-1
          k=B%ColInd(kk)
          temp(k)=temp(k)+ajj*B%Val(kk)
        END DO
      END DO
      DO j=C%RowPtr(i),C%RowPtr(i+1)-1
        C%Val(j)=temp(C%ColInd(j))
        temp(C%ColInd(j))=0.
      END DO
    END DO
  END SUBROUTINE SparseMult
  !
  ! SYMBOLIC MATRIX + MATRIX  -->  MatC = MatA + MatB
  SUBROUTINE SymbolicAdd(MatC,MatA,MatB)
    TYPE(CSR_Matrix_T), INTENT(IN)  :: MatA
    TYPE(CSR_Matrix_T), INTENT(IN)  :: MatB
    TYPE(CSR_Matrix_T), INTENT(OUT) :: MatC
    !
    INTEGER :: i,ii,j,jj,k,kk
    INTEGER, ALLOCATABLE :: TmpCol(:)
    INTEGER, ALLOCATABLE :: PermVec(:)
    INTEGER, ALLOCATABLE :: Indizes(:)
    INTEGER :: ColLen
    INTEGER :: currentlength
    INTEGER :: sameCnt
    !
    IF (.NOT.((MatA%m==MatB%m).AND.(MatA%n==MatB%n))) THEN
      WRITE(*,*) 'Wrong Matrix Dim'
      WRITE(*,*) 'A: ',MatA%m,MatA%n
      WRITE(*,*) 'B: ',MatB%m,MatB%n
      CALL FinishMPI()
      STOP 'STOP'
    END IF
    !
    CALL New_CSR(MatC,MatA%m,MatA%n)
    !
    DO i=1,MatC%m
      currentlength=(MatA%RowPtr(i+1)-MatA%RowPtr(i))+(MatB%RowPtr(i+1)-MatB%RowPtr(i))
      ALLOCATE(Indizes(currentlength))
      Indizes=0
      sameCnt=0
      DO ii=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        DO j=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
          IF (MatA%ColInd(ii)==MatB%ColInd(j)) THEN
            sameCnt=sameCnt+1
            Indizes(sameCnt)=j
          END IF
        END DO
      END DO
      MatC%RowPtr(i+1)=MatC%RowPtr(i)+currentlength-sameCnt
      DEALLOCATE(Indizes)
    END DO
    !
    ! Allocate ColInd
    ALLOCATE(MatC%ColInd(MatC%RowPtr(MatC%m+1)-1))
    MatC%ColInd=0
    !
    kk=1
    DO i=1,MatC%m
      k=1
      currentlength=MatC%RowPtr(i+1)-MatC%RowPtr(i)
      sameCnt=0
      DO ii=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        DO j=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
          IF (MatA%ColInd(ii)==MatB%ColInd(j)) THEN
            sameCnt=sameCnt+1
          END IF
        END DO
      END DO
      !
      ALLOCATE(TmpCol(currentlength+sameCnt))
      TmpCol=0
      ALLOCATE(PermVec(currentlength+sameCnt))
      PermVec=0
      !
      DO jj=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        TmpCol(k)=MatA%ColInd(jj)
        k=k+1
      END DO
      DO jj=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
        TmpCol(k)=MatB%ColInd(jj)
        k=k+1
      END DO
      !
      CALL unirnk(TmpCol,PermVec,ColLen)
      DO j=1,ColLen
        MatC%ColInd(kk)=TmpCol(PermVec(j))
        kk=kk+1
      END DO
      !
      DEALLOCATE(TmpCol)
      DEALLOCATE(PermVec)
    END DO
    !
    ALLOCATE(MatC%Val(MatC%RowPtr(MatC%m+1)-1))
    MatC%Val=ZERO
    MatC%nnz=MatC%RowPtr(MatC%m+1)-1
  END SUBROUTINE SymbolicAdd
  !
  !
  ! SPARSE  MATRIX + MATRIX  -->  MatC = MatA (Sub_Add) MatB 
  SUBROUTINE SparseAdd(MatC,MatA,MatB,Sub_Add)
    TYPE(CSR_Matrix_T), INTENT(IN)    :: MatA, MatB
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: MatC
    CHARACTER, INTENT(IN) :: Sub_Add
    !
    ! Temp variables
    INTEGER, ALLOCATABLE :: tmpCol(:)
    INTEGER, ALLOCATABLE :: perm(:)
    REAL(RealKind), ALLOCATABLE :: tmpBVal(:)
    REAL(RealKind), ALLOCATABLE :: tmpVal(:)
    INTEGER :: lenColIndA, lenColIndB
    INTEGER :: i
    !
    IF ((Sub_Add=='-').OR.(Sub_Add=='+')) THEN
      IF (Sub_Add=='-') THEN
        ALLOCATE(tmpBVal(MatB%RowPtr(MatB%m+1)-1))
        tmpBVal=-MatB%Val
      ELSE
        tmpBVal=MatB%Val
      END IF
    ELSE
      WRITE(*,*) 'Sub_Add should be "-" or "+" !!'
      CALL FinishMPI()
      STOP 'STOP'
    END IF
    !
    DO i=1,MatC%m
      lenColIndA=MatA%RowPtr(i+1)-MatA%RowPtr(i)
      lenColIndB=MatB%RowPtr(i+1)-MatB%RowPtr(i)
      ALLOCATE(tmpCol(lenColIndA+lenColIndB))
      tmpCol=0
      ALLOCATE(tmpVal(lenColIndA+lenColIndB))
      tmpVal=ZERO
      ALLOCATE(Perm(lenColIndA+lenColIndB))
      Perm=0
      tmpCol=(/MatA%ColInd(MatA%RowPtr(i):MatA%RowPtr(i+1)-1), &
      &        MatB%ColInd(MatB%RowPtr(i):MatB%RowPtr(i+1)-1) /)
      tmpVal=(/MatA%Val(MatA%RowPtr(i):MatA%RowPtr(i+1)-1),    &
      &        tmpBVal(MatB%RowPtr(i):MatB%RowPtr(i+1)-1) /)
      !
      CALL sort2(tmpCol,Perm)
      tmpVal=tmpVal(Perm)
      CALL CompressList(tmpCol,tmpVal)
      MatC%Val(MatC%RowPtr(i):MatC%RowPtr(i+1)-1)=tmpVal
      !
      DEALLOCATE(tmpCol)
      DEALLOCATE(tmpVal)
      DEALLOCATE(Perm)
    END DO
    DEALLOCATE(tmpBVal)
  END SUBROUTINE SparseAdd
  !
  !
  ! SORT ALGORITHM FOR SYMBOLIC MATRIX*MATRIX CALC
  SUBROUTINE Sort(v)
    INTEGER :: v(:)
    !
    INTEGER :: i,j,temp
    DO i=1,SIZE(v)
      DO j=i+1,SIZE(v)
        IF (v(i)>v(j)) THEN
          temp=v(i)
          v(i)=v(j)
          v(j)=temp
        END IF
      END DO
    END DO
  END SUBROUTINE Sort
  !
  !
  ! SORT ALGORITHM FOR MATRIX+MATRIX
  SUBROUTINE Sort2(v,perm)
    INTEGER :: v(:)
    INTEGER :: perm(:)
    !
    INTEGER :: i,j,temp1,temp2
    DO i=1,SIZE(v)
      perm(i)=i
    END DO
    !
    DO i=1,SIZE(v)
      DO j=i+1,SIZE(v)
        IF (v(i)>v(j)) THEN
          ! Feldelemente sortieren
          temp1=v(i)
          v(i)=v(j)
          v(j)=temp1
          ! Permutationsvektor erzeugen
          temp2=perm(i)
          perm(i)=perm(j)
          perm(j)=temp2
        END IF
      END DO
    END DO
  END SUBROUTINE Sort2
  !
  !
  SUBROUTINE CompressList(ColInd,Val)
    INTEGER, ALLOCATABLE :: ColInd(:)
    REAL(RealKind), ALLOCATABLE, OPTIONAL :: Val(:)
    !
    INTEGER :: i,j,iList,MemberCol
    INTEGER :: TempListCol(SIZE(ColInd))
    REAL(RealKind) :: MemberVal
    REAL(RealKind) :: TempListVal(SIZE(ColInd))
    LOGICAL :: Insert
    !
    TempListVal=ZERO
    iList=0
    !
    S1:DO i=1,SIZE(ColInd)
      MemberCol=ColInd(i)
      IF (PRESENT(Val)) MemberVal=Val(i)
      Insert=.TRUE.
      S2:DO j=1,iList
        IF (MemberCol==TempListCol(j)) THEN
          Insert=.FALSE.
          IF (PRESENT(Val)) TempListVal(iList)=TempListVal(iList)+MemberVal
          EXIT S2
        END IF
      END DO S2
      IF (Insert) THEN
        iList=iList+1
        TempListCol(iList)=MemberCol
        IF (PRESENT(Val)) TempListVal(iList)=MemberVal
      END IF
    END DO S1
    DEALLOCATE(ColInd)
    IF (PRESENT(Val)) DEALLOCATE(Val)
    ALLOCATE(ColInd(1:iList))
    IF (PRESENT(Val)) ALLOCATE(Val(1:iList))
    ColInd=TempListCol(1:iList)
    IF (PRESENT(Val)) Val=TempListVal(1:iList)
  END SUBROUTINE CompressList
  !
  !
  ! SPARSE JACOBIMATRIX CALC
  SUBROUTINE JacobiMatrix(gMat,aMat,rVec,yVec,jMat)
    !
    ! jMat = gMat*Dr*aMat*invDy;
    !
    TYPE(CSR_Matrix_T), INTENT(IN) :: gMat
    TYPE(CSR_Matrix_T), INTENT(IN) :: aMat
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: jMat
    REAL(RealKind), INTENT(IN) :: rVec(aMat%m)
    REAL(RealKind), INTENT(IN) :: yVec(aMat%n)
    !
    INTEGER :: i,j,jj,k,kk
    REAL(RealKind) :: ajj
    REAL(RealKind) :: temp(MAX(gMat%m,gMat%n,aMat%n))
    !
    temp=ZERO
    !jMat%m=aMat%n
    !jMat%n=jMat%m
    !
    jMat%Val=ZERO
    !print*, 'nID=', MPI_ID, 'SumJacColInd=',SUM(jMat%ColInd), 'sumf0=', sum(abs(rVec)),&
    !&           'sumY=',SUM(ABS(yVec)), 'sumAvals=', SUM(ABS(aMat%val)), 'sumBATvals=',SUM(ABS(gMat%val))
    !
    DO i=1,gMat%m
      DO jj=gMat%RowPtr(i),gMat%RowPtr(i+1)-1
        j=gMat%ColInd(jj)
        ajj=gMat%Val(jj)*rVec(j)
        DO kk=aMat%RowPtr(j),aMat%RowPtr(j+1)-1
          k=aMat%ColInd(kk)
          temp(k)=temp(k)+ajj*aMat%Val(kk)/yVec(k)
        END DO
      END DO
      DO j=jMat%RowPtr(i),jMat%RowPtr(i+1)-1
        jMat%Val(j)=temp(jMat%ColInd(j))
        temp(jMat%ColInd(j))=ZERO
      END DO
    END DO
  END SUBROUTINE JacobiMatrix
  !
  !
  ! SPARSE MITER CALCULATION_CLASSIC
  SUBROUTINE Miter_Classic(gMat,aMat,rVec,yVec,h,g,miter)
    !
    ! miter = Id - h*g*gMat*Dr*aMat*invDy;
    !
    TYPE(CSR_Matrix_T), INTENT(IN) :: gMat
    TYPE(CSR_Matrix_T), INTENT(IN) :: aMat
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: miter
    REAL(RealKind), INTENT(IN) :: rVec(:)
    REAL(RealKind), INTENT(IN) :: yVec(:)
    REAL(RealKind), INTENT(IN) :: h, g
    !
    INTEGER :: i,ii,j,jj,k,kk
    REAL(RealKind) :: ajj
    REAL(RealKind) :: hg
    REAL(RealKind) :: temp(MAX(gMat%m,gMat%n,aMat%n))
    !
    Miter%Val = ZERO
    temp      = ZERO
    hg        = h*g
    !print*, 'debug :: hg=', hg
    !
    DO i = 1 , gMat%m

      DO jj = gMat%RowPtr(i) , gMat%RowPtr(i+1)-1
        j   = gMat%ColInd(jj)
        ajj = gMat%Val(jj) * rVec(j) 
        
        ! calculate the derivative dY'/dY
        DO kk = aMat%RowPtr(j) , aMat%RowPtr(j+1)-1
          k       = aMat%ColInd(kk)
          temp(k) = temp(k) + ajj * aMat%Val(kk)/yVec(k)
        END DO
        
        IF ( combustion ) THEN
        ! calculate the drivative dY'/dT 

        ! calculate the drivative dT'/dY

        ! calculate the drivative dT'/dT

        END IF

      END DO

      ! update the matrix entrie
      ii = 1
      DO jj = miter%RowPtr(i) , miter%RowPtr(i+1)-1
        j = miter%ColInd(jj)
        IF ( i==j ) THEN
          miter%Val(jj) = ONE - (hg)*temp(j)
        ELSE
          miter%Val(jj) =     - (hg)*temp(j)
        END IF
        ii = ii + 1
        temp(j) = ZERO
      END DO
    END DO
  END SUBROUTINE Miter_Classic
  !
  !
  SUBROUTINE SymbolicExtendedMatrix(A,BAT,MiterOut)
    TYPE(CSR_Matrix_T), INTENT(IN) :: A, BAT
    TYPE(CSR_Matrix_T), INTENT(OUT) :: MiterOut
    TYPE(CSR_Matrix_T) :: Miter
    TYPE(CSR_Matrix_T) :: Id
    INTEGER :: LenColInd
    INTEGER, ALLOCATABLE :: tmpColInd(:)
    !
    INTEGER :: i,j
    !
    CALL New_CSR(Miter,A%m+BAT%m,A%n+BAT%n)
    !
    DO i=1,A%m
      Miter%RowPtr(i+1)=Miter%RowPtr(i)+(A%RowPtr(i+1)-A%RowPtr(i))
    END DO
    DO i=1,BAT%m
      Miter%RowPtr(A%m+i+1)=Miter%RowPtr(A%m+i)+(BAT%RowPtr(i+1)-BAT%RowPtr(i))
    END DO
    !
    ALLOCATE(Miter%ColInd(Miter%RowPtr(Miter%m+1)-1))
    Miter%ColInd=0
    j=1
    DO i=1,A%m
      LenColInd=A%RowPtr(i+1)-A%RowPtr(i)
      ALLOCATE(tmpColInd(LenColInd))
      tmpColInd=(/A%m+(A%ColInd(A%RowPtr(i):(A%RowPtr(i+1)-1)))/)
      CALL sort(tmpColInd)
      Miter%ColInd(Miter%RowPtr(i):Miter%RowPtr(i+1)-1)=tmpColInd
      DEALLOCATE(tmpColInd)
    END DO
    !
    DO i=1,BAT%m
      LenColInd=BAT%RowPtr(i+1)-BAT%RowPtr(i)
      ALLOCATE(tmpColInd(LenColInd))
      tmpColInd=(/BAT%ColInd(BAT%RowPtr(i):(BAT%RowPtr(i+1)-1))/)
      CALL sort(tmpColInd)
      Miter%ColInd(Miter%RowPtr(A%m+i):Miter%RowPtr(A%m+i+1)-1)=tmpColInd
      DEALLOCATE(tmpColInd)
    END DO
    !
    ALLOCATE(Miter%Val(Miter%RowPtr(Miter%m+1)-1))
    Miter%Val=ZERO
    !
    ! build sparse unit matrix
    CALL SparseID(Id,A%m+BAT%m)
    !
    DO i=1,Id%m
      Id%RowPtr(i+1)=Id%RowPtr(i)+1
      Id%ColInd(i)=i
    END DO
    CALL SymbolicAdd(MiterOut,Id,Miter)
  END SUBROUTINE SymbolicExtendedMatrix
  !
  !
  SUBROUTINE BuildSymbolicClassicMatrix(CL,Jac,RowGamma)
    TYPE(CSR_Matrix_T), INTENT(OUT) :: CL
    TYPE(CSR_Matrix_T), INTENT(IN)  :: Jac
    REAL(RealKind),     INTENT(IN)  :: RowGamma

    TYPE(CSR_Matrix_T) :: Id
    TYPE(CSR_Matrix_T) :: CL0
    INTEGER, ALLOCATABLE :: tmpColInd(:)
    !
    INTEGER :: i, j, jj, m, n, nnz
    !
    !------------------------------------------------------------------------------
    ! --- Set Matrix dimensions and nonzeros 
    !------------------------------------------------------------------------------
    !
    !
    IF ( combustion ) THEN
      m    = Jac%n + 1               ! nummber of rows
      n    = Jac%n + 1               ! number of coloumns
      nnz  = Jac%RowPtr(Jac%m+1)-1      &  ! nonzeros of Jacobian
      &       + Jac%m + Jac%n + 1       ! Dc,U^T and Dr,~K and X (down right)
    ELSE
      m    = Jac%n               ! nummber of rows
      n    = Jac%n               ! number of coloumns
      nnz  = Jac%RowPtr(Jac%m+1)-1        ! nonzeros of Jacobian
    END IF
    CALL New_CSR ( CL0 , m , n , nnz )
    CALL SparseID( Id  , n )
    ALLOCATE(CL0%DiagPtr(m))
    CL0%DiagPtr  = -12

    IF ( combustion ) THEN
      !
      DO i = 1 , m - 1
        CL0%RowPtr(i+1) = CL0%RowPtr(i) + (Jac%RowPtr(i+1)-Jac%RowPtr(i)) + 1

        CL0%ColInd( CL0%RowPtr(i):CL0%RowPtr(i+1)-1 ) =                           & 
        &    (/ Jac%ColInd(Jac%RowPtr(i):(Jac%RowPtr(i+1)-1)) , n /)
      END DO
      
      ! set last row (full row)
      CL0%RowPtr(m+1) = CL0%RowPtr(m) + n
      CL0%ColInd(CL0%RowPtr(m):CL0%RowPtr(m+1)-1) = [( i , i = 1 , n )]

      CALL SymbolicAdd(CL,Id,CL0)
      ALLOCATE(CL%RowVectorPtr(m-1))
      ALLOCATE(CL%ColVectorPtr(m-1))
      CL%RowVectorPtr = -ONE
      CL%ColVectorPtr = -ONE
      CL%Val          = ONE
    ELSE
      CALL SymbolicAdd(CL,Id,Jac)
    END IF
    !
    ALLOCATE(CL%DiagPtr(m))
    CL%DiagPtr  = -12
    ! get diagonal pointer
    DO i = 1 , n
      IF ( combustion .AND. i < n ) CL%ColVectorPtr(i)  = CL%RowPtr(i+1) - 1
      DO jj = CL%RowPtr(i) , CL%RowPtr(i+1)-1 
        j = CL%ColInd(jj)
        IF ( i == j ) CL%DiagPtr(i) = jj
      END DO
    END DO
    IF ( combustion ) CL%RowVectorPtr = [( i , i = 1 , n )]
  END SUBROUTINE BuildSymbolicClassicMatrix
  !
  !
  SUBROUTINE BuildSymbolicExtendedMatrix(EX,A,BAT,g)
    TYPE(CSR_Matrix_T), INTENT(OUT) :: EX
    TYPE(CSR_Matrix_T), INTENT(IN) :: A, BAT
    REAL(RealKind),     INTENT(IN) :: g

    INTEGER, ALLOCATABLE :: tmpColInd(:)
    !
    INTEGER :: i, j, dimBig, nnzBig
    !
    !------------------------------------------------------------------------------
    ! --- Set big Matrix dimensions and nonzeros 
    !------------------------------------------------------------------------------
    !
    IF ( combustion ) THEN
      DimBig  = A%m   + BAT%m+1               ! nummber of rows
      nnzBig  = A%nnz + BAT%nnz      &        ! nonzeros of alpha and (beta-alpha)^T
      &       + 2*A%n + 2*BAT%n + 1           ! Dc,U^T and Dr,~K and X (down right)
    ELSE
      DimBig  = A%m   + BAT%m                 ! nummber of rows
      nnzBig  = A%nnz + BAT%nnz      &        ! nonzeros of alpha and (beta-alpha)^T
      &       + A%n   + BAT%n                 ! Dr, Dc
    END IF
    
    CALL New_CSR( EX , DimBig , DimBig , nnzBig )
    
    ALLOCATE(EX%DiagPtr(DimBig))          ! entire diagonal
    ALLOCATE(EX%DiagPtr_R(A%m))         ! reaction rates
    ALLOCATE(EX%DiagPtr_C(BAT%m))       ! species concentrations
    EX%DiagPtr        = -13
    EX%DiagPtr_R      = -13
    EX%DiagPtr_C      = -13

    ! Allocate a full row vector and a full column vector
    IF ( combustion ) THEN
      ALLOCATE(EX%ColVectorPtr(A%m))    ! | Vector
      ALLOCATE(EX%RowVectorPtr(BAT%m))  ! _ Vector
      EX%ColVectorPtr = -13
      EX%RowVectorPtr = -13
    END IF
    !
    !---------------------------------------------------------------------------------
    ! --- Set Row Pointer, Coloum Index and 'Values'(=1.0d0 for Dr, Dc, ~K, U^T and X)
    !---------------------------------------------------------------------------------
    !
    !                    FIRST PART
    !          _                                  _ 
    !         |    Diag_1    |   g*alpha    | g*~K |
    !         |--------------+--------------+------|
    ! miter = |              |              |      |
    !         !--------------+--------------+------|
    !         |_             |              |     _|
    !
    IF ( combustion ) THEN

      DO i = 1 , A%m 

        EX%RowPtr(i+1)  = EX%RowPtr(i) + A%RowPtr(i+1) - A%RowPtr(i) + 2

        EX%ColInd(EX%RowPtr(i):EX%RowPtr(i+1)-1) =                                 & 
        &       (/  i  , A%m + A%ColInd(A%RowPtr(i):(A%RowPtr(i+1)-1))    , DimBig  /)

        EX%Val(EX%RowPtr(i):EX%RowPtr(i+1)-1)    =                                 & 
        &       (/ ONE , g*A%Val(A%ColInd(A%RowPtr(i):(A%RowPtr(i+1)-1))) ,  ONE /)
     
        ! set pointers for better access 
        EX%DiagPtr(i)       = EX%RowPtr(i)
        EX%DiagPtr_R(i)     = EX%RowPtr(i)
        EX%ColVectorPtr(i)  = EX%RowPtr(i+1) - 1

      END DO

    ELSE

      DO i = 1 , A%m 

        EX%RowPtr(i+1)  = EX%RowPtr(i) + A%RowPtr(i+1) - A%RowPtr(i) + 1

        EX%ColInd(EX%RowPtr(i):EX%RowPtr(i+1)-1)  =                        & 
        &       (/  i  , A%m + A%ColInd(A%RowPtr(i):(A%RowPtr(i+1)-1))    /)

        EX%Val(EX%RowPtr(i):EX%RowPtr(i+1)-1)     =                        & 
        &       (/ ONE , g*A%Val(A%ColInd(A%RowPtr(i):(A%RowPtr(i+1)-1))) /)

        ! set pointers for better access 
        EX%DiagPtr(i)   = EX%RowPtr(i)
        EX%DiagPtr_R(i) = EX%RowPtr(i)

      END DO

    END IF
    !
    !                    Second PART
    !          _                                  _ 
    !         |   Diag_1_nR  |   g*alpha    | g*~K |
    !         |--------------+--------------+------|
    ! miter = |   BAT_Mat    |  Diag_1_nS   |      |
    !         !--------------+--------------+------|
    !         |_             |              |     _|
    ! same for combustion and atmopheric stuff
    DO i=1,BAT%m

      EX%RowPtr(A%m+i+1)  = EX%RowPtr(A%m+i) + BAT%RowPtr(i+1) - BAT%RowPtr(i) + 1

      EX%ColInd(EX%RowPtr(A%m+i):EX%RowPtr(A%m+i+1)-1)  =                 &
      &         (/ BAT%ColInd(BAT%RowPtr(i):(BAT%RowPtr(i+1)-1)) , A%m+i /)

      EX%Val(EX%RowPtr(A%m+i):EX%RowPtr(A%m+i+1)-1)     =                 &
      &         (/ BAT%Val(BAT%RowPtr(i):(BAT%RowPtr(i+1)-1))    ,  ONE  /)
     
      EX%DiagPtr(i+A%m)   = EX%RowPtr(A%m+i+1)-1
      EX%DiagPtr_C(i)     = EX%RowPtr(A%m+i+1)-1

    END DO
    !                    Third PART
    !          _                                  _ 
    !         |   Diag_1_nR  |   g*alpha    | g*~K |
    !         |--------------+--------------+------|
    ! miter = |   BAT_Mat    |  Diag_1_nS   |      |
    !         !--------------+--------------+------|
    !         |_             |    -U^TD_c   |   1 _|
    !
    IF ( combustion ) THEN

      EX%RowPtr(DimBig+1)=EX%RowPtr(DimBig)+A%n+1
      
      EX%ColInd(EX%RowPtr(DimBig):EX%RowPtr(DimBig+1)-1) = [( i , i = BAT%n+1 , BAT%n+A%n+1 )]
      EX%Val(EX%RowPtr(DimBig):EX%RowPtr(DimBig+1)-1)    = -99.d0
      EX%RowVectorPtr(:) = [( i , i = BAT%n+1 , BAT%n+A%n+1 )]
      EX%DiagPtr(DimBig) = EX%RowPtr(DimBig+1)-1

    END IF
    !call printsparse(EX , '*')
  !
  END SUBROUTINE BuildSymbolicExtendedMatrix
  !
  !
  !
  SUBROUTINE SetLUvaluesEX(LU,m,n,invrVec,hyVec,KVec,UVec,X,constantValues)
    !
    ! Set values to block matrix
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: LU  
    INTEGER :: m, n
    REAL(RealKind), INTENT(IN) :: invrVec(:), KVec(:)
    REAL(RealKind), INTENT(IN) :: hyVec(:), UVec(:)
    REAL(RealKind), INTENT(IN) :: constantValues(:)
    REAL(RealKind), INTENT(IN) :: X

    integer :: i, ip
    !
    !
    !LU%Val = ZERO
    LU%Val = constantValues
    !do i=1,size(constantValues)
    !   write(9999,*), 'Debug___    constant values = ', constantValues(i)
    !end do
    !stop
    !
    !DO i=1,n
    !  LU%Val( LU%DiagPtr_R(i) ) = invrVec(i)
    !END DO
    !DO i=1,m
    !  LU%Val( LU%DiagPtr_C(i) ) = hyVec(i)
    !END DO
    IF ( combustion ) THEN
      !          _                                  _ 
      !         |   D_r^(-1)   |   g*alpha    | g*~K |
      !         |--------------+--------------+------|
      ! miter = |(beta-alpha)^T|    1/h*D_c   |      |
      !         !--------------+--------------+------|
      !         |_             | -1/h*U^TD_c  |   X _|
      !
      ! HIER HIER HIER mit schleife probierten
      LU%Val( LU%DiagPtr_R      ) = invrVec
      LU%Val( LU%DiagPtr_C      ) = hyVec  
      LU%Val( LU%ColVectorPtr   ) = KVec
      LU%Val( LU%RowVectorPtr   ) = UVec
      LU%Val( LU%DiagPtr(m+n+1) ) = X
    ELSE
      !          _                            _
      !         |   D_r^(-1)   |   g*alpha    |
      !         |--------------+--------------|
      ! miter = |(beta-alpha)^T|    1/h*D_c   |
      !         |_             |             _|
      !
      !LU%Val( LU%DiagPtr_P(  1:n  ) )  = invrVec
      !LU%Val( LU%DiagPtr_P(n+1:n+m) )  = hyVec
      LU%Val( LU%DiagPtr_R )  = invrVec
      LU%Val( LU%DiagPtr_C )  = hyVec
      !DO i = 1 , n
      !  ip = LU%Permu(i)
      !  LU%Val( LU%DiagPtr(ip) ) = invrVec(i)
      !END DO
      !!
      !DO  i = 1 , m
      !  ip = LU%Permu(n+i)
      !  LU%Val( LU%DiagPtr(ip) ) = hyVec(i)
      !END DO
    END IF
  END SUBROUTINE SetLUvaluesEX 
  !
  !
  SUBROUTINE SetLUvaluesCL(LU,A,Permu,UVec,KVec)
    !
    ! Set values to block matrix
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: LU
    TYPE(CSR_Matrix_T), INTENT(IN)    :: A
    INTEGER,            INTENT(IN)    :: Permu(:)
    REAL(RealKind),     INTENT(IN)    :: UVec(:), KVec(:)
    !          _ 
    INTEGER :: i
    !
    LU%Val = ZERO
    !print*, ' debug :: permu = ', permu
    !stop

    !IF ( combustion ) THEN
    !  ! fehlt noch
    !ELSE
    !  LU%Val( Permu ) = A%Val
    !END IF
    !
    DO i=1,A%RowPtr(A%m+1)-1
      LU%Val(Permu(i))=A%Val(i)
    END DO
    !
  END SUBROUTINE SetLUvaluesCL 
  !
  ! Permutes the values in Miter and writes it to LU structur of Miter
  ! Permutation vector is generated in this routine
  SUBROUTINE Get_LU_Permutaion(Permutation,LU,A,m,n)
    ! INOUT
    TYPE(CSR_Matrix_T) :: LU
    ! IN
    TYPE(CSR_Matrix_T) :: A
    INTEGER            :: m, n, nnzA
    ! OUT
    INTEGER, ALLOCATABLE :: Permutation(:)
    ! TEMP
    INTEGER :: i,ip,j,jj,jp,jjp,jp1
 
    nnzA=A%RowPtr(A%m+1)-1
  
   
    IF (.NOT.ALLOCATED(Permutation)) ALLOCATE(Permutation(nnzA))
    Permutation = -14

    LU%Val = ZERO

    DO i = 1 , A%n  
      ip = LU%Permu(i)
      DO jj = A%RowPtr(i) , A%RowPtr(i+1)-1
        j  = A%ColInd(jj)
        jp = LU%Permu(A%ColInd(jj))
        DO jjp = LU%RowPtr(ip) , LU%RowPtr(ip+1)-1
          jp1 = LU%ColInd(jjP)
          IF ( jp1 == jp ) THEN
            LU%Val(jjP) = A%Val(jj)
            Permutation(jj) = jjP
          END IF  
        END DO  
      END DO  
    END DO

  END SUBROUTINE Get_LU_Permutaion
  !
  !
  SUBROUTINE Miter_Extended(CSR_Mat,m,n,invrVec,hyVec)
    ! Set values to block matrix
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: CSR_Mat
    ! miter = | invDiagrVec     ~~~~   |
    !         |   ~~~~~~    DiagyVec/h |
    !
    INTEGER :: m, n
    REAL(RealKind), INTENT(IN) :: invrVec(:)
    REAL(RealKind), INTENT(IN) :: hyVec(:)
    !TEMP
    !
    !
    CSR_Mat%Val(CSR_Mat%RowPtr(1:n))=invrVec(:)
    CSR_Mat%Val(CSR_Mat%RowPtr(n+2:n+m+1)-1)=hyVec(:)
    !
  END SUBROUTINE Miter_Extended
  !
  !
  ! PRINT SPARSE MATRIX (compressed Row format)
  SUBROUTINE PrintRhs(Rhs,FileName)
    REAL(RealKind) :: Rhs(:)
    CHARACTER(*), OPTIONAL :: FileName
    !
    INTEGER :: i
    !
    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName))//'.Rhs',STATUS='UNKNOWN')
    !
    DO i=1,SIZE(Rhs)
      WRITE(99,'(1X,E23.14)') Rhs(i)
    END DO
  END SUBROUTINE PrintRhs
  !
  !
  ! PRINT SPARSE MATRIX (compressed Row format)
  SUBROUTINE PrintSparseMatrix(A,FileName)
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    CHARACTER(*), OPTIONAL :: FileName
    !
    INTEGER :: i,j,jj
    !
    !
    WRITE(*,*) 'Sparse Matrix: ', FileName
    WRITE(*,*) 'dimension:     ', A%m,A%n
    WRITE(*,*) 'nonzeros:      ', SIZE(A%ColInd)
    ! 
    DO i=1,A%m
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        WRITE(*,'(1X,I5,1X,I5,10X,E23.14)') i,j,A%Val(jj)
      END DO
    END DO
  END SUBROUTINE PrintSparseMatrix
  !
  !
  !
  SUBROUTINE WriteSparseMatrix(A,FileName,nr,ns)
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    CHARACTER(*), OPTIONAL :: FileName
    INTEGER, OPTIONAL :: nr,ns
    !
    INTEGER :: i,j,jj
    !

    ! 
    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName))//'.SparseMat',STATUS='UNKNOWN')

    WRITE(99,*) '###########################################################'
    WRITE(99,*) '##############  Sparse Matrix  Matlab input ###############'
    WRITE(99,*) '###########################################################'
    WRITE(99,*) 
    WRITE(99,*) '###########################################################'
    WRITE(99,*) '#     Name:        ' , ADJUSTL(FileName)
    WRITE(99,*) '#'  
    WRITE(99,*) '#     Dimension:   ' , A%m,' x ',A%n
    WRITE(99,*) '#     Nonzeros:    ' , A%nnz
    WRITE(99,*) '#     nreac, nspc: ' , nr,' , ',ns
    WRITE(99,*) '###########################################################'
    WRITE(99,*)

    WRITE(99,*) '###########################################################'
    WRITE(99,*) '# Sparse Matrix in row/column - index format'
    WRITE(99,*) '###########################################################'
    WRITE(99,*) 'MATRIX'
    DO i=1,A%m
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        WRITE(99,'(1X,I12,10X,I12,10X,E23.14)') i,j,A%Val(jj)
      END DO
    END DO
    WRITE(99,*)
    WRITE(99,*)

    IF (ALLOCATED(A%DiagPtr)) THEN
      WRITE(99,*) '###########################################################'
      WRITE(99,*) '# Array for diagonal entries'
      WRITE(99,*) '###########################################################'
      WRITE(99,*) 'DIAG_PTR'
      DO i=1,SIZE(A%DiagPtr)
        WRITE(99,'(1X,I12,10X,I12)') i,A%DiagPtr(i)
      END DO
      WRITE(99,*)
      WRITE(99,*)
    END IF


    IF ( EXTENDED ) THEN
      IF (ALLOCATED(A%DiagPtr_R)) THEN
        WRITE(99,*) '###########################################################'
        WRITE(99,*) '# Array for reaction rate diagonal entries  '
        WRITE(99,*) '#                        (only EXTENDED case) '
        WRITE(99,*) '###########################################################'
        WRITE(99,*) 'DIAG_R_PTR'
        DO i=1,SIZE(A%DiagPtr_R)
          WRITE(99,'(1X,I12,10X,I12)') i,A%DiagPtr_R(i)
        END DO
        WRITE(99,*)
        WRITE(99,*)
      END IF

      IF (ALLOCATED(A%DiagPtr_C)) THEN
        WRITE(99,*) '###########################################################'
        WRITE(99,*) '# Array for species concentration diagonal  '
        WRITE(99,*) '#                        (only EXTENDED case) '
        WRITE(99,*) '###########################################################'
        WRITE(99,*) 'DIAG_C_PTR'
        DO i=1,SIZE(A%DiagPtr_C)
          WRITE(99,'(1X,I12,10X,I12)') i,A%DiagPtr_C(i)
        END DO
        WRITE(99,*)
        WRITE(99,*)
      END IF
    END IF

    IF (ALLOCATED(A%Permu)) THEN
      WRITE(99,*) '###########################################################'
      WRITE(99,*) '# Array for LU permutation and inverse permutation '
      WRITE(99,*) '###########################################################'
      WRITE(99,*) 'DIAG_PERMUTATION'
      DO i=1,SIZE(A%Permu)
        WRITE(99,'(1X,I12,10X,I12)') i,A%Permu(i), A%InvPer(i)
      END DO
      WRITE(99,*)
      WRITE(99,*)
    END IF

    CLOSE(99)
  END SUBROUTINE WriteSparseMatrix
  !
  !
  ! PRINT SPARSE MATRIX (compressed Row format)
  SUBROUTINE PrintSparse2(n,rowind,colind,val,FileName)
    INTEGER :: rowind(:)
    INTEGER :: colind(:)
    REAL(RealKind) :: val(:)
    INTEGER :: n
    CHARACTER(*), OPTIONAL :: FileName
    !
    INTEGER :: i
    !
    WRITE(*,*) 'Print Matrix: ',FileName
    WRITE(*,*) 'Dim: ',n,' x ',n , 'nnz: ',SIZE(ColInd)
    WRITE(*,*)
    ! 
    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName))//'.SparseMat',STATUS='UNKNOWN')
    DO i=1,SIZE(RowInD)
      IF (FileName=='*') THEN     !nur Hauptdiagonale ausgeben
        WRITE(*,'(1X,I5,1X,I5,10X,E23.14)') RowInd(i),ColInd(i),val(i)
      ELSE
        WRITE(99,'(1X,I5,1X,I5,10X,E23.14)') RowInd(i),ColInd(i),val(i)
      END IF
    END DO
    CLOSE(99)
  END SUBROUTINE PrintSparse2
  !
  !
  !
  !
  ! PRINT SPARSE MATRIX (RowInd/ColInd format)
  SUBROUTINE PrintMatrix(A)
    TYPE(SpRowIndColInd_T), INTENT(IN) :: A
    !
    INTEGER :: j
    !
    WRITE(*,*) '    i ','    j   ','           Wert  '
    WRITE(*,*) '------------------------------------------------'
    DO j=1,SIZE(A%ColInd)
      WRITE(*,'(1X,I5,1X,I5,10X,E12.4)') A%RowInd(j),A%ColInd(j),A%Val(j)
    END DO
  END SUBROUTINE PrintMatrix
  !
  !
  ! convert compressed row format to rowindex column index format
  SUBROUTINE CompRowToIndRow(MatIn,MatOut)
    ! only n by n matrices
    TYPE(CSR_Matrix_T), INTENT(IN) :: MatIn
    TYPE(SpRowIndColInd_T), INTENT(OUT) :: MatOut
    !
    INTEGER :: i,j,jj
    MatOut%m=MatIn%m
    MatOut%n=MatIn%RowPtr(MatIn%m+1)-1
    !
    ALLOCATE(MatOut%RowInd(MatOut%n))
    MatOut%RowInd=0
    ALLOCATE(MatOut%ColInd(MatOut%n))
    MatOut%ColInd=0
    ALLOCATE(MatOut%Val(MatOut%n))
    MatOut%Val=ZERO
    j=1
    DO i=1,MatIn%m
      DO jj=MatIn%RowPtr(i),MatIn%RowPtr(i+1)-1
        MatOut%RowInd(j)=i
        MatOut%ColInd(j)=MatIn%ColInd(jj)
        MatOut%Val(j)=MatIn%Val(jj)
        j=j+1
      END DO
    END DO
  END SUBROUTINE CompRowToIndRow
  !
  !
  ! Matrix*Vector1+Vector2 (rhs)
  SUBROUTINE DiffConcDt(DcDt,A,Vec1)
    !IN
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    REAL(RealKind), INTENT(IN) :: Vec1(:)
    !OUT
    REAL(RealKind), INTENT(OUT) :: DcDt(A%m)
    REAL(RealKind) :: Tmp
    INTEGER :: i,jj
    !
    DO i=1,A%m
      !Tmp=ZERO
      !DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      !  Tmp=Tmp+A%Val(jj)*Vec1(A%ColInd(jj))
      !END DO
      !DcDt(i)=Tmp
      DcDt(i) = SUM(  A%Val(A%RowPtr(i):A%RowPtr(i+1)-1 )         &
      &               * Vec1(A%ColInd(A%RowPtr(i):A%RowPtr(i+1)-1)))
    END DO
  END SUBROUTINE DiffConcDt
  !
  !
  ! Matrix*Vector1+Vector2 (rhs)
  SUBROUTINE MatVecMult(Rhs,A,Vec1,Vec2)
    !IN
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    REAL(RealKind), INTENT(IN) :: Vec1(:), Vec2(:)
    !OUT
    REAL(RealKind), INTENT(OUT) :: Rhs(A%m)
    !TEMP
    REAL(RealKind) :: Tmp
    INTEGER :: i,jj
    !
    DO i=1,A%m
      Rhs(i) = SUM( A%Val(A%RowPtr(i):A%RowPtr(i+1)-1)              &
      &             * Vec1(A%ColInd(A%RowPtr(i):A%RowPtr(i+1)-1)) ) &
      &         + Vec2(i)
    END DO
    !write(343,*) '-----------',SIZE(A%val)
    !do i=1,SIZE(A%val)
    ! write(343,*) A%val(i)
    !end do 
    !write(343,*) '-----------',SIZE(vec1)
    !do i=1,SIZE(vec1)
    ! write(343,*) vec1(i)
    !end do 
    ! write(343,*) '-----------',SIZE(vec2)
    !do i=1,SIZE(vec2)
    ! write(343,*) vec2(i)
    !end do
    !stop
  END SUBROUTINE MatVecMult
  !
  !
  ! Print ordinary matrix or vector
  SUBROUTINE PM(M)
    REAL(RealKind), INTENT(IN) :: M(:,:)
    INTEGER :: i
    !
    DO i=1,SIZE(M,1)
        WRITE(*,*) M(i,:)
    END DO
  END SUBROUTINE PM
  !
  !
  ! Print ordinary  vector
  SUBROUTINE PV(V)
    REAL(RealKind), INTENT(IN) :: V(:)
    INTEGER :: i
    !
    DO i=1,SIZE(V)
        WRITE(*,*) i, V(i)
    END DO
  END SUBROUTINE PV
  !
  !
  SUBROUTINE WriteSpRowIndColInd_T(Mat,rhs)
    TYPE(sprowindcolind_t) :: Mat
    REAL(RealKind) :: rhs(:)
    !
    INTEGER :: i
    !
    !
    199 format(I6,4X,I6,4X,E18.12,4X,E18.12)
    !OPEN(UNIT=91,FILE='ChemieMat'//'.spmat',STATUS='UNKNOWN')
    !
    WRITE(*,*) '  '
    WRITE(*,*) '  i        j          val               rhs   '
    WRITE(*,*) '=====================================================================>'
    WRITE(*,*) ' '
    !
    DO i=1,SIZE(Mat%RowInd)
      WRITE(*,199) Mat%RowInd(i),Mat%ColInd(i),Mat%Val(i),rhs(i)
    END DO 
    !CLOSE(91)  
  END SUBROUTINE WriteSpRowIndColInd_T
  !
  !
  SUBROUTINE PermuToInvPer(InvPer,Permu)
    INTEGER, ALLOCATABLE :: InvPer(:)
    INTEGER :: Permu(:)
    !
    INTEGER :: i
    !
    IF (.NOT.ALLOCATED(InvPer)) ALLOCATE(InvPer(SIZE(Permu)))
    !
    DO i=1,SIZE(Permu)
      InvPer(Permu(i))=i
    END DO  
  END SUBROUTINE PermuToInvPer
  !
  !
  SUBROUTINE PrintPerm(Permu,InvPermu,FileName)
    INTEGER :: Permu(:), InvPermu(:)
    CHARACTER(*), OPTIONAL :: FileName
    !
    INTEGER :: i
    !
    WRITE(*,*) 'Print Permu: ',FileName
    WRITE(*,*) 'Dim: ',SIZE(Permu)
    WRITE(*,*)
    ! 
    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName))//'.Permu',STATUS='UNKNOWN')
    DO i=1,SIZE(Permu)
      IF (FileName=='*') THEN     !nur Hauptdiagonale ausgeben
        WRITE(*,'(1X,I5,1X,I5)') Permu(i),InvPermu(i)
      ELSE
        WRITE(99,'(1X,I5,1X,I5)') Permu(i),InvPermu(i)
      END IF
    END DO
    CLOSE(99)

  END SUBROUTINE PrintPerm
  !
  !
  !
  SUBROUTINE SparseLU(A)
    TYPE(CSR_Matrix_T) :: A
    !
    REAL(RealKind) :: w(A%n)
    REAL(RealKind) :: alpha
    INTEGER :: i,j,jj,kk
    !INTEGER :: OPCOUNT

    !OPCOUNT=0
    !
    DO i=1,A%n
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        w(A%ColInd(jj))=A%Val(jj)
      END DO
      DO jj=A%RowPtr(i),A%DiagPtr(i)-1
        j=A%ColInd(jj)
        alpha=w(j)/A%Val(A%DiagPtr(j))
        !OPCOUNT=OPCOUNT+1
        w(j)=alpha
        DO kk=A%DiagPtr(j)+1,A%RowPtr(j+1)-1
          w(A%ColInd(kk))=w(A%ColInd(kk))-alpha*A%Val(kk)
          !OPCOUNT=OPCount+1
        END DO
      END DO
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        A%Val(jj)=w(A%ColInd(jj))
      END DO
    END DO
    !print*, 'opcount=', opcount
    !stop
  END SUBROUTINE SparseLU
  !
  !
  SUBROUTINE SolveSparse(LU,rhs)
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: LU
    REAL(RealKind) :: Rhs(:)
    !
    INTEGER :: i,jj
    REAL(RealKind) :: b(LU%n)
    REAL(RealKind) :: DiagPtr(LU%n)
    !
    !DiagPtr(:)=(/ LU%DiagPtr , LU%DiagPtr2  /)
    !--  Permutation of right hand side
    DO i=1,LU%n
      b(LU%Permu(i)) = Rhs(i)      
    END DO
    !
    !--  L-solve
    DO i=2,LU%n
      DO jj=LU%RowPtr(i),LU%DiagPtr(i)-1
        b(i)=b(i)-LU%Val(jj)*b(LU%ColInd(jj))
      END DO
    END DO
    !
    !--  U-solve
    DO i=LU%n,1,-1
      DO jj=LU%DiagPtr(i)+1,LU%RowPtr(i+1)-1
        b(i)=b(i)-LU%Val(jj)*b(LU%ColInd(jj))
      END DO
      b(i)=b(i)/LU%Val(LU%DiagPtr(i))
    END DO
    !
    !--  Back-Permutation of solution
    DO i=1,LU%n
      Rhs(LU%InvPer(i)) = b(i)       
    END DO
  END SUBROUTINE SolveSparse
END MODULE Sparse_Mod
