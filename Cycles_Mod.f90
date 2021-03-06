! WRITE(*,'(A)') '  Breakpoint:   '; READ(*,*)
MODULE Cycles_Mod
  USE Kind_Mod,   ONLY: dp
  USE Sparse_Mod, ONLY: CSR_Matrix_T, SortVec, SortVecDesc,     &
                      & SortVecDesc2, SpRowIndColInd_T,         &
                      & CSR_to_SpRowIndColInd, Free_Matrix_CSR, &
                      & New_CSR, Copy_CSR, CompressIntegerArray,&
                      & PrintSparseMatrix

  USE Reac_Mod,    ONLY: y_name
  USE Control_Mod, ONLY: List

  IMPLICIT NONE


  TYPE Node_T
    INTEGER :: Id, cnt=0
    INTEGER,  ALLOCATABLE :: v(:)
    TYPE(Node_T), POINTER :: next=>NULL()
  END TYPE Node_T

  TYPE Tarjan_T
    INTEGER :: end
    INTEGER, ALLOCATABLE :: index(:)
    INTEGER, ALLOCATABLE :: lowLink(:)
    LOGICAL, ALLOCATABLE :: onStack(:)
  END TYPE Tarjan_T

  TYPE(Node_T),   POINTER, PRIVATE :: cycFirst, cyc, sccFirst, scc_List
  TYPE(Tarjan_T),          PRIVATE :: Tar
  TYPE(List),              PRIVATE :: Stapel
  TYPE(List), ALLOCATABLE, PRIVATE :: Tar_SCC(:), blocked_List(:), Graph(:)
  INTEGER,    ALLOCATABLE, PRIVATE :: Tar_stack(:)
  INTEGER,                 PRIVATE :: Tar_idx, cntSCC
  LOGICAL,    ALLOCATABLE, PRIVATE :: blocked_Vertex(:)

  INTEGER :: lvl
  INTEGER, PARAMETER :: maxlen = 2

  CONTAINS
 
  ! This subroutine computes the elementary circuits (cycles) of a direced graph
  ! acording to Johnsons algorithm (see: https://doi.org/10.1137/0204007) 
  !
  ! INPUT: (IA,JA) if there is a edge from vertex I to vertex J
  !        length(IA) = length(JA) = nnz(A)
  !
  FUNCTION Find_Elem_Circuits(A,FAMS) RESULT(Cyclic_Set_Out)
  !SUBROUTINE Find_Elem_Circuits(A,SpcList)
    USE Control_Mod,   ONLY: Families_T, BSP, OutputPath

    TYPE(List), ALLOCATABLE :: Cyclic_Set_Out(:) 
    ! IN:
    TYPE(CSR_Matrix_T)   :: A, sub_Ak
    TYPE(Families_T)     :: FAMS(:)
    INTEGER, ALLOCATABLE :: SpcList(:)

    TYPE(SpRowIndColInd_T) :: sub_A
    TYPE(List), ALLOCATABLE :: scc(:) !, sub_Ak(:)
    INTEGER :: n, s, least_node, nspc
    INTEGER :: i, j, k, ispc , cnt, nnz, iFam, nFam
    INTEGER :: iC
    INTEGER, ALLOCATABLE :: comp_nodes(:)
    INTEGER, ALLOCATABLE :: iCyc_red(:)
    LOGICAL :: dummy
    TYPE(Node_T), POINTER :: printCYC

    TYPE(List), ALLOCATABLE :: Cyclic_Set(:)
    INTEGER                 :: nCycles=0, nCycles_red=0

     777 FORMAT(10X,A)

    IF ( A%m /= A%n ) THEN
      WRITE(*,*) '   Input matrix is no square matrix   '
      RETURN
    END IF

    WRITE(*,777) REPEAT('*',39)
    WRITE(*,777) '********** Structure Analysis *********'
    WRITE(*,777) REPEAT('*',39)
    WRITE(*,*)

    n    = A%n
    nnz  = A%nnz

    sub_A = CSR_to_SpRowIndColInd(A)

    ! Initialize global arrays
    ALLOCATE( blocked_Vertex(n) , blocked_List(n) , Stapel%List(0) )
    blocked_Vertex = .FALSE.

    DO s=1,n; ALLOCATE( blocked_List(s)%list(0) ); END DO
    blocked_List%len = 0

    Stapel%len  = 0
    Stapel%List = 0

    ALLOCATE(Cyc); NULLIFY(Cyc%next)

    cycFirst => Cyc
    printCYC => Cyc
    
    ! diese schleife später nur über ausgewählte knoten s
    s = 1
    nFam = SIZE(FAMS)

    DO iFam = 1 , nFam

      SpcList = [ FAMS(iFam)%Index ]
      nspc = SIZE(SpcList)

      DO j = 1 , nspc

        ! Generate Subgraph 
        IF ( j > 1 ) CALL Purge_Vertex( sub_A , SpcList(j-1) ) 

        ! find strong connected components (Tarjan algorithm)
        CALL Tarjan( scc , ispc , SpcList(j) , sub_A )

        IF ( scc(ispc)%len > 1 ) THEN
          
          comp_nodes = scc(ispc)%List
          least_node = SpcList(j)

          sub_Ak = Generate_Submatrix( comp_nodes , sub_A )

          blocked_Vertex = .FALSE.
          DO i = 1,scc(ispc)%len
            IF (ALLOCATED(blocked_List(i)%List)) THEN
              DEALLOCATE(blocked_List(i)%List)
              ALLOCATE(blocked_List(i)%List(0))
            END IF
            blocked_List(i)%len = 0
          END DO

          s = least_node

          dummy = Circuit( s, s, sub_Ak )

          DEALLOCATE(scc)
          s = s + 1
          CALL Progress(iFam,nFam,j,nspc)
        ELSE
          CALL Progress(iFam,nFam,nspc,nspc)
          EXIT
        END IF
      END DO

    END DO ! 1,SIZE(FAM)
    
    cntSCC = 0
    printCYC => cycFirst
    DO WHILE (ASSOCIATED(printCYC))
      IF ( printCYC%cnt > 1 ) cntSCC = cntSCC + 1
      printCYC => printCYC%next
    END DO
    nCycles = cntSCC
    ALLOCATE(Cyclic_Set(nCycles))

    ! write paths to file and save it in a allcatable array
    OPEN(UNIT=99,FILE=OutputPath//'all_species_cycles_'//TRIM(BSP)//'.txt',STATUS='UNKNOWN')
    WRITE(99,*)
    WRITE(99,'(A,*(I0))') '   Number of all found cycles:   ',nCycles
    WRITE(99,'(A)') '   Cycle Length:       species: '
    printCYC => cycFirst
    i = 0
    DO 
      IF ( ASSOCIATED(PrintCYC%next) ) THEN
        IF ( printCYC%cnt > 1 ) THEN

          i = i + 1

          Cyclic_Set(i)%len  = printCYC%cnt+1
          Cyclic_Set(i)%List = [printCYC%v]

          ! zum anschauen
          WRITE(99,'(10X,I2,8X,*(A))') Cyclic_Set(i)%len,                      &
          & (TRIM(y_name(Cyclic_Set(i)%List(j)))//' -> ' ,j=1,Cyclic_Set(i)%len-1), &
          &  TRIM(y_name(Cyclic_Set(i)%List(Cyclic_Set(i)%len)))


        END IF
      ELSE
        EXIT
      END IF
      printCYC => printCYC%next
    END DO
    CLOSE(99)

    
    ! write paths to file and save it in a allcatable array
    cntSCC = 0
    ALLOCATE(iCyc_red(0))
    DO iC = 1, nCycles

      iFam = Find_Familie(Cyclic_Set(iC)%List(1))

      IF ( iFam > 0 ) THEN
        cnt = 0
        DO i = 2 , Cyclic_Set(iC)%len-1
          IF ( Is_Spc_In_Cycle( Cyclic_Set(iC)%List(i)             &
            &                 , FAMS(iFam)%Index                   &
            &                 , SIZE(FAMS(iFam)%Index) ) > 0 )  THEN
            cnt = cnt + 1
          END IF
        END DO

        IF ( cnt == Cyclic_Set(iC)%len-2 ) THEN 
          cntSCC = cntSCC + 1
          iCyc_red = [iCyc_red , iC]
        END IF

      END IF

    END DO
    
    WRITE(*,*); WRITE(*,*) 
    WRITE(*,'(10X,A,I0)') '    Number of cycles = ',cntSCC
    WRITE(*,*); WRITE(*,*) 

    nCycles_red = cntSCC
    ALLOCATE(Cyclic_Set_Out(nCycles_red))
    !OPEN(UNIT=98,FILE=OutputPath//'family_species_cycles_'//TRIM(BSP)//'.txt',STATUS='UNKNOWN')
    !WRITE(98,*)
    !WRITE(98,'(A,*(I0))') '   Number of cycles afer sorting out:   ',nCycles_red
    !WRITE(98,'(A)') '   Cycle Length:       species: '
    DO i = 1, nCycles_red
      iC = iCyc_red(i)
      Cyclic_Set_Out(i)%len = Cyclic_Set(iC)%len
      Cyclic_Set_Out(i)%List = [Cyclic_Set(iC)%List]
      !WRITE(98,'(10X,I2,8X,*(A))') Cyclic_Set_Out(i)%len,                      &
      !& (TRIM(y_name(Cyclic_Set_Out(i)%List(k)))//' -> ' ,k=1,Cyclic_Set_Out(i)%len-1), &
      !&  TRIM(y_name(Cyclic_Set_Out(i)%List(Cyclic_Set_Out(i)%len)))
      ! matlab output
      !WRITE(97,'(*(I0,2X))') (Cyclic_Set_Out(i)%List(j) , j=1,Cyclic_Set_Out(i)%len)
    END DO

    CLOSE(98)


    CONTAINS

    FUNCTION Find_Familie(iSpc) RESULT(iFa)
      INTEGER :: iFa
      INTEGER :: iSpc
      INTEGER :: ii, jj

      iFa=0
      DO ii = 1 , SIZE(FAMS)
        DO jj = 1 , SIZE(FAMS(ii)%Index)
          IF ( FAMS(ii)%Index(jj) == iSpc ) THEN
            iFa = ii
            RETURN
          END IF
        END DO
      END DO
    END FUNCTION Find_Familie

    FUNCTION Is_Spc_In_Cycle(iSpc,List,n) RESULT(iRow)
      INTEGER             :: iRow
      INTEGER, INTENT(IN) :: List(:)
      INTEGER, INTENT(IN) :: iSpc, n
      INTEGER :: i

      iRow=0
      DO i=1,n
        IF ( List(i)==iSpc ) THEN
          iRow=i  
          RETURN
        END IF
      END DO
    END FUNCTION Is_Spc_In_Cycle

  END FUNCTION Find_Elem_Circuits

  
  RECURSIVE SUBROUTINE Unblock(u)
    INTEGER :: u, w

    blocked_Vertex(u) = .FALSE.

    DO w = 1 , blocked_List(u)%len
      IF (blocked_Vertex(w)) CALL Unblock(w)
    END DO

    IF (ALLOCATED(blocked_List(u)%List)) THEN
      DEALLOCATE(blocked_List(u)%List)
      ALLOCATE(blocked_List(u)%List(0))
    END IF
    blocked_List(u)%len = 0
  END SUBROUTINE Unblock


  RECURSIVE LOGICAL FUNCTION Circuit(v,s,Ak) RESULT(f)
    ! INPUT:
    INTEGER    :: v, s
    TYPE(CSR_Matrix_T) :: Ak
    ! TEMP:
    INTEGER :: w, ww, j
    INTEGER, ALLOCATABLE :: Bnode(:)
    LOGICAL :: found_v

    LOGICAL :: dbg=.FALSE. !, inFAM

    f = .FALSE.
    !inFAM = .FALSE.
    lvl = lvl + 1
    !WRITE(*,*) ' level = ', lvl

    Stapel%len  = Stapel%len + 1
    Stapel%List = [Stapel%List , v]
    blocked_Vertex(v) = .TRUE.
    
    !IF (dbg) THEN
    !  333 FORMAT(A,*(I0,1X))
    !  334 FORMAT(A,*(L0,1X))
    !  WRITE(*,333) ' Stapel len   = ',Stapel%len 
    !  WRITE(*,333) ' Stapel List  = ',Stapel%List
    !  !WRITE(*,334) ' blocked vert = ',blocked_Vertex 
    !  WRITE(*,*) ''
    !END IF

    
    DO ww = Ak%RowPtr(v),Ak%RowPtr(v+1)-1
      w = Ak%ColInd(ww)

      !DO j = 1,nallFAM;  IF (allFAM(j)==w) inFAM=.TRUE.;  END DO
      
      ! found new cycle if w == s
      IF ( w == s ) THEN
        f = .TRUE.
        Cyc%v = [Stapel%List , s]
        Cyc%cnt = SIZE(Cyc%v) - 1

        !WRITE(*,'(A,I3,A,*(I0,1X))') '   Kreis der Länge: ',SIZE([Stapel%List , s]),' mit den Knoten => ',[Stapel%List , s]
        ALLOCATE(Cyc%next); Cyc => Cyc%next

      ! if no new cycle is found check if vertex w is
      ! blocked jet, if it is not blocked continue 
      ! ciruit with node w

      ! AN DER STELLE DIE FAMILIENLISTE DURCHGEHEN UND NUR WEITERMACHEN
      ! WENN w ALS FAMILIENMITGLIED DEKLARIERT WURDE, AUßERDEM WEITER 
      ! UNTEN NOCHMAL SCHAUEN OB ES NOCHMAL NÖTIG IST DIE LISTE ZU DURCHLAUFEN
      ELSE IF ( .NOT.blocked_Vertex(w) .AND. Stapel%len <= maxlen ) THEN
        f = Circuit(w,s,Ak)
        !inFAM = .FALSE.
      END IF
    END DO

    ! if a new cycle was found unblock node v
    IF ( f .OR. Stapel%len > maxlen ) THEN
      CALL Unblock(v)
    ! otherwise check neighbours of node v and block them
    ! if there not in the blocked_List jet
    ELSE
      DO ww = Ak%RowPtr(v),Ak%RowPtr(v+1)-1
        w = Ak%ColInd(ww)
        found_v = (findF(blocked_List(w)%List , v) /= 0)

        IF ( .NOT.found_v ) THEN
          blocked_List(w)%List = [blocked_List(w)%List , v]
          !IF (dbg) WRITE(*,333) ' blocked List  = ',blocked_list(w)%List
        END IF

      END DO
    END IF

    Stapel%len  = Stapel%len - 1
    Stapel%List = Stapel%List(1:Stapel%len)

    !IF (dbg) THEN
      !WRITE(*,*) '  OUT'
      !WRITE(*,333) ' Stapel len   = ',SIZE(Stapel%LIst) 
      !WRITE(*,333) ' Stapel List  = ',Stapel%List
      !WRITE(*,334) ' blocked vert = ',blocked_Vertex 
      !WRITE(*,*) ''
      !WRITE(*,*) '      NEXT    '
      !READ(*,*)
      !WRITE(*,*) ''
    !END IF
  END FUNCTION Circuit

  FUNCTION CSR_to_AdjList(CSR) RESULT(Adj)
    ! INPUT:
    TYPE(CSR_Matrix_T)  :: CSR
    ! OUTPUT:
    TYPE(List), ALLOCATABLE :: Adj(:)
    ! TEMP:
    INTEGER :: i 
  
    ALLOCATE(Adj(CSR%m))
  
    DO i = 1,CSR%m
      Adj(i)%List = [CSR%ColInd(CSR%RowPtr(i):CSR%RowPtr(i+1)-1)]
      Adj(i)%len  = CSR%RowPtr(i+1)-CSR%RowPtr(i)
    END DO
  
  END FUNCTION CSR_to_AdjList


  SUBROUTINE Tarjan(scc,ispc,Target_Spc,A)
    ! OUT:
    TYPE(List), ALLOCATABLE :: scc(:)
    INTEGER :: ispc
    ! IN:
    INTEGER :: Target_Spc
    TYPE(SpRowIndColInd_T)  :: A
    ! TEMP:
    INTEGER :: i , j , n, nSCC, tmpV, least_node
    INTEGER, ALLOCATABLE :: nVertexSCC(:), q(:)
    
    n = A%n
    least_node = MINVAL([A%RowInd , A%ColInd])
    
    ALLOCATE(scc_List); NULLIFY(scc_List%next)
    cntSCC   = 0
    sccFirst => scc_List

    ALLOCATE(Tar%index(n), Tar%lowLink(n), Tar%onStack(n))
    Tar_idx     = 0
    Tar%index   = 0
    Tar%lowLink = 0
    Tar%onStack = .FALSE.
    
    ! loop through all nodes
    DO i = 1,n
      IF (Tar%index(i)==0) CALL StrongConnect(i,A)
    END DO

    nSCC = scc_List%id-1
    
    IF (ALLOCATED(scc)) DEALLOCATE(scc)
    ALLOCATE( scc(nSCC), nVertexSCC(nSCC), q(nSCC) )
    nVertexSCC = 0;  q = 0

    scc_List => sccFirst

    ! generate count vector
    DO i = 1, nSCC
      nVertexSCC(i) = SIZE(scc_List%v)

      ! hier wäre es von vorteil sowohl knoten UND kanten 
      ! der zusammenhangskomponenten zu speichern
      ALLOCATE( scc(i)%List(nVertexSCC(i)) )

      scc(i)%len  = nVertexSCC(i)
      scc(i)%List = scc_List%v
      scc_List   => scc_List%next
      !WRITE(*,'(I2,A,I2,A,*(I2,1X))'), i , ' len = ', scc(i)%len , ' List = ', scc(i)%List
      DO j = 1,nVertexSCC(i)
        IF (scc(i)%List(j) == Target_Spc) ispc = i
      END DO
      
    END DO


    ! sort by decreasing size
    !CALL SortVecDesc2(nVertexSCC,q)
    !scc = scc(q)
    
    DEALLOCATE(Tar%index, Tar%lowLink, Tar%onStack)
  END SUBROUTINE Tarjan


  RECURSIVE SUBROUTINE StrongConnect(i,A)
      
    TYPE(SpRowIndColInd_T) :: A
    INTEGER :: i

    INTEGER, ALLOCATABLE :: out_Edges(:)
    INTEGER, ALLOCATABLE :: theSCC(:)
    INTEGER :: n, k, j, nSCC

    LOGICAL :: dbg=.false.
  
    Tar_idx = Tar_idx + 1

    ! visit the node
    Tar%index(i)   = Tar_idx
    Tar%lowLink(i) = Tar_idx
    Tar%onStack(i) = .TRUE.


    IF (.NOT.ALLOCATED(Tar_stack)) ALLOCATE(Tar_stack(0))

    Tar_stack = [i , Tar_stack]

    IF (dbg) THEN
      WRITE(*,*) ''
      WRITE(*,'(A,I0,A,*(I3))') '   Tar%index(',i,') = ',Tar%index(i)
      WRITE(*,'(A,I0,A,*(I3))') ' Tar%lowLink(',i,') = ',Tar%lowLink(i)
      WRITE(*,'(A,I0,A,*(L))')  ' Tar%onStack(',i,') =  ',Tar%onStack(i)
      WRITE(*,'(A,*(I2,1X))')      '       Tar_stack = ',Tar_stack
    END IF

    IF ( A%n < i ) THEN
      WRITE(*,*) ' Permit operation on rectangular matrices. '
      STOP ' StrongComponents '
    END IF

    ALLOCATE(out_Edges(0))

    DO j=1,SIZE(A%ColInd)
      IF ( A%ColInd(j)==i ) THEN
        out_Edges = [out_Edges , A%RowInd(j)]
      END IF
    END DO

    n = SIZE(out_Edges)
    IF (dbg) WRITE(*,'(A,*(I3,1X))') '  ---> out_Edges = ', out_Edges
    
    DO k = 1 , n
      j = out_Edges(k)

      IF (dbg) THEN
        WRITE(*,'(A,I0,A,*(I3))') '   Tar%index(',j,') = ',Tar%index(j)
        WRITE(*,'(A,I0,A,*(L3))') ' Tar%onStack(',j,') = ',Tar%onStack(j)
      END IF

      ! if not visited, visit
      IF ( Tar%index(j) == 0 ) THEN

        CALL StrongConnect(j,A)
        ! carry back lowlink, if lower
        Tar%lowLink(i) = MINVAL( Tar%lowLink([i,j]) )

      END IF
      IF ( Tar%onStack(j) ) THEN
      
        IF (dbg) WRITE(*,'(2(A,I0),A,*(I3))') &
          & ' MIN(( Tar%lowLink(',i,') , Tar%index(',j,') ))', &
          & Tar%lowLink(i) ,Tar%index(j)
        ! carry back index, if lower
        Tar%lowLink(i) = MINVAL( [Tar%lowLink(i) , Tar%index(j)] )
     
      END IF
      IF (dbg) WRITE(*,'(A,I0,A,*(I3))') &
        &' Tar%lowLink(',i,') = ',Tar%lowLink(i)
    END DO
  
    IF (dbg) WRITE(*,'(A,I0,A,I2,A,*(I3))') &
      &' Tar%lowLink(',i,') = Tar%index(',i,') ',Tar%lowLink(i) ,Tar%index(i)
    

    IF ( Tar%lowLink(i) == Tar%index(i) ) THEN
      ! label a new SCC
      cntSCC = cntSCC + 1
      theSCC = [ Tar_stack( 1 : findL(Tar_stack,i) ) ]
      nSCC   = SIZE(theSCC)
      CALL SortVec(theSCC)
      Tar_stack = [ Tar_stack(nSCC+1:) ]
      Tar%onStack(theSCC) = .FALSE.

      ALLOCATE(scc_List%v(nSCC))
      scc_List%Id = cntSCC;    scc_List%v  = theSCC

      ALLOCATE(scc_List%next); scc_List => scc_List%next
      scc_List%Id = cntSCC + 1
    END IF

    IF (dbg) THEN
      WRITE(*,*) ''
      WRITE(*,*) '     next step     '
      WRITE(*,*) ''
      !READ(*,*)
    END IF

  
  END SUBROUTINE StrongConnect
  
  ! functoin finds the last index j for entry a_j == i
  FUNCTION findL(Array,i) RESULT(Idx)
    INTEGER :: Idx
    INTEGER :: Array(:)
    INTEGER :: i
    INTEGER :: j
    DO j = 1,SIZE(Array)
      IF (Array(j)==i) Idx = j !; EXIT
    END DO
  END FUNCTION findL

  FUNCTION findF(Array,i) RESULT(Idx)
    INTEGER :: Idx
    INTEGER :: Array(:)
    INTEGER :: i
    INTEGER :: j
    Idx = 0
    DO j = 1,SIZE(Array)
      IF (Array(j)==i) THEN
        Idx = j 
        RETURN
      END IF
    END DO
  END FUNCTION findF


  FUNCTION find2(Arr1,Arr2,i) RESULT(Idx)
    INTEGER, ALLOCATABLE :: Idx(:)
    INTEGER :: Arr1(:), Arr2(:)
    INTEGER :: i
    INTEGER :: j
    ALLOCATE(Idx(0))
    DO j = 1,SIZE(Arr1)
      IF (Arr1(j)==i.OR.Arr2(j)==i) THEN
        Idx = [Idx , j]
      END IF
    END DO
  END FUNCTION find2

  SUBROUTINE CompressIntegerArray2(Array,Array2)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: Array(:),Array2(:)
    INTEGER, ALLOCATABLE :: tmpArray(:,:)
    
    INTEGER :: i, cnt
    
    ALLOCATE(tmpArray(COUNT(Array/=0),2))

    cnt = 0
    DO i=1,SIZE(Array)
      IF (Array(i)/=0) THEN
        cnt = cnt + 1
        tmpArray(cnt,1) = Array(i)
        tmpArray(cnt,2) = Array2(i)
      END IF
    END DO
    Array  = [tmpArray(:,1)]
    Array2 = [tmpArray(:,2)]
  END SUBROUTINE CompressIntegerArray2

  SUBROUTINE Purge_Vertex(A,s)
    TYPE(SpRowIndColInd_T), INTENT(INOUT) :: A
    INTEGER, INTENT(IN) :: s

    INTEGER :: i, j

    DO i = 1,SIZE(A%RowInd)
      IF (A%RowInd(i)==s .OR. A%ColInd(i)==s) THEN
        A%RowInd(i) = 0
        A%ColInd(i) = 0
      END IF
    END DO
    CALL CompressIntegerArray2(A%RowInd,A%ColInd)

  END SUBROUTINE Purge_Vertex

  FUNCTION Generate_Submatrix(vlist,A) RESULT(subA)
    TYPE(CSR_Matrix_T)     :: subA
    TYPE(SpRowIndColInd_T) :: A
    INTEGER :: vlist(:)

    INTEGER :: i, j1, j2, nnz, nlist, cnt
    INTEGER, ALLOCATABLE :: fixI(:)

    nnz = SIZE(A%RowInd)
    nlist = SIZE(vlist)

    !print*, ' start gen matrix : ',nnz, nlist
    !WRITE(*,'(A,*(I3))') ' start mit vlist  = ',vlist2
    !WRITE(*,'(A,*(I3))') ' start mit RowInd = ',A%RowInd
    !WRITE(*,'(A,*(I3))') ' start mit ColInd = ',A%ColInd

    ALLOCATE(fixI(nnz))

    cnt  = 0
    fixI = 0

    DO i =1,nnz
      DO j1=1,nlist
        IF (A%RowInd(i)==vlist(j1)) THEN
          DO j2 = 1,nlist
            IF (A%ColInd(i)==vlist(j2)) THEN
              cnt = cnt + 1
              fixI(cnt) = i
            END IF
          END DO
        END IF
      END DO
    END DO
    CALL CompressIntegerArray(fixI);  nnz = SIZE(fixI)

    subA = New_CSR( A%m , A%n , nnz , &
                  & A%RowInd(fixI) ,  &
                  & A%ColInd(fixI) ,  &
                  & A%val(fixI)       )

    !WRITE(*,'(A,*(I3))') ' start mit RowInd = ',A%RowInd(fixI)
    !WRITE(*,'(A,*(I3))') ' start mit ColInd = ',A%ColInd(fixI)
    !WRITE(*,*) ' new A rowptr = ',subA%RowPtr
    !WRITE(*,*) ' new A colind = ',subA%ColInd

  END FUNCTION Generate_Submatrix

  SUBROUTINE Progress(i,ii,j,k)
    INTEGER(4)  :: i,ii,j,k 
    ! print the progress bar.
    WRITE(*,'(A1,14X,4(A,I0),A,$)') char(13),'Family :: (',i,'/',ii,')    Node :: (',j,'/',k,')  processed.'
  END SUBROUTINE Progress

END MODULE Cycles_Mod
