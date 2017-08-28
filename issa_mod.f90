MODULE issa
  
  USE Kind_Mod

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE ISSA_nue(nue,cyc,RS)
    USE mo_reac,     ONLY: neq, nspc, y_name
    USE mo_control,  ONLY: List, Chain
    USE Sparse_Mod,  ONLY: CSR_Matrix_T, New_CSR
    USE ChemSys_Mod, ONLY: ReactionStruct_T
    ! IN:
    TYPE(CSR_Matrix_T)      :: nue
    TYPE(List), ALLOCATABLE :: cyc(:)
    TYPE(ReactionStruct_T), ALLOCATABLE :: RS(:)
    ! OUT:
    TYPE(CSR_Matrix_T), ALLOCATABLE :: nue_plus(:), nue_minus(:)
    ! TEMP:
    INTEGER :: i, ii, iR, iE, iP, j,  k, n_c, len_c, nE, nP, nnz_k
    INTEGER,     ALLOCATABLE :: tsi(:)
    TYPE(Chain), ALLOCATABLE :: R_k(:)

    n_c = SIZE(cyc)

    write(*,*) '_____________________________________________'
    write(*,*)
    write(*,*) '  Testing ISSA stuff: building nue+ and nue- '
    write(*,*) '_____________________________________________'
    write(*,*)

    ! build nue_plus and nue_minus (see atmospheric env. 39 (2005) page: 4343)
    ALLOCATE( nue_plus(n_c) , nue_minus(n_c) , R_k(n_c) )

    LOOP_OVER_ALL_CYCLES: DO k = 1,n_c ! Kette Nummer: k
      
      ! First we need to detect the reactions of each species cycle
      len_c = cyc(k)%len
      tsi   = cyc(k)%List

      ALLOCATE( R_k(k)%sName(len_c-1,2), R_k(k)%sIdx(len_c-1,2), R_k(k)%rIdx(len_c-1) )

      !write(*,*) REPEAT('-',32)//'  Durchsuche Kette: ',k,'  der Länge: ',len_c
      !write(*,'(5X,*(A))')  ( TRIM(y_name(cyc(k)%List(j)))//' -> ' , j=1,len_c-1 ) ,&
      !&                       TRIM(y_name(cyc(k)%List(len_c)))
      !write(*,*) 

      DO ii = 1,len_c-1
        ALLOCATE( R_k(k)%rIdx(ii)%List(0) )
        DO iR = 1,neq

          nE = SIZE(RS(iR)%Educt)
          nP = SIZE(RS(iR)%Product)

          DO iE = 1,nE
            IF ( tsi(ii) == RS(iR)%Educt(iE)%iSpecies ) THEN
              R_k(k)%sName(ii,1) = RS(iR)%Educt(iE)%Species
              R_k(k)%sIdx(ii,1)  = RS(iR)%Educt(iE)%iSpecies
            
              DO iP = 1,nP
                IF ( tsi(ii+1) == RS(iR)%Product(iP)%iSpecies ) THEN
                  R_k(k)%sName(ii,2)   = TRIM(RS(iR)%Product(iP)%Species)
                  R_k(k)%sIdx(ii,2)    = RS(iR)%Product(iP)%iSpecies
                  R_k(k)%rIdx(ii)%List = [R_k(k)%rIdx(ii)%List , iR]
                  !write(*,*) '  in Reaktion :: ',iR,'  mit  ',TRIM(RS(iR)%Line1),  &
                  !&          '  mit Vorgänger: ',TRIM(RS(iR)%Educt(iE)%Species),   &
                  !&          '  und Nachfolger:',TRIM(RS(iR)%Product(iP)%Species)
                END IF 
              END DO

            END IF
          END DO

        END DO
      END DO
      DEALLOCATE(tsi)

      nnz_k = 0
      DO j = 1,len_c-1
        nnz_k = nnz_k + SIZE(R_k(k)%rIdx(j)%List)
      END DO
      write(*,*) '  alle Ketten aufgebaut => nun nichtnullen zählen: k,size: ',k, nnz_k

      ! build the matrices nue+ and nue-
      nue_plus(k)  = New_CSR( nspc, neq, nnz_k )
      nue_minus(k) = New_CSR( nspc, neq, nnz_k )

    END DO LOOP_OVER_ALL_CYCLES

    CALL ShowChain(R_k)

  END SUBROUTINE ISSA_nue

  SUBROUTINE ShowChain(C1)
    USE mo_control, ONLY: Chain
    TYPE(Chain) :: C1(:)
    INTEGER :: k, ii, iR, j
   
    write(*,*)
    write(*,*) REPEAT('-',50)//' Print the Reaction-Chains'

    DO k=1,SIZE(C1)
      write(*,*) '   Kette Nummer: ',k, ' der Länge: ',SIZE(C1(k)%sIdx,1)
      DO ii=1,SIZE(C1(k)%sName(:,:),1)
        write(*,*) '   [' , TRIM(C1(k)%sName(ii,1))//' , '//TRIM(C1(k)%sName(ii,2))//']  mit Reaktionen:  ' , & 
        &          '   [' , (  C1(k)%rIdx(ii)%List(j) ,j=1,SIZE(C1(k)%rIdx(ii)%List) ) ,']'
        !write(*,*) '   [' , C1(k)%sIdx(1),' , 'C1(k)%sIdx(2) , ']'
      END DO
      write(*,*)
    END DO


  END SUBROUTINE ShowChain


END MODULE issa
