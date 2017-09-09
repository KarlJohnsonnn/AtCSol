MODULE issa
  
  USE Kind_Mod
  USE Sparse_Mod, ONLY: SpRowIndColInd_T, CSR_Matrix_T

  IMPLICIT NONE

  ! OUT of ISSA_nue:
  !TYPE(SpRowIndColInd_T), ALLOCATABLE :: nue_p(:), nue_m(:)
  INTEGER, ALLOCATABLE :: S_imp(:), R_imp(:)    ! Sets of importent Species/Reactions
  INTEGER              :: len_S_imp=0, len_R_imp=0

  TYPE, EXTENDS(CSR_Matrix_T) :: CSR_Matrix_T_issa
    REAl(dp), ALLOCATABLE :: nvc(:)
  END TYPE CSR_Matrix_T_issa

  TYPE(CSR_Matrix_T), ALLOCATABLE :: nue_p(:), nue_m(:)
  !TYPE(CSR_Matrix_T_issa), ALLOCATABLE :: nue1_p(:), nue1_m(:)

  CONTAINS  


  SUBROUTINE ISSA_nue(nue,cyc,RS)
    USE mo_reac,     ONLY: neq, nspc, y_name
    USE mo_control,  ONLY: List, Chain, ZERO
    USE Sparse_Mod,  ONLY: CSR_Matrix_T, New_CSR, CSR_to_SpRowIndColInd,  &
    &                   SpRowIndColInd_to_CSR, PrintSparse2, PrintSparse3
    USE ChemSys_Mod, ONLY: ReactionStruct_T
    ! IN:
    TYPE(CSR_Matrix_T)      :: nue
    TYPE(List), ALLOCATABLE :: cyc(:)
    TYPE(ReactionStruct_T), ALLOCATABLE :: RS(:)
    ! TEMP:
    INTEGER :: i, ii, iR, iE, iP, j, jj,  k, n_c, len_c, nE, nP, nnz_k
    INTEGER,     ALLOCATABLE :: tsi(:), perm(:)
    INTEGER,     ALLOCATABLE :: ri_p(:), ci_p(:), ri_m(:), ci_m(:)
    REAL(dp),    ALLOCATABLE :: v_p(:), v_m(:)
    TYPE(Chain), ALLOCATABLE :: R_k(:)
    TYPE(SpRowIndColInd_T)   :: sp_nue

    ! Debugging Formats

    n_c = SIZE(cyc)

    write(*,*) '_____________________________________________'
    write(*,*)
    write(*,*) '  Testing ISSA stuff: building nue+ and nue- '
    write(*,*) '_____________________________________________'
    write(*,*)

    ! build nue_p and nue_m (see atmospheric env. 39 (2005) page: 4343)
    ALLOCATE( nue_p(0:n_c) , nue_m(0:n_c) , R_k(0:n_c) )

    sp_nue = CSR_to_SpRowIndColInd(nue)

    LOOP_OVER_ALL_CYCLES: DO k = 1,n_c ! Kette Nummer: k
      
      ! First we need to detect the reactions of each species cycle
      len_c = cyc(k)%len
      tsi   = cyc(k)%List

      ALLOCATE( ri_p(0), ci_p(0), v_p(0), ri_m(0), ci_m(0), v_m(0) )

      !ALLOCATE( R_k(k)%sName(len_c-1,2), R_k(k)%sIdx(len_c-1,2), R_k(k)%rIdx(len_c-1) )

      !write(*,*) REPEAT('-',32)//'  Browse Cycle: ',k,'  with length: ',len_c
      !write(*,'(5X,*(A))')  ( TRIM(y_name(cyc(k)%List(j)))//' -> ' , j=1,len_c-1 ) ,&
      !&                       TRIM(y_name(cyc(k)%List(len_c)))
      !write(*,*) 

      DO ii = 1,len_c-1
        !ALLOCATE( R_k(k)%rIdx(ii)%List(0), R_k(k)%rIdx(ii)%ListE(0), R_k(k)%rIdx(ii)%ListP(0) ) 
        DO iR = 1,neq

          nE = SIZE(RS(iR)%Educt)
          nP = SIZE(RS(iR)%Product)

          DO iE = 1,nE
            IF ( tsi(ii) == RS(iR)%Educt(iE)%iSpecies ) THEN
              !R_k(k)%sName(ii,1) = RS(iR)%Educt(iE)%Species
              !R_k(k)%sIdx(ii,1)  = RS(iR)%Educt(iE)%iSpecies
            
              DO iP = 1,nP
                IF ( tsi(ii+1) == RS(iR)%Product(iP)%iSpecies ) THEN
                  !R_k(k)%sName(ii,2)   = TRIM(RS(iR)%Product(iP)%Species)
                  !R_k(k)%sIdx(ii,2)    = RS(iR)%Product(iP)%iSpecies
                  !R_k(k)%rIdx(ii)%List = [R_k(k)%rIdx(ii)%List , iR]
                  !R_k(k)%rIdx(ii)%ListE = [-RS(iR)%Educt(iE)%Koeff,  R_k(k)%rIdx(ii)%ListE]
                  !R_k(k)%rIdx(ii)%ListP = [RS(iR)%Product(iP)%Koeff, R_k(k)%rIdx(ii)%ListP]
                  !write(*,'(A14,I0,A)',    ADVANCE='NO') '  Reaction :: ',iR,':  '//TRIM(RS(iR)%Line1)
                  !write(*,'(A9,F4.1,1X,A)',ADVANCE='NO') '  Educt: ', RS(iR)%Educt(iP)%Koeff,TRIM(RS(iR)%Educt(iE)%Species)
                  !write(*,'(A15,F4.1,1X,A)')             '  and Product: ', RS(iR)%Product(iP)%Koeff,TRIM(RS(iR)%Product(iP)%Species)

                  !write(*,'(A,I0,A3,2X,*(I0,5X))') '  spc in (beta-alpha)_(iR,j) = ',iR,' | ',nue%ColInd(nue%RowPtr(iR):nue%RowPtr(iR+1)-1)
                  !write(*,'(A,I0,A3,*(F4.1,2X))')  '  val in (beta-alpha)_(iR,j) = ',iR,' | ',nue%Val(nue%RowPtr(iR):nue%RowPtr(iR+1)-1)

                  DO jj = nue%RowPtr(iR),nue%RowPtr(iR+1)-1
                    IF (nue%ColInd(jj)==RS(iR)%Educt(iE)%iSpecies .OR. nue%ColInd(jj)==RS(iR)%Product(iP)%iSpecies) THEN
                      !write(*,'(A,I0,A1,I0,A4,F5.1)') '  spc in (beta-alpha)_(',iR,',',nue%ColInd(jj),') = ',nue%Val(jj)
                      IF ( nue%Val(jj) > ZERO ) THEN
                        ci_p = [ci_p, iR]
                        ri_p = [ri_p, nue%ColInd(jj)]
                        v_p  = [v_p, ABS(nue%Val(jj))]
                      ELSE
                        ci_m = [ci_m, iR]
                        ri_m = [ri_m, nue%ColInd(jj)]
                        v_m  = [v_m, ABS(nue%Val(jj))]
                      END IF
                    END IF
                  END DO

                END IF 
              END DO

            END IF
          END DO

        END DO
      END DO
      DEALLOCATE(tsi)

      CALL Progress3(k,n_c)

      ! build the sparse matrices nue_plus_ijk and nue_minus_ijk
      CALL SortVecAsc_I(ri_p,perm)
      ci_p = ci_p(perm)
      v_p  = v_p(perm)
      nnz_k = SIZE(ri_p)
      nue_p(k) = New_CSR(nue%n,nue%m,nnz_k,ri_p,ci_p,v_p)
      DEALLOCATE(ri_p,ci_p,v_p,perm)

      CALL SortVecAsc_I(ri_m,perm)
      ci_m = ci_m(perm)
      v_m  = v_m(perm)
      nnz_k = SIZE(ri_m)
      nue_m(k) = New_CSR(nue%n,nue%m,nnz_k,ri_m,ci_m,v_m)
      DEALLOCATE(ri_m,ci_m,v_m,perm)

      !CALL PrintSparse3(nue_p(k)%m, nue_p(k)%n, nue_p(k)%RowPtr, nue_p(k)%ColInd , nue_p(k)%Val , 'nue_p')
      !CALL PrintSparse3(nue_m(k)%m, nue_m(k)%n, nue_m(k)%RowPtr, nue_m(k)%ColInd , nue_m(k)%Val , 'nue_m')

      !write(*,*) '  Cycle: ',k, '  with nnz_p = ',nue_p(k)%nnz, '  and nnz_m = ',nue_m(k)%nnz

    END DO LOOP_OVER_ALL_CYCLES
    WRITE(*,*)

    !CALL ShowChain(R_k)

  END SUBROUTINE ISSA_nue

  SUBROUTINE ISSA_iter(Rate_t,Conc,t,h)

    USE Rates_Mod,  ONLY: ReactionRatesAndDerivative_ChemKin, ReactionRates_Tropos
    USE mo_control, ONLY: Teq, ZERO, List
    USE mo_reac,    ONLY: neq, nspc
    USE Sparse_Mod, ONLY: DAX_sparse

    ! IN:
    REAL(dp) :: Rate_t(neq)
    REAL(dp) :: Conc(nspc)
    REAL(dp) :: t, h

    ! TEMP:
    REAL(dp), DIMENSION(neq) :: Rate_th, DRatedT, avg_Rate
    INTEGER                  :: i, j, k, ii, jj, n_c, nnz, n_Simp
    REAL(dp)                 :: v, f_ijk, g_ijk
    INTEGER,    ALLOCATABLE  :: F_ik(:), G_ik(:)
    TYPE(List), ALLOCATABLE  :: tF(:), tG(:)
    INTEGER,    ALLOCATABLE  :: perm(:)

    REAL(dp)                 :: sum_p(nspc), sum_m(nspc)

    ! debugging format
    101 FORMAT( 3(A,I0), 2(A,Es10.3,2X) )
    102 FORMAT( A,I0, 2(A,I0,A,Es10.3,2X) )
    
    n_c       = SIZE(nue_p)-1 ! because nue_p index from 0 to n_c
    len_S_imp = SIZE(S_imp)

    IF (Teq) THEN
      CALL ReactionRatesAndDerivative_ChemKin( t+h , Conc , Rate_th , DRatedT )
    ELSE
      CALL ReactionRates_Tropos( t+h , Conc , Rate_th )
    END IF
    avg_Rate = ABS(Rate_th-Rate_t)/h

    !write(*,'(5X,A,I0,A,Es10.2)') ( 'avg_Rate(',j,') = ',avg_Rate(j) , j=1,neq ) 

    write(*,*)
    write(*,*) '----- Calculating normal. valu. coeffs -----'
    write(*,*)

    LOOP_OVER_ALL_CYCLES: DO k = 1,n_c
      
      !---------------------------------------------------------------------------------
      ! (a) "For the actual group of important species (index set S_imp) the valuation 
      !      coefficients f_ijk, g_ijk are calculated. At the start S_imp contains the 
      !      target species only." (See Atmosph. Env. 39 (2005) page 4343)
      !---------------------------------------------------------------------------------

      sum_p = DAX_sparse( nue_p(k) , avg_Rate )
      sum_m = DAX_sparse( nue_m(k) , avg_Rate )

      !DO i =1,nspc    ! debug test
      !  write(*,102) '  Cycle: ',k,'.  sum(nue_p_(',i,',:)*r(:)) = ', sum_p(i) ,' and  sum(nue_m_(',i,',:)*r(:)) = ', sum_m(i)
      !END DO
      !write(*,*)

      

!   ALT muss vlt komplett weg
!
!      CALC__f_ijk: DO jj = 1,nue_p(k)%nnz
!
!
!      ! ACHTUNG i muss wieder her!!!!
!
!        !i = nue_p(k)%RowInd(jj)
!        !ALLOCATE( tF(ii) )
!
!        DO ii = 1,n_Simp
!          IF ( S_imp(ii) == i ) THEN
!            j = nue_p(k)%ColInd(jj)
!            v = nue_p(k)%Val(jj)
!            f_ijk = v*avg_Rate(j) / sum_p
!            !write(*,101) '  Cycle: ',k,'.  nue_p_(',i,',',j,') = ', v,',   f = ',f_ijk
!            !nue_p(k)%nvc(jj) = f_ijk
!            !tF(ii)%ListE = [tF(ii)%ListE , f_ijk]
!          END IF
!        END DO
!        !CALL SortVecAsc_R( tF(ii)%ListE , perm )
!        !write(*,*) '  f_ijk_perm = ', ( tF(ii)%ListE , ii=1,SIZE(tf) )
!        !DEALLOCATE(tF)
!
!      END DO CALC__f_ijk
!      write(*,*)
!
!      DO jj = 1,nue_m(k)%nnz
!      ! ACHTUNG i muss wieder her!!!!
!        !i = nue_m(k)%RowInd(jj)
!        DO ii = 1,n_Simp
!          IF ( i == S_imp(ii) ) THEN
!            j = nue_m(k)%ColInd(jj)
!            v = nue_m(k)%Val(jj)
!            g_ijk = v*avg_Rate(j) / sum_m
!            !write(*,101) '  Cycle: ',k,'.  nue_m_(',i,',',j,') = ', v,',   g = ',g_ijk
!            !nue_m(k)%nvc(jj) = g_ijk
!          END IF
!        END DO
!      END DO
!      !write(*,*)

      !-------------------------------------------------------------------------------------
      ! (b) "The maximum member groups of redundant reactions (index sets F_ik, G_ik) with
      !      the property [SUM__(j in F_ik) f_ijk] < eps and [SUM__(j in G_ik) g_ijk] < eps
      !      are determined. Especially reactions with f_ijk = 0, g_ijk = 0 are always part 
      !      of F_ik, G_ik respectivley. The treshold value eps with 0 <= eps <= 1 controls
      !      the reduction intensity." (See Atmosph. Env. 39 (2005) page 4343)
      !-------------------------------------------------------------------------------------
      



    END DO LOOP_OVER_ALL_CYCLES

    stop 'issa_mod'
    
  END SUBROUTINE ISSA_iter

  SUBROUTINE SortVecAsc_R(v,q)
    USE mo_control, ONLY: big

    REAL(dp), INTENT(INOUT) :: v(:)
    INTEGER,  ALLOCATABLE, OPTIONAL :: q(:)
    !
    INTEGER :: i,n,im(1)
    REAL(dp), ALLOCATABLE :: tv(:)
   
    n = SIZE(v)
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

END MODULE issa
