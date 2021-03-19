MODULE Lumping_Mod
  
  USE Kind_Mod
  USE Sparse_Mod, ONLY: CSR_Matrix_T
  USE Control_Mod
  USE qsort_c_module
  USE Chemsys_Mod

  IMPLICIT NONE

! ##########################################################################
! ###                                                                    ###                  
! ###       Systematic reduction of complex tropospheric chemical        ###
! ###   mechanisms, Part II: Lumping using a time-scale based approach   ###
! ###                                                                    ###                  
! ##########################################################################
!
!  lumping technique according to Whitehouse 2004 
!  (best used for mechanisms big as MCM)
!
!  Whitehouse, L. E., Tomlin, A. S., and Pilling, M. J.: 
!  Systematic reduction of complex tropospheric chemical mechanisms, 
!  Part II: Lumping using a time-scale based approach, 
!  Atmos. Chem. Phys., 4, 2057–2081, 
!  https://doi.org/10.5194/acp-4-2057-2004, 2004. 
!
!

  TYPE linked_list
     TYPE(linked_list), POINTER   :: next => NULL()
     INTEGER                      :: id = -1
     CHARACTER(LenName)           :: id_char = ''
   CONTAINS
     PROCEDURE :: put  => put_ll
     PROCEDURE :: free => free_ll
  END TYPE linked_list

  TYPE lumping_group
     TYPE (lumping_group), POINTER    :: next => NULL()
     TYPE (linked_list), ALLOCATABLE  :: reactions(:)
     TYPE (linked_list)               :: spc
     INTEGER                          :: nspc = 1
     INTEGER                          :: id = -1
   CONTAINS
     PROCEDURE :: put  => put_lumping_group
     PROCEDURE :: free => free_lumping_group
     PROCEDURE :: add  => lumping_group_add
  END TYPE lumping_group

  TYPE NewReac_Data
     CHARACTER(LenLine)    :: Line3Template = ''
     CHARACTER(LenLine)    :: Line4 = ''
     TYPE(Duct_T), ALLOCATABLE  :: Educt(:)
     TYPE(Duct_T), ALLOCATABLE  :: Product(:)
     CHARACTER(LenType)    :: Type,  TypeConstant
  END TYPE NewReac_Data

  CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!                              ALGORITHM
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  !******************************
  !******************************
  !
  !  REMARKS:
  !
  ! currently ONLY gas phase, especially no equilibrium reactions (only reactions that occur in A and B)
  !
  !******************************
  !******************************
  
  SUBROUTINE lump_System(tau)
    ! input:
    !
    !  lifetimes of species at given times (tau(i,j) : lifetime of species i at timepoint j)
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: tau
    !
    ! end input


    !  lumped system matrices
    !TYPE(CSR_Matrix_T), INTENT(OUT) :: lumped_A, lumped_B
    !TYPE(CSR_Matrix_T) :: lumped_A, lumped_B
    TYPE(NewReac_Data), ALLOCATABLE :: l_ReactionSystem(:)

    ! lumping groups storage
    TYPE(lumping_group), TARGET  :: first_group 
    TYPE(lumping_group), POINTER :: current_group

    ! vector to test if species/reaction is already lumped
    LOGICAL, ALLOCATABLE :: Spc_lumped(:), reac_lumped(:)

    LOGICAL :: iSpc_lumpable, jSpc_lumpable, same_lifetimes, same_reactypes, equivs_present, similar_k
    INTEGER :: timestep, nSpc, nReac, iSpc, jSpc, n_iSpcReacs, n_jSpcReacs, i, j, nTau
    INTEGER, ALLOCATABLE :: iSpcReacs(:), jSpcReacs(:), reac_equivs(:), iSpc_reac_types(:), jSpc_reac_types(:), Reacs_Array(:)
    
    nReac = A%m
    nSpc  = A%n
    nTau  = SIZE(tau,2)
    
    ALLOCATE(reac_lumped(nReac),                &
           & Spc_lumped(nSpc))                   

    Spc_lumped  = .FALSE.
    Reac_lumped = .FALSE.

    current_group=>first_group

    DO iSpc=1,nSpc
      IF ( .NOT. Spc_lumped(iSpc) ) THEN
        Spc_lumped(iSpc) = .TRUE.

        ! find reactions in which iSpc occurs as educt
        CALL Spc_involving_reacs(iSpcReacs,n_iSpcReacs,iSpc)
     
        ! create new lumping group with iSpc as initial species
        ! every species will be in one lumping group (possibly containing only this species)
        CALL NextLumpingGroup(current_group,iSpc,iSpcReacs)

        ! check if one or more iSpcReacs are already lumped elsewhere -> species iSpc cannot be lumped!
        iSpc_lumpable = .TRUE.
        IF ( n_iSpcReacs == 0 ) THEN
          iSpc_lumpable = .FALSE. ! don't lump what exists only as product (or do??)
        ELSE
          DO i=1,n_iSpcReacs
            IF ( reac_lumped(iSpcReacs(i)) ) THEN
              iSpc_lumpable=.FALSE.
            END IF
          END DO
        END IF

        IF ( iSpc_lumpable ) THEN

          ! declare all reactions involving iSpc as educt as lumped
          DO i=1,n_iSpcReacs
            reac_lumped(iSpcReacs(i))=.TRUE.
          END DO

          ! find reaction types of species i (occuring numbers of educts of iSpcReacs)
          CALL find_reaction_types(iSpc_reac_types,iSpcReacs)

          ! check all species with greater indices than iSpc, 
          ! smaller ones have already been checked before
          DO jSpc=iSpc+1,nSpc
            same_lifetimes=.TRUE.
            DO timestep=1,nTau
            !~~~~~~~~~~~~~~~~~~~~~
            ! maybe vectorize this loop if faster
            ! 1/(machine_eps+J_i,i) in tau calculation?
            ! think where eps_tau could come from, pass as argument, create constant?
            !~~~~~~~~~~~~~~~~~~~~~
              IF (ABS(tau(iSpc,timestep)-tau(jSpc,timestep))>eps_tau) THEN
                same_lifetimes=.FALSE.
              END IF
            END DO  

            IF (same_lifetimes) THEN
              ! find reactions in which iSpc occurs as educt
              CALL Spc_involving_reacs(jSpcReacs,n_jSpcReacs,jSpc)
              ! check if one or more jSpcReacs are already lumped elsewhere -> species jSpc cannot be lumped!
              jSpc_lumpable = .TRUE.
              IF ( n_jSpcReacs == 0 ) THEN
                jSpc_lumpable = .FALSE.
              ELSE
                DO j=1,n_jSpcReacs
                  IF ( reac_lumped(jSpcReacs(j)) ) THEN
                    jSpc_lumpable=.FALSE.
                  END IF
                END DO
              END IF
              
              IF (jSpc_lumpable) THEN
                ! find reaction types of species i (occuring numbers of educts of iSpcReacs)
                CALL find_reaction_types(jSpc_reac_types,jSpcReacs)

                ! check for equal reaction types of iSpc and jSpc
                CALL QsortC_int(iSpc_reac_types)
                CALL QsortC_int(jSpc_reac_types)
                same_reactypes=.FALSE.
                IF (ALL(iSpc_reac_types == jSpc_reac_types)) THEN
                  same_reactypes=.TRUE.
                END IF

                ! check for same size too (so no reactions are left over)
                IF (same_reactypes .AND. n_iSpcReacs==n_jSpcReacs) THEN
                  ALLOCATE(reac_equivs(n_jSpcReacs))
                  equivs_present=.TRUE.
                  ! now check for reaction equivalents (reactions with same reactants excluding iSpc and jSpc)
                  CALL check_reac_equivs(equivs_present,reac_equivs,iSpcReacs,jSpcReacs,iSpc,jSpc)
                  IF (equivs_present) THEN
                          !TO BE IMPLEMENTED: check rate constant errors
                    similar_k=.TRUE.
                    CALL compare_k(similar_k,reac_equivs,iSpcReacs,jSpcReacs)
                    IF (similar_k) THEN
                      CALL current_group%add(jSpc,reac_equivs,jSpcReacs)
                      ! declare all jSpcReacs as lumped
                      DO j=1,n_jSpcReacs
                        reac_lumped(jSpcReacs(j))=.TRUE.
                      END DO
                      Spc_lumped(jSpc) = .TRUE.
                    END IF
                  END IF ! equivs_present
                  DEALLOCATE(reac_equivs)
                END IF ! same_reactypes and number of reacs
                DEALLOCATE(jSpc_reac_types)
              END IF ! jSpc_lumpable
              DEALLOCATE(jSpcReacs)
            END IF ! same_lifetimes
          END DO ! jSpc loop
          ! de-lump reactions if no species was lumped with iSpc
          IF ( current_group%nspc == 1 ) THEN
            DO i=1,n_iSpcReacs
              reac_lumped(iSpcReacs(i))=.FALSE.
            END DO
          END IF
          DEALLOCATE(iSpc_reac_types)
        END IF ! iSpc_lumpable
        DEALLOCATE(iSpcReacs)
      END IF ! iSpc not lumped
    END DO ! iSpc loop
    CALL WriteLumpingGroups(first_group)

    CALL build_LumpedReacSys(l_ReactionSystem, first_group, reac_lumped)

    CALL Print_LumpedSysFile(l_ReactionSystem, first_group, 'LUMPING/'//TRIM(BSP)//'_lumped.sys')

  END SUBROUTINE lump_System



! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!                       ALGORITHM-SUBROUTINES
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE Print_LumpedSysFile(l_RS, first_group, FileName)
    ! input:
    !
    TYPE(NewReac_Data), DIMENSION(:), INTENT(IN) :: l_RS
    !
    TYPE(lumping_group), TARGET, INTENT(IN)      :: first_group
    !
    CHARACTER(*), INTENT(IN)                     :: FileName
    !
    ! end input

    INTEGER :: i, j, nR, nEducts, nProducts
    REAL(dp) :: Koeff
    CHARACTER(lenName) :: Species, Koeff_Str, len_Str
    TYPE(lumping_group), POINTER :: group
    TYPE(linked_list), POINTER :: Spc_list
    
    group=>first_group
    
    OPEN(UNIT=989,FILE=ADJUSTL(TRIM(FileName)),STATUS='UNKNOWN')
  
    WRITE(989,'(A)') '# ================= '//TRIM(FileName)//' ================='
    WRITE(989,'(A)') '# = Please copy the data into your sys-file for ='
    WRITE(989,'(A)') '# =============== chemical input. ==============='
    WRITE(989,'(A)') '#'
    WRITE(989,'(A)') '#  ===================   Unit options   ======================'
    WRITE(989,'(A)') ''
    WRITE(989,'(A)') 'UNIT GAS    0   #    Gas phase units     (0 = molec/cm3, 1 = mol/m3)'
    WRITE(989,'(A)') 'UNIT AQUA   0   #    Aqueous phase units (0 = mol/l)'
    WRITE(989,'(A)') ''
    WRITE(989,'(A)') '#'
    WRITE(989,'(A)') '# This is a lumped system. Several original species were lumped.'
    WRITE(989,'(A)') '# The reactions below consist of these lumped species.' 
    WRITE(989,'(A)') '#'
    WRITE(989,'(A)') '# The lumped species have names and ingredients of the original system as follows:'
    DO WHILE ( ASSOCIATED(group) )
      IF ( group%nspc > 1 ) THEN
        WRITE(989,'(A)') '# '
        WRITE(989,'(A)') '# '//ADJUSTL(group%spc%id_char)
        WRITE(989,'(A)') '# contains '

        Spc_list => group%spc
        WRITE(989,'(A)', ADVANCE ='NO') '# '//TRIM(ADJUSTL(Spc_list%id_char))
        DO WHILE ( ASSOCIATED(Spc_list%next) )
          Spc_list => Spc_list%next
          WRITE(989,'(A)', ADVANCE='NO') ', '//TRIM(ADJUSTL(Spc_list%id_char))
        END DO
        WRITE(989,'(A)') ''
      END IF
      group=>group%next
    END DO
    WRITE(989,'(A)') '#'
    WRITE(989,'(A)') '#'
    WRITE(989,'(A)') '#'



    nR = SIZE(l_RS)
    DO i=1,nR
      WRITE (989,*)
      
      WRITE (989,'(A)') 'CLASS: '//TRIM(l_RS(i)%Type)

      nEducts   = SIZE(l_RS(i)%Educt)
      nProducts = SIZE(l_RS(i)%Product)
      DO j=1,nEducts
        Koeff = l_RS(i)%Educt(j)%Koeff
        Species = l_RS(i)%Educt(j)%Species
        IF ( ANINT(Koeff)-Koeff < milli ) THEN
          WRITE (Koeff_Str,'(F4.1)') Koeff
        ELSE
          WRITE (Koeff_Str,'(F5.3)') Koeff
        END IF
        WRITE (len_Str,'(I3)') LEN_TRIM(ADJUSTL(Koeff_Str))+LEN_TRIM(Species)+1
        WRITE (989,'(A'//TRIM(len_Str)//')',ADVANCE='NO') TRIM(ADJUSTL(Koeff_Str))//' '//TRIM(Species)
        IF ( j<nEducts ) THEN
          WRITE (989,'(A3)',ADVANCE='NO') ' + '
        ELSE 
          WRITE (989,'(A3)',ADVANCE='NO') ' = '
        END IF
      END DO
      
      DO j=1,nProducts
        Koeff = l_RS(i)%Product(j)%Koeff
        Species = l_RS(i)%Product(j)%Species
        IF ( ANINT(Koeff)-Koeff < milli ) THEN
          WRITE (Koeff_Str,'(F4.1)') Koeff
        ELSE
          WRITE (Koeff_Str,'(F5.3)') Koeff
        END IF
        WRITE (len_Str,*) LEN_TRIM(ADJUSTL(Koeff_Str))+LEN_TRIM(Species)+1
        WRITE (989,'(A'//TRIM(len_Str)//')',ADVANCE='NO') TRIM(ADJUSTL(Koeff_Str))//' '//TRIM(Species)
        IF ( j<nProducts ) THEN
          WRITE (989,'(A3)',ADVANCE='NO') ' + '
        END IF
      END DO
    
      WRITE (989,*) 
      WRITE (989,'(A)') TRIM(l_RS(i)%Line3Template)
      WRITE (989,'(A)') TRIM(l_RS(i)%Line4)
    END DO

    CLOSE(989)    
  END SUBROUTINE Print_LumpedSysFile


  SUBROUTINE build_LumpedReacSys(l_RS, first_group, reac_lumped)
    ! input:
    !
    TYPE(NewReac_Data), ALLOCATABLE :: l_RS(:)
    !
    TYPE(lumping_group), TARGET :: first_group
    !
    LOGICAL, DIMENSION(:) :: reac_lumped
    !
    ! end input

    CHARACTER(LenName), ALLOCATABLE :: Spc_Old2New(:)
    TYPE(lumping_group), POINTER :: group
    INTEGER :: i, j, k, l, iR, kR, kE, NewReac, nLumpedReacs, nRemainingReacs, nNewReacs,          &
             & nLumpedSpc, nRemainingSpc, nNewSpc, nEducts, nProducts, counter
    
    INTEGER, ALLOCATABLE :: iSpc_Old2New(:), Spc_array(:), Reacs_Array(:)

    group=>first_group

    ! calculate numbers of lumped/old/new species/reactions
    nRemainingReacs = COUNT(.NOT. reac_lumped)
    IF ( group%nspc == 1 ) THEN
      nRemainingSpc = 1
      nLumpedSpc    = 0
      nLumpedReacs  = 0
    ELSE
      nRemainingSpc = 0
      nLumpedSpc    = 1
      nLumpedReacs  = SIZE(group%reactions)
    END IF
    DO WHILE ( ASSOCIATED(group%next) )
      IF ( group%nspc == 1 ) THEN
        nRemainingSpc = nRemainingSpc + 1
      ELSE
        nLumpedSpc    = nLumpedSpc + 1
        nLumpedReacs  = nLumpedReacs + SIZE(group%reactions)
      END IF
      group=>group%next
    END DO
    nNewSpc = nRemainingSpc + nLumpedSpc
    nNewReacs = nRemainingReacs + nLumpedReacs

    ALLOCATE( Spc_Old2New(A%n+SIZE(ListNonReac2)), iSpc_Old2New(A%n+SIZE(ListNonReac2)), l_RS(nNewReacs) )
    
    ! create mask: Spc_Old2New(old_species_id) = new_species_id
    Spc_Old2New = ''
    group=>first_group
    DO i=1,nNewSpc
      CALL LL2Array(group%spc,Spc_array)
      DO j=1,SIZE(Spc_array)
        iSpc_Old2New(Spc_array(j))=i
      END DO
      IF ( group%nspc > 1 ) THEN
        CALL NewNames(group)
      END IF
      group=>group%next
    END DO
    ! add passive species
    DO i=1,SIZE(ListNonReac2)
      iSpc_Old2New(nNewSpc+i) = nNewSpc + i
    END DO
    
    ! collect all lumped reactions
    NewReac = 1
    group=>first_group
    DO i=1,nNewSpc
      IF ( group%nspc > 1 ) THEN
        DO j=1,SIZE(group%reactions)
          CALL LL2Array(group%reactions(j),Reacs_Array)
          CALL Reac2NewReac(Reacs_array)
        END DO
      END IF
      group=>group%next
    END DO
    ! collect all remaining reactions
    DO i=1,SIZE(reac_lumped)
      IF ( .NOT. reac_lumped(i) ) THEN
        CALL Reac2NewReac((/i/))
      END IF
    END DO

   CONTAINS

     SUBROUTINE Reac2NewReac(Reacs_old)
       INTEGER, DIMENSION(:) :: Reacs_old

       INTEGER :: iR_old
       TYPE(Duct_T) :: current

       iR_old = Reacs_old(1)
       
       ! calculate number of educts and products
       nEducts = SIZE(ReactionSystem(iR_old)%Educt)
       nProducts = 0
       DO k=1,SIZE(Reacs_old)
         kR = Reacs_old(k)
         nProducts = nProducts + SIZE(ReactionSystem(kR)%Product)
       END DO
          
       ALLOCATE( l_RS(NewReac)%Educt(nEducts) ) 
       ALLOCATE( l_RS(NewReac)%Product(nProducts) ) 

       ! write new educts and products, 
       ! lumped species gets the name of the first species in lumping group
       ! educts
       DO k=1,nEducts
         current = ReactionSystem(iR_old)%Educt(k)
         l_RS(NewReac)%Educt(k)          = current
         ! correct id to new species id
         l_RS(NewReac)%Educt(k)%iSpecies = iSpc_Old2New(current%iSpecies)
       END DO
       ! products
       counter = 1
       DO k=1,SIZE(Reacs_old)
         kR=Reacs_old(k)
         DO l=LBOUND(ReactionSystem(kR)%Product,1),UBOUND(ReactionSystem(kR)%Product,1)
           current = ReactionSystem(kR)%Product(l)
           l_RS(NewReac)%Product(counter) = current
           ! concern name of species
           IF ( Spc_Old2New(current%iSpecies) /= '' ) THEN ! duct has a new name
             l_RS(NewReac)%Product(counter)%Species  = Spc_Old2New(current%iSpecies)
           ELSE ! old name is maintained
             l_RS(NewReac)%Product(counter)%Species = current%Species
           END IF
           ! correct id to new species id
           l_RS(NewReac)%Product(counter)%iSpecies = iSpc_Old2New(current%iSpecies)
           counter=counter+1
         END DO
       END DO

       ! collect remaining properties
       l_RS(NewReac)%Type           = ReactionSystem(iR_old)%Type
       l_RS(NewReac)%TypeConstant   = ReactionSystem(iR_old)%TypeConstant
       l_RS(NewReac)%Line3Template  = ReactionSystem(iR_old)%Line3
       l_RS(NewReac)%Line4          = ReactionSystem(iR_old)%Line4

       NewReac = NewReac + 1

     END SUBROUTINE

     SUBROUTINE NewNames(group)
       TYPE(lumping_group), POINTER, INTENT(IN) :: group
       
       INTEGER, ALLOCATABLE :: Spc_array(:)
       INTEGER :: i, iR, iSpc
       CHARACTER(lenName) :: NewName

       CALL LL2Array(group%spc,Spc_array)
       
       ! the new lumped species gets the name of the first species in the lumping group
       iSpc = Spc_array(1)

       ! set iR to a reaction involving the first spc
       iR = group%reactions(1)%id

       DO i=1,SIZE(ReactionSystem(iR)%Educt)
         IF ( ReactionSystem(iR)%Educt(i)%iSpecies == iSpc ) THEN
           NewName = ReactionSystem(iR)%Educt(i)%Species
         END IF
       END DO

       DO i=2,SIZE(Spc_array)
         Spc_Old2New(Spc_array(i)) = NewName 
       END DO

     END SUBROUTINE NewNames

  END SUBROUTINE build_LumpedReacSys


  SUBROUTINE compare_k(similar_k, reac_equivs, iSpcReacs, jSpcReacs)
    ! input:
    !
    LOGICAL, INTENT(INOUT) :: similar_k
    !
    INTEGER, DIMENSION(:), INTENT(IN) :: reac_equivs, iSpcReacs, jSpcReacs
    !
    ! end input
 
    INTEGER :: i, reac, iR, jR, nReacs
    TYPE(ReactionStruct_T) :: iReac, jReac

    nReacs = SIZE(iSpcReacs)
    DO reac=1,nReacs
      iR = iSpcReacs(reac)
      jR = jSpcReacs(reac)

      iReac = ReactionSystem(iR)
      jReac = ReactionSystem(jR)
      
      IF (iReac%TypeConstant == jReac%TypeConstant .AND. iReac%Type == jReac%Type ) THEN
        DO i=1,SIZE(iReac%Constants)
          IF ( .NOT. ABS( ( iReac%Constants(i)-jReac%Constants(i) ) / ( iReac%Constants(i)+eps ) ) < eps_k ) THEN
            similar_k = .FALSE.
          END IF
        END DO
      ELSE
        similar_k = .FALSE.
      END IF

      IF (.NOT. similar_k) EXIT
    END DO

  END SUBROUTINE compare_k


 ! SUBROUTINE build_LumpedSys(lumped_A, lumped_B, first_group, Spc_lumped, reac_lumped)
 !   ! input:
 !   !
 !   TYPE(CSR_Matrix_T) :: lumped_A, lumped_B
 !   !
 !   TYPE(lumping_group), TARGET, INTENT(IN) :: first_group
 !   !
 !   LOGICAL, ALLOCATABLE, INTENT(IN) :: Spc_lumped(:), reac_lumped(:)
 !   !
 !   ! end input
!
!    TYPE(Species_T), ALLOCATABLE, TARGET :: NewListGas2(:)
!    INTEGER, ALLOCATABLE :: lA_RowPtr(:), lA_ColInd(:), lA_val(:),                      &
!                          & lB_RowPtr(:), lB_ColInd(:), lB_val(:),                      &
!                          & Spc_Old2New(:), remainingReacs(:), groupSpcs(:), eqReacs(:)
!    TYPE(lumping_group), POINTER :: group
!    INTEGER :: i, j, k, l, iR, counter, nLumpedSpc, nRemainingSpc, nLumpedReacs, nRemainingReacs,  &
!             & nNewSpc, nNewReacs, nLumpedEducts, nLumpedProducts, newReac, lA_Col, lB_Col,        &
!             & nEducts_iR, nProducts_iR, educt_iR, product_iR
!
!    group => first_group
!
!    IF ( COUNT(.NOT. Spc_lumped) > 0 ) THEN
!      WRITE (*,*) 'ERROR: some species were skipped and not put into a lumping group'
!      STOP
!    END IF
!
!    ! count lumping species
!    nNewSpc = 1 ! first group
!    DO WHILE ( ASSOCIATED(group%next) )
!      nNewSpc = nNewSpc + 1
!      group => group%next
!    END DO
!    group => first_group
!
!    nLumpedReacs = 0
!    nLumpedEducts = 0
!    nLumpedProducts = 0
!    ! count lumped reactions and educts
!    DO i=1,nNewSpc
!      nLumpedReacs = nLumpedReacs + SIZE(group%reactions)
!      DO j=1,SIZE(group%reactions)
!      ! only count educts if this reaction will be a lumped one or will remain unmodified in the new system
!      ! (and not if its lumped elsewhere)
!        IF ( (group%nspc>1) .OR. (.NOT. reac_lumped(group%reactions(j)%id) ) ) THEN
!          nLumpedEducts = nLumpedEducts + A%RowPtr(group%reactions(j)%id+1) - A%RowPtr(group%reactions(j)%id)
!        END IF
!      END DO
!
 !     IF ( ASSOCIATED(group%next) ) THEN
 !       group=>group%next
 !     END IF
 !   END DO
 !   group => first_group
 !   ! count remaining reactions and educts
 !   nRemainingReacs = COUNT(.NOT. reac_lumped)
 !   ALLOCATE(remainingReacs(nRemainingReacs))
 !   counter = 1
 !   DO i=1,A%m
 !     IF ( .NOT. reac_lumped(i) ) THEN
 !       remainingReacs(counter) = i
 !       counter = counter + 1
 !       nLumpedEducts = nLumpedEducts + A%RowPtr(i+1) - A%RowPtr(i)
 !     END IF
 !   END DO
!
 !   nNewReacs = nRemainingReacs + nLumpedReacs
 !   ALLOCATE(NewListGas2(nNewSpc),                                                    &
 !          & lA_RowPtr(nNewReacs+1), lA_ColInd(nLumpedEducts), lA_val(nLumpedEducts), &
 !          & lB_RowPtr(nNewReacs+1), lB_ColInd(B%nnz), lB_val(B%nnz))
!
 !   ALLOCATE(Spc_Old2New(A%n))
!
 !   ! old_id now participates in new_id  ( Spc_Old2New(old_id) = new_id )
 !   Spc_Old2New = 0
 !   DO i=1,nNewSpc
 !     CALL LL2Array(group%spc, groupSpcs)
 !     DO j=1,SIZE(groupSpcs)
 !       Spc_Old2New(groupSpcs(j)) = i
 !     END DO
 !     IF ( ASSOCIATED(group%next) ) THEN
 !       group=>group%next
 !     END IF
 !   END DO
 !   group => first_group
!
 !   ! initialize counters for new reactions (newReac) and ColInd/val for A and B
 !   newReac = 1
 !   lA_Col  = 1
 !   lB_Col  = 1
!
 !   ! build lumped_A vectors
 !   counter = 1
 !   DO i=1,nNewSpc
 !     IF ( group%nspc==1 ) THEN ! add this one as one of the new species
 !       NewListGas2(i) = ListGas2(group%spc%id)
 !       DO j=1,SIZE(group%reactions)
 !         iR = group%reactions(j)%id
 !         ! collect reactions that are not lumped elsewhere to maintain them in the new system
 !         IF (.NOT. reac_lumped(iR)) THEN
 !           
 !           ! fill A vectors
!
 !           nEducts_iR = A%RowPtr(iR+1)-A%RowPtr(iR)
 !           lA_RowPtr(newReac+1)=lA_RowPtr(newReac)+nEducts_iR
 !           DO k=0,nEducts_iR-1
 !             lA_ColInd(lA_Col) = Spc_Old2New(A%ColInd(A%RowPtr(iR)+k))
 !             lA_val(lA_Col) = A%val(A%RowPtr(iR)+k)
 !             lA_Col = lA_Col + 1
 !           END DO
!
 !           ! fill B vectors
!
 !           nProducts_iR = B%RowPtr(iR+1)-B%RowPtr(iR)
 !           lB_RowPtr(newReac+1)=lB_RowPtr(newReac)+nProducts_iR
 !           DO k=0,nProducts_iR-1
 !               !~~~~~~~~~~~~~~~~~~~~~~~~
 !               !
 !               ! some values may occur multiple times (2 species that got lumped together are in the products of two iRs 
 !               ! -> same new id, will occur twice in ColInd, is that a problem or does it just get added?)
 !               !
 !               !~~~~~~~~~~~~~~~~~~~~~~~~
 !             lB_ColInd(lB_Col) = Spc_Old2New(B%ColInd(B%RowPtr(iR)+k))
 !             lB_val(lB_Col) = B%val(B%RowPtr(iR)+k)
 !             lB_Col = lB_Col + 1
 !           END DO
!
 !           newReac = newReac + 1
 !         END IF
 !       END DO
 !     ELSE ! create new lumped species
 !       NewListGas2(i) = ListGas2(group%spc%id) ! same properties in Species_T assumed for all members of the lumping group
 !       NewListGas2(i)%Species = 'Lu'//NewListGas2(i)%Species  ! rename the new species
 !       counter = counter + 1
!
 !       DO j=1,SIZE(group%reactions)
 !         iR = group%reactions(j)%id
 !         
 !         ! fill A vectors
 !         
 !         nEducts_iR = A%RowPtr(iR+1)-A%RowPtr(iR)
 !         lA_RowPtr(newReac+1)=lA_RowPtr(newReac)+nEducts_iR
 !         DO k=0,nEducts_iR-1
 !           educt_iR = A%ColInd(A%RowPtr(iR)+k)
 !           
 !           lA_ColInd(lA_Col) = Spc_Old2New(educt_iR)
 !           lA_val(lA_Col) = A%val(A%RowPtr(iR)+k) ! lumped Spc gets same stoichometric coefficient as old species as they are all
 !                                                  ! the same (or 1 ??)
 !           
 !           lA_Col = lA_Col + 1
 !         END DO
!
 !         ! fill B vectors
!
 !         CALL LL2Array(group%reactions(j), eqReacs)
 !         lB_RowPtr(newReac+1) = lB_RowPtr(newReac)
 !         DO k=1,SIZE(eqReacs)
 !           iR = eqReacs(k)
 !           nProducts_iR = B%RowPtr(iR+1)-B%RowPtr(iR)
 !           lB_RowPtr(newReac+1)=lB_RowPtr(newReac+1)+nProducts_iR
 !           DO l=0,nProducts_iR-1
 !             ! product_iR is the new id for the old product
 !             product_iR = Spc_Old2New(B%ColInd(B%RowPtr(iR)+l))
 !             !~~~~~~~~~~~~~~~~~~~~~~~~
 !             !
 !             ! if product_iR == i another species of the lumping group occured as product
 !             ! -> problematic???
 !             !
 !             ! EVEN if product_iR /= i product_iR could be a newly lumped species and so
 !             ! for every member of a lumping group the reactions that produce it
 !             ! need to be known (always, so maybe collect them somewhere, could be extracted from B)
 !             !
 !             !~~~~~~~~~~~~~~~~~~~~~~~~
 !             
 !             lB_ColInd(lB_Col) = product_iR
 !             lB_val(lB_Col) = B%val(B%RowPtr(iR)+l) !THIS IS SIGMA 
 !             
 !             lB_Col = lB_Col + 1
 !           END DO
 !         END DO
!
 !         newReac = newReac + 1
 !       END DO
 !     END IF
 !     IF ( ASSOCIATED(group%next) ) THEN
 !       group=>group%next
 !     END IF
 !   END DO
!
 !   !build lumped_A matrix
 !   lumped_A = New_CSR(nNewReacs,nNewSpc)
 !   ALLOCATE(lumped_A%ColInd(nLumpedEducts), lumped_A%val(nLumpedEducts))
 !   lumped_A%RowPtr = lA_RowPtr
 !   lumped_A%ColInd = lA_ColInd
 !   lumped_A%val = lA_val
!
 !   !build lumped_B matrix
 !   lumped_B = New_CSR(nNewReacs,nNewSpc)
 !   ALLOCATE(lumped_B%ColInd(B%nnz), lumped_B%val(B%nnz))
 !   lumped_B%RowPtr = lB_RowPtr
 !   lumped_B%ColInd = lB_ColInd
 !   lumped_B%val = lB_val
!
 ! END SUBROUTINE build_LumpedSys


  SUBROUTINE Spc_involving_reacs(SpcReacs,nSpcReacs,Spc)
    ! input:
    !
    !  species ID
    INTEGER, INTENT(IN) :: Spc
    !
    ! end input

    ! output:
    !
    !  indices of reactions that involve Spc as educt
    INTEGER, ALLOCATABLE, INTENT(OUT) :: SpcReacs(:)
    !
    !  number of reactions (=size(SpcReacs))
    INTEGER :: nSpcReacs
    !
    ! end output


    INTEGER :: iReac, nReac, SpcReac
        
    nReac=A%m
    nSpcReacs=0

    ! find number of reactions in which iSpc participates
    ! to be able to allocate vector of iSpc involving reactions
    DO iReac=1,nReac
      ! test if iSpc occurs in educts of reaction iReac
      IF ( ANY( A%ColInd(A%RowPtr(iReac):A%RowPtr(iReac+1)-1) == Spc ) ) THEN
        nSpcReacs=nSpcReacs+1
      END IF
    END DO

    IF (nSpcReacs > 0) THEN
      ALLOCATE(SpcReacs(nSpcReacs))
      ! find reactions in which iSpc participates
      SpcReac=1
      DO iReac=1,nReac
        ! test if iSpc occurs in educts of reaction iReac
        IF( ANY( A%ColInd(A%RowPtr(iReac):A%RowPtr(iReac+1)-1) == Spc ) ) THEN
          SpcReacs(SpcReac)=iReac
          SpcReac=SpcReac+1
        END IF
      END DO
    ELSE
      ALLOCATE(SpcReacs(1))
      SpcReacs(1)=0
    END IF

  END SUBROUTINE Spc_involving_reacs



  SUBROUTINE find_reaction_types(types,SpcReacs) 
    ! input:
    !
    INTEGER, DIMENSION(:), INTENT(IN) :: SpcReacs
    !
    ! end input
    
    ! output:
    INTEGER, ALLOCATABLE, INTENT(OUT) :: types(:)
    ! end output
        
    LOGICAL, ALLOCATABLE :: types_bool(:)
    INTEGER :: iT, iR, iR_type, nSpcReacs, maxtype, n_types
       
    ! first find maximum reaction type 
    nSpcReacs=SIZE(SpcReacs,1)
    DO iR=1,nSpcReacs
      ! check how many educts are present in reaction iR
      iR_type=A%RowPtr(iR+1)-A%RowPtr(iR)
      ! check if type greater than maxtype, if so change maxtype
      IF (iR_type>maxtype) THEN
        maxtype=iR_type
      END IF
    END DO

    ALLOCATE(types_bool(maxtype))
    types_bool = .FALSE.
    
    ! then collect types and count them 
    n_types=0
    DO iR=1,nSpcReacs
      ! check how many educts are present in reaction iR
      iR_type=A%RowPtr(iR+1)-A%RowPtr(iR)
      IF (types_bool(iR_type) .EQV. .FALSE.) THEN
        ! set types-entry to true for this type if not already
        types_bool(iR_type)=.TRUE.
        n_types=n_types+1
      END IF
    END DO
    ALLOCATE(types(n_types))
    
    ! finally save types in a vector simply containing all types
    iT=0
    DO iR=1,maxtype
      IF (types_bool(iR)) THEN
        iT=iT+1
        types(iT)=iR
      END IF
    END DO

  END SUBROUTINE find_reaction_types

  SUBROUTINE check_reac_equivs(equivs_present,reac_equivs,iSpcReacs,jSpcReacs,iSpc,jSpc) 
    ! input:
    !
    INTEGER, DIMENSION(:), INTENT(IN) :: iSpcReacs, jSpcReacs
    INTEGER, INTENT(IN) :: iSpc, jSpc
    !
    ! end input
    
    ! output:
    !
    !  logical variable to see if every iSpcReaction has an equivalent jSpcReaction
    !  (meaning educts are the same excluding iSpc and jSpc)
    LOGICAL, INTENT(OUT) :: equivs_present
    !
    INTEGER, ALLOCATABLE, INTENT(OUT) :: reac_equivs(:)
    ! end output
        
    INTEGER :: nSpcReacs, iR, jR, iE, jE, iSpcReac, jSpcReac, i_educt, j_educt
    INTEGER, ALLOCATABLE :: i_educts(:), i_eductsCoeff(:), j_educts(:), j_eductsCoeff(:)
    LOGICAL :: found_iR_equiv, educts_equal
    LOGICAL, DIMENSION(SIZE(jSpcReacs)) :: reac_available
        
    nSpcReacs=SIZE(jSpcReacs)
    ALLOCATE(reac_equivs(nSpcReacs))
    reac_equivs=ZERO
    reac_available=.TRUE.
    equivs_present=.TRUE.
    DO iR=1,nSpcReacs
      found_iR_equiv=.FALSE.  

      iSpcReac=iSpcReacs(iR)
      
      ! allocate number of educts for iSpcReac(iR)
      ALLOCATE(i_educts(A%RowPtr(iSpcReac+1)-A%RowPtr(iSpcReac)))
      ALLOCATE(i_eductsCoeff(A%RowPtr(iSpcReac+1)-A%RowPtr(iSpcReac)))
      
      ! determine educts of reaction iSpcReac (excluding iSpc)
      DO iE=1,SIZE(i_educts)
        i_educt=A%ColInd(A%RowPtr(iSpcReac)+iE-1)
     
        i_eductsCoeff(iE) = A%val(A%RowPtr(iSpcReac)+iE-1) 
        IF (i_educt==iSpc) THEN
          i_educts(iE)=0
        ELSE 
          i_educts(iE)=i_educt
        END IF
      END DO
      
      ! search for an equivalent reaction to iSpcReac
      DO jR=1,nSpcReacs
        IF (reac_available(jR) .AND. (.NOT. found_iR_equiv)) THEN
          jSpcReac=jSpcReacs(jR)
    
          ! allocate number of educts for jSpcReac(jR)
          ALLOCATE(j_educts(A%RowPtr(jSpcReac+1)-A%RowPtr(jSpcReac)))
          ALLOCATE(j_eductsCoeff(A%RowPtr(jSpcReac+1)-A%RowPtr(jSpcReac)))
        
          ! determine educts of reaction jSpcReac (excluding jSpc)
          DO jE=1,SIZE(j_educts)
            j_educt=A%ColInd(A%RowPtr(jSpcReac)+jE-1)

            j_eductsCoeff(jE) = A%val(A%RowPtr(jSpcReac)+jE-1) 
            IF (j_educt==jSpc) THEN
              j_educts(jE)=0
            ELSE 
              j_educts(jE)=j_educt
            END IF
          END DO

          ! compare educts and stoichometric coefficients
          educts_equal = .FALSE.
          IF ( SIZE(i_educts) == SIZE(j_educts) ) THEN
            DO jE=1,SIZE(j_educts)
              educts_equal = .FALSE.
              DO iE=1,SIZE(i_educts)
                IF ( (j_educts(jE) == i_educts(iE)) .AND. (j_eductsCoeff(jE) == i_eductsCoeff(iE)) ) THEN
                  educts_equal = .TRUE.
                END IF
              END DO
              IF (.NOT. educts_equal) THEN
                EXIT 
              END IF
            END DO
          END IF

          IF (educts_equal) THEN
            ! if same, collect the equivalent jR to iR, block jR for next tests 
            ! and EXIT searching for an iR equivalent
            reac_equivs(iR)=jR
            reac_available(jR)=.FALSE.
            found_iR_equiv=.TRUE.
          END IF
          DEALLOCATE(j_educts)
          DEALLOCATE(j_eductsCoeff)
        END IF
      END DO ! jR
      IF ( .NOT. found_iR_equiv ) THEN
        ! if this is false, no equivalent for iR was found,
        ! so iSpc and jSpc cannot be lumped
        equivs_present=.FALSE.
        EXIT 
      END IF
      DEALLOCATE(i_educts)
      DEALLOCATE(i_eductsCoeff)
    END DO ! iR

  END SUBROUTINE check_reac_equivs



! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!                          TYPE-PROCEDURES
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  SUBROUTINE NextLumpingGroup(group, Spc, SpcReacs)
    INTEGER :: Spc
    INTEGER, DIMENSION(:) :: SpcReacs
    TYPE(lumping_group), POINTER :: group

    TYPE(lumping_group), POINTER :: temp => NULL()
    LOGICAL :: found_spc, reached_end

    found_spc = .FALSE.
    reached_end = .FALSE.   
    ! put a new group
    CALL group%put(Spc,SpcReacs)

    ! then point to that group (important) and catch some errors
    temp => group
    DO WHILE ( (.NOT. found_spc) .AND. (.NOT. reached_end) ) 
      IF ( temp%spc%id == Spc ) THEN
        group=>temp
        found_spc = .TRUE.
      ELSE IF ( ASSOCIATED(temp%next) ) THEN
        temp => temp%next
      ELSE
        reached_end = .TRUE. ! should not happen because Spc was put in group
      END IF
    END DO
    IF (.NOT. found_spc) THEN
      WRITE (*,*) 'ERROR: lumping group put did not work properly'
      STOP
    END IF
  END SUBROUTINE NextLumpingGroup


  RECURSIVE SUBROUTINE put_ll(list,id,given_id_char)
    CLASS(linked_list) , INTENT(INOUT) :: list
    INTEGER                            :: id
    CHARACTER(LenName) , OPTIONAL      :: given_id_char

    CHARACTER(LenName) :: id_char
    
    IF ( PRESENT(given_id_char) ) THEN
      id_char = given_id_char
    ELSE
      id_char = ''
    END IF

    IF ( .NOT. ASSOCIATED(list%next) .AND. list%id .EQ. -1 ) THEN
      list%id=id
      list%id_char=id_char
    ELSE
      IF ( .NOT. ASSOCIATED(list%next) ) ALLOCATE(list%next)
      CALL put_ll(list%next,id,id_char)
    END IF
  END SUBROUTINE put_ll 

  RECURSIVE SUBROUTINE free_ll(list)
    CLASS(linked_list), INTENT(inout) :: list
    IF (ASSOCIATED(list%next)) THEN
       CALL free_ll(list%next)
       DEALLOCATE(list%next)
    END IF
    list%next => NULL()
  END SUBROUTINE free_ll

  FUNCTION get_Spc_char(Spc, SpcReac) RESULT(Spc_char)
    INTEGER, INTENT(IN) :: Spc
    INTEGER, INTENT(IN), OPTIONAL :: SpcReac

    CHARACTER(LenName) :: Spc_char

    INTEGER :: some_reac,i

    Spc_Char=''

    IF ( PRESENT(SpcReac) .AND. SpcReac>0 .AND. ALLOCATED(ReactionSystem(SpcReac)%Educt) ) THEN
      DO i=LBOUND(ReactionSystem(SpcReac)%Educt,1),UBOUND(ReactionSystem(SpcReac)%Educt,1)
        IF ( ReactionSystem(SpcReac)%Educt(i)%iSpecies == Spc ) THEN
          Spc_char = ReactionSystem(SpcReac)%Educt(i)%Species
        END IF
      END DO
    ELSE
      ! find where Spc appears as product
      some_reac = FINDLOC(B%ColInd,Spc,DIM=1)
      ! find which reaction this is
      IF ( some_reac > 0 ) THEN
        some_reac = MINLOC(B%RowPtr,MASK=B%RowPtr>some_reac,DIM=1) - 1
      END IF

      IF ( ALLOCATED(ReactionSystem(some_reac)%Product) .AND. some_reac > 0) THEN
        DO i=LBOUND(ReactionSystem(some_reac)%Product,1),UBOUND(ReactionSystem(some_reac)%Product,1)
          IF ( ReactionSystem(some_reac)%Product(i)%iSpecies == Spc ) THEN
            Spc_char = ReactionSystem(some_reac)%Product(i)%Species
          END IF
        END DO
      ELSE 
        ! find where Spc appears as educt
        some_reac = FINDLOC(A%ColInd,Spc,DIM=1)
        ! find which reaction this is
        IF ( some_reac > 0 ) THEN
          some_reac = MINLOC(A%RowPtr,MASK=A%RowPtr>some_reac,DIM=1) - 1
        END IF  

        IF ( ALLOCATED(ReactionSystem(some_reac)%Educt) .AND. some_reac > 0) THEN
          DO i=LBOUND(ReactionSystem(some_reac)%Educt,1),UBOUND(ReactionSystem(some_reac)%Educt,1)
            IF ( ReactionSystem(some_reac)%Educt(i)%iSpecies == Spc ) THEN
              Spc_char = ReactionSystem(some_reac)%Educt(i)%Species
            END IF
          END DO
        END IF
      END IF
    END IF

    IF ( Spc_Char == '' ) THEN
      WRITE(*,*) 'ERROR: could not find char name of species ',Spc
      STOP
    END IF

  END FUNCTION get_Spc_char
  
  RECURSIVE SUBROUTINE put_lumping_group(group,Spc,SpcReacs,given_id)
    CLASS(lumping_group), INTENT(INOUT) :: group
    INTEGER, DIMENSION(:), INTENT(IN)   :: SpcReacs
    INTEGER, INTENT(IN)                 :: Spc
    INTEGER, OPTIONAL                   :: given_id

    INTEGER :: iReac, id, i, SpcReac
    CHARACTER(LenName) :: Spc_char = ''

    IF (PRESENT(given_id)) THEN
      id=given_id
    ELSE
      id=1
    END IF

    ! if there is no next group and the group itself is empty (<=> id=-1) then fill it with values
    IF ( (.NOT. ASSOCIATED(group%next)) .AND. (group%id == -1) ) THEN
      ALLOCATE(group%reactions(SIZE(SpcReacs)))
      DO iReac=1,SIZE(SpcReacs)
       CALL group%reactions(iReac)%put(id=SpcReacs(iReac))
      END DO
      
      Spc_char = get_Spc_char(Spc,SpcReacs(1))
      
      CALL group%spc%put(Spc,Spc_char)
      group%id=id
    ELSE ! else go on searching for an empty group
      IF ( .NOT. ASSOCIATED(group%next) ) ALLOCATE(group%next)
      id=group%id+1
      CALL put_lumping_group(group%next,Spc,SpcReacs,id)
    END IF
  END SUBROUTINE put_lumping_group
 
  RECURSIVE SUBROUTINE free_lumping_group(group)
    CLASS(lumping_group), INTENT(INOUT) :: group
    INTEGER :: i

    IF (ASSOCIATED(group%next)) THEN
      CALL free_lumping_group(group%next)
      DEALLOCATE(group%next)
    END IF
    group%next => NULL()
    
    DO i=1,SIZE(group%reactions)
      CALL group%reactions(i)%free() 
    END DO
    DEALLOCATE(group%reactions)
    CALL group%spc%free()
  END SUBROUTINE free_lumping_group

  SUBROUTINE lumping_group_add(group, Spc, reac_equivs, SpcReacs)
    ! input:
    !
    !  lumping group that is to be extended by Spc and its reactions
    CLASS(lumping_group) :: group
    !
    !  array to assign equivalent reactions 
    !  (reac_equivs(i)=j <=> j-th reaction of Spc is equivalent to i-th reaction of first Spc of group)
    INTEGER, DIMENSION(:) :: reac_equivs
    !  
    !  Spc id and vector of reaction id's
    INTEGER :: Spc
    INTEGER, DIMENSION(:) :: SpcReacs
    !
    ! end input

    INTEGER :: i
    CHARACTER(LenName) :: Spc_char

    Spc_char = get_Spc_char(Spc,SpcReacs(1))

    CALL group%spc%put(Spc,Spc_char)

    DO i=1,SIZE(reac_equivs)
      CALL group%reactions(i)%put(SpcReacs(reac_equivs(i)))
    END DO

    group%nspc = group%nspc + 1
    
  END SUBROUTINE lumping_group_add


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!                               WRITING
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




  SUBROUTINE WriteLumpingGroups(first_group)

    ! input:
    !
    TYPE(lumping_group) :: first_group

    INTEGER :: nLumpSpc, nLumpReacs, group_count

    TYPE(linked_list) :: unlumped_Spc

    nLumpSpc = 0
    nLumpReacs = A%m
    group_count = 0

    WRITE (*,*) ''
    WRITE (*,*) ' _________________________ '
    WRITE (*,*) '|                         |'
    WRITE (*,*) '|         LUMPING         |'
    WRITE (*,*) '|_________________________|'
    WRITE (*,*) ''
    WRITE (*,*) ''

    CALL WriteLumpingGroups_rec(first_group,nLumpSpc,nLumpReacs,group_count,unlumped_Spc)

  END SUBROUTINE WriteLumpingGroups

  RECURSIVE SUBROUTINE WriteLumpingGroups_rec(group,nLumpSpc,nLumpReacs,group_count,unlumped_Spc)
   
    ! input:
    !
    TYPE(lumping_group), TARGET :: group
    !
    INTEGER :: nLumpSpc, nLumpReacs, group_count
    !
    TYPE(linked_list), TARGET :: unlumped_Spc
    
    TYPE(linked_list), POINTER :: temp_list
    INTEGER, ALLOCATABLE :: Spcs_array(:), Spcs_offset(:), reacs_array(:)
    CHARACTER, ALLOCATABLE :: Spcs_unlumped(:)
    
    INTEGER :: i, j, iR, len_Spcs, n_Spcs
    CHARACTER(LEN=6) :: i_str, iR_str, id_str

    ! count this group as it will be one lumped species
    nLumpSpc = nLumpSpc + 1

    CALL LL2Array(group%spc,Spcs_array)
    n_Spcs=SIZE(Spcs_array)
    
    IF ( n_Spcs > 1 ) THEN
      group_count = group_count + 1
      
      ! count how many reacs are lumped (n_Spcs species are lumped to one -> (n_Spcs-1)*n_SpcReacs reactions are omitted
      nLumpReacs = nLumpReacs - (n_Spcs - 1)*SIZE(group%reactions)

      id_str=''
      WRITE (id_str,'(I6)') group_count
      WRITE (*,*) ''
      WRITE (*,*) ''
      WRITE (*,*) id_str//'. lumping group '
      WRITE (*,*) ''
      WRITE (*,*) '      -number of species:',group%nspc
      WRITE (*,'(A)',ADVANCE='NO') '       -species: '
      temp_list=>group%spc
      WRITE (*,'(A)',ADVANCE='NO') TRIM(ADJUSTL(temp_list%id_char))
      DO WHILE ( ASSOCIATED(temp_list%next) )
        temp_list=>temp_list%next
        WRITE (*,'(A)',ADVANCE='NO') ', '//TRIM(ADJUSTL(temp_list%id_char))
      END DO
      WRITE (*,*) ''
      WRITE (*,*) '      -reactions:'
      DO i=1,SIZE(group%reactions)
        i_str=''
        ! convert i to string
        WRITE (i_str,'(I6)') i
        CALL LL2Array(group%reactions(i), reacs_array)

        WRITE (*,'(A)',ADVANCE='NO') '    '//i_str//'. eq class: '
        DO iR=1,SIZE(reacs_array)
          iR_str=''
          WRITE (iR_str,'(I6)') reacs_array(iR)
          WRITE (*,'(A)',ADVANCE='NO') TRIM(iR_str)
        END DO
        WRITE (*,*) '' ! reset last advance=no
        !CALL WriteLinkedList(group%reactions(i))
      END DO
    ELSE !n_Spcs=1
     ! add Spc to unlumped species
     CALL unlumped_Spc%put(group%spc%id,group%spc%id_char)
    END IF

    IF (ASSOCIATED(group%next)) THEN
      CALL WriteLumpingGroups_rec(group%next, nLumpSpc, nLumpReacs, group_count, unlumped_Spc)
    ELSE
      ! END
      WRITE (*,*) ''
      WRITE (*,*) ''
      WRITE (*,*) ''
      WRITE (*,*) '     Unlumped Species:'
      temp_list=>unlumped_Spc
      WRITE (*,'(A)',ADVANCE='NO') '        '//TRIM(ADJUSTL(temp_list%id_char))
      DO WHILE ( ASSOCIATED(temp_list%next) )
        temp_list=>temp_list%next
        WRITE (*,'(A)',ADVANCE='NO') ', '//TRIM(ADJUSTL(temp_list%id_char))
      END DO
      WRITE (*,*) ''
      WRITE (*,*) ''
      WRITE (*,*) ''
      WRITE (*,*) ''
      WRITE (*,*) ' ______________________________________________________________'
      WRITE (*,*) '|                                                              |'
      WRITE (*,*) '|  Lumped',A%n,     'species and',A%m,       'reactions        |'
      WRITE (*,*) '|  to    ',nLumpSpc,'species and',nLumpReacs,'reactions.       |'
      WRITE (*,*) '|______________________________________________________________|'
    END IF
  END SUBROUTINE WriteLumpingGroups_rec

  !SUBROUTINE GasSpc_LL2String(Spcs, Spc_LL)
  !  CHARACTER, ALLOCATABLE :: Spcs(:)
  !  TYPE(linked_list) :: Spc_LL
  ! 
  !  INTEGER :: len_Spcs, i, j, n_Spcs
  !  INTEGER, ALLOCATABLE :: Spcs_array(:), Spcs_offset(:), reacs_array(:)
  ! 
  !  len_Spcs = 0
  !  CALL LL2Array(Spc_LL,Spcs_array)
  !  n_Spcs=SIZE(Spcs_array)
  !  ALLOCATE(Spcs_offset(n_Spcs+1))
  !  DO i=1,n_Spcs
  !    Spcs_offset(i)=len_Spcs+1
  !    ! add length of species name plus comma plus space
  !    len_Spcs=len_Spcs+LEN(TRIM(ListGas2(Spcs_array(i))%Species))+2
  !  END DO
  !  ! remove last comma and space
  !  !len_Spcs=len_Spcs-2
  !  Spcs_offset(n_Spcs+1)=len_Spcs+1
  !  ALLOCATE(Spcs(len_Spcs))
 ! 
 !   DO i=1,n_Spcs
 !     DO j=0,Spcs_offset(i+1)-3-Spcs_offset(i)
 !       Spcs(Spcs_offset(i)+j) = ListGas2(Spcs_array(i))%Species(j+1:j+1) ! j+1 to overcome initial blank space
 !     END DO
 !     Spcs(Spcs_offset(i+1)-2:Spcs_offset(i+1)-2) = ','
 !     Spcs(Spcs_offset(i+1)-1:Spcs_offset(i+1)-1) = ' '
 !   END DO
 !   ! delete last comma
 !   Spcs(len_Spcs-1)=' '
!
!  END SUBROUTINE GasSpc_LL2String

  RECURSIVE SUBROUTINE WriteLinkedList(list)
    TYPE(linked_list) :: list

    INTEGER :: i

    WRITE (*,*) list%id
    IF (ASSOCIATED(list%next)) THEN
      CALL WriteLinkedList(list%next)
    END IF
  END SUBROUTINE WriteLinkedList



! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!                                UTIL
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  SUBROUTINE LL2Array(list, list_array)
    TYPE(linked_list), TARGET, INTENT(IN) :: list
    INTEGER, ALLOCATABLE, INTENT(OUT) :: list_array(:)

    TYPE(linked_list), POINTER :: temp
    INTEGER :: i,len

    len=1
    temp=>list
    DO WHILE (ASSOCIATED(temp%next))
      len=len+1
      temp=>temp%next
    END DO
    ALLOCATE(list_array(len))
    temp=>list

    DO i=1,len
      list_array(i)=temp%id
      temp=>temp%next
    END DO
  END SUBROUTINE LL2Array

END MODULE Lumping_Mod
