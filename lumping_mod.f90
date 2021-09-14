MODULE Lumping_Mod
  
  USE Kind_Mod
  USE Sparse_Mod, ONLY: CSR_Matrix_T
  USE Control_Mod
  USE qsort_c_module
  USE Chemsys_Mod
  USE Reac_Mod, ONLY: y_name
  USE IO_Mod

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
     REAL(dp)                     :: weight = -1.0_dp
   CONTAINS
     PROCEDURE :: put  => put_ll
     PROCEDURE :: free => free_ll
  END TYPE linked_list

  TYPE linked_tree
     TYPE(linked_tree)  , ALLOCATABLE :: child(:)
     INTEGER                          :: id = -1
     CHARACTER(LenName)               :: id_char = ''
     INTEGER            , ALLOCATABLE :: weight(:)
     LOGICAL            , ALLOCATABLE :: relevant(:)
   CONTAINS
     PROCEDURE :: put  => put_tree
     PROCEDURE :: free => free_tree
     PROCEDURE :: cut  => cut_tree
  END TYPE linked_tree

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
     CHARACTER(LenLine)         :: Line3Template = ''
     CHARACTER(LenLine)         :: Line3 = ''
     CHARACTER(LenLine)         :: Line4 = ''
     TYPE(Duct_T), ALLOCATABLE  :: Educt(:)
     TYPE(Duct_T), ALLOCATABLE  :: Product(:)
     CHARACTER(LenType)         :: Type,  TypeConstant
     REAL(dp),     ALLOCATABLE  :: Constants(:)            
  END TYPE NewReac_Data

  TYPE Emis_Family
     CHARACTER(LenType)               :: Name = ''
     REAL(dp)                         :: Emis_val = mONE
     REAL(dp),           ALLOCATABLE  :: Portions(:)
     INTEGER,            ALLOCATABLE  :: Members(:)
     CHARACTER(LenName), ALLOCATABLE  :: Members_char(:)
  END TYPE Emis_Family

  ! positions of record in files
  INTEGER, ALLOCATABLE :: FluxPositions(:), ConcPositions(:), Positions(:)
  ! status of files
  LOGICAL :: fluxfile_exists, concfile_exists
  ! dummy counter
  INTEGER :: counter1 = 0, counter2 = 0

  CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!                              ALGORITHM
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

  SUBROUTINE lump_System(tau, ConcMatrix)
    ! input:
    !
    !  lifetimes of species at given times (tau(i,j) : lifetime of species i at timepoint j)
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: tau
    !  concentrations of species at given times (ConcMatrix(i,j) : concentration of species i at timepoint j)
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: ConcMatrix
    !
    ! end input
    

    ! lumped ReactionSystem
    TYPE(NewReac_Data), ALLOCATABLE :: l_ReactionSystem(:)

    ! lumping groups storage
    TYPE(lumping_group), TARGET  :: first_group 
    TYPE(lumping_group), POINTER :: current_group

    ! array of emission families (needed for calculating lumped stoichiometric coefficients)
    TYPE(Emis_Family), ALLOCATABLE :: EmisFams(:)
    ! tree of sources
    TYPE(linked_tree)              :: tree

    ! vector to test if species/reaction is already lumped
    LOGICAL, ALLOCATABLE :: Spc_lumped(:), reac_lumped(:)

    LOGICAL :: iSpc_lumpable, jSpc_lumpable, same_lifetimes, same_reactypes, equivs_present, similar_k, constant_sigma
    INTEGER :: timestep, nSpc, nReac, iSpc, jSpc, n_iSpcReacs, n_jSpcReacs, i, j, nLumpSpc, nLumpReacs, dummy
    INTEGER, ALLOCATABLE :: iSpcReacs(:), jSpcReacs(:), reac_equivs(:), Preserve_Spc(:), &
                          & iSpc_reac_types(:), jSpc_reac_types(:), Reacs_Array(:)
    REAL(dp), ALLOCATABLE :: time_flux(:), dt_flux(:)
    REAL(dp) :: Conc_ratio = 0.0_dp
    
    ! counters (overall and succesful for each lumping step)
    INTEGER :: cLumpablei  = 0      &
    &        , csLumpablei = 0      &
    &        , cTau        = 0      &
    &        , csTau       = 0      &
    &        , cLumpablej  = 0      &
    &        , csLumpablej = 0      &
    &        , cTypes      = 0      &
    &        , csTypes     = 0      &
    &        , cEquivs     = 0      &
    &        , csEquivs    = 0      &
    &        , cK          = 0      &
    &        , csK         = 0      &
    &        , cSigma      = 0      &
    &        , csSigma     = 0      &
    &        , cdummy      = 0

   !PERMTEST
    INTEGER, ALLOCATABLE :: Perm(:), SortVec(:)
    REAL(dp), ALLOCATABLE :: Rates(:), ratios(:)
    INTEGER        :: io_stat = 0, somespc = 0
    CHARACTER(200) :: io_msg  = ''
    CHARACTER(LenLine) :: TESTEST  = ''
    REAL(dp) :: sum1 = 0,sum2 = 0,sum3 = 0,sum4 = 0,sum5 = 0  
    LOGICAL :: testbool 

    CALL Logo3 



    ! get positions of record out of meta files
    FluxFile     = 'REDUCTION/flux_'//TRIM(BSP)//'.dat'
    FluxMetaFile = 'REDUCTION/fluxmeta_'//TRIM(BSP)//'.dat'   
    ConcFile     = 'LUMPING/conc_'//TRIM(BSP)//'.dat'
    ConcMetaFile = 'LUMPING/concmeta_'//TRIM(BSP)//'.dat'   

    CALL CollectPositions(FluxPositions, fluxfile_exists, iStpFlux, TimeFluxRead, FluxMetaUnit, FluxMetaFile, 2)
    CALL CollectPositions(ConcPositions, concfile_exists, iStpConc, TimeConcRead, ConcMetaUnit, ConcMetaFile, 0)

    !LumpingControlFile    = 'LUMPING/'//TRIM(BSP)//'.ctrl'
    current_group=>first_group
    nReac = A%m
    nSpc  = A%n
    cdummy = 0
    
    ALLOCATE(reac_lumped(nReac),                &
           & Spc_lumped(nSpc))                   

    Spc_lumped  = .FALSE.
    Reac_lumped = .FALSE.

    ! preserve species that are to be preserved
    Preserve_Spc = Read_Preserve_Spc(LumpingControlFile)
    DO i=1,SIZE(Preserve_Spc)
      iSpc = Preserve_Spc(i)
      Spc_lumped(iSpc) = .TRUE.
      CALL NextLumpingGroup(current_group,iSpc,iSpcReacs,n_iSpcReacs)
    END DO

    ! scan system and collect lumpable species in lumping groups
    DO iSpc=1,nSpc
      IF ( .NOT. Spc_lumped(iSpc) ) THEN
        CALL CounterNextStep(cdummy,cLumpablei)

        Spc_lumped(iSpc) = .TRUE.

        ! create new lumping group with iSpc as initial species and find iSpcReacs
        ! every species will be in one lumping group (possibly containing only this species)
        CALL NextLumpingGroup(current_group,iSpc,iSpcReacs,n_iSpcReacs)

        CALL check_lumpable(iSpc_lumpable,iSpc,iSpcReacs,n_iSpcReacs)

        IF ( iSpc_lumpable ) THEN
          CALL CounterNextStep(csLumpablei,cdummy)

          ! declare all reactions involving iSpc as educt as lumped
          ! this will be reversed if no species is lumped with iSpc
          DO i=1,n_iSpcReacs
            reac_lumped(iSpcReacs(i))=.TRUE.
          END DO

          ! find reaction types of species i (occuring numbers of educts of iSpcReacs)
          CALL find_reaction_types(iSpc_reac_types,iSpcReacs)

          ! check all species with greater indices than iSpc, 
          ! smaller ones have already been checked before
          DO jSpc=iSpc+1,nSpc
            CALL CounterNextStep(cdummy,cTau)

            CALL compare_lifetimes(same_lifetimes,tau,iSpc,jSpc)

            IF ( same_lifetimes ) THEN
              CALL CounterNextStep(csTau,cLumpablej)
              
              ! find reactions in which iSpc occurs as educt
              CALL Spc_involving_reacs(jSpcReacs,n_jSpcReacs,jSpc)
              
              CALL check_lumpable(jSpc_lumpable,jSpc,jSpcReacs,n_jSpcReacs)
              
              IF ( jSpc_lumpable ) THEN
                CALL CounterNextStep(csLumpablej,cTypes)
                
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
                  CALL CounterNextStep(csTypes, cEquivs)

                  ALLOCATE(reac_equivs(n_jSpcReacs))
                  ! now check for reaction equivalents (reactions with same reactants excluding iSpc and jSpc)
                  CALL check_reac_equivs(equivs_present,reac_equivs,iSpcReacs,jSpcReacs,iSpc,jSpc)

                  IF (equivs_present) THEN
                    CALL CounterNextStep(csEquivs,cK)
                    
                    CALL compare_k(similar_k,reac_equivs,iSpcReacs,jSpcReacs)
                    
                    IF (similar_k) THEN
                      CALL CounterNextStep(csK, cSigma)

                      CALL find_sigma(constant_sigma, Conc_ratio, ConcMatrix, iSpc, jSpc)

                      IF (constant_sigma) THEN
                        CALL CounterNextStep(csSigma,cdummy)
                      
                        ! LUMP jSpc to iSpc
                        CALL current_group%add(jSpc,reac_equivs,jSpcReacs,Conc_ratio)
                      
                        ! declare all jSpcReacs as lumped
                        DO j=1,n_jSpcReacs
                          reac_lumped(jSpcReacs(j))=.TRUE.
                        END DO

                        Spc_lumped(jSpc) = .TRUE.
                      
                      END IF
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
          ! otherwise correct the concentration ratios to sigmas
          IF ( current_group%nspc == 1 ) THEN
            DO i=1,n_iSpcReacs
              reac_lumped(iSpcReacs(i))=.FALSE.
            END DO
          ELSE
            CALL Ratios2Sigmas(current_group)
          END IF
          DEALLOCATE(iSpc_reac_types)
        END IF ! iSpc_lumpable
        DEALLOCATE(iSpcReacs)
      END IF ! iSpc not lumped
    END DO ! iSpc loop

    ! count how many species/reactions remain (n_Spcs species are lumped to one -> (n_Spcs-1)*n_SpcReacs reactions are omitted
    nLumpSpc  = 0
    nLumpReacs = A%m
    current_group=>first_group
    DO WHILE (ASSOCIATED(current_group))
      nLumpSpc = nLumpSpc + 1
      nLumpReacs = nLumpReacs - (current_group%nspc - 1)*SIZE(current_group%reactions)
      current_group=>current_group%next
    END DO 

    EmisFams = Read_Emis_Families(LumpingControlFile)

    CALL build_LumpedReacSys(l_ReactionSystem, first_group, reac_lumped)
         
    CALL Print_LumpedSysFile(l_ReactionSystem, first_group, 'LUMPING/'//TRIM(BSP)//'_lumped.sys')

    !CALL WriteLumpingGroups(first_group)
    CALL WriteCounters()
    WRITE(*,*) '          Number of ducts with coefficient < nano  = ', counter2

!    CALL WriteTree(tree)
!!somespc = PositionSpeciesAll('OCC(C)OO')
!somespc = PositionSpeciesAll('OCCC(O)(C)C(=O)OON(=O)=O')
!IF (somespc == 0) THEN 
!  WRITE(*,*) 'lumping spc not existing'
!ELSE
!WRITE(*,*) 'BUILD CALL INI'
!CALL BuildEmisTree(somespc,tree, (/ 0 /), (/ 0 /))
!!CALL WriteTree(tree)
!!CALL CleanEmisTree(tree)
!WRITE(*,*) '2MATLAB CALL INI'
!CALL EmisTreeToMatlab(tree)
!somespc = PositionSpeciesAll('OCCC(O)C(=O)OON(=O)=O')
!IF (somespc == 0) THEN 
!  WRITE(*,*) 'lumping spc not existent'
!  STOP
!END IF
!CALL tree%free()
!WRITE(*,*) 'BUILD CALL INI 2'
!CALL BuildEmisTree(somespc,tree, (/ 0 /), (/ 0 /))
!!CALL WriteTree(tree)
!!CALL CleanEmisTree(tree)
!WRITE(*,*) '2MATLAB CALL INI 2'
!CALL EmisTreeToMatlab(tree)
!    ! CALL GroupTrees2Matlab(315,first_group)
!END IF
    CALL end_lumping

    
    CALL CONTROL_SUBR(l_ReactionSystem)

   CONTAINS

    SUBROUTINE WriteCounters()

      WRITE(*,*) ''
      WRITE(*,*) '          Trials Analysis for each step of lumping:'
      WRITE(*,*) ''
      WRITE(*,*) ''
      WRITE(*,626) '              Step 1: new lumping species           - accepted: ', csLumpablei, ' of ', cLumpablei, '   (',100.0*REAL(csLumpablei)/REAL(cLumpablei), '%)'
      WRITE(*,626) '              Step 2: life time comparison          - accepted: ', csTau,       ' of ', cTau, '   (',100*REAL(csTau)/REAL(cTau), '%)'
      WRITE(*,626) '              Step 3: possible member of new group  - accepted: ', csLumpablej, ' of ', cLumpablej, '   (',100*REAL(csLumpablej)/REAL(cLumpablej), '%)'
      WRITE(*,626) '              Step 4: reaction type comparison      - accepted: ', csTypes,     ' of ', cTypes, '   (',100*REAL(csTypes)/REAL(cTypes), '%)'
      WRITE(*,626) '              Step 5: equivalent reaction pairs     - accepted: ', csEquivs,    ' of ', cEquivs, '   (',100*REAL(csEquivs)/REAL(cEquivs), '%)'
      WRITE(*,626) '              Step 6: similar rate constants        - accepted: ', csK,         ' of ', cK, '   (',100*REAL(csK)/REAL(cK), '%)'
      WRITE(*,626) '              Step 7: constant sigmas               - accepted: ', csSigma,     ' of ', cSigma, '   (',100*REAL(csSigma)/REAL(cSigma), '%)'
      WRITE(*,*) ''
      WRITE(*,*) ''

      626 FORMAT(A64,I8,A4,I8,A4,F5.1,A2)
    END SUBROUTINE WriteCounters

    SUBROUTINE CounterNextStep(succesful,new)
      INTEGER :: succesful, new

      succesful = succesful + 1
      new       = new + 1
    END SUBROUTINE CounterNextStep

    SUBROUTINE check_lumpable(Spc_lumpable, Spc, SpcReacs, n_SpcReacs)
      LOGICAL :: Spc_lumpable
      INTEGER :: Spc, n_SpcReacs
      INTEGER, ALLOCATABLE :: SpcReacs(:)
     
      INTEGER :: i

      Spc_lumpable = .TRUE.

      IF ( n_SpcReacs == 0 ) THEN
        Spc_lumpable = .FALSE. ! don't lump what exists only as product (or do??)
      ELSE
        DO i=1,n_SpcReacs
          IF ( reac_lumped(SpcReacs(i)) ) THEN
            Spc_lumpable=.FALSE.
          END IF
        END DO
      END IF
    END SUBROUTINE check_lumpable
  
    SUBROUTINE end_lumping()
  
      WRITE (*,*) ''
      WRITE (*,*) ''
      WRITE (*,*) '           ______________________________________________________________'
      WRITE (*,*) '          |                                                              |'
      WRITE (*,*) '          |  Lumped',A%n,     'species and',A%m,       'reactions        |'
      WRITE (*,*) '          |  to    ',nLumpSpc,'species and',nLumpReacs,'reactions.       |'
      WRITE (*,*) '          |______________________________________________________________|'
      WRITE (*,*) ''
  
    END SUBROUTINE end_lumping

  END SUBROUTINE lump_System

  SUBROUTINE CONTROL_SUBR(l_ReactionSystem)
    TYPE(NewReac_Data), DIMENSION(:) :: l_ReactionSystem(:)

    INTEGER :: i, j, k
    CHARACTER(LenName) :: char1, char2

    DO i=1,SIZE(y_name)
      IF ( PositionSpeciesAll(y_name(i))/=i ) THEN
        WRITE(*,*) 'ERRORI',i,y_name(i),PositionSpeciesAll(y_name(i))
      END IF
      DO j=1,SIZE(ReactionSystem)
        DO k=LBOUND(ReactionSystem(j)%Educt,1),UBOUND(ReactionSystem(j)%Educt,1)
          IF ( (ReactionSystem(j)%Educt(k)%iSpecies == i .AND.                                       &
             & TRIM(ADJUSTL(ReactionSystem(j)%Educt(k)%Species)) /= TRIM(ADJUSTL(y_name(i)))) .OR.   &
             & (ReactionSystem(j)%Educt(k)%iSpecies /= i .AND.                                       &
             & TRIM(ADJUSTL(ReactionSystem(j)%Educt(k)%Species)) == TRIM(ADJUSTL(y_name(i))))        ) THEN
            WRITE(*,*) ''
            WRITE(*,*) 'Reaction Educt analysis:'
            WRITE(*,*) 'Reaction ',j 
            WRITE(*,*) 'Species: i=',i
            WRITE(*,*) 'y_name=', y_name(i)
            WRITE(*,*) 'Species=',ReactionSystem(j)%Educt(k)%Species
            WRITE(*,*) 'iSpecies=', ReactionSystem(j)%Educt(k)%iSpecies
            WRITE(*,*) 'PosSpecAll(Species)=', PositionSpeciesAll(ReactionSystem(j)%Educt(k)%Species)
          END IF
        END DO
        DO k=LBOUND(ReactionSystem(j)%Product,1),UBOUND(ReactionSystem(j)%Product,1)
          IF ( (ReactionSystem(j)%Product(k)%iSpecies == i .AND.                                       &
             & TRIM(ADJUSTL(ReactionSystem(j)%Product(k)%Species)) /= TRIM(ADJUSTL(y_name(i)))) .OR.   &
             & (ReactionSystem(j)%Product(k)%iSpecies /= i .AND.                                       &
             & TRIM(ADJUSTL(ReactionSystem(j)%Product(k)%Species)) == TRIM(ADJUSTL(y_name(i))))        ) THEN
            WRITE(*,*) ''
            WRITE(*,*) 'Reaction Product analysis:'
            WRITE(*,*) 'Reaction ',j 
            WRITE(*,*) 'Species: i=',i
            WRITE(*,*) 'y_name=', y_name(i)
            WRITE(*,*) 'Species=',ReactionSystem(j)%Product(k)%Species
            WRITE(*,*) 'iSpecies=', ReactionSystem(j)%Product(k)%iSpecies
            WRITE(*,*) 'PosSpecAll(Species)=', PositionSpeciesAll(ReactionSystem(j)%Product(k)%Species)
          END IF
        END DO
      END DO
    END DO


    DO i=1,SIZE(l_ReactionSystem)
      DO j=LBOUND(l_ReactionSystem(i)%Educt,1),UBOUND(l_ReactionSystem(i)%Educt,1)
        IF ( l_ReactionSystem(i)%Educt(j)%iSpecies<1 .OR. TRIM(ADJUSTL(l_ReactionSystem(i)%Educt(j)%Species))=='' ) THEN
          WRITE(*,*) ''
          WRITE(*,*) 'lumped Reaction ',i,':'
          WRITE(*,*) 'Educt ',j,':',l_ReactionSystem(i)%Educt(j)%iSpecies,l_ReactionSystem(i)%Educt(j)%Species
        END IF
      END DO
      DO j=LBOUND(l_ReactionSystem(i)%Product,1),UBOUND(l_ReactionSystem(i)%Product,1)
        IF ( l_ReactionSystem(i)%Product(j)%iSpecies<1 .OR. TRIM(ADJUSTL(l_ReactionSystem(i)%Product(j)%Species))=='' ) THEN
          WRITE(*,*) ''
          WRITE(*,*) 'lumped Reaction ',i,':'
          WRITE(*,*) 'Product ',j,':',l_ReactionSystem(i)%Product(j)%iSpecies,l_ReactionSystem(i)%Product(j)%Species
        END IF
      END DO
    END DO

  END SUBROUTINE

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!                       ALGORITHM SUBROUTINES
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE Ratios2Sigmas(group)
    TYPE(lumping_group), TARGET :: group

    REAL(dp) :: sum_ratios
    TYPE(linked_list), POINTER :: temp

    sum_ratios = 0.0_dp
    temp=>group%spc
    DO
      sum_ratios = sum_ratios + temp%weight
      IF (.NOT. ASSOCIATED(temp%next)) EXIT
      temp=>temp%next
    END DO

    temp=>group%spc
    DO
      temp%weight = temp%weight / sum_ratios
      IF (.NOT. ASSOCIATED(temp%next)) EXIT
      temp=>temp%next
    END DO
  END SUBROUTINE Ratios2Sigmas

  SUBROUTINE find_sigma(constant_sigma, Conc_ratio, ConcMatrix, iSpc, jSpc)
    LOGICAL :: constant_sigma
    REAL(dp) :: Conc_ratio
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: ConcMatrix
    INTEGER :: iSpc, jSpc

    REAL(dp) :: expected, varcoeff

    ! calculate expected value of the ratio conc(jSpc)/conc(iSpc)
    ! care with +eps!!!!!!!!!!! 1e-16 could be much in mol!
    expected = SUM( ConcMatrix(jSpc,2:)/(ConcMatrix(iSpc,2:)+eps) )/(SIZE(ConcMatrix,2)-1)

    ! calculate empirical variation coefficient
    ! sum of squared differences
    varcoeff = SUM(( ConcMatrix(jSpc,2:)/(ConcMatrix(iSpc,2:)+eps)-expected )**2)
    ! unbiased standard dev
    varcoeff = varcoeff / (SIZE(ConcMatrix,2)-2)
    ! variation coefficient (~ procentual deviation)
    varcoeff = SQRT(varcoeff)/expected

    IF ( varcoeff<.5) THEN
      constant_sigma = .TRUE.
    ELSE
      constant_sigma = .FALSE.
    END IF
    Conc_ratio = expected

  END SUBROUTINE find_sigma

  SUBROUTINE GroupTrees2Matlab(group_id, first_group, clean_given)
    INTEGER :: group_id
    TYPE(lumping_group), TARGET :: first_group
    LOGICAL, OPTIONAL :: clean_given

    TYPE(lumping_group), POINTER :: temp
    TYPE(linked_list), POINTER :: spc
    INTEGER :: i
    TYPE(linked_tree) :: tree
    LOGICAL :: clean

    IF ( PRESENT(clean_given) ) THEN
      clean = clean_given
    ELSE
      clean = .FALSE.
    END IF
  
    temp=>first_group
    i=0
    ! find group_id-th lumping group (group = more than one species)
    DO 
      IF (temp%nspc>1) i=i+1
      IF (i==group_id) EXIT
      temp=>temp%next
    END DO

    spc=>temp%spc
    DO i=1,temp%nspc
WRITE(*,*) '1'
      CALL tree%free()
WRITE(*,*) '2'
      !CALL BuildEmisTree(spc%id,tree,0)
WRITE(*,*) '3'
      IF (clean) CALL CleanEmisTree(tree)
WRITE(*,*) '4'
      CALL EmisTreeToMatlab(tree)

WRITE(*,*) '5'
      spc => spc%next
WRITE(*,*) '6'
    END DO

  END SUBROUTINE GroupTrees2Matlab

  FUNCTION CheckIfLumped(Spc1,Spc2,first_group) RESULT(LumpedTogether)
    INTEGER, INTENT(IN) :: Spc1, Spc2
    TYPE(lumping_group), TARGET, INTENT(IN) :: first_group

    LOGICAL :: LumpedTogether

    TYPE(lumping_group), POINTER :: group
    INTEGER, ALLOCATABLE :: groupspc(:)

    group=>first_group
    LumpedTogether = .FALSE.
    DO
      CALL LL2Array(group%spc, groupspc)
      IF ( ANY( groupspc==Spc1 ) .AND. ANY( groupspc==Spc2 ) ) THEN 
        LumpedTogether = .TRUE.
        EXIT
      END IF

      IF (ASSOCIATED(group%next)) THEN
        group=>group%next
      ELSE
        EXIT
      END IF
    END DO
  END FUNCTION CheckIfLumped

  RECURSIVE SUBROUTINE CheckTreeAllocs(tree, level_given)
    TYPE(linked_tree), TARGET :: tree
    INTEGER, OPTIONAL :: level_given
    
    INTEGER :: i, level

    IF ( PRESENT(level_given) ) THEN
      level = level_given
    ELSE
      level = 1
    END IF
    
    IF (ALLOCATED(tree%child) .AND. (.NOT. ALLOCATED(tree%weight))) THEN
      WRITE(*,*) 'Child but not weight allocated at level', level, ' for id ', tree%id
    ELSE IF (ALLOCATED(tree%weight) .AND. (.NOT. ALLOCATED(tree%child))) THEN
      WRITE(*,*) 'Weight but not child allocated at level ', level, ' for id ', tree%id
    END IF

    IF ( ALLOCATED(tree%child) ) THEN
IF (LBOUND(tree%child,1)/=1) WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ERROR WRONG LBOUND'
      DO i=1,SIZE(tree%child)
        CALL CheckTreeAllocs(tree%child(i), level+1)
      END DO
    END IF
  END SUBROUTINE CheckTreeAllocs

  SUBROUTINE cut_tree(tree,pos)
    CLASS(linked_tree), INTENT(INOUT) :: tree
    TYPE(linked_tree), ALLOCATABLE :: new_child(:)
    INTEGER, ALLOCATABLE :: new_weight(:)
    INTEGER :: pos, nChild, nWeight, i

WRITE(*,*) 'in cut check 1'
CALL CheckTreeAllocs(tree)

    nChild=SIZE(tree%child)
    nWeight=SIZE(tree%weight)
    IF ( nChild/=nWeight ) THEN
      WRITE(*,*) 'WARNING, cut_tree: size of child does not match size of weight in ',tree%id, '. No cutting done.'
      RETURN
    END IF
    IF (pos>nChild .OR. pos<1) THEN
      WRITE(*,*) 'WARNING, cut_tree: pos does not match child size of tree of ', tree%id,'. No cutting done.'
      RETURN
    END IF

    ALLOCATE( new_child(nChild-1), new_weight(nChild-1) )
    DO i=1,nChild
      IF ( i/=pos ) THEN
        IF (i>pos) THEN
          new_child(i-1)  = tree%child(i)
          new_weight(i-1) = tree%weight(i)
        ELSE
          new_child(i)  = tree%child(i)
          new_weight(i) = tree%weight(i)
        END IF
      END IF
    END DO
    !IF (pos==1) THEN
    !  new_child  = tree%child(pos+1:nChild)
    !  new_weight = tree%weight(pos+1:nChild)
    !ELSE IF (pos==nChild) THEN
    !  new_child  = tree%child(1:pos-1)
    !  new_weight = tree%weight(1:pos-1)
    !ELSE
    !  new_child(1:pos-1)       = tree%child(1:pos-1)
    !  new_child(pos:nChild-1)  = tree%child(pos+1:nChild)
    !  new_weight(1:pos-1)      = tree%weight(1:pos-1)
    !  new_weight(pos:nChild-1) = tree%weight(pos+1:nChild)
    !END IF

DO i=1,nChild-1
WRITE(*,*) 'in cut check news1',i
CALL CheckTreeAllocs(new_child(i))
END DO
WRITE(*,*) 'in cut check 2'
CALL CheckTreeAllocs(tree)

    DO i=1,nChild
      CALL tree%child(i)%free()
    END DO

DO i=1,nChild-1
WRITE(*,*) 'in cut check news2',i
CALL CheckTreeAllocs(new_child(i))
END DO
WRITE(*,*) 'in cut check 3'
CALL CheckTreeAllocs(tree)
    
    DEALLOCATE( tree%child, tree%weight )
    ALLOCATE( tree%child(nChild-1), tree%weight(nChild-1) )

WRITE(*,*) 'in cut check 4'
CALL CheckTreeAllocs(tree)

    tree%child  = new_child
    tree%weight = new_weight

DO i=1,nChild-1
WRITE(*,*) 'in cut check news3',i
CALL CheckTreeAllocs(new_child(i))
END DO
WRITE(*,*) 'in cut check 5'
CALL CheckTreeAllocs(tree)

WRITE(*,*) 'in cut check 6'
CALL CheckTreeAllocs(tree)
WRITE(*,*) '811', SIZE(tree%child), SIZE(tree%weight)
  END SUBROUTINE cut_tree
  
  RECURSIVE SUBROUTINE CleanEmisTree(tree)
    TYPE(linked_tree), INTENT(INOUT) :: tree

    REAL(dp), ALLOCATABLE :: ratios(:), dummy(:)
    REAL(dp) :: tol
    INTEGER, ALLOCATABLE :: Perm(:)
    INTEGER :: i, j, pos


    IF (( tree%id==-1                  ) .OR. &
       &( .NOT. ALLOCATED(tree%child)  ) .OR. &
       &( .NOT. ALLOCATED(tree%weight) )      ) RETURN

    IF ( .NOT. fluxfile_exists ) THEN 
      WRITE(*,*) 'Missing flux file in CleanEmisTree'
      RETURN
    END IF

    IF (ALLOCATED(tree%relevant)) DEALLOCATE(tree%relevant)
    ALLOCATE(tree%relevant(SIZE(tree%child)))
    tree%relevant = .TRUE.

    ! get ratios of mass moved by reactions producing tree%id
    ratios = MassMoveRatios(tree%weight)
    CALL SortVecAscReal(ratios,Perm)

!WRITE(*,*) ''
!WRITE(*,*) 'weight of ',y_name(tree%id)
!WRITE(*,*) tree%weight
!WRITE(*,*) 'ratios of ', y_name(tree%id)
!WRITE(*,*) ratios
!WRITE(*,*) ''

    ! tolerance to omit all other branches
    tol = .9

    ! find number of branches to keep
    DO i=1,SIZE(ratios)-1
      IF ( SUM(ratios(i+1:))>tol ) THEN
        !pos=Perm(i)
        ! some branches have been cut -> new positions in array, find new position
        !DO j=1,i-1
        !  IF ( Perm(j)<Perm(i) ) pos=pos-1
        !END DO
        !CALL tree%cut(pos)
        tree%relevant(Perm(i)) = .FALSE.
      ELSE
        EXIT
      END IF
    END DO

    DO i=1,SIZE(tree%child)
      !CALL CleanEmisTree(tree%child(i))
      IF ( tree%relevant(i) ) THEN
        CALL CleanEmisTree(tree%child(i))
      END IF
    END DO
  END SUBROUTINE CleanEmisTree

  FUNCTION MassMoveRatios(reacs) RESULT(ratios)
    INTEGER, DIMENSION(:) :: reacs
    
    REAL(dp), ALLOCATABLE :: ratios(:), Rates(:)
    INTEGER :: i, nReacs, io_stat = 0, dummy
    REAL(dp) :: total
    CHARACTER(200) :: io_msg  = ''

    IF ( .NOT. fluxfile_exists ) THEN
      WRITE(*,*) 'Missing flux file in MassMoveRatios'
      RETURN
    END IF

    nReacs = SIZE(reacs)

    IF (ALLOCATED(ratios)) DEALLOCATE(ratios)
    ALLOCATE(Rates(nr),ratios(nReacs))
    CALL OpenFile_rStream(FluxUnit,FluxFile)
    DO i=1,iStpFlux
      READ( FluxUnit, POS=Positions(i), IOSTAT=io_stat, IOMSG=io_msg) Rates
      IF ( io_stat>0 ) WRITE(*,'(10X,A,I0,A)') '   ERROR reading fluxes.dat :: ',io_stat,'  '//TRIM(io_msg)
      IF ( io_stat<0 ) WRITE(*,'(10X,A,I0,A)') '   WARNING (?) reading fluxes.dat :: ',io_stat,'  '//TRIM(io_msg)
      ratios = ratios + ABS(Rates(reacs))
    END DO
    CLOSE(FluxUnit)

    ! check for NaN or large number
    DO i=1,nReacs
      IF ( ratios(i)/=ratios(i) ) THEN
        WRITE(*,*) 'WARNING: NaN read from flux file.'
      ELSE IF (ratios(i) > giga**11) THEN
        WRITE(*,*) 'WARNING: Very large number (>1e99) read from flux file.'
      END IF
    END DO

    total = SUM(ratios)
    !WRITE(*,*) 'portions', ratios
    !WRITE(*,*) 'total', total
    ratios = ratios/total

  END FUNCTION MassMoveRatios

  RECURSIVE SUBROUTINE CollectNodesEdges(tree,starts,ends,weights,mask)
    TYPE(linked_tree) :: tree
    TYPE(linked_list), TARGET :: starts, ends, weights, mask

    INTEGER :: i,j, eductid
    CHARACTER(LenName) :: eductsname = '', eductname = ''
    TYPE(linked_list), POINTER :: tmpstarts, tmpends, tmpweights
    LOGICAL :: edge_exists

    IF ( ALLOCATED(tree%child) ) THEN
      IF ( .NOT. ALLOCATED(tree%relevant) ) THEN
        ALLOCATE(tree%relevant(SIZE(tree%child)))
        tree%relevant = .TRUE.
      END IF
      DO i=1,SIZE(tree%child)
        IF ( tree%relevant(i) ) THEN
          eductname  = ''
          eductsname = ''
          edge_exists = .FALSE.
          tmpstarts  => starts
          tmpends    => ends
          tmpweights => weights
 
          ! find out if edge already exists
          DO 
            IF ( tmpstarts%id  == tree%id          .AND. &
               & tmpends%id    == tree%child(i)%id .AND. &
               & tmpweights%id == tree%weight(i)         ) THEN
            
              edge_exists = .TRUE.
              EXIT
            END IF
            IF ( ASSOCIATED(tmpstarts%next)  .AND. &
                 ASSOCIATED(tmpends%next)    .AND. &
                 ASSOCIATED(tmpweights%next)       ) THEN
                 tmpstarts  => tmpstarts%next
                 tmpends    => tmpends%next
                 tmpweights => tmpweights%next
            ELSE
              EXIT
            END IF
          END DO

          IF ( .NOT. edge_exists ) THEN
            ! find weight (=other educts)
            IF ( SIZE(ReactionSystem(tree%weight(i))%Educt)>1 ) THEN
              DO j=LBOUND(ReactionSystem(tree%weight(i))%Educt,1),UBOUND(ReactionSystem(tree%weight(i))%Educt,1)
                eductid   = ReactionSystem(tree%weight(i))%Educt(j)%iSpecies
                eductname = ReactionSystem(tree%weight(i))%Educt(j)%Species
                IF (eductid/=tree%child(i)%id) eductsname = TRIM(ADJUSTL(eductsname))//TRIM(ADJUSTL(eductname))
              END DO
            ELSE
              eductsname = 'nothing'
            END IF

            ! collect start, end of edge with reaction and reacting educts as weight
            CALL starts%put(tree%id)
            CALL ends%put(tree%child(i)%id)
            CALL weights%put(tree%weight(i),eductsname)
          END IF
        END IF
      END DO

      DO i=1,SIZE(tree%child)
        IF ( tree%relevant(i) ) THEN
          IF ( SearchLinkedList(mask,id=tree%child(i)%id)==0 ) CALL mask%put(tree%child(i)%id)
          CALL CollectNodesEdges(tree%child(i),starts,ends,weights,mask)
        END IF
      END DO
    END IF
  END SUBROUTINE CollectNodesEdges

  RECURSIVE SUBROUTINE BuildEmisTree(Spc, tree, parentspc, parentreac)
    ! input:
    INTEGER, INTENT(IN) :: Spc
    !
    INTEGER, DIMENSION(:) :: parentreac, parentspc
    !
    ! output:
    TYPE(linked_tree) :: tree

    TYPE(linked_list) :: children
    INTEGER :: current, current_reac, i, j, iEduct
    CHARACTER(LenName) :: Educt

    LOGICAL :: testest


    WRITE(*,*) SIZE(parentspc)

    tree%id = Spc
    tree%id_char = y_name(Spc)

    ! find reactions that produce Spc
    current = 0
    DO WHILE ( current<SIZE(B%ColInd) )

      IF ( FINDLOC(B%ColInd(current+1:),Spc,DIM=1)/=ZERO ) THEN ! another reaction with Spc as product was found
        ! refresh current index
        current = current + FINDLOC(B%ColInd(current+1:),Spc,DIM=1)
        ! find current reaction
        current_reac = MINLOC(B%RowPtr,MASK=B%RowPtr>current,DIM=1) - 1
 
        ! add educts that are C-species as children of Spc (parent in chemical chronology)
        DO i=LBOUND(ReactionSystem(current_reac)%Educt,1),UBOUND(ReactionSystem(current_reac)%Educt,1)
              
          iEduct = ReactionSystem(current_reac)%Educt(i)%iSpecies
          Educt  = ReactionSystem(current_reac)%Educt(i)%Species

          ! check if it could be a child
          IF (      INDEX(Educt,'C')  > 0       .AND. &
             &  TRIM(ADJUSTL(Educt)) /='CO'     .AND. &
             &  TRIM(ADJUSTL(Educt)) /='CL'     .AND. &
             &  FINDLOC(parentspc,iEduct,DIM=1) == 0        ) THEN

            ! current_reac serves as a marker (reaction id) to find the corresponding rate constant later
            CALL tree%put(   ReactionSystem(current_reac)%Educt(i)%iSpecies &
                         & , current_reac                                   &
                         & , ReactionSystem(current_reac)%Educt(i)%Species  ) 
            CALL children%put(iEduct)
          ELSE IF (FINDLOC(parentspc,iEduct,DIM=1)>0 .AND. FINDLOC(parentspc,iEduct,DIM=1)<SIZE(parentspc)) THEN
          !ELSE IF (FINDLOC(parentspc,iEduct,DIM=1)>0) THEN
            testest=.FALSE.
            WRITE(*,*) 'iEduct=',TRIM(ADJUSTL(y_name(iEduct))),', Spc=',TRIM(ADJUSTL(y_name(Spc)))
            WRITE(*,*) 'PRE spc:'
            DO j=1,SIZE(parentspc)
              WRITE(*,*) TRIM(ADJUSTL(y_name(parentspc(j))))
            END DO
            WRITE(*,*) ''
            DO j=1,SIZE(parentspc)
              IF (testest) THEN
                WRITE(*,'(A4,I5,A2,A)',ADVANCE='NO') ' <=(',parentreac(j-1),') ', TRIM(ADJUSTL(y_name(parentspc(j))))
              END IF
              IF (parentspc(j) == iEduct) THEN
                testest = .TRUE.
                WRITE(*,*) TRIM(ADJUSTL(y_name(iEduct)))
              END IF
            END DO
            WRITE(*,'(A4,I5,A2,A)',ADVANCE='NO') ' <=(',parentreac(SIZE(parentreac)),') ', TRIM(ADJUSTL(y_name(Spc)))     
            WRITE(*,'(A4,I5,A2,A)',ADVANCE='NO') ' <=(',current_reac,') ', TRIM(ADJUSTL(y_name(iEduct)))     
            WRITE(*,*) ''
            WRITE(*,*) ''
            
            STOP

          END IF
        END DO
      ELSE ! Spc did not appear as product again
        EXIT
      END IF
    END DO
    CALL children%free()

    ! recursively find source emissions of the added children
    IF ( ALLOCATED(tree%child) ) THEN 
      DO i=1,SIZE(tree%child)
        CALL BuildEmisTree(tree%child(i)%id, tree%child(i), [parentspc, Spc], [parentreac, tree%weight(i)])
      END DO
    END IF
  END SUBROUTINE BuildEmisTree

!  RECURSIVE SUBROUTINE BuildEmisTree(Spc, tree, parent)
!    ! input:
!    INTEGER, INTENT(IN) :: Spc
!    !
!    INTEGER, INTENT(IN) :: parent
!    !
!    ! output:
!    TYPE(linked_tree) :: tree
!
!    TYPE(linked_list) :: children
!    INTEGER :: current, current_reac, i, iEduct
!    CHARACTER(LenName) :: Educt
!
!    LOGICAL :: testest
!    
!
!    tree%id = Spc
!    tree%id_char = y_name(Spc)
!
!    ! find reactions that produce Spc
!    current = 0
!    DO WHILE ( current<SIZE(B%ColInd) )
!
!      IF ( FINDLOC(B%ColInd(current+1:),Spc,DIM=1)/=ZERO ) THEN ! another reaction with Spc as product was found
!        ! refresh current index
!        current = current + FINDLOC(B%ColInd(current+1:),Spc,DIM=1)
!        ! find current reaction
!        current_reac = MINLOC(B%RowPtr,MASK=B%RowPtr>current,DIM=1) - 1
!
!        ! TEST IF REAC CONTAINS SPC AS PRODUCT
!        !testest=.FALSE.
!        !DO i=LBOUND(ReactionSystem(current_reac)%Product,1),UBOUND(ReactionSystem(current_reac)%Product,1)
!        !  IF (ReactionSystem(current_reac)%Product(i)%iSpecies==Spc) testest=.TRUE.
!        !END DO
!        !IF (.NOT. testest) THEN 
!        !  WRITE(*,*) 'ERROR WRONG FUNCTION FOUND'
!        !ELSE 
!        !  WRITE(*,*) 'RIGHT FUNCTION FOUND'
!        !END IF
! 
!        ! add educts that are C-species as children of Spc (parent in chemical chronology)
!        DO i=LBOUND(ReactionSystem(current_reac)%Educt,1),UBOUND(ReactionSystem(current_reac)%Educt,1)
!              
!          iEduct = ReactionSystem(current_reac)%Educt(i)%iSpecies
!          Educt  = ReactionSystem(current_reac)%Educt(i)%Species
!
!          ! check if it could be a child
!          IF (      INDEX(Educt,'C')  > 0       .AND. &
!             &  TRIM(ADJUSTL(Educt)) /='CO'     .AND. &
!             &  TRIM(ADJUSTL(Educt)) /='CL'     .AND. &
!             &              iEduct   /= parent        ) THEN
!
!!WRITE(*,*) TRIM(ADJUSTL(y_name(Spc))),' SOURCE: ',ReactionSystem(current_reac)%Educt(i)%iSpecies, TRIM(ADJUSTL(ReactionSystem(current_reac)%Educt(i)%Species)), ' IN REACTION ',current_reac       
!WRITE(*,*) '1'
!CALL WriteChildren(tree)
!
!            ! current_reac serves as a marker (reaction id) to find the corresponding rate constant
!            CALL tree%put(   ReactionSystem(current_reac)%Educt(i)%iSpecies &
!                         & , current_reac                                   &
!                         & , ReactionSystem(current_reac)%Educt(i)%Species  ) 
!            CALL children%put(iEduct)
!WRITE(*,*) '2'
!CALL WriteChildren(tree)
!          END IF
!        END DO
!      ELSE ! Spc did not appear as product again
!        EXIT
!      END IF
!    END DO
!    CALL children%free()
!
!    ! recursively find source emissions of the added children
!    IF ( ALLOCATED(tree%child) ) THEN 
!      DO i=1,SIZE(tree%child)
!        CALL BuildEmisTree(tree%child(i)%id, tree%child(i), Spc)
!      END DO
!    END IF
!  END SUBROUTINE BuildEmisTree

  SUBROUTINE WriteChildren(tree)
    TYPE(linked_tree) :: tree

    INTEGER, ALLOCATABLE :: children(:)
    INTEGER :: i

    IF (ALLOCATED(tree%child)) THEN
      ALLOCATE(children(SIZE(tree%child)))
      DO i=1,SIZE(tree%child)
        children(i) = tree%child(i)%id
      END DO

      WRITE(*,*) 'Children of ',tree%id,': ',children
    ELSE
      WRITE(*,*) 'Children of ',tree%id,': none'
    END IF
  END SUBROUTINE WriteChildren

  SUBROUTINE compare_lifetimes(same_lifetimes, tau, iSpc, jSpc)
    LOGICAL :: same_lifetimes
    REAL(dp), DIMENSION(:,:) :: tau
    INTEGER :: iSpc,jSpc

    INTEGER :: timestep,nTau

    nTau=SIZE(tau,2)
      
    same_lifetimes=.TRUE.
    DO timestep=1,nTau
    !~~~~~~~~~~~~~~~~~~~~~
    ! maybe vectorize this loop if faster
    !~~~~~~~~~~~~~~~~~~~~~
      IF (ABS((tau(iSpc,timestep)-tau(jSpc,timestep))/(tau(iSpc,timestep)+100*eps))>eps_tau) THEN
        same_lifetimes=.FALSE.
      END IF
    END DO  
  END SUBROUTINE compare_lifetimes


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
    INTEGER :: i, j, iR, kR, kE, NewReac, nLumpedReacs, nRemainingReacs, nNewReacs,          &
             & nLumpedSpc, nRemainingSpc, nNewSpc, nEducts, nProducts
    
    INTEGER , ALLOCATABLE :: iSpc_Old2New(:), Spc_array(:), Reacs_Array(:)
    REAL(dp), ALLOCATABLE :: sigma(:), group_sigma(:)

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

    ALLOCATE( Spc_Old2New(A%n+SIZE(ListNonReac2))   &
            , iSpc_Old2New(A%n+SIZE(ListNonReac2))  &
            , l_RS(nNewReacs)                       &
            , sigma(A%n+SIZE(ListNonReac2))         )
    
    ! create mask: Spc_Old2New(old_species_id) = new_species_id
    Spc_Old2New = ''
    sigma = ONE
    group=>first_group
    DO i=1,nNewSpc
      CALL LL2Array(group%spc,Spc_array)
      DO j=1,SIZE(Spc_array)
        iSpc_Old2New(Spc_array(j))=i
      END DO
      IF ( group%nspc > 1 ) THEN
        CALL NewProperties(group, Spc_Old2New, sigma)
      END IF
      group=>group%next
    END DO

    IF (MINVAL(sigma)<0 .OR. MAXVAL(sigma)>1) THEN 
      WRITE(*,*) 'Sigma contains values out of [0,1]. Abort!'
      STOP
    END IF

    ! add passive species
    DO i=1,SIZE(ListNonReac2)
      iSpc_Old2New(A%n+i) = nNewSpc + i
      !Spc_Old2New not changed, passive species wont be lumped
    END DO
    
    ! collect all lumped reactions
    NewReac = 1
    group=>first_group
    DO i=1,nNewSpc
      IF ( group%nspc > 1 ) THEN
        CALL LL2Array_weight(group%spc,group_sigma)
        DO j=1,SIZE(group%reactions)
          CALL LL2Array(group%reactions(j),Reacs_Array)
          CALL LumpReacs(Reacs_array,group_sigma)
        END DO
      END IF
      group=>group%next
    END DO
    ! collect all remaining reactions
    DO i=1,SIZE(reac_lumped)
      IF ( .NOT. reac_lumped(i) ) THEN
        CALL LumpReacs((/i/),(/ONE/))
      END IF
    END DO

   CONTAINS

     SUBROUTINE LumpReacs(Reacs_old, group_sigma)
       INTEGER, DIMENSION(:) :: Reacs_old
       REAL(dp), DIMENSION(:) :: group_sigma

       INTEGER :: pos, k, l, iR_old, currentColon, currentSpace, counter
       TYPE(Duct_T) :: current
       CHARACTER(LenLine) :: Line3, Line3Template
       TYPE(linked_list) :: products

       iR_old = Reacs_old(1)
       ALLOCATE(l_RS(NewReac)%Constants(SIZE(ReactionSystem(iR_old)%Constants)))
       
       ! calculate number of educts and products
       ! meanwhile collect the reaction constants and use the mean values for new constants
       nEducts = SIZE(ReactionSystem(iR_old)%Educt)
       nProducts = 0
       l_RS(NewReac)%Constants = ZERO
       DO k=1,SIZE(Reacs_old)
         kR = Reacs_old(k)
         l_RS(NewReac)%Constants = l_RS(NewReac)%Constants + ReactionSystem(kR)%Constants
         
         DO l=LBOUND(ReactionSystem(kR)%Product,1),UBOUND(ReactionSystem(kR)%Product,1)
           ! check if product already appeared
           pos = SearchLinkedList(products,id=iSpc_Old2New(ReactionSystem(kR)%Product(l)%iSpecies))
           IF (pos == 0) THEN
             CALL products%put(iSpc_Old2New(ReactionSystem(kR)%Product(l)%iSpecies))
             nProducts = nProducts + 1
           END IF
         END DO
       END DO
       l_RS(NewReac)%Constants = l_RS(NewReac)%Constants / SIZE(Reacs_old)
          
       ALLOCATE( l_RS(NewReac)%Educt(nEducts) ) 
       ALLOCATE( l_RS(NewReac)%Product(nProducts) )

       ! write new educts and products, 
       ! lumped species gets the name of the first species in lumping group
       ! educts:
       DO k=1,nEducts
         current = ReactionSystem(iR_old)%Educt(k)
         l_RS(NewReac)%Educt(k)          = current
         ! correct id to new species id
         l_RS(NewReac)%Educt(k)%iSpecies = iSpc_Old2New(current%iSpecies)
         ! note: no sigma correction needad as the lumped species gets sigma=1 as educt and other educts can't be lumped spc
       END DO
       ! products:
       CALL products%free()
       counter = 1
       DO k=1,SIZE(Reacs_old)
         kR=Reacs_old(k)
         DO l=LBOUND(ReactionSystem(kR)%Product,1),UBOUND(ReactionSystem(kR)%Product,1)

           current = ReactionSystem(kR)%Product(l)
           ! search for earlier appearance of current product species
           pos = SearchLinkedList(products,id=iSpc_Old2New(current%iSpecies))
           
           IF ( pos == 0 ) THEN ! new product species
        
             CALL products%put(iSpc_Old2New(current%iSpecies))

             l_RS(NewReac)%Product(counter) = current

             ! apply sigma for all products of one lumped species
             l_RS(NewReac)%Product(counter)%Koeff = group_sigma(k) * l_RS(NewReac)%Product(counter)%Koeff
             
             ! concern lumping status of product
             IF ( Spc_Old2New(current%iSpecies) /= '' ) THEN ! duct has a new name and so was lumped
               ! apply new name
               l_RS(NewReac)%Product(counter)%Species  = Spc_Old2New(current%iSpecies)
               ! correct Koeff to Koeff*sigma (not by this lumping step, but because current product was lumped itself)
               l_RS(NewReac)%Product(counter)%Koeff = sigma(current%iSpecies)*l_RS(NewReac)%Product(counter)%Koeff
             END IF
             ! correct id to new species id
             l_RS(NewReac)%Product(counter)%iSpecies = iSpc_Old2New(current%iSpecies)
             counter=counter+1

           ELSE ! product species already appeared as product, increase stoichometric coeff
        
             ! correct Koeff to Koeff*sigma (sigma=1 if spc was not lumped)
             l_RS(NewReac)%Product(pos)%Koeff = l_RS(NewReac)%Product(pos)%Koeff + group_sigma(k)*sigma(current%iSpecies)*current%Koeff

           END IF

         END DO
       END DO

       ! assemble new line 3
       Line3Template = ReactionSystem(iR_old)%Line3
       currentColon = INDEX(Line3Template,':')
       ! maintain type of constant (before first colon)
       Line3 = Line3Template(1:currentColon)
       DO k=1,SIZE(l_RS(NewReac)%Constants)
         ! isolate next parameter name (between whitespace after a number and a colon)
         IF ( k==1 ) THEN
           currentSpace = currentColon + 1
         ELSE
           currentSpace = currentColon + INDEX(ADJUSTL(Line3Template(currentColon+1:)),' ')
           currentSpace = currentSpace + LEN_TRIM(Line3Template(currentColon+1:)) - LEN_TRIM(ADJUSTL(Line3Template(currentColon+1:)))
         END IF
         currentColon = currentColon + INDEX(Line3Template(currentColon+1:),':')

         ! add parameter name and parameter value to Line3
         Line3 = TRIM(ADJUSTL(Line3))//'  '//TRIM(ADJUSTL(Line3Template(currentSpace:currentColon)))
         Line3 = TRIM(ADJUSTL(Line3))//' '//TRIM(ADJUSTL(Real2CharExp(l_RS(NewReac)%Constants(k))))
       END DO

       ! collect remaining properties
       l_RS(NewReac)%Type           = ReactionSystem(iR_old)%Type
       l_RS(NewReac)%TypeConstant   = ReactionSystem(iR_old)%TypeConstant
       l_RS(NewReac)%Line3Template  = ReactionSystem(iR_old)%Line3
       l_RS(NewReac)%Line4          = ReactionSystem(iR_old)%Line4
       l_RS(NewReac)%Line3          = Line3

       NewReac = NewReac + 1

     END SUBROUTINE LumpReacs

     SUBROUTINE NewProperties(group, Spc_Old2New, sigma)
       TYPE(lumping_group), POINTER, INTENT(IN) :: group
       CHARACTER(LenName) , DIMENSION(:) :: Spc_Old2New
       REAL(dp)           , DIMENSION(:) :: sigma
       
       INTEGER, ALLOCATABLE :: Spc_array(:)
       REAL(dp), ALLOCATABLE :: sigma_array(:)
       INTEGER :: i, iR, iSpc
       CHARACTER(lenName) :: NewName

       CALL LL2Array(group%spc,Spc_array)
       CALL LL2Array_weight(group%spc,sigma_array)

       ! the new lumped species gets the name of the first species in the lumping group
       iSpc = Spc_array(1)
       NewName = y_name(iSpc)

       DO i=1,SIZE(Spc_array)
         Spc_Old2New(Spc_array(i)) = NewName 
         sigma(Spc_array(i))       = sigma_array(i) 
       END DO

     END SUBROUTINE NewProperties

  END SUBROUTINE build_LumpedReacSys


  SUBROUTINE compare_k(similar_k, reac_equivs, iSpcReacs, jSpcReacs)
    ! input:
    !
    LOGICAL, INTENT(OUT) :: similar_k
    !
    INTEGER, DIMENSION(:), INTENT(IN) :: reac_equivs, iSpcReacs, jSpcReacs
    !
    ! end input
 
    INTEGER :: i, reac, iR, jR, nReacs
    TYPE(ReactionStruct_T) :: iReac, jReac

    similar_k = .TRUE.

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


  SUBROUTINE put_tree(tree,id,weight,given_id_char)
    CLASS(linked_tree) :: tree
    
    INTEGER, INTENT(IN) :: id, weight
    CHARACTER(LenName), OPTIONAL :: given_id_char

    TYPE(linked_tree), ALLOCATABLE :: temp(:)
    INTEGER :: nChildren, i

    IF ( .NOT. ALLOCATED(tree%child) ) THEN ! first child
      ALLOCATE(tree%child(1))
      IF (ALLOCATED(tree%weight)) DEALLOCATE(tree%weight)
      ALLOCATE(tree%weight(1))
      tree%child(1)%id = id
      tree%weight(1) = weight
      IF ( PRESENT(given_id_char) ) THEN
        tree%child(1)%id_char = given_id_char
      END IF
    ELSE ! add next child
      ! increase size of child array by one
      nChildren = SIZE(tree%child)
      ALLOCATE(temp(nChildren))
      temp = tree%child
      DO i=1,nChildren
        CALL tree%child(i)%free()
      END DO
      DEALLOCATE(tree%child)
      ALLOCATE(tree%child(nChildren+1))
      tree%child(1:nChildren) = temp(:)

      ! apply desired values to (new) last child
      tree%child(nChildren+1)%id = id
      tree%weight = [ tree%weight , weight ]
      IF ( PRESENT(given_id_char) ) THEN
        tree%child(nChildren+1)%id_char = given_id_char
      END IF
      DO i=1,nChildren
        CALL temp(i)%free()
      END DO
      DEALLOCATE(temp)
    END IF
  END SUBROUTINE put_tree 

  RECURSIVE SUBROUTINE free_tree(tree)
    CLASS(linked_tree) :: tree

    INTEGER :: i
    
    IF ( ALLOCATED(tree%child) ) THEN
      DO i=1,SIZE(tree%child)
        CALL free_tree(tree%child(i))
      END DO
      DEALLOCATE(tree%child)
    END IF
    IF ( ALLOCATED(tree%weight) ) DEALLOCATE( tree%weight )
    IF ( ALLOCATED(tree%relevant) ) DEALLOCATE( tree%relevant )
    !tree%id=-1
    !tree%id_char=''

  END SUBROUTINE free_tree
  
  RECURSIVE SUBROUTINE WriteTree(tree,step_given)
    TYPE(linked_tree) :: tree
    INTEGER, OPTIONAL :: step_given

    INTEGER :: i, step

    IF (PRESENT(step_given)) THEN
      step = step_given
    ELSE
      step=1
    END IF
    
    WRITE (*,'(I5)',ADVANCE='NO') tree%id
    IF ( tree%id_char /= '' ) THEN
      WRITE (*,'(A)',ADVANCE = 'NO') '/'//TRIM(ADJUSTL(tree%id_char))
    END IF

    IF ( ALLOCATED(tree%child) ) THEN
      WRITE (*,'(A1,I1,A)',ADVANCE='NO') '(',step,' '
      DO i=1,SIZE(tree%child)
        CALL WriteTree(tree%child(i),step+1)
      END DO  
      WRITE (*,'(A,I2,A)',ADVANCE='NO') ' ',step,')'
    END IF
    IF ( step/=1 ) THEN
      WRITE (*,'(A2)',ADVANCE='NO') ', '
    ELSE
      WRITE(*,*) ''
    END IF
  END SUBROUTINE WriteTree


  RECURSIVE SUBROUTINE FindLeafs(tree,leafs)
    TYPE(linked_tree) :: tree
    TYPE(linked_list) :: leafs

    INTEGER :: i

    IF (ALLOCATED(tree%child)) THEN
      DO i=1,SIZE(tree%child)
        CALL FindLeafs(tree%child(i),leafs)
      END DO
    ELSE
      CALL leafs%put(tree%id)
    END IF
  END SUBROUTINE FindLeafs


  SUBROUTINE NextLumpingGroup(group, Spc, SpcReacs, n_SpcReacs)
    INTEGER :: Spc, n_SpcReacs
    INTEGER, ALLOCATABLE :: SpcReacs(:)
    TYPE(lumping_group), POINTER :: group

    TYPE(lumping_group), POINTER :: temp => NULL()
    LOGICAL :: found_spc, reached_end

    found_spc = .FALSE.
    reached_end = .FALSE.   

    ! find reactions in which iSpc occurs as educt
    CALL Spc_involving_reacs(SpcReacs,n_SpcReacs,Spc)  

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


  RECURSIVE SUBROUTINE put_ll(list, id, id_char, weight)
    CLASS(linked_list) , INTENT(INOUT) :: list
    INTEGER                            :: id
    CHARACTER(LenName) , OPTIONAL      :: id_char
    REAL(dp)           , OPTIONAL      :: weight
   
    IF ( .NOT. ASSOCIATED(list%next) .AND. list%id .EQ. -1 ) THEN
      list%id=id
      IF (PRESENT(id_char)) list%id_char = id_char
      IF (PRESENT(weight))  list%weight  = weight
    ELSE
      IF ( .NOT. ASSOCIATED(list%next) ) ALLOCATE(list%next)

      IF (PRESENT(id_char) .AND. PRESENT(weight)) THEN
        CALL put_ll(list%next,id,id_char,weight)
      ELSE IF (PRESENT(id_char)) THEN
        CALL put_ll(list%next,id,id_char=id_char)
      ELSE IF (PRESENT(weight)) THEN
        CALL put_ll(list%next,id,weight=weight)
      ELSE
        CALL put_ll(list%next,id)
      END IF
    END IF
  END SUBROUTINE put_ll

  RECURSIVE SUBROUTINE free_ll(list)
    CLASS(linked_list), INTENT(inout) :: list
    IF (ASSOCIATED(list%next)) THEN
       CALL free_ll(list%next)
       DEALLOCATE(list%next)
    END IF
    list%next => NULL()
    list%id=-1
    list%id_char=''
    list%weight=-1.0_dp
  END SUBROUTINE free_ll

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
      
      !Spc_char = get_Spc_char(Spc,SpcReacs(1))
      Spc_char = y_name(Spc)

      CALL group%spc%put(Spc,Spc_char,ONE)
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

  SUBROUTINE lumping_group_add(group, Spc, reac_equivs, SpcReacs, Conc_ratio)
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
    !  ratio of concentrations jSpc/iSpc
    REAL(dp) :: Conc_ratio
    !
    ! end input

    INTEGER :: i
    CHARACTER(LenName) :: Spc_char

    !Spc_char = get_Spc_char(Spc,SpcReacs(1))
    Spc_char = y_name(Spc)
    CALL group%spc%put(Spc,Spc_char,Conc_ratio)

    DO i=1,SIZE(reac_equivs)
      CALL group%reactions(i)%put(SpcReacs(reac_equivs(i)))
    END DO

    group%nspc = group%nspc + 1
    
  END SUBROUTINE lumping_group_add



! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!                               I/O
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
        WRITE(989,'(A)') '# '//TRIM(ADJUSTL(group%spc%id_char))
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
        Koeff_Str = Real2CharDec(Koeff)
        WRITE (len_Str,'(I3)') LEN_TRIM(ADJUSTL(Koeff_Str))+LEN_TRIM(Species)+1
        WRITE (989,'(A'//TRIM(len_Str)//')',ADVANCE='NO') TRIM(ADJUSTL(Koeff_Str))//' '//TRIM(Species)
        IF ( j<nEducts ) THEN
          WRITE (989,'(A3)',ADVANCE='NO') ' + '
        ELSE 
          WRITE (989,'(A3)',ADVANCE='NO') ' = '
        END IF
      END DO
     
!counter2 = 0 
      DO j=1,nProducts
        Koeff = l_RS(i)%Product(j)%Koeff
        Species = l_RS(i)%Product(j)%Species
        Koeff_Str = Real2CharDec(Koeff)
        WRITE (len_Str,*) LEN_TRIM(ADJUSTL(Koeff_Str))+LEN_TRIM(Species)+1
        WRITE (989,'(A'//TRIM(len_Str)//')',ADVANCE='NO') TRIM(ADJUSTL(Koeff_Str))//' '//TRIM(Species)
        IF ( j<nProducts ) THEN
          WRITE (989,'(A3)',ADVANCE='NO') ' + '
        END IF
      END DO
    
      WRITE (989,*) 
      WRITE (989,'(A)') TRIM(l_RS(i)%Line3)
      WRITE (989,'(A)') TRIM(l_RS(i)%Line4)
    END DO

    CLOSE(989)    
  END SUBROUTINE Print_LumpedSysFile

  SUBROUTINE EmisTreeToMatlab(tree)
    TYPE(linked_tree), TARGET :: tree

    CHARACTER(LenName) :: Spc_name, num
    TYPE(linked_list), TARGET  :: starts, ends, weights, mask, leafs
    TYPE(linked_list), POINTER :: temp
    INTEGER :: treeID

    CALL mask%put(tree%id)
    CALL CollectNodesEdges(tree,starts,ends,weights,mask)
    CALL FindLeafs(tree,leafs)

    OPEN(UNIT=989,FILE='LUMPING/EmisTree_'//TRIM(ADJUSTL(tree%id_char))//'.txt',STATUS='UNKNOWN')
    
    ! write start points of edges
    temp=>starts
    WRITE (989,'(A)',ADVANCE='NO') 't = ['
    DO
      treeID = SearchLinkedList(mask,id=temp%id)
      WRITE(num,*) treeID 
      WRITE(989,'(A)',ADVANCE='NO') TRIM(ADJUSTL(num))//' '
      IF ( ASSOCIATED(temp%next) ) THEN
        temp=>temp%next
      ELSE
        EXIT
      END IF
    END DO
    WRITE(989,'(A)',ADVANCE='NO') '];'
    WRITE(989,'(A)') ''

    ! write end points of edges
    temp=>ends
    WRITE (989,'(A)',ADVANCE='NO') 's = ['
    DO 
      treeID = SearchLinkedList(mask,id=temp%id)
      WRITE(num,*) treeID 
      WRITE(989,'(A)',ADVANCE='NO') TRIM(ADJUSTL(num))//' '
      IF ( ASSOCIATED(temp%next) ) THEN
        temp=>temp%next
      ELSE
        EXIT
      END IF
    END DO
    WRITE(989,'(A)',ADVANCE='NO') '];'
    WRITE(989,'(A)') ''

    ! write matlab commands
    WRITE(989,'(A)') 'G=digraph(s,t);'

    ! write NodeLabel
    temp=>mask
    WRITE (989,'(A)',ADVANCE='NO') 'NodeLabels = {'
    DO
      IF (temp%id == tree%id) THEN
        num = 'Start: '//TRIM(ADJUSTL(y_name(temp%id)))
      ELSE IF (SearchLinkedList(leafs,id=temp%id) > 0) THEN
        num = 'End: '//TRIM(ADJUSTL(y_name(temp%id)))
      ELSE
        num = y_name(temp%id)
      END IF
      WRITE(989,'(A)',ADVANCE='NO') ' '''//TRIM(ADJUSTL(num))//''' '
      IF ( ASSOCIATED(temp%next) ) THEN
        temp=>temp%next
      ELSE
        EXIT
      END IF
    END DO
    WRITE(989,'(A)',ADVANCE='NO') '};'
    WRITE(989,'(A)') ''
    
    ! write EdgeLabel
    temp=>weights
    WRITE (989,'(A)',ADVANCE='NO') 'EdgeLabelsUnsorted = {'
    DO 
      WRITE(num,*) temp%id
      WRITE(989,'(A)',ADVANCE='NO') ' '''//TRIM(ADJUSTL(num))//'/'//TRIM(ADJUSTL(temp%id_char))//''' '
      IF ( ASSOCIATED(temp%next) ) THEN
        temp=>temp%next
      ELSE
        EXIT
      END IF
    END DO
    WRITE(989,'(A)',ADVANCE='NO') '};'
    WRITE(989,'(A)') ''
    WRITE(989,'(A)') '[~,Perm] = sort(s);'
    WRITE(989,'(A)') 'EdgeLabels=EdgeLabelsUnsorted(Perm);'
    WRITE(989,'(A)') 'plot(G,''NodeLabel'',NodeLabels,''EdgeLabel'',EdgeLabels);'


    CLOSE(989)

  END SUBROUTINE EmisTreeToMatlab


  FUNCTION Read_Emis_Families(FileName) RESULT(EmisFams)
    USE Reac_Mod, ONLY: y_name
    USE ChemSys_Mod,  ONLY: PositionSpeciesAll
    ! OUT:
    TYPE(Emis_Family), ALLOCATABLE :: EmisFams(:)
    ! IN:
    CHARACTER(*) :: FileName
    ! TEMP:
    INTEGER :: iPos, j, nFams, current, FirstBlank, SecondBlank, temp_n
    REAL(dp) :: NewPortion
    CHARACTER(LenLine) :: Line, temp
    CHARACTER(LenName) :: NewMember
    CHARACTER(LenName), ALLOCATABLE :: temp_arr(:)

    OPEN(UNIT=99,FILE=TRIM(ADJUSTL(FileName)),STATUS='UNKNOWN')
    REWIND(99)

    iPos = FindSection(99,'EMIS_FAMILIES')
    IF ( iPos==-1 ) THEN
      ALLOCATE(EmisFams(1))
      RETURN
    END IF
    
    ! count families
    nFams = 0
    DO
      iPos = FindSection(99,'FAMILY')
      IF ( iPos /= -1 ) THEN
        nFams = nFams + 1
      ELSE
        EXIT
      END IF
    END DO   
    ALLOCATE(EmisFams(nFams))

    ! walk to the first family
    REWIND(99)
    iPos = FindSection(99,'EMIS_FAMILIES')
    current = 0
    ALLOCATE(temp_arr(0))
    DO
      READ(99,'(A)') Line; Line = TRIM(ADJUSTL(Line))
      IF ( Line == 'END_EMIS_FAMILIES' ) THEN
        EXIT
      ELSE IF ( Line(1:6) == 'FAMILY' ) THEN ! go to next family
        
        current = current + 1
     
        ! FirstBlank is between FAMILY and Name, SecondBlank is between Name and Emission value
        FirstBlank  = INDEX(Line,' ')
        SecondBlank = FirstBlank + INDEX(ADJUSTL(Line(FirstBlank:)),' ')
     
        ! collect name and emission value of family
        EmisFams(current)%Name = TRIM(ADJUSTL(Line(FirstBlank:SecondBlank)))
        temp = TRIM(ADJUSTL(Line(SecondBlank:)))
        READ(temp,*) EmisFams(current)%Emis_val
        ALLOCATE(EmisFams(current)%Members(0),EmisFams(current)%Members_char(0),EmisFams(current)%Portions(0))
      
      ELSE IF ( Line /= '' ) THEN ! add new family member
        
        ! FirstBlank is between family member (species) and emission portion
        FirstBlank = INDEX(Line,' ')
        
        NewMember = TRIM(ADJUSTL(Line(:FirstBlank)))
        temp = TRIM(ADJUSTL(Line(FirstBlank:)))
        READ(temp,*) NewPortion

        ! increase size of members_char by one, [ , ]-syntax doesnt work
        ! optimize maybe
        temp_n = SIZE(EmisFams(current)%Members_char)
        DEALLOCATE(temp_arr); ALLOCATE(temp_arr(temp_n+1))
        temp_arr(1:temp_n) = EmisFams(current)%Members_char
        temp_arr(temp_n+1) = NewMember
        DEALLOCATE(EmisFams(current)%Members_char); ALLOCATE(EmisFams(current)%Members_char(temp_n+1))

        EmisFams(current)%Members_char = temp_arr
        EmisFams(current)%Members      = [ EmisFams(current)%Members     , PositionSpeciesAll(NewMember) ]
        EmisFams(current)%Portions     = [ EmisFams(current)%Portions    , NewPortion                    ]
      
      END IF
    END DO
    
   ! DO j=1,nFams
   !   WRITE(*,*) 'Family Name ',EmisFams(j)%Name
   !   WRITE(*,*) 'Emis val: ',Real2CharExp(EmisFams(j)%Emis_val)
   !   WRITE(*,*) 'Members char: ',EmisFams(j)%Members_char
   !   WRITE(*,*) 'Members: ',EmisFams(j)%Members
   !   WRITE(*,*) 'Portions: ',EmisFams(j)%Portions
   !   WRITE(*,*) ''
   !   WRITE(*,*) ''
   !   WRITE(*,*) ''
   ! END DO

    CLOSE(99)

  END FUNCTION Read_Emis_Families


  FUNCTION Read_Preserve_Spc(FileName) RESULT(Idx)
    USE ChemSys_Mod,  ONLY: PositionSpeciesAll
    USE Control_Mod,  ONLY: LenLine
    USE Reac_Mod,     ONLY: y_name
    ! This routine will read the species that are not to be lumped.
    ! OUT:
    INTEGER, ALLOCATABLE :: Idx(:)
    ! IN:
    CHARACTER(*)  :: FileName
    ! TEMP:
    INTEGER       :: iPos, j
    CHARACTER(LenLine) :: Line
     
    OPEN(UNIT=99,FILE=TRIM(ADJUSTL(FileName)),STATUS='UNKNOWN')
    REWIND(99)

    ALLOCATE(Idx(0))
    iPos = FindSection(99,'PRESERVE')
    IF ( iPos /= -1 ) THEN
      DO
        READ(99,'(A)') Line;  Line = ADJUSTL(Line)
        IF (TRIM(Line) == 'END_PRESERVE') EXIT
        IF (TRIM(Line) == '' ) CYCLE
        iPos = PositionSpeciesAll(Line)
        IF ( iPos > 0 ) Idx = [ Idx , iPos ]
      END DO
    END IF
    CLOSE(99)
    
    WRITE(*,*)
    WRITE(*,777) 'Species that are preserved while lumping:'
    WRITE(*,*)
    IF ( SIZE(Idx)==0 ) THEN
      WRITE(*,*) '    none'
    ELSE
      DO j=1,SIZE(Idx)
        WRITE(*,'(10X,A,I0,A,I6,5X,A)') '    Preserved Spc ',j,' = ',Idx(j),TRIM(y_name(Idx(j)))
      END DO
    END IF
    WRITE(*,*); WRITE(*,*)

    
    777 FORMAT(10X,A)
  END FUNCTION Read_Preserve_Spc

  SUBROUTINE WriteLumpingTimes()
    CHARACTER(8) ::  unit1  = '',   &
    &                unit2  = '',   &
    &                unit3  = ''
   
    CALL ConvertTime(TimeFluxRead,unit1)
    CALL ConvertTime(TimeConcRead,unit2)
    CALL ConvertTime(TimeLumping,unit3)
    WRITE(*,*)
    WRITE(*,'(11X,A,1X,F10.4,A)') 'Time reading flux-dataset = ', TimeFluxRead, unit1
    WRITE(*,'(11X,A,1X,F10.4,A)') 'Time reading conc-dataset = ', TimeConcRead, unit2
    WRITE(*,'(11X,A,1X,F10.4,A)') 'Time lumping procedure    = ', TimeLumping, unit3
    WRITE(*,*);
  END SUBROUTINE WriteLumpingTimes

  SUBROUTINE WriteLumpingGroups(first_group)

    ! input:
    !
    TYPE(lumping_group) :: first_group

    INTEGER :: group_count

    TYPE(linked_list) :: unlumped_Spc

    group_count = 0

    CALL WriteLumpingGroups_rec(first_group,group_count,unlumped_Spc)

  END SUBROUTINE WriteLumpingGroups

  RECURSIVE SUBROUTINE WriteLumpingGroups_rec(group,group_count,unlumped_Spc)
   
    ! input:
    !
    TYPE(lumping_group), TARGET :: group
    !
    INTEGER :: group_count
    !
    TYPE(linked_list), TARGET :: unlumped_Spc
    
    TYPE(linked_list), POINTER :: temp_list
    INTEGER, ALLOCATABLE :: Spcs_offset(:), reacs_array(:)
    CHARACTER, ALLOCATABLE :: Spcs_unlumped(:)
    REAL(dp), ALLOCATABLE :: sigmas(:)
    
    INTEGER :: i, j, iR, len_Spcs, n_Spcs
    CHARACTER(LEN=6) :: i_str, iR_str, id_str

    CALL LL2Array_weight(group%spc,sigmas)
    
    IF ( group%nspc > 1 ) THEN
      group_count = group_count + 1
      
      id_str=''
      WRITE (id_str,'(I6)') group_count
      WRITE (*,*) ''
      WRITE (*,*) ''
      WRITE (*,*) id_str//'. lumping group '
      WRITE (*,*) ''
      WRITE (*,*) '      -number of species:',group%nspc
      WRITE (*,*) '      -sigmas:', sigmas
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
      END DO
    ELSE !n_Spcs=1
     ! add Spc to unlumped species
     CALL unlumped_Spc%put(group%spc%id,group%spc%id_char)
    END IF

    IF (ASSOCIATED(group%next)) THEN
      CALL WriteLumpingGroups_rec(group%next, group_count, unlumped_Spc)
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
    END IF
  END SUBROUTINE WriteLumpingGroups_rec

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

  SUBROUTINE LL2Array_weight(list, list_array)
    TYPE(linked_list), TARGET, INTENT(IN) :: list
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: list_array(:)

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
      list_array(i)=temp%weight
      temp=>temp%next
    END DO
  END SUBROUTINE LL2Array_weight

  FUNCTION SearchLinkedList(list,id,id_char) RESULT(pos)
    TYPE(linked_list),  INTENT(IN), TARGET    :: list
    INTEGER,            INTENT(IN), OPTIONAL  :: id
    CHARACTER(LenName), INTENT(IN), OPTIONAL  :: id_char

    INTEGER :: pos

    TYPE(linked_list), POINTER :: temp

    temp=>list
    pos=1
    DO
     IF ( PRESENT(id) ) THEN
       IF ( temp%id==id ) THEN
         EXIT
       END IF
     END IF
     IF ( PRESENT(id_char) ) THEN
       IF ( temp%id_char == id_char ) THEN
         EXIT
       END IF
     END IF
     IF ( ASSOCIATED(temp%next) ) THEN
       temp=>temp%next
       pos = pos+1
     ELSE
       pos=0
       EXIT
     END IF
    END DO

  END FUNCTION SearchLinkedList

  SUBROUTINE Logo3()

    WRITE (*,*) ''
    WRITE (*,*) '                  _________________________ '
    WRITE (*,*) '                 |                         |'
    WRITE (*,*) '                 |         LUMPING         |'
    WRITE (*,*) '                 |_________________________|'
    WRITE (*,*) ''
    WRITE (*,*) ''
  
  END SUBROUTINE Logo3

  FUNCTION Real2CharExp(x) RESULT(xChar)
    REAL(dp) :: x
    CHARACTER(LenLine) :: xChar

    INTEGER :: magofmag, decim, aux
    CHARACTER(LenLine) :: nChars, relevant_digits, exp_char
    REAL(dp) :: magnitude, temp

    IF ( x/=0 ) THEN
      magnitude = LOG(ABS(x))/LOG(TEN)
      IF ( INT(magnitude) /= 0 ) THEN
        magofmag = INT(LOG(ABS(magnitude))/LOG(TEN))
      ELSE
        magofmag = 0
      END IF
    ELSE
      magnitude = 0
      magofmag  = 0
    END IF
    WRITE (exp_char,*) magofmag + 1
    temp = x
    temp = temp*(TEN**(mONE*CEILING(magnitude)))
    decim = 0
    ! find number of relevant digits after comma (decim)
    DO WHILE ( ABS(temp-ANINT(temp))>nano .AND. decim<16 )
      temp = temp*10
      decim=decim+1
    END DO
    decim = MAX(decim,1)

    IF (x<ZERO) THEN
      aux=6
    ELSE
      aux=5
    END IF
    WRITE (nChars,*) decim + magofmag + aux
    WRITE (relevant_digits,*) decim
    
    WRITE (xChar,'(E'//TRIM(ADJUSTL(nChars))//'.'//TRIM(ADJUSTL(relevant_digits))//'E'//TRIM(ADJUSTL(exp_char))//')') x

  END FUNCTION Real2CharExp

  FUNCTION Real2CharDec(x) RESULT(xChar)
    REAL(dp) :: x
    CHARACTER(LenLine) :: xChar

    INTEGER :: mag, decim, aux
    CHARACTER(LenLine) :: nChars, relevant_digits
    REAL(dp) :: magnitude, temp

    IF ( x/=0 ) THEN
      magnitude = LOG(ABS(x))/LOG(TEN)
    ELSE
      magnitude = 0
    END IF
    temp = x
    decim = 0
    ! find number of relevant digits after comma (decim)
    DO WHILE ( ABS(temp-ANINT(temp))>micro .AND. decim<16 )
      temp = temp*10
      decim=decim+1
    END DO
    decim = MAX(decim,1)
    mag = INT(magnitude) + 1

    IF (x<ZERO) THEN
      aux=2
    ELSE
      aux=1
    END IF
    WRITE (nChars,*) decim + MAX(mag,1) + aux
    WRITE (relevant_digits,*) decim

    WRITE (xChar,'(F'//TRIM(ADJUSTL(nChars))//'.'//TRIM(ADJUSTL(relevant_digits))//')') x

    !WRITE(*,*) 'x=',x
    !WRITE(*,*) 'xChar=',xChar
    !WRITE(*,*) ''

  END FUNCTION Real2CharDec

  FUNCTION FindSection(iUnit,SectionName,iLine_given) RESULT(iLine)
    USE Control_Mod, ONLY: lenLine
    ! OUT:
    INTEGER      :: iLine
    ! IN: 
    INTEGER      :: iUnit
    CHARACTER(*) :: SectionName
    INTEGER, OPTIONAL :: iLine_given
    ! TEMP: 
    CHARACTER(LenLine) :: Line
    INTEGER      :: io_err
   
    iLine = 0
    IF ( PRESENT(iLine_given) ) THEN ! walk to iLine_given to start searching procedure from there
      DO WHILE ( iLine < iLine_given )
        READ(iUnit,'(A)',IOSTAT=io_err) Line
        iLine = iLine + 1
      END DO
    END IF
    DO
      READ(iUnit,'(A)',IOSTAT=io_err) Line
      IF (io_err==0) THEN
        IF ( INDEX(TRIM(ADJUSTL(Line)),SectionName) /= 0 )  EXIT
        iLine = iLine + 1
      ELSEIF (io_err>0.OR.io_err<0) THEN
        iLine = -1
        EXIT
      END IF
    END DO
     777 FORMAT(10X,A)
  END FUNCTION FindSection

END MODULE Lumping_Mod
