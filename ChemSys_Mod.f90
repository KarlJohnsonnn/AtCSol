!===================================================================================!
!                                                                                   !
!                     Module for reading the chemical system                        !
!                 file .sys and building the coeficient matrices                    !
!                                                                                   !
!===================================================================================!
MODULE Chemsys_Mod
  !
  USE Kind_Mod
  USE Meteo_Mod
  USE Sparse_Mod
  USE String_Mod
  USE LexicalStringSort
  USE hashtbl
  USE UniRnk_Mod
  USE Control_Mod
  USE Reac_Mod
  USE InputTool_Mod
  USE NetCDF_Mod
  !
  IMPLICIT NONE
  !
  !
  INTEGER, PARAMETER ::  maxLENinActDuct=9
  ! 
  TYPE Duct_T
    CHARACTER(LenName) :: Species=''
    CHARACTER(LenType) :: Type
    REAL(dp)           :: Koeff
    INTEGER            :: iSpecies=0
  END TYPE Duct_T

  TYPE Special_T
    INTEGER                         :: nVariables = 0
    INTEGER,            ALLOCATABLE :: iVariables(:)
    CHARACTER(LenName), ALLOCATABLE :: cVariables(:)
    CHARACTER(LenLine)              :: Formula = ''
    LOGICAL                         :: Temp = .FALSE.
  END TYPE Special_T
  !
  ! LIST FORM
  TYPE Reaction_T
    CHARACTER(LenType)        :: Type, TypeConstant
    CHARACTER(LenName)        :: Comment
    CHARACTER(LenLine)        :: Line1, Line2, Line3, Line4
    CHARACTER(LenName)        :: Factor
    TYPE(Duct_T),     POINTER :: Educt(:)=>NULL(), Product(:)=>NULL()
    REAL(dp),     ALLOCATABLE :: Constants(:)
    TYPE(Duct_T),     POINTER :: InActEduct(:)=>NULL(), InActProduct(:)=>NULL()
    TYPE(Special_T)           :: Special
    INTEGER                   :: nInActEd=0, nInActPro=0
    TYPE(Reaction_T), POINTER :: Next=>NULL()
  END TYPE Reaction_T
  !
  ! ARRAY FORM
  TYPE ReactionStruct_T
    CHARACTER(LenType)              :: Type,  TypeConstant
    CHARACTER(LenLine)              :: Line1='' , Line2='' , Line3='', Line4=''
    LOGICAL                         :: bR = .FALSE. , brX = .FALSE. 
    CHARACTER(LenName)              :: Factor = ''
    CHARACTER(LenName)              :: Comment = ''
    CHARACTER(2)                    :: direction = ''
    REAL(dp)                        :: SumAqCoef     
    TYPE(Special_T)                 :: Special
    TYPE(Duct_T)  ,     ALLOCATABLE :: Educt(:), Product(:)
    REAL(dp),           ALLOCATABLE :: Constants(:)
    REAL(dp),           ALLOCATABLE :: LowConst(:), HighConst(:), TroeConst(:) ! combustion press dep reactions
    REAL(dp),           ALLOCATABLE :: InActEduct(:), InActProduct(:)
    INTEGER                         :: nInActEd = 0, nInActPro = 0, nActEd = 0, nActPro = 0
    INTEGER                         :: nConst = 0
    INTEGER                         :: HenrySpc = 0
    LOGICAL                         :: TB = .FALSE. , TBextra=.FALSE.
    INTEGER,            ALLOCATABLE :: TBidx(:)
    CHARACTER(LenName), ALLOCATABLE :: TBspc(:)
    REAL(dp),           ALLOCATABLE :: TBalpha(:)
    CHARACTER(LenName), ALLOCATABLE :: InActEductSpc(:), InActProductSpc(:)
  END TYPE ReactionStruct_T
  !
  !
  TYPE ListReaction_T
    TYPE(Reaction_T), POINTER :: Start=>NULL()
    TYPE(Reaction_T), POINTER :: End=>NULL()
    INTEGER :: LenList=0
  END TYPE ListReaction_T
  !
  TYPE Species_T
    CHARACTER(LenName) :: Species=''
    LOGICAL            :: isHenry=.FALSE.
    REAL(dp)     :: Hf=0.0d0, Gf=0.0d0, Cp=0.0d0
  END TYPE Species_T

  

  
  TYPE(Reaction_T), POINTER   :: System
  TYPE(ListReaction_T), SAVE  :: ListRGas, ListRHenry, ListRAqua,        &
  &                              ListRDiss, ListRSolid, ListRPartic,     &
  &                              ListRMicro
  !
  TYPE(hash_tbl_sll)          :: ListAqua, ListGas, ListSolid,           &
  &                              ListPartic, ListNonReac, ListAtoms
  TYPE(hash_tbl_sll)          :: ListFamilies
  !
  TYPE(Species_T), ALLOCATABLE, TARGET :: ListAqua2(:), ListGas2(:),     &
  &                               ListSolid2(:), ListPartic2(:), &
  &                               ListNonReac2(:)
  INTEGER :: InputUnit=10
  INTEGER, PARAMETER :: MaxEduct=10
  INTEGER, PARAMETER :: MaxProduct=10
  !
  CHARACTER(33), PARAMETER :: SetSpecies='ABCDEFGHIJKLMNOPQRSTUVWXYZapsc[]()=+*'
  CHARACTER(14), PARAMETER :: SetConstants='ABCDEFGKINMOR/'
  CHARACTER(12), PARAMETER :: SetExponent='0123456789+-'

  TYPE Element_T
    CHARACTER(5) :: Element=''
  END TYPE Element_T

  TYPE(Element_T) :: Elements(11)=(/       &
  &                    Element_T('(')      &
  &                    ,Element_T(')')     &
  &                    ,Element_T('exp')   &
  &                    ,Element_T('+')     &
  &                    ,Element_T('-')     &
  &                    ,Element_T('*')     &
  &                    ,Element_T('/')     &
  &                    ,Element_T('**')    &
  &                    ,Element_T('abs')   &
  &                    ,Element_T('sqrt')  &
  &                    ,Element_T('log')   /)
 
  
  INTEGER :: nsr                        ! # activ species + all Reactions

 
  !
  CHARACTER(20) :: Filename
  CHARACTER(20) :: IniName
  !
  REAL(dp), PARAMETER :: RGas=8.3145d0
  REAL(dp), PARAMETER :: TRef=280.0d0 !298.15d0
  !
  TYPE(Reaction_T), POINTER :: Current
  TYPE(ReactionStruct_T), ALLOCATABLE :: ReactionSystem(:)
  TYPE(ListReaction_T), ALLOCATABLE :: CompleteReactionList(:)
  !
  !
  REAL(dp), ALLOCATABLE :: Emis(:)          & ! emission values
  !&                            , InitValAct(:)    & ! initial values activ spc
  &                      , InitValInAct(:)    ! initial values inactiv spc
  !
  !
  CHARACTER(LenName), ALLOCATABLE :: RO2spcG(:) , RO2spcA(:)
  INTEGER, ALLOCATABLE :: RO2idxG(:) , RO2idxA(:)
  !
  !
  !REAL(dp) :: aH2O
  !
  REAL(dp), ALLOCATABLE :: sumBAT(:)         ! sum_j=1,n_s (b_ij-a_ij),  i el. N_R

  INTEGER :: fNumber = 0
  !
  CONTAINS
  ! ------------------------------------
  ! -----------SUBROUTINEN--------------
  ! ------------------------------------
  !
  SUBROUTINE SortReactionList(ReacStructOut,ReacStructIn)
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStructIn(:)
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStructOut(:)
    !
    INTEGER :: i
    CHARACTER(20), ALLOCATABLE :: ReacTypeSorted(:)
    INTEGER, ALLOCATABLE :: iReacTypeSorted(:)
    !
    ! sort the reaction list --> TypeConstant
    ALLOCATE(ReacTypeSorted(neq))
    ALLOCATE(iReacTypeSorted(neq))  
    DO i=1,neq
      ReacTypeSorted(i)=ReacStructIn(i)%Type
    END DO
    CALL StringSort(ReacTypeSorted,iReacTypeSorted)
    ALLOCATE(ReacStructOut(neq))
    !
    DO i=1,SIZE(ReacStructIN)
      ReacStructOut(i)=ReacStructIn(iReacTypeSorted(i))
    END DO
    DEALLOCATE(ReacStructIn)
    DEALLOCATE(iReacTypeSorted)
  END SUBROUTINE SortReactionList
  !
  !
  SUBROUTINE ReadSpecies(Out)
    LOGICAL :: Out
    !
    CHARACTER(LenName) :: Species
    CHARACTER(LenType) :: Type
    INTEGER :: Pos

    READ(InputUnit,'(a100)',END=1) Species

    DO
      Pos = SCAN( Species , "'" )
      IF ( Pos > 0 ) THEN
        Species(Pos:) = Species(Pos+1:)
      ELSE
        EXIT
      END IF
    END DO
    IF ( Species /= '' ) CALL InsertSpecies(Species,Type)

    Out = .FALSE.
    GO TO 999

  1 CONTINUE
   
    Out = .TRUE.
999 CONTINUE
  END SUBROUTINE ReadSpecies


  SUBROUTINE ReadReaction_neu(Out)
    LOGICAL :: Out
    !
    INTEGER :: iLine,PosColon,Pos,is
    CHARACTER(LenLine) :: LocLine
    CHARACTER(LenLine) :: Line(1:4)
    CHARACTER(20) :: CLASS
    CHARACTER(40) :: TypeR
    INTEGER :: idxFAC
    
    iLine = 0
    DO

      READ( InputUnit , '(A400)' , IOSTAT=is ) LocLine
      idxFAC = INDEX(LocLine,'$')

      IF ( idxFAC > 0 ) THEN
        SELECT CASE (TRIM(LocLine(idxFAC:)))
          CASE ('$H2','$O2N2','$M','$O2','$N2','$H2O','$O2O2','$aH2O','$+M','$(+M)','$RO2','$RO2aq')
            IF ( TRIM(LocLine(idxFAC:)) == '$H2'    ) nr_FAC_H2    = nr_FAC_H2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2N2'  ) nr_FAC_O2N2  = nr_FAC_O2N2  + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$M'     ) nr_FAC_M     = nr_FAC_M     + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2'    ) nr_FAC_O2    = nr_FAC_O2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$N2'    ) nr_FAC_N2    = nr_FAC_N2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$H2O'   ) nr_FAC_H2O   = nr_FAC_H2O   + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2O2'  ) nr_FAC_O2O2  = nr_FAC_O2O2  + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$RO2'   ) nr_FAC_RO2   = nr_FAC_RO2   + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$RO2aq' ) nr_FAC_RO2aq = nr_FAC_RO2aq + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$aH2O'  ) nr_FAC_aH2O  = nr_FAC_aH2O  + 2
            nr_FACTOR = nr_FACTOR + 1
          CASE DEFAULT
            WRITE(*,*) '  Unknown FACTOR:  ', TRIM(LocLine(idxFAC:)), '   at Line:  ???'
        END SELECT
      END IF

      IF ( ABS(is) > 0 ) EXIT
      
      ! if no comment or blank line then
      IF ( ADJUSTL(LocLine(1:1)) /= '#'       .AND.  &
      &    ADJUSTL(LocLine(1:7)) /= 'COMMENT' .AND.  &
      &    LEN(TRIM(LocLine)) > 0 ) THEN
        iLine = iLine + 1
        Line(iLine) = LocLine
        IF ( iLine == 4 ) EXIT
      END IF

    END DO

    IF ( iLine >= 3 ) THEN
      Pos = SCAN(Line(1),'#')
      IF ( Pos > 0 ) Line(1) = Line(1)(:Pos-1)
      
      ! if new reaction line starts go one line back
      IF ( INDEX(Line(4),'CLASS') > 0 ) BACKSPACE(InputUnit)
      
      ! read the reaction TYPE
      PosColon = Index(Line(1),':')
      CLASS    = ADJUSTL(Line(1)(PosColon+1:))
      
      ! count the number of each reaction type
      SELECT CASE (CLASS)
        
        CASE ('GAS')      ! gaseous phase reactions

          nr_gas = nr_gas + 1
          CALL InsertReaction( ListRGas , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('PHOTO','PHOTO2','PHOTO3','PHOTAB','PHOTABC','PHOTMCM')
              IF ( TypeR == 'PHOTAB'  ) reac_par(iPHOTAB)%n_reac  = reac_par(iPHOTAB)%n_reac + 1
              IF ( TypeR == 'PHOTABC' ) reac_par(iPHOTABC)%n_reac = reac_par(iPHOTABC)%n_reac + 1
              IF ( TypeR == 'PHOTMCM' ) reac_par(iPHOTMCM)%n_reac = reac_par(iPHOTMCM)%n_reac + 1
              IF ( TypeR == 'PHOTO'   ) reac_par(iPHOTO)%n_reac   = reac_par(iPHOTO)%n_reac + 1
              IF ( TypeR == 'PHOTO2'  ) reac_par(iPHOTO2)%n_reac  = reac_par(iPHOTO2)%n_reac + 1
              IF ( TypeR == 'PHOTO3'  ) reac_par(iPHOTO3)%n_reac  = reac_par(iPHOTO3)%n_reac + 1
            CASE ('CONST')
              reac_par(iCONST)%n_reac = reac_par(iCONST)%n_reac + 1
            CASE ('TEMP','TEMP1','TEMP2','TEMP3','TEMP4')
              IF ( TypeR == 'TEMP' )  reac_par(iTEMP)%n_reac  = reac_par(iTEMP)%n_reac + 1
              IF ( TypeR == 'TEMP1' ) reac_par(iTEMP1)%n_reac = reac_par(iTEMP1)%n_reac + 1
              IF ( TypeR == 'TEMP2' ) reac_par(iTEMP2)%n_reac = reac_par(iTEMP2)%n_reac + 1
              IF ( TypeR == 'TEMP3' ) reac_par(iTEMP3)%n_reac = reac_par(iTEMP3)%n_reac + 1
              IF ( TypeR == 'TEMP4' ) reac_par(iTEMP4)%n_reac = reac_par(iTEMP4)%n_reac + 1
            CASE ('TROE','TROEF','TROEQ','TROEQF','TROEXP','TROEMCM')
              IF ( TypeR == 'TROE'    ) reac_par(iTROE)%n_reac    = reac_par(iTROE)%n_reac + 1
              IF ( TypeR == 'TROEF'   ) reac_par(iTROEF)%n_reac   = reac_par(iTROEF)%n_reac + 1 
              IF ( TypeR == 'TROEQ'   ) reac_par(iTROEQ)%n_reac   = reac_par(iTROEQ)%n_reac + 1
              IF ( TypeR == 'TROEQF'  ) reac_par(iTROEQF)%n_reac  = reac_par(iTROEQF)%n_reac + 1 
              IF ( TypeR == 'TROEXP'  ) reac_par(iTROEXP)%n_reac  = reac_par(iTROEXP)%n_reac + 1 
              IF ( TypeR == 'TROEMCM' ) reac_par(iTROEMCM)%n_reac = reac_par(iTROEMCM)%n_reac + 1
            CASE ('SPEC1','SPEC2','SPEC3','SPEC4','SPEC1MCM',  &
            &     'SPEC2MCM','SPEC3MCM','SPEC4MCM','SPEC5MCM', &
            &     'SPEC6MCM','SPEC7MCM','SPEC8MCM','SPEC9MCM'  )
              nr_G_spec = nr_G_spec + 1
              IF ( TypeR == 'SPEC1' ) reac_par(iSPEC1)%n_reac = reac_par(iSPEC1)%n_reac + 1
              IF ( TypeR == 'SPEC2' ) reac_par(iSPEC2)%n_reac = reac_par(iSPEC2)%n_reac + 1
              IF ( TypeR == 'SPEC3' ) reac_par(iSPEC3)%n_reac = reac_par(iSPEC3)%n_reac + 1
              IF ( TypeR == 'SPEC4' ) reac_par(iSPEC4)%n_reac = reac_par(iSPEC4)%n_reac + 1
              IF ( TypeR == 'SPEC1MCM' ) reac_par(iSPEC1MCM)%n_reac = reac_par(iSPEC1MCM)%n_reac + 1
              IF ( TypeR == 'SPEC2MCM' ) reac_par(iSPEC2MCM)%n_reac = reac_par(iSPEC2MCM)%n_reac + 1
              IF ( TypeR == 'SPEC3MCM' ) reac_par(iSPEC3MCM)%n_reac = reac_par(iSPEC3MCM)%n_reac + 1
              IF ( TypeR == 'SPEC4MCM' ) reac_par(iSPEC4MCM)%n_reac = reac_par(iSPEC4MCM)%n_reac + 1
              IF ( TypeR == 'SPEC5MCM' ) reac_par(iSPEC5MCM)%n_reac = reac_par(iSPEC5MCM)%n_reac + 1
              IF ( TypeR == 'SPEC6MCM' ) reac_par(iSPEC6MCM)%n_reac = reac_par(iSPEC6MCM)%n_reac + 1
              IF ( TypeR == 'SPEC7MCM' ) reac_par(iSPEC7MCM)%n_reac = reac_par(iSPEC7MCM)%n_reac + 1
              IF ( TypeR == 'SPEC8MCM' ) reac_par(iSPEC8MCM)%n_reac = reac_par(iSPEC8MCM)%n_reac + 1
              IF ( TypeR == 'SPEC9MCM' ) reac_par(iSPEC9MCM)%n_reac = reac_par(iSPEC9MCM)%n_reac + 1
            CASE ('S4H2O')
              reac_par(iS4H2O)%n_reac = reac_par(iS4H2O)%n_reac + 1 
            CASE ('T1H2O')
              reac_par(iT1H2O)%n_reac = reac_par(iT1H2O)%n_reac + 1 
            CASE ('SPECIAL')
              reac_par(iSPECIAL)%n_reac = reac_par(iSPECIAL)%n_reac + 1 
            CASE ('HOM1')
              reac_par(iHOM1)%n_reac = reac_par(iHOM1)%n_reac + 1 
            CASE DEFAULT
              WRITE(*,*) '  Unknown gaseous reaction: ', TypeR
          END SELECT

        CASE ('HENRY')        ! phase transfer pseudo-reactions 

          nr_henry = nr_henry + 1
          CALL InsertReaction( ListRHenry , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('TEMP','TEMP1','TEMP2','TEMP3')
              IF ( TypeR == 'TEMP' )  reac_par(iTEMP)%n_reac  = reac_par(iTEMP)%n_reac + 1
              IF ( TypeR == 'TEMP1' ) reac_par(iTEMP1)%n_reac = reac_par(iTEMP1)%n_reac + 1
              IF ( TypeR == 'TEMP2' ) reac_par(iTEMP2)%n_reac = reac_par(iTEMP2)%n_reac + 1
              IF ( TypeR == 'TEMP3' ) reac_par(iTEMP3)%n_reac = reac_par(iTEMP3)%n_reac + 1
            CASE ('CONST')
              reac_par(iCONST)%n_reac = reac_par(iCONST)%n_reac + 1
            CASE ('SPECIAL')
              reac_par(iSPECIAL)%n_reac = reac_par(iSPECIAL)%n_reac + 1 
            CASE DEFAULT
              WRITE(*,*) '  Unknown phase transfer reaction: ', TypeR
          END SELECT

        CASE ('AQUA')         ! aquatic phase reactions 

          nr_aqua = nr_aqua + 1
          CALL InsertReaction( ListRAqua , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('PHOTAB','PHOTABC','PHOTMCM')
              IF ( TypeR == 'PHOTAB'  ) reac_par(iPHOTAB)%n_reac  = reac_par(iPHOTAB)%n_reac + 1
              IF ( TypeR == 'PHOTABC' ) reac_par(iPHOTABC)%n_reac = reac_par(iPHOTABC)%n_reac + 1
              IF ( TypeR == 'PHOTMCM' ) reac_par(iPHOTMCM)%n_reac = reac_par(iPHOTMCM)%n_reac + 1
            CASE ('CONST')
              reac_par(iCONST)%n_reac = reac_par(iCONST)%n_reac + 1
            CASE ('TEMP','Temp1''TEMP2','TEMP3','TEMP4')
              IF ( TypeR == 'TEMP' )  reac_par(iTEMP)%n_reac  = reac_par(iTEMP)%n_reac + 1
              IF ( TypeR == 'TEMP1' ) reac_par(iTEMP1)%n_reac = reac_par(iTEMP1)%n_reac + 1
              IF ( TypeR == 'TEMP2' ) reac_par(iTEMP2)%n_reac = reac_par(iTEMP2)%n_reac + 1
              IF ( TypeR == 'TEMP3' ) reac_par(iTEMP3)%n_reac = reac_par(iTEMP3)%n_reac + 1
              IF ( TypeR == 'TEMP4' ) reac_par(iTEMP4)%n_reac = reac_par(iTEMP4)%n_reac + 1
            CASE ('ASPEC1','ASPEC2','ASPEC3','ASPEC4')
              IF ( TypeR == 'ASPEC1' ) reac_par(iASPEC1)%n_reac = reac_par(iASPEC1)%n_reac + 1
              IF ( TypeR == 'ASPEC2' ) reac_par(iASPEC2)%n_reac = reac_par(iASPEC2)%n_reac + 1
              IF ( TypeR == 'ASPEC3' ) reac_par(iASPEC3)%n_reac = reac_par(iASPEC3)%n_reac + 1
              !IF ( TypeR == 'ASPEC4' ) reac_par(iASPEC4)%n_reac = reac_par(iASPEC4)%n_reac + 1
            CASE ('SPECIAL')
              reac_par(iSPECIAL)%n_reac = reac_par(iSPECIAL)%n_reac + 1 
            CASE DEFAULT
              WRITE(*,*) '  Unknown aqueous reaction: ', TypeR
          END SELECT

        CASE ('DISS')        ! fast aquatic phase equil. reactions 

          nr_diss = nr_diss + 1
          CALL InsertReaction( ListRDiss , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('DCONST','DTEMP','DTEMP2','DTEMP3','DTEMP4','DTEMP5','MESKHIDZE')
              IF ( TypeR == 'DCONST' ) reac_par(iDCONST)%n_reac = reac_par(iDCONST)%n_reac + 1
              IF ( TypeR == 'DTEMP'  ) reac_par(iDTEMP)%n_reac  = reac_par(iDTEMP)%n_reac + 1
              IF ( TypeR == 'DTEMP2' ) reac_par(iDTEMP2)%n_reac = reac_par(iDTEMP2)%n_reac + 1
              IF ( TypeR == 'DTEMP3' ) reac_par(iDTEMP3)%n_reac = reac_par(iDTEMP3)%n_reac + 1
              IF ( TypeR == 'DTEMP4' ) reac_par(iDTEMP4)%n_reac = reac_par(iDTEMP4)%n_reac + 1
              IF ( TypeR == 'DTEMP5' ) reac_par(iDTEMP5)%n_reac = reac_par(iDTEMP5)%n_reac + 1
              IF ( TypeR == 'MESKHIDZE' ) reac_par(iMESKHIDZE)%n_reac  = reac_par(iMESKHIDZE)%n_reac + 1
            CASE ('SPECIAL')
              reac_par(iSPECIAL)%n_reac = reac_par(iSPECIAL)%n_reac + 1 
            CASE DEFAULT
              WRITE(*,*) '  Unknown dissociation reaction: ', TypeR
          END SELECT

        CASE ('SOLID')

          nr_solid = nr_solid + 1
          CALL InsertReaction( ListRSolid , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('EQUI')
              reac_par(iEQUI)%n_reac = reac_par(iEQUI)%n_reac + 1
            CASE ('DTEMP3')
              reac_par(iDTEMP3)%n_reac = reac_par(iDTEMP3)%n_reac + 1
            CASE ('SPECIAL')
              reac_par(iSPECIAL)%n_reac = reac_par(iSPECIAL)%n_reac + 1 
            CASE DEFAULT
              WRITE(*,*) '  Unknown solid reaction: ', TypeR
          END SELECT

        CASE ('PARTI')

          nr_parti = nr_parti + 1
          CALL InsertReaction( ListRPartic , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('SPECIAL')
              reac_par(iSPECIAL)%n_reac = reac_par(iSPECIAL)%n_reac + 1 
            CASE DEFAULT
              WRITE(*,*) '  Unknown particle reaction: ', TypeR
          END SELECT

        CASE ('MICROPHYS')

          nr_micphys = nr_micphys + 1
          CALL InsertReaction( ListRMicro , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('SPECIAL')
              reac_par(iSPECIAL)%n_reac = reac_par(iSPECIAL)%n_reac + 1 
            CASE DEFAULT
              WRITE(*,*) '  Unknown microphysic reaction: ', TypeR
          END SELECT

        CASE DEFAULT
          WRITE(*,*) '  Unknown reaction CLASS: ', CLASS
          STOP
      END SELECT

      Out = .FALSE.
    ELSE
      Out = .TRUE.
    END IF

  END SUBROUTINE ReadReaction_neu




  !
  SUBROUTINE ReadReaction(Out)
    LOGICAL :: Out
    !
    INTEGER :: iLine,PosColon,Pos,is
    CHARACTER(LenLine) :: LocLine
    CHARACTER(LenLine) :: Line(1:4)
    CHARACTER(20) :: CLASS
    CHARACTER(40) :: TypeR
    INTEGER :: idxFAC
    
    iLine = 0
    DO

      READ( InputUnit , '(A400)' , IOSTAT=is ) LocLine
      idxFAC = INDEX(LocLine,'$')

      IF ( idxFAC > 0 ) THEN
        SELECT CASE (TRIM(LocLine(idxFAC:)))
          CASE ('$H2','$O2N2','$M','$O2','$N2','$H2O','$O2O2','$aH2O','$+M','$(+M)','$RO2','$RO2aq')
            IF ( TRIM(LocLine(idxFAC:)) == '$H2'    ) nr_FAC_H2    = nr_FAC_H2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2N2'  ) nr_FAC_O2N2  = nr_FAC_O2N2  + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$M'     ) nr_FAC_M     = nr_FAC_M     + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2'    ) nr_FAC_O2    = nr_FAC_O2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$N2'    ) nr_FAC_N2    = nr_FAC_N2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$H2O'   ) nr_FAC_H2O   = nr_FAC_H2O   + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2O2'  ) nr_FAC_O2O2  = nr_FAC_O2O2  + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$RO2'   ) nr_FAC_RO2   = nr_FAC_RO2   + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$RO2aq' ) nr_FAC_RO2aq = nr_FAC_RO2aq + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$aH2O'  ) nr_FAC_aH2O  = nr_FAC_aH2O  + 2
            nr_FACTOR = nr_FACTOR + 1
          CASE DEFAULT
            WRITE(*,*) '  Unknown FACTOR:  ', TRIM(LocLine(idxFAC:)), '   at Line:  ???'
        END SELECT
      END IF

      IF ( ABS(is) > 0 ) EXIT
      
      ! if no comment or blank line then
      IF ( ADJUSTL(LocLine(1:1)) /= '#'       .AND.  &
      &    ADJUSTL(LocLine(1:7)) /= 'COMMENT' .AND.  &
      &    LEN(TRIM(LocLine)) > 0 ) THEN
        iLine = iLine + 1
        Line(iLine) = LocLine
        IF ( iLine == 4 ) EXIT
      END IF

    END DO

    IF ( iLine >= 3 ) THEN
      Pos = SCAN(Line(1),'#')
      IF ( Pos > 0 ) Line(1) = Line(1)(:Pos-1)
      
      ! if new reaction line starts go one line back
      IF ( INDEX(Line(4),'CLASS') > 0 ) BACKSPACE(InputUnit)
      
      ! read the reaction TYPE
      PosColon = Index(Line(1),':')
      CLASS    = ADJUSTL(Line(1)(PosColon+1:))
      
      ! count the number of each reaction type
      SELECT CASE (CLASS)
        
        CASE ('GAS')      ! gaseous phase reactions

          nr_gas = nr_gas + 1
          CALL InsertReaction( ListRGas , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('PHOTO','PHOTO2','PHOTO3','PHOTAB','PHOTABC','PHOTMCM')
              nr_G_photo = nr_G_photo + 1
              IF ( TypeR == 'PHOTAB'  ) nr_PHOTab  = nr_PHOTab  + 1
              IF ( TypeR == 'PHOTABC' ) nr_PHOTabc = nr_PHOTabc + 1
              IF ( TypeR == 'PHOTMCM' ) nr_PHOTmcm = nr_PHOTmcm + 1
              IF ( TypeR == 'PHOTO'   ) nr_PHOTOkpp  = nr_PHOTOkpp  + 1
              IF ( TypeR == 'PHOTO2'  ) nr_PHOTO2kpp = nr_PHOTO2kpp + 1
              IF ( TypeR == 'PHOTO3'  ) nr_PHOTO3kpp = nr_PHOTO3kpp + 1
            CASE ('CONST')
              nr_G_const = nr_G_const + 1
              nr_CONST   = nr_CONST   + 1
            CASE ('TEMP','TEMP1','TEMP2','TEMP3','TEMP4')
              nr_G_temp  = nr_G_temp + 1
              IF ( TypeR == 'TEMP' )  nr_TEMP  = nr_TEMP  + 1
              IF ( TypeR == 'TEMP1' ) nr_TEMP1 = nr_TEMP1 + 1
              IF ( TypeR == 'TEMP2' ) nr_TEMP2 = nr_TEMP2 + 1
              IF ( TypeR == 'TEMP3' ) nr_TEMP3 = nr_TEMP3 + 1
              IF ( TypeR == 'TEMP4' ) nr_TEMP4 = nr_TEMP4 + 1
            CASE ('TROE','TROEF','TROEQ','TROEQF','TROEXP','TROEMCM')
              nr_G_troe = nr_G_troe + 1
              IF ( TypeR == 'TROE'    ) nr_TROE    = nr_TROE    + 1
              IF ( TypeR == 'TROEF'   ) nr_TROEf   = nr_TROEf   + 1
              IF ( TypeR == 'TROEQ'   ) nr_TROEq   = nr_TROEq   + 1
              IF ( TypeR == 'TROEQF'  ) nr_TROEqf  = nr_TROEqf  + 1
              IF ( TypeR == 'TROEXP'  ) nr_TROExp  = nr_TROExp  + 1
              IF ( TypeR == 'TROEMCM' ) nr_TROEmcm = nr_TROEmcm + 1
            CASE ('SPEC1','SPEC2','SPEC3','SPEC4','SPEC1MCM',  &
            &     'SPEC2MCM','SPEC3MCM','SPEC4MCM','SPEC5MCM', &
            &     'SPEC6MCM','SPEC7MCM','SPEC8MCM','SPEC9MCM'  )
              nr_G_spec = nr_G_spec + 1
              IF ( TypeR == 'SPEC1' ) nr_SPEC1 = nr_SPEC1 + 1
              IF ( TypeR == 'SPEC2' ) nr_SPEC2 = nr_SPEC2 + 1
              IF ( TypeR == 'SPEC3' ) nr_SPEC3 = nr_SPEC3 + 1
              IF ( TypeR == 'SPEC4' ) nr_SPEC4 = nr_SPEC4 + 1
              IF ( TypeR == 'SPEC1MCM' ) nr_SPEC1mcm = nr_SPEC1mcm + 1
              IF ( TypeR == 'SPEC2MCM' ) nr_SPEC2mcm = nr_SPEC2mcm + 1
              IF ( TypeR == 'SPEC3MCM' ) nr_SPEC3mcm = nr_SPEC3mcm + 1
              IF ( TypeR == 'SPEC4MCM' ) nr_SPEC4mcm = nr_SPEC4mcm + 1
              IF ( TypeR == 'SPEC5MCM' ) nr_SPEC5mcm = nr_SPEC5mcm + 1
              IF ( TypeR == 'SPEC6MCM' ) nr_SPEC6mcm = nr_SPEC6mcm + 1
              IF ( TypeR == 'SPEC7MCM' ) nr_SPEC7mcm = nr_SPEC7mcm + 1
              IF ( TypeR == 'SPEC8MCM' ) nr_SPEC8mcm = nr_SPEC8mcm + 1
              IF ( TypeR == 'SPEC9MCM' ) nr_SPEC9mcm = nr_SPEC9mcm + 1
            CASE ('S4H2O')
              nr_S4H2O = nr_S4H2O + 1
            CASE ('T1H2O')
              nr_T1H2O = nr_T1H2O + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_G_special = nr_G_special + 1
            CASE ('HOM1')
              nr_HOM1 = nr_HOM1 + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown gaseous reaction: ', TypeR
          END SELECT

        CASE ('HENRY')        ! phase transfer pseudo-reactions 

          nr_henry = nr_henry + 1
          CALL InsertReaction( ListRHenry , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('TEMP','TEMP1','TEMP2','TEMP3')
              IF ( TypeR == 'TEMP' )  nr_TEMP  = nr_TEMP  + 1
              IF ( TypeR == 'TEMP1' ) nr_TEMP1 = nr_TEMP1 + 1
              IF ( TypeR == 'TEMP2' ) nr_TEMP2 = nr_TEMP2 + 1
              IF ( TypeR == 'TEMP3' ) nr_TEMP3 = nr_TEMP3 + 1
            CASE ('CONST')
              nr_CONST = nr_CONST + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_H_special = nr_H_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown phase transfer reaction: ', TypeR
          END SELECT

        CASE ('AQUA')         ! aquatic phase reactions 

          nr_aqua = nr_aqua + 1
          CALL InsertReaction( ListRAqua , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('PHOTAB','PHOTABC','PHOTMCM')
              nr_A_photo = nr_A_photo + 1
              IF ( TypeR == 'PHOTAB'  ) nr_PHOTab  = nr_PHOTab  + 1
              IF ( TypeR == 'PHOTABC' ) nr_PHOTabc = nr_PHOTabc + 1
              IF ( TypeR == 'PHOTMCM' ) nr_PHOTmcm = nr_PHOTmcm + 1
            CASE ('CONST')
              nr_A_const = nr_A_const + 1
              nr_CONST   = nr_CONST   + 1
            CASE ('TEMP','Temp1''TEMP2','TEMP3','TEMP4')
              nr_A_temp  = nr_A_temp + 1
              IF ( TypeR == 'TEMP' )  nr_TEMP  = nr_TEMP  + 1
              IF ( TypeR == 'TEMP1' ) nr_TEMP1 = nr_TEMP1 + 1
              IF ( TypeR == 'TEMP2' ) nr_TEMP2 = nr_TEMP2 + 1
              IF ( TypeR == 'TEMP3' ) nr_TEMP3 = nr_TEMP3 + 1
              IF ( TypeR == 'TEMP4' ) nr_TEMP4 = nr_TEMP4 + 1
            CASE ('ASPEC1','ASPEC2','ASPEC3','ASPEC4')
              nr_A_spec  = nr_A_spec + 1
              IF ( TypeR == 'ASPEC1' ) nr_ASPEC1 = nr_ASPEC1 + 1
              IF ( TypeR == 'ASPEC2' ) nr_ASPEC2 = nr_ASPEC2 + 1
              IF ( TypeR == 'ASPEC3' ) nr_ASPEC3 = nr_ASPEC3 + 1
              IF ( TypeR == 'ASPEC4' ) nr_ASPEC4 = nr_ASPEC4 + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_A_special = nr_A_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown aqueous reaction: ', TypeR
          END SELECT

        CASE ('DISS')        ! fast aquatic phase equil. reactions 

          nr_diss = nr_diss + 1
          CALL InsertReaction( ListRDiss , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('DCONST','DTEMP','DTEMP2','DTEMP3','DTEMP4','DTEMP5','MESKHIDZE')
              IF ( TypeR == 'DCONST'    ) nr_DCONST = nr_DCONST + 1
              IF ( TypeR == 'DTEMP'     ) nr_DTEMP  = nr_DTEMP  + 1
              IF ( TypeR == 'DTEMP2'    ) nr_DTEMP2 = nr_DTEMP2 + 1
              IF ( TypeR == 'DTEMP3'    ) nr_DTEMP3 = nr_DTEMP3 + 1
              IF ( TypeR == 'DTEMP4'    ) nr_DTEMP4 = nr_DTEMP4 + 1
              IF ( TypeR == 'DTEMP5'    ) nr_DTEMP5 = nr_DTEMP5 + 1
              IF ( TypeR == 'MESKHIDZE' ) nr_Meskhidze = nr_Meskhidze + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_D_special = nr_D_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown dissociation reaction: ', TypeR
          END SELECT

        CASE ('SOLID')

          nr_solid = nr_solid + 1
          CALL InsertReaction( ListRSolid , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('EQUI')
              nr_S_equi = nr_S_equi + 1
            CASE ('DTEMP3')
              nr_S_temp = nr_S_temp + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_S_special = nr_S_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown solid reaction: ', TypeR
          END SELECT

        CASE ('PARTI')

          nr_parti = nr_parti + 1
          CALL InsertReaction( ListRPartic , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_P_special = nr_P_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown particle reaction: ', TypeR
          END SELECT

        CASE ('MICROPHYS')

          nr_micphys = nr_micphys + 1
          CALL InsertReaction( ListRMicro , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_M_special = nr_M_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown microphysic reaction: ', TypeR
          END SELECT

        CASE DEFAULT
          WRITE(*,*) '  Unknown reaction CLASS: ', CLASS
          STOP
      END SELECT

      Out = .FALSE.
    ELSE
      Out = .TRUE.
    END IF

  END SUBROUTINE ReadReaction
  !
  !
  SUBROUTINE CompressParty(Ducts,Perm,Len)
    INTEGER, ALLOCATABLE :: Ducts(:)
    INTEGER, ALLOCATABLE :: Perm(:)
    INTEGER :: Len
    !
    ! sort ColInd and Val for acc column indx
    Perm=0
    Len=0
    CALL unirnk(Ducts,Perm,Len)
    Ducts=Ducts(Perm)
    CALL CompressList(Ducts)
  END SUBROUTINE CompressParty
  !
  !
  SUBROUTINE PrintSpecies(ListName,Unit,phs)
    TYPE(Species_T) :: ListName(:)
    INTEGER :: Unit
    INTEGER, ALLOCATABLE :: spclist(:)
    CHARACTER(1), OPTIONAL :: phs
    !
    INTEGER :: i

    DO i=1,SIZE(ListName)
      WRITE(Unit,*) "'"//TRIM(ListName(i)%Species)//"'"
    END DO
  END SUBROUTINE PrintSpecies
  !
  !
  SUBROUTINE SpcIdx(ListName,idx)
    TYPE(Species_T) :: ListName(:)
    INTEGER :: idx
    !
    INTEGER :: i
    !
    DO i=1,SIZE(ListName)
      IF (i==idx) THEN
        WRITE(*,*) "'"//TRIM(ListName(idx)%Species)//"'"
      END IF
    END DO
  END SUBROUTINE SpcIdx
  !
  !
  SUBROUTINE PrintHeadSpecies(Filename,Unit)
    INTEGER :: Unit
    CHARACTER(*) :: Filename
    !
    CHARACTER(8) :: Date
    CHARACTER(10) :: Time
    INTEGER(8) :: Value(8)
    !
    CALL DATE_AND_TIME(Date,Time,VALUES=Value)
    !
    WRITE(Unit,*) ' ==========================================================='
    WRITE(Unit,*) ' ========  0-dim Simulation of chemical mechanisms  ========'
    WRITE(Unit,*) ' ========     Output -  Chemical Reaction Data      ========'
    WRITE(Unit,*) ' ==========================================================='
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' Created:             ', Date(7:8),'.',Date(5:6),'.',Date(1:4)
    WRITE(Unit,*) ' Chemical Mechanism:  ', TRIM(ADJUSTL(FileName))
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' =================     Units         ======================='
    WRITE(Unit,*) ''
    IF (UnitGas==0) THEN
      WRITE(Unit,*) ' Gas Phase Units:     molec/cm3'
    ELSE
      WRITE(Unit,*) ' Gas Phase Units:     mol/m3'
    END IF
    IF (UnitAqua==0) THEN
      WRITE(Unit,*) ' Aqueous Phase Units: mol/l'
    END IF
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' =================    Numbers        ======================='
    WRITE(Unit,*) ''
    WRITE(Unit,*) ns_GAS &
                 +ns_AQUA &
                 +ns_PARTI &
                 +ns_KAT &
                 +ns_SOLID,  '     Number of Species'
    WRITE(Unit,*) ns_GAS,    '     No. of gaseous species'
    WRITE(Unit,*) ns_AQUA,   '     No. of aqueous species'
    WRITE(Unit,*) ns_PARTI, '     No. of particular species'
    WRITE(Unit,*) ns_SOLID,  '     No. of solid   species'
    WRITE(Unit,*) ns_KAT,'     Number of Non-reactive Species '
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' =================   Species Names   ======================='
    WRITE(Unit,*) ''
  END SUBROUTINE PrintHeadSpecies
  !
  !
  SUBROUTINE PrintFinalReactions(Unit)
    INTEGER :: Unit
    !
    WRITE(Unit,*) ''
    WRITE(Unit,*) ''
    WRITE(Unit,*) '========================================================='
    WRITE(Unit,*) '========              End  TAPE2                 ========'
    WRITE(Unit,*) '========     M3TRAS:  Chemical Reaction Data     ========'
    WRITE(Unit,*) '========================================================='
  END SUBROUTINE PrintFinalReactions
  !
  !
  SUBROUTINE PrintHeadReactions(Unit)
    INTEGER :: Unit
  
    nr    = nr_gas  + 2*nr_henry + nr_aqua   &
    &     + 2*nr_diss + nr_solid + nr_parti &
    &     + nr_micphys

    nr_D_Temp = nr_DTEMP  + nr_DTEMP2 &
    &         + nr_DTEMP3 + nr_DTEMP5
    
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' ================   Description of Reactions   =============='
    WRITE(Unit,*) ''
    WRITE(Unit,*) nr,          '        NREAK  : Number of Reactions'
    WRITE(Unit,*) nr_gas,      '        NGAS   : Gas phase reactions'
    WRITE(Unit,*) nr_G_photo,  '           Gaseous PHOTO - type reactions'
    WRITE(Unit,*) nr_G_const,  '           Gaseous CONST - type reactions'
    WRITE(Unit,*) nr_G_temp,   '           Gaseous TEMP - type reactions'
    WRITE(Unit,*) nr_SimpTB,   '           Gaseous Simple three-body - type reactions'
    WRITE(Unit,*) nr_G_lind,   '           Gaseous Lindemann - type reactions'
    WRITE(Unit,*) nr_G_troe,   '           Gaseous TROE - type reactions'
    WRITE(Unit,*) nr_G_spec,   '           Gaseous SPEC - type reactions'
    WRITE(Unit,*) nr_G_special,'           Gaseous SPECIAL formula - type reactions'
    WRITE(Unit,*) nr_henry,    '        NHENRY : Henry Equilib. reactions'
    WRITE(Unit,*) nr_diss,     '        NDISS  : Dissociation reactions'
    WRITE(Unit,*) nr_DCONST,   '           Aqueous DCONST - type reactions'
    WRITE(Unit,*) nr_D_TEMP,   '           Aqueous DTEMP - type reactions'
    WRITE(Unit,*) nr_D_special,'           Aqueous SPECIAL formula reactions'
    WRITE(Unit,*) nr_aqua,     '        NAQUA  : Aquatic Equilib. reactions'
    WRITE(Unit,*) nr_A_photo,  '           Aqueous PHOTO - type reactions'
    WRITE(Unit,*) nr_A_const,  '           Aqueous CONST - type reactions'
    WRITE(Unit,*) nr_A_temp,   '           Aqueous TEMP - type reactions'
    WRITE(Unit,*) nr_A_spec,   '           Aqueous SPEC - type reactions'
    WRITE(Unit,*) nr_A_special,'           Aqueous SPECIAL formula - type reactions'
    WRITE(Unit,*) nr_parti,   '        NPARTI : Particulare reactions   '
    WRITE(Unit,*) nr_P_special,'           Partic SPECIAL formula - type reactions'
    WRITE(Unit,*) nr_solid,    '        NSOLID : Solid Equilib. reactions'
    WRITE(Unit,*) nr_S_temp,   '           Solid DTEMP3 - type reactions'
    WRITE(Unit,*) nr_S_equi,   '           Solid EQUI - type reactions'
    WRITE(Unit,*) nr_S_spec,   '           Solid SPEC - type reactions'
    WRITE(Unit,*) nr_S_special,'           Solid SPECIAL formula - type reactions'
    WRITE(Unit,*) nr_micphys,  '        NMICRO : Microphysical reactions'
    WRITE(Unit,*) nr_M_special,'           Micro SPECIAL formula - type reactions'
    WRITE(Unit,*)
    WRITE(Unit,*) ' ======================  Reactions   ========================'
    WRITE(Unit,*) ''
  END SUBROUTINE PrintHeadReactions
  !
  SUBROUTINE Print_SysFile(RS,IndexSet,NewName)
    TYPE(ReactionStruct_T), INTENT(IN) :: RS(:)
    INTEGER,      OPTIONAL, INTENT(IN) :: IndexSet(:)
    CHARACTER(*), OPTIONAL             :: NewName
    
    INTEGER :: iR, j
    INTEGER :: nR
    LOGICAL :: done = .FALSE.
    
    OPEN(UNIT=989,FILE=ADJUSTL(TRIM(NewName)),STATUS='UNKNOWN')
  
    WRITE(989,'(A)') '# ================= '//TRIM(NewName)//' ================='
    WRITE(989,'(A)') '# = Please copy the data into your sys-file for ='
    WRITE(989,'(A)') '# =============== chemical input. ==============='
    WRITE(989,'(A)') '#'
    WRITE(989,'(A)') '#  ===================   Unit options   ======================'
    WRITE(989,'(A)') ''
    WRITE(989,'(A)') 'UNIT GAS    0   #    Gas phase units     (0 = molec/cm3, 1 = mol/m3)'
    WRITE(989,'(A)') 'UNIT AQUA   0   #    Aqueous phase units (0 = mol/l)'
    WRITE(989,'(A)') ''
    WRITE(989,'(A)') '#'

 
    IF (PRESENT(IndexSet)) THEN
      nR = SIZE(IndexSet)
      j = 0

      ! Print IndexSet of reactions
      DO 
        j = j + 1
        iR = IndexSet(j)

        IF ( .NOT. RS(iR)%bR ) THEN
          ! if direction is --> then
          CALL PrintReac(iR)
        ELSE
          ! for equilibrium reactions we need to check if the correct one is printet
          ! by correct one i mean the reaction which occurs in the original system
          ! not the reversed reaction string, becaus the evaluation of the rate 
          ! constant would be garbage if the reversed reaction string is printet
          IF ( iR-1 /= IndexSet(j-1) ) THEN
            CALL PrintReac(iR-1)
          END IF
        END IF

        IF ( j >= nR ) EXIT
      END DO
    ELSE
      nR = SIZE(RS)

      ! Print all reaction in RS list
      iR = 0
      DO 
        iR = iR + 1
        CALL PrintReac(iR)
        IF ( iR >= nR ) EXIT
        IF ( MAXVAL(INDEX(RS(iR)%Type,['DISS','HENR'])) > 0 ) iR = iR + 1
      END DO
    
    END IF
    CLOSE(989)

      CONTAINS
        
      SUBROUTINE PrintReac(iR)
        INTEGER :: iR
        WRITE(989,'(A)') 'CLASS: '//TRIM(RS(iR)%Type)
        WRITE(989,'(A)') TRIM(RS(iR)%Line1)
        IF ( TRIM(RS(iR)%Line3) /= '' ) WRITE(989,'(A)') TRIM(RS(iR)%Line3)
        !IF ( TRIM(RS(i)%Special%Formula) /= '' ) WRITE(989,'(A,I0)') 'SPECIAL: '//TRIM(RS(i)%Special%Formula)//';  ',RS(i)%Special%nVariables
        IF ( TRIM(RS(iR)%Factor)  /= '' ) WRITE(989,'(A)') 'FACTOR: '//TRIM(RS(iR)%Factor)
        WRITE(989,'(A)') 
      END SUBROUTINE PrintReac
  END SUBROUTINE Print_SysFile

  FUNCTION RemoveSpaces(String) RESULT(StringOut)
    ! replaces multiple spaces in string by one space
    CHARACTER(*), INTENT(IN) :: String
    CHARACTER(LEN(String))   :: StringOUT

    INTEGER :: i

    StringOut = TRIM(String)

    DO
      i = INDEX(TRIM(StringOut),'  ')
      IF ( i == 0 ) EXIT
      StringOut = TRIM(ADJUSTL(StringOut(:i-1))) &
      &           //' '//                        &
      &           TRIM(ADJUSTL(StringOut(i+1:)))
    END DO

  END FUNCTION RemoveSpaces

  !
  SUBROUTINE Print_ChemFile(RS,File,Unit,CK)
    ! IN:
    TYPE(ReactionStruct_T), ALLOCATABLE :: RS(:)
    CHARACTER(*) :: File
    INTEGER      :: Unit
    LOGICAL      :: CK
    ! TEMP:
    INTEGER      :: io_stat
    INTEGER      :: i,j,m,iR
    INTEGER      :: nEduct,nProd
    TYPE(Duct_T) :: ActiveEduct(30)
    TYPE(Duct_T) :: ActiveProduct(30)
    !
    INTEGER      :: nnzA, nnzB
    !
    INTEGER,  ALLOCATABLE :: tmpIdx(:)
    REAL(dp), ALLOCATABLE :: tmpVal(:)
    INTEGER,  ALLOCATABLE :: permutation(:)
    INTEGER :: newLen

    
     
    !-----------------------------------------------------------------------
    ! --- Build the reaction system
    !-----------------------------------------------------------------------
    IF ( .NOT.CK ) THEN
      CALL AllListsToArray( RS            &
      &                   , ListRGas    , ListRHenry  &
      &                   , ListRAqua   , ListRDiss   &
      &                   , ListRSolid  , ListRPartic &
      &                   , ListRMicro  )
    END IF
    

    OPEN(unit=Unit, file=File, status='replace', action='write', access='sequential', iostat=io_stat)
    IF ( io_stat /= 0 ) WRITE(*,*) '  ERROR creating chem-file :: ',io_stat
    REWIND(ChemUnit)
    !----------------------------------------------------------------
    ! ---  build the coeficient matrices and write .chem
    CALL PrintHeadSpecies ( File , Unit )
    
    IF ( ns_GAS   > 0 ) CALL PrintSpecies( ListGas2     , Unit )
    IF ( ns_AQUA  > 0 ) CALL PrintSpecies( ListAqua2    , Unit )
    IF ( ns_SOLID > 0 ) CALL PrintSpecies( ListSolid2   , Unit )
    IF ( ns_PARTI > 0 ) CALL PrintSpecies( ListPartic2  , Unit )
    IF ( ns_KAT   > 0 ) CALL PrintSpecies( ListNonReac2 , Unit )
     
    CALL PrintHeadReactions( Unit )


    !-----------------------------------------------------------------------
    ! --- print reactions and build A, B and (B-A) structure
    !-----------------------------------------------------------------------
   
    ! set matrix dimensions
    A%m  = neq;   A%n  = nspc
    B%m  = neq;   B%n  = nspc
    BA%m = neq;   BA%n = nspc
    
    ! Standart alloc
    ALLOCATE(A%RowPtr(A%m+1),B%RowPtr(B%m+1),BA%RowPtr(BA%m+1))
    A%RowPtr  = 0;    A%RowPtr(1)  = 1
    B%RowPtr  = 0;    B%RowPtr(1)  = 1
    BA%RowPtr = 0;    BA%RowPtr(1) = 1
    
    DO iR=1,neq
      ! count activ educts in reaction iR
      nEduct = 0
      !print*, 'DEBUG::chemsys    sizeRSe,p= ',iR,SIZE(RS(iR)%Educt),SIZE(RS(iR)%Product)
      !print*, 'DEBUG::chemsys    reaktion = ',TRIM(RS(iR)%Line1)
      DO i=1,SIZE(RS(iR)%Educt)
        SELECT CASE(RS(iR)%Educt(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','Micro','GAS')
            nEduct = nEduct + 1
            ActiveEduct(nEduct) = RS(iR)%Educt(i)
            !print*, 'debug::chemssys   ActiveEduct(nEduct)=',ActiveEduct(nEduct)
        END SELECT
      END DO
      ! count activ products in reaction iR
      nProd = 0
      DO i=1,SIZE(RS(iR)%Product)
        SELECT CASE(RS(iR)%Product(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','Micro','GAS')
            nProd = nProd + 1
            ActiveProduct(nProd) = RS(iR)%Product(i)
            !print*, 'debug::chemssys   ActiveProduct(nProd)=',ActiveProduct(nProd)
        END SELECT
      END DO
   
      !iR = iR + 1
      WRITE(Unit,*)
      WRITE(Unit,'(A,I6,A)') '#-----------', iR ,'. Reaction ----------- '
    
      WRITE(Unit,*) TRIM(RS(iR)%Type)//'   '//TRIM(RS(iR)%TypeConstant)
     
      WRITE(Unit,*) SIZE(RS(iR)%Educt), SIZE(RS(iR)%Product), nEduct, nProd
      
      ! Educt Matrix A
      IF ( nEduct > 1 ) THEN
        ALLOCATE( tmpIdx(nEduct), tmpVal(nEduct) )
        tmpIdx = 0;  tmpVal = ZERO
        DO j=1,nEduct
          tmpIdx(j) = PositionSpeciesAll(ActiveEduct(j)%Species)
          tmpVal(j) = ActiveEduct(j)%Koeff
        END DO
        CALL CompressList(tmpIdx,tmpVal)
        A%RowPtr(iR+1) = A%RowPtr(iR) + SIZE(tmpIdx)
        DEALLOCATE( tmpIdx , tmpVal )
      ELSE
        A%RowPtr(iR+1) = A%RowPtr(iR) + nEduct
      END IF

      ! Product Matrix B
      IF (nProd>1) THEN
        ALLOCATE( tmpIdx(nProd), tmpVal(nProd) )
        tmpIdx = 0;  tmpVal = ZERO
        DO j=1,nProd
          tmpIdx(j) = PositionSpeciesAll(ActiveProduct(j)%Species)
          tmpVal(j) = ActiveProduct(j)%Koeff
        END DO
        CALL CompressList(tmpIdx,tmpVal)
        B%RowPtr(iR+1) = B%RowPtr(iR) + SIZE(tmpIdx)
        DEALLOCATE( tmpIdx , tmpVal )
      ELSE
        B%RowPtr(iR+1) = B%RowPtr(iR) + nProd
      END IF

      ! ----------------------------------------------------
      ! SpeziesIndx Edukt=> 1:#Edukt von Reaktion iR
      ! SpeziesIndx Produkt=> #Edukt+1:#Edukt+#Produkt von Reaktion iR
      ! #aktiver Stoffe der Reaktion
      WRITE(Unit,*) (PositionSpeciesAll(RS(iR)%Educt(i)%Species),  &
      &             i=1,SIZE(RS(iR)%Educt)),                       &
      &             (PositionSpeciesAll(RS(iR)%Product(i)%Species),&
      &             i=1,SIZE(RS(iR)%Product)),                     &
      &             nEduct+nProd
      ! 
      !----------------------------------------------------
      ! Tupel: (SpeziesIndex,-Koeffzien) für alle aktiven Edukte (links)
      ! Tupel: (SpeziesIndex,+Koeffzien) für alle aktiven Produkte (rechts)
      WRITE(Unit,'(*(7X,I5,3X,F6.3))', ADVANCE='NO')                     &
      &                   ( PositionSpeciesAll(ActiveEduct(i)%Species)   &
      &                  ,  -ActiveEduct(i)%Koeff,i=1,nEduct )   &
      &                  ,( PositionSpeciesAll(ActiveProduct(i)%Species) &
      &                  ,   ActiveProduct(i)%Koeff,i=1,nProd )
      WRITE(Unit,*)
      !
      IF (RS(iR)%TypeConstant=='SPECIAL') THEN
        WRITE(Unit,*) TRIM(RS(iR)%Line3)
      ELSE
        WRITE(Unit,*) SIZE(RS(iR)%Constants), RS(iR)%Constants
      END IF
      !
      ! #Reaktionskonstanten, Reaktionskonstanten 1:#

      IF (RS(iR)%Factor /= '') WRITE(Unit,*) 'FACTOR:  ',RS(iR)%Factor

      SELECT CASE (RS(iR)%Factor)
        CASE ('$RO2');   hasRO2   = .TRUE.
        CASE ('$RO2aq'); hasRO2aq = .TRUE.
      END SELECT
      
      IF (CK) WRITE(Unit,*) 'EXTRA1:  ',ADJUSTL(TRIM(RS(iR)%Line2))
      IF (CK) WRITE(Unit,*) 'EXTRA2:  ',ADJUSTL(TRIM(RS(iR)%Line3))
    END DO
   
    ! loop again to set ColInd and Val on A and B
    nnzA = 0
    nnzB = 0
  
    ALLOCATE( A%ColInd(A%RowPtr(A%m+1)-1) , A%Val(A%RowPtr(A%m+1)-1) &
    &       , B%ColInd(B%RowPtr(B%m+1)-1) , B%Val(B%RowPtr(B%m+1)-1) )
    A%ColInd = 0; A%Val = ZERO
    B%ColInd = 0; B%Val = ZERO
    
    ALLOCATE(sumBAT(A%m)); sumBAT = ZERO
    !
    DO iR = 1,neq
      nEduct = 0
      DO i=1,SIZE(RS(iR)%Educt)
        SELECT CASE(RS(iR)%Educt(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','Micro','GAS')
            nEduct = nEduct + 1
            ActiveEduct(nEduct) = RS(iR)%Educt(i)
        END SELECT
      END DO
      nProd = 0
      DO i=1,SIZE(RS(iR)%Product)
        SELECT CASE(RS(iR)%Product(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','Micro','GAS')
            nProd = nProd + 1
            ActiveProduct(nProd) = RS(iR)%Product(i)
        END SELECT
      END DO
      
      ! set ColInd and Val on A and B
      IF (nEduct>1) THEN
        ALLOCATE(tmpIdx(nEduct),tmpVal(nEduct),permutation(nEduct))
        tmpIdx = 0;  tmpVal = ZERO;   permutation = 0
        DO j=1,nEduct
          tmpIdx(j) = PositionSpeciesAll(ActiveEduct(j)%Species)
          tmpVal(j) = ActiveEduct(j)%Koeff
        END DO
        
        ! sort ColInd and Val for acc column indx
        CALL unirnk(tmpIdx,permutation,newLen)
        tmpIdx = tmpIdx(permutation); tmpVal = tmpVal(permutation)
        CALL CompressList(tmpIdx,tmpVal)
        !

        DO m=1,newLen
          nnzA = nnzA + 1
          A%ColInd(nnzA) = tmpIdx(m)
          A%Val(nnzA) = tmpVal(m)
          sumBAT(iR)  = sumBAT(iR) - tmpVal(m)
        END DO
        DEALLOCATE(tmpIdx,tmpVal,permutation)
      ELSE
        ! reactions with only one educt
        DO m=1,nEduct
          nnzA = nnzA + 1
          A%ColInd(nnzA) = PositionSpeciesAll(ActiveEduct(m)%Species)
          A%Val(nnzA) = ActiveEduct(m)%Koeff
          sumBAT(iR)  = sumBAT(iR) - ActiveEduct(m)%Koeff
        END DO
      END IF
      !
      IF (nProd>1) THEN
        ALLOCATE(tmpIdx(nProd),tmpVal(nProd),permutation(nProd))
        tmpIdx = 0;  tmpVal = ZERO;   permutation = 0
        DO j=1,nProd
          tmpIdx(j) = PositionSpeciesAll(ActiveProduct(j)%Species)
          tmpVal(j) = ActiveProduct(j)%Koeff
        END DO
        CALL unirnk(tmpIdx,permutation,newLen)
        tmpIdx = tmpIdx(permutation); tmpVal = tmpVal(permutation)
        CALL CompressList(tmpIdx,tmpVal)
        !
        DO m=1,newLen
          nnzB = nnzB+1
          B%ColInd(nnzB) = tmpIdx(m)
          B%Val(nnzB) = tmpVal(m)
          sumBAT(iR)  = sumBAT(iR) + tmpVal(m)
        END DO
        DEALLOCATE(tmpIdx,tmpVal,permutation)
      ELSE
        DO m=1,nProd
          nnzB = nnzB + 1
          B%ColInd(nnzB) = PositionSpeciesAll(ActiveProduct(m)%Species)
          B%Val(nnzB) = ActiveProduct(m)%Koeff
          sumBAT(iR)  = sumBAT(iR) + ActiveProduct(m)%Koeff
        END DO
      END IF
    END DO

    A%nnz = nnzA
    B%nnz = nnzB

    CALL PrintFinalReactions( Unit )
    CLOSE(ChemUnit)
  END SUBROUTINE Print_ChemFile

  SUBROUTINE Setup_SpeciesOrder(A)
    
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    INTEGER :: iR, j, jj
    INTEGER :: nnz, cnt
    INTEGER, ALLOCATABLE :: tmpFO1(:), tmpFO2(:)
    INTEGER, ALLOCATABLE :: tmpSO1(:), tmpSO2(:)
    INTEGER, ALLOCATABLE :: tmpHO1(:), tmpHO2(:)
    REAL(dp), ALLOCATABLE :: atmpHO(:)
    REAL(dp), PARAMETER :: big = -99999999999999.d0

    nnz = A%nnz

    ALLOCATE( tmpFO1(nnz), tmpFO2(nnz), &
            & tmpSO1(nnz), tmpSO2(nnz), &
            & tmpHO1(nnz), tmpHO2(nnz), &
            & atmpHO(nnz)               )

    tmpFO1  = 0;      tmpFO2  = 0
    tmpSO1  = 0;      tmpSO2  = 0
    tmpHO1  = 0;      tmpHO2  = 0
    atmpHO  = big

    cnt = 0
    DO iR = 1 , neq
      DO jj = A%RowPtr(iR) , A%RowPtr(iR+1)-1
        cnt = cnt + 1
        IF      (A%Val(jj) == ONE) THEN
          tmpFO1(cnt) = iR
          tmpFO2(cnt) = A%ColInd(jj)
        ELSE IF (A%Val(jj) == TWO) THEN
          tmpSO1(cnt) = iR
          tmpSO2(cnt) = A%ColInd(jj)
        ELSE
          tmpHO1(cnt) = iR
          tmpHO2(cnt) = A%ColInd(jj)
          atmpHO(cnt) = A%Val(jj)
        END IF
      END DO
    END DO

    CALL CompressIntegerArray(tmpFO1); CALL CompressIntegerArray(tmpFO2)
    CALL CompressIntegerArray(tmpSO1); CALL CompressIntegerArray(tmpSO2)
    CALL CompressIntegerArray(tmpHO1); CALL CompressIntegerArray(tmpHO2)
    CALL CompressDoubleArray(atmpHO)
    nFirst_order  = SIZE(tmpFO1)
    nSecond_order = SIZE(tmpSO1)
    nHigher_order = SIZE(tmpHO1)

    ALLOCATE(iFO(nFirst_order,2))
    ALLOCATE(iSO(nSecond_order,2))
    ALLOCATE(iHO(nHigher_order,2))
    ALLOCATE(aHO(nHigher_order))

    iFO(:,1) = tmpFO1; iFO(:,2) = tmpFO2
    iSO(:,1) = tmpSO1; iSO(:,2) = tmpSO2
    iHO(:,1) = tmpHO1; iHO(:,2) = tmpHO2
    aHO      = atmpHO

  END SUBROUTINE Setup_SpeciesOrder


  SUBROUTINE InputChemicalData(InitFileName,DataFileName,MeteoFileName)
    CHARACTER(*) :: InitFileName, DataFileName, MeteoFileName

    INTEGER :: iSpc, i, iFr, iPos, lb, ub, id
    INTEGER :: io_stat
    REAL(dp) :: pi43, LWC
    CHARACTER(60) :: string = ''

    INTEGER, ALLOCATABLE :: katPhase(:,:)
    !
    ! for pH Start
    REAL(dp) :: kappa

    id = MPI_id + 1
    !
    pi43=4.0d0/3.0d0*PI
    !
    ! this is for mass transfer (accom , diffus term)
    ALLOCATE( InitValAct(nspc) , InitValKat(ns_KAT) ,  y_emi(nspc) , y_depos(nspc) )
    ALLOCATE( henry_diff(nspc+ns_KAT), henry_accom(nspc+ns_KAT) )

    InitValAct = 1.0e-20_dp
    InitValKat = 1.0e-20_dp
    y_emi   = ZERO
    y_depos = ZERO
    henry_diff  = ZERO;   henry_diff(1:ns_GAS)  = 5.0e-6_dp
    henry_accom = ZERO;   henry_accom(1:ns_GAS) = 5.0e-5_dp
    !
    !
    !=========================================================================
    !===  Read and save species names in lists
    !=========================================================================
    !
    !-- Open .chem-file and skip the first 24 lines (head)
    OPEN(unit=ChemUnit, file=ChemFile, status='old', action='read', access='sequential', iostat=io_stat)
    IF ( io_stat /= 0 ) WRITE(*,*) '  ERROR opening chem-file :: ',io_stat
    REWIND(ChemUnit)
    
    DO i=1,24;  READ(ChemUnit,*);  END DO         ! skip the first 24 lines

    !---  set indices for pH and water dissociation 
    Hp_ind   = 0
    OHm_ind  = 0
    aH2O_ind = 0
    SO4mm_ind = 0
    HSO4m_ind = 0
    Temp_ind = nspc + 1
    
    !=========================================================================
    !===  Read and Split Species Names and Initital Values
    !=========================================================================
    
    ALLOCATE(y_name(nspc+ns_KAT))

    ! read gaseous phase species
    IF (ns_GAS>0) THEN
      bGs(1) = 1
      bGs(2) = ns_GAS
      iGs = [(i, i=bGs(1),bGs(2))]

      DO iSpc = bGs(1),bGs(2)
        READ(ChemUnit,*,IOSTAT=io_stat)  y_name(iSpc)
      END DO

      InitValAct(iGs) = 1.e-20_dp
      CALL Read_INI_file( InitFileName , InitValAct, InitValKat , 'GAS' , 'INITIAL' )

      N2O2    = N2 + O2        ! passive species N2+O2
      N2O2H2O = N2O2 + H2O     ! passive species N2+O2+H2Oc

      ! fractions of passive species
      mN2  = N2/N2O2H2O
      mO2  = O2/N2O2H2O
      mH2O = H2O/N2O2H2O

      !---  Read gase phase emission values
      CALL Read_INI_file( InitFileName , y_emi , InitValKat , 'GAS' , 'EMISS' )
      CALL Read_INI_file( InitFileName , y_depos , InitValKat , 'GAS' , 'DEPOS' )

    END IF

    ! read aqueous phase species
    IF (ns_AQUA>0) THEN
      CALL Read_AFrac( AFrac, InitFileName )
      CALL Read_Modes( Mode , InitFileName )

      bAs(1) = ns_GAS + 1
      bAs(2) = ns_GAS + ns_AQUA
      iAs = [(i, i=bAs(1),bAs(2))]

      DO iSpc = bAs(1),bAs(2)
        READ(ChemUnit,*,IOSTAT=io_stat)  y_name(iSpc)
        IF ( y_name(iSpc) == 'Hp'        ) Hp_ind    = iSpc
        IF ( y_name(iSpc) == 'OHm'       ) OHm_ind   = iSpc
        IF ( y_name(iSpc) == 'SO4mm_ind' ) SO4mm_ind = iSpc
        IF ( y_name(iSpc) == 'HSO4m_ind' ) HSO4m_ind = iSpc
      END DO

      !=== Set  Chemical DATA  
      !===  (Molar Mass, Charges, Densities) 
      ALLOCATE( Charge(ns_AQUA)   &      ! charge of ions
      &       , SolubInd(ns_AQUA) &      ! solubility
      &       , MolMass(ns_AQUA)  &      ! molar mass of species
      &       , SpcDens(ns_AQUA)  &      ! density of species
      &       , OrgIndex(ns_AQUA) &      ! carbon atoms
      &       , CC(ns_AQUA)       &      ! compound class
      &       , ActIndex(ns_AQUA) )      ! index for calculation of activity coefficient

      !---  Set default values
      Charge   = ZERO
      SolubInd = ZERO
      MolMass  = ZERO
      SpcDens  = ONE
      OrgIndex = ZERO
      CC       = '  '
      ActIndex = ZERO
      
      !---  Determine Charges of Ions
      DO iSpc=1,ns_AQUA
        string = y_name(ns_GAS+iSpc)
        i = 1
        DO ! loop for cations
          iPos = INDEX(string(i:),'p')
          IF ( iPos <= 0 )  EXIT 
          Charge(iSpc) = Charge(iSpc) + 1
          i = i + iPos
        END DO
        i = 1
        DO ! loop for anions
          iPos = INDEX(string(i:),'m')
          IF ( iPos <= 0 )  EXIT 
          Charge(iSpc) = Charge(iSpc) - 1
          i = i + iPos
        END DO
      END DO

      InitValAct(iAs) = 1.e-16_dp
      CALL Read_INI_file( InitFileName , InitValAct, InitValKat , 'AQUA' , 'INITIAL' )


      DO i = 1 , SIZE(AFrac%Species)
        iPos = PositionSpeciesAll(AFrac%Species(i))
        IF ( iPos > 0 ) THEN
          InitValAct(iPos) = Mode%Number(id) * kilo       &  ! [#/m3]
          &                * AFrac%Frac1(i)               &  ! [g/g]
          &                * Mode%Radius(id)**3 * pi34    &  ! [m3]
          &                * Mode%Density(id)             &  ! [kg/m3]
          &                / AFrac%MolMass(i) * mol2Part     ! 1/[kg/mol] * [molec/mol]
        END IF
      END DO


      !---  Initial pH by charge balance 
      IF ( pHSet .AND. (Hp_ind*OHm_ind)>0 )  THEN
        Kappa = pHValue( InitValAct(bAs(1):bAs(2)) )
        IF ( Kappa > ZERO )  THEN
          InitValAct(Hp_ind) = Kappa
        ELSE 
          InitValAct(OHm_ind) = InitValAct(OHm_ind) + InitValAct(Hp_ind) - Kappa
        END IF
      END IF
    END IF  

    ! read solid phase species
    IF (ns_SOLID>0) THEN

      bSs(1) = ns_GAS + ns_AQUA + 1
      bSs(2) = ns_GAS + ns_AQUA + ns_SOLID
      iSs    = [(i, i=bSs(1),bSs(2))]
    
      DO iSpc = bSs(1),bSs(2)
        READ(ChemUnit,*,IOSTAT=io_stat)  y_name(iSpc) 
      END DO

      ! no solid input jet
    END IF

    ! read particle phase species
    IF (ns_PARTI>0) THEN

      bPs(1) = ns_GAS + ns_AQUA + ns_SOLID + 1
      bPs(2) = ns_GAS + ns_AQUA + ns_SOLID + ns_PARTI
      iPs    = [(i, i=bPs(1),bPs(2))]

      DO iSpc = bPs(1),bPs(2)
        READ(ChemUnit,*,IOSTAT=io_stat)  y_name(iSpc) 
      END DO
      
      ! no parti input jet
    END IF

    ! read catalytic phase species
    IF (ns_KAT>0) THEN

      lb = ns_GAS + ns_AQUA + ns_SOLID + ns_PARTI + 1
      ub = ns_GAS + ns_AQUA + ns_SOLID + ns_PARTI + ns_KAT
    
      DO iSpc = lb,ub
        READ(ChemUnit,*,IOSTAT=io_stat)  y_name(iSpc)
        henry_diff(iSpc)  = 5.0e-6_dp
        henry_accom(iSpc) = 5.0e-5_dp
        IF ( y_name(iSpc) == '[aH2O]' ) aH2O_ind = iSpc-(ns_AQUA+ns_GAS+ns_SOLID+ns_PARTI)
      END DO
      
      IF ( aH2O_ind > 0 ) THEN
        LWC  = pseudoLWC(tBegin)
        InitValKat(aH2O_ind) = InitValKat(aH2O_ind) * LWC * mol2part    ! convert aH2O [mol/L] to [molec/cm3]
      END IF
      
    END IF
  
    IF ( MPI_master .AND. ns_AQUA>0) THEN
      IF (hp_ind==0)   WRITE(*,*) '   ReadChem...Warning: Hp  not in mechanism!' 
      IF (ohm_ind==0)  WRITE(*,*) '   ReadChem...Warning: OHm  not in mechanism!' 
      IF (ah2o_ind==0) WRITE(*,*) '   ReadChem...Warning: aH2O  not in mechanism!' 
    END IF

    REWIND(ChemUnit);  CLOSE(ChemUnit)
    
    !
    !--- Read thermodynamic data,....
    CALL Read_SpeciesData( henry_diff, henry_accom, DataFileName )


    !
    !--- Boundaries for reaction phases
    bGr = [ 1                                    , nr_gas ]
    bHr = [ nr_gas+1                             , nr_gas+nr_henry ]
    bAr = [ nr_gas+nr_henry+1                    , nr_gas+nr_henry+nr_liquid ]
    bSr = [ nr_gas+nr_henry+nr_liquid+1          , nr_gas+nr_henry+nr_liquid+nr_solid ]
    bPr = [ nr_gas+nr_henry+nr_liquid+nr_solid+1 , nr_gas+nr_henry+nr_liquid+nr_solid+nr_parti ]
    
    iGr = [(i, i=bGr(1),bGr(2))]
    iHr = [(i, i=bHr(1),bHr(2))]
    iAr = [(i, i=bAr(1),bAr(2))]
    iSr = [(i, i=bSr(1),bSr(2))]
    iPr = [(i, i=bPr(1),bPr(2))]
  END SUBROUTINE InputChemicalData



  SUBROUTINE Read_SpeciesData(y_acc,y_diff,FileName)
    REAL(dp) :: y_acc(:) , y_diff(:) 
    CHARACTER(*) :: FileName
    !
    !
    CHARACTER(100) :: SpeciesName
    INTEGER :: iPos, i
    LOGICAL :: Back=.FALSE.
    REAL(dp) :: mm, alpha, dg, org, c1
    REAL(dp) :: nue
    CHARACTER(10) :: ro2d
    CHARACTER(10) :: c2
    INTEGER :: slash
    !INTEGER, ALLOCATABLE :: allRO2(:)
    !CHARACTER(100), ALLOCATABLE :: allRO2name(:)
    INTEGER :: ic1
   
    CALL OpenIniFile(FileName)
    !
    i=0
    !-----------------------------------------------------------
    ! --- Gas Phase thermodynamic data
    !
    GAS: DO 
      CALL LineFile( Back                         &
      &            , Start1 = 'BEGIN_DATAGAS'     &
      &            , End    = 'END_DATAGAS'       &
      &            , Name1  = SpeciesName         &
      &            , R1 = mm, R2 = alpha, R3 = dg )
      IF (Back)   EXIT
      !
      iPos = PositionSpeciesAll(SpeciesName)
      IF ( iPos > 0) THEN
        IF (alpha==ZERO .AND. dg==ZERO) CYCLE GAS 
        !
        y_acc(iPos)  = 1.0e-12_dp/(3.0_dp*dg)
        nue = SQRT(8.0e+03_dp*8.313_dp/Pi/mm)
        y_diff(iPos) = 4.0e-06_dp/(3.0_dp*alpha*nue)
        !
      END IF
    END DO GAS
    CALL RewindFile
    CALL ClearIniFile

    !
    !-----------------------------------------------------------
    ! --- Aqua Phase thermodynamic data
    !
    !!hier wird warhscneinlich das flasche abgespeichert
    IF ( ns_AQUA>0 ) THEN
      AQUA: DO 
        i = i + 1
        CALL LineFile( Back                     &
        &            , Start1 = 'BEGIN_DATAQUA' &
        &            , End    = 'END_DATAQUA'   &
        &            , Name1  = SpeciesName     &
        &            , R1 = mm, R2 = alpha      &
        &            , R3 = dg, R4 = org        &
        &            , Name2 = c2, rLen = 4     )
        IF (Back)  EXIT
        !
        iPos = PositionSpeciesAll(SpeciesName)
        IF ( iPos>0 ) THEN 
          IF (iPos<ns_GAS+ns_AQUA ) THEN
            iPos = iPos - ns_GAS
            MolMass(iPos)  = mm
            Charge(iPos)   = alpha
            SolubInd(iPos) = dg
            OrgIndex(iPos) = org
            CC(iPos) = ADJUSTL(TRIM(c2))
          ELSE
            ! do not save passive species infos?
          END IF
        END IF
      END DO AQUA
      CALL RewindFile
      CALL ClearIniFile
    END IF
    !
    !-----------------------------------------------------------
    ! --- Gas Phase RO2
    !
    IF ( hasRO2 ) THEN
      CALL LineFile( Back, Start1='BEGIN_DATARO2',    &
      &              End='END_DATARO2',               &
      &              Name1=SpeciesName,               &
      &              R1=c1)

      ALLOCATE(RO2(INT(c1)));  RO2 = 0 
      i = 0
      DO
        CALL LineFile( Back, Start1='BEGIN_DATARO2',  &
        &              End='END_DATARO2',             &
        &              Name1=SpeciesName)
     
        IF ( Back ) EXIT
        slash = INDEX(SpeciesName,'_')
        IF ( slash>0 ) SpeciesName(slash:slash)='/'
        
        IF (PositionSpeciesAll(SpeciesName) > 0) THEN
          i = i + 1
          RO2(i) = PositionSpeciesAll(SpeciesName)
        END IF
      END DO
      CALL RewindFile
      CALL ClearIniFile
      CALL CompressIntegerArray(RO2);   nRO2 = SIZE(RO2)
     
    END IF
    !
    !-----------------------------------------------------------
    ! --- Aqua Phase RO2
    !
    IF ( hasRO2aq ) THEN
      CALL LineFile( Back, Start1='BEGIN_DATARO2aq',  &
      &              End='END_DATARO2aq',             &
      &              Name1=SpeciesName,               &
      &              R1=c1)
      
      ALLOCATE(RO2aq(INT(c1)));      RO2aq = 0
      
      i=0
      DO
        CALL LineFile( Back, Start1='BEGIN_DATARO2aq',  &
        &              End='END_DATARO2aq',             &
        &              Name1=SpeciesName)
        IF (Back) EXIT
        IF (PositionSpeciesALL(SpeciesName)>0) THEN
          i = i + 1
          RO2aq(i)=PositionSpeciesAll(SpeciesName)
        END IF
      END DO
      CALL CompressIntegerArray(RO2aq);   nRO2aq = SIZE(RO2aq)

    END IF
    CALL CloseIniFile
      
    !  WRITE(333,*) ' nRO2=',SIZE(RO2)
    !DO i=1,SIZE(RO2)
    !  WRITE(333,*) i, RO2(i)
    !END DO 
    !WRITE(333,*) ' nRO2aq=',SIZE(RO2aq)
    !DO i=1,SIZE(RO2aq)
    !  WRITE(333,*) i, RO2aq(i)
    !END DO 
  END SUBROUTINE Read_SpeciesData
  !


  SUBROUTINE Read_INI_file(FileName,Activ,Passiv,Phase,section)
    CHARACTER(*) :: FileName
    !
    REAL(dp)      :: Activ(:)
    REAL(dp)      :: Passiv(:)
    CHARACTER(*)  :: Phase, section
    !
    INTEGER :: iPos
    REAL(dp) :: c1
    CHARACTER(LenName) :: SpeciesName
    LOGICAL :: Back=.FALSE.

    ! read initial values
    CALL OpenIniFile(FileName)
    DO
      CALL LineFile( Back,                       &
      &              Start1 = 'BEGIN_'//Phase,   &
      &              Start2 = 'BEGIN_'//section,   &
      &              End    = 'END_'//section,     &
      &              Name1  = SpeciesName,       &
      &              R1     = c1                 )
      IF (Back)   EXIT
      
      iPos = PositionSpeciesAll(SpeciesName)
      !print*, ' debug :: '//Phase//' spc :: ', TRIM(SpeciesName), c1, back, iPos
      IF (iPos>0) THEN
        SpeciesName = ADJUSTL(SpeciesName)

        IF (SpeciesName(1:1)=='['.AND. LEN_TRIM(SpeciesName)<maxLENinActDuct) THEN
          Passiv(iPos-nspc) = c1
          IF (Phase=='GAS') THEN
            IF (TRIM(SpeciesName)=='[H2O]') H2O = c1
            IF (TRIM(SpeciesName)=='[N2]')  N2  = c1
            IF (TRIM(SpeciesName)=='[O2]')  O2  = c1
            ns_G_KAT = ns_G_KAT + 1
          ELSE IF (Phase=='AQUA') THEN
            ns_A_KAT = ns_A_KAT + 1
          END IF
        ELSE
          !iPos=PositionSpeciesAll(SpeciesName)
          Activ(iPos) = c1
        END IF
    
      END IF
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_INI_file

  !
  SUBROUTINE Read_Diag(DiagSpc_Pos,DiagSpc_Phase,FileName)
    INTEGER, ALLOCATABLE :: DiagSpc_Pos(:)
    CHARACTER(1), ALLOCATABLE :: DiagSpc_Phase(:)
    CHARACTER(*) :: FileName
    !
    INTEGER :: iPos, cnt, lb, ub
    CHARACTER(50) :: SpeciesName
    LOGICAL :: Back
    !
    ! Read initial values of aqua spc
    CALL OpenIniFile(FileName)
    cnt=0
    DO
      CALL LineFile( Back,                &
      &              Start1='BEGIN_DIAG', &
      &              End   ='END_DIAG',   &
      &              Name1 =SpeciesName   )
      IF (Back)   EXIT
      !
      IF ( ADJUSTL(SpeciesName(1:1)) /= '#' .AND. &
        & PositionSpeciesAll(SpeciesName) > 0    ) cnt = cnt + 1
    END DO
    CALL CloseIniFile
    !
    ALLOCATE(DiagSpc_Pos(cnt),DiagSpc_Phase(cnt))
    ALLOCATE(iNcdfGas(0),iNcdfAqua(0),iNcdfSolid(0),iNcdfParti(0))
    DiagSpc_Pos   = 0
    DiagSpc_Phase = '-'
    !
    CALL OpenIniFile(FileName)
    cnt=0
    DO
      CALL LineFile( Back,                &
      &              Start1='BEGIN_DIAG', &
      &              End   ='END_DIAG',   &
      &              Name1 =SpeciesName   )
      IF (Back)   EXIT
      !
      IF (ADJUSTL(SpeciesName(1:1))/='#') THEN
        iPos = PositionSpeciesAll(SpeciesName)
        IF ( iPos > 0 ) THEN
          cnt = cnt + 1
          DiagSpc_Pos(cnt) = iPos
          IF      ( bGs(1) <= iPos .AND. iPos <= bGs(2) ) THEN
            DiagSpc_Phase(cnt)='g'
            iNcdfGas = [iNcdfGas,iPos]
            nNcdfGas = nNcdfGas + 1
          ELSE IF ( bAs(1) <= iPos .AND. iPos <= bAs(2) ) THEN
            DiagSpc_Phase(cnt)='a'
            iNcdfAqua = [iNcdfAqua,iPos]
            nNcdfAqua = nNcdfAqua + 1
          ELSE IF ( bSs(1) <= iPos .AND. iPos <= bSs(2) ) THEN
            DiagSpc_Phase(cnt)='s'
            iNcdfSolid = [iNcdfSolid,iPos]
            nNcdfSolid = nNcdfSolid + 1
          ELSE IF ( bPs(1) <= iPos .AND. iPos <= bPs(2) ) THEN
            DiagSpc_Phase(cnt)='p'
            iNcdfParti = [iNcdfParti,iPos]
            nNcdfParti = nNcdfParti + 1
          ELSE
            DiagSpc_Phase(cnt)='k'
          END IF
        END IF
      END IF
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_Diag


  SUBROUTINE Read_AFrac(AFrac,FileName) 
    CHARACTER(*)  :: FileName
    TYPE(AFrac_T) :: AFrac

    INTEGER       :: cnt
    REAL(dp)      :: c1,c2,c3,c4
    CHARACTER(20) :: SpeciesName
    LOGICAL       :: Back=.FALSE.
    
    ! Read AFRAC values
    CALL OpenIniFile(FileName)
    cnt=0
    DO
      CALL LineFile( Back,                      &
      &              Start1 ='BEGIN_AQUA',      &
      &              Start2 ='BEGIN_AFRAC',     &
      &              End    ='END_AFRAC',       &
      &              Name1  =SpeciesName,       &
      &              R1=c1, R2=c2, R3=c3, R4=c4 )
      IF (Back)   EXIT
      cnt=cnt+1
    END DO
    CALL CloseIniFile
    ALLOCATE(AFrac%Species(cnt),AFrac%MolMass(cnt) &
    &       ,AFrac%Charge(cnt), AFrac%SolubInd(cnt)&
    &       ,AFrac%Frac1(cnt))
    !
    CALL OpenIniFile(FileName)
    cnt=0
    DO
      CALL LineFile( Back,                      &
      &              Start1 ='BEGIN_AQUA',      &
      &              Start2 ='BEGIN_AFRAC',     &
      &              End    ='END_AFRAC',       &
      &              Name1  =SpeciesName,       &
      &              R1=c1, R2=c2, R3=c3, R4=c4 )
      IF (Back)   EXIT
      !
      cnt=cnt+1
      AFrac%Species(cnt) = TRIM(SpeciesName)
      AFrac%MolMass(cnt) = c1
      AFrac%Charge(cnt)  = INT(c2)
      AFrac%SolubInd(cnt)= c3
      AFrac%Frac1(cnt)   = c4
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_AFrac

  SUBROUTINE Read_Modes(Mode,FileName)
    CHARACTER(*) :: FileName
    TYPE(Modes_T) :: Mode
    !
    INTEGER :: cnt
    REAL(dp) :: c1,c2,c3
    LOGICAL :: Back
    REAL(dp) :: LWC
    !
    ! Read SPEK values
    CALL OpenIniFile(FileName)
    CALL LineFile( Back,                 &
    &              Start1 ='BEGIN_AQUA', &
    &              Start2 ='BEGIN_SPEK', &
    &              End    ='END_SPEK',   &
    &              R1=c1                 )
    !
    nFrac = INT(c1)
    ALLOCATE( Mode%Radius(nFrac)  , Mode%Number(nFrac)    &
    &       , Mode%Density(nFrac) , Mode%wetRadius(nFrac) ) 
    !
    LWC = pseudoLWC(tBegin)
    cnt = 0
    REWIND(InputUnit_Initials)
    CALL ClearIniFile()
    DO 
      CALL LineFile( Back,                 &
      &              Start1 ='BEGIN_AQUA', &
      &              Start2 ='BEGIN_SPEK', &
      &              End    ='END_SPEK',   &
      &              R1=c1, R2=c2, R3=c3   )
      IF (Back)   EXIT
      !
      ! HIER KÖNNTE MAL WIEDER WAS PASSIEREN
      IF (c1 < ONE) THEN
        cnt = cnt + 1
        Mode%Radius(cnt) = REAL(c1,KIND=dp)
        Mode%Number(cnt) = REAL(c2,KIND=dp)
        Mode%Density(cnt) = REAL(c3,KIND=dp)
        Mode%wetRadius(cnt) = (Pi34*LWC/Mode%Number(cnt))**(rTHREE)*0.1_dp
      END IF
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_Modes


  SUBROUTINE InputDatFile(FileName)
    CHARACTER(*) :: FileName
    !
    !
    INTEGER :: i
    REAL(dp) :: c1
    CHARACTER(80) :: SpeciesName
    CHARACTER(5) :: ro2d
    LOGICAL :: Back
    !
    CALL OpenIniFile(FileName)
    CALL LineFile( Back,                   &
    &              Start1='BEGIN_DATARO2', &
    &              End   ='END_DATARO2',   &
    &              Name1 =ro2d, R1=c1      )
    !
    i=0
    DO
      CALL LineFile( Back,                   &
      &              Start1='BEGIN_DATARO2', &
      &              End   ='END_DATARO2',   &
      &              Name1 =SpeciesName      )
      IF (Back) EXIT
      !
      IF (PositionSpecies(SpeciesName)>0) i=i+1
    END DO
    IF (i>0) THEN
      ALLOCATE(RO2spcG(i))
      RO2spcG=''
      ALLOCATE(RO2idxG(i))
      RO2idxG=0
      CALL RewindFile
      CALL ClearIniFile
      c1=0
      CALL LineFile( Back, Start1='BEGIN_DATARO2',  &
      &              End='END_DATARO2',             &
      &              Name1=ro2d,                    &
      &              R1=c1)
      !
      i=0
      DO
        CALL LineFile( Back, Start1='BEGIN_DATARO2',  &
        &              End='END_DATARO2',             &
        &              Name1=SpeciesName)
        IF (Back) EXIT
        IF (PositionSpecies(SpeciesName)>0) THEN
          i=i+1
          RO2spcG(i)=SpeciesName
          RO2idxG(i)=PositionSpecies(SpeciesName)
        END IF
      END DO
    END IF
    CLOSE(InputUnit) 
  END SUBROUTINE InputDatFile


  !
  SUBROUTINE InsertReaction(List,Line,TypeR)
    TYPE(ListReaction_T) :: List
    CHARACTER(*) :: Line(1:4)
    CHARACTER(*) :: TypeR
    !
    INTEGER :: PosColon,PosEqual,PosFactor,PosSpecial
    CHARACTER(LenLine) :: Left,Right
    TYPE(Reaction_T), POINTER :: Reaction
    CHARACTER(LenName), ALLOCATABLE :: Ducts(:)

    INTEGER :: i
    !
    IF (ASSOCIATED(List%Start)) THEN
      ALLOCATE(List%End%Next)
      List%End => List%End%Next
    ELSE
      ALLOCATE(List%End)
      List%Start => List%End
    END IF
    List%LenList = List%LenList + 1
    Reaction => List%End
    PosColon = Index(Line(1),':')
    Reaction%Type  = ADJUSTL(Line(1)(PosColon+1:))
    Reaction%Line1 = TRIM(Line(2))
    Reaction%Line3 = TRIM(Line(3))
    PosEqual = Index(Reaction%Line1,' = ')
    IF (PosEqual==0) THEN
      WRITE(*,*); WRITE(*,*)
      WRITE(*,'(10X,A,I0,A)') 'ERROR: Missing seperator " = " in reaction ',fNumber,':  '//TRIM(Line(2)) 
      WRITE(*,*); WRITE(*,*)
      STOP 
    ELSE
      fNumber = fNumber + 1
    END IF
    
    ! extract educts
    Left = Reaction%Line1(:PosEqual-1)
    CALL ExtractSpecies( Left, Reaction%Educt,     &
    &                    Reaction%InActEduct,      &
    &                    Reaction%nInActEd         )  ! geändert
    
    ! extract products
    Right = Reaction%Line1(PosEqual+3:)
    CALL ExtractSpecies( Right, Reaction%Product,  &
    &                    Reaction%InActProduct,    &
    &                    Reaction%nInActPro        )! geändert


    IF ( INDEX(Line(3),'SPECIAL') > 0 ) THEN
      Ducts = [Reaction%Educt(:)%Species , Reaction%Product(:)%Species]
      CALL ExtractConstants(Line(3),Ducts,Reaction%Constants,Reaction%TypeConstant,Reaction%Special)
    ELSE
      ! extract constants
      CALL ExtractConstants(Line(3),Ducts,Reaction%Constants,Reaction%TypeConstant)
    END IF
    Reaction%Line2  = Line(3)
    PosFactor = INDEX(Line(4),'FACTOR: ')

    IF (PosFactor > 0) THEN
      Reaction%Factor = TRIM(Line(4)(PosFactor+8:)) 
    ELSE
      Reaction%Factor = ''
    END IF
    !
    TypeR = ADJUSTL(Reaction%TypeConstant)
  END SUBROUTINE InsertReaction
  !
  !
  SUBROUTINE ReadUnits
    !
    INTEGER :: Pos
    CHARACTER(LenLine) :: LocLine
    INTEGER :: ios
    CHARACTER(400) :: iom
    LOGICAL :: comments
    !
    REWIND(InputUnit)
    DO 
      READ(InputUnit,'(A400)',IOSTAT=ios,IOMSG=iom) LocLine
      IF ( ios /= 0 ) THEN
        WRITE(*,*) ' Error Reading Units'
        WRITE(*,*) ' Error Message  ::  ',TRIM(iom)
      ELSE
        comments = (ADJUSTL(LocLine(1:1)) == '#'        .OR. &
        &           ADJUSTL(LocLine(1:7)) == 'COMMENT'  .OR. &
        &           LEN_TRIM(LocLine)     == 0 )
        IF ( comments ) CYCLE

        Pos = INDEX(LocLine,'GAS')
        IF ( Pos > 0 ) READ(LocLine(Pos+3:),*) UnitGas 
        Pos = INDEX(LocLine,'AQUA')
        IF ( Pos > 0 ) THEN
          READ(LocLine(Pos+4:),*) UnitAQUA 
          EXIT
        END IF
      END IF
    END DO

  END SUBROUTINE ReadUnits
  !
  !
  SUBROUTINE ExtractConstants(String,Ducts,Constants,Type,SpecialForm)
    CHARACTER(*) :: String
    REAL(dp), ALLOCATABLE :: Constants(:)
    !REAL(dp), POINTER :: Constants(:)
    CHARACTER(*) :: Type
    !
    INTEGER :: PosColon,PosName,PosComment,PosSemiColon
    INTEGER :: i,PosNum1,PosNum2,PosNum3,NumNum,PosElem
    CHARACTER(4)  :: NameNumNum
    CHARACTER(10) :: DummyString
    CHARACTER(LEN(String)) :: LocString
    CHARACTER(LEN(String)) :: LocString1
    CHARACTER(LEN(String)) :: NameConstant
    REAL(dp) :: Dummy
    INTEGER :: is

    ! this is for the new special formula input
    CHARACTER(LenLine)   :: SString
    TYPE(Special_T), OPTIONAL :: SpecialForm
    CHARACTER(LenName), ALLOCATABLE :: Ducts(:)
    INTEGER, ALLOCATABLE :: iSortedDucts(:)
    INTEGER :: nvs, cnt, idxDuct, lenDuct
    !
    LocString  = String
    String     = ''
    PosColon   = INDEX(LocString,':')
    Type       = LocString(:PosColon-1)
    LocString  = ADJUSTL(LocString(PosColon+1:))
    PosComment = INDEX(LocString,'#')

    IF ( PosComment > 0 ) LocString = LocString(:PosComment-1)
    
    LocString1 = LocString

    IF ( Type /= 'SPECIAL' ) THEN
      ALLOCATE(Constants(0))
      DO
        PosColon = INDEX(LocString1,':')
        IF ( PosColon > 0 ) THEN
          LocString1 = ADJUSTL(LocString1(PosColon+1:))
          READ( LocString1 , * , IOSTAT=is ) Dummy, DummyString
          PosName   = INDEX(LocString1,TRIM(DummyString))
          Constants = [Constants , Dummy]
          IF ( PosName > 0 ) LocString1 = LocString1(PosName:)
        ELSE
          EXIT
        END IF
      END DO

    ELSE
      ! special rate constatn formula
      IF ( PosColon==0 ) THEN
        WRITE(*,*) '  Missing seperator ":" for TypeConstant SPECIAL'
        STOP
      END IF

      PosSemiColon = Index(LocString,';')
      IF ( PosSemiColon==0 ) THEN
        WRITE(*,*) '  Missing seperator ";" for TypeConstant SPECIAL'
        STOP
      END IF

      ! save formula
      SString = ADJUSTL(LocString1(:PosSemiColon-1))

      ! read number of variables in formula
      READ( LocString1(PosSemiColon+1:) , * , IOSTAT=is ) nvs

      SpecialForm%nVariables = nvs      
      SpecialForm%Formula  = ADJUSTL(SString)
      IF (INDEX(SString,'TEMP')>0) SpecialForm%Temp = .TRUE.

      ALLOCATE(SpecialForm%cVariables(nvs))
      ALLOCATE(iSortedDucts(SIZE(Ducts)) )

      CALL StringSort(Ducts,iSortedDucts)

      cnt = 0
      DO i = SIZE(Ducts), 1 , -1
        idxDuct = INDEX(SString,TRIM(Ducts(i)))

        ! if a species if katalytic ( on both sides of the reaction equation )
        IF ( i > 1 ) THEN
          IF ( Ducts(i)==Ducts(i-1) ) CYCLE
        END IF

        IF ( idxDuct > 0 ) THEN
          cnt = cnt + 1
          SpecialForm%cVariables(cnt) = ADJUSTL(Ducts(i))

          DO
            idxDuct = INDEX(SString,TRIM(Ducts(i)))
            IF ( idxDuct == 0 ) EXIT
            lenDuct = LEN_TRIM(Ducts(i))
            SString(idxDuct:idxDuct+lenDuct) = REPEAT('_',lenDuct)
          END DO

        END IF
      END DO

      IF ( SpecialForm%Temp ) THEN
        SpecialForm%cVariables(cnt+1) = 'TEMP'
      ELSE
        WRITE(*,*) '      Warning: Missing temperature "TEMP" in SPECIAL Formula :: '&
        &           ,TRIM(SpecialForm%Formula)
      END IF

    END IF
  END SUBROUTINE ExtractConstants
  !
  !
  SUBROUTINE ExtractSpecies(String,Duct,InActDuct,NumInActDucts)
    ! IN:
    CHARACTER(*) :: String
    ! OUT:
    TYPE(Duct_T), POINTER :: Duct(:)
    TYPE(Duct_T), POINTER :: InActDuct(:)
    INTEGER :: NumInActDucts
    !
    INTEGER :: PosMinus,PosPlus,NumSpec,PosSpecies
    REAL(dp) :: PreFac !NumberSpecies
    CHARACTER(LenLine) :: Species
    CHARACTER(LEN(String)) :: LocString
    INTEGER :: sbL, sbR
    INTEGER :: ios
    LOGICAL :: dummy=.FALSE.

    !
    LocString = String
    NumSpec   = 1

    !WRITE(*,*) ' Current sys String :: ',TRIM(String)

    ! count species
    DO
     PosPlus  = INDEX(LocString,' + ')
     PosMinus = INDEX(LocString,' - ')

     IF ( PosPlus>0 .AND. (PosMinus==0 .OR. PosMinus>PosPlus) ) THEN
       LocString = LocString(PosPlus+3:)
       NumSpec   = NumSpec + 1
     END IF

     IF ( PosMinus>0 .AND. (PosPlus==0 .OR. PosPlus>PosMinus) ) THEN
       LocString = LocString(PosMinus+3:)
       NumSpec   = NumSpec + 1
     END IF

     IF ( PosPlus==0 .AND. PosMinus==0 ) EXIT
    END DO
   
    ALLOCATE( Duct(NumSpec) , InActDuct(NumSpec) )
    LocString = String
    NumSpec   = 0
    DO
      PosPlus  = INDEX(LocString,' + ')
      PosMinus = INDEX(LocString,' - ')
      IF ( PosMinus>0 .AND. PosMinus<PosPlus ) THEN
        PreFac = mONE
      ELSE
        PreFac = ONE
      END IF
      IF      ( PosPlus>0  .AND. (PosPlus<PosMinus.OR.PosMinus==0) ) THEN
        Species   = ADJUSTL(LocString(:PosPlus-1))
        LocString = LocString(PosPlus+3:)
      ELSE IF ( PosMinus>0 .AND. (PosMinus<PosPlus.OR.PosPlus==0) ) THEN
        Species   = ADJUSTL(LocString(:PosMinus-1))
        LocString = LocString(PosMinus+3:)
      ELSE
        Species = ADJUSTL(LocString)
      END IF

      PosSpecies = SCAN(Species,SetSpecies)
      NumSpec = NumSpec + 1


      ! check if there is a dummy species e.g. 0.000 with no species
      ! or species: (dummy)
      dummy = SCAN(TRIM(ADJUSTL(Species)) , SetSpecies) == 0 .OR. &
            & INDEX(TRIM(ADJUSTL(Species)) , '(dummy)') /= 0 

      !WRITE(*,*), ' species  ::: ', TRIM(Species)//'   is dummy = ', dummy

      IF (PosSpecies==1) THEN           
        sbL = INDEX(TRIM(Species),'[');  sbR = INDEX(TRIM(Species),']')

        ! check if species if passive
        IF ( sbL==1 .AND. LEN_TRIM(Species)==sbR-sbL+1 ) THEN
          ! works if there's just one InActEduct
          InActDuct(1)%Koeff   = PreFac
          InActDuct(1)%Species = Species
          NumInActDucts        = NumInActDucts + 1
        END IF
        Duct(NumSpec)%Koeff   = PreFac
        Duct(NumSpec)%Species = Species
      ELSE

        IF (.NOT.dummy) THEN
          READ( Species(1:PosSpecies-1) ,*, IOSTAT=ios) Duct(NumSpec)%Koeff
          IF ( ios == 0 ) THEN
            Duct(NumSpec)%Koeff   = PreFac * Duct(NumSpec)%Koeff
            Duct(NumSpec)%Species = Species(PosSpecies:)
          ELSE 
            WRITE(*,*) ' Error reading species and coefficients :: ', Species
            WRITE(*,*) ' IOSTAT = ', ios
          END IF
        ELSE
          Duct(NumSpec)%Koeff   = ZERO
          Duct(NumSpec)%Species = '(dummy)'
          !WRITE(*,*) ' ?? Dummy Species = ',TRIM(Species)
        END IF
      END IF

      ! Syntax check for missing seperators, e.g. "+CO2" insted of "+ CO2" or "CO2+" insted of "CO2 +"
      IF  ( INDEX(TRIM(ADJUSTL(Duct(NumSpec)%Species)),' ',.TRUE.) > 0 ) THEN
        WRITE(*,*); WRITE(*,*)
        WRITE(*,777) 'ERROR: Missing white space in reaction substring: '//TRIM(String)
        WRITE(*,777) '       Species = '//TRIM(Duct(NumSpec)%Species)
        WRITE(*,777) '       Check syntax in '//TRIM(SysFile)//'.sys!'
        WRITE(*,*); WRITE(*,*)
        CALL FinishMPI()
        STOP 
      END IF

      ! if no dummy species was found then add new species to hash table
      CALL InsertSpecies( Duct(NumSpec)%Species, Duct(NumSpec)%Type )

      ! if there are no further species exit the loop
      IF ( PosPlus==0 .AND. PosMinus==0 ) EXIT

    END DO
    777 FORMAT(10X,A)
  END SUBROUTINE ExtractSpecies
  
  !
  SUBROUTINE InsertSpecies(Species,Type)
    CHARACTER(*) :: Species
    CHARACTER(*) :: Type

    INTEGER :: len

    len = LEN_TRIM(Species)

    IF (Species(1:1)=='p') THEN
      CALL InsertHash( ListPartic , TRIM(ADJUSTL(Species)) , ns_PARTI)
      Type = 'Partic'
    ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
      CALL InsertHash( ListAqua   , TRIM(ADJUSTL(Species)) , ns_AQUA)
      Type = 'Aqua'
    ELSE IF (Species(1:1)=='s') THEN
      CALL InsertHash( ListSolid  , TRIM(ADJUSTL(Species)) , ns_SOLID)
      Type = 'Solid'
    ELSE IF ( Species(1:1)=='[' .AND. Species(len:len)==']' &
            &                   .AND. len<maxLENinActDuct   ) THEN
      CALL InsertHash( ListNonReac, TRIM(ADJUSTL(Species)) , ns_KAT)
      Type = 'Inert'
    ELSE IF ( MAXVAL(INDEX(Species(1:1) , ['(','0'])) > 0 ) THEN  ! dummy species
      Type = 'Dummy'
    ELSE
      CALL InsertHash( ListGas ,TRIM(ADJUSTL(Species)) ,     ns_GAS)
      Type = 'Gas'
    END IF
  END SUBROUTINE InsertSpecies
  !
  !
  FUNCTION PositionListSpecies(Species)
    TYPE(Species_T), POINTER :: PositionListSpecies
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpecies
    !
    PositionSpecies=0
    PositionListSpecies=>NULL()
    IF (Species(1:1)=='p') THEN
      PositionSpecies=GetHash(ListPartic,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListPartic2(PositionSpecies)
    ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
      PositionSpecies=GetHash(ListAqua,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListAqua2(PositionSpecies)
    ELSE IF (Species(1:1)=='s') THEN
      PositionSpecies=GetHash(ListSolid,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListSolid2(PositionSpecies)
    ELSE IF (Species(1:1)=='['.AND.Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
      PositionSpecies=GetHash(ListNonReac,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListNonReac2(PositionSpecies)
    ELSE
      PositionSpecies=GetHash(ListGas,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListGas2(PositionSpecies)
    END IF
  END FUNCTION PositionListSpecies
  !
  !
  FUNCTION PositionSpeciesGas(Species)
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpeciesGas
    ! 
    PositionSpeciesGas=0
    PositionSpeciesGas=GetHash(ListGas,TRIM(ADJUSTL(Species)))
  END FUNCTION
  !
  !
  FUNCTION PositionSpeciesCK(Species)
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpeciesCK
    ! 
    PositionSpeciesCK=0
    PositionSpeciesCK=GetHash(ListGas,TRIM(ADJUSTL(Species)))
  END FUNCTION PositionSpeciesCK
  !
  !
  !
  FUNCTION PositionSpecies(Species)
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpecies
    ! 
    PositionSpecies=0
    IF (Species(1:1)=='p') THEN
      PositionSpecies=GetHash(ListPartic,TRIM(ADJUSTL(Species)))     
    ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
      PositionSpecies=GetHash(ListAqua,TRIM(ADJUSTL(Species)))       &
      &               + ns_GAS
    ELSE IF (Species(1:1)=='s') THEN
      PositionSpecies=GetHash(ListSolid,TRIM(ADJUSTL(Species)))      
    ELSE IF (Species(1:1)=='['.AND.                                  &
    &        Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
      PositionSpecies=GetHash(ListNonReac,TRIM(ADJUSTL(Species)))  
    ELSE
      PositionSpecies=GetHash(ListGas,TRIM(ADJUSTL(Species)))
    END IF
  END FUNCTION PositionSpecies

  !
  FUNCTION PositionSpeciesAll(Species) RESULT(Pos)
    CHARACTER(*) :: Species
    !
    INTEGER :: Pos
    INTEGER :: len
   
    Pos = 0
    len = LEN_TRIM(Species)

    !WRITE(*,*) ' number ob species = ',species, len
  
    ! Combustion system
    IF ( Teq ) THEN
      Pos = -1
      Pos = GetHash(ListGas,TRIM(ADJUSTL(Species)))
    ELSE

    ! tropospheric system
      IF ( Species(1:1) == 'p' ) THEN
        Pos = GetHash(ListPartic,TRIM(ADJUSTL(Species))) 
        IF (Pos>0) Pos = Pos + ns_GAS + ns_AQUA + ns_SOLID         
      
      ! AQUA 
      ELSE IF ( Species(1:1)=='a'.OR.SCAN(Species,'pm')>0 ) THEN
        Pos = GetHash(ListAqua,TRIM(ADJUSTL(Species)))
        IF (Pos>0) Pos = Pos + ns_GAS
 
      ! SOLID
      ELSE IF ( Species(1:1)=='s' ) THEN
        Pos = GetHash(ListSolid,TRIM(ADJUSTL(Species)))      
        IF (Pos>0) Pos = Pos + ns_GAS + ns_AQUA         

      ! NonReac
      ELSE IF ( Species(1:1) == '[' .AND. Species(len:len) == ']' .AND. &
      &        len < maxLENinActDuct) THEN
        Pos = GetHash(ListNonReac,TRIM(ADJUSTL(Species)))    
        IF (Pos>0) Pos = Pos + ns_GAS + ns_AQUA + ns_SOLID + ns_PARTI
        
      ! GAS
      ELSE
        Pos = GetHash(ListGas,TRIM(ADJUSTL(Species)))
      END IF
    END IF
  END FUNCTION PositionSpeciesAll
  !
  FUNCTION PositionAtom(Atom) RESULT(Pos)
    CHARACTER(*) :: Atom
    !
    INTEGER :: Pos
    ! 
    Pos = 0
    Pos = GetHash(ListAtoms,TRIM(ADJUSTL(Atom)))
  END FUNCTION PositionAtom
  !
  !
  !
  SUBROUTINE OpenFile(FileName,Type)
    CHARACTER(*) :: Filename
    CHARACTER(*) :: Type
    !
    LOGICAL :: ExistFile
    !
    INQUIRE(FILE=TRIM(Filename)//'.'//TRIM(Type),EXIST=ExistFile)
    IF (ExistFile) THEN
      OPEN(UNIT=InputUnit,FILE=TRIM(Filename)//'.'//TRIM(Type),STATUS='UNKNOWN')
    END IF
  END SUBROUTINE OpenFile
  !
  !
  SUBROUTINE CloseFile(FileName,Type)
    CHARACTER(*) :: Filename
    CHARACTER(*) :: Type
    !
    LOGICAL :: ExistFile
    !
    INQUIRE(FILE=TRIM(Filename)//'.'//TRIM(Type),EXIST=ExistFile)
    IF (ExistFile) CLOSE(UNIT=InputUnit)
  END SUBROUTINE CloseFile
  !
  !
  SUBROUTINE ReadThermoData(FileName)
    CHARACTER(*) :: Filename
    !
    CHARACTER(LenLine) :: LocLine
    CHARACTER(LenName) :: Name
    !
    INTEGER :: is,iLine,i
    REAL(8) :: Hf,Gf,Cp
    TYPE(Species_T), POINTER :: Species
    TYPE(Reaction_T), POINTER :: Current
    !
    CALL OpenFile(FileName,'dat')
    iLine=0
    DO
      READ(InputUnit,'(a400)',IOSTAT=is) LocLine
      IF (ABS(is)>0) THEN
        EXIT
      END IF
      IF (ADJUSTL(LocLine(1:1))/='#'.AND.LEN(TRIM(LocLine))>0) THEN
        READ(LocLine,*) Name
        IF (PositionSpecies(Name)>0) THEN
          iLine=iLine+1
        END IF
      END IF
    END DO
    REWIND(InputUnit)
    DO
      READ(InputUnit,'(a400)',IOSTAT=is) LocLine
      IF (ABS(is)>0) THEN
        EXIT
      END IF
      IF (ADJUSTL(LocLine(1:1))/='#'.AND.LEN(TRIM(LocLine))>0) THEN
        READ(LocLine,*) Name,Hf,Gf,Cp
        Species=>PositionListSpecies(Name)
        IF (ASSOCIATED(Species)) THEN
          Species%Hf=Hf
          Species%Gf=Gf
          Species%Cp=Cp
        END IF
      END IF
    END DO
    CLOSE(InputUnit)
    !
    Current=>ListRSolid%Start
    DO
      IF (ASSOCIATED(Current)) THEN
        IF (Current%TypeConstant=='DTEMP3') THEN
          Hf=0.0d0
          Gf=0.0d0
          Cp=0.0d0
          DO i=1,SIZE(Current%Educt)
            Species=>PositionListSpecies(Current%Educt(i)%Species)
            Hf=Hf-Current%Educt(i)%Koeff*Species%Hf
            Gf=Gf-Current%Educt(i)%Koeff*Species%Gf
            Cp=Cp-Current%Educt(i)%Koeff*Species%Cp
          END DO
          DO i=1,SIZE(Current%Product)
            Species=>PositionListSpecies(Current%Product(i)%Species)
            Hf=Hf+Current%Product(i)%Koeff*Species%Hf
            Gf=Gf+Current%Product(i)%Koeff*Species%Gf
            Cp=Cp+Current%Product(i)%Koeff*Species%Cp
          END DO
          Hf=Hf*1000.0d0
          Gf=Gf*1000.0d0
          IF (ALLOCATED(Current%Constants)) THEN
            DEALLOCATE(Current%Constants)
          END IF
          ALLOCATE(Current%Constants(3))
          Current%Constants(1)=EXP(-Gf/(RGas*TRef))
          Current%Constants(2)=Cp/RGas
          Current%Constants(3)=-Hf/RGas+Cp*TRef/RGas
        END IF
        Current=>Current%Next
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE ReadThermoData
  !
  !
  SUBROUTINE ReadSystem(FileName)
    CHARACTER(*) :: Filename
    !
    LOGICAL :: Out

    FileName = FileName(:INDEX(FileName,'.')-1)
    !
    CALL InitHashTable(ListAqua,100)
    CALL InitHashTable(ListGas,100)
    CALL InitHashTable(ListSolid,100)
    CALL InitHashTable(ListPartic,100)
    CALL InitHashTable(ListNonReac,100)
    CALL OpenFile(FileName,'spc')
    DO
      CALL ReadSpecies(Out);  IF (Out) EXIT
    END DO
    CALL CloseFile(FileName,'spc')
    CALL OpenFile(FileName,'sys')
    CALL ReadUnits
    DO
      CALL ReadReaction(Out); IF (Out) EXIT
    END DO
    CALL CloseFile(FileName,'sys')
    ALLOCATE(ListGas2(ns_GAS))
    CALL HashTableToList(ListGas,ListGas2)
    CALL SortList(ListGas2)
    CALL ListToHashTable(ListGas2,ListGas)
    ALLOCATE(ListAqua2(ns_AQUA))
    CALL HashTableToList(ListAqua,ListAqua2)
    CALL SortList(ListAqua2)
    CALL ListToHashTable(ListAqua2,ListAqua)
    ALLOCATE(ListSolid2(ns_SOLID))
    CALL HashTableToList(ListSolid,ListSolid2)
    CALL SortList(ListSolid2)
    CALL ListToHashTable(ListSolid2,ListSolid)
    ALLOCATE(ListPartic2(ns_PARTI))
    CALL HashTableToList(ListPartic,ListPartic2)
    CALL SortList(ListPartic2)
    CALL ListToHashTable(ListPartic2,ListPartic)
    ALLOCATE(ListNonReac2(ns_KAT))
    CALL HashTableToList(ListNonReac,ListNonReac2)
    CALL SortList(ListNonReac2)
    CALL ListToHashTable(ListNonReac2,ListNonReac)

  END SUBROUTINE ReadSystem
  !
  !
  SUBROUTINE SortList(List)
    TYPE(Species_T) :: List(:)
    !
    TYPE(Species_T) :: Temp
    INTEGER :: i,j
    !
    DO i=1,SIZE(List)
      DO j=1,SIZE(List)-i
        IF (List(j+1)%Species<List(j)%Species) THEN
          Temp=List(j+1)
          List(j+1)=List(j)
          List(j)=Temp
        END IF
      END DO
    END DO
    DO i=1,SIZE(List)
      IF (List(i)%Species=='OHm') THEN
        Temp=List(i)
        IF (i==SIZE(List)) EXIT
        DO j=i,SIZE(List)-1
          List(j)=List(j+1)
        END DO
        List(SIZE(List))=Temp
        EXIT
      END IF
    END DO
    DO i=1,SIZE(List)
      IF (List(i)%Species=='Hp') THEN
        Temp=List(i)
        IF (i==SIZE(List)) EXIT
        DO j=i,SIZE(List)-1
          List(j)=List(j+1)
        END DO
        List(SIZE(List))=Temp
        EXIT
      END IF
    END DO
  END SUBROUTINE SortList
  !
  !
  SUBROUTINE ListToHashTable(List,Table)
    TYPE(Species_T) :: List(:)
    TYPE(hash_tbl_sll) :: Table
    !
    INTEGER :: i
    DO i=1,SIZE(List)
      CALL Table%put(TRIM(ADJUSTL(List(i)%Species)),i)
    END DO
  END SUBROUTINE ListToHashTable
  !
  !
  SUBROUTINE HashTableToList(Table,List)
    TYPE(hash_tbl_sll) :: Table
    TYPE(Species_T) :: List(:)
    !
    INTEGER :: i,j
    TYPE(sllist), POINTER :: child => NULL()
    !
    DO i=LBOUND(table%vec,dim=1), UBOUND(table%vec,dim=1)
      IF (ALLOCATED(table%vec(i)%key)) THEN
        DO j=1,SIZE(table%vec(i)%key)
          List(table%vec(i)%Val)%Species(j:j)=table%vec(i)%key(j)
        END DO
      END IF
      Child=>table%vec(i)%Child
      DO
        IF (ASSOCIATED(Child)) THEN
          DO j=1,SIZE(Child%key)
            List(Child%Val)%Species(j:j)=Child%key(j)
          END DO
          Child=>Child%Child
        ELSE
          EXIT
        END IF
      END DO
    END DO
  END SUBROUTINE HashTableToList
  !
  !
  SUBROUTINE ReactionListToArray(ReacList,ReacArr)
    TYPE(ListReaction_T) :: ReacList
    TYPE(ReactionStruct_T), ALLOCATABLE, INTENT(OUT) :: ReacArr(:)
    INTEGER :: i, j, TmpArraySize
    INTEGER :: ListLen=0
    !
    Current=>ReacList%Start
    DO 
      IF (.NOT.ASSOCIATED(Current)) EXIT
      ListLen=ListLen+1
      Current=>Current%Next
    END DO
    !
    ALLOCATE(ReacArr(ListLen))
    !
    Current=>ReacList%Start
    i=1
    !
    DO
      IF (ASSOCIATED(Current)) THEN
        ReacArr(i)%Type=Current%Type
        ReacArr(i)%TypeConstant=Current%TypeConstant
        ReacArr(i)%Line1=Current%Line1
        ReacArr(i)%Line2=Current%Line2
        ReacArr(i)%Factor=Current%Factor
        !
        TmpArraySize=SIZE(Current%Educt)
        ALLOCATE(ReacArr(i)%Educt(TmpArraySize))
        DO j=1,TmpArraySize
          ReacArr(i)%Educt(j)%Species=Current%Educt(j)%Species
          ReacArr(i)%Educt(j)%Type=Current%Educt(j)%Type
          ReacArr(i)%Educt(j)%Koeff=Current%Educt(j)%Koeff
        END DO
        !
        TmpArraySize=SIZE(Current%Product)
        ALLOCATE(ReacArr(i)%Product(TmpArraySize))
        DO j=1,TmpArraySize
          ReacArr(i)%Product(j)%Species=Current%Product(j)%Species
          ReacArr(i)%Product(j)%Type=Current%Product(j)%Type
          ReacArr(i)%Product(j)%Koeff=Current%Product(j)%Koeff
        END DO
        !
        ALLOCATE(ReacArr(i)%Constants(SIZE(Current%Constants)))
        ReacArr(i)%Constants=Current%Constants
        !
        TmpArraySize=SIZE(Current%Educt)
        ALLOCATE(ReacArr(i)%InActEduct(TmpArraySize))
        ALLOCATE(ReacArr(i)%InActEductSpc(TmpArraySize))
        DO j=1,TmpArraySize
          ReacArr(i)%InActEduct(j)=Current%InActEduct(j)%Koeff
          ReacArr(i)%InActEductSpc(j)=Current%InActEduct(j)%Species
        END DO
        !
        TmpArraySize=SIZE(Current%InActProduct)
        ALLOCATE(ReacArr(i)%InActProduct(TmpArraySize))
        DO j=1,TmpArraySize
          ReacArr(i)%InActProduct(j)=Current%InActProduct(j)%Koeff    
        END DO
        !    
        ReacArr(i)%nInActEd=Current%nInActEd
        ReacArr(i)%nInActPro=Current%nInActPro
      ELSE
        EXIT
      END IF
      Current=>Current%Next
      i=i+1
    END DO
  END SUBROUTINE ReactionListToArray
  !
  !
  !
  SUBROUTINE AllListsToArray(ReacStruct,LGas,LHenry,LAqua,LDiss,LSolid,LPartic,LMicro)
    ! OUT:
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStruct(:)
    ! IN:
    TYPE(ListReaction_T) :: LGas, LHenry, LAqua, LDiss  &
    &                     , LSolid, LPartic, LMicro
    
    INTEGER :: i, j, iList, iEq
    INTEGER :: nList

    INTEGER :: icnt(47), icntFAC(10), iHen

    ! #Reactions
    neq  = nr_gas + 2*nr_henry + nr_aqua + 2*nr_diss &
    &    + nr_solid + nr_parti + nr_micphys
    nr   = neq
    
    ! #Spezies berechnen
    nspc = ns_GAS + ns_AQUA + ns_SOLID + ns_PARTI
    ns   = nspc
    
    nr_liquid = nr_aqua + 2*nr_diss

    !-----------------------------------------------------------------------
    ! --- Set logicals
    !-----------------------------------------------------------------------
		IF ( ns_GAS   > 0 ) hasGasSpc   = .TRUE.; nPhases = nPhases + 1
		IF ( ns_AQUA  > 0 ) hasAquaSpc  = .TRUE.; nPhases = nPhases + 1
		IF ( ns_SOLID > 0 ) hasSolidSpc = .TRUE.; nPhases = nPhases + 1
		IF ( ns_PARTI > 0 ) hasPartiSpc = .TRUE.; nPhases = nPhases + 1

		IF ( nr_gas    > 0 ) hasGasReac    = .TRUE.
		IF ( nr_aqua   > 0 ) hasAquaReac   = .TRUE.
		IF ( nr_henry  > 0 ) hasHenryReac  = .TRUE.
		IF ( nr_solid  > 0 ) hasSolidReac  = .TRUE.
		IF ( nr_parti  > 0 ) hasPartiReac  = .TRUE.
		IF ( nr_diss   > 0 ) hasDissReac   = .TRUE.
		IF ( nr_liquid > 0 ) hasLiquidReac = .TRUE.

    hasPhotoReac  = (nr_G_photo+nr_A_photo) > 0
		hasFactorReac = nr_FACTOR > 0
    
    nList = 0
    IF (ASSOCIATED(LGas%Start))    nList=nList+1
    IF (ASSOCIATED(LHenry%Start))  nList=nList+1
    IF (ASSOCIATED(LAqua%Start))   nList=nList+1
    IF (ASSOCIATED(LDiss%Start))   nList=nList+1
    IF (ASSOCIATED(LSolid%Start))  nList=nList+1
    IF (ASSOCIATED(LPartic%Start)) nList=nList+1
    IF (ASSOCIATED(LMicro%Start))  nList=nList+1
    ALLOCATE(CompleteReactionList(nList))
    nList = 0
    IF (ASSOCIATED(LGas%Start))   THEN
      nList=nList+1; CompleteReactionList(nList)=LGas
    END IF
    IF (ASSOCIATED(LHenry%Start)) THEN
      !ALLOCATE(isHenry(2*nr_henry,2)); isHenry = 0
      nList=nList+1; CompleteReactionList(nList)=LHenry
    END IF
    IF (ASSOCIATED(LAqua%Start))  THEN
      nList=nList+1; CompleteReactionList(nList)=LAqua
    END IF
    IF (ASSOCIATED(LDiss%Start))  THEN
      nList=nList+1; CompleteReactionList(nList)=LDiss
    END IF
    IF (ASSOCIATED(LSolid%Start))  THEN
      nList=nList+1; CompleteReactionList(nList)=LSolid
    END IF
    IF (ASSOCIATED(LPartic%Start)) THEN
      nList=nList+1; CompleteReactionList(nList)=LPartic
    END IF
    IF (ASSOCIATED(LMicro%Start)) THEN
      nList=nList+1; CompleteReactionList(nList)=LMicro
    END IF

    ALLOCATE( ReacStruct(neq) )
  
    CALL AllocateRTarrays
 
    i=1
    iHen = 0
    icnt = 0
    icntFAC = 0

    DO iList=1,nList
      Current => CompleteReactionList(iList)%Start
      DO WHILE (ASSOCIATED(Current)) 
        ReacStruct(i)%Type   = Current%Type
        ReacStruct(i)%Line1  = Current%Line1
        ReacStruct(i)%Line2  = Current%Line2
        ReacStruct(i)%Line3  = Current%Line3
        ReacStruct(i)%Factor = Current%Factor
        ReacStruct(i)%TypeConstant = Current%TypeConstant
         
        CALL Setup_iFACTOR( i, icntFAC, Current%Factor)

        ! forward direaction
        CALL Setup_ReacParameter( i , icnt ,           &
                                & Current%Type ,       &
                                & Current%TypeConstant,&
                                & Current%Constants,   &
                                & TRIM(Current%Line1)  )
        !
        ReacStruct(i)%nActEd  = SIZE(Current%Educt)
        ReacStruct(i)%nActPro = SIZE(Current%Product)
        ALLOCATE( ReacStruct(i)%Educt(ReacStruct(i)%nActEd),   &
                & ReacStruct(i)%Product(ReacStruct(i)%nActPro) )
            
        ! aus A + B -> .5 C + .5 D + .5 C + .5 D ===> A + B -> C + D

        DO j = 1 , ReacStruct(i)%nActEd
          !write(*,*) 'curr educ = ',i,j,TRIM(Current%Educt(j)%Species)
          ReacStruct(i)%Educt(j)%Species  = Current%Educt(j)%Species
          ReacStruct(i)%Educt(j)%Type     = Current%Educt(j)%Type
          ReacStruct(i)%Educt(j)%Koeff    = Current%Educt(j)%Koeff
          ReacStruct(i)%Educt(j)%iSpecies = PositionSpeciesAll(Current%Educt(j)%Species)
        END DO
        !CALL RemoveDuplicateSpecies(iCol, iVal)
            
        !DO j = 1 , ReacStruct(i)%nActPro
        !  !write(*,*) 'curr prod = ',i,j,TRIM(Current%Product(j)%Species)
        !  ReacStruct(i)%Product(j)%Species  = Current%Product(j)%Species
        !  ReacStruct(i)%Product(j)%Type     = Current%Product(j)%Type
        !  ReacStruct(i)%Product(j)%Koeff    = Current%Product(j)%Koeff
        !  ReacStruct(i)%Product(j)%iSpecies = PositionSpeciesAll(Current%Product(j)%Species)
        !END DO
        ReacStruct(i)%Product = CleanUpDucts(Current%Product)
        ReacStruct(i)%nActPro = SIZE(ReacStruct(i)%Product)
        !IF (i == 103 .OR. i == 102) THEN
        !  WRITE(*,*) 
        !  WRITE(*,*) ' nducts = ', ReacStruct(i)%nActPro
        !  DO j = 1 , ReacStruct(i)%nActPro
        !    !write(*,*) 'curr prod = ',i,j,TRIM(Current%Product(j)%Species)
        !   WRITE(*,*) ReacStruct(i)%Product(j)%Species 
        !    WRITE(*,*) ReacStruct(i)%Product(j)%Type    
        !    WRITE(*,*) ReacStruct(i)%Product(j)%Koeff   
        !    WRITE(*,*) ReacStruct(i)%Product(j)%iSpecies
        !  END DO
        !END IF
        
        IF ( ReacStruct(i)%TypeConstant == 'SPECIAL' ) THEN
          j = SIZE(Current%Special%cVariables)
          IF (Current%Special%Temp) THEN
            ALLOCATE(ReacStruct(i)%Special%cVariables(j),ReacStruct(i)%Special%iVariables(j-1))
          ELSE
            ALLOCATE(ReacStruct(i)%Special%cVariables(j),ReacStruct(i)%Special%iVariables(j))
          END IF

          ReacStruct(i)%Special%nVariables = j
          ReacStruct(i)%Special%Formula  = Current%Special%Formula
          ReacStruct(i)%Special%Temp     = Current%Special%Temp
          DO j = 1,ReacStruct(i)%Special%nVariables
            ReacStruct(i)%Special%cVariables(j) = Current%Special%cVariables(j)
            IF ( Current%Special%cVariables(j) /= 'TEMP' ) THEN 
              ReacStruct(i)%Special%iVariables(j) = PositionSpeciesAll(Current%Special%cVariables(j))
            END IF
          END DO
        ELSE
          ReacStruct(i)%nConst  = SIZE(Current%Constants)
          ALLOCATE(ReacStruct(i)%Constants(ReacStruct(i)%nConst))
          ReacStruct(i)%Constants = Current%Constants
        END IF
        !
        ! DAS HIER MUSS ANDERS WERDEN: INDIZES SO ABSPEICHER DASS SIE AUF INDEX DER TEMP3 REAKTIONEN ZEIGEN
        ! iR%iHENRY(icnt(29),1) müsste auf icnt(7) zeigen
        IF ( Current%Type == 'HENRY' ) THEN
          ReacStruct(i)%direction = 'GA'
          ReacStruct(i)%HenrySpc  = PositionSpeciesAll(ReacStruct(i)%Educt(1)%Species)

          icnt(29) = icnt(29) + 1
          iR%iHENRY(icnt(29),1) = i
          iR%iHENRY(icnt(29),2) = PositionSpeciesAll(ReacStruct(i)%Educt(1)%Species)
          iR%iHENRY(icnt(29),3) = i + 1
          iR%iHENRY(icnt(29),4) = PositionSpeciesAll(ReacStruct(i)%Product(1)%Species)
        END IF
        !
        !
        ReacStruct(i)%nInActEd  = Current%nInActEd
        ReacStruct(i)%nInActPro = Current%nInActPro
        ALLOCATE( ReacStruct(i)%InActEduct(Current%nInActEd),      &
                & ReacStruct(i)%InActEductSpc(Current%nInActEd),   & 
                & ReacStruct(i)%InActProduct(Current%nInActPro),   &
                & ReacStruct(i)%InActProductSpc(Current%nInActPro) )

        DO j = 1 , Current%nInActEd
          ReacStruct(i)%InActEduct(j)    = Current%InActEduct(j)%Koeff
          ReacStruct(i)%InActEductSpc(j) = Current%InActEduct(j)%Species
          nFirst_orderKAT = nFirst_orderKAT + 1
        END DO
        DO j = 1 , Current%nInActPro
          ReacStruct(i)%InActProduct(j)    = Current%InActProduct(j)%Koeff    
          ReacStruct(i)%InActProductSpc(j) = Current%InActProduct(j)%Species
        END DO
        
        ReacStruct(i)%SumAqCoef = SUM(Current%Educt%Koeff) - ONE
       
        IF (ReacStruct(i)%nInActEd > 0 ) THEN
          IF (TRIM(ADJUSTL(ReacStruct(i)%InActEductSpc(1)))=='[aH2O]') THEN 
            ReacStruct(i)%SumAqCoef = ReacStruct(i)%SumAqCoef + ONE
          END IF
        END IF
        IF ( ReacStruct(i)%Type=='AQUA'.OR. ReacStruct(i)%Type=='DISS' ) THEN
          IF ( ReacStruct(i)%SumAqCoef > ZERO ) nr_HOaqua = nr_HOaqua + 1
        END IF
        !
        ! for equilibrium reactions save <-- direction
        SELECT CASE (Current%Type)
          CASE ('DISS','HENRY')
            i=i+1
            iEq = INDEX(Current%Line1,' = ')
            ReacStruct(i)%Type   = Current%Type
            ReacStruct(i)%Line1  = TRIM(Current%Line1(iEq+3:))//' = '//TRIM(Current%Line1(:iEq))
            ReacStruct(i)%Line2  = 'reverse reaction'
            ReacStruct(i)%bR     = .TRUE.
            ReacStruct(i)%Line3  = Current%Line3
            ReacStruct(i)%Factor = Current%Factor
            ReacStruct(i)%TypeConstant = Current%TypeConstant
            CALL Setup_iFACTOR(i,icntFAC,Current%Factor)
           
            ReacStruct(i)%nActEd  = SIZE(Current%Product)
            ReacStruct(i)%nActPro = SIZE(Current%Educt)
            ALLOCATE( ReacStruct(i)%Educt(ReacStruct(i)%nActEd),   &
                    & ReacStruct(i)%Product(ReacStruct(i)%nActPro) )

            DO j=1,ReacStruct(i)%nActEd
              ReacStruct(i)%Educt(j)%Species  = Current%Product(j)%Species
              ReacStruct(i)%Educt(j)%Type     = Current%Product(j)%Type
              ReacStruct(i)%Educt(j)%Koeff    = Current%Product(j)%Koeff
              ReacStruct(i)%Educt(j)%iSpecies = PositionSpeciesAll(Current%Product(j)%Species)
            END DO
            DO j=1,ReacStruct(i)%nActPro
              ReacStruct(i)%Product(j)%Species  = Current%Educt(j)%Species
              ReacStruct(i)%Product(j)%Type     = Current%Educt(j)%Type
              ReacStruct(i)%Product(j)%Koeff    = Current%Educt(j)%Koeff
              ReacStruct(i)%Product(j)%iSpecies = PositionSpeciesAll(Current%Educt(j)%Species)
            END DO
            
            IF ( ReacStruct(i)%TypeConstant == 'SPECIAL' ) THEN
              ReacStruct(i)%Special%nVariables = ReacStruct(i-1)%Special%nVariables
              ReacStruct(i)%Special%Formula    = ReacStruct(i-1)%Special%Formula
              ReacStruct(i)%Special%Temp       = ReacStruct(i-1)%Special%Temp
              ReacStruct(i)%Special%cVariables = ReacStruct(i-1)%Special%cVariables
              ReacStruct(i)%Special%iVariables = ReacStruct(i-1)%Special%iVariables
            ELSE
              ReacStruct(i)%nConst  = SIZE(Current%Constants)
              ALLOCATE(ReacStruct(i)%Constants(ReacStruct(i)%nConst))
              ReacStruct(i)%Constants = Current%Constants
            END IF
            
            IF ( Current%Type == 'HENRY' ) THEN
              ReacStruct(i)%direction = 'AG'
              ReacStruct(i)%HenrySpc  = PositionSpeciesAll(ReacStruct(i)%Product(1)%Species)
            END IF
            !
            ReacStruct(i)%nInActEd  = Current%nInActPro
            ReacStruct(i)%nInActPro = Current%nInActEd
            ALLOCATE( ReacStruct(i)%InActEduct(Current%nInActPro),    &
                    & ReacStruct(i)%InActEductSpc(Current%nInActPro), & 
                    & ReacStruct(i)%InActProduct(Current%nInActEd),   &
                    & ReacStruct(i)%InActProductSpc(Current%nInActEd) )

            DO j = 1 , Current%nInActPro
              ReacStruct(i)%InActEduct(j)    = Current%InActProduct(j)%Koeff
              ReacStruct(i)%InActEductSpc(j) = Current%InActProduct(j)%Species
              nFirst_orderKAT = nFirst_orderKAT + 1
            END DO
            DO j = 1 , Current%nInActEd
              ReacStruct(i)%InActProduct(j)    = Current%InActEduct(j)%Koeff    
              ReacStruct(i)%InActProductSpc(j) = Current%InActEduct(j)%Species
            END DO
            
            ReacStruct(i)%SumAqCoef = SUM(Current%Product%Koeff) - ONE
            
            IF (ReacStruct(i)%nInActEd > 0 ) THEN
              IF (TRIM(ADJUSTL(ReacStruct(i)%InActEductSpc(1)))=='[aH2O]') THEN 
                ReacStruct(i)%SumAqCoef = ReacStruct(i)%SumAqCoef + ONE
              END IF
            END IF
            IF ( ReacStruct(i)%Type=='AQUA'.OR. ReacStruct(i)%Type=='DISS' ) THEN
              IF ( ReacStruct(i)%SumAqCoef > ZERO ) nr_HOaqua = nr_HOaqua + 1
            END IF
            !
        END SELECT
        !
        Current=>Current%Next
        i=i+1
      END DO
    END DO
    
    ! build the array for mass action products of inactive species (katalytic?)
    ALLOCATE( iFO_kat(nfirst_orderKAT,2) )
    nFirst_orderKAT = 0

    ! counting the aquatic reactions with more than one educt
    ALLOCATE( iR%iHOaqua(nr_HOaqua), iR%HOaqua(nr_HOaqua) )
    nr_HOaqua = 0

    DO i=1,neq
      
      IF ( ReacStruct(i)%nInActEd > 0 ) THEN
        nFirst_orderKAT = nFirst_orderKAT + 1
        iFO_kat(nfirst_orderKAT,1) = i
        iFO_kat(nfirst_orderKAT,2) = PositionSpeciesAll(ReacStruct(i)%InActEductSpc(1)) - nspc
      END IF
      
      IF ( ReacStruct(i)%Type=='AQUA'.OR. ReacStruct(i)%Type=='DISS' ) THEN
        IF ( ReacStruct(i)%SumAqCoef > ZERO ) THEN
          nr_HOaqua = nr_HOaqua + 1
          iR%iHOaqua(nr_HOaqua) = i
          iR%HOaqua(nr_HOaqua)  = ReacStruct(i)%SumAqCoef
        END IF
      END IF

    END DO
    
    !print*, ' nHenry = ',nHenry, nHenryga, nHenryag
    !stop

    CONTAINS 

      FUNCTION CleanUpDucts(DuctsIN) RESULT(Ducts)
        TYPE(Duct_T), INTENT(IN)  :: DuctsIN(:)
        TYPE(Duct_T), ALLOCATABLE :: Ducts(:)

        INTEGER :: i, j, n, Len, nDupes
        INTEGER,        ALLOCATABLE :: tmp_iSpc(:) 
        CHARACTER(10),  ALLOCATABLE :: tmp_Type(:)
        CHARACTER(100), ALLOCATABLE :: tmp_cSpc(:)
        REAL(dp),       ALLOCATABLE :: tmp_Koeff(:)
        INTEGER,        ALLOCATABLE :: Perm(:)

        
        n = SIZE(DuctsIN)
        
        IF ( n > 1 ) THEN
        
          ALLOCATE(tmp_cSpc(n), tmp_Type(n), tmp_Koeff(n), tmp_iSpc(n))

          DO i = 1 , n 
            tmp_cSpc(i)  = DuctsIN(i)%Species
            tmp_Type(i)  = DuctsIN(i)%Type
            tmp_Koeff(i) = DuctsIN(i)%Koeff
            tmp_iSpc(i)  = PositionSpeciesAll(DuctsIN(i)%Species)
          END DO 

          CALL SortVecAsc2(Tmp_iSpc,Perm)

          DO i = 1 , n-1
            IF ( tmp_iSpc(i) == tmp_iSpc(i+1) ) THEN
              tmp_Koeff(i+1)  = tmp_Koeff(i+1) + tmp_Koeff(i)
              tmp_iSpc(i) = -1              
            END IF
          END DO

          nDupes = COUNT(tmp_iSpc==-1)

          IF ( nDupes == 0 ) THEN

            ALLOCATE(Ducts(n))
            DO i = 1 , n
              Ducts(i)%iSpecies = tmp_iSpc(i)
              Ducts(i)%Species  = tmp_cSpc(i)(:)
              Ducts(i)%Type     = tmp_Type(i)(:)
              Ducts(i)%Koeff    = tmp_Koeff(i)
            END DO 

          ELSE 

            ALLOCATE(Ducts(n-nDupes))
            Len = 0
            
            DO i = 1 , n
              IF ( tmp_iSpc(i) > 0 ) THEN
                Len = Len + 1
                Ducts(Len)%iSpecies = tmp_iSpc(i)
                Ducts(Len)%Species  = tmp_cSpc(i)(:)
                Ducts(Len)%Type     = tmp_Type(i)(:)
                Ducts(Len)%Koeff    = tmp_Koeff(i)
              END IF
            END DO
              
          END IF

        ELSE
          ALLOCATE(Ducts(1))
          Ducts%iSpecies = PositionSpeciesAll(DuctsIN(1)%Species)
          Ducts%Species  = DuctsIN%Species
          Ducts%Type     = DuctsIN%Type
          Ducts%Koeff    = DuctsIN%Koeff
        END IF
      END FUNCTION CleanUpDucts
  END SUBROUTINE AllListsToArray
  !
  SUBROUTINE Setup_iFACTOR(iReac,icntFAC,Factor)
    CHARACTER(*), INTENT(IN) :: Factor
    INTEGER,      INTENT(IN) :: iReac
    INTEGER,      INTENT(INOUT) :: icntFAC(:)

    SELECT CASE (Factor)
      CASE ('$H2');   icntFAC(1)=icntFAC(1)+1; iR%iFAC_H2(icntFAC(1))=iReac
      CASE ('$O2N2'); icntFAC(2)=icntFAC(2)+1; iR%iFAC_O2N2(icntFAC(2))=iReac
      CASE ('$M');    icntFAC(3)=icntFAC(3)+1; iR%iFAC_M(icntFAC(3))=iReac
      CASE ('$O2');   icntFAC(4)=icntFAC(4)+1; iR%iFAC_O2(icntFAC(4))=iReac
      CASE ('$N2');   icntFAC(5)=icntFAC(5)+1; iR%iFAC_N2(icntFAC(5))=iReac
      CASE ('$H2O');  icntFAC(6)=icntFAC(6)+1; iR%iFAC_H2O(icntFAC(6))=iReac
      CASE ('$RO2');  icntFAC(7)=icntFAC(7)+1; iR%iFAC_RO2(icntFAC(7))=iReac; hasRO2=.TRUE.
      CASE ('$O2O2'); icntFAC(8)=icntFAC(8)+1; iR%iFAC_O2O2(icntFAC(8))=iReac
      CASE ('$aH2O'); icntFAC(9)=icntFAC(9)+1; iR%iFAC_aH2O(icntFAC(9))=iReac
      CASE ('$RO2aq'); icntFAC(10)=icntFAC(10)+1; iR%iFAC_RO2aq(icntFAC(10))=iReac; hasRO2aq=.TRUE.
    END SELECT
  END SUBROUTINE Setup_iFACTOR


  SUBROUTINE Setup_ReacParameter_neu(iReac,icnt,Typ,TypeR,C,Line1)
  REAL(dp), INTENT(IN) :: C(:)
  CHARACTER(*),   INTENT(IN) :: Typ
  CHARACTER(*),   INTENT(IN) :: TypeR
  CHARACTER(*),   INTENT(IN) :: Line1
  INTEGER :: iReac, idx
  INTEGER   :: icnt(:)

  idx = 0 

!  SELECT CASE ( TRIM(TypeR) )
!    CASE ('PHOTABC'); idx = iPHOTABC
!    CASE ('PHOTMCM'); idx = iPHOTMCM
!    CASE ('PHOTAB');  idx = iPHOTAB
!    CASE ('CONST'); idx = iCONST
!    CASE ('TEMP');  idx = iTEMP
!    CASE ('TEMP1'); idx = iTEMP1
!    CASE ('TEMP2'); idx = iTEMP2
!    CASE ('TEMP3'); idx = iTEMP3
!    CASE ('TEMP4'); idx = iTEMP4
!    CASE ('TROE');    idx = iTROE
!    CASE ('TROEF');   idx = iTROEF
!    CASE ('TROEQ');   idx = iTROEQ
!    CASE ('TROEQF');  idx = iTROEQF
!    CASE ('TROEXP');  idx = iTROEXP
!    CASE ('TROEMCM'); idx = iTROEMCM
!    CASE ('SPEC1'); idx = iSPEC1
!    CASE ('SPEC2'); idx = iSPEC2
!    CASE ('SPEC3'); idx = iSPEC3
!    CASE ('SPEC4'); idx = iSPEC4
!    CASE ('SPEC1MCM'); idx = iSPEC1MCM
!    CASE ('SPEC2MCM'); idx = iSPEC2MCM
!    CASE ('SPEC3MCM'); idx = iSPEC3MCM
!    CASE ('SPEC4MCM'); idx = iSPEC4MCM
!    CASE ('SPEC5MCM'); idx = iSPEC5MCM
!    CASE ('SPEC6MCM'); idx = iSPEC6MCM
!    CASE ('SPEC7MCM'); idx = iSPEC7MCM
!    CASE ('SPEC8MCM'); idx = iSPEC8MCM
!    CASE ('SPEC9MCM'); idx = iSPEC9MCM
!    CASE ('S4H2O');  idx = iS4H2O
!    CASE ('T1H2O');  idx = iT1H2O
!    CASE ('ASPEC1'); idx = iASPEC1
!    CASE ('ASPEC2'); idx = iASPEC2
!    CASE ('ASPEC3'); idx = iASPEC3
!    CASE ('ASPEC4'); idx = iASPEC4
!    CASE ('DCONST'); idx = iDCONST
!    CASE ('DTEMP');  idx = iDTEMP
!    CASE ('DTEMP2'); idx = iDTEMP2
!    CASE ('DTEMP3'); idx = iDTEMP3
!    CASE ('DTEMP4'); idx = iDTEMP4
!    CASE ('DTEMP5'); idx = iDTEMP5
!    CASE ('MESKHIDZE'); idx = iMESKHIDZE
!    CASE ('PHOTO');   idx = iPHOTO
!    CASE ('PHOTO2');  idx = iPHOTO2 
!    CASE ('PHOTO3');  idx = iPHOTO3
!    CASE ('SPECIAL'); idx = iSPECIAL
!    CASE ('HOM1');    idx = iHOM1
!    CASE DEFAULT
!      WRITE(*,*) ''
!      WRITE(*,*) ' Reaction Type unknown:  ',TRIM(TypeR),'  --> check input file'
!      WRITE(*,*) ''
!  END SELECT

  IF ( idx > 0 ) THEN
    IF ( SIZE(C) < reac_par(idx)%n_par ) CALL ErrorMSG(iReac,TRIM(Line1))
    icnt(idx) = icnt(idx) + 1
    reac_par(idx)%act = .TRUE.
    !reac_val(idx)%iR( icnt(idx)     ) = iReac 
    reac_val(idx)%vR( icnt(idx) , : ) = C
  END IF



  CONTAINS

    SUBROUTINE ErrorMSG(iR,Line)
      INTEGER      :: iR
      CHARACTER(*) :: Line
      WRITE(*,*);  WRITE(*,*)
      WRITE(*,*) '  ERROR -- > check parameter in reaction: ', iR ,'  ::  '//Line
      WRITE(*,*);  WRITE(*,*)
      STOP
    END SUBROUTINE ErrorMSG

END SUBROUTINE Setup_ReacParameter_neu

  SUBROUTINE Setup_ReacParameter(iReac,icnt,Typ,TypeR,C,Line1)
    REAL(dp), INTENT(IN) :: C(:)
    CHARACTER(*),   INTENT(IN) :: Typ
    CHARACTER(*),   INTENT(IN) :: TypeR
    CHARACTER(*),   INTENT(IN) :: Line1
    INTEGER,        INTENT(IN) :: iReac
    INTEGER,        INTENT(INOUT) :: icnt(47)
    !WRITE(*,*) ' iReac =', iReac, TRIM(ReactionSystem(iReac)%Line1), '    Consts = ', C

    SELECT CASE ( TRIM(TypeR) )
      CASE ('PHOTABC')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(3)=icnt(3)+1; iR%iPHOTabc(icnt(3))=iReac; iR%PHOTabc(icnt(3),:)=C 

        !IF ( SIZE(C) < reac_par(ind_PHOTABC)%n_par ) CALL ErrorMSG(iReac,Line1)
        !reac_par(ind_PHOTABC)%n_reac = reac_par(ind_PHOTABC)%n_reac + 1
        !iR%iPHOTabc( reac_par(ind_PHOTABC)%n_reac     ) = iReac 
        !iR%PHOTabc ( reac_par(ind_PHOTABC)%n_reac , : ) = C 
      CASE ('PHOTMCM')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(4)=icnt(4)+1; iR%iPHOTmcm(icnt(4))=iReac; iR%PHOTmcm(icnt(4),:)=C 
      CASE ('PHOTAB')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(2)=icnt(2)+1; iR%iPHOTab(icnt(2))=iReac;  iR%PHOTab(icnt(2),:)=C
      CASE ('CONST')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
        icnt(1)=icnt(1)+1; iR%iCONST(icnt(1))=iReac;   iR%CONST(icnt(1))=C(1)
      CASE ('TEMP')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(46)=icnt(46)+1; iR%iTEMP1(icnt(46))=iReac; iR%TEMP1(icnt(46),:)=C 
      CASE ('TEMP1')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(5)=icnt(5)+1; iR%iTEMP1(icnt(5))=iReac; iR%TEMP1(icnt(5),:)=C 
      CASE ('TEMP2')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(6)=icnt(6)+1; iR%iTEMP2(icnt(6))=iReac; iR%TEMP2(icnt(6),:)=C 
      CASE ('TEMP3')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(7)=icnt(7)+1; iR%iTEMP3(icnt(7))=iReac; iR%TEMP3(icnt(7),:)=C 
      CASE ('TEMP4')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(8)=icnt(8)+1; iR%iTEMP4(icnt(8))=iReac; iR%TEMP4(icnt(8),:)=C 
      CASE ('TROE')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(9)=icnt(9)+1; iR%iTROE(icnt(9))=iReac;  iR%TROE(icnt(9),:)=C 
      CASE ('TROEF')
        IF ( SIZE(C)<5 ) CALL ErrorMSG(iReac,Line1)
        icnt(10)=icnt(10)+1; iR%iTROEf(icnt(10))=iReac;   iR%TROEf(icnt(10),:)=C  
      CASE ('TROEQ')
        IF ( SIZE(C)<6 ) CALL ErrorMSG(iReac,Line1)
        icnt(11)=icnt(11)+1; iR%iTROEq(icnt(11))=iReac;   iR%TROEq(icnt(11),:)=C 
      CASE ('TROEQF')
        IF ( SIZE(C)<7 ) CALL ErrorMSG(iReac,Line1)
        icnt(12)=icnt(12)+1; iR%iTROEqf(icnt(12))=iReac;  iR%TROEqf(icnt(12),:)=C 
      CASE ('TROEXP')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(13)=icnt(13)+1; iR%iTROExp(icnt(13))=iReac;  iR%TROExp(icnt(13),:)=C 
      CASE ('TROEMCM')
        IF ( SIZE(C)<10 ) CALL ErrorMSG(iReac,Line1)
        icnt(14)=icnt(14)+1; iR%iTROEmcm(icnt(14))=iReac; iR%TROEmcm(icnt(14),:)=C 
      CASE ('SPEC1')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(15)=icnt(15)+1; iR%iSPEC1(icnt(15))=iReac;   iR%SPEC1(icnt(15),:)=C 
      CASE ('SPEC2')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(16)=icnt(16)+1; iR%iSPEC2(icnt(16))=iReac;   iR%SPEC2(icnt(16),:)=C 
      CASE ('SPEC3')
        IF ( SIZE(C)<6 ) CALL ErrorMSG(iReac,Line1)
        icnt(17)=icnt(17)+1; iR%iSPEC3(icnt(17))=iReac;   iR%SPEC3(icnt(17),:)=C 
      CASE ('SPEC4')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(18)=icnt(18)+1; iR%iSPEC4(icnt(18))=iReac;   iR%SPEC4(icnt(18),:)=C 
      CASE ('SPEC1MCM')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(19)=icnt(19)+1; iR%iSPEC1mcm(icnt(19))=iReac; iR%SPEC1mcm(icnt(19),:)=C 
      CASE ('SPEC2MCM')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(20)=icnt(20)+1; iR%iSPEC2mcm(icnt(20))=iReac; iR%SPEC2mcm(icnt(20),:)=C 
      CASE ('SPEC3MCM')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(21)=icnt(21)+1; iR%iSPEC3mcm(icnt(21))=iReac; iR%SPEC3mcm(icnt(21),:)=C 
      CASE ('SPEC4MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(22)=icnt(22)+1; iR%iSPEC4mcm(icnt(22))=iReac; iR%SPEC4mcm(icnt(22),:)=C 
      CASE ('SPEC5MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(23)=icnt(23)+1; iR%iSPEC5mcm(icnt(23))=iReac; iR%SPEC5mcm(icnt(23),:)=C 
      CASE ('SPEC6MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(24)=icnt(24)+1; iR%iSPEC6mcm(icnt(24))=iReac; iR%SPEC6mcm(icnt(24),:)=C 
      CASE ('SPEC7MCM')
        IF ( SIZE(C)<6 ) CALL ErrorMSG(iReac,Line1)
        icnt(25)=icnt(25)+1; iR%iSPEC7mcm(icnt(25))=iReac; iR%SPEC7mcm(icnt(25),:)=C 
      CASE ('SPEC8MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(26)=icnt(26)+1; iR%iSPEC8mcm(icnt(26))=iReac; iR%SPEC8mcm(icnt(26),:)=C 
      CASE ('SPEC9MCM')
        IF ( SIZE(C)<10 ) CALL ErrorMSG(iReac,Line1)
        icnt(47)=icnt(47)+1; iR%iSPEC9mcm(icnt(47))=iReac; iR%SPEC9mcm(icnt(47),:)=C 
      CASE ('S4H2O')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(27)=icnt(27)+1; iR%iS4H2O(icnt(27))=iReac;  iR%S4H2O(icnt(27),:)=C 
      CASE ('T1H2O')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(28)=icnt(28)+1; iR%iT1H2O(icnt(28))=iReac;  iR%T1H2O(icnt(28),:)=C 
      CASE ('ASPEC1')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(30)=icnt(30)+1; iR%iASPEC1(icnt(30))=iReac; iR%ASPEC1(icnt(30),:)=C 
      CASE ('ASPEC2')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(31)=icnt(31)+1; iR%iASPEC2(icnt(31))=iReac; iR%ASPEC2(icnt(31),:)=C 
      CASE ('ASPEC3')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(32)=icnt(32)+1; iR%iASPEC3(icnt(32))=iReac; iR%ASPEC3(icnt(32),:)=C 
      CASE ('ASPEC4')
        !IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac)
        !icnt(33)=icnt(33)+1; iR%iASPEC4(icnt(33))=iReac; iR%ASPEC4(icnt(33),:)=C(1)
      CASE ('DCONST')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(34)=icnt(34)+1
        iR%iDCONST(icnt(34),1)=iReac
        iR%iDCONST(icnt(34),2)=iReac+1
        iR%DCONST(icnt(34),:)=C 
      CASE ('DTEMP')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(35)=icnt(35)+1; iR%iDTEMP(icnt(35),1)=iReac
        iR%iDTEMP(icnt(35),2)=iReac+1; iR%DTEMP(icnt(35),:)=C 
      CASE ('DTEMP2')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(36)=icnt(36)+1; iR%iDTEMP2(icnt(36),1)=iReac
        iR%iDTEMP2(icnt(36),2)=iReac+1 ; iR%DTEMP2(icnt(36),:)=C
      CASE ('DTEMP3')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(37)=icnt(37)+1; iR%iDTEMP3(icnt(37),1)=iReac
        iR%iDTEMP3(icnt(37),2)=iReac+1; iR%DTEMP3(icnt(37),:)=C
      CASE ('DTEMP4')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(38)=icnt(38)+1; iR%iDTEMP4(icnt(38),1)=iReac
        iR%iDTEMP4(icnt(38),2)=iReac+1; iR%DTEMP4(icnt(38),:)=C
      CASE ('DTEMP5')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(39)=icnt(39)+1; iR%iDTEMP5(icnt(39),1)=iReac
        iR%iDTEMP5(icnt(39),2)=iReac+1; iR%DTEMP5(icnt(39),:)=C
      CASE ('MESKHIDZE')
        IF ( SIZE(C)<7 ) CALL ErrorMSG(iReac,Line1)
        icnt(40)=icnt(40)+1; iR%iMeskhidze(icnt(40),1)=iReac
        iR%iMeskhidze(icnt(40),2)=iReac+1; iR%Meskhidze(icnt(40),:)=C 
      CASE ('PHOTO')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
         icnt(41)=icnt(41)+1; iR%iPHOTOkpp(icnt(41))=iReac; iR%PHOTOkpp(icnt(41))=C(1)
      CASE ('PHOTO2')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
        icnt(42)=icnt(42)+1; iR%iPHOTO2kpp(icnt(42))=iReac; iR%PHOTO2kpp(icnt(42))=C(1)
      CASE ('PHOTO3')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
        icnt(43)=icnt(43)+1; iR%iPHOTO3kpp(icnt(43))=iReac; iR%PHOTO3kpp(icnt(43))=C(1)
      CASE ('SPECIAL')
        !IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(44)=icnt(44)+1; iR%iSPECIAL(icnt(44))=iReac
      CASE ('HOM1')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(45)=icnt(45)+1; iR%iHOM1(icnt(45))=iReac; iR%HOM1(icnt(45),:)=C 
      CASE DEFAULT
        WRITE(*,*) ''
        WRITE(*,*) ' Reaction Type unknown:  ',TRIM(TypeR),'  --> check input file'
        WRITE(*,*) ''
    END SELECT

    CONTAINS

      SUBROUTINE ErrorMSG(iR,Line)
        INTEGER      :: iR
        CHARACTER(*) :: Line
        WRITE(*,*);  WRITE(*,*)
        WRITE(*,*) '  ERROR -- > check parameter in reaction: ', iR ,'  ::  '//Line
        WRITE(*,*);  WRITE(*,*)
        STOP
      END SUBROUTINE ErrorMSG

  END SUBROUTINE Setup_ReacParameter
  
  
  SUBROUTINE CheckConstants(RS)
    TYPE(ReactionStruct_T) :: RS(:)     ! reaction system
    CHARACTER(15) :: Const_T            ! constant type for reaction i
    !
    INTEGER :: i,j
    !

    DO i = 1,SIZE(RS)
      Const_T = RS(i)%TypeConstant
      DO j=1,nReacTypes  ! 37 comes from mo_reac def_par
        IF ( TRIM(reac_par(j)%name_type) == ADJUSTL(TRIM(Const_T)) ) EXIT
      END DO
      IF ( SIZE(RS(i)%Constants) /= reac_par(j)%n_par ) THEN
        WRITE(*,*) 'ERROR: Wrong number of constants:'
        WRITE(*,*) '----->  reaction:     ',i, '   ', TRIM(RS(i)%Line1)
        WRITE(*,*) '----->  desired #consts: ', reac_par(j)%n_par, j
        WRITE(*,*) '----->  actual  #consts: ', SIZE(RS(i)%Constants)
        WRITE(*,*) '       Check sys-file for syntax errors!'
        CALL FinishMPI();  STOP
      END IF

    END DO
  END SUBROUTINE CheckConstants


  SUBROUTINE AllocateRTarrays()
    INTEGER :: iType

    ! allocate index arrays and parameter arrays for the new vectorized version
    ALLOCATE( iR%iCONST(nr_CONST), iR%CONST(nr_CONST))

    ALLOCATE( iR%iPHOTabc(nr_PHOTabc), iR%iPHOTab(nr_PHOTab), iR%iPHOTmcm(nr_PHOTmcm))
    ALLOCATE( iR%PHOTabc(nr_PHOTabc,3), iR%PHOTab(nr_PHOTab,2), iR%PHOTmcm(nr_PHOTmcm,3))

    ALLOCATE( iR%iTEMP(nr_TEMP),  iR%iTEMP1(nr_TEMP1),  iR%iTEMP2(nr_TEMP2),  iR%iTEMP3(nr_TEMP3),  iR%iTEMP4(nr_TEMP4))
    ALLOCATE( iR%TEMP(nr_TEMP,3), iR%TEMP1(nr_TEMP1,2), iR%TEMP2(nr_TEMP2,2), iR%TEMP3(nr_TEMP3,2), iR%TEMP4(nr_TEMP4,2))

    ALLOCATE( iR%iTROE(nr_TROE),     iR%iTROEf(nr_TROEf),   iR%iTROEq(nr_TROEq),   &
            & iR%iTROEqf(nr_TROEqf), iR%iTROExp(nr_TROExp), iR%iTROEmcm(nr_TROEmcm))
    ALLOCATE( iR%TROE(nr_TROE,4),     iR%TROEf(nr_TROEf,5),   iR%TROEq(nr_TROEq,6),   &
            & iR%TROEqf(nr_TROEqf,7), iR%TROExp(nr_TROExp,5), iR%TROEmcm(nr_TROEmcm,10))

    ALLOCATE( iR%iSPEC1(nr_SPEC1), iR%iSPEC2(nr_SPEC2), &
            & iR%iSPEC3(nr_SPEC3), iR%iSPEC4(nr_SPEC4))
    ALLOCATE( iR%SPEC1(nr_SPEC1,2), iR%SPEC2(nr_SPEC2,2), iR%SPEC3(nr_SPEC3,6), iR%SPEC4(nr_SPEC4,4))

    ALLOCATE( iR%iSPEC1mcm(nr_SPEC1mcm), iR%iSPEC2mcm(nr_SPEC2mcm),&
            & iR%iSPEC3mcm(nr_SPEC3mcm), iR%iSPEC4mcm(nr_SPEC4mcm),&
            & iR%iSPEC5mcm(nr_SPEC5mcm), iR%iSPEC6mcm(nr_SPEC6mcm),&
            & iR%iSPEC7mcm(nr_SPEC7mcm), iR%iSPEC8mcm(nr_SPEC8mcm),&
            & iR%iSPEC9mcm(nr_SPEC9mcm) )
    ALLOCATE( iR%SPEC1mcm(nr_SPEC1mcm,3), iR%SPEC2mcm(nr_SPEC2mcm,3),&
            & iR%SPEC3mcm(nr_SPEC3mcm,2), iR%SPEC4mcm(nr_SPEC4mcm,4),&
            & iR%SPEC5mcm(nr_SPEC5mcm,4), iR%SPEC6mcm(nr_SPEC6mcm,4),&
            & iR%SPEC7mcm(nr_SPEC7mcm,6), iR%SPEC8mcm(nr_SPEC8mcm,4),&
            & iR%SPEC9mcm(nr_SPEC9mcm,10) )

    ALLOCATE( iR%iS4H2O(nr_S4H2O), iR%iT1H2O(nr_T1H2O))
    ALLOCATE( iR%S4H2O(nr_S4H2O,4), iR%T1H2O(nr_T1H2O,2))

    ALLOCATE( iR%iASPEC1(nr_ASPEC1), iR%iASPEC2(nr_ASPEC2),& 
            & iR%iASPEC3(nr_ASPEC3), iR%iASPEC4(nr_ASPEC4) )
    ALLOCATE( iR%ASPEC1(nr_ASPEC1,2), iR%ASPEC2(nr_ASPEC2,3),&
            & iR%ASPEC3(nr_ASPEC3,2), iR%ASPEC4(nr_ASPEC4,3) ) 

    ALLOCATE( iR%iDTEMP(nr_DTEMP,2),   iR%iDTEMP2(nr_DTEMP2,2),                           &
            & iR%iDTEMP3(nr_DTEMP3,2), iR%iDTEMP4(nr_DTEMP4,2), iR%iDTEMP5(nr_DTEMP5,2) )
    ALLOCATE( iR%DTEMP(nr_DTEMP,3),    iR%DTEMP2(nr_DTEMP2,4),                           &
            & iR%DTEMP3(nr_DTEMP3,4),  iR%DTEMP4(nr_DTEMP4,4),  iR%DTEMP5(nr_DTEMP5,3) )

    ALLOCATE( iR%iDCONST(nr_DCONST,2), iR%iMeskhidze(nr_Meskhidze,2) )
    ALLOCATE( iR%DCONST(nr_DCONST,2),  iR%Meskhidze(nr_Meskhidze,7) )
    
    ALLOCATE( iR%iFAC_H2(nr_FAC_H2), iR%iFAC_O2N2(nr_FAC_O2N2), iR%iFAC_M(nr_FAC_M),        &
            & iR%iFAC_O2(nr_FAC_O2), iR%iFAC_N2(nr_FAC_N2), iR%iFAC_H2O(nr_FAC_H2O),        &
            & iR%iFAC_RO2(nr_FAC_RO2), iR%iFAC_O2O2(nr_FAC_O2O2), iR%iFAC_aH2O(nr_FAC_aH2O),&
            & iR%iFAC_RO2aq(nr_FAC_RO2aq) )

    !  dim1 = iReac, dim2= iSpc, dim3 = iR_g->a / iR_a->g
    ALLOCATE( iR%iHENRY(nr_henry,4) )

    ALLOCATE( iR%iPHOTOkpp(nr_PHOTOkpp), iR%iPHOTO2kpp(nr_PHOTO2kpp), &
            & iR%iPHOTO3kpp(nr_PHOTO3kpp),  iR%PHOTOkpp(nr_PHOTOkpp), &
            & iR%PHOTO2kpp(nr_PHOTO2kpp),  iR%PHOTO3kpp(nr_PHOTO3kpp) )

    ALLOCATE( iR%iSpecial(nr_special) ) 

    ALLOCATE( iR%iHOM1(nr_HOM1) , iR%HOM1(nr_HOM1,3) ) 

  END SUBROUTINE AllocateRTarrays


  SUBROUTINE SearchReactions(Species)
    CHARACTER(*) :: Species
    CHARACTER(80) :: tmpSpc
    INTEGER :: iR, jD, uPath
    INTEGER :: cRcnt, pRcnt

    uPath = 13
    tmpSpc = TRIM(ADJUSTL(Species))
    cRcnt = 0
    pRcnt = 0

    OPEN ( UNIT=uPath , FILE='REACTION_PATHS/'//TRIM(tmpSpc)//'_path.txt' , STATUS='REPLACE' )
    WRITE(uPath,*) ' ********************************************************************************************'
    WRITE(uPath,*) '  '
    WRITE(uPath,*) '  Chemical Mechanism ::              ', TRIM(BSP)
    WRITE(uPath,*) '     System contains ::              ', neq , ' reactions'
    WRITE(uPath,*) '                                     ', nspc, ' species'
    WRITE(uPath,*) '  '
    WRITE(uPath,*) '  All reactions including species :: ', TRIM(tmpSpc)
    WRITE(uPath,*) '  '

    DO iR = 1 , neq
      ! Check educts
      DO jD = 1 , ReactionSystem(iR)%nActEd
        IF (TRIM(ReactionSystem(iR)%Educt(jD)%Species) == TRIM(tmpSpc))   cRcnt = cRcnt + 1
      END DO
      ! Check products
      DO jD = 1 , ReactionSystem(iR)%nActPro
        IF (TRIM(ReactionSystem(iR)%Product(jD)%Species) == TRIM(tmpSpc)) pRcnt = pRcnt + 1
      END DO
    END DO

    WRITE(uPath,*) '    + Number of Reactions where ',TRIM(tmpSpc),' is involved: ', cRcnt+pRcnt
    WRITE(uPath,*) '        - Number of consuming Reactions: ', cRcnt
    WRITE(uPath,*) '        - Number of producing Reactions: ', pRcnt

    DO iR = 1 , neq
      ! Check educts
      DO jD = 1 , ReactionSystem(iR)%nActEd
        IF (TRIM(ReactionSystem(iR)%Educt(jD)%Species) == TRIM(tmpSpc)) THEN
          CALL PrintReaction(iR,uPath)
        END IF
      END DO
      ! Check products
      DO jD = 1 , ReactionSystem(iR)%nActPro
        IF (TRIM(ReactionSystem(iR)%Product(jD)%Species) == TRIM(tmpSpc)) THEN
          CALL PrintReaction(iR,uPath)
        END IF
      END DO
    END DO

    CLOSE( UNIT=13 )

    WRITE(*,*) '  All reactions containing ',TRIM(tmpSpc), &
    &          ' saved in REACTION_PATHs/'//TRIM(tmpSpc)//'_path.txt'
  END SUBROUTINE SearchReactions

  SUBROUTINE PrintReaction(iR,Unit)
    INTEGER :: iR
    INTEGER :: Unit

    WRITE(Unit,*) ''
    WRITE(Unit,*) ' ********************************************************************************************'
    WRITE(Unit,*) '  Reaction Number   :: ', iR
    WRITE(Unit,*) '  Reaction Class    :: ', TRIM(ReactionSystem(iR)%Type)
    WRITE(Unit,*) '  Constant Type     :: ', TRIM(ReactionSystem(iR)%TypeConstant)
    WRITE(Unit,*) '  Reaction          :: ', TRIM(ReactionSystem(iR)%Line1)
    WRITE(Unit,*) '  Order of Reaction :: ', INT(SUM(ReactionSystem(iR)%Educt%Koeff))
    WRITE(Unit,*) '  Factor            :: ', TRIM(ReactionSystem(iR)%Factor)
    WRITE(Unit,*) '  Constants         :: ', ReactionSystem(iR)%Constants
    WRITE(Unit,*) ' ********************************************************************************************'
    WRITE(Unit,*) ''
  END SUBROUTINE PrintReaction

  SUBROUTINE Permute_Species(perm)
    INTEGER, ALLOCATABLE :: perm(:)
    INTEGER :: j
    INTEGER :: par(4) = [ 99999999,99999998,99999997,99999996 ]
    INTEGER, ALLOCATABLE :: p1(:), p2(:), p3(:), p4(:)
    INTEGER, ALLOCATABLE :: idx(:)

    ! Henry species in the middle
    ALLOCATE(idx(nspc))
    idx(1:ns_GAS)  = par(1)
    idx(iR%iHENRY(:,2)) = par(2)
    idx(ns_GAS+1:) = par(4)
    idx(iR%iHENRY(:,4)) = par(3)

    DO j=1,nspc
      IF     ( idx(j) /= par(1) ) THEN
        p1 = [ p1 , j ]
      ELSEIF ( idx(j) /= par(2) ) THEN
        p2 = [ p2 , j ]
      ELSEIF ( idx(j) /= par(3) ) THEN
        p3 = [ p3 , j ]
      ELSEIF ( idx(j) /= par(4) ) THEN
        p4 = [ p4 , j ]
      END IF
    END DO
    perm = [p1, p2, p3, p4]

  END SUBROUTINE Permute_Species


END MODULE Chemsys_Mod
