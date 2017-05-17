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
  USE mo_unirnk
  USE mo_control
  USE mo_reac
  USE InputTool_Mod
  USE NetCDF_Mod
  !USE Functions_Mod
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: LenLine=405
  INTEGER, PARAMETER :: LenName=100
  !
  INTEGER, PARAMETER ::  maxLENinActDuct=9
  ! 
  TYPE Duct_T
    CHARACTER(LenName) :: Species=''
    CHARACTER(20)      :: Type
    REAL(dp)           :: Koeff
    INTEGER            :: iSpecies = 0
  END TYPE Duct_T
  !
  ! LIST FORM
  TYPE Reaction_T
    CHARACTER(20)      :: Type, TypeConstant
    CHARACTER(LenLine) :: Line1, Line2, Line3
    CHARACTER(LenName) :: Factor
    TYPE(Duct_T)  , POINTER   :: Educt(:)=>NULL(), Product(:)=>NULL()
    REAL(dp), POINTER         :: Constants(:)=>NULL()
    TYPE(Duct_T)  , POINTER   :: InActEduct(:)=>NULL(), InActProduct(:)=>NULL()
    INTEGER                   :: nInActEd=0, nInActPro=0
    TYPE(Reaction_T), POINTER :: Next=>NULL()
  END TYPE Reaction_T
  !
  ! ARRAY FORM
  TYPE ReactionStruct_T
    CHARACTER(20)       :: Type,  TypeConstant
    CHARACTER(LenLine)  :: Line1='' , Line2='' , Line3=''
    LOGICAL             :: bR = .FALSE. , brX = .FALSE. 
    CHARACTER(LenName)  :: Factor
    CHARACTER(2)        :: direction
    REAL(dp)      :: SumAqCoef     
    TYPE(Duct_T)  , ALLOCATABLE :: Educt(:), Product(:)
    REAL(dp), ALLOCATABLE       :: Constants(:)
    REAL(dp), ALLOCATABLE       :: LowConst(:), HighConst(:), TroeConst(:) ! combustion press dep reactions
    REAL(dp), ALLOCATABLE       :: InActEduct(:), InActProduct(:)
    INTEGER                     :: nInActEd = 0, nInActPro = 0, nActEd = 0, nActPro = 0
    INTEGER                     :: nConst = 0
    INTEGER                     :: HenrySpc = 0
    LOGICAL                     :: TB = .FALSE. , TBextra=.FALSE.
    INTEGER, ALLOCATABLE        :: TBidx(:)
    CHARACTER(100), ALLOCATABLE :: TBspc(:)
    REAL(dp), ALLOCATABLE       :: TBalpha(:)
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
    REAL(dp)     :: Hf=0.0d0, Gf=0.0d0, Cp=0.0d0
  END TYPE Species_T
  !
  !
  TYPE Element_T
    CHARACTER(5) :: Element=''
  END TYPE Element_T
  !
  !
  TYPE AFRAC_T
    CHARACTER(LenName) :: Species=''
    REAL(dp) :: MolMass           ! [g/mol]
    INTEGER        :: Charge            ! ladung (+,-,++,--,...)
    REAL(dp) :: SolubInd          ! Löslichkeitsindex
    REAL(dp) :: Frac1             ! [g/g]
  END TYPE AFRAC_T
  !
  TYPE(AFRAC_T), ALLOCATABLE :: InitAFrac(:)
  !
  !
  TYPE SPEK_T
    REAL(dp) :: Radius            ! [m]   radius partivle
    REAL(dp) :: wetRadius         ! [m]   radius droplett
    REAL(dp) :: Number            ! [#/cm3]
    REAL(dp) :: Density           ! [kg/m3]
  END TYPE SPEK_T
  !
  TYPE(SPEK_T), ALLOCATABLE :: SPEK(:)
  !
  !
  TYPE NReacType_T
    INTEGER :: GasPhoto, GasPhotAB, GasPhotABC, GasPhotMCM        &
    &        , GasConst, Temp, Temp1, Temp2, Temp3, Troe, Troef   &
    &        , TroeQ, Spec1, Spec2, Spec3, Spec4, Spec1MCM        &
    &        , Spec2MCM, Spec3MCM, SPec4MCM, Spec5MCM, Spec6MCM   &
    &        , Spec7MCM, Spec8MCM, S4H2O, Henry, AquaPhoto        &
    &        , AquaPhotAB, AquaPhotABC, AquaPhotMCM, AquaConst    &
    &        , AquaTemp, AquaTemp1, AquaTemp2, AquaTemp3, Special & 
    &        , DTemp, DTemp1, DTemp2, DTemp3, DTemp4, DTemp5      &
    &        , Meskhidze, Equi, SolidSpecial, Parti, Microphys, SolidDTemp3
  END TYPE NReacType_T
  !
  TYPE(NReacType_T) :: NTypes
  !
  TYPE(Element_T) :: Elements(11)=(/                   &
  &                                 Element_T('(')     &
  &                                ,Element_T(')')     &
  &                                ,Element_T('exp')   &
  &                                ,Element_T('+')     &
  &                                ,Element_T('-')     &
  &                                ,Element_T('*')     &
  &                                ,Element_T('/')     &  
  &                                ,Element_T('**')    &
  &                                ,Element_T('abs')   &
  &                                ,Element_T('sqrt')  &
  &                                ,Element_T('log')   &
  &                                /)
  !
  !
  TYPE(Reaction_T), POINTER   :: System
  TYPE(ListReaction_T), SAVE  :: ListRGas, ListRHenry, ListRAqua,        &
  &                              ListRDiss, ListRSolid, ListRPartic,     &
  &                              ListRMicro
  !
  TYPE(hash_tbl_sll)          :: ListAqua, ListGas, ListSolid,           &
  &                              ListPartic, ListNonReac
  !
  TYPE(Species_T), ALLOCATABLE, TARGET :: ListAqua2(:), ListGas2(:),     &
  &                                       ListSolid2(:), ListPartic2(:), &
  &                                       ListNonReac2(:)
  INTEGER :: InputUnit=10
  INTEGER, PARAMETER :: MaxEduct=10
  INTEGER, PARAMETER :: MaxProduct=10
  !
  CHARACTER(33), PARAMETER :: SetSpecies='ABCDEFGHIJKLMNOPQRSTUVWXYZapsc[]()=+*'
  CHARACTER(14), PARAMETER :: SetConstants='ABCDEFGKINMOR/'
  CHARACTER(12), PARAMETER :: SetExponent='0123456789+-'
  !
  INTEGER ::  NumberSpeciesMicro=0               
  !
  INTEGER ::  NumberReactionsPartic=0            &
  &         , NumberReactionsMicro=0             
  !
  INTEGER :: nsr                        ! # activ species + all Reactions
  !
  INTEGER :: UnitGas=0
  INTEGER :: UnitAqua=0
  !
  CHARACTER(20) :: Filename !='Salt'
  CHARACTER(20) :: IniName  !='Salt'
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
  &                            , InitValInAct(:)    ! initial values inactiv spc
  !
  !
  CHARACTER(LenName), ALLOCATABLE :: RO2spcG(:) , RO2spcA(:)
  INTEGER, ALLOCATABLE :: RO2idxG(:) , RO2idxA(:)
  !
  !
  !REAL(dp) :: aH2O
  !
  REAL(dp), ALLOCATABLE :: sumBAT(:)         ! sum_j=1,n_s (b_ij-a_ij),  i el. N_R
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
    CHARACTER(100) :: Species
    CHARACTER(20) :: Type
    INTEGER :: Pos
    !
    READ(InputUnit,'(a100)',END=1) Species
    DO
      Pos=SCAN(Species,"'")
      IF (Pos>0) THEN
        Species(Pos:)=Species(Pos+1:)
      ELSE
        EXIT
      END IF
    END DO
    IF (Species/='') THEN
      CALL InsertSpecies(Species,Type)
    END IF
    Out=.FALSE.
    GO TO 999
  1 CONTINUE
    Out=.TRUE.
999 CONTINUE
  END SUBROUTINE ReadSpecies
  !
  SUBROUTINE ReadReaction(Out)
    LOGICAL :: Out
    !
    INTEGER :: iLine,PosColon,Pos,is
    CHARACTER(LenLine) :: LocLine
    CHARACTER(LenLine) :: Line(1:4)
    CHARACTER(20) :: Type
    CHARACTER(40) :: TypeR
    INTEGER :: idxFAC
    !
    !
    !
    iLine=0
    DO
      READ(InputUnit,'(a400)',IOSTAT=is) LocLine
      idxFAC = INDEX(LocLine,'$')
      IF ( idxFAC > 0 ) THEN
        SELECT CASE (TRIM(LocLine(idxFAC:)))
          CASE ('$H2','$O2N2','$M','$O2','$N2','$H2O','$O2O2','$aH2O','$+M','$(+M)','$RO2','$RO2aq')
            IF (TRIM(LocLine(idxFAC:))=='$H2')   nFAC_H2   = nFAC_H2 + 1
            IF (TRIM(LocLine(idxFAC:))=='$O2N2') nFAC_O2N2 = nFAC_O2N2 + 1
            IF (TRIM(LocLine(idxFAC:))=='$M')    nFAC_M    = nFAC_M + 1
            IF (TRIM(LocLine(idxFAC:))=='$O2')   nFAC_O2   = nFAC_O2 + 1
            IF (TRIM(LocLine(idxFAC:))=='$N2')   nFAC_N2   = nFAC_N2 + 1
            IF (TRIM(LocLine(idxFAC:))=='$H2O')  nFAC_H2O  = nFAC_H2O + 1
            IF (TRIM(LocLine(idxFAC:))=='$O2O2') nFAC_O2O2 = nFAC_O2O2 + 1
            IF (TRIM(LocLine(idxFAC:))=='$aH2O') nFAC_aH2O = nFAC_aH2O + 2
            IF (TRIM(LocLine(idxFAC:))=='$RO2')  nFAC_RO2  = nFAC_RO2 + 1
            IF (TRIM(LocLine(idxFAC:))=='$RO2aq') nFAC_RO2aq = nFAC_RO2aq + 1
            nFACTOR = nFACTOR + 1
        END SELECT
      END IF
      IF (ABS(is)>0) THEN
        EXIT
      END IF
      IF (ADJUSTL(LocLine(1:1))/='#'.AND.        &
      &   ADJUSTL(LocLine(1:7))/='COMMENT'.AND.  &
      &   LEN(TRIM(LocLine))>0) THEN
        iLine=iLine+1
        Line(iLine)=LocLine
        IF (iLine==4) THEN
          EXIT
        END IF
      END IF
    END DO
    IF (iLine>=3) THEN
      Pos=SCAN(Line(1),'#')
      IF (Pos>0) THEN
        Line(1)=Line(1)(:Pos-1)
      END IF
      !
      IF (INDEX(Line(4),'CLASS')>0) THEN
        BACKSPACE(InputUnit)
      END IF
      PosColon=Index(Line(1),':')
      Type=ADJUSTL(Line(1)(PosColon+1:))
      ! 
      !
      SELECT CASE (Type)
        CASE ('GAS')
          nreakgas=nreakgas+1
          CALL InsertReaction(ListRGas,Line,TypeR)
          !print*, 'gasTypeR = ', TypeR
          SELECT CASE (TypeR)
            CASE ('PHOTO','PHOTO2','PHOTO3') ! kpp style photolytic reactions
              IF (TypeR=='PHOTO')  nPHOTOkpp  = nPHOTOkpp  + 1
              IF (TypeR=='PHOTO2') nPHOTO2kpp = nPHOTO2kpp + 1
              IF (TypeR=='PHOTO3') nPHOTO3kpp = nPHOTO3kpp + 1
              nreakgphoto = nreakgphoto+1
            CASE ('PHOTAB','PHOTABC','PHOTMCM')
              IF (TypeR=='PHOTAB')   nPHOTab  = nPHOTab  + 1
              IF (TypeR=='PHOTABC')  nPHOTabc = nPHOTabc + 1
              IF (TypeR=='PHOTMCM')  nPHOTmcm = nPHOTmcm + 1
              nreakgphoto = nreakgphoto+1
            CASE ('CONST')
              nCONST = nCONST + 1
              nreakgconst = nreakgconst+1
            CASE ('TEMP1','TEMP2','TEMP3')
              IF (TypeR=='TEMP1') nTEMP1 = nTEMP1 + 1
              IF (TypeR=='TEMP2') nTEMP2 = nTEMP2 + 1
              IF (TypeR=='TEMP3') nTEMP3 = nTEMP3 + 1
              IF (TypeR=='TEMP4') nTEMP4 = nTEMP4 + 1
              nreakgtemp = nreakgtemp + 1
            CASE ('TROE','TROEF','TROEQ','TROEQF','TROEXP','TROEMCM')
              IF (TypeR=='TROE')    nTROE    = nTROE    + 1
              IF (TypeR=='TROEF')   nTROEf   = nTROEf   + 1
              IF (TypeR=='TROEQ')   nTROEq   = nTROEq   + 1
              IF (TypeR=='TROEQF')  nTROEqf  = nTROEqf  + 1
              IF (TypeR=='TROEXP')  nTROExp  = nTROExp  + 1
              IF (TypeR=='TROEMCM') nTROEmcm = nTROEmcm + 1
              nreakgtroe = nreakgtroe + 1
            CASE ('SPEC1','SPEC2','SPEC3','SPEC4','SPEC1MCM','SPEC2MCM',            &
            &     'SPEC3MCM','SPEC4MCM','SPEC5MCM','SPEC6MCM','SPEC7MCM','SPEC8MCM')
              IF (TypeR=='SPEC1') nSPEC1 = nSPEC1 + 1
              IF (TypeR=='SPEC2') nSPEC2 = nSPEC2 + 1
              IF (TypeR=='SPEC3') nSPEC3 = nSPEC3 + 1
              IF (TypeR=='SPEC4') nSPEC4 = nSPEC4 + 1
              IF (TypeR=='SPEC1MCM') nSPEC1mcm = nSPEC1mcm + 1
              IF (TypeR=='SPEC2MCM') nSPEC2mcm = nSPEC2mcm + 1
              IF (TypeR=='SPEC3MCM') nSPEC3mcm = nSPEC3mcm + 1
              IF (TypeR=='SPEC4MCM') nSPEC4mcm = nSPEC4mcm + 1
              IF (TypeR=='SPEC5MCM') nSPEC5mcm = nSPEC5mcm + 1
              IF (TypeR=='SPEC6MCM') nSPEC6mcm = nSPEC6mcm + 1
              IF (TypeR=='SPEC7MCM') nSPEC7mcm = nSPEC7mcm + 1
              IF (TypeR=='SPEC8MCM') nSPEC8mcm = nSPEC8mcm + 1
              nreakgspec = nreakgspec + 1
            CASE ('S4H2O')
              nS4H2O = nS4H2O + 1
            CASE ('T1H2O')
              nT1H2O = nT1H2O + 1
          END SELECT
        CASE ('HENRY')
          nreakhenry = nreakhenry + 1
          nHENRY = nHENRY + 1
          CALL InsertReaction(ListRHenry,Line,TypeR)
          IF (TypeR=='TEMP1') nTEMP1 = nTEMP1 + 1
          IF (TypeR=='TEMP2') nTEMP2 = nTEMP2 + 1
          IF (TypeR=='TEMP3') nTEMP3 = nTEMP3 + 1
          IF (TypeR=='CONST') nCONST = nCONST + 1
        CASE ('AQUA')
          nreakaqua = nreakaqua+1
          CALL InsertReaction(ListRAqua,Line,TypeR)
          SELECT CASE (TypeR)
            CASE ('PHOTAB','PHOTABC','PHOTMCM')
              IF (TypeR=='PHOTAB')   nPHOTab  = nPHOTab  + 1
              IF (TypeR=='PHOTABC')  nPHOTabc = nPHOTabc + 1
              IF (TypeR=='PHOTMCM')  nPHOTmcm = nPHOTmcm + 1
              nreakaphoto = nreakaphoto + 1
            CASE ('CONST')
              nCONST = nCONST + 1
              nreakaconst = nreakaconst + 1
            CASE ('TEMP','Temp1''TEMP2','TEMP3')
              IF (TypeR=='TEMP1') nTEMP1 = nTEMP1 + 1
              IF (TypeR=='TEMP2') nTEMP2 = nTEMP2 + 1
              IF (TypeR=='TEMP3') nTEMP3 = nTEMP3 + 1
              IF (TypeR=='TEMP4') nTEMP4 = nTEMP4 + 1
              nreakatemp = nreakatemp + 1
            CASE ('ASPEC1','ASPEC2','ASPEC3','ASPEC4')
              IF (TypeR=='ASPEC1') nASPEC1 = nASPEC1 + 1
              IF (TypeR=='ASPEC2') nASPEC2 = nASPEC2 + 1
              IF (TypeR=='ASPEC3') nASPEC3 = nASPEC3 + 1
              IF (TypeR=='ASPEC4') nASPEC4 = nASPEC4 + 1
            !CASE ('SPECIAL')
            !  NTypes%Special=NTypes%Special+1
            !  nreakaspec=nreakaspec+1
          END SELECT
        CASE ('DISS')
          nreakdissoc=nreakdissoc+1
          CALL InsertReaction(ListRDiss,Line,TypeR)
          SELECT CASE (TypeR)
            CASE ('DCONST','DTEMP','DTEMP2','DTEMP3','DTEMP4','DTEMP5','MESKHIDZE')
              IF (TypeR=='DCONST')   nDCONST = nDCONST + 1
              IF (TypeR=='DTEMP')    nDTEMP  = nDTEMP  + 1
              IF (TypeR=='DTEMP2')   nDTEMP2 = nDTEMP2 + 1
              IF (TypeR=='DTEMP3')   nDTEMP3 = nDTEMP3 + 1
              IF (TypeR=='DTEMP4')   nDTEMP4 = nDTEMP4 + 1
              IF (TypeR=='DTEMP5')   nDTEMP5 = nDTEMP5 + 1
              IF (TypeR=='MESKHIDZE') nMeskhidze = nMeskhidze + 1
          END SELECT
        !CASE ('SOLID')
        !  nreaksolid = nreaksolid + 1
        !  CALL InsertReaction(ListRSolid,Line,TypeR)
        !  SELECT CASE (TypeR)
        !    CASE ('EQUI')
        !      NTypes%Equi=NTypes%Equi+1
        !      nreaksolidEqui=nreaksolidEqui+1
        !    CASE ('DTEMP3')
        !      NTypes%SolidDTemp3=NTypes%SolidDTemp3+1
        !      nreaksolidtemp=nreaksolidtemp+1
        !    CASE ('SPECIAL')
        !      NTypes%SolidSpecial=NTypes%SolidSpecial+1
        !      nreaksolidspec=nreaksolidspec+1
        !  END SELECT
        !CASE ('PARTI')
        !  NTypes%Parti=NTypes%Parti+1
        !  NumberReactionsPartic=NumberReactionsPartic+1
        !  CALL InsertReaction(ListRPartic,Line,TypeR)
        !CASE ('MICROPHYS')
        !  NTypes%Microphys=NTypes%Microphys+1
        !  NumberReactionsMicro=NumberReactionsMicro+1
        !  CALL InsertReaction(ListRMicro,Line,TypeR)
      CASE DEFAULT
        WRITE(*,*) '  Unknown reaction type: ', Type
      END SELECT
      Out=.FALSE.
    ELSE
      Out=.TRUE.
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
  SUBROUTINE PrintSpecies(ListName,Unit)
    TYPE(Species_T) :: ListName(:)
    INTEGER :: Unit
    !
    INTEGER :: i
    !
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
    WRITE(Unit,*) ' Created:             ',Date(7:8),'.',Date(5:6),'.',Date(1:4)
    WRITE(Unit,*) ' Chemical Mechanism:  '                                     &
    &             , TRIM(ADJUSTL(FileName(INDEX(FileName,'/',.TRUE.)+1:))),'.chem'
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
    WRITE(Unit,*) ntGAS &
                 +ntAQUA &
                 +ntPart &
                 +ntKAT &
                 +ntSOLID,  '     Number of Species'
    WRITE(Unit,*) ntGAS,    '     No. of gaseous species'
    WRITE(Unit,*) ntAQUA,   '     No. of aqueous species'
    WRITE(Unit,*) ntPart, '     No. of particular species'
    WRITE(Unit,*) ntSOLID,  '     No. of solid   species'
    WRITE(Unit,*) ntKAT,'     Number of Non-reactive Species '
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
    !
    !
    nreak=nreakgas               &
    &    +nreakhenry             &
    &    +nreakaqua              &
    &    +nreakdissoc            &
    &    +nreaksolid             &
    &    +NumberReactionsPartic  &
    &    +NumberReactionsMicro
    !
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' ================   Description of Reactions   =============='
    WRITE(Unit,*) ''
    WRITE(Unit,*) nreak,          '        NREAK   : Number of Reactions'
    WRITE(Unit,*) nreakgas,       '        NGAS   : Gas phase reactions'
    WRITE(Unit,*) nreakgphoto,    '           Gaseous PHOTO - type reactions'
    WRITE(Unit,*) nreakgconst,    '           Gaseous CONST - type reactions'
    WRITE(Unit,*) nreakgtemp,     '           Gaseous TEMP - type reactions'
    WRITE(Unit,*) nreakSimpTB,    '           Gaseous Simple three-body - type reactions'
    WRITE(Unit,*) nreakglind,     '           Gaseous Lindemann - type reactions'
    WRITE(Unit,*) nreakgtroe,     '           Gaseous TROE - type reactions'
    WRITE(Unit,*) nreakgspec,     '           Gaseous SPECIAL - type reactions'
    WRITE(Unit,*) nreakhenry,     '        NHENRY : Henry Equilib. reactions'
    WRITE(Unit,*) nreakdissoc,    '        NDISS  : Dissociation reactions'
    WRITE(Unit,*) nreakaqua,      '        NAQUA  : Aquatic Equilib. reactions'
    WRITE(Unit,*) nreakaphoto,    '           Aqueous PHOTO - type reactions'
    WRITE(Unit,*) nreakaconst,    '           Aqueous CONST - type reactions'
    WRITE(Unit,*) nreakatemp,     '           Aqueous TEMP - type reactions'
    WRITE(Unit,*) nreakaspec,     '           Aqueous SPECIAL - type reactions'
    WRITE(Unit,*) NumberReactionsPartic,      '        NPARTI  : Particulare reactions   '
    WRITE(Unit,*) nreaksolid,     '        NSOLID  : Solid Equilib. reactions'
    WRITE(Unit,*) nreaksolidtemp, '           Solid DTEMP3 - type reactions'
    WRITE(Unit,*) nreaksolidEqui, '           Solid EQUI - type reactions'
    WRITE(Unit,*) nreaksolidspec, '           Solid SPECIAL - type reactions'
    WRITE(Unit,*) NumberReactionsMicro,       '        NMICRO  : Microphysical reactions'
    WRITE(Unit,*)
    WRITE(Unit,*) ' ======================  Reactions   ========================'
    WRITE(Unit,*) ''
  END SUBROUTINE PrintHeadReactions
  !
  SUBROUTINE PrintSortedReactions(ReacStruct)
    TYPE(ReactionStruct_T), INTENT(IN) :: ReacStruct(:)
    INTEGER :: Unit=989
    !
    INTEGER :: i
    !
    OPEN(UNIT=989,FILE=ADJUSTL(TRIM(ChemFile))//'_sorted.sys',STATUS='UNKNOWN')
    !
    WRITE(Unit,*) '# ================= Sorted  '//TRIM(BSP)//'.sys ================='
    WRITE(Unit,*) '# = Please copy the data into your sys-file for ='
    WRITE(Unit,*) '# =============== chemical input. ==============='
    WRITE(Unit,*) '#'
    WRITE(Unit,*) '#  ===================   Unit options   ======================'
    WRITE(Unit,*) ''
    WRITE(Unit,*) 'UNIT GAS    0   #    Gas phase units     (0 = molec/cm3, 1 = mol/m3)'
    WRITE(Unit,*) 'UNIT AQUA   0   #    Aqueous phase units (0 = mol/l)'
    WRITE(Unit,*) 'UNIT AQUA   0   #    Aqueous phase units (0 = mol/l)'
    WRITE(Unit,*) ''
    WRITE(Unit,*) '#  ===================   Gas Phase      ======================'
    WRITE(Unit,*) '#'
    !
    i=0
    DO 
      i=i+1
      WRITE(Unit,*) 'CLASS: ' ,TRIM(ReacStruct(i)%Type)
      WRITE(Unit,*) TRIM(ReacStruct(i)%Line1)
      WRITE(Unit,*) TRIM(ReacStruct(i)%Line3)
      WRITE(Unit,*) 'FACTOR: ',TRIM(ReacStruct(i)%Factor)
      WRITE(Unit,*) ''
      IF (i>=SIZE(ReacStruct)) EXIT
      IF (ReacStruct(i)%Type=='DISS'.OR.ReacStruct(i)%Type=='HENRY') i=i+1
    END DO
    CLOSE(989)
  END SUBROUTINE PrintSortedReactions
  !
  SUBROUTINE PrintReactions(ReacStruct,Unit,CK)
    !TYPE(ListReaction_T) :: List
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStruct(:)
    INTEGER :: Unit
    LOGICAL, OPTIONAL :: CK
    !
    INTEGER :: i,k,l,m,iLoop
    INTEGER, SAVE :: NumberReaction=0
    INTEGER :: NumActiveEduct,NumActiveProduct
    TYPE(Duct_T) :: ActiveEduct(30)
    TYPE(Duct_T) :: ActiveProduct(30)
    !
    INTEGER :: ValCntA=1
    INTEGER :: ValCntB=1
    INTEGER :: EductCnt=0
    INTEGER :: ProductCnt=0
    !
    INTEGER, ALLOCATABLE :: tmpColA(:),tmpColB(:)
    REAL(dp), ALLOCATABLE :: tmpValA(:),tmpValB(:)
    INTEGER, ALLOCATABLE :: APermVec(:)
    INTEGER :: AColLen
    INTEGER :: cntEQ1, cntEQ2, cntNEQ1 
    !
    !
    ! set matrix dimensions
    A%m=neq
    A%n=nspc
    B%m=neq
    B%n=nspc
    BA%m=neq
    BA%n=nspc
    !
    ! Standart alloc
    ALLOCATE(A%RowPtr(A%m+1))
    A%RowPtr=0
    A%RowPtr(1)=1
    ALLOCATE(B%RowPtr(B%m+1))
    B%RowPtr=0
    B%RowPtr(1)=1
    ALLOCATE(BA%RowPtr(BA%m+1))
    BA%RowPtr=0
    BA%RowPtr(1)=1
    
    !
    DO iLoop=1,neq
      ! count activ educts in reaction iLoop
      NumActiveEduct=0
      !print*, 'DEBUG::chemsys    sizeRSe,p= ',iloop,SIZE(ReactionSystem(iLoop)%Educt),SIZE(ReactionSystem(iLoop)%Product)
      !print*, 'DEBUG::chemsys    reaktion = ',TRIM(ReactionSystem(iLoop)%Line1)
      DO i=1,SIZE(ReactionSystem(iLoop)%Educt)
        SELECT CASE(ReactionSystem(iLoop)%Educt(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','GAS')
            NumActiveEduct=NumActiveEduct+1
            ActiveEduct(NumActiveEduct)=ReactionSystem(iLoop)%Educt(i)
            !print*, 'debug::chemssys   ActiveEduct(NumActiveEduct)=',ActiveEduct(NumActiveEduct)
        END SELECT
      END DO
      ! count activ products in reaction iLoop
      NumActiveProduct=0
      DO i=1,SIZE(ReactionSystem(iLoop)%Product)
        SELECT CASE(ReactionSystem(iLoop)%Product(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','GAS')
            NumActiveProduct=NumActiveProduct+1
            ActiveProduct(NumActiveProduct)=ReactionSystem(iLoop)%Product(i)
            !print*, 'debug::chemssys   ActiveProduct(NumActiveProduct)=',ActiveProduct(NumActiveProduct)
        END SELECT
      END DO
      !
      NumberReaction=NumberReaction+1
      WRITE(Unit,'(a12,i5,a23)')                                                  &
      &             '#-----------',NumberReaction,'. Reaction ----------- '
      !
      WRITE(Unit,*) TRIM(ReactionSystem(iLoop)%Type)//'   '//                   &
      &             TRIM(ReactionSystem(iLoop)%TypeConstant)
      !
      WRITE(Unit,*) SIZE(ReactionSystem(iLoop)%Educt),                          &
      &             SIZE(ReactionSystem(iLoop)%Product),                        &
      &             NumActiveEduct,NumActiveProduct
      !
      ! set RowPtr of A and B
      ! Educt Matrix A
      IF (NumActiveEduct>1) THEN
        ALLOCATE(tmpColA(NumActiveEduct))
        ALLOCATE(tmpValA(NumActiveEduct))
        tmpColA=0
        tmpValA=0.0d0
        DO l=1,NumActiveEduct
          tmpColA(l)=PositionSpeciesAll(ActiveEduct(l)%Species)
          tmpValA(l)=ActiveEduct(l)%Koeff
        END DO
        CALL CompressList(tmpColA,tmpValA)
        A%RowPtr(NumberReaction+1)=A%RowPtr(NumberReaction)+SIZE(tmpColA)
      ELSE
        A%RowPtr(NumberReaction+1)=A%RowPtr(NumberReaction)+NumActiveEduct
      END IF
      ! Product Matrix B
      IF (NumActiveProduct>1) THEN
        ALLOCATE(tmpColB(NumActiveProduct))
        ALLOCATE(tmpValB(NumActiveProduct))
        tmpColB=0
        tmpValB=0.0d0
        DO l=1,NumActiveProduct
          tmpColB(l)=PositionSpeciesAll(ActiveProduct(l)%Species)
          tmpValB(l)=ActiveProduct(l)%Koeff
        END DO
        CALL CompressList(tmpColB,tmpValB)
        B%RowPtr(NumberReaction+1)=B%RowPtr(NumberReaction)+SIZE(tmpColB)
      ELSE
        B%RowPtr(NumberReaction+1)=B%RowPtr(NumberReaction)+NumActiveProduct
      END IF
      IF (ALLOCATED(tmpColA)) DEALLOCATE(tmpColA)
      IF (ALLOCATED(tmpValA)) DEALLOCATE(tmpValA)
      IF (ALLOCATED(tmpColB)) DEALLOCATE(tmpColB)
      IF (ALLOCATED(tmpValB)) DEALLOCATE(tmpValB)
      ! ----------------------------------------------------
      ! SpeziesIndx Edukt=> 1:#Edukt von Reaktion NumberReaction
      ! SpeziesIndx Produkt=> #Edukt+1:#Edukt+#Produkt von Reaktion NumberReaction
      ! #aktiver Stoffe der Reaktion
      WRITE(Unit,*) (PositionSpeciesAll(ReactionSystem(iLoop)%Educt(i)%Species),   &
      &             i=1,SIZE(ReactionSystem(iLoop)%Educt)),                     &
      &             (PositionSpeciesAll(ReactionSystem(iLoop)%Product(i)%Species), &
      &             i=1,SIZE(ReactionSystem(iLoop)%Product)),                   &
      &             NumActiveEduct+NumActiveProduct
      ! 
      !----------------------------------------------------
      ! Tupel: (SpeziesIndex,-Koeffzien) für alle aktiven Edukte (links)
      ! Tupel: (SpeziesIndex,+Koeffzien) für alle aktiven Produkte (rechts)
      !WRITE(Unit,'(7X,I5,3X,F6.3)',ADVANCE='NO')                                  &
      !&                   ( PositionSpeciesAll(ActiveEduct(i)%Species)            &
      !&                  ,  -ActiveEduct(i)%Koeff,i=1,NumActiveEduct     )        
      !WRITE(Unit,'(7X,I5,3X,F6.3)',ADVANCE='NO')                                  &
      !&                  ,( PositionSpeciesAll(ActiveProduct(i)%Species)          &
      !&                  ,   ActiveProduct(i)%Koeff,i=1,NumActiveProduct )
      WRITE(Unit,*)                                 &
      &                   ( PositionSpeciesAll(ActiveEduct(i)%Species)            &
      &                  ,  -ActiveEduct(i)%Koeff,i=1,NumActiveEduct     )        &
      &                  ,( PositionSpeciesAll(ActiveProduct(i)%Species)          &
      &                  ,   ActiveProduct(i)%Koeff,i=1,NumActiveProduct )
      !
      IF (ReactionSystem(iLoop)%TypeConstant=='SPECIAL') THEN
        WRITE(Unit,*) ReactionSystem(iLoop)%                                    &
        &             Line2(:LEN(TRIM(ReactionSystem(iLoop)%Line2))-1)
      END IF
      !
      ! #Reaktionskonstanten, Reaktionskonstanten 1:#
      WRITE(Unit,*) SIZE(ReactionSystem(iLoop)%Constants),                      &
                    ReactionSystem(iLoop)%Constants
      WRITE(Unit,*) 'FACTOR:  ',ReactionSystem(iLoop)%Factor

      SELECT CASE (ReactionSystem(iLoop)%Factor)
      !  CASE ('$H2','$O2N2','$M','$O2','$N2','$H2O','aH2O','$+M','$(+M)')
      !    nFACTOR = nFACTOR + 1
        CASE ('$RO2')
          hasRO2  = .TRUE.
      !    nFACTOR = nFACTOR + 1
        CASE ('$RO2aq')
          hasRO2aq = .TRUE.
      !    nFACTOR  = nFACTOR + 1
      END SELECT

      !
      IF (PRESENT(CK)) WRITE(Unit,*) 'EXTRA1:  ',ADJUSTL(TRIM(ReactionSystem(iLoop)%Line2))
      IF (PRESENT(CK)) WRITE(Unit,*) 'EXTRA2:  ',ADJUSTL(TRIM(ReactionSystem(iLoop)%Line3))
    END DO
    !
    ! loop again to set ColInd and Val on A and B
    NumberReaction=0
    ValCntA=0
    ValCntB=0
    EductCnt=0
    ProductCnt=0
    k=1
    !
    ! standart alloc 
    ALLOCATE(A%ColInd(A%RowPtr(A%m+1)-1))
    ALLOCATE(A%Val(A%RowPtr(A%m+1)-1))
    A%ColInd=0
    A%Val=ZERO
    !
    ALLOCATE(B%ColInd(B%RowPtr(B%m+1)-1))
    ALLOCATE(B%Val(B%RowPtr(B%m+1)-1))
    B%ColInd=0
    B%Val=ZERO
    !
    ALLOCATE(sumBAT(A%m))
    sumBAT=ZERO
    !
    DO iLoop=1,neq
      NumActiveEduct=0
      DO i=1,SIZE(ReactionSystem(iLoop)%Educt)
        SELECT CASE(ReactionSystem(iLoop)%Educt(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','GAS')
            NumActiveEduct=NumActiveEduct+1
            ActiveEduct(NumActiveEduct)=ReactionSystem(iLoop)%Educt(i)
        END SELECT
      END DO
      NumActiveProduct=0
      DO i=1,SIZE(ReactionSystem(iLoop)%Product)
        SELECT CASE(ReactionSystem(iLoop)%Product(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','GAS')
            NumActiveProduct=NumActiveProduct+1
            ActiveProduct(NumActiveProduct)=ReactionSystem(iLoop)%Product(i)
        END SELECT
      END DO
      !
      !
      NumberReaction=NumberReaction+1
      ! set ColInd and Val on A and B
      IF (NumActiveEduct>1) THEN
        ALLOCATE(tmpColA(NumActiveEduct))
        ALLOCATE(tmpValA(NumActiveEduct))
        ALLOCATE(APermVec(NumActiveEduct))
        tmpColA=0
        tmpValA=ZERO
        APermVec=0
        DO l=1,NumActiveEduct
          tmpColA(l)=PositionSpeciesAll(ActiveEduct(l)%Species)
          tmpValA(l)=ActiveEduct(l)%Koeff
        END DO
        !
        ! sort ColInd and Val for acc column indx
        CALL unirnk(tmpColA,APermVec,AColLen)
        tmpColA=tmpColA(APermVec)
        tmpValA=tmpValA(APermVec)
        CALL CompressList(tmpColA,tmpValA)
        !

        DO m=1,AColLen
          ValCntA=ValCntA+1
          A%ColInd(ValCntA)=tmpColA(m)
          A%Val(ValCntA)=tmpValA(m)
          
          ! this is for the TempX  backward reaction 
          sumBAT(iLoop)=sumBAT(iLoop)-tmpValA(m)
        END DO
        DEALLOCATE(APermVec)
      ELSE
        ! reactions with only one educt
        DO i=1,NumActiveEduct
          ValCntA=ValCntA+1
          A%ColInd(ValCntA)=PositionSpeciesAll(ActiveEduct(i)%Species)
          A%Val(ValCntA)=ActiveEduct(i)%Koeff
          sumBAT(iLoop)=sumBAT(iLoop)-ActiveEduct(i)%Koeff
        END DO
      END IF
      !
      IF (NumActiveProduct>1) THEN
         ALLOCATE(tmpValB(NumActiveProduct))
         ALLOCATE(tmpColB(NumActiveProduct))
         ALLOCATE(APermVec(NumActiveProduct))
        tmpColB=0
        tmpValB=ZERO
        APermVec=0
        DO l=1,NumActiveProduct
          tmpColB(l)=PositionSpeciesAll(ActiveProduct(l)%Species)
          tmpValB(l)=ActiveProduct(l)%Koeff
        END DO
        CALL unirnk(tmpColB,APermVec,AColLen)
        tmpColB=tmpColB(APermVec)
        tmpValB=tmpValB(APermVec)
        CALL CompressList(tmpColB,tmpValB)
        !
        DO m=1,AColLen
          ValCntB=ValCntB+1
          B%ColInd(ValCntB)=tmpColB(m)
          B%Val(ValCntB)=tmpValB(m)
          !
          ! this is for the TempX  backward reaction 
          sumBAT(iLoop)=sumBAT(iLoop)+tmpValB(m)
        END DO
        DEALLOCATE(APermVec)
      ELSE
        !IF (NumbActiveProduct>0) THEN
        DO i=1,NumActiveProduct
          ValCntB=ValCntB+1
          B%ColInd(ValCntB) = PositionSpeciesAll(ActiveProduct(i)%Species)
          B%Val(ValCntB) = ActiveProduct(i)%Koeff
          sumBAT(iLoop)  = sumBAT(iLoop) + ActiveProduct(i)%Koeff
        END DO
        !END IF
      END IF
      IF (ALLOCATED(tmpColB)) DEALLOCATE(tmpColB)
      IF (ALLOCATED(tmpValB)) DEALLOCATE(tmpValB)
      IF (ALLOCATED(tmpColA)) DEALLOCATE(tmpColA)
      IF (ALLOCATED(tmpValA)) DEALLOCATE(tmpValA)
    END DO

    A%nnz=A%RowPtr(A%m+1)-1
    B%nnz=B%RowPtr(B%m+1)-1
  END SUBROUTINE PrintReactions

  SUBROUTINE GatherSpeciesOrder(A)
    
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    INTEGER :: iR, j, jj
    INTEGER :: nnz, cnt
    INTEGER, ALLOCATABLE :: tmpFO1(:), tmpFO2(:)
    INTEGER, ALLOCATABLE :: tmpSO1(:), tmpSO2(:)
    INTEGER, ALLOCATABLE :: tmpHO1(:), tmpHO2(:)
    REAL(dp), ALLOCATABLE :: atmpHO(:)
    REAL(dp), PARAMETER :: big = -99999999999999.d0

    nnz = A%nnz

    ALLOCATE(tmpFO1(nnz),tmpFO2(nnz),&
            &tmpSO1(nnz),tmpSO2(nnz),&
            &tmpHO1(nnz),tmpHO2(nnz),&
            &atmpHO(nnz))

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

    ALLOCATE(first_order(nFirst_order,2))
    ALLOCATE(second_order(nSecond_order,2))
    ALLOCATE(higher_order(nHigher_order,2))
    ALLOCATE(ahigher_order(nHigher_order))

    first_order(:,1)  = tmpFO1; first_order(:,2)  = tmpFO2
    second_order(:,1) = tmpSO1; second_order(:,2) = tmpSO2
    higher_order(:,1) = tmpHO1; higher_order(:,2) = tmpHO2
    ahigher_order     = atmpHO

  END SUBROUTINE GatherSpeciesOrder


  !
  ! Read Chemical Data (initial values and Emisions)
  SUBROUTINE InputChemicalData(InitFileName,DataFileName,MeteoFileName)
    CHARACTER(*) :: InitFileName, DataFileName, MeteoFileName
     !
    !REAL(dp), ALLOCATABLE :: GasInAct(:), AqInAct(:)
    !
    INTEGER, PARAMETER :: GasRateUnit=0 ! ???
    INTEGER :: jt, i, iPos
    REAL(dp) :: pi43, LWC
    CHARACTER(60) :: string = ''
    !
    ! for pH Start
    REAL(dp) :: kappa
    !
    pi43=4.0d0/3.0d0*PI
    !
    ! this is for mass transfer (accom , diffus term)
    ALLOCATE( henry_diff(nspc+ntkat), henry_accom(nspc+ntkat) )
    henry_diff = ZERO;    henry_accom = ZERO
    henry_diff(1:ntGas)=5.0d-6
    henry_accom(1:ntGas)=5.0d-5
    !
    !
    !=========================================================================
    !===  Read and save species names in lists
    !=========================================================================
    !
    !-- Open .chem-file and skip the first 24 lines (head)
    OPEN(UNIT=89,FILE=ADJUSTL(TRIM(ChemFile))//'.chem',STATUS='UNKNOWN')
    REWIND(89)
    DO i=1,24
      READ(89,*)
    END DO
    !---  set indices for pH and water dissociation 
    Hp_ind =0
    OHm_ind=0
    aH2O_ind=0
    Temp_ind=nspc+1
    nDIM=nspc
    !
    !--- allocate space for a list of all species (incl. kat spc)
    ALLOCATE(y_name(nspc+ntkat))
    !
    DO jt=1,ntGas
      READ(89,*)  y_name(jt)
    END DO
    DO jt=1,ntAqua
      READ(89,*)  y_name(ntGas+jt)
      IF ( y_name(ntGas+jt)=='Hp' )   Hp_ind=jt+ntGas
      IF ( y_name(ntGas+jt)=='OHm' ) OHm_ind=jt+ntGas
    END DO
    DO jt=1,ntkat
      READ(89,*)  y_name(ntGas+ntAqua+jt)
      IF ( y_name(ntGas+ntAqua+jt)=='[aH2O]' ) aH2O_ind=jt 
    END DO
    !
    IF (MPI_ID==0.AND.ntAqua>=1) THEN
      IF (hp_ind==0)   WRITE(*,*) 'ReadChem...Warning: Hp  not in mechanism!' 
      IF (ohm_ind==0)  WRITE(*,*) 'ReadChem...Warning: OHm  not in mechanism!' 
      IF (ah2o_ind==0) WRITE(*,*) 'ReadChem...Warning: aH2O  not in mechanism!' 
    END IF
    CLOSE(89)
    !
    DO i=1,ntkat
      IF ( y_name(nspc+i)/='[aH2O]' )  THEN
        henry_diff(nspc+i)=5.d-6
        henry_accom(nspc+i)=5.d-5
      END IF
    END DO
    !
    !=========================================================================
    !===  Split Species Names
    !=========================================================================
    !
    IF (ntGas >=1) ALLOCATE (GasName(ntGas))          ! gas phase species
    IF (ntAqua>=1) ALLOCATE (AquaName(ntAqua))        ! aqueous phase species
    IF (ntkat >=1) ALLOCATE (PassName(ntkat))         ! passive phase species
    !
    FORALL ( jt=1:ntGas ) GasName(jt)  = ADJUSTL(y_name(jt))
    FORALL ( jt=1:ntAqua) AquaName(jt) = ADJUSTL(y_name(ntGas+jt))
    FORALL ( jt=1:ntKat ) PassName(jt) = ADJUSTL(y_name(ntGas+ntAqua+jt))
    !
    !=========================================================================
    !===  Set  Chemical DATA  
    !===  (Molar Mass, Charges, Densities) 
    !=========================================================================
    !
    IF ( ntAqua>=1 ) THEN
      !---  Allocate arrays
      ALLOCATE ( Charge(ntAqua),   &      ! charge of ions
               & SolubInd(ntAqua), &      ! solubility
               & MolMass(ntAqua),  &      ! molar mass of species
               & SpcDens(ntAqua),  &      ! density of species
               & OrgIndex(ntAqua), &      ! carbon atoms
               & CC(ntAqua),       &      ! compound class
               & ActIndex(ntAqua)  )      ! index for calculation of activity coefficient
  
      !---  Set default values
      Charge   = ZERO
      SolubInd = ZERO
      MolMass  = ZERO
      SpcDens  = ONE
      OrgIndex = ZERO
      CC       = '  '
      ActIndex = ZERO
      !
      !--------------------------------------------------------------
      !---  Determine Charges of Ions
      DO jt=1,ntAqua
        string=AquaName(jt)
        i=1
        DO ! loop for cations
          iPos=INDEX(string(i:),'p')
          IF (iPos<=0)  EXIT 
          Charge(jt)=Charge(jt) + 1
          i=i+iPos
        END DO
        i=1
        DO ! loop for anions
          iPos=INDEX(string(i:),'m')
          IF (iPos<=0)  EXIT 
          Charge(jt)=Charge(jt) - 1
          i=i+iPos
        END DO
      END DO
    END IF
    !
    !=========================================================================
    !===  Read Initial Values from file   
    !=========================================================================
    !---------------------------------
    ! Units for Chemistry
    ALLOCATE( InitValAct(nspc) , y_e(nspc) )
    ALLOCATE( InitValKat(ntKat) )
    InitValAct=1.0d-20
    InitValKat=1.0d-20
    y_e=ZERO
    !
    !---  Read gase phase inititals
    !---------------------------------
    CALL Read_EMISS( InitFileName , y_e(:) )
    CALL Read_GASini( InitFileName , InitValAct(1:nspc), InitValKat(:) )
    ntKat = ntKatGas
    !
    !--- Read thermodynamic data,....
    CALL Read_SpeciesData(henry_diff,henry_accom,DataFileName)
    !---------------------------------
    !
    !---  Read/Calculate aqua phase initials
    !---------------------------------
    IF ( ntAqua>=1 ) THEN
      !
      CALL Read_AQUAini( InitValAct(ntGas+1:), InitValKat(:), InitFileName )
      ntKat = ntKat + ntKatAqua
      CALL Read_AFRAC( InitAFrac, InitFileName )
      CALL Read_SPEK( SPEK, InitFileName)
      
      LWC  = pseudoLWC(tAnf)
     
      IF ( ntKatAqua>=1 ) THEN
        DO i=1,ntKat
          IF (y_name(nspc+i)=='[aH2O]') THEN
            InitValKat(i) = InitValKat(i)*LWC*mol2part    ! convert aH2O [mol/L] to [molec/cm3]
          END IF
        END DO
      END IF

      !
      !-----------------------------------------
      ! calculate initial aqueus concentrations 
      !-----------------------------------------
      !InitValAct(ntGAS+1:)=1.0d-16     ! notwendig? (willi 10.5.)
      DO i=1,SIZE(InitAFrac)
        iPos=PositionSpeciesAll(InitAFrac(i)%Species)
        IF (iPos>0) THEN
          !
          InitValAct(iPos) = SPEK(1)%Number * kilo          &  ! [#/m3]
          &                * InitAFrac(i)%Frac1              &  ! [g/g]
          &                * (pi43*(SPEK(1)%Radius)**3)      &  ! [m3]
          &                * SPEK(1)%Density                 &  ! [kg/m3]
          &                / (InitAFrac(i)%MolMass)             ! 1/[kg/mol] 
          InitValAct(iPos) = InitValAct(iPos) * mol2Part        ! *[molec/mol]
        END IF
      END DO
    END IF
    
    !======================================================
    !---  Compute pH value and number of ions
    !---  Initial pH by charge balance 
    IF ( pHSet>=1.AND.ntAqua>0 )  THEN
      Kappa = pHValue(InitValAct(ntGas+1:))
      IF ( Kappa > ZERO )  THEN
        InitValAct(Hp_ind)  = Kappa
      ELSE 
        InitValAct(OHm_ind) = InitValAct(OHm_ind) + InitValAct(Hp_ind) - Kappa
      END IF
    END IF
    !
  END SUBROUTINE InputChemicalData
  !
  !
  !
  SUBROUTINE Read_SpeciesData(y_acc,y_diff,FileName)
    REAL(dp) :: y_acc(:) , y_diff(:) 
    CHARACTER(*) :: FileName
    !
    !
    CHARACTER(100) :: SpeciesName
    INTEGER :: iPos, i
    LOGICAL :: Back=.FALSE.
    REAL(dp) :: mm, alpha, dg, c1
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
      CALL LineFile( Back,                                     &
      &              Start1='BEGIN_DATAGAS',                   &
      &              End   ='END_DATAGAS',                     &
      &              Name1=SpeciesName, R1=mm, R2=alpha, R3=dg )
      IF (Back)   EXIT
      !
      iPos=PositionSpeciesAll(SpeciesName)
      IF (iPos>0) THEN
        IF (alpha==0.0d0 .AND. dg==0.0d0) CYCLE GAS 
        !
        y_acc(iPos)=1.0d-12/(3.0d0*dg)
        nue=SQRT(8.0d+03*8.313d0/Pi/mm)
        y_diff(iPos)=4.0d-06/(3.0d0*alpha*nue)
        !
      END IF
    END DO GAS
    CALL RewindFile
    CALL ClearIniFile
    !
    !-----------------------------------------------------------
    ! --- Aqua Phase thermodynamic data
    !
    IF ( ntAQUA>0 ) THEN
      AQUA: DO 
        i=i+1
        CALL LineFile( Back,                                     &
        &              Start1='BEGIN_DATAQUA',                   &
        &              End   ='END_DATAQUA',                     &
        &              Name1=SpeciesName, R1=mm, R2=alpha, R3=dg ,Name2=c2)
        IF (Back)  EXIT
        !
        iPos=PositionSpeciesAll(SpeciesName)
        IF ( iPos>0 .AND. iPos<ntGAS+ntAQUA) THEN
          IF ( alpha==0.0d0 .AND. dg==0.0d0 ) CYCLE AQUA
          !print*, ' wir sind nie hier'
          CC(iPos-ntGas)=ADJUSTL(TRIM(c2))
          y_acc(iPos)=1.0d-12/(3.0d0*dg)
          nue=SQRT(8.0d+03*8.313d0/Pi/mm)
          y_diff(iPos)=4.0d-06/(3.0d0*alpha*nue)
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

      ALLOCATE(RO2(INT(c1)))
      RO2 = 0 
      i   = 0
      DO
        CALL LineFile( Back, Start1='BEGIN_DATARO2',  &
        &              End='END_DATARO2',             &
        &              Name1=SpeciesName)
     
        IF ( Back ) EXIT
        slash=INDEX(SpeciesName,'_')
        IF ( slash>0 ) THEN
          SpeciesName(slash:slash)='/'
        END IF
        IF (PositionSpeciesAll(SpeciesName)>0) THEN
          i=i+1
          RO2(i) = PositionSpeciesAll(SpeciesName)
        END IF
      END DO
      CALL RewindFile
      CALL ClearIniFile
      
      CALL CompressIntegerArray(RO2);   nRO2 = SIZE(RO2)
     
    END IF
    !
    !stop 'chemsysmod'
    !-----------------------------------------------------------
    ! --- Aqua Phase RO2
    !
    IF ( hasRO2aq ) THEN
      CALL LineFile( Back, Start1='BEGIN_DATARO2aq',    &
      &              End='END_DATARO2aq',               &
      &              Name1=SpeciesName,               &
      &              R1=c1)
      !
      i=0

      ALLOCATE(RO2aq(INT(c1)))
      RO2aq = 0
      
      i=0
      DO
        CALL LineFile( Back, Start1='BEGIN_DATARO2aq',  &
        &              End='END_DATARO2aq',             &
        &              Name1=SpeciesName)
        IF (Back) EXIT
        IF (PositionSpeciesALL(SpeciesName)>0) THEN
          i=i+1
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
  SUBROUTINE Read_EMISS(FileName,Emi)
    CHARACTER(*) :: FileName
    REAL(dp) :: Emi(:)
    !
    INTEGER :: iPos
    REAL(dp) :: c1
    CHARACTER(20) :: SpeciesName
    LOGICAL :: Back=.FALSE.
    !
    ! read emission values
    Emi(:)=ZERO
    CALL OpenIniFile(FileName)
    DO
      CALL LineFile( Back,                    &
      &              Start1='BEGIN_GAS',      &
      &              Start2='BEGIN_EMISS',    &
      &              End   ='END_EMISS',      &
      &              Name1=SpeciesName, R1=c1 )
      IF (Back)   EXIT
      !
      iPos=PositionSpeciesALL(SpeciesName)
      IF (iPos>0) Emi(iPos)=c1
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_EMISS
  !
  !
  SUBROUTINE Read_GASini(FileName,GASact,InAct)
    CHARACTER(*) :: FileName
    !
    REAL(dp) :: GASact(:)
    REAL(dp) :: InAct(:)
    !
    INTEGER :: iPos
    REAL(dp) :: c1
    CHARACTER(20) :: SpeciesName
    LOGICAL :: Back=.FALSE.
    !
    !
    GASact = 1.0d-20 
    InAct  = 1.0d-20
    ! read initial values
    CALL OpenIniFile(FileName)
    DO
      CALL LineFile( Back,                       &
      &              Start1 = 'BEGIN_GAS',       &
      &              Start2 = 'BEGIN_INITIAL',   &
      &              End    = 'END_INITIAL',     &
      &              Name1  = SpeciesName,       &
      &              R1     = c1                 )
      IF (Back)   EXIT
      !
      iPos = PositionSpeciesAll(SpeciesName)
      !print*, ' debug :: gas spc :: ', TRIM(SpeciesName), c1, back, iPos
      IF (iPos>0) THEN
        SpeciesName = ADJUSTL(SpeciesName)
        IF (SpeciesName(1:1)=='['.AND. LEN(TRIM(SpeciesName))<maxLENinActDuct) THEN
          InAct(iPos-nspc) = c1
          ntKatGas = ntKatGas + 1
        ELSE
          !iPos=PositionSpeciesAll(SpeciesName)
          GASact(iPos) = c1
       END IF
      END IF
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_GASini
  !
  !
  SUBROUTINE Read_Diag(DiagSpc,DiagPhase,FileName)
    INTEGER, ALLOCATABLE :: DiagSpc(:)
    CHARACTER(1), ALLOCATABLE :: DiagPhase(:)
    CHARACTER(*) :: FileName
    !
    INTEGER :: iPos, cnt
    CHARACTER(50) :: SpeciesName
    LOGICAL :: Back
    !
    !

    ! Read initial values of aqua spc
    CALL OpenIniFile(FileName)
    cnt=0
    DO
      CALL LineFile( Back,                &
      &              Start1='BEGIN_DIAG', &
      &              End='END_DIAG',      &
      &              Name1=SpeciesName    )
      IF (Back)   EXIT
      !
      IF (ADJUSTL(SpeciesName(1:1))/='#') THEN
        iPos=PositionSpeciesAll(SpeciesName)
        IF (iPos>0) THEN
          cnt=cnt+1
        END IF
      END IF
    END DO
    CALL CloseIniFile
    !
    ALLOCATE(DiagSpc(cnt))
    ALLOCATE(DiagPhase(cnt))
    DiagSpc=0
    DiagPhase='e'
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
        iPos=PositionSpeciesAll(TRIM(ADJUSTL(SpeciesName)))
        IF (iPos>0) THEN
          cnt=cnt+1
          !print*, cnt, SpeciesName, ipos
          DiagSpc(cnt)=iPos
          IF (iPos<=ntGas) THEN
            DiagPhase(cnt)='g'
          ELSE IF (iPos>ntGas.AND.iPos<=ntGas+ntAqua) THEN
            DiagPhase(cnt)='a'
          END IF
        END IF
      END IF
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_Diag
  !
  SUBROUTINE Read_AQUAini(AQUAact,AquaInAct,FileName)
    REAL(dp) :: AQUAact(:)
    REAL(dp) :: AquaInAct(:)
    CHARACTER(*) :: FileName
    !
    INTEGER :: iPos
    REAL(dp) :: c1
    CHARACTER(20) :: SpeciesName
    LOGICAL :: Back
    !
    ! Read initial values of aqua spc
    AQUAact(:)=1.0d-16
    CALL OpenIniFile(FileName)
    DO
      CALL LineFile( Back,                    &
      &              Start1 = 'BEGIN_AQUA',   &
      &              Start2 = 'BEGIN_INITIAL',&
      &              End    = 'END_INITIAL',  &
      &              Name1  = SpeciesName,    &
      &              R1     = c1              )
      IF (Back)   EXIT
      
      iPos = PositionSpecies(SpeciesName)
      SpeciesName = ADJUSTL(SpeciesName)

      !print*, ' debug :: gas spc :: ', TRIM(SpeciesName), c1, back, iPos
      IF (iPos>0) THEN
        IF (SpeciesName(1:1)=='['.AND. LEN(TRIM(SpeciesName))<maxLENinActDuct) THEN
          AquaInAct(iPos)=c1
          ntKatAqua = ntKatAqua + 1
        ELSE
          AquaAct(iPos) = c1
       END IF
      END IF
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_AQUAini
  !
  !
  SUBROUTINE Read_AFRAC(AFrac,FileName) 
    CHARACTER(*) :: FileName
    TYPE(AFRAC_T), ALLOCATABLE :: AFrac(:)
    !
    INTEGER :: cnt
    REAL(dp) :: c1,c2,c3,c4
    CHARACTER(20) :: SpeciesName
    LOGICAL :: Back=.FALSE.
    !
    !
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
    ALLOCATE(AFrac(cnt)) 
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
      AFrac(cnt)%Species=SpeciesName
      AFrac(cnt)%MolMass=c1
      AFrac(cnt)%Charge=INT(c2)
      AFrac(cnt)%SolubInd=c3
      AFrac(cnt)%Frac1=c4
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_AFRAC
  !  
  !
  SUBROUTINE Read_SPEK(SPEK,FileName)
    CHARACTER(*) :: FileName
    TYPE(SPEK_T), ALLOCATABLE :: SPEK(:)
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
    ALLOCATE(SPEK(INT(c1))) 
    !
    LWC=pseudoLWC(tAnf)
    cnt=1
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
      IF (c1<1.0d0) THEN
        SPEK(1)%Radius=REAL(c1,KIND=dp)
        SPEK(1)%Number=REAL(c2,KIND=dp)
        SPEK(1)%wetRadius=(3.0d0/4.0d0/PI*LWC/SPEK(1)%Number)**(1.0d0/3.0d0)
        SPEK(1)%Density=REAL(c3,KIND=dp)
      END IF
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_SPEK
  !
  !
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
    CALL LineFile( Back, Start1='BEGIN_DATARO2',    &
    &              End='END_DATARO2',               &
    &              Name1=ro2d,                      &
    &              R1=c1)
    !
    i=0
    DO
      CALL LineFile( Back, Start1='BEGIN_DATARO2',  &
      &              End='END_DATARO2',             &
      &              Name1=SpeciesName)
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
  !
  SUBROUTINE InsertReaction(List,Line,TypeR)
    TYPE(ListReaction_T) :: List
    CHARACTER(*) :: Line(1:4)
    CHARACTER(*)  :: TypeR
    !
    INTEGER :: PosColon,PosEqual
    CHARACTER(LenLine) :: Left,Right
    TYPE(Reaction_T), POINTER :: Reaction
    !
    IF (ASSOCIATED(List%Start)) THEN
      ALLOCATE(List%End%Next)
      List%End=>List%End%Next
    ELSE
      ALLOCATE(List%End)
      List%Start=>List%End
    END IF
    List%LenList=List%LenList+1
    Reaction=>List%End
    PosColon=Index(Line(1),':')
    Reaction%Type=ADJUSTL(Line(1)(PosColon+1:))
    Reaction%Line1=Line(2)
    Reaction%Line3=TRIM(Line(3))                        ! save 3rd Line of reaction
    PosEqual=Index(Reaction%Line1,' = ')
    IF (PosEqual==0)  STOP 'Trennzeichen Falsch'
    !
    Left=Reaction%Line1(1:PosEqual-1)
    CALL ExtractSpecies( Left, Reaction%Educt,     &
    &                    Reaction%InActEduct,      &
    &                    Reaction%nInActEd)  ! geändert
    Right=Reaction%Line1(PosEqual+3:)
    CALL ExtractSpecies( Right, Reaction%Product,  &
    &                    Reaction%InActProduct,    &
    &                    Reaction%nInActPro)! geändert
    CALL ExtractConstants(Line(3),Reaction%Constants,Reaction%TypeConstant)
    Reaction%Line2=Line(3)
    Reaction%Factor=TRIM(Line(4)(9:)) 
    IF (Reaction%Factor(1:1)/='$') Reaction%Factor=' NONE '
    !
    TypeR=Reaction%TypeConstant
  END SUBROUTINE InsertReaction
  !
  !
  SUBROUTINE ReadUnits
     !
     INTEGER :: iLine,Pos
     CHARACTER(LenLine) :: LocLine
     CHARACTER(LenLine) :: Line(3)
     !
     iLine=0
     BACKSPACE(InputUnit)
     DO
       READ(InputUnit,'(a200)',END=1) LocLine
       IF (ADJUSTL(LocLine(1:1))/='#'.AND.         &
       &   ADJUSTL(LocLine(1:7))/='COMMENT'.AND.   &
       &   LEN(TRIM(LocLine))>0) THEN
         iLine=iLine+1
         Line(iLine)=LocLine
         IF (iLine==3) THEN
           Pos=INDEX(Line(1),'GAS')
           IF (Pos>0) THEN
             Line(1)=Line(1)(Pos+3:)
             READ(Line(1),*) UnitGas
           END IF
           Pos=INDEX(Line(2),'AQUA')
           IF (Pos>0) THEN
             Line(2)=Line(2)(Pos+4:)
             READ(Line(2),*) UnitAqua
           END IF
           EXIT
         END IF
       END IF
     END DO
  1  CONTINUE
  END SUBROUTINE ReadUnits
  !
  !
  SUBROUTINE ExtractConstants(String,Constants,Type)
    CHARACTER(*) :: String
    REAL(dp), POINTER :: Constants(:)
    CHARACTER(*) :: Type
    !
    INTEGER :: NumColon,PosColon,PosName,PosComment
    INTEGER :: i,PosNum1,PosNum2,PosNum3,NumNum,PosElem
    CHARACTER(4)  :: NameNumNum
    CHARACTER(10) :: DummyString
    CHARACTER(LEN(String)) :: LocString
    CHARACTER(LEN(String)) :: LocString1
    CHARACTER(LEN(String)) :: NameConstant
    REAL(dp) :: Dummy
    INTEGER :: is
    !
    LocString=String
    String=''
    PosColon=Index(LocString,':')
    Type=LocString(:PosColon-1)
    LocString=ADJUSTL(LocString(PosColon+1:))
    PosComment=INDEX(LocString,'#')
    IF (PosComment>0) LocString=LocString(:PosComment-1)
    !
    NumColon=0
    LocString1=LocString
    IF (Type/='SPECIAL') THEN
      DO
        PosColon=Index(LocString1,':')
        IF (PosColon>0) THEN
          LocString1=ADJUSTL(LocString1(PosColon+1:))
          READ(LocString1,*,IOSTAT=is) Dummy,DummyString
          NumColon=NumColon+1
          PosName=INDEX(LocString1,TRIM(DummyString))
          IF (PosName>0) THEN
            LocString1=LocString1(PosName:)
          END IF
        ELSE
          EXIT
        END IF
      END DO
      ALLOCATE(Constants(NumColon))
      NumColon=0
      DO
        PosColon=Index(LocString,':')
        IF (PosColon>0) THEN
          LocString=ADJUSTL(LocString(PosColon+1:))
          NumColon=NumColon+1
          READ(LocString,*,IOSTAT=is) Constants(NumColon),DummyString
          PosName=INDEX(LocString,TRIM(DummyString))
          IF (PosName>0) LocString=LocString(PosName:)
          !
        ELSE
          EXIT
        END IF
      END DO
    ELSE
      NumNum=0
      DO
        PosNum1=SCAN(LocString1,'1234567890.')
        PosNum2=LEN(LocString1(MAX(PosNum1,1):))
        DO i=1,SIZE(Elements)
          PosElem=INDEX(LocString1(MAX(PosNum1,1):),TRIM(Elements(i)%Element))
          IF (PosElem>0) THEN
            PosNum2=MIN(PosNum2,PosElem-1)
          END IF
        END DO
        PosNum2=PosNum2+PosNum1-1
        IF (PosNum2>0) THEN
          IF ( (LocString1(PosNum2:PosNum2)=='e'.OR.      &
          &     LocString1(PosNum2:PosNum2)=='E').AND.    &
          &    (LocString1(PosNum2+1:PosNum2+1)=='-'.OR.  &
          &     LocString1(PosNum2+1:PosNum2+1)=='+')) THEN
            PosNum3=PosNum2+2
            PosNum2=LEN(LocString(MAX(PosNum1,1):))
            DO i=1,SIZE(Elements)
              PosElem=INDEX(LocString1(MAX(PosNum3,1):),TRIM(Elements(i)%Element))
              IF (PosElem>0) PosNum2=MIN(PosNum2,PosElem-1)
            END DO
            PosNum2=PosNum2+PosNum3-1
          ELSE
            PosNum2=PosNum2-PosNum1+1
          END IF
        ELSE
          PosNum2=PosNum2-PosNum1+1
        END IF
        PosNum2=PosNum2+PosNum1-1
        IF (PosNum1>0) THEN
          LocString1=LocString1(PosNum2+1:)
          NumNum=NumNum+1
        ELSE
          EXIT
        END IF
      END DO
      ALLOCATE(Constants(NumNum))
      NumNum=0
      DO
        PosNum1=SCAN(LocString,'1234567890.')
        PosNum2=LEN(LocString(MAX(PosNum1,1):))
        DO i=1,SIZE(Elements)
          PosElem=INDEX(LocString(MAX(PosNum1,1):),TRIM(Elements(i)%Element))
          IF (PosElem>0) PosNum2=MIN(PosNum2,PosElem-1)
        END DO
        PosNum2=PosNum2+PosNum1-1
        IF (PosNum2>0) THEN
          IF ( (LocString(PosNum2:PosNum2)=='e'.OR.      &
          &     LocString(PosNum2:PosNum2)=='E').AND.    &
          &    (LocString(PosNum2+1:PosNum2+1)=='-'.OR.  &
          &     LocString(PosNum2+1:PosNum2+1)=='+')) THEN
            PosNum3=PosNum2+2
            PosNum2=LEN(LocString(MAX(PosNum1,1):))
            DO i=1,SIZE(Elements)
              PosElem=INDEX(LocString(MAX(PosNum3,1):),TRIM(Elements(i)%Element))
              IF (PosElem>0) PosNum2=MIN(PosNum2,PosElem-1)
            END DO
            PosNum2=PosNum2+PosNum3-1
          ELSE
            PosNum2=PosNum2-PosNum1+1
          END IF
        ELSE
          PosNum2=PosNum2-PosNum1+1
        END IF
        PosNum2=PosNum2+PosNum1-1
        IF (PosNum1>0) THEN
          NameConstant=LocString(PosNum1:PosNum2)
          NumNum=NumNum+1
          READ(NameConstant,*) Constants(NumNum)
          WRITE(NameNumNum,'(I2)') NumNum
          String=TRIM(String)//LocString(:PosNum1-1)//'$'//TRIM(ADJUSTL(NameNumNum))
          LocString=LocString(PosNum2+1:)
        ELSE
          String=TRIM(String)//TRIM(LocString)
          EXIT
        END IF
      END DO
    END IF
  END SUBROUTINE ExtractConstants
  !
  !
  SUBROUTINE ExtractSpecies(String,Duct,InActDuct,NumInActDucts)
    CHARACTER(*) :: String
    TYPE(Duct_T), POINTER :: Duct(:)
    TYPE(Duct_T), POINTER :: InActDuct(:)
    INTEGER :: NumInActDucts
    !
    INTEGER :: PosMinus,PosPlus,NumSpec,PosSpecies
    REAL(dp) :: PreFac !NumberSpecies
    CHARACTER(LenLine) :: Species
    CHARACTER(LEN(String)) :: LocString
    INTEGER :: sbL, sbR
    !
    LocString=String
    NumSpec=1
    DO
     PosPlus=INDEX(LocString,' + ')
     PosMinus=INDEX(LocString,' - ')
      IF (PosMinus==0.OR.PosMinus>PosPlus) THEN
        IF      (PosPlus>0) THEN
          LocString=LocString(PosPlus+3:)
          NumSpec=NumSpec+1
        END IF
      END IF
      IF (PosPlus==0.OR.PosPlus>PosMinus) THEN
        IF (PosMinus>0) THEN
          LocString=LocString(PosMinus+3:)
          NumSpec=NumSpec+1
        END IF
      END IF
      IF (PosPlus==0.AND.PosMinus==0) EXIT
    END DO
    !
    !
    ALLOCATE(Duct(NumSpec))
    ALLOCATE(InActDuct(NumSpec))
    LocString=String
    NumSpec=0
    DO
      PosPlus =INDEX(LocString,' + ')
      PosMinus=INDEX(LocString,' - ')
      IF (PosMinus>0.AND.PosMinus<PosPlus) THEN
        PreFac=-1.0d0
      ELSE
        PreFac=1.0d0
      END IF
      IF (PosPlus>0.AND.(PosPlus<PosMinus.OR.PosMinus==0)) THEN
        Species=ADJUSTL(LocString(:PosPlus-1))
        LocString=LocString(PosPlus+3:)
      ELSE IF (PosMinus>0.AND.(PosMinus<PosPlus.OR.PosPlus==0)) THEN
        Species=ADJUSTL(LocString(:PosMinus-1))
        LocString=LocString(PosMinus+3:)
      ELSE
        Species=ADJUSTL(LocString)
      END IF
      !
      ! check syntax for missing white space 
      !e.g.:  NO2 = O3PX +NO -->  missing white space between + and NO
      IF (INDEX(TRIM(Species),' +')>0.OR.INDEX(TRIM(Species),' -')>0) THEN
        WRITE(*,*) 'Missing white space: ', TRIM(Species)
        WRITE(*,*) 'Check syntax in .sys!'
        CALL FinishMPI()
        STOP 
      END IF
      !
      PosSpecies=SCAN(Species,SetSpecies)
      NumSpec=NumSpec+1
      IF (PosSpecies==1) THEN           
        sbL = INDEX(TRIM(Species),'[')
        sbR = INDEX(TRIM(Species),']')
        IF (sbL==1 .AND. LEN(TRIM(Species))==sbR-sbL+1) THEN
          ! works if there's just one InActEduct
          InActDuct(1)%Koeff=PreFac
          InActDuct(1)%Species=Species
          NumInActDucts=NumInActDucts+1
        END IF
        Duct(NumSpec)%Koeff=PreFac
        Duct(NumSpec)%Species=Species
      ELSE
        !print*, 'spc(1:pos)::: ', Species
        READ(Species(1:PosSpecies-1),*) Duct(NumSpec)%Koeff
        Duct(NumSpec)%Koeff=PreFac*Duct(NumSpec)%Koeff
        Duct(NumSpec)%Species=Species(PosSpecies:)
      END IF
      CALL InsertSpecies(Duct(NumSpec)%Species,Duct(NumSpec)%Type)
      IF (PosPlus==0.AND.PosMinus==0) EXIT
    END DO
  END SUBROUTINE ExtractSpecies
  !
  !
  SUBROUTINE InsertSpecies(Species,Type)
    CHARACTER(*) :: Species
    CHARACTER(*) :: Type
    !
    IF (Species(1:1)=='p') THEN
      CALL InsertHash(ListPartic,TRIM(ADJUSTL(Species)),ntPart)
      Type='Partic'
    ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
      CALL InsertHash(ListAqua,TRIM(ADJUSTL(Species)),ntAQUA)
      Type='Aqua'
    ELSE IF (Species(1:1)=='s') THEN
      CALL InsertHash(ListSolid,TRIM(ADJUSTL(Species)),ntSOLID)
      Type='Solid'
    ELSE IF (Species(1:1)=='['.AND.LEN(TRIM(Species))<maxLENinActDuct.AND. &
    &        Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']')       THEN
      CALL InsertHash(ListNonReac,TRIM(ADJUSTL(Species)),ntKAT)
      Type='Inert'
    ELSE IF (Species(1:1)=='(') THEN
    ELSE
      CALL InsertHash(ListGas,TRIM(ADJUSTL(Species)),ntGAS)
      Type='Gas'
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
      &               + ntGAS
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
  !
  FUNCTION PositionSpeciesAll(Species)
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpeciesAll
    !
    ! 
    PositionSpeciesAll=0
    !
    ! Combustion system
    IF ( TempEq ) THEN
      PositionSpeciesAll=-1
      PositionSpeciesAll=GetHash(ListGas,TRIM(ADJUSTL(Species)))
    ELSE
    ! tropospheric system
      IF (Species(1:1)=='p') THEN
        PositionSpeciesAll=GetHash(ListPartic,TRIM(ADJUSTL(Species))) 
        IF (PositionSpeciesAll>0) THEN
          PositionSpeciesAll=PositionSpeciesAll+ntGAS+ntAQUA+ntSOLID         
        END IF
      ! 
      ! AQUA 
      ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
        PositionSpeciesAll=GetHash(ListAqua,TRIM(ADJUSTL(Species)))
        IF (PositionSpeciesAll>0) THEN
          PositionSpeciesAll=PositionSpeciesAll+ntGAS
        END IF
      !
      ! SOLID
      ELSE IF (Species(1:1)=='s') THEN
        PositionSpeciesAll=GetHash(ListSolid,TRIM(ADJUSTL(Species)))      
        IF (PositionSpeciesAll>0) THEN
          PositionSpeciesAll=PositionSpeciesAll+ntGAS+ntAQUA         
        END IF
      !
      ! NonReac
      ELSE IF (Species(1:1)=='['.AND.                                  &
      &        Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']'.AND.&
      &        LEN(TRIM(Species))<maxLENinActDuct) THEN
        PositionSpeciesAll=GetHash(ListNonReac,TRIM(ADJUSTL(Species)))    
        IF (PositionSpeciesAll>0) THEN
          !PositionSpeciesAll=PositionSpeciesAll+ntGAS+ntAQUA+ntSOLID+ntPart
          PositionSpeciesAll=PositionSpeciesAll+ntGAS+ntAQUA+ntSOLID+ntPart
        END IF
      ! GAS
      ELSE
        PositionSpeciesAll=GetHash(ListGas,TRIM(ADJUSTL(Species)))
      END IF
    END IF
  END FUNCTION PositionSpeciesAll
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
          IF (ASSOCIATED(Current%Constants)) THEN
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
    !
    CALL InitHashTable(ListAqua,100)
    CALL InitHashTable(ListGas,100)
    CALL InitHashTable(ListSolid,100)
    CALL InitHashTable(ListPartic,100)
    CALL InitHashTable(ListNonReac,100)
    CALL OpenFile(FileName,'spc')
    DO
      CALL ReadSpecies(Out)
      IF (Out) EXIT
    END DO
    CALL CloseFile(FileName,'spc')
    CALL OpenFile(FileName,'sys')
    CALL ReadUnits
    DO
      CALL ReadReaction(Out)
      IF (Out) EXIT
    END DO
    CALL CloseFile(FileName,'sys')
    ALLOCATE(ListGas2(ntGAS))
    CALL HashTableToList(ListGas,ListGas2)
    CALL SortList(ListGas2)
    CALL ListToHashTable(ListGas2,ListGas)
    ALLOCATE(ListAqua2(ntAQUA))
    CALL HashTableToList(ListAqua,ListAqua2)
    CALL SortList(ListAqua2)
    CALL ListToHashTable(ListAqua2,ListAqua)
    ALLOCATE(ListSolid2(ntSOLID))
    CALL HashTableToList(ListSolid,ListSolid2)
    CALL SortList(ListSolid2)
    CALL ListToHashTable(ListSolid2,ListSolid)
    ALLOCATE(ListPartic2(ntPart))
    CALL HashTableToList(ListPartic,ListPartic2)
    CALL SortList(ListPartic2)
    CALL ListToHashTable(ListPartic2,ListPartic)
    ALLOCATE(ListNonReac2(ntKAT))
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
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStruct(:)
    TYPE(ListReaction_T) , OPTIONAL :: LGas, LHenry, LAqua, LDiss,                   &
    &                                  LSolid, LPartic, LMicro
    INTEGER :: i, j, iList, iEq
    INTEGER :: nList

    INTEGER :: icnt(47), icntFAC(10), iHen
            
    !
    ! #Reactions
    neq= nreakgas     &
    &  +2*nreakhenry &
    &  +nreakaqua    &
    &  +2*nreakdissoc  
    !
    ! #Spezies berechnen
    nspc= ntGAS+ntAQUA    
    !
    !
    !
    IF (PRESENT(LGas))    nList=nList+1
    IF (PRESENT(LHenry))  nList=nList+1
    IF (PRESENT(LAqua))   nList=nList+1
    IF (PRESENT(LDiss))   nList=nList+1
    IF (PRESENT(LSolid))  nList=nList+1
    IF (PRESENT(LPartic)) nList=nList+1
    IF (PRESENT(LMicro))  nList=nList+1
    ALLOCATE(CompleteReactionList(nList))
    nList=0
    IF (PRESENT(LGas)) THEN
      nList=nList+1
      CompleteReactionList(nList)=LGas
    END IF
    IF (PRESENT(LHenry)) THEN
      nList=nList+1
      CompleteReactionList(nList)=LHenry
    END IF
    IF (PRESENT(LAqua)) THEN 
      nList=nList+1
      CompleteReactionList(nList)=LAqua
    END IF
    IF (PRESENT(LDiss)) THEN 
      nList=nList+1
      CompleteReactionList(nList)=LDiss
    END IF
    IF (PRESENT(LSolid)) THEN
      nList=nList+1
      CompleteReactionList(nList)=LSolid
    END IF
    IF (PRESENT(LPartic)) THEN
      nList=nList+1
      CompleteReactionList(nList)=LPartic
    END IF
    IF (PRESENT(LMicro)) THEN
      nList=nList+1
      CompleteReactionList(nList)=LMicro
    END IF

    ALLOCATE( ReacStruct(neq) )
    !
    CALL AllocateRTarrays
    !
    i=1
    iHen = 0
    icnt = 0
    icntFAC = 0
    !
    DO iList=1,nList
      Current=>CompleteReactionList(iList)%Start
      DO WHILE (ASSOCIATED(Current)) 
        ReacStruct(i)%Type   = Current%Type
        ReacStruct(i)%Line1  = Current%Line1
        ReacStruct(i)%Line2  = Current%Line2
        ReacStruct(i)%Line3  = Current%Line3
        ReacStruct(i)%Factor = Current%Factor
        ReacStruct(i)%TypeConstant = Current%TypeConstant
         
        CALL GatherReacFACTOR(i,icntFAC,Current%Factor)

        ! forward direaction
        CALL GatherReacArrays( i , icnt ,           &
                             & Current%Type ,       &
                             & Current%TypeConstant,&
                             & Current%Constants    )
        !
        ReacStruct(i)%nActEd  = SIZE(Current%Educt)
        ReacStruct(i)%nActPro = SIZE(Current%Product)
        ALLOCATE( ReacStruct(i)%Educt(ReacStruct(i)%nActEd),   &
                & ReacStruct(i)%Product(ReacStruct(i)%nActPro) )

        DO j = 1 , ReacStruct(i)%nActEd
          ReacStruct(i)%Educt(j)%Species  = Current%Educt(j)%Species
          ReacStruct(i)%Educt(j)%Type     = Current%Educt(j)%Type
          ReacStruct(i)%Educt(j)%Koeff    = Current%Educt(j)%Koeff
          ReacStruct(i)%Educt(j)%iSpecies = PositionSpeciesAll(Current%Educt(j)%Species)
        END DO
        DO j = 1 , ReacStruct(i)%nActPro
          ReacStruct(i)%Product(j)%Species  = Current%Product(j)%Species
          ReacStruct(i)%Product(j)%Type     = Current%Product(j)%Type
          ReacStruct(i)%Product(j)%Koeff    = Current%Product(j)%Koeff
          ReacStruct(i)%Product(j)%iSpecies = PositionSpeciesAll(Current%Product(j)%Species)
        END DO
        
        ReacStruct(i)%nConst  = SIZE(Current%Constants)
        ALLOCATE(ReacStruct(i)%Constants(ReacStruct(i)%nConst))
        ReacStruct(i)%Constants = Current%Constants
        !
        IF ( Current%Type == 'HENRY' ) THEN
          ReacStruct(i)%direction = 'GA'
          ReacStruct(i)%HenrySpc  = PositionSpeciesAll(ReacStruct(i)%Educt(1)%Species)

          icnt(29) = icnt(29) + 1
          RTind2%iHENRY(icnt(29),1) = i
          RTind2%iHENRY(icnt(29),2) = PositionSpeciesAll(ReacStruct(i)%Educt(1)%Species)
          RTind2%iHENRY(icnt(29),3) = i + 1
          RTind2%iHENRY(icnt(29),4) = PositionSpeciesAll(ReacStruct(i)%Product(1)%Species)
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
          IF ( ReacStruct(i)%SumAqCoef > ZERO ) nHOaqua = nHOaqua + 1
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
            CALL GatherReacFACTOR(i,icntFAC,Current%Factor)
           
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
            
            ReacStruct(i)%nConst = SIZE(Current%Constants)
            ALLOCATE(ReacStruct(i)%Constants(ReacStruct(i)%nConst))
            ReacStruct(i)%Constants = Current%Constants
            !
            IF ( Current%Type == 'HENRY' ) THEN
              ReacStruct(i)%direction = 'AG'
              ReacStruct(i)%HenrySpc  = PositionSpeciesAll(ReacStruct(i)%Product(1)%Species)
            END IF
            !
            ReacStruct(i)%nInActEd  = Current%nInActPro
            ReacStruct(i)%nInActPro = Current%nInActEd
            ALLOCATE( ReacStruct(i)%InActEduct(Current%nInActPro),   &
                    & ReacStruct(i)%InActEductSpc(Current%nInActPro),& 
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
              IF ( ReacStruct(i)%SumAqCoef > ZERO ) nHOaqua = nHOaqua + 1
            END IF
            !
        END SELECT
        !
        Current=>Current%Next
        i=i+1
      END DO
    END DO
    
    ! build the array for mass action products of inactive species (katalytic?)
    ALLOCATE( first_orderKAT(nfirst_orderKAT,2) )
    nFirst_orderKAT = 0

    ! counting the aquatic reactions with more than one educt
    ALLOCATE( RTind2%iHOaqua(nHOaqua), RTpar2%HOaqua(nHOaqua) )
    nHOaqua = 0

    DO i=1,neq
      
      IF ( ReacStruct(i)%nInActEd > 0 ) THEN
        nFirst_orderKAT = nFirst_orderKAT + 1
        first_orderKAT(nfirst_orderKAT,1) = i
        first_orderKAT(nfirst_orderKAT,2) = PositionSpeciesAll(ReacStruct(i)%InActEductSpc(1)) - nspc
      END IF
      
      IF ( ReacStruct(i)%Type=='AQUA'.OR. ReacStruct(i)%Type=='DISS' ) THEN
        IF ( ReacStruct(i)%SumAqCoef > ZERO ) THEN
          nHOaqua = nHOaqua + 1
          RTind2%iHOaqua(nHOaqua) = i
          RTpar2%HOaqua(nHOaqua)  = ReacStruct(i)%SumAqCoef
        END IF
      END IF

    END DO
    
    !print*, ' nHenry = ',nHenry, nHenryga, nHenryag
    !stop
  END SUBROUTINE AllListsToArray
  !
  SUBROUTINE GatherReacFACTOR(iR,icntFAC,Factor)
    CHARACTER(*), INTENT(IN) :: Factor
    INTEGER,      INTENT(IN) :: iR
    INTEGER,      INTENT(INOUT) :: icntFAC(:)

    SELECT CASE (Factor)
      CASE ('$H2')
        icntFAC(1) = icntFAC(1) + 1
        RTind2%iFAC_H2(icntFAC(1)) = iR
      CASE ('$O2N2')
        icntFAC(2) = icntFAC(2) + 1
        RTind2%iFAC_O2N2(icntFAC(2)) = iR
      CASE ('$M')
        icntFAC(3) = icntFAC(3) + 1
        RTind2%iFAC_M(icntFAC(3)) = iR
      CASE ('$O2')
        icntFAC(4) = icntFAC(4) + 1
        RTind2%iFAC_O2(icntFAC(4)) = iR
      CASE ('$N2')
        icntFAC(5) = icntFAC(5) + 1
        RTind2%iFAC_N2(icntFAC(5)) = iR
      CASE ('$H2O')
        icntFAC(6) = icntFAC(6) + 1
        RTind2%iFAC_H2O(icntFAC(6)) = iR
      CASE ('$RO2')
        icntFAC(7) = icntFAC(7) + 1
        RTind2%iFAC_RO2(icntFAC(7)) = iR
        hasRO2  = .TRUE.
      CASE ('$O2O2')
        icntFAC(8) = icntFAC(8) + 1
        RTind2%iFAC_O2O2(icntFAC(8)) = iR
      CASE ('$aH2O')
        icntFAC(9) = icntFAC(9) + 1
        RTind2%iFAC_aH2O(icntFAC(9)) = iR
      CASE ('$RO2aq')
        icntFAC(10) = icntFAC(10) + 1
        RTind2%iFAC_RO2aq(icntFAC(10)) = iR
        hasRO2aq = .TRUE.
    END SELECT
  END SUBROUTINE GatherReacFACTOR

  SUBROUTINE GatherReacArrays(iR,icnt,Typ,TypeR,C)
    REAL(dp), INTENT(IN) :: C(:)
    CHARACTER(*),   INTENT(IN) :: Typ
    CHARACTER(*),   INTENT(IN) :: TypeR
    INTEGER,        INTENT(IN) :: iR
    INTEGER,        INTENT(INOUT) :: icnt(47)

    SELECT CASE ( TypeR )
      CASE ('CONST')
        icnt(1) = icnt(1) + 1
        RTind2%iCONST(icnt(1))   = iR
        RTpar2%CONST(icnt(1)) = C(1)
      CASE ('PHOTAB') 
        icnt(2) = icnt(2) + 1
        RTind2%iPHOTab(icnt(2)) = iR
        RTpar2%PHOTab(icnt(2),:) = C
      CASE ('PHOTABC')
        icnt(3) = icnt(3) + 1
        RTind2%iPHOTabc(icnt(3)) = iR
        RTpar2%PHOTabc(icnt(3),:) = C 
      CASE ('PHOTMCM') 
        icnt(4) = icnt(4) + 1
        RTind2%iPHOTmcm(icnt(4)) = iR 
        RTpar2%PHOTmcm(icnt(4),:) = C 
      CASE ('TEMP1') 
        icnt(5) = icnt(5) + 1
        RTind2%iTEMP1(icnt(5)) = iR 
        RTpar2%TEMP1(icnt(5),:) = C 
      CASE ('TEMP2') 
        icnt(6) = icnt(6) + 1
        RTind2%iTEMP2(icnt(6)) = iR 
        RTpar2%TEMP2(icnt(6),:) = C 
      CASE ('TEMP3')
        icnt(7) = icnt(7) + 1
        RTind2%iTEMP3(icnt(7)) = iR 
        RTpar2%TEMP3(icnt(7),:) = C 
      CASE ('TEMP4') 
        icnt(8) = icnt(8) + 1
        RTind2%iTEMP4(icnt(8)) = iR 
        RTpar2%TEMP4(icnt(8),:) = C 
      CASE ('TROE')  
        icnt(9) = icnt(9) + 1
        RTind2%iTROE(icnt(9)) = iR 
        RTpar2%TROE(icnt(9),:) = C 
      CASE ('TROEF') 
        icnt(10) = icnt(10) + 1
        RTind2%iTROEf(icnt(10)) = iR  
        RTpar2%TROEf(icnt(10),:) = C  
      CASE ('TROEQ') 
        icnt(11) = icnt(11) + 1
        RTind2%iTROEq(icnt(11)) = iR 
        RTpar2%TROEq(icnt(11),:) = C 
      CASE ('TROEQF')
        icnt(12) = icnt(12) + 1
        RTind2%iTROEqf(icnt(12)) = iR 
        RTpar2%TROEqf(icnt(12),:) = C 
      CASE ('TROEXP')
        icnt(13) = icnt(13) + 1
        RTind2%iTROExp(icnt(13)) = iR 
        RTpar2%TROExp(icnt(13),:) = C 
      CASE ('TROEMCM') 
        icnt(14) = icnt(14) + 1
        RTind2%iTROEmcm(icnt(14)) = iR 
        RTpar2%TROEmcm(icnt(14),:) = C 
      CASE ('SPEC1') 
        icnt(15) = icnt(15) + 1
        RTind2%iSPEC1(icnt(15)) = iR 
        RTpar2%SPEC1(icnt(15),:) = C 
      CASE ('SPEC2') 
        icnt(16) = icnt(16) + 1
        RTind2%iSPEC2(icnt(16)) = iR 
        RTpar2%SPEC2(icnt(16),:) = C 
      CASE ('SPEC3') 
        icnt(17) = icnt(17) + 1
        RTind2%iSPEC3(icnt(17)) = iR 
        RTpar2%SPEC3(icnt(17),:) = C 
      CASE ('SPEC4') 
        icnt(18) = icnt(18) + 1
        RTind2%iSPEC4(icnt(18)) = iR 
        RTpar2%SPEC4(icnt(18),:) = C 
      CASE ('SPEC1MCM')
        icnt(19) = icnt(19) + 1
        RTind2%iSPEC1mcm(icnt(19)) = iR 
        RTpar2%SPEC1mcm(icnt(19),:) = C 
      CASE ('SPEC2MCM') 
        icnt(20) = icnt(20) + 1
        RTind2%iSPEC2mcm(icnt(20)) = iR 
        RTpar2%SPEC2mcm(icnt(20),:) = C 
      CASE ('SPEC3MCM') 
        icnt(21) = icnt(21) + 1
        RTind2%iSPEC3mcm(icnt(21)) = iR 
        RTpar2%SPEC3mcm(icnt(21),:) = C 
      CASE ('SPEC4MCM')
        icnt(22) = icnt(22) + 1
        RTind2%iSPEC4mcm(icnt(22)) = iR 
        RTpar2%SPEC4mcm(icnt(22),:) = C 
      CASE ('SPEC5MCM') 
        icnt(23) = icnt(23) + 1
        RTind2%iSPEC5mcm(icnt(23)) = iR 
        RTpar2%SPEC5mcm(icnt(23),:) = C 
      CASE ('SPEC6MCM') 
        icnt(24) = icnt(24) + 1
        RTind2%iSPEC6mcm(icnt(24)) = iR 
        RTpar2%SPEC6mcm(icnt(24),:) = C 
      CASE ('SPEC7MCM') 
        icnt(25) = icnt(25) + 1
        RTind2%iSPEC7mcm(icnt(25)) = iR 
        RTpar2%SPEC7mcm(icnt(25),:) = C 
      CASE ('SPEC8MCM') 
        icnt(26) = icnt(26) + 1
        RTind2%iSPEC8mcm(icnt(26)) = iR 
        RTpar2%SPEC8mcm(icnt(26),:) = C 
      CASE ('S4H2O') 
        icnt(27) = icnt(27) + 1
        RTind2%iS4H2O(icnt(27)) = iR 
        RTpar2%S4H2O(icnt(27),:) = C 
      CASE ('T1H2O') 
        icnt(28) = icnt(28) + 1
        RTind2%iT1H2O(icnt(28)) = iR 
        RTpar2%T1H2O(icnt(28),:) = C 
      CASE ('ASPEC1') 
        icnt(30) = icnt(30) + 1
        RTind2%iASPEC1(icnt(30)) = iR 
        RTpar2%ASPEC1(icnt(30),:) = C 
      CASE ('ASPEC2') 
        icnt(31) = icnt(31) + 1
        RTind2%iASPEC2(icnt(31)) = iR 
        RTpar2%ASPEC2(icnt(31),:) = C 
      CASE ('ASPEC3') 
        icnt(32) = icnt(32) + 1
        RTind2%iASPEC3(icnt(32)) = iR 
        RTpar2%ASPEC3(icnt(32),:) = C 
      CASE ('ASPEC4') 
        icnt(33) = icnt(33) + 1
        RTind2%iASPEC4(icnt(33)) = iR 
        RTpar2%ASPEC4(icnt(33),:) = C(1)
      CASE ('DCONST') 
        icnt(34) = icnt(34) + 1
        RTind2%iDCONST(icnt(34),1) = iR
        RTind2%iDCONST(icnt(34),2) = iR + 1
        RTpar2%DCONST(icnt(34),:)  = C 
      CASE ('DTEMP')  
        icnt(35) = icnt(35) + 1
        RTind2%iDTEMP(icnt(35),1) = iR
        RTind2%iDTEMP(icnt(35),2) = iR + 1
        RTpar2%DTEMP(icnt(35),:)  = C 
      CASE ('DTEMP2') 
        icnt(36) = icnt(36) + 1
        RTind2%iDTEMP2(icnt(36),1)  = iR
        RTind2%iDTEMP2(icnt(36),2)  = iR + 1 
        RTpar2%DTEMP2(icnt(36),:) = C 
      CASE ('DTEMP3') 
        icnt(37) = icnt(37) + 1
        RTind2%iDTEMP3(icnt(37),1) = iR
        RTind2%iDTEMP3(icnt(37),2) = iR + 1
        RTpar2%DTEMP3(icnt(37),:)  = C 
      CASE ('DTEMP4') 
        icnt(38) = icnt(38) + 1
        RTind2%iDTEMP4(icnt(38),1) = iR
        RTind2%iDTEMP4(icnt(38),2) = iR + 1
        RTpar2%DTEMP4(icnt(38),:)  = C 
      CASE ('DTEMP5') 
        icnt(39) = icnt(39) + 1
        RTind2%iDTEMP5(icnt(39),1) = iR
        RTind2%iDTEMP5(icnt(39),2) = iR + 1
        RTpar2%DTEMP5(icnt(39),:)  = C 
      CASE ('MESKHIDZE') 
        icnt(40) = icnt(40) + 1
        RTind2%iMeskhidze(icnt(40),1)  = iR
        RTind2%iMeskhidze(icnt(40),2)  = iR + 1
        RTpar2%Meskhidze(icnt(40),:) = C 
      CASE ('PHOTO') 
        icnt(41) = icnt(41) + 1
        RTind2%iPHOTOkpp(icnt(41)) = iR
        RTpar2%PHOTOkpp(icnt(41))  = C(1)
      CASE ('PHOTO2') 
        icnt(42) = icnt(42) + 1
        RTind2%iPHOTO2kpp(icnt(42)) = iR
        RTpar2%PHOTO2kpp(icnt(42))  = C(1)
      CASE ('PHOTO3') 
        icnt(43) = icnt(43) + 1
        RTind2%iPHOTO3kpp(icnt(43)) = iR
        RTpar2%PHOTO3kpp(icnt(43))  = C(1)
      CASE DEFAULT
        WRITE(*,*) ''
        WRITE(*,*) ' Reaction Type unknown:  ',TypeR,'  --> check input file'
        WRITE(*,*) ''
    END SELECT

  END SUBROUTINE GatherReacArrays
  !
  !
  SUBROUTINE CheckConstants(RS)
    TYPE(ReactionStruct_T) :: RS(:)     ! reaction system
    CHARACTER(15) :: Const_T            ! constant type for reaction i
    !
    INTEGER :: i,j
    !

    DO i=1,SIZE(RS)
      Const_T=RS(i)%TypeConstant
      DO j=1,nr_reac 
        IF (TRIM(var_par(j)%str_type)==ADJUSTL(TRIM(Const_T))) EXIT
      END DO
      IF (SIZE(RS(i)%Constants)/=var_par(j)%anzp) THEN
        WRITE(*,*) 'ERROR: Wrong number of constants in reaction: ',i  
        WRITE(*,*) '----->  reaction:     ',i, '   ', TRIM(RS(i)%Line1)
        WRITE(*,*) '----->  soll #consts: ', var_par(j)%anzp, j
        WRITE(*,*) '----->  ist  #consts: ', SIZE(RS(i)%Constants)
        WRITE(*,*) '       Check sys-file for syntax errors!'
        CALL FinishMPI()
        STOP 'STOP'
      END IF
    END DO
  END SUBROUTINE CheckConstants

  SUBROUTINE CompressIntegerArray(Array)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: Array(:)
    INTEGER, ALLOCATABLE :: tmpArray(:)
    
    INTEGER :: i, N, cnt, M
    
    N = COUNT(Array/=0)
    ALLOCATE(tmpArray(N))

    cnt = 0
    DO i=1,SIZE(Array)
      IF (Array(i)/=0) THEN
        cnt = cnt + 1
        !tmpArray(cnt) = i
        tmpArray(cnt) = Array(i)
      END IF
    END DO
    DEALLOCATE(Array)
    ALLOCATE(Array(N))
    Array = tmpArray
  END SUBROUTINE CompressIntegerArray


  SUBROUTINE CompressDoubleArray(Array)
    REAL(dp), ALLOCATABLE, INTENT(INOUT) :: Array(:)
    REAL(dp), ALLOCATABLE :: tmpArray(:)
    
    INTEGER :: i, N, cnt
    REAL(dp), PARAMETER :: big = -99999999999999.d0
    
    N = COUNT( Array /= big )
    ALLOCATE(tmpArray(N))

    cnt = 0
    DO i=1,SIZE(Array)
      IF ( Array(i) /= big ) THEN
        cnt = cnt + 1
        tmpArray(cnt) = Array(i)
      END IF
    END DO
    DEALLOCATE(Array)
    ALLOCATE(Array(N))
    Array = tmpArray
  END SUBROUTINE CompressDoubleArray

  SUBROUTINE AllocateRTarrays()

    ! allocate index arrays and parameter arrays for the new vectorized version
    ALLOCATE( RTind2%iCONST(nCONST), RTpar2%CONST(nCONST))

    ALLOCATE( RTind2%iPHOTabc(nPHOTabc), RTind2%iPHOTab(nPHOTab), RTind2%iPHOTmcm(nPHOTmcm))
    ALLOCATE( RTpar2%PHOTabc(nPHOTabc,3), RTpar2%PHOTab(nPHOTab,2), RTpar2%PHOTmcm(nPHOTmcm,3))

    ALLOCATE( RTind2%iTEMP1(nTEMP1), RTind2%iTEMP2(nTEMP2), RTind2%iTEMP3(nTEMP3), RTind2%iTEMP4(nTEMP4))
    ALLOCATE( RTpar2%TEMP1(nTEMP1,2), RTpar2%TEMP2(nTEMP2,2), RTpar2%TEMP3(nTEMP3,2), RTpar2%TEMP4(nTEMP4,2))

    ALLOCATE( RTind2%iTROE(nTROE),     RTind2%iTROEf(nTROEf),   RTind2%iTROEq(nTROEq),   &
            & RTind2%iTROEqf(nTROEqf), RTind2%iTROExp(nTROExp), RTind2%iTROEmcm(nTROEmcm))
    ALLOCATE( RTpar2%TROE(nTROE,4),     RTpar2%TROEf(nTROEf,5),   RTpar2%TROEq(nTROEq,6),   &
            & RTpar2%TROEqf(nTROEqf,7), RTpar2%TROExp(nTROExp,5), RTpar2%TROEmcm(nTROEmcm,10))

    ALLOCATE( RTind2%iSPEC1(nSPEC1), RTind2%iSPEC2(nSPEC2), RTind2%iSPEC3(nSPEC3), RTind2%iSPEC4(nSPEC4))
    ALLOCATE( RTpar2%SPEC1(nSPEC1,2), RTpar2%SPEC2(nSPEC2,2), RTpar2%SPEC3(nSPEC3,6), RTpar2%SPEC4(nSPEC4,4))

    ALLOCATE( RTind2%iSPEC1mcm(nSPEC1mcm), RTind2%iSPEC2mcm(nSPEC2mcm),&
            & RTind2%iSPEC3mcm(nSPEC3mcm), RTind2%iSPEC4mcm(nSPEC4mcm),&
            & RTind2%iSPEC5mcm(nSPEC5mcm), RTind2%iSPEC6mcm(nSPEC6mcm),&
            & RTind2%iSPEC7mcm(nSPEC7mcm), RTind2%iSPEC8mcm(nSPEC8mcm) )
    ALLOCATE( RTpar2%SPEC1mcm(nSPEC1mcm,3), RTpar2%SPEC2mcm(nSPEC2mcm,3),&
            & RTpar2%SPEC3mcm(nSPEC3mcm,2), RTpar2%SPEC4mcm(nSPEC4mcm,4),&
            & RTpar2%SPEC5mcm(nSPEC5mcm,4), RTpar2%SPEC6mcm(nSPEC6mcm,4),&
            & RTpar2%SPEC7mcm(nSPEC7mcm,6), RTpar2%SPEC8mcm(nSPEC8mcm,4) )

    ALLOCATE( RTind2%iS4H2O(nS4H2O), RTind2%iT1H2O(nT1H2O))
    ALLOCATE( RTpar2%S4H2O(nS4H2O,4), RTpar2%T1H2O(nT1H2O,2))

    ALLOCATE( RTind2%iASPEC1(nASPEC1), RTind2%iASPEC2(nASPEC2),& 
            & RTind2%iASPEC3(nASPEC3), RTind2%iASPEC4(nASPEC4) )
    ALLOCATE( RTpar2%ASPEC1(nASPEC1,2), RTpar2%ASPEC2(nASPEC2,3),&
            & RTpar2%ASPEC3(nASPEC3,2), RTpar2%ASPEC4(nASPEC4,3) ) 

    ALLOCATE( RTind2%iDTEMP(nDTEMP,2),   RTind2%iDTEMP2(nDTEMP2,2),                           &
            & RTind2%iDTEMP3(nDTEMP3,2), RTind2%iDTEMP4(nDTEMP4,2), RTind2%iDTEMP5(nDTEMP5,2) )
    ALLOCATE( RTpar2%DTEMP(nDTEMP,3),    RTpar2%DTEMP2(nDTEMP2,4),                           &
            & RTpar2%DTEMP3(nDTEMP3,4),  RTpar2%DTEMP4(nDTEMP4,4),  RTpar2%DTEMP5(nDTEMP5,3) )

    ALLOCATE( RTind2%iDCONST(nDCONST,2), RTind2%iMeskhidze(nMeskhidze,2) )
    ALLOCATE( RTpar2%DCONST(nDCONST,2),  RTpar2%Meskhidze(nMeskhidze,7) )
    
    ALLOCATE( RTind2%iFAC_H2(nFAC_H2), RTind2%iFAC_O2N2(nFAC_O2N2), RTind2%iFAC_M(nFAC_M),        &
            & RTind2%iFAC_O2(nFAC_O2), RTind2%iFAC_N2(nFAC_N2), RTind2%iFAC_H2O(nFAC_H2O),        &
            & RTind2%iFAC_RO2(nFAC_RO2), RTind2%iFAC_O2O2(nFAC_O2O2), RTind2%iFAC_aH2O(nFAC_aH2O),&
            & RTind2%iFAC_RO2aq(nFAC_RO2aq) )

    !  dim1 = iReac, dim2= iSpc, dim3 = iR_g->a / iR_a->g
    ALLOCATE( RTind2%iHENRY(nHENRY,4) )

    ALLOCATE( RTind2%iPHOTOkpp(nPHOTOkpp), RTind2%iPHOTO2kpp(nPHOTO2kpp), RTind2%iPHOTO3kpp(nPHOTO3kpp),& 
            & RTpar2%PHOTOkpp(nPHOTOkpp),  RTpar2%PHOTO2kpp(nPHOTO2kpp),  RTpar2%PHOTO3kpp(nPHOTO3kpp)  )
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
END MODULE Chemsys_Mod
