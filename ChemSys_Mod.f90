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
    REAL(RealKind)     :: Koeff
  END TYPE Duct_T
  !
  ! LIST FORM
  TYPE Reaction_T
    CHARACTER(20)       :: Type, TypeConstant
    CHARACTER(LenLine) :: Line1, Line2, Line3
    CHARACTER(LenName) :: Factor
    TYPE(Duct_T)  , POINTER   :: Educt(:)=>NULL(),                  &
    &                            Product(:)=>NULL()
    REAL(RealKind), POINTER   :: Constants(:)=>NULL()
    TYPE(Duct_T)  , POINTER   :: InActEduct(:)=>NULL(),             &
    &                            InActProduct(:)=>NULL()
    INTEGER                   :: nInActEd=0, nInActPro=0
    TYPE(Reaction_T), POINTER :: Next=>NULL()
  END TYPE Reaction_T
  !
  ! ARRAY FORM
  TYPE ReactionStruct_T
    CHARACTER(20)       :: Type,  TypeConstant
    CHARACTER(LenLine)  :: Line1           &
    &                    , Line2=''        &  ! Line2 = BackReaction if nessessary
    &                    , Line3=''
    CHARACTER(LenName)  :: Factor
    CHARACTER(2)        :: direction
    REAL(RealKind)      :: SumAqCoef     
    TYPE(Duct_T)  , ALLOCATABLE     :: Educt(:), Product(:)
    REAL(RealKind), ALLOCATABLE     :: Constants(:)
    REAL(RealKind), ALLOCATABLE     :: LowConst(:), HighConst(:)
    REAL(RealKind), ALLOCATABLE     :: TroeConst(:)
    REAL(RealKind), ALLOCATABLE     :: InActEduct(:), InActProduct(:)
    INTEGER                         :: nInActEd=0, nInActPro=0
    INTEGER                         :: nActEd=0  , nActPro=0
    INTEGER                         :: NumConst=0
    INTEGER                         :: HenrySpc=0
    INTEGER, ALLOCATABLE            :: TB(:)
    REAL(RealKind), ALLOCATABLE     :: TBalpha(:)
    CHARACTER(LenName), ALLOCATABLE :: InActEductSpc(:)
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
    REAL(RealKind)     :: Hf=0.0d0, Gf=0.0d0, Cp=0.0d0
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
    REAL(RealKind) :: MolMass           ! [g/mol]
    INTEGER        :: Charge            ! ladung (+,-,++,--,...)
    REAL(RealKind) :: SolubInd          ! Löslichkeitsindex
    REAL(RealKind) :: Frac1             ! [g/g]
  END TYPE AFRAC_T
  !
  TYPE(AFRAC_T), ALLOCATABLE :: InitAFrac(:)
  !
  !
  TYPE SPEK_T
    REAL(RealKind) :: Radius            ! [m]   radius partivle
    REAL(RealKind) :: wetRadius         ! [m]   radius droplett
    REAL(RealKind) :: Number            ! [#/cm3]
    REAL(RealKind) :: Density           ! [kg/m3]
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
  INTEGER :: NSactNR                        ! # activ species + all Reactions
  !
  INTEGER :: UnitGas=0
  INTEGER :: UnitAqua=0
  !
  CHARACTER(20) :: Filename !='Salt'
  CHARACTER(20) :: IniName  !='Salt'
  !
  REAL(RealKind), PARAMETER :: RGas=8.3145d0
  REAL(RealKind), PARAMETER :: TRef=280.0d0 !298.15d0
  !
  TYPE(Reaction_T), POINTER :: Current
  TYPE(ReactionStruct_T), ALLOCATABLE :: ReactionSystem(:)
  TYPE(ListReaction_T), ALLOCATABLE :: CompleteReactionList(:)
  !
  !
  REAL(RealKind), ALLOCATABLE :: Emis(:)          & ! emission values
  !&                            , InitValAct(:)    & ! initial values activ spc
  &                            , InitValInAct(:)    ! initial values inactiv spc
  !
  !
  CHARACTER(LenName), ALLOCATABLE :: RO2spcG(:) , RO2spcA(:)
  INTEGER, ALLOCATABLE :: RO2idxG(:) , RO2idxA(:)
  !
  !
  REAL(RealKind) :: aH2O
  !
  REAL(RealKind), ALLOCATABLE :: sumBAT(:)         ! sum_j=1,n_s (b_ij-a_ij),  i el. N_R
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
  !========
  SUBROUTINE mkArray(Arr,Typen)
    INTEGER :: Arr(47)
    TYPE(NReacType_T) :: Typen
    !
    !
    Arr(1)=Typen%GasPhoto
    Arr(2)=Typen%GasPhotAB
    Arr(3)=Typen%GasPhotABC
    Arr(4)=Typen%GasPhotMCM
    Arr(5)=Typen%GasConst
    Arr(6)=Typen%Temp
    Arr(7)=Typen%Temp1
    Arr(8)=Typen%Temp2
    Arr(9)=Typen%Temp3
    Arr(10)=Typen%Troe
    Arr(11)=Typen%Troef
    Arr(12)=Typen%TroeQ
    Arr(13)=Typen%Spec1
    Arr(14)=Typen%Spec2
    Arr(15)=Typen%Spec3
    Arr(16)=Typen%Spec4
    Arr(17)=Typen%Spec1MCM
    Arr(18)=Typen%Spec2MCM
    Arr(19)=Typen%Spec3MCM
    Arr(20)=Typen%SPec4MCM
    Arr(21)=Typen%Spec5MCM
    Arr(22)=Typen%Spec6MCM
    Arr(23)=Typen%Spec7MCM
    Arr(24)=Typen%Spec8MCM
    Arr(25)=Typen%S4H2O
    Arr(26)=Typen%Henry
    Arr(27)=Typen%AquaPhoto
    Arr(28)=Typen%AquaPhotAB
    Arr(29)=Typen%AquaPhotABC
    Arr(30)=Typen%AquaPhotMCM 
    Arr(31)=Typen%AquaConst
    Arr(32)=Typen%AquaTemp
    Arr(33)=Typen%AquaTemp1
    Arr(34)=Typen%AquaTemp2
    Arr(35)=Typen%AquaTemp3
    Arr(36)=Typen%Special
    Arr(37)=Typen%DTemp
    Arr(38)=Typen%DTemp1
    Arr(39)=Typen%DTemp2
    Arr(40)=Typen%DTemp3
    Arr(41)=Typen%DTemp4
    Arr(42)=Typen%DTemp5
    Arr(43)=Typen%Meskhidze
    Arr(44)=Typen%Equi
    Arr(45)=Typen%SolidSpecial
    Arr(46)=Typen%Parti
    Arr(47)=Typen%Microphys
  END SUBROUTINE mkArray
  !
  !
  SUBROUTINE InitNReacType()
    NTypes%GasPhoto=0
    NTypes%GasPhotAB=0
    NTypes%GasPhotABC=0
    NTypes%GasPhotMCM=0
    NTypes%GasConst=0
    NTypes%Temp=0
    NTypes%Temp1=0
    NTypes%Temp2=0
    NTypes%Temp3=0
    NTypes%Troe=0
    NTypes%Troef=0
    NTypes%TroeQ=0
    NTypes%Spec1=0
    NTypes%Spec2=0
    NTypes%Spec3=0
    NTypes%Spec4=0
    NTypes%Spec1MCM=0
    NTypes%Spec2MCM=0
    NTypes%Spec3MCM=0
    NTypes%SPec4MCM=0
    NTypes%Spec5MCM=0
    NTypes%Spec6MCM=0
    NTypes%Spec7MCM=0
    NTypes%Spec8MCM=0
    NTypes%S4H2O=0
    NTypes%Henry=0
    NTypes%AquaPhoto=0
    NTypes%AquaPhotAB=0
    NTypes%AquaPhotABC=0
    NTypes%AquaPhotMCM =0
    NTypes%AquaConst=0
    NTypes%AquaTemp=0
    NTypes%AquaTemp1=0
    NTypes%AquaTemp2=0
    NTypes%AquaTemp3=0
    NTypes%Special=0
    NTypes%DTemp=0
    NTypes%DTemp1=0
    NTypes%DTemp2=0
    NTypes%DTemp3=0
    NTypes%DTemp4=0
    NTypes%DTemp5=0
    NTypes%Meskhidze=0
    NTypes%Equi=0
    NTypes%SolidSpecial=0
    NTypes%Parti=0
    NTypes%Microphys=0
    NTypes%SolidDTemp3=0
  END SUBROUTINE InitNReacType
  !
  !
  SUBROUTINE ReadReaction(Out)
    LOGICAL :: Out
    !
    INTEGER :: iLine,PosColon,Pos,is
    CHARACTER(LenLine) :: LocLine
    CHARACTER(LenLine) :: Line(1:4)
    CHARACTER(20) :: Type
    CHARACTER(40) :: TypeR
    !
    !
    !
    iLine=0
    DO
      READ(InputUnit,'(a400)',IOSTAT=is) LocLine
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
          SELECT CASE (TypeR)
            !CASE ('PHOTO','PHOTMCM')
            CASE ('PHOTO','PHOTAB','PHOTABC','PHOTMCM')
              IF (TypeR=='PHOTO')    NTypes%GasPhoto=NTypes%GasPhoto+1
              IF (TypeR=='PHOTAB')   NTypes%GasPhotAB=NTypes%GasPhotAB+1
              IF (TypeR=='PHOTABC')  NTypes%GasPhotABC=NTypes%GasPhotABC+1
              IF (TypeR=='PHOTMCM')  NTypes%GasPhotMCM=NTypes%GasPhotMCM+1
              nreakgphoto=nreakgphoto+1
            CASE ('CONST')
              NTypes%GasConst=NTypes%GasConst+1
              nreakgconst=nreakgconst+1
            CASE ('TEMP','TEMP1','TEMP2','TEMP3')
              IF (TypeR=='TEMP')  NTypes%Temp=NTypes%Temp+1
              IF (TypeR=='TEMP1') NTypes%Temp1=NTypes%Temp1+1
              IF (TypeR=='TEMP2') NTypes%Temp2=NTypes%Temp2+1
              IF (TypeR=='TEMP3') NTypes%Temp3=NTypes%Temp3+1
              nreakgtemp=nreakgtemp+1
            CASE ('TROE','TROEF','TROEQ')
              IF (TypeR=='TROE')  NTypes%Troe=NTypes%Troe+1
              IF (TypeR=='TROEF') NTypes%TroeF=NTypes%TroeF+1
              IF (TypeR=='TROEQ') NTypes%TRoeQ=NTypes%TroeQ+1
              nreakgtroe=nreakgtroe+1
            CASE ('SPEC1','SPEC2','SPEC3','SPEC4','SPEC1MCM','SPEC2MCM',            &
            &     'SPEC3MCM','SPEC4MCM','SPEC5MCM','SPEC6MCM','SPEC7MCM','SPEC8MCM')
              IF (TypeR=='SPEC1') NTypes%Spec1=NTypes%Spec1+1
              IF (TypeR=='SPEC2') NTypes%Spec2=NTypes%Spec2+1
              IF (TypeR=='SPEC3') NTypes%Spec3=NTypes%Spec3+1
              IF (TypeR=='SPEC4') NTypes%Spec4=NTypes%Spec4+1
              IF (TypeR=='SPEC1MCM') NTypes%Spec1MCM=NTypes%Spec1MCM+1
              IF (TypeR=='SPEC2MCM') NTypes%Spec2MCM=NTypes%Spec2MCM+1
              IF (TypeR=='SPEC3MCM') NTypes%Spec3MCM=NTypes%Spec3MCM+1
              IF (TypeR=='SPEC4MCM') NTypes%Spec4MCM=NTypes%Spec4MCM+1
              IF (TypeR=='SPEC5MCM') NTypes%Spec5MCM=NTypes%Spec5MCM+1
              IF (TypeR=='SPEC6MCM') NTypes%Spec6MCM=NTypes%Spec6MCM+1
              IF (TypeR=='SPEC7MCM') NTypes%Spec7MCM=NTypes%Spec7MCM+1
              IF (TypeR=='SPEC8MCM') NTypes%Spec8MCM=NTypes%Spec8MCM+1
              nreakgspec=nreakgspec+1
            CASE ('S4H2O')
              NTypes%S4H2O=NTypes%S4H2O+1
              nreakgspec=nreakgspec+1
          END SELECT
        CASE ('HENRY')
          nreakhenry=nreakhenry+1
          NTypes%Henry=NTypes%Henry+1
          CALL InsertReaction(ListRHenry,Line,TypeR)
        CASE ('AQUA')
          nreakaqua=nreakaqua+1
          CALL InsertReaction(ListRAqua,Line,TypeR)
          SELECT CASE (TypeR)
            !CASE ('PHOTO','PHOTMCM')
            CASE ('PHOTO','PHOTAB','PHOTABC','PHOTMCM')
              IF (TypeR=='PHOTO')   NTypes%AquaPhoto=NTypes%AquaPhoto+1
              IF (TypeR=='PHOTAB')  NTypes%AquaPhotAB=NTypes%AquaPhotAB+1
              IF (TypeR=='PHOTABC') NTypes%AquaPhotABC=NTypes%AquaPhotABC+1
              IF (TypeR=='PHOTMCM') NTypes%AquaPhotMCM=NTypes%AquaPhotMCM+1
              nreakaphoto=nreakaphoto+1
            CASE ('CONST')
              NTypes%AquaConst=NTypes%AquaConst+1
              nreakaconst=nreakaconst+1
            CASE ('TEMP','Temp1''TEMP2','TEMP3')
              IF (TypeR=='TEMP')   NTypes%AquaTemp=NTypes%AquaTemp+1
              IF (TypeR=='TEMP1')   NTypes%AquaTemp1=NTypes%AquaTemp1+1
              IF (TypeR=='TEMP2')   NTypes%AquaTemp2=NTypes%AquaTemp2+1
              IF (TypeR=='TEMP3')   NTypes%AquaTemp3=NTypes%AquaTemp3+1
              nreakatemp=nreakatemp+1
            CASE ('SPECIAL')
              NTypes%Special=NTypes%Special+1
              nreakaspec=nreakaspec+1
          END SELECT
        CASE ('DISS')
          nreakdissoc=nreakdissoc+1
          CALL InsertReaction(ListRDiss,Line,TypeR)
          SELECT CASE (TypeR)
            CASE ('DTEMP','DTEMP1','DTEMP2','DTEMP3','DTEMP4','DTEMP5','MESKHIDZE')
              IF (TypeR=='DTEMP')   NTypes%DTemp=NTypes%DTemp+1
              IF (TypeR=='DTEMP1')   NTypes%DTemp1=NTypes%DTemp1+1
              IF (TypeR=='DTEMP2')   NTypes%DTemp2=NTypes%DTemp2+1
              IF (TypeR=='DTEMP3')   NTypes%DTemp3=NTypes%DTemp3+1
              IF (TypeR=='DTEMP4')   NTypes%DTemp4=NTypes%DTemp4+1
              IF (TypeR=='DTEMP5')   NTypes%DTemp5=NTypes%DTemp5+1
              IF (TypeR=='MESKHIDZE')   NTypes%Meskhidze=NTypes%Meskhidze+1
          END SELECT
        CASE ('SOLID')
          nreaksolid=nreaksolid+1
          CALL InsertReaction(ListRSolid,Line,TypeR)
          SELECT CASE (TypeR)
            CASE ('EQUI')
              NTypes%Equi=NTypes%Equi+1
              nreaksolidEqui=nreaksolidEqui+1
            CASE ('DTEMP3')
              NTypes%SolidDTemp3=NTypes%SolidDTemp3+1
              nreaksolidtemp=nreaksolidtemp+1
            CASE ('SPECIAL')
              NTypes%SolidSpecial=NTypes%SolidSpecial+1
              nreaksolidspec=nreaksolidspec+1
          END SELECT
        CASE ('PARTI')
          NTypes%Parti=NTypes%Parti+1
          NumberReactionsPartic=NumberReactionsPartic+1
          CALL InsertReaction(ListRPartic,Line,TypeR)
        CASE ('MICROPHYS')
          NTypes%Microphys=NTypes%Microphys+1
          NumberReactionsMicro=NumberReactionsMicro+1
          CALL InsertReaction(ListRMicro,Line,TypeR)
      END SELECT
      Out=.FALSE.
    ELSE
      Out=.TRUE.
    END IF
  END SUBROUTINE ReadReaction
  !========
  !
  SUBROUTINE mkParty(ReacStruct)
    TYPE(ReactionStruct_T) :: ReacStruct(:)
    !
    INTEGER, ALLOCATABLE :: DuctsPHOTOreac(:)     &
    &                     , DuctsCONSTreac(:)     &
    &                     , DuctsTEMPreac(:)      &
    &                     , DuctsTROEreac(:)      &
    &                     , DuctsSPECIALreac(:)   &
    &                     , DuctsHENRYreac(:)     &
    &                     , DuctsAPHOTOreac(:)    &
    &                     , DuctsACONSTreac(:)    & 
    &                     , DuctsATEMPreac(:)     & 
    &                     , DuctsASPECIALreac(:)  & 
    &                     , DuctsDISSreac(:)      & 
    &                     , DuctsSEQUIreac(:)     & 
    &                     , DuctsSSPECIALreac(:)  & 
    &                     , DuctsSTEMPreac(:)     &
    &                     , DuctsPARTIreac(:)     &
    &                     , DuctsMICPHYSreac(:)     
    !
    INTEGER :: cntPHOTO=0, cntCONST=0,  cntTEMP=0      &
    &        , cntTROE=0 , cntSPEC=0                   &
    &        , cntHENRY=0, cntAPHOTO=0, cntACONST=0    &
    &        , cntATEMP=0, cntASPEC=0,  cntDISS=0      &
    &        , cntSEQUI=0, cntSTEMP=0,  cntSSPEC=0     &
    &        , cntPARTI=0, cntMICPHYS=0
    !
    INTEGER, ALLOCATABLE :: PermVec(:)
    INTEGER :: i, j, ColLen
    !
    ALLOCATE(DuctsPHOTOreac(2*nspc))
    ALLOCATE(DuctsCONSTreac(3*nspc))
    ALLOCATE(DuctsTEMPreac(5*nspc))
    ALLOCATE(DuctsTROEreac(2*nspc))
    ALLOCATE(DuctsSPECIALreac(2*nspc))
    ALLOCATE(DuctsHENRYreac(2*nspc))
    ALLOCATE(DuctsAPHOTOreac(2*nspc))
    ALLOCATE(DuctsACONSTreac(2*nspc))
    ALLOCATE(DuctsATEMPreac(2*nspc))
    ALLOCATE(DuctsASPECIALreac(2*nspc))
    ALLOCATE(DuctsDISSreac(2*nspc))
    ALLOCATE(DuctsSEQUIreac(2*nspc))
    ALLOCATE(DuctsSSPECIALreac(2*nspc))
    ALLOCATE(DuctsSTEMPreac(2*nspc))
    ALLOCATE(DuctsPARTIreac(2*nspc))
    ALLOCATE(DuctsMICPHYSreac(2*nspc))
    ALLOCATE(PermVec(2*nspc))
    DuctsPHOTOreac=0
    DuctsCONSTreac=0
    DuctsTEMPreac=0
    DuctsTROEreac=0
    DuctsSPECIALreac=0
    DuctsHENRYreac=0
    DuctsAPHOTOreac=0
    DuctsACONSTreac=0
    DuctsATEMPreac=0
    DuctsASPECIALreac=0
    DuctsDISSreac=0
    DuctsSEQUIreac=0
    DuctsSSPECIALreac=0
    DuctsSTEMPreac=0
    DuctsPARTIreac=0
    DuctsMICPHYSreac=0
    PermVec=0
    !
    DO i=1,neq
      SELECT CASE (ReacStruct(i)%Type)
        CASE ('GAS')
          SELECT CASE (ReacStruct(i)%TypeConstant)
            CASE ('PHOTO','PHOTAB','PHOTABC','PHOTMCM')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntPHOTO=cntPHOTO+1
                DuctsPHOTOreac(cntPHOTO)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
            CASE ('CONST')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntCONST=cntCONST+1
                DuctsCONSTreac(cntCONST)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
            CASE ('TEMP','TEMP1','TEMP2')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntTEMP=cntTEMP+1
                DuctsTEMPreac(cntTEMP)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
            CASE ('TROE','TROEF','TROEQ')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntTROE=cntTROE+1
                DuctsTROEreac(cntTROE)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
            CASE ('SPEC1','SPEC2','SPEC3','SPEC4','SPEC1MCM','SPEC2MCM',            &
            &     'SPEC3MCM','SPEC4MCM','SPEC5MCM','SPEC6MCM','SPEC7MCM','SPEC8MCM')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntSPEC=cntSPEC+1
                DuctsSPECIALreac(cntSPEC)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
            CASE ('SPECIAL')
          END SELECT
        CASE ('HENRY')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntHENRY=cntHENRY+1
                DuctsHENRYreac(cntHENRY)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
        CASE ('AQUA')
          SELECT CASE (ReacStruct(i)%TypeConstant)
            CASE ('PHOTO','PHOTAB','PHOTABC','PHOTMCM')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntAPHOTO=cntAPHOTO+1
                DuctsAPHOTOreac(cntAPHOTO)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
            CASE ('CONST')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntACONST=cntACONST+1
                DuctsACONSTreac(cntACONST)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
            CASE ('TEMP','TEMP1','TEMP3')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntATEMP=cntATEMP+1
                DuctsATEMPreac(cntATEMP)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
            CASE ('SPECIAL')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntASPEC=cntASPEC+1
                DuctsASPECIALreac(cntASPEC)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
          END SELECT
        CASE ('DISS')
          SELECT CASE (ReacStruct(i)%TypeConstant)
            CASE ('DTEMP','DTEMP2','DTEMP3','DTEMP4','DTEMP5','MESKHIDZE')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntDISS=cntDISS+1
                DuctsDISSreac(cntDISS)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
          END SELECT
        CASE ('SOLID')
          SELECT CASE (ReacStruct(i)%TypeConstant)
            CASE ('EQUI')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntSEQUI=cntSEQUI+1
                DuctsSEQUIreac(cntSEQUI)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
            CASE ('DTEMP3')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntSTEMP=cntSTEMP+1
                DuctsSTEMPreac(cntSTEMP)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
            CASE ('SPECIAL')
              DO j=1,SIZE(ReacStruct(i)%Educt)
                cntSSPEC=cntSSPEC+1
                DuctsSSPECIALreac(cntSSPEC)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
                &                                                         %Species)
              END DO
          END SELECT
        CASE ('PARTI')
          DO j=1,SIZE(ReacStruct(i)%Educt)
            cntPARTI=cntPARTI+1
            DuctsPARTIreac(cntPARTI)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
            &                                                         %Species)
          END DO
        CASE ('MICROPHYS')
          DO j=1,SIZE(ReacStruct(i)%Educt)
            cntMICPHYS=cntMICPHYS+1
            DuctsMICPHYSreac(cntMICPHYS)=PositionSpeciesAll(ReacStruct(i)%Educt(j)  &
            &                                                         %Species)
          END DO
      END SELECT
    END DO
    CALL CompressParty(DuctsPHOTOreac,PermVec,ColLen)
    CALL CompressParty(DuctsCONSTreac,PermVec,ColLen)
    CALL CompressParty(DuctsTEMPreac,PermVec,ColLen)
    CALL CompressParty(DuctsTROEreac,PermVec,ColLen)
    CALL CompressParty(DuctsSPECIALreac,PermVec,ColLen)
    CALL CompressParty(DuctsHENRYreac,PermVec,ColLen)
    CALL CompressParty(DuctsACONSTreac,PermVec,ColLen)
    CALL CompressParty(DuctsAPHOTOreac,PermVec,ColLen)
    CALL CompressParty(DuctsATEMPreac,PermVec,ColLen)
    CALL CompressParty(DuctsASPECIALreac,PermVec,ColLen)
    CALL CompressParty(DuctsDISSreac,PermVec,ColLen)
    CALL CompressParty(DuctsSEQUIreac,PermVec,ColLen)
    CALL CompressParty(DuctsSTEMPreac,PermVec,ColLen)
    CALL CompressParty(DuctsSSPECIALreac,PermVec,ColLen)
  END SUBROUTINE mkParty
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
    WRITE(Unit,*) ' ========                 MODMEP / TESTVERSION      ========'
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
    WRITE(Unit,*) nreak,            '        NREAK   : Number of Reactions'
    WRITE(Unit,*) nreakgas,         '        NGAS   : Gas phase reactions'
    WRITE(Unit,*) nreakgphoto,    '           Gaseous PHOTO - type reactions'
    WRITE(Unit,*) nreakgconst,    '           Gaseous CONST - type reactions'
    WRITE(Unit,*) nreakgtemp,     '           Gaseous TEMP - type reactions'
    WRITE(Unit,*) nreakgtroe,     '           Gaseous TROE - type reactions'
    WRITE(Unit,*) nreakgspec,  '           Gaseous SPECIAL - type reactions'
    WRITE(Unit,*) nreakhenry,       '        NHENRY : Henry Equilib. reactions'
    WRITE(Unit,*) nreakdissoc,        '        NDISS  : Dissociation reactions'
    WRITE(Unit,*) nreakaqua,        '        NAQUA  : Aquatic Equilib. reactions'
    WRITE(Unit,*) nreakaphoto,   '           Aqueous PHOTO - type reactions'
    WRITE(Unit,*) nreakaconst,   '           Aqueous CONST - type reactions'
    WRITE(Unit,*) nreakatemp,    '           Aqueous TEMP - type reactions'
    WRITE(Unit,*) nreakaspec, '           Aqueous SPECIAL - type reactions'
    WRITE(Unit,*) NumberReactionsPartic,      '        NPARTI  : Particulare reactions   '
    WRITE(Unit,*) nreaksolid,       '        NSOLID  : Solid Equilib. reactions'
    WRITE(Unit,*) nreaksolidtemp, '           Solid DTEMP3 - type reactions'
    WRITE(Unit,*) nreaksolidEqui,   '           Solid EQUI - type reactions'
    WRITE(Unit,*) nreaksolidspec,'           Solid SPECIAL - type reactions'
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
    REAL(RealKind), ALLOCATABLE :: tmpValA(:),tmpValB(:)
    INTEGER, ALLOCATABLE :: APermVec(:)
    INTEGER :: AColLen
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
      DO i=1,SIZE(ReactionSystem(iLoop)%Educt)
        SELECT CASE(ReactionSystem(iLoop)%Educt(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','GAS')
            NumActiveEduct=NumActiveEduct+1
            ActiveEduct(NumActiveEduct)=ReactionSystem(iLoop)%Educt(i)
        END SELECT
      END DO
      ! count activ products in reaction iLoop
      NumActiveProduct=0
      DO i=1,SIZE(ReactionSystem(iLoop)%Product)
        SELECT CASE(ReactionSystem(iLoop)%Product(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','GAS')
            NumActiveProduct=NumActiveProduct+1
            ActiveProduct(NumActiveProduct)=ReactionSystem(iLoop)%Product(i)
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
      IF (ReactionSystem(iLoop)%Factor=='$RO2') hasRO2=.TRUE.
      IF (ReactionSystem(iLoop)%Factor=='$RO2aq') hasRO2aq=.TRUE.
      !
      IF (PRESENT(CK)) WRITE(Unit,*) 'FACTOR:  ',ADJUSTL(TRIM(ReactionSystem(iLoop)%Line2))
      IF (PRESENT(CK)) WRITE(Unit,*) 'FACTOR:  ',ADJUSTL(TRIM(ReactionSystem(iLoop)%Line3))
    END DO
    !
    ! loop again to set ColInd and Val on A and B
    NumberReaction=0
    ValCntA=0
    ValCntB=1
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
    sumBAT=0.0d0
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
          !
          ! this is for the TempX  backward reaction 
          sumBAT(iLoop)=sumBAT(iLoop)-tmpValA(m)
        END DO
        DEALLOCATE(APermVec)
      ELSE
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
          B%ColInd(ValCntB)=tmpColB(m)
          B%Val(ValCntB)=tmpValB(m)
          !
          ! this is for the TempX  backward reaction 
          sumBAT(iLoop)=sumBAT(iLoop)+tmpValB(m)
          !
          ValCntB=ValCntB+1
        END DO
        DEALLOCATE(APermVec)
      ELSE
        DO i=1,NumActiveProduct
          B%ColInd(ValCntB)=PositionSpeciesAll(ActiveProduct(i)%Species)
          B%Val(ValCntB)=ActiveProduct(i)%Koeff
          sumBAT(iLoop)=sumBAT(iLoop)+ActiveProduct(i)%Koeff
          ValCntB=ValCntB+1
        END DO
      END IF
      IF (ALLOCATED(tmpColB)) DEALLOCATE(tmpColB)
      IF (ALLOCATED(tmpValB)) DEALLOCATE(tmpValB)
      IF (ALLOCATED(tmpColA)) DEALLOCATE(tmpColA)
      IF (ALLOCATED(tmpValA)) DEALLOCATE(tmpValA)
    END DO
  END SUBROUTINE PrintReactions
  !
  !
  ! Read Chemical Data (initial values and Emisions)
  SUBROUTINE InputChemicalData(InitFileName,DataFileName,MeteoFileName)
    CHARACTER(*) :: InitFileName, DataFileName, MeteoFileName
     !
    !REAL(RealKind), ALLOCATABLE :: GasInAct(:), AqInAct(:)
    !
    INTEGER, PARAMETER :: GasRateUnit=0 ! ???
    INTEGER :: jt, i, iPos
    REAL(RealKind) :: pi43, LWC
    CHARACTER(60) :: string = ''
    !
    ! for pH Start
    REAL(RealKind) :: kappa
    REAL(RealKind), EXTERNAL :: pHValue
    !
    pi43=4.0d0/3.0d0*PI
    !
    ! this is for mass transfer (accom , diffus term)
    ALLOCATE( y_c1(nspc+ntkat), y_c2(nspc+ntkat) )
    y_c1(1:ntGas)=5.0d-6
    y_c2(1:ntGas)=5.0d-5
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
        y_c1(nspc+i)=5.d-6
        y_c2(nspc+i)=5.d-5
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
    FORALL (jt=1:ntGas) GasName(jt)=ADJUSTL(y_name(jt))
    FORALL (jt=1:ntAqua) AquaName(jt)=ADJUSTL(y_name(ntGas+jt))
    FORALL (jt=1:ntKat) PassName(jt)=ADJUSTL(y_name(ntGas+ntAqua+jt))
    !
    !=========================================================================
    !===  Set  Chemical DATA  
    !===  (Molar Mass, Charges, Densities) 
    !=========================================================================
    !
    IF ( ntAqua>=1 ) THEN
      !---  Allocate arrays
      ALLOCATE (Charge(ntAqua))             ! charge of ions
      ALLOCATE (SolubInd(ntAqua))           ! solubility
      ALLOCATE (MolMass(ntAqua))            ! molar mass of species
      ALLOCATE (SpcDens(ntAqua))            ! density of species
      ALLOCATE (OrgIndex(ntAqua))           ! carbon atoms
      ALLOCATE (CC(ntAqua))                 ! compound class
      ALLOCATE (ActIndex(ntAqua))           ! index for calculation of activity coefficient
  
      !---  Set default values
      Charge(:)=0.d0
      SolubInd(:)=0.d0
      MolMass(:)=0.d0
      SpcDens(:)=1.d0
      OrgIndex(:)=0.d0
      CC(:)='  '
      ActIndex(:)=0.d0
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
    !
    !--- Read thermodynamic data,....
    CALL Read_SpeciesData(y_c1,y_c2,DataFileName)
    !---------------------------------
    !
    !---  Read/Calculate aqua phase initials
    !---------------------------------
    IF ( ntAqua>=1 ) THEN
      !
      CALL Read_AQUAini( InitValAct(ntGas+1:), InitValKat(:), InitFileName )
      CALL Read_AFRAC( InitAFrac, InitFileName )
      CALL Read_SPEK( SPEK, InitFileName)
      !
      !
      LWC=pseudoLWC(tAnf)
      aH2O=aH2OmolperL*LWC*mol2part 
      !
      IF ( ntKat>=1 ) THEN
        DO i=1,ntKat
          IF (y_name(nspc+i)=='[aH2O]') THEN
            InitValKat(i)=InitValKat(i)*LWC*mol2part    ! convert fluessigwasser to molec/cm3
            aH2O = InitValKat(i)
          END IF
        END DO
      END IF
      !
      !-----------------------------------------
      ! calculate initial aqueus concentrations 
      !-----------------------------------------
      InitValAct(ntGAS+1:)=1.0d-16
      DO i=1,SIZE(InitAFrac)
        iPos=PositionSpeciesAll(InitAFrac(i)%Species)
        IF (iPos>0) THEN
          !
          InitValAct(iPos) = SPEK(1)%Number*1.0d+03          &  ! [#/m3]
          &                * InitAFrac(i)%Frac1              &  ! [g/g]
          &                * (pi43*(SPEK(1)%Radius)**3)      &  ! [m3]
          &                * SPEK(1)%Density                 &  ! [kg/m3]
          &                / (InitAFrac(i)%MolMass)             ! 1/[kg/mol] 
          InitValAct(iPos) = InitValAct(iPos) * mol2Part        ! *[molec/mol]
        END IF
      END DO
    END IF
    !
    !----------------------------------
    ! initial values of passive species 
    !----------------------------------
    DO i=1,SIZE(InitValKat)
      IF (ListNonReac2(i)%Species/='[aH2O]') THEN
        InitValKat(i)=InitValKat(i)*GasFac
      END IF
    END DO
    !======================================================
    !---  Compute pH value and number of ions
    !---  Initial pH by charge balance 
    IF ( pHSet>=1.AND.ntAqua>0 )  THEN
      Kappa = pHValue(InitValAct(ntGas+1:))
      IF ( Kappa>0.d0 )  THEN
        InitValAct(Hp_ind)=Kappa
      ELSE 
        InitValAct(OHm_ind)=InitValAct(OHm_ind)+InitValAct(Hp_ind)-Kappa
      END IF
    END IF
    !
  END SUBROUTINE InputChemicalData
  !
  !
  !
  SUBROUTINE Read_SpeciesData(y_acc,y_diff,FileName)
    REAL(RealKind) :: y_acc(:) , y_diff(:) 
    CHARACTER(*) :: FileName
    !
    !
    CHARACTER(240) :: line
    CHARACTER(100) :: SpeciesName
    INTEGER :: iPos, i
    LOGICAL :: Back=.FALSE.
    REAL(RealKind) :: mm, alpha, dg, c1
    REAL(RealKind) :: nue
    CHARACTER(10) :: ro2d
    CHARACTER(10) :: c2
    INTEGER :: slash
    INTEGER, ALLOCATABLE :: allRO2(:)
    CHARACTER(100), ALLOCATABLE :: allRO2name(:)
   
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
    IF ( ntGAS>0 ) THEN
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
      !
      ALLOCATE(allRO2(c1))
      ALLOCATE(allRO2name(c1))
      allRO2(:)=0
      allRO2name(:)='dummy'
      !
      !
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
        IF ( Back ) EXIT
        slash=INDEX(SpeciesName,'_')
        IF ( slash>0 ) THEN
          SpeciesName(slash:slash)='/'
        END IF
        IF (PositionSpeciesAll(SpeciesName)>0) THEN
          i=i+1
          allRO2name(i)=SpeciesName
          allRO2(i)=PositionSpeciesAll(SpeciesName)
          !print*, 'debug:: ro2 :',i,PositionSpeciesAll(SpeciesName), SpeciesName
        END IF
      END DO
      CALL RewindFile
      CALL ClearIniFile
    END IF
    !
    c1=c1-COUNT(allRO2==0,1)
    ALLOCATE(RO2spcG(c1))
    ALLOCATE(RO2(c1))
    RO2spcG(:)=allRO2name(1:c1)
    RO2(:)=allRO2(1:c1)
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
      DO
        CALL LineFile( Back, Start1='BEGIN_DATARO2aq',  &
        &              End='END_DATARO2aq',             &
        &              Name1=SpeciesName)
        IF (Back) EXIT
        IF (PositionSpeciesAll(SpeciesName)>0) i=i+1
      END DO
      IF (i>0) THEN
        ALLOCATE(RO2spcA(i))
        RO2spcA=''
        ALLOCATE(RO2aq(i))
        RO2aq=0
        CALL RewindFile
        CALL ClearIniFile
        c1=0
        CALL LineFile( Back, Start1='BEGIN_DATARO2aq',  &
        &              End='END_DATARO2aq',             &
        &              Name1=ro2d,                    &
        &              R1=c1)
        !
        i=0
        DO
          CALL LineFile( Back, Start1='BEGIN_DATARO2aq',  &
          &              End='END_DATARO2aq',             &
          &              Name1=SpeciesName)
          IF (Back) EXIT
          IF (PositionSpeciesALL(SpeciesName)>0) THEN
            i=i+1
            RO2spcA(i)=SpeciesName
            RO2aq(i)=PositionSpeciesAll(SpeciesName)
          END IF
        END DO
      END IF
    END IF
    CALL CloseIniFile
      
      WRITE(333,*) ' nRO2=',SIZE(RO2)
    DO i=1,SIZE(RO2)
      WRITE(333,*) i, RO2(i)
    END DO 
      WRITE(333,*) ' nRO2aq=',SIZE(RO2aq)
    DO i=1,SIZE(RO2aq)
      WRITE(333,*) i, RO2aq(i)
    END DO 
  END SUBROUTINE Read_SpeciesData
  !
  SUBROUTINE Read_EMISS(FileName,Emi)
    CHARACTER(*) :: FileName
    REAL(RealKind) :: Emi(:)
    !
    INTEGER :: iPos
    REAL(RealKind) :: c1
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
    REAL(RealKind) :: GASact(:)
    REAL(RealKind) :: InAct(:)
    !
    INTEGER :: iPos
    REAL(RealKind) :: c1
    CHARACTER(20) :: SpeciesName
    LOGICAL :: Back=.FALSE.
    !
    !
    GASact=1.0d-20
    InAct=1.0d-20
    ! read initial values
    CALL OpenIniFile(FileName)
    DO
      CALL LineFile( Back,                     &
      &              Start1='BEGIN_GAS',       &
      &              Start2='BEGIN_INITIAL',   &
      &              End   ='END_INITIAL',     &
      &              Name1 =SpeciesName, R1=c1 )
      IF (Back)   EXIT
      !
      iPos=PositionSpeciesAll(SpeciesName)
      IF (iPos>0) THEN
        IF (SpeciesName(1:1)=='[') THEN
          InAct(iPos-nspc)=c1
        ELSE
          iPos=PositionSpeciesAll(SpeciesName)
          GASact(iPos)=c1
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
    REAL(RealKind) :: AQUAact(:)
    REAL(RealKind) :: AquaInAct(:)
    CHARACTER(*) :: FileName
    !
    INTEGER :: iPos
    REAL(RealKind) :: c1
    CHARACTER(20) :: SpeciesName
    LOGICAL :: Back
    !
    ! Read initial values of aqua spc
    CALL OpenIniFile(FileName)
    DO
      CALL LineFile( Back,                    &
      &              Start1='BEGIN_AQUA',     &
      &              Start2='BEGIN_INITIAL',  &
      &              End   ='END_INITIAL',    &
      &              Name1=SpeciesName,R1=c1  )
      IF (Back)   EXIT
      !
      iPos=PositionSpecies(SpeciesName)
      IF (iPos>0) THEN
        IF (SpeciesName(1:1)=='[') THEN
          AquaInAct(iPos)=c1
        ELSE
          AquaAct(iPos)=c1
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
    REAL(RealKind) :: c1,c2,c3,c4
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
    REAL(RealKind) :: c1,c2,c3
    LOGICAL :: Back
    REAL(RealKind) :: LWC
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
        SPEK(1)%Radius=REAL(c1,KIND=RealKind)
        SPEK(1)%Number=REAL(c2,KIND=RealKind)
        SPEK(1)%wetRadius=(3.0d0/4.0d0/PI*LWC/SPEK(1)%Number)**(1.0d0/3.0d0)
        SPEK(1)%Density=REAL(c3,KIND=RealKind)
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
    REAL(RealKind) :: c1
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
    REAL(RealKind), POINTER :: Constants(:)
    CHARACTER(*) :: Type
    !
    INTEGER :: NumColon,PosColon,PosName,PosComment
    INTEGER :: i,PosNum1,PosNum2,PosNum3,NumNum,PosElem
    CHARACTER(4)  :: NameNumNum
    CHARACTER(10) :: DummyString
    CHARACTER(LEN(String)) :: LocString
    CHARACTER(LEN(String)) :: LocString1
    CHARACTER(LEN(String)) :: NameConstant
    REAL(RealKind) :: Dummy
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
    REAL(RealKind) :: PreFac !NumberSpecies
    CHARACTER(LenLine) :: Species
    CHARACTER(LEN(String)) :: LocString
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
        IF (Species(1:1)=='['.AND.LEN(TRIM(Species))<maxLENinActDuct) THEN
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
    ! PARTIC
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
        PositionSpeciesAll=PositionSpeciesAll+ntGAS+ntAQUA+ntSOLID+ntPart
      END IF
    ! GAS
    ELSE
      PositionSpeciesAll=GetHash(ListGas,TRIM(ADJUSTL(Species)))
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
    INTEGER :: i, j, iList, TmpArraySize
    INTEGER :: nList
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
      idxHp=PositionSpeciesAll('Hp')
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
    !
    ALLOCATE(ReacStruct(neq))
    ALLOCATE(ReactionTypes(neq))
    !
    i=1
    !
    DO iList=1,nList
      Current=>CompleteReactionList(iList)%Start
      DO WHILE (ASSOCIATED(Current)) 
        ReacStruct(i)%Type=Current%Type
        ReactionTypes(i)=Current%Type
        ReacStruct(i)%TypeConstant=Current%TypeConstant
        ReacStruct(i)%Line1=Current%Line1
        ReacStruct(i)%Line2=Current%Line2
        ReacStruct(i)%Line3=Current%Line3
        ReacStruct(i)%Factor=Current%Factor
        !
        TmpArraySize=SIZE(Current%Educt)
        ReacStruct(i)%nActEd=tmpArraySize
        ALLOCATE(ReacStruct(i)%Educt(TmpArraySize))
        DO j=1,TmpArraySize
          ReacStruct(i)%Educt(j)%Species=Current%Educt(j)%Species
          ReacStruct(i)%Educt(j)%Type=Current%Educt(j)%Type
          ReacStruct(i)%Educt(j)%Koeff=Current%Educt(j)%Koeff
        END DO
        !
        TmpArraySize=SIZE(Current%Product)
        ReacStruct(i)%nActPro=TmpArraySize
        ALLOCATE(ReacStruct(i)%Product(TmpArraySize))
        DO j=1,TmpArraySize
          ReacStruct(i)%Product(j)%Species=Current%Product(j)%Species
          ReacStruct(i)%Product(j)%Type=Current%Product(j)%Type
          ReacStruct(i)%Product(j)%Koeff=Current%Product(j)%Koeff
        END DO
        !
        ReacStruct(i)%NumConst=SIZE(Current%Constants)
        ALLOCATE(ReacStruct(i)%Constants(ReacStruct(i)%NumConst))
        ReacStruct(i)%Constants=Current%Constants
        !
        IF (Current%Type=='HENRY') THEN
          ReacStruct(i)%direction='GA'
          ReacStruct(i)%HenrySpc=PositionSpeciesAll(ReacStruct(i)%Educt(1)%Species)
          !ReacStruct(i)%Constants(1)=ReacStruct(i)%Constants(1)
        END IF
        !
        !
        !
        TmpArraySize=SIZE(Current%Educt)
        ALLOCATE(ReacStruct(i)%InActEduct(TmpArraySize))
        ALLOCATE(ReacStruct(i)%InActEductSpc(TmpArraySize))
        DO j=1,TmpArraySize
          ReacStruct(i)%InActEduct(j)=Current%InActEduct(j)%Koeff
          ReacStruct(i)%InActEductSpc(j)=Current%InActEduct(j)%Species
        END DO
        !
        TmpArraySize=SIZE(Current%InActProduct)
        ALLOCATE(ReacStruct(i)%InActProduct(TmpArraySize))
        DO j=1,TmpArraySize
          ReacStruct(i)%InActProduct(j)=Current%InActProduct(j)%Koeff    
        END DO
        !    
        ReacStruct(i)%nInActEd=Current%nInActEd
        ReacStruct(i)%nInActPro=Current%nInActPro
        !
        ReacStruct(i)%SumAqCoef=SUM(Current%Educt(:)%Koeff)-1.0d0
        !
        !IF (TRIM(ADJUSTL(ReacStruct(i)%InActEductSpc(1)))=='[aH2O]') THEN 
        !  ReacStruct(i)%SumAqCoef=ReacStruct(i)%SumAqCoef+1
        !END IF
        ! for equilibrium reactions save <-- direction
        SELECT CASE (Current%Type)
          CASE ('DISS','HENRY')
            i=i+1
            ReacStruct(i)%Type=Current%Type
            ReactionTypes(i)=Current%Type
            ReacStruct(i)%TypeConstant=Current%TypeConstant
            ReacStruct(i)%Line1=Current%Line1
            ReacStruct(i)%Line2='BackReaction'
            ReacStruct(i)%Line3=Current%Line3
            ReacStruct(i)%Factor=Current%Factor
            !
            TmpArraySize=SIZE(Current%Product)
            ReacStruct(i)%nActEd=TmpArraySize
            ALLOCATE(ReacStruct(i)%Educt(TmpArraySize))
            DO j=1,TmpArraySize
              ReacStruct(i)%Educt(j)%Species=Current%Product(j)%Species
              ReacStruct(i)%Educt(j)%Type=Current%Product(j)%Type
              ReacStruct(i)%Educt(j)%Koeff=Current%Product(j)%Koeff
            END DO
            !
            TmpArraySize=SIZE(Current%Educt)
            ReacStruct(i)%nActPro=TmpArraySize
            ALLOCATE(ReacStruct(i)%Product(TmpArraySize))
            DO j=1,TmpArraySize
              ReacStruct(i)%Product(j)%Species=Current%Educt(j)%Species
              ReacStruct(i)%Product(j)%Type=Current%Educt(j)%Type
              ReacStruct(i)%Product(j)%Koeff=Current%Educt(j)%Koeff
            END DO
            !
            ReacStruct(i)%NumConst=SIZE(Current%Constants)
            ALLOCATE(ReacStruct(i)%Constants(ReacStruct(i)%NumConst))
            ReacStruct(i)%Constants=Current%Constants
            !
            IF (Current%Type=='HENRY') THEN
              ReacStruct(i)%direction='AG'
              ReacStruct(i)%HenrySpc=PositionSpeciesAll(ReacStruct(i)%Product(1)%Species)
            END IF
            !
            !
            TmpArraySize=SIZE(Current%InActProduct)
            ALLOCATE(ReacStruct(i)%InActEduct(TmpArraySize))
            ALLOCATE(ReacStruct(i)%InActEductSpc(TmpArraySize))
            DO j=1,TmpArraySize
              ReacStruct(i)%InActEduct(j)=Current%InActProduct(j)%Koeff
              ReacStruct(i)%InActEductSpc(j)=Current%InActProduct(j)%Species
            END DO
            !
            TmpArraySize=SIZE(Current%InActEduct)
            ALLOCATE(ReacStruct(i)%InActProduct(TmpArraySize))
            DO j=1,TmpArraySize
              ReacStruct(i)%InActProduct(j)=Current%InActEduct(j)%Koeff    
            END DO
            !    
            ReacStruct(i)%nInActEd=Current%nInActPro
            ReacStruct(i)%nInActPro=Current%nInActEd
            !
            !
            ReacStruct(i)%SumAqCoef=SUM(Current%Product(:)%Koeff)-1.0d0
            !
            IF (TRIM(ADJUSTL(ReacStruct(i)%InActEductSpc(1)))=='[aH2O]') THEN 
              ReacStruct(i)%SumAqCoef=ReacStruct(i)%SumAqCoef+1
            END IF
            !
        END SELECT
        !
        Current=>Current%Next
        i=i+1
      END DO
    END DO
  END SUBROUTINE AllListsToArray
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
END MODULE Chemsys_Mod
