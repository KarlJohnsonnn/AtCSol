MODULE mo_ckinput
  !
  USE hashtbl
  USE Kind_Mod
  USE Meteo_Mod
  USE ChemSys_Mod, ONLY: ReactionStruct_t, Duct_T, ListAqua, ListGas        &
  &                    , ListSolid, ListPartic, ListNonReac, ReactionSystem &
  &                    , UnitGas, UnitAqua, ListGas2, Species_T             &
  &                    , ListToHashTable, HashTableToList, SortList         &
  &                    , PositionSpeciesGas
  !
  USE mo_reac, ONLY: ntGas, ntAqua, ntSolid, ntPart, ntKat, neq, nspc, nReak&
  &                    , nreakgas,  nreakgconst, nreakgphoto, nreakgspec    &
  &                    , nreakgtemp, nreakgtroe, nDIM                       &
  &                    , lowA,lowB,lowC,lowD,lowE,lowF,lowG                 &
  &                    , highA,highB,highC,highD,highE,highF,highG
  USE NetCDF_Mod
  !
  IMPLICIT NONE
  !
  !
  INTEGER, PRIVATE, PARAMETER :: LenLine=405
  INTEGER, PRIVATE, PARAMETER :: LenName=100
  !
  CHARACTER(51), PRIVATE, PARAMETER :: CKSetSpecies='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz()'
  CHARACTER(14), PRIVATE, PARAMETER :: SetConstants='ABCDEFGKINMOR/'
  CHARACTER(11), PRIVATE, PARAMETER :: CKSetNumber='0123456789.'
  !
  LOGICAL :: Comment
  !
CONTAINS
  !
  !
  SUBROUTINE OpenFile(nUnit,FileName,Type)
    INTEGER :: nUnit
    CHARACTER(*) :: Filename
    CHARACTER(*) :: Type
    !
    LOGICAL :: ExistFile
    !
    INQUIRE(FILE=TRIM(Filename)//'.'//TRIM(Type),EXIST=ExistFile)
    IF (ExistFile) THEN
      OPEN(UNIT=nUnit,FILE=TRIM(Filename)//'.'//TRIM(Type),STATUS='UNKNOWN')
    END IF
  END SUBROUTINE OpenFile
  !
  !
  SUBROUTINE CloseFile(nUnit,FileName,Type)
    INTEGER :: nUnit
    CHARACTER(*) :: Filename
    CHARACTER(*) :: Type
    !
    LOGICAL :: ExistFile
    !
    INQUIRE(FILE=TRIM(Filename)//'.'//TRIM(Type),EXIST=ExistFile)
    IF (ExistFile) CLOSE(UNIT=nUnit)
  END SUBROUTINE CloseFile
  !
  !
  SUBROUTINE Read_Thermodata(ThermSwitchTemp,DatThermo,UnitThermo,nSpc)
    CHARACTER(*) :: DatThermo
    INTEGER :: UnitThermo, nSpc
    
    CHARACTER(80) :: iLine   !i =1,2,3,4
    CHARACTER(18) :: Species(nSpc)
    CHARACTER(18) :: SpeciesInfo(nspc)
    CHARACTER(6)  :: RefDataCode(nSpc)
    CHARACTER(2)  :: Atoms(nSpc,4)
    INTEGER       :: nAtoms(nSpc,4)
    REAL(8)       :: TempRange(nSpc,2)
    REAL(8)       :: MolMass(nSpc)
    CHARACTER(1)  :: Phase(nSpc)
    REAL(8)       :: ta,tb,tc,td,te,tf,tg
    REAL(8)       :: H0_29815R(nspc)
    REAL(8),ALLOCATABLE :: ThermSwitchTemp(:)
    !
    INTEGER :: i,j,n
    INTEGER :: idxWhiteSpace
    INTEGER :: ALLOC_ERR
    CHARACTER(18) :: tSpcName
    LOGICAL :: THERMO=.FALSE.
    !
    ! FORMATS
    1 FORMAT(A18,A6,4(A2,I3),A1,2F10.3,F13.7,A1,I1) 
    2 FORMAT(5E15.8,I1)
    !
    !
    CALL OpenFile(UnitThermo,DatThermo(1:INDEX(DatThermo,'.')-1),'dat')
    ALLOCATE(lowA(nspc),lowB(nspc),lowC(nspc),lowD(nspc),lowE(nspc),lowF(nspc),lowG(nspc),STAT=ALLOC_ERR)
    ALLOCATE(highA(nspc),highB(nspc),highC(nspc),highD(nspc),highE(nspc),highF(nspc),highG(nspc),STAT=ALLOC_ERR)
    ALLOCATE(ThermSwitchTemp(nspc))
    !
    i=0
    DO
      !
      READ(UnitThermo,'(A80)') iLine
      IF (iLine(1:3)=='END'.OR.iLine(1:3)=='end') EXIT
      !
      IF (THERMO.OR.TRIM(ADJUSTL(iLine))=='thermo'.OR.TRIM(ADJUSTL(iLine))=='THERMO') THEN
        IF (.NOT.THERMO) THERMO=.TRUE.
        SELECT CASE (iLine(80:80))
          CASE ('1')
            !i=i+1
            READ(iLine,1) tSpcName
            ! hier statt inkrem. i speziesnummer suchen
            tSpcName=TRIM(ADJUSTL(tSpcName))
            i=PositionSpeciesGas(tSpcName(1:INDEX(tSpcName,' ')-1))
            !
            READ(iLine,1) tSpcName,RefDataCode(i),(Atoms(i,j),    &
            &             nAtoms(i,j),j=1,4),Phase(i),ta,tb,tc,n
            idxWhiteSpace=INDEX(tSpcName,' ')
            !
            WRITE(Species(i),*) tSpcName(1:idxWhiteSpace-1)
            WRITE(SpeciesInfo(i),*) tSpcName(idxWhiteSpace:)
            !
            TempRange(i,1)=REAL(ta,KIND=RealKind)
            TempRange(i,2)=REAL(tb,KIND=RealKind)
            ThermSwitchTemp(i)=TempRange(i,2)
            MolMass(i)=REAL(tc,KIND=RealKind)
          CASE ('2')
            READ(iLine,2) ta,tb,tc,td,te,n
            highA(i)=REAL(ta,KIND=RealKind)
            highB(i)=REAL(tb,KIND=RealKind)
            highC(i)=REAL(tc,KIND=RealKind)
            highD(i)=REAL(td,KIND=RealKind)
            highE(i)=REAL(te,KIND=RealKind)
          CASE ('3')
            READ(iLine,2) tf,tg,ta,tb,tc,n
            highF(i)=REAL(tf,KIND=RealKind)
            highG(i)=REAL(tg,KIND=RealKind)
            lowA(i)=REAL(ta,KIND=RealKind)
            lowB(i)=REAL(tb,KIND=RealKind)
            lowC(i)=REAL(tc,KIND=RealKind)
          CASE ('4')
            READ(iLine,2) td,te,tf,tg,ta,n
            lowD(i)=REAL(td,KIND=RealKind)
            lowE(i)=REAL(te,KIND=RealKind)
            lowF(i)=REAL(tf,KIND=RealKind)
            lowG(i)=REAL(tg,KIND=RealKind)
            H0_29815R(i)=REAL(ta,KIND=RealKind)
          CASE DEFAULT
            CONTINUE
        END SELECT
      END IF
    END DO
    CALL CloseFile(UnitThermo,DatThermo(1:INDEX(DatThermo,'.')-1),'dat')
  END SUBROUTINE Read_Thermodata
  !
  !
  SUBROUTINE Read_Elements(DataReac,UnitReac)
    CHARACTER(*) :: DataReac
    INTEGER :: UnitReac
    
    INTEGER, PARAMETER :: nMaxElements=130
    CHARACTER(200) :: iLine   !i =1,2,3,4
    CHARACTER(2)   :: tElem(nMaxElements)
    CHARACTER(2),ALLOCATABLE   :: Elements(:)
    !
    INTEGER :: i
    INTEGER :: iWS
    !
    CALL OpenFile(UnitReac,DataReac(1:LEN(DataReac)-4),'sys')
    !
    i=0
    DO
      READ(UnitReac,'(A200)') iLine
      !
      ! check if next line is end
      IF (iLine(1:3)=='END'.OR.iLine(1:3)=='end') THEN
        CALL CutCArray(Elements,tElem,i)
        EXIT
      END IF
      !
      IF (iLine(1:4)=='ELEM'.OR.iLine(1:4)=='elem') THEN
        i=0
        READ(UnitReac,'(A200)') iLine
        !
        DO
          i=i+1
          iWS=INDEX(ADJUSTL(iLine),' ')
          WRITE(tElem(i),'(A2)') ADJUSTL(iLine(:iWS-1))
          iLine=ADJUSTL(iLine(iWS:))
          !
          IF (iLine=='') EXIT
        END DO
      END IF
    END DO
    CALL CloseFile(UnitReac,DataReac(1:LEN(DataReac)-4),'sys')
  END SUBROUTINE Read_Elements
  !
  !
  SUBROUTINE Read_Species(DataReac,UnitReac)
    CHARACTER(*) :: DataReac
    INTEGER :: UnitReac
    !
    INTEGER, PARAMETER :: nMaxElements=130
    CHARACTER(200)     :: iLine   !i =1,2,3,4
    !
    CHARACTER(200)     :: headline  
    INTEGER :: i
    INTEGER :: iWS
    !
    CHARACTER(10) :: phase
    !
    !
    ! Init Hash Tables
    ntGas=0
    ntAqua=0
    ntSolid=0
    ntPart=0
    ntkat=0
    CALL Init_HashTbls(ListGas)
    !
    CALL OpenFile(UnitReac,TRIM(DataReac),'sys')
    !
    i=0
    CALL FindSection(UnitReac,'species',headline)
    DO
    i=i+1
      READ(UnitReac,'(A200)') iLine
      IF (iLine(1:3)=='END'.OR.iLine(1:3)=='end') EXIT
      !
      DO
        IF (ADJUSTL(iLine(1:3))=='end'.OR.iLine=='') EXIT
        iLine=ADJUSTL(iLine)
        iWS=INDEX(iLine,' ')
        CALL InsertSpecies(iLine(1:iWS-1),phase)
        iLine=iLine(iWS:)
      END DO
    END DO
    CALL CloseFile(UnitReac,TRIM(DataReac),'sys')
    !
    nspc=ntGas+ntAqua+ntSolid+ntPart+ntkat

    !
    !print*, 'Total number of species = ', nSpc
  END SUBROUTINE Read_Species
  !
  SUBROUTINE Read_Reaction(DataReac,UnitReac)
    CHARACTER(*) :: DataReac
    INTEGER :: UnitReac
    !
    CHARACTER(200)     :: iLine   !i =1,2,3,4
    CHARACTER(200)     :: locString
    CHARACTER(200)     :: headline  
    CHARACTER(200)     :: ductstr
    CHARACTER(200)     :: dummyString
    !
    INTEGER :: i
    INTEGER :: iReac
    INTEGER :: iWS
    INTEGER :: iNxtSpc
    INTEGER :: iM, iKlammerM
    !
    INTEGER :: iKl,iKr
    CHARACTER(10) :: auxiliary
    !
    INTEGER :: nduct, nEducts, nProducts
    INTEGER :: fPosPlus , fPosEq, fPosFw
    INTEGER :: idxDucts(6)!,idxDuctsP(6)
    !REAL(RealKind), ALLOCATABLE :: valDucts(:)
    INTEGER, ALLOCATABLE :: SPCind(:,:)
    !
    !INTEGEr :: TableNspc
    CHARACTER(15) :: units
    REAL(RealKind) :: tmpReal
    !
    LOGICAL :: bR
    !
    CALL OpenFile(UnitReac,TRIM(DataReac),'sys')
    REWIND UnitReac
    !
    CALL FindSection(UnitReac,'reactions',headline)
    !
    ! unit
    iWS=INDEX(headline,' ')
    headline=TRIM(ADJUSTL(headline(iWs:)))
    SELECT CASE (headline)
      CASE ('CAL/MOLE','cal/mole')
        units='CAL/MOLE'
      CASE ('KCAL/MOLE','kcal/mole')
        units='KCAL/MOLE'
      CASE ('JOULES/MOLE','joules/mole')
        units='JOULES/MOLE'
      CASE DEFAULT
        units='CAL/MOLE'
        WRITE(*,*) ' Using cal/mole '
    END SELECT
    !
    ! first count reations, ...
    iReac=0
    DO
      ! read next line
      READ(UnitReac,'(A200)') iLine
      ! adjust lef
      iLine=ADJUSTL(iLine)
      ! exit cond
      IF ( MAX(INDEX(iLine,'END'),INDEX(iLine,'end'))>0 )  EXIT
      !
      !print*, 'debug::mockin   iline=',iline
      ! if no exit iReac++
      !          reversible reactions
      IF ( MAX(INDEX(iLine,' = '),INDEX(iLine,' <=> '))>0 )  iReac=iReac+2
      !          irreversible reactions
      IF ( INDEX(iLine,' => ')>0 )  iReac=iReac+1
    END DO
    REWIND UnitReac
    nReak=iReac
    !print*, 'debug::mockin   nreak=',nreak
    !CALL PrintHashTableCK(ListGas,TableNspc)
    !IF (TableNspc>nSpc) STOP ' More Species in reactin system than declared in species section! '
    !
    ! ALLOCATE reaction system structur
    ALLOCATE(ReactionSystem(nReak))
    ! ALLOCATE matrices alpha, beta, (beta-alpha)^T
    !alpha%m=nReak
    !alpha%n=nSpc
    !beta%m =nReak
    !beta%n =nSpc
    !ALLOCATE(alpha%RowPtr(nReak+1),beta%RowPtr(nReak+1))
    !alpha%RowPtr(1)=1
    !beta%RowPtr(1) =1
    !
    ! allocate temp array for spc indices (max 6 educts, 6 products)
    ALLOCATE(SPCind(nReak,12))
    SPCind=-1
    !
    !
    CALL FindSection(UnitReac,'reactions',headline)
    iReac     = 0
    fPosPlus  = 0
    fPosEq    = 0
    fPosFw    = 0
    iLine       = ''
    dummyString = ''
    LocString   = ''
    DO
      ! exit cond
      IF ( MAX(INDEX(iLine,'END'),INDEX(iLine,'end'))>0 ) EXIT
      !IF (iLine(1:3)=='END'.OR.iLine(1:3)=='end') EXIT
      !
      IF ( TRIM(ADJUSTL(LocString))==dummyString) THEN
        ! read next line
        READ(UnitReac,'(A200)') iLine
      ELSE
        iLine=TRIM(ADJUSTL(LocString))
      END IF
      !
      dummyString=TRIM(ADJUSTL(iLine))
      !print*, 'DEBUG::mo_ckinput   vor skip  iLine=',TRIM(iLine)
      !
      CALL SkipLines(UnitReac,iLine)
      !print*, 'DEBUG::mo_ckinput   nachskip  iLine=',TRIM(iLine)
      !
      ! adjust left and cut off comments
      IF (INDEX(iLine,'!')>0)  THEN
        iLine=iLine(1:INDEX(iLine,'!')-1)
      END IF
      iLine=ADJUSTL(iLine)
      !
      !
      ! if no exit iReac++
      !          reversible reactions
      fPosEq=MAX(INDEX(iLine,' = '),INDEX(iLine,' <=> '))
      fPosFw=INDEX(iLine,' => ')
      !
      ! IF REACTION LINE
      IF (fPosEq>0.OR.fPosFW>0)  THEN
        iReac=iReac+1
        ReactionSystem(iReac)%Type='GAS'    ! immer gas bei chemkin?
        !
        ! if  =   or   <=>   reaction
        IF (fPosEq>0) THEN
          bR=.TRUE.
          ReactionSystem(iReac+1)%Type='GAS'    ! immer gas bei chemkin?
        END IF
        !
        ! FIRST LINE OF REACTION  ---  Extract Arrhenius Coeff 
        ALLOCATE(ReactionSystem(iReac)%Constants(3))
        IF (bR) ALLOCATE(ReactionSystem(iReac+1)%Constants(3))
        LocString=iLine
        DO i=3,1,-1
          LocString=ADJUSTR(LocString)
          iWS=INDEX(LocString,' ',.TRUE.)        ! find first whitespace <--
          READ(LocString(iWS:),'(E10.4)') tmpReal
          ReactionSystem(iReac)%Constants(i)=REAL(tmpReal,KIND=RealKind)
          IF (bR) ReactionSystem(iReac+1)%Constants(i)=REAL(tmpReal,KIND=RealKind)
          LocString=LocString(:iWS)
        END DO
        !
        ! left over string = reaction
        LocString=TRIM(ADJUSTL(LocString))
        !
        IF (fPosEq>0) THEN
          ReactionSystem(iReac)%Line1=TRIM(ADJUSTL(LocString(1:fPosEq-1)))//' => '// &
          &                           TRIM(ADJUSTL(LocString(fPosEq+4:)))   ! save reaction string
          !print*, 'DEBUG::mo_ckinput  hinreaktion = ',iReac,TRIM(ReactionSystem(iReac)%Line1)
        ELSE
          ReactionSystem(iReac)%Line1=LocString               ! save reaction string
        END IF
        !
        nduct=0
        !
        !
        !p
        IF (bR) THEN
          ! save reaction string
          ReactionSystem(iReac+1)%Line1=TRIM(ADJUSTL(LocString(fPosEq+4:)))//' => '// &
          &                             TRIM(ADJUSTL(LocString(1:fPosEq-1)))
          !print*, 'DEBUG::mo_ckinput  rückrealinks= ',TRIM(ADJUSTL(LocString(fPosEq+4:)))
          !print*, 'DEBUG::mo_ckinput  rückrearecht= ',TRIM(ADJUSTL(LocString(1:fPosEq-1)))
          !print*, 'DEBUG::mo_ckinput  rückreaktion= ',iReac+1,TRIM(ReactionSystem(iReac+1)%Line1)
        END IF
        ! extract the constant type by checking the reaction if +m or (+M) appears, and cut M off
        !print*, 'DEBUG::mo_ckinput getkonst  vor locstr= ',locstring
        CALL GetConstantType(iM,iKlammerM,LocString,ReactionSystem(iReac)%TypeConstant,ReactionSystem(iReac)%Factor)
        !print*, 'DEBUG::mo_ckinput getkonst nach locstr= ',locstring
        !
        ! count educts and products
        CALL CountDooku(nEducts,nProducts,LocString)
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
        !
        !print*, ' debug::mockin   iR, nedu,npro, factor=', iReac,neducts,nproducts,ReactionSystem(iReac)%Factor

        !print*, 'DEBUG::mo_ck..    ned,npr=',neducts,nproducts,locstring
        ALLOCATE(ReactionSystem(iReac)%Educt(nEducts))
        IF (bR) ALLOCATE(ReactionSystem(iReac+1)%Educt(nProducts))
        ALLOCATE(ReactionSystem(iReac)%Product(nProducts))
        IF (bR) ALLOCATE(ReactionSystem(iReac+1)%Product(nEducts))
        !
        !
        IF (bR) THEN
          ! save reaction string
          ReactionSystem(iReac+1)%TypeConstant=ReactionSystem(iReac)%TypeConstant
          ReactionSystem(iReac+1)%Factor=ReactionSystem(iReac)%Factor
        END IF
        !
        fPosEq=INDEX(LocString,'=')
        ductStr=ADJUSTL(LocString(1:fPosEq-2))
        !
!----------------------------------------------------------------------------------------------
! 
!       build row pointer array for stoechiometric coefficients of matrix alpha
!       (left side of reaction
!
!----------------------------------------------------------------------------------------------
        !ALLOCATE(idxDucts(6))
        !ALLOCATE(idxDuctsP(6))
        !
        DO i=1,nEducts
          fPosPlus=INDEX(ductStr,'+')
          !
          CALL ExtractKoeff(iNxtSpc,tmpReal,ductStr)
          !
          IF (fPosPlus>0) THEN
            ReactionSystem(iReac)%Educt(i)%Species=TRIM(ADJUSTL(ductStr(iNxtSpc:fPosPlus-1)))
            print*, 'debug::mock..       spc  =  ', TRIM(ADJUSTL(ductStr(iNxtSpc:fPosPlus-1)))
            IF (bR) ReactionSystem(iReac+1)%Product(i)%Species=TRIM(ADJUSTL(ductStr(iNxtSpc:fPosPlus-1)))
          ELSE
            ReactionSystem(iReac)%Educt(i)%Species=TRIM(ADJUSTL(ductStr(iNxtSpc:)))
            print*, 'debug::mock.. last  spc  =  ', TRIM(ADJUSTL(ductStr(iNxtSpc:)))
            IF (bR) ReactionSystem(iReac+1)%Product(i)%Species=TRIM(ADJUSTL(ductStr(iNxtSpc:)))
          END IF
          !
          ! forward
          SPCind(iReac,i)=PositionSpeciesGas(ReactionSystem(iReac)%Educt(i)%Species)
          ReactionSystem(iReac)%Educt(i)%Type='GAS'                ! immer gas in chemkin?
          ReactionSystem(iReac)%Educt(i)%Koeff=REAL(tmpReal,KIND=RealKind)
          !
          idxDucts(i)=PositionSpeciesGas(ReactionSystem(iReac)%Educt(i)%Species)
          !
          ! backward
          IF (bR) THEN
            SPCind(iReac+1,6+i)=PositionSpeciesGas(ReactionSystem(iReac)%Educt(i)%Species)
            ReactionSystem(iReac+1)%Product(i)%Type='GAS'                ! immer gas in chemkin?
            ReactionSystem(iReac+1)%Product(i)%Koeff=REAL(tmpReal,KIND=RealKind)
          END IF
          !
          ductStr=TRIM(ADJUSTL(ductStr(fPosPlus+1:)))  
        END DO
        !
        !stop
        ! HIER NOCH COLIND FÜR DIE STÖCH MATRIZEN BESTIMMEN EINFACH MIT SPCIND ARBEITEN
        !
        !IF (nEducts>1) THEN
        !  alpha%RowPtr(iReac+1)=alpha%RowPtr(iReac)+nEducts
        !ELSE
        !  alpha%RowPtr(iReac+1)=alpha%RowPtr(iReac)+1
        !END IF
        !
        ! PRODUCT SIDE OF REACTION
        !
        !
        ductStr=ADJUSTL(LocString(fPosEq+2:))
        DO i=1,nProducts
          fPosPlus=INDEX(ductStr,'+')
          !
          CALL ExtractKoeff(iNxtSpc,tmpReal,ductStr)
          !
          IF (fPosPlus>0) THEN
            ReactionSystem(iReac)%Product(i)%Species=TRIM(ADJUSTL(ductStr(iNxtSpc:fPosPlus-1)))
            IF (bR) ReactionSystem(iReac+1)%Educt(i)%Species=TRIM(ADJUSTL(ductStr(iNxtSpc:fPosPlus-1)))
          ELSE
            ReactionSystem(iReac)%Product(i)%Species=TRIM(ADJUSTL(ductStr(iNxtSpc:)))
            IF (bR) ReactionSystem(iReac+1)%Educt(i)%Species=TRIM(ADJUSTL(ductStr(iNxtSpc:)))
          END IF
          ! forward
          SPCind(iReac,6+i)=PositionSpeciesGas(ReactionSystem(iReac)%Product(i)%Species)
          ReactionSystem(iReac)%Product(i)%Type='GAS'                
          ReactionSystem(iReac)%Product(i)%Koeff=REAL(tmpReal,KIND=RealKind)
          ! backward
          IF (bR) THEN
            SPCind(iReac+1,i)=PositionSpeciesGas(ReactionSystem(iReac)%Product(i)%Species)
            ReactionSystem(iReac+1)%Educt(i)%Type='GAS'                
            ReactionSystem(iReac+1)%Educt(i)%Koeff=REAL(tmpReal,KIND=RealKind)
          END IF
          !
          ductStr=TRIM(ADJUSTL(ductStr(fPosPlus+1:)))  
        END DO
        !
        !
        !IF (nProducts>1) THEN
        !  beta%RowPtr(iReac+1)=beta%RowPtr(iReac)+nProducts
        !  IF (bR) alpha%RowPtr(iReac+2)=alpha%RowPtr(iReac+1)+nProducts
        !  IF (bR) beta%RowPtr(iReac+2)=Beta%RowPtr(iReac+1)+nEducts
        !ELSE
        !  beta%RowPtr(iReac+1)=beta%RowPtr(iReac)+1
        !  IF (bR) beta%RowPtr(iReac+2)=beta%RowPtr(iReac+1)+1
        !  IF (bR) alpha%RowPtr(iReac+2)=alpha%RowPtr(iReac+1)+1
        !END IF
        !
        !
        IF (bR) THEN
          !
          iLine=''
          ! skip empty lines and comment lines

          !print*, 'DEBUG::mo_ck.. vor    iR,locstr= ',iReac,locstring
          !CALL SkipLines(UnitReac,iLine)
          LocString=ADJUSTL(iLine)
          !print*, 'DEBUG::mo_ck.. nach   iR,locstr= ',iReac,locstring
          !
          iKl=INDEX(LocString,'/')
          iKr=INDEX(LocString,'/',.TRUE.)
          !
          auxiliary=ADJUSTL(LocString(1:iKl-1))
          !
          ! if there are extra coefs for backreaction
          IF (auxiliary=='rev'.OR.auxiliary=='REV') THEN
            !ReactionSystem(iReac+1)%Factor=auxiliary
            ReactionSystem(iReac+1)%Line3=auxiliary
            IF (.NOT.ALLOCATED(ReactionSystem(iReac+1)%Constants)) ALLOCATE(ReactionSystem(iReac+1)%Constants(3))
            LocString=ADJUSTL(LocString(ikl+1:ikr-1))
            DO i=1,3
              iWS=INDEX(LocString,' ')
              READ(LocString(1:iWS-1),'(E10.4)') tmpReal
              ReactionSystem(iReac+1)%Constants(i)=REAL(tmpReal,KIND=RealKind)
              LocString=ADJUSTL(LocString(iWs+1:))
            END DO
          END IF
        END IF
        !
        ! print backreaction
        IF (bR) THEN
          ReactionSystem(iReac+1)%Line2='BackReaction'
        END IF
      END IF ! Reaction => or ( = or <=> )
      !
      !==================================================
      ! IF AN    +m, (+m)     OR   +M, (+M) is involved
      !==================================================
      ! handle additional parameters for reverse reaction etc.
      !
      IF (ReactionSystem(iReac)%Factor=='$(+M)') THEN
        !
        IF (auxiliary/='rev'.OR.auxiliary/='REV') THEN
          ! hier sind wir schon in der richtigen zeile und es muss nichts geskippt werden  
        ELSE
          iLine=''
        END IF
        ! skip empty lines and comment lines
        CALL SkipLines(UnitReac,iLine)
        LocString=ADJUSTL(iLine)
        !
        iKl=INDEX(LocString,'/')
        iKr=INDEX(LocString,'/',.TRUE.)
        IF (iKl>0) THEN
          auxiliary=ADJUSTL(LocString(1:iKl-1))
          LocString=ADJUSTL(LocString(iKl+1:iKr-1))
          ReactionSystem(iReac)%Line3=TRIM(auxiliary)
          IF (bR) ReactionSystem(iReac+1)%Line3=TRIM(auxiliary)
        ELSE
          ReactionSystem(iReac)%Line3=''
          IF (bR) ReactionSystem(iReac+1)%Line3=''
        END IF
        !
        !
        IF (ReactionSystem(iReac)%Line3=='low'.OR.ReactionSystem(iReac)%Line3=='LOW') THEN
          ALLOCATE(ReactionSystem(iReac)%LowConst(3))
          ALLOCATE(ReactionSystem(iReac+1)%LowConst(3))
          DO i=1,3
            iWS=INDEX(LocString,' ')
            READ(LocString(1:iWS-1),'(E10.4)') tmpReal
            ReactionSystem(iReac)%LowConst(i)=REAL(tmpReal,KIND=RealKind)
            ReactionSystem(iReac+1)%LowConst(i)=REAL(tmpReal,KIND=RealKind)
            LocString=ADJUSTL(LocString(iWs+1:))
          END DO
        ELSEIF (ReactionSystem(iReac)%Line3=='high'.OR.ReactionSystem(iReac)%Line3=='HIGH') THEN
          ALLOCATE(ReactionSystem(iReac)%HighConst(3))
          DO i=1,3
            iWS=INDEX(LocString,' ')
            READ(LocString(1:iWS-1),'(E10.4)') tmpReal
            ReactionSystem(iReac)%HighConst(i)=REAL(tmpReal,KIND=RealKind)
            ReactionSystem(iReac+1)%HighConst(i)=REAL(tmpReal,KIND=RealKind)
            LocString=ADJUSTL(LocString(iWs+1:))
          END DO
        END IF
        !
        iLine=''
        ! skip empty lines and comment lines
        CALL SkipLines(UnitReac,iLine)
        LocString=ADJUSTL(iLine)
        iKl=INDEX(LocString,'/')
        iKr=INDEX(LocString,'/',.TRUE.)
        !
        auxiliary=ADJUSTL(LocString(1:iKl-1))
        !
        ! extract troe parameters 
        IF (MAXVAL(INDEX(LocString,(/'troe','TROE'/)))>0) THEN
          ALLOCATE(ReactionSystem(iReac)%TroeConst(4))
          ALLOCATE(ReactionSystem(iReac+1)%TroeConst(4))
          LocString=ADJUSTL(LocString(iKl+1:iKr-1))
          DO i=1,4
            iWS=INDEX(LocString,' ')
            READ(LocString(1:iWS-1),'(E10.4)') tmpReal
            ReactionSystem(iReac)%TroeConst(i)=REAL(tmpReal,KIND=RealKind)
            ReactionSystem(iReac+1)%TroeConst(i)=REAL(tmpReal,KIND=RealKind)
            LocString=ADJUSTL(LocString(iWs+1:))
          END DO
        END IF
      END IF
      !   CASE IF +m OR +M is involved in reaction iReac
      IF (ReactionSystem(iReac)%Factor=='$+M'.OR.ReactionSystem(iReac)%Factor=='$(+M)') THEN
        iLine=''          ! clear actual line 
        !
        !print*, 'debug::mockinp    vor line= ',locString
        ! skip empty lines and comment lines and read next line
        CALL SkipLines(UnitReac,iLine)
        !
        LocString=ADJUSTL(iLine)
        !print*, 'debug::mockinp   nach line= ',locString
        !
        iKl=INDEX(LocString,'/')
        iKr=INDEX(LocString,'/',.TRUE.)
        !
        ! extract 3rd body species and alpha values
        IF (iKL>0) THEN
          CALL ComputeThirdBody(  ReactionSystem(iReac)%TB      &
          &                     , ReactionSystem(iReac)%TBspc   &
          &                     , ReactionSystem(iReac)%TBalpha &
          &                     , LocString)
          IF (bR) THEN
            ALLOCATE(ReactionSystem(iReac+1)%TB(SIZE(ReactionSystem(iReac)%TB)))
            ALLOCATE(ReactionSystem(iReac+1)%TBspc(SIZE(ReactionSystem(iReac)%TB)))
            ALLOCATE(ReactionSystem(iReac+1)%TBalpha(SIZE(ReactionSystem(iReac)%TB)))
            ReactionSystem(iReac+1)%TB=ReactionSystem(iReac)%TB
            ReactionSystem(iReac+1)%TBspc=ReactionSystem(iReac)%TBspc
            ReactionSystem(iReac+1)%TBalpha=ReactionSystem(iReac)%TBalpha
          END IF
        END IF
        !
      END IF
      IF (bR) iReac=iReac+1
      IF (bR) bR=.FALSE.
    END DO                      ! next reaction
    CALL CloseFile(UnitReac,DataReac,'sys')
    !
    neq=nReak
    nreakgas=nReak
    ALLOCATE(ListGas2(ntGAS))
    CALL HashTableToList(ListGas,ListGas2)
    CALL SortList(ListGas2)
    CALL ListToHashTable(ListGas2,ListGas)
    !
  END SUBROUTINE Read_Reaction
  !
  !
  SUBROUTINE SkipLines(UnitR,Line)
    INTEGER, INTENT(INOUT) :: UnitR
    CHARACTER(*) :: Line
    !
    ! skip empty lines and comment lines
    DO WHILE (Line(1:1)=='!'.OR.Line=='')
      READ(UnitR,'(A200)') Line
      Line=ADJUSTL(Line)
    END DO
    !
  END SUBROUTINE SkipLines
  !
  !
  SUBROUTINE NewLine(UnitReac,LocString,i,j)
    INTEGER :: UnitReac
    CHARACTER(200) :: iLine
    CHARACTER(200) :: LocString
    INTEGER, OPTIONAL :: i,j
    ! read next line
    ! extract troe parameters 
    READ(UnitReac,'(A200)') iLine
    LocString=ADJUSTL(iLine)
    IF (PRESENT(i).AND.PRESENT(j)) THEN
      i=INDEX(LocString,'/')
      j=INDEX(LocString,'/',.TRUE.)
    END IF
  END SUBROUTINE NewLine
  !
  !
  SUBROUTINE ComputeThirdBody(indM,spcM,aM,Line)
    ! out:
    INTEGER, ALLOCATABLE        :: indM(:)
    CHARACTER(*), ALLOCATABLE   :: spcM(:)
    REAL(RealKind), ALLOCATABLE :: aM(:)
    ! in:
    CHARACTER(*), INTENT(IN) :: Line
    ! temp:
    CHARACTER(LEN(Line)) :: locLine
    INTEGER :: ind, kl1, kl2
    REAL(RealKind) :: tmp
    ind=0
    locLine=Line
    DO
      IF (locLine=='') EXIT
      kl1=INDEX(locLine,'/')
      locLine=ADJUSTL(locLine(kl1+1:))
      kl1=INDEX(locLine,'/')
      locLine=ADJUSTL(locLine(kl1+1:))
      ind=ind+1
    END DO
    !
    ALLOCATE(indM(ind))
    ALLOCATE(spcM(ind))
    ALLOCATE(aM(ind))
    locLine=Line
    ind=0
    DO
      IF (locLine=='') EXIT
      ind=ind+1
      kl1=INDEX(locLine,'/')
      kl2=INDEX(locLine(kl1+1:),'/')
      !print*, 'DEBUG::ckINput_____________'
      !print*, 'DEBUG::ckINpt                  ind=    ', ind
      !print*, 'DEBUG::ckINput             locline=    ', locline
      !print*, 'DEBUG::ckINput              species=   ', locLine(1:kl1-1)
      !print*, 'DEBUG::ckINput  posspc(unsortiert) =   ', PositionSpeciesGas(locLine(1:kl1-1))
      !print*, 'DEBUG::ckINput        locline(kl1:)=   ', locline(kl1+1:)
      !print*, 'DEBUG::ckINput_____________'
      !
      indM(ind)=PositionSpeciesGas(locLine(1:kl1-1))
      spcM(ind)=locLine(1:kl1-1)
      READ(locLine(kl1+1:kl1+kl2-1),'(E12.6)') tmp
      aM(ind)=REAL(tmp,KIND=RealKind)
      !
      locLine=ADJUSTL(locLine(kl1+kl2+1:))
    END DO
  END SUBROUTINE ComputeThirdBody
  !
  !
  SUBROUTINE GetConstantType(idxM,idxM2,String,ConstType,Factor)
    INTEGER, INTENT(OUT) :: idxM,idxM2
    CHARACTER(*),INTENT(INOUT) :: String
    CHARACTER(*),INTENT(INOUT) :: ConstType
    CHARACTER(*),INTENT(INOUT) :: Factor
    !
    INTEGER :: i
    !
    ! save constant type and cut off +m or (+m) (no longer needed)
    DO i=1,2 
    !
      idxM=MAXVAL(INDEX(String,(/'+M','+m'/)))
      idxM2=MAXVAL(INDEX(String,(/'(+M)','(+m)'/)))
      !
      IF (idxM>0) THEN
        IF (idxM2>0) THEN 
          ConstType='PRESSX'
          Factor='$(+M)'
          String(idxM2:idxM2+3)='    '
          nreakgtroe=nreakgtroe+1
        ELSE
          ConstType='TEMPX'
          Factor='$+M'
          String(idxM:idxM+2)='  '
          nreakgtemp=nreakgtemp+1
        END IF
      ELSE
        ! normal arrhenius without 3rd body
        ConstType='TEMPX'
        Factor='None'
        nreakgtemp=nreakgtemp+1
      END IF
    END DO
    !
  END SUBROUTINE GetConstantType
  !
  !
  SUBROUTINE ExtractKoeff(NextSpecies,Koeff,String)
    CHARACTER(*), INTENT(IN) :: String
    !
    INTEGER :: NextSpecies
    REAL(RealKind) :: Koeff
    !
    ! check if spc has coef
    IF (SCAN(String(1:1),CKSetNumber)>0) THEN
      NextSpecies=SCAN(String,CKSetSpecies)            ! find the next letter
      READ(String(1:NextSpecies-1),'(E12.6)') Koeff
    ELSE
      NextSpecies=1
      Koeff=1.0d0
    END IF
  END SUBROUTINE ExtractKoeff
  !
  !
  SUBROUTINE CountDooku(nE,nP,inLine)
    CHARACTER(*), INTENT(IN) :: inLine
    INTEGER :: nE, nP
    !
    CHARACTER(LEN(inLine)) :: Line
    INTEGER :: PosPlus,PosEq
    !INTEGER :: tmp
    !
    INTEGER :: edIndex(6),prIndex(6)
    !
    nE=0
    nP=0
    Line=inLine
    DO
      Line=ADJUSTL(Line)
      PosPlus=INDEX(Line,'+')
      PosEq  =INDEX(Line,'=')
      !
      IF (PosEq>0) THEN
        nE=nE+1
        IF (PosPlus>0.AND.PosPlus<PosEq) THEN
          edIndex(nE)=PositionSpeciesGas(Line(:PosPlus-1))
          Line=ADJUSTL(Line(PosPlus+1:))
        ELSE
          edIndex(nE)=PositionSpeciesGas(Line(:PosEq-1))
          Line=ADJUSTL(Line(PosEq+1:))
        END IF
      ELSE
        EXIT
      END IF
    END DO
    !
    DO
      PosPlus=INDEX(Line,'+')
      nP=nP+1
      IF (PosPlus>0) THEN
        prIndex(nP)=PositionSpeciesGas(Line(:PosPlus-1))
        Line=ADJUSTL(Line(PosPlus+1:))
      ELSE
        prIndex(nP)=PositionSpeciesGas(TRIM(Line(:)))
        EXIT
      END IF
    END DO
    prIndex(nP)=PositionSpeciesGas(Line(:))
  END SUBROUTINE CountDooku
  !
  !
  SUBROUTINE FindSection(UnitReac,SecName,dummy)
    CHARACTER(*) :: SecName
    INTEGER :: UnitReac
    !
    CHARACTER(200) :: dummy
    !
    dummy=ADJUSTL(dummy)
    !
    DO WHILE (dummy(1:LEN(SecName))/=SecName)
      READ(UnitReac,'(A200)') dummy
    END DO
  END SUBROUTINE
  !
  SUBROUTINE CutCArray(OutArr,InArr,nCut)
    CHARACTER(2)               :: InArr(:)
    CHARACTER(2), ALLOCATABLE  :: OutArr(:)
    INTEGER :: nCut
    !
    ALLOCATE(OutArr(nCut))
    OutArr=InArr(1:nCut)
    !
  END SUBROUTINE CutCArray
  !
  !
  SUBROUTINE Init_HashTbls(L1,L2,L3,L4,L5)
    TYPE(hash_tbl_sll)           :: L1
    TYPE(hash_tbl_sll), OPTIONAL :: L2, L3, L4, L5
    !
    CALL InitHashTable(L1,100)
    IF (PRESENT(L2)) CALL InitHashTable(L2,100)
    IF (PRESENT(L3)) CALL InitHashTable(L3,100)
    IF (PRESENT(L4)) CALL InitHashTable(L4,100)
    IF (PRESENT(L5)) CALL InitHashTable(L5,100)
    !
  END SUBROUTINE Init_HashTbls
  !
  !
  SUBROUTINE InsertSpecies(Species,Type)
    CHARACTER(*) :: Species
    CHARACTER(*) :: Type
    !
    IF (Species(1:1)=='p') THEN
      CALL InsertHash(ListPartic,TRIM(ADJUSTL(Species)),ntPart)
      Type='Partic'
    !ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
    !  CALL InsertHash(ListAqua,TRIM(ADJUSTL(Species)),ntAqua)
    !  Type='Aqua'
    ELSE IF (Species(1:1)=='s') THEN
      CALL InsertHash(ListSolid,TRIM(ADJUSTL(Species)),ntSolid)
      Type='Solid'
    ELSE IF (Species(1:1)=='['.AND.LEN(TRIM(Species))<10.AND. &
      &      Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
      CALL InsertHash(ListNonReac,TRIM(ADJUSTL(Species)),ntkat)
      Type='Inert'
    ELSE IF (Species(1:1)=='(') THEN
    ELSE
      CALL InsertHash(ListGas,TRIM(ADJUSTL(Species)),ntGas)
      Type='Gas'
    END IF
  END SUBROUTINE InsertSpecies
  !
  !
  SUBROUTINE GetSpeciesNames(FileName,Species)
    CHARACTER(60), ALLOCATABLE :: Species(:)
    CHARACTER(*)               :: FileName
    !
    INTEGER :: i
    !
    !-- Open .chem-file and skip the first 24 lines (head)
    OPEN(UNIT=89,FILE=ADJUSTL(TRIM(FileName))//'.chem',STATUS='UNKNOWN')
    REWIND(89)
    DO i=1,24
      READ(89,*)
    END DO
    !--- allocate space for a list of all species (incl. kat spc)
    ALLOCATE(Species(nspc))
    !
    DO i=1,ntGas
      READ(89,*)  Species(i)
    END DO
    CLOSE(89)
    !
  END SUBROUTINE GetSpeciesNames
  !
  !
  SUBROUTINE PrintReactionSystem(ReacSys)
    TYPE(ReactionStruct_T) :: ReacSys(:)
    !
    INTEGER ::i,j
    !
    WRITE(*,*) ''
    WRITE(*,*) '  Reaction System  '
    WRITE(*,*) ''
    !
    DO i=1,SIZE(ReacSys)
      WRITE(*,*) '+---------------+'
      WRITE(*,'(A12,I4,A2)') ' | REACTION: ',i,' |'
      WRITE(*,*) '+---------------+-------------------------------------------------'
      WRITE(*,*) '| CLASS:   ',TRIM(ReacSys(i)%Type)
      WRITE(*,*) '|     ',TRIM(ReacSys(i)%Line1)
      IF(TRIM(ReacSys(i)%Line2)/='') WRITE(*,*) '|         *',TRIM(ReacSys(i)%Line2)//'*'
      WRITE(*,*) '| FACTOR:  ',ReacSys(i)%Factor
      !
      WRITE(*,*) '| educts:   ','  koeff ','  Species      ',' Hashtbl index '
      WRITE(*,*) '|--------     ----------------------------------------------------'
      WRITE(*,'(A2,12X,F6.3,4X,A15,3X,I4)') (' |',ReacSys(i)%Educt(j)%Koeff, ReacSys(i)%Educt(j)%Species , &
      &                            PositionSpeciesGas(ReacSys(i)%Educt(j)%Species),         &
      &                            j=1,SIZE(ReacSys(i)%Educt))
      WRITE(*,*) '| products: '
      WRITE(*,*) '|----------   -----------------------------------------------------'
      WRITE(*,'(A2,12X,F6.3,4X,A15,3X,I4)') (' |',ReacSys(i)%Product(j)%Koeff, ReacSys(i)%Product(j)%Species, &
      &                            PositionSpeciesGas(ReacSys(i)%Product(j)%Species),          &
      &                            j=1,SIZE(ReacSys(i)%Product))
      WRITE(*,*) '|'
      IF (ReacSys(i)%Factor=='$+M') THEN
        WRITE(*,*) '| Arrhenius Coefs A , b , E :'
        WRITE(*,'(A2,10X,E10.4)') (' |', ReacSys(i)%Constants(j),j=1,SIZE(ReacSys(i)%Constants))
        WRITE(*,*) '|'
        WRITE(*,*) '| 3rd body species with alpha: '
        WRITE(*,'(A2,5X,I4,F6.2)') (' |', ReacSys(i)%TB(j),ReacSys(i)%TBalpha(j),j=1,SIZE(ReacSys(i)%TB))
      ELSE IF (ReacSys(i)%Factor=='$(+M)') THEN
        IF (ReacSys(i)%Line3=='low'.OR.ReacSys(i)%Line3=='LOW') THEN
          WRITE(*,*) '| fall-off reactions:'
          WRITE(*,*) '| (LOW pressure parameter): ','(HIGH pressure parameter): '
          WRITE(*,'(A2,10X,E10.4,10X,E10.4)') (' |', ReacSys(i)%LowConst(j), ReacSys(i)%Constants(j)    &
          &                                     ,j=1,SIZE(ReacSys(i)%LowConst) )
          WRITE(*,*) '|'
        ELSE IF (ReacSys(i)%Line3=='high'.OR.ReacSys(i)%Line3=='HIGH') THEN
          WRITE(*,*) '| fall-off reactions:'
          WRITE(*,*) '| (HIGH pressure parameter): ','(LOW pressure parameter): '
          WRITE(*,'(A2,10X,E10.4,10X,E10.4)') (' |', ReacSys(i)%HighConst(j),ReacSys(i)%Constants(j)   &
          &                                    ,j=1,SIZE(ReacSys(i)%HighConst))
          WRITE(*,*) '|'
        END IF
        WRITE(*,*) '| Troe Parameter for pressur dependen reaction: '
        WRITE(*,'(A2,5X,E10.4)') (' |', ReacSys(i)%TroeConst(j),j=1,SIZE(ReacSys(i)%TroeConst))
        WRITE(*,*) '|'
        WRITE(*,*) '| 3rd body species with alpha: '
        WRITE(*,'(A2,5X,I4,F6.2)') (' |', ReacSys(i)%TB(j),ReacSys(i)%TBalpha(j),j=1,SIZE(ReacSys(i)%TB))
        WRITE(*,*) '|'
      ELSE
        WRITE(*,*) '| Arrhenius Coefs A , b , E :'
        WRITE(*,'(A2,10X,E10.4)') (' |', ReacSys(i)%Constants(j) ,j=1,SIZE(ReacSys(i)%Constants))
        WRITE(*,*) '|'
      END IF
      WRITE(*,*) '|-----------------------------------------------------------------'
      WRITE(*,*) ''
    END DO
  END SUBROUTINE PrintReactionSystem
  !
  !
  SUBROUTINE PrintHeadReactionsCK(Unit)
    INTEGER :: Unit
    !
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
    WRITE(Unit,*)
    WRITE(Unit,*) ' ======================  Reactions   ========================'
    WRITE(Unit,*) ''
  END SUBROUTINE PrintHeadReactionsCK
  !
END MODULE mo_ckinput
