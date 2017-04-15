MODULE mo_ckinput
  !
  USE hashtbl
  USE Kind_Mod
  USE Meteo_Mod
  USE ChemSys_Mod, ONLY: ReactionStruct_t, Duct_T, ListAqua, ListGas        &
  &                    , ListSolid, ListPartic, ListNonReac, ReactionSystem &
  &                    , UnitGas, UnitAqua, ListGas2, Species_T             &
  &                    , ListToHashTable, HashTableToList, SortList         &
  &                    , PositionSpeciesGas,PositionSpeciesAll
  !
  USE mo_reac,     ONLY: ntGas, ntAqua, ntSolid, ntPart, ntKat, neq, nspc, nReak&
  &                    , nreakgas,  nreakgconst, nreakgphoto, nreakgspec    &
  &                    , nreakgtemp, nreakgtroe, nDIM                       &
  &                    , lowA,lowB,lowC,lowD,lowE,lowF,lowG                 &
  &                    , highA,highB,highC,highD,highE,highF,highG
  USE mo_control,  ONLY: Rcal
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
    CHARACTER(30) :: Species(nSpc)
    CHARACTER(30) :: SpeciesInfo(nspc)
    CHARACTER(6)  :: RefDataCode(nSpc)
    CHARACTER(2)  :: Atoms(nSpc,4)
    INTEGER       :: nAtoms(nSpc,4)
    REAL(8)       :: TempRange(nSpc,3)
    REAL(8)       :: MolMass(nSpc)
    CHARACTER(1)  :: Phase(nSpc)
    REAL(8)       :: ta,tb,tc,td,te,tf,tg
    REAL(8)       :: H0_29815R(nspc)
    REAL(8),ALLOCATABLE :: ThermSwitchTemp(:)
    !
    INTEGER :: i,j,n, cnt
    INTEGER :: idxWhiteSpace
    INTEGER :: ALLOC_ERR
    CHARACTER(30) :: tSpcName
    LOGICAL :: THERMO=.FALSE.

    REAL(RealKind) :: thin1, thin2, thin3
    REAL(RealKind) :: ThermoIntervall(3)


    CHARACTER(1) :: rNumber='-'
    CHARACTER(4) :: dstr='----'
    !
    ! FORMATS
    1 FORMAT(A18,A6,4(A2,I3),A1,2F10.3,F13.7,A1,I1) 
    2 FORMAT(5(E15.8),I1)
    4 FORMAT(4(E15.8),I1)
    !
    !
    CALL OpenFile(UnitThermo,DatThermo(1:INDEX(DatThermo,'.')-1),'dat')
    ALLOCATE(lowA(nspc),lowB(nspc),lowC(nspc),lowD(nspc),lowE(nspc),lowF(nspc),lowG(nspc),STAT=ALLOC_ERR)
    ALLOCATE(highA(nspc),highB(nspc),highC(nspc),highD(nspc),highE(nspc),highF(nspc),highG(nspc),STAT=ALLOC_ERR)
    ALLOCATE(ThermSwitchTemp(nspc))
    lowA = -9999999999.0d0
    REWIND(UnitThermo)
    !
    i=0
    cnt=0
    ! find thermo and read 3 temp values
    DO
      READ(UnitThermo,*) iLine
      IF ( MAXVAL(INDEX(iLine,(/'THERMO','thermo'/))) > 0 ) THEN
        READ(UnitThermo,*) thin1, thin2, thin3
        ThermoIntervall(1)  = REAL(thin1,RealKind)
        ThermoIntervall(2)  = REAL(thin2,RealKind)
        ThermoIntervall(3)  = REAL(thin3,RealKind)
        EXIT
      END IF
    END DO

    DO
      !
      READ(UnitThermo,'(A80)') iLine
      iLine = ADJUSTR(iLine)

      rNumber = '-'
      IF ( SCAN(iLine(80:80),'1234') > 0 ) THEN
        rNumber = iLine(80:80)
      END IF

      IF ( MAXVAL(INDEX(iLine,['END','end'])) > 0 ) EXIT
      !
      SELECT CASE (rNumber)
        CASE ('1')
          READ(iLine,1) tSpcName
          i=PositionSpeciesAll(tSpcName(1:INDEX(tSpcName,' ')-1))

          IF ( i <= 0 ) THEN
            READ(UnitThermo,*) iLine
            READ(UnitThermo,*) iLine
            READ(UnitThermo,*) iLine
          ELSE
            cnt=cnt+1
            !
            READ(iLine,1) tSpcName,RefDataCode(i),(Atoms(i,j),  &
            &       nAtoms(i,j),j=1,4),Phase(i),ta,tb,tc,n
            idxWhiteSpace=INDEX(tSpcName,' ')
            !
            WRITE(Species(i),*) tSpcName(1:idxWhiteSpace-1)
            WRITE(SpeciesInfo(i),*) tSpcName(idxWhiteSpace:)
            !
            TempRange(i,1)=REAL(ta,KIND=RealKind)
            TempRange(i,2)=REAL(tb,KIND=RealKind)
            TempRange(i,3)=REAL(tc,KIND=RealKind)
            ThermSwitchTemp(i)=TempRange(i,3)
            MolMass(i)=REAL(tc,KIND=RealKind)
          END IF
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
        !READ(iLine,4) td,te,tf,tg,ta,n
        READ(iLine,4) td,te,tf,tg,n
        lowD(i)=REAL(td,KIND=RealKind)
        lowE(i)=REAL(te,KIND=RealKind)
        lowF(i)=REAL(tf,KIND=RealKind)
        lowG(i)=REAL(tg,KIND=RealKind)
       ! H0_29815R(i)=REAL(ta,KIND=RealKind)
       CASE DEFAULT
        CONTINUE
      END SELECT
    END DO
    CALL CloseFile(UnitThermo,DatThermo(1:INDEX(DatThermo,'.')-1),'dat')
    IF ( cnt < nspc ) THEN
      WRITE(*,*) '    Some species are missing in thermodynamic data (*.dat)'
      DO i=1,nspc
        IF (lowA(i)==-9999999999.0d0) THEN
          WRITE(*,*) ' Species ',i,' = ',TRIM(y_name(i)),' is missing!'
        END IF
      END DO
      STOP 'mo_ckinput'
    END IF
  END SUBROUTINE Read_Thermodata
  !
  !
  SUBROUTINE Read_Elements(DataReac,UnitReac)
    CHARACTER(*) :: DataReac
    INTEGER :: UnitReac
    
    INTEGER, PARAMETER :: nMaxElements=130
    CHARACTER(80) :: iLine   !i =1,2,3,4
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
      READ(UnitReac,'(A80)') iLine
      !
      ! check if next line is end
      IF (iLine(1:3)=='END'.OR.iLine(1:3)=='end') THEN
        CALL CutCArray(Elements,tElem,i)
        EXIT
      END IF
      !
      IF (iLine(1:4)=='ELEM'.OR.iLine(1:4)=='elem') THEN
        i=0
        READ(UnitReac,'(A80)') iLine
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
    CHARACTER(100)     :: iLine   !i =1,2,3,4
    !
    CHARACTER(100)     :: headline  
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
      READ(UnitReac,'(A100)') iLine
      iLine = ADJUSTL(iLine)
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
    CHARACTER(80) :: iLine   !i =1,2,3,4
    CHARACTER(80) :: locString
    CHARACTER(80) :: headline  
    CHARACTER(80) :: ductstr
    CHARACTER(80) :: dummyString
    CHARACTER(2)  :: EqSign
    !
    INTEGER :: i
    INTEGER :: iReac
    INTEGER :: iWS
    INTEGER :: iNxtSpc
    INTEGER :: iM, iKlammerM
    INTEGER :: io_err
    INTEGER :: nRowThirdBodys
    !
    INTEGER :: iKl,iKr
    CHARACTER(10) :: auxiliary
    !
    INTEGER :: nduct, nEducts, nProducts
    INTEGER :: fPosPlus , fPosEq, fPosFw
    INTEGER :: idxDuctE(6),idxDuctP(6)
    REAL(RealKind) :: KoefDuctE(6), KoefDuctP(6)
    CHARACTER(20)  :: NamesDuctE(6),NamesDuctP(6)
    REAL(RealKind), ALLOCATABLE :: valDucts(:)
    !
    !INTEGEr :: TableNspc
    CHARACTER(15) :: units
    REAL(RealKind) :: tmpReal,tmpReal2,tmpReal3
    !
    LOGICAL :: bR
    LOGICAL :: nxtReac
    LOGICAL :: ende
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
    END SELECT
    !
    ! first count reations, ...
    iReac = 0
    DO
      ! read next line
      READ(UnitReac,'(A80)') iLine
     
      iLine=ADJUSTL(iLine)

      ! exit cond
      IF ( MAXVAL(INDEX(iLine,['END','end']))>0 )  EXIT

      IF   ( INDEX(iLine,'=>')>0 ) THEN
        iReac = iReac+1    ! irreversible reactions
        IF ( INDEX(iLine,'<=>')>0) THEN
          iReac = iReac+1  ! reversible reactions with arrows
        END IF
      END IF
      !
    END DO
    REWIND UnitReac
    nReak=iReac
    !print*, ' ges reak = ', ireac
    !
    ! ALLOCATE reaction system structur
    ALLOCATE(ReactionSystem(nReak))
    CALL FindSection(UnitReac,'reactions',headline)
    nxtReac   = .FALSE.
    ende      = .FALSE.
    iReac     = 0
    fPosPlus  = 0
    fPosEq    = 0
    fPosFw    = 0
    iLine       = ''
    dummyString = ''
    LocString   = ''

    ! gather first reaction line
    CALL NextLine(UnitReac,iLine,ende,nxtReac)

    READ_REACTION_MECHANISM: DO
      
      IF ( ende ) EXIT READ_REACTION_MECHANISM  ! exit reading mechnism
         
      REACTION_LINE: IF (nxtReac)  THEN

        !print*,
        !print*, 'DEBUG:: vor  reaktion lesen iLine = ',ireac,TRIM(iLine)

        iReac = iReac + 1
        
        ReactionSystem(iReac)%Type='GAS'    ! all gaseous
        
        fPosEq = INDEX(iLine,'<=')
        ! if   <=>   reaction
        IF (fPosEq>0) THEN
          bR=.TRUE.
          ReactionSystem(iReac+1)%Type = 'GAS'    ! all gaseous
          ReactionSystem(iReac+1)%bR   = .TRUE.    ! all gaseous
        ELSE 
          fPosFw = INDEX(iLine,'=>')
          IF (fPosFw>0) bR = .FALSE.
        END IF
        
        ! FIRST LINE OF REACTION  ---  Extract Arrhenius Coeff 
        ALLOCATE(ReactionSystem(iReac)%Constants(3))
        IF (bR) ALLOCATE(ReactionSystem(iReac+1)%Constants(3))

        LocString = ADJUSTL(iLine)

        DO i=3,1,-1
          LocString = ADJUSTR(LocString)
          iWS = INDEX(LocString,' ',.TRUE.)        ! find first whitespace <--
          READ(LocString(iWS:),*,IOSTAT=io_err) tmpReal

          IF (io_err == 0 ) THEN
            ReactionSystem(iReac)%Constants(i)   = REAL(tmpReal,KIND=RealKind)
            IF (bR) &
            ReactionSystem(iReac+1)%Constants(i) = REAL(tmpReal,KIND=RealKind)
          ELSE
            WRITE(*,*) ''
            WRITE(*,*) '   Something went wrong reading the Arrhenius parameter'
            CALL  PrintError(io_err,iReac,LocString)
          END IF
          LocString = TRIM(LocString(:iWS))
        END DO

        ! save reaction string
        IF (fPosEq>0) THEN
          fPosEq = INDEX(LocString,'<=')
          ReactionSystem(iReac)%Line1   = TRIM(ADJUSTL(LocString( 1:fPosEq-1) ))//' => '// &
          &                               TRIM(ADJUSTL(LocString( fPosEq+3:  )))
          ReactionSystem(iReac+1)%Line1 = TRIM(ADJUSTL(LocString( fPosEq+3:  )))//' => '// &
          &                               TRIM(ADJUSTL(LocString( 1:fPosEq-1 )))
        ELSE
          ReactionSystem(iReac)%Line1   = TRIM(ADJUSTL(LocString))
        END IF
     
        ! extract the constant type by checking the reaction if +m or (+M) appears, and cut M off
        CALL GetConstantType( iM , iKlammerM , LocString ,         &
          &                   ReactionSystem(iReac)%TypeConstant , &
          &                   ReactionSystem(iReac)%Factor         )

        IF (bR) THEN
          ! save reaction string
          ReactionSystem(iReac+1)%TypeConstant = &
          &                   ReactionSystem(iReac)%TypeConstant
          ReactionSystem(iReac+1)%Factor  = &
          &                   ReactionSystem(iReac)%Factor
        END IF

        !print*, 'DEBUG::mo_ckinput getkonst nach locstr= ',locstring

        IF      ( INDEX(iLine,'<=>')>0 ) THEN
          EqSign(1:2) = '<='
        ELSE IF ( INDEX(iLine,'=>')>0 ) THEN
          EqSign(1:2) = '=>'
        END IF
     
        ! count educts and products
        ! get spc numbers and stoechiom coefs
        CALL CountDooku( nEducts   , idxDuctE , KoefDuctE , NamesDuctE ,&
                       & nProducts , idxDuctP , KoefDuctP , NamesDuctP ,&
                       & LocString , TRIM(EqSign))

        ! place species, stoecho coefs and type in reactionsystem struct
        ALLOCATE(ReactionSystem(iReac)%Educt(nEducts))
        IF (bR) ALLOCATE(ReactionSystem(iReac+1)%Product(nEducts))
        
        DO i=1,nEducts
          ReactionSystem(iReac)%Educt(i)%Type   = 'GAS'                ! immer gas in chemkin?
          ReactionSystem(iReac)%Educt(i)%Species= NamesDuctE(i)
          ReactionSystem(iReac)%Educt(i)%Koeff  = KoefDuctE(i)
          IF (bR) THEN
            ReactionSystem(iReac+1)%Product(i)%Type   = 'GAS'                ! immer gas in chemkin?
            ReactionSystem(iReac+1)%Product(i)%Species= NamesDuctE(i)
            ReactionSystem(iReac+1)%Product(i)%Koeff  = KoefDuctE(i)
          END IF
        END DO

        ALLOCATE(ReactionSystem(iReac)%Product(nProducts))
        IF (bR) ALLOCATE(ReactionSystem(iReac+1)%Educt(nProducts))

        DO i=1,nProducts
          ReactionSystem(iReac)%Product(i)%Type   = 'GAS'                ! immer gas in chemkin?
          ReactionSystem(iReac)%Product(i)%Species= NamesDuctP(i)
          ReactionSystem(iReac)%Product(i)%Koeff  = KoefDuctP(i)
          IF (bR) THEN
            ReactionSystem(iReac+1)%Educt(i)%Type   = 'GAS'                ! immer gas in chemkin?
            ReactionSystem(iReac+1)%Educt(i)%Species= NamesDuctP(i)
            ReactionSystem(iReac+1)%Educt(i)%Koeff  = KoefDuctP(i)
          END IF
        END DO

        ! read next line after reaction line (specific reaction parameter)
        CALL NextLine(UnitReac,iLine,ende,nxtReac)
        IF ( nxtReac ) EXIT REACTION_LINE  ! => no extra parameter, next reaction
        IF ( ende ) EXIT READ_REACTION_MECHANISM  ! exit reading mechnism
        
        
        EXPLICITE_REVERSE_REACTION_COEF: IF (bR .AND.  &
          .NOT.ReactionSystem(iReac)%Factor == '$(+M)' ) THEN

          ! skip empty lines and comment lines
          IF ( MAXVAL(INDEX(iLine,['DUPLICATE','duplicate']))>0 ) THEN
            CALL NextLine(UnitReac,iLine,ende,nxtReac)
            IF ( nxtReac ) EXIT EXPLICITE_REVERSE_REACTION_COEF
            IF ( ende ) EXIT READ_REACTION_MECHANISM  ! exit reading mechnism
          END IF

          IF ( MAXVAL(INDEX(iLine,['REV','rev']))>0 ) THEN

            ReactionSystem(iReac+1)%bRexp = .TRUE.
            LocString = ADJUSTL(iLine)
            iKl = INDEX(LocString,'/')
            iKr = INDEX(LocString,'/',.TRUE.)
            
            IF (.NOT.ALLOCATED(ReactionSystem(iReac+1)%Constants)) &
              ALLOCATE(ReactionSystem(iReac+1)%Constants(3))
            
            ! get rid of 'REV' or 'rev' and slashes 
            LocString = ADJUSTL(LocString(ikl+1:ikr-1))

            READ(LocString,*,IOSTAT=io_err) tmpReal, tmpReal2, tmpReal3

            IF ( io_err == 0 ) THEN
              ReactionSystem(iReac+1)%Constants(1) = REAL(tmpReal, KIND=RealKind)
              ReactionSystem(iReac+1)%Constants(2) = REAL(tmpReal2,KIND=RealKind)
              ReactionSystem(iReac+1)%Constants(3) = REAL(tmpReal3,KIND=RealKind)
            ELSE
              CALL  PrintError(io_err,iReac,LocString)
            END IF
            IF (.NOT.nxtReac) CALL NextLine(UnitReac,iLine,ende,nxtReac)
            IF ( ende ) EXIT READ_REACTION_MECHANISM  ! exit reading mechnism

          ELSE
            EXIT EXPLICITE_REVERSE_REACTION_COEF
          END IF
        END IF EXPLICITE_REVERSE_REACTION_COEF
        
        ! print backreaction
        IF (bR) THEN
          ReactionSystem(iReac+1)%Line2='BackReaction'
          ReactionSystem(iReac+1)%bR=.TRUE.
        END IF
      END IF REACTION_LINE ! Reaction => or ( = or <=> )
      !
      !print*, 'DEBUG:: nach REACTION_LINE  iLine = ',TRIM(iLine)

      !==================================================
      ! IF AN    +m, (+m)     OR   +M, (+M) is involved
      !==================================================
      ! handle additional parameters for reverse reaction etc.
      IF (INDEX(iLine,'/')>0) THEN
        ! probably thrid body spc or further parameter
      ELSE
        IF (.NOT.nxtReac) CALL NextLine(UnitReac,iLine,ende,nxtReac)
      END IF
      IF ( ende ) EXIT READ_REACTION_MECHANISM  ! exit reading mechnism
      !print*, 'DEBUG:: nach reaktion lesen iLine = ',ireac,TRIM(iLine)
      !IF ( nxtReac ) EXIT REACTION_LINE  ! => no extra parameter, next reaction

      !
      THIRD_BODY_M: IF ( ReactionSystem(iReac)%Factor == '$(+M)' ) THEN

        !print*, ' debug :: THIRD_BODY_M :: ',LocString
        LocString = ADJUSTL(iLine)

        iKl = INDEX(LocString,'/')
        iKr = INDEX(LocString,'/',.TRUE.)
        IF ( iKl > 0 ) THEN
          auxiliary = ADJUSTL(LocString(1:iKl-1))
          LocString = ADJUSTL(LocString(iKl+1:iKr-1))
          ReactionSystem(iReac)%Line3 = TRIM(auxiliary)
          IF (bR) ReactionSystem(iReac+1)%Line3 = TRIM(auxiliary)
        ELSE
          ReactionSystem(iReac)%Line3 = ''
          IF (bR) ReactionSystem(iReac+1)%Line3 = ''
        END IF
        !
        !
        IF (ReactionSystem(iReac)%Line3=='low'.OR.  &
        &   ReactionSystem(iReac)%Line3=='LOW'      ) THEN

          ALLOCATE(ReactionSystem(iReac)%LowConst(3))
          ALLOCATE(ReactionSystem(iReac+1)%LowConst(3))

          READ(LocString,*,IOSTAT=io_err) tmpReal, tmpReal2, tmpReal3

          IF ( io_err == 0 ) THEN
            ReactionSystem(iReac)%LowConst(1) = REAL(tmpReal, KIND=RealKind)
            ReactionSystem(iReac)%LowConst(2) = REAL(tmpReal2,KIND=RealKind)
            ReactionSystem(iReac)%LowConst(3) = REAL(tmpReal3,KIND=RealKind)
            ReactionSystem(iReac+1)%LowConst  = ReactionSystem(iReac)%LowConst
          ELSE
            CALL  PrintError(io_err,iReac,LocString)
          END IF

        ELSEIF (ReactionSystem(iReac)%Line3=='high'.OR. &
        &       ReactionSystem(iReac)%Line3=='HIGH'     ) THEN
          ALLOCATE(ReactionSystem(iReac)%HighConst(3))
          ALLOCATE(ReactionSystem(iReac+1)%HighConst(3))

          READ(LocString,*,IOSTAT=io_err) tmpReal, tmpReal2, tmpReal3

          IF ( io_err == 0 ) THEN
            ReactionSystem(iReac)%HighConst(1) = REAL(tmpReal, KIND=RealKind)
            ReactionSystem(iReac)%HighConst(2) = REAL(tmpReal2,KIND=RealKind)
            ReactionSystem(iReac)%HighConst(3) = REAL(tmpReal3,KIND=RealKind)
            ReactionSystem(iReac+1)%HighConst  = ReactionSystem(iReac)%HighConst
          ELSE
            CALL  PrintError(io_err,iReac,LocString)
          END IF

        END IF
        
        ! skip empty lines and comment lines
        IF (.NOT.nxtReac) CALL NextLine(UnitReac,iLine,ende,nxtReac)
        IF ( nxtReac ) EXIT THIRD_BODY_M  ! => no extra parameter, next reaction
        IF ( ende ) EXIT READ_REACTION_MECHANISM  ! exit reading mechnism
        
        LocString = ADJUSTL(iLine)
        ! extract troe parameters 
        IF (MAXVAL(INDEX(LocString,['troe','TROE']))>0) THEN
          ALLOCATE(ReactionSystem(iReac)%TroeConst(4))
          ALLOCATE(ReactionSystem(iReac+1)%TroeConst(4))

          LocString=ADJUSTL(iLine)
          iKl=INDEX(LocString,'/')
          iKr=INDEX(LocString,'/',.TRUE.)

          READ(LocString(iKl:iKr),*,IOSTAT=io_err) tmpReal, tmpReal2, tmpReal3

          IF ( io_err == 0 ) THEN
            ReactionSystem(iReac)%TroeConst(1) = REAL(tmpReal, KIND=RealKind)
            ReactionSystem(iReac)%TroeConst(2) = REAL(tmpReal2,KIND=RealKind)
            ReactionSystem(iReac)%TroeConst(3) = REAL(tmpReal3,KIND=RealKind)
            ReactionSystem(iReac+1)%TroeConst  = ReactionSystem(iReac)%TroeConst
          ELSE
            CALL  PrintError(io_err,iReac,LocString)
          END IF

          IF (.NOT.nxtReac) CALL NextLine(UnitReac,iLine,ende,nxtReac)
          IF ( ende ) EXIT READ_REACTION_MECHANISM  ! exit reading mechnism
        END IF
      END IF THIRD_BODY_M
      !print*, 'DEBUG:: nach third body M         = ',ireac,TRIM(iLine)

      !   CASE IF +m OR +M is involved in reaction iReac
      THIRD_BODY_M2: IF ( ReactionSystem(iReac)%Factor == '$+M' .OR. &
        &                 ReactionSystem(iReac)%Factor == '$(+M)'    ) THEN

        nRowThirdBodys = 0

        IF ( nxtReac ) THEN
          ! all spezies gleichmaessig
          ReactionSystem(iReac)%TBextra = .TRUE.
          IF (bR) ReactionSystem(iReac+1)%TBextra = .TRUE.
        ELSE
          EXTRACT_THIRD_BODY_SPC: DO
  
            LocString=ADJUSTL(iLine)
            
            iKl=INDEX(LocString,'/')
            iKr=INDEX(LocString,'/',.TRUE.)
           
            ! extract 3rd body species and alpha values
            IF ( iKl>0 ) THEN
              
              nRowThirdBodys = nRowThirdBodys + 1
  
              !print*, 'DEBUG:: in schleife 3rd body iline= ',nRowThirdBodys,TRIM(LocString)
              
              ReactionSystem(iReac)%TBextra = .TRUE.
              CALL ExtractThirdBodys(  ReactionSystem(iReac)%TBidx   &
              &                      , ReactionSystem(iReac)%TBspc   &
              &                      , ReactionSystem(iReac)%TBalpha &
              &                      , LocString, nRowThirdBodys)
            END IF
  
            ! skip empty lines and comment lines
            IF (.NOT.nxtReac) CALL NextLine(UnitReac,iLine,ende,nxtReac)
            IF ( nxtReac .OR. ende ) EXIT EXTRACT_THIRD_BODY_SPC  ! => no extra parameter, next reaction
    
          END DO EXTRACT_THIRD_BODY_SPC
          IF (bR) THEN
            ReactionSystem(iReac+1)%TBextra = .TRUE.
            nRowThirdBodys = SIZE(ReactionSystem(iReac)%TBidx)
            ALLOCATE(ReactionSystem(iReac+1)%TBidx(nRowThirdBodys), &
            &        ReactionSystem(iReac+1)%TBspc(nRowThirdBodys), &
            &        ReactionSystem(iReac+1)%TBalpha(nRowThirdBodys))
            ReactionSystem(iReac+1)%TBidx   = ReactionSystem(iReac)%TBidx
            ReactionSystem(iReac+1)%TBspc   = ReactionSystem(iReac)%TBspc
            ReactionSystem(iReac+1)%TBalpha = ReactionSystem(iReac)%TBalpha
          END IF
        END IF

      END IF THIRD_BODY_M2
      !print*, 'DEBUG:: nach third body M2        = ',ireac,TRIM(iLine)

      IF (bR) THEN
        iReac = iReac + 1
        bR    = .FALSE.
      END IF
      
      !print*, ' LOOP number done => next   ',iReac
      !read(*,*)

    END DO READ_REACTION_MECHANISM                      ! next reaction
    CALL CloseFile(UnitReac,DataReac,'sys')
    !
    neq       = nReak
    nreakgas  = nReak
    ALLOCATE(ListGas2(ntGAS))
    CALL HashTableToList(ListGas,ListGas2)
    CALL SortList(ListGas2)
    CALL ListToHashTable(ListGas2,ListGas)
    !
  END SUBROUTINE Read_Reaction
  !
  !
  SUBROUTINE NextLine(UnitR,Line,ende,nxtReac)
    INTEGER, INTENT(INOUT) :: UnitR
    CHARACTER(*) :: Line
    LOGICAL :: ende
    LOGICAL :: nxtReac
    !
    ende = .FALSE.
    nxtReac = .FALSE.
    Line = ''
    
    ! skip empty lines and comment lines
    DO 
      IF ( Line(1:1)=='!' .OR. Line=='' ) THEN
        READ(UnitR,'(A80)') Line
        Line = ADJUSTL(Line)
      ELSE
        IF ( MAXVAL(INDEX(Line,['END','end']))>0 ) ende = .TRUE.
        IF ( MAXVAL(INDEX(Line,['<=','=>']))>0) nxtReac = .TRUE.
        EXIT
      END IF
    END DO
    !
  END SUBROUTINE NextLine
  !
  !
  SUBROUTINE NewLine(UnitReac,LocString,i,j)
    INTEGER :: UnitReac
    CHARACTER(80) :: iLine
    CHARACTER(80) :: LocString
    INTEGER, OPTIONAL :: i,j
    ! read next line
    ! extract troe parameters 
    READ(UnitReac,'(A80)') iLine
    LocString=ADJUSTL(iLine)
    IF (PRESENT(i).AND.PRESENT(j)) THEN
      i=INDEX(LocString,'/')
      j=INDEX(LocString,'/',.TRUE.)
    END IF
  END SUBROUTINE NewLine
  !
  !
  SUBROUTINE ExtractThirdBodys(indM,spcM,aM,Line,iRow)
    ! out:
    INTEGER,        ALLOCATABLE, INTENT(INOUT) :: indM(:)
    CHARACTER(*),   ALLOCATABLE, INTENT(INOUT) :: spcM(:)
    REAL(RealKind), ALLOCATABLE, INTENT(INOUT) :: aM(:)
    ! in:
    CHARACTER(*), INTENT(IN) :: Line
    INTEGER,      INTENT(IN) :: iRow
    ! temp:
    INTEGER,        ALLOCATABLE :: tindM(:)
    CHARACTER(20),  ALLOCATABLE :: tspcM(:)
    REAL(RealKind), ALLOCATABLE :: taM(:)

    CHARACTER(LEN(Line)) :: locLine
    INTEGER :: nSpc, kl1, kl2, io_err
    INTEGER :: presentSpc
    REAL(RealKind) :: tmp

    locLine = Line

    ! if there where already 3rd body in line befor
    IF ( iRow > 1 ) THEN
      nSpc  = SIZE(indM) 
      ALLOCATE (tindM(nSpc),tspcM(nSpc),taM(nSpc))
      tindM = indM 
      tspcM = spcM
      taM   = aM
    ELSE
      nSpc = 0
    END IF

    DO
      IF (locLine=='') EXIT
      ! two slashes = one 3rd body spc
      kl1=INDEX(locLine,'/')  ; locLine=ADJUSTL(locLine(kl1+1:))
      kl1=INDEX(locLine,'/')  ; locLine=ADJUSTL(locLine(kl1+1:))
      nSpc=nSpc+1
    END DO

    IF(ALLOCATED(indM)) DEALLOCATE(indM)
    IF(ALLOCATED(spcM)) DEALLOCATE(spcM)
    IF(ALLOCATED(aM))   DEALLOCATE(aM)

    ALLOCATE(indM(nSpc),spcM(nSpc),aM(nSpc))

    locLine=Line
    IF ( iRow > 1 ) THEN
      nSpc = SIZE(tindM)
      indM(1:nSpc) = tindM
      spcM(1:nSpc) = tspcM
      aM(1:nSpc)   = taM
    ELSE
      nSpc=0
    END IF

    DO
      IF ( locLine == '' ) EXIT
      nSpc  = nSpc + 1
      kl1   = INDEX(locLine,'/')
      kl2   = INDEX(locLine(kl1+1:),'/')
      
      indM(nSpc) = -1
      spcM(nSpc) = ADJUSTL(TRIM(locLine(1:kl1-1)))
      READ(locLine(kl1+1:kl1+kl2-1),*,IOSTAT=io_err) tmp

      IF ( io_err == 0 ) THEN
        aM(nSpc) = REAL(tmp,KIND=RealKind)
        locLine  = ADJUSTL(locLine(kl1+kl2+1:))
      ELSE
        WRITE(*,*) '   Error reading thirdbody species!'
        WRITE(*,*) '      Error code:      ',io_err
        WRITE(*,*) '      String failed:   ',locLine
        STOP '  Reading chemical system '
      END IF

    END DO
  END SUBROUTINE ExtractThirdBodys

!***************************************************************************************************
!***************************************************************************************************
!
!         FUNCTIONS FOR  CONVERTING MOLE FRACTION; MASS FRACTION; MOLAR CONCENTRATION 
  
  FUNCTION MassFr_to_MoleFr(MassFr) RESULT(MoleFr)
    REAL(RealKind), ALLOCATABLE :: MoleFr(:)
    REAL(RealKind), INTENT(IN)  :: MassFr(:)     ! Mass fraction 
    !TEMP
    REAL(RealKind) :: W     ! mean molecular weight of a mixture

    W = MeanMolecularWeight( MassFr=MassFr )
    MoleFr   = MassFr * W * rMW

  END FUNCTION MassFr_to_MoleFr

  
  FUNCTION MassFr_To_MoleConc(MassFr,rho) RESULT(MoleConc)
    REAL(RealKind), ALLOCATABLE :: MoleConc(:) ! Concentration  [mol/cm3]
    REAL(RealKind), INTENT(IN)  :: MassFr(:)   ! Mass fraction  [g/g]
    REAL(RealKind), INTENT(IN)  :: rho         ! Mass density   [g/cm3] 
  
    MoleConc = rho * MassFr * rMW 

  END FUNCTION MassFr_To_MoleConc

  
  FUNCTION MoleFr_to_MassFr(MoleFr) RESULT(MassFr)
    REAL(RealKind), ALLOCATABLE :: MassFr(:)  ! Mass fraction [g/g]
    REAL(RealKind), INTENT(IN)  :: MoleFr(:)  ! Mole fraction  [mol/mol]
    REAL(RealKind) :: W   ! mean molecular weight of a mixture
   
    W = MeanMolecularWeight( MoleFr=MoleFr ) 
    MassFr  = MW * MoleFr / W

  END FUNCTION MoleFr_to_MassFr


  FUNCTION MoleConc_to_MassFr(MoleConc) RESULT(MassFr)
    REAL(RealKind), ALLOCATABLE :: MassFr(:)  ! Mass fraction [g/g]
    REAL(RealKind), INTENT(IN)  :: MoleConc(:)  ! Mole fraction  [mol/cm3]
    REAL(RealKind) :: W   ! mean molecular weight of a mixture
   
    MassFr  = MW * MoleConc / SUM( MoleConc * MW )

  END FUNCTION MoleConc_to_MassFr


  FUNCTION MoleFr_To_MoleConc(MoleFr,rho,Press,Temp) RESULT(MoleConc)
    REAL(RealKind), ALLOCATABLE :: MoleConc(:)      ! Mole concentration
    REAL(RealKind), INTENT(IN)  :: MoleFr(:)        ! Mole fraction 
    REAL(RealKind), INTENT(IN), OPTIONAL  :: rho    ! [g/cm3]
    REAL(RealKind), INTENT(IN), OPTIONAL  :: Press  ! [dyn/cm2]
    REAL(RealKind), INTENT(IN), OPTIONAL  :: Temp   ! [K]

    REAL(RealKind)              :: W

    IF (PRESENT(rho)) THEN
      W = MeanMolecularWeight( MoleFr=MoleFr )
      MoleConc = MoleFr * rho / W
    ELSE IF (PRESENT(Press).AND.PRESENT(Temp)) THEN
      MoleConc = MoleFr * Press / (Rerg * SUM(MoleFr) * Temp) 
    ELSE
      WRITE(*,*) '  The function needs either the density value rho or pressure+temperature value! '
      STOP
    END IF

  END FUNCTION MoleFr_To_MoleConc
  
!
!***************************************************************************************************
!***************************************************************************************************

  FUNCTION Pressure(MoleConc,Temp) RESULT(P)
    REAL(RealKind) :: P           ! Pressure in [Pa]
    REAL(RealKind), INTENT(IN)  :: MoleConc(:) ! in [mol/cm3] 
    REAL(RealKind), INTENT(IN)  :: Temp        ! in [K]
                        
    P = SUM( MoleConc ) * Temp * Rerg * dyncm2_to_Pa ! Rerg in [erg/mol/K]
 
  END FUNCTION Pressure
  

  FUNCTION Density(MoleConc) RESULT(rho)
    REAL(RealKind) :: rho          ! in [kg/cm3]
    REAL(RealKind), INTENT(IN)  :: MoleConc(:)  ! in [mol/cm3]
    
    rho = kilo * SUM( MoleConc * MW ) 

  END FUNCTION Density
  

  FUNCTION MeanMolecularWeight(MassFr,MoleFr,MoleConc) RESULT(W)
    REAL(RealKind) :: W
    REAL(RealKind), INTENT(IN), OPTIONAL :: MassFr(:)
    REAL(RealKind), INTENT(IN), OPTIONAL :: MoleFr(:)
    REAL(RealKind), INTENT(IN), OPTIONAL :: MoleConc(:)

    IF      ( PRESENT(MassFr) ) THEN
      W = ONE / SUM( MassFr * rMW )
    ELSE IF ( PRESENT(MoleFr) ) THEN
      W = SUM( MoleFr * MW )
    ELSE IF ( PRESENT(MoleConc) ) THEN
      W = SUM( MoleConc * MW ) / SUM( MoleConc )
    ELSE
      WRITE(*,*) '  The function needs either mass fraction or mole fraction or mole concentration! '
      STOP
    END IF
    
  END FUNCTION MeanMolecularWeight

  !
  !is the mass average mixture specific  heat at constant volume,
  SUBROUTINE MassAveMixSpecHeat(cvmixture,dUdT,MassFr,MoleConc)
    !IN
    REAL(RealKind) :: dUdT(:)                 ! in [-]
    REAL(RealKind), OPTIONAL :: MassFr(:)     ! in [g/g]
    REAL(RealKind), OPTIONAL :: MoleConc(:)   ! in [mol/cm3]
    !OUT
    REAL(RealKind) :: cvmixture         ! in [J/kg/K]
    !TEMP
    REAL(RealKind) :: ravgConc          ! in [cm3/mol][mol/g]

    IF (PRESENT(MassFr)) THEN
      cvmixture = kilo * R * SUM( MassFr * dUdT * rMW )
    ELSE IF (PRESENT(MoleConc)) THEN
      ravgConc  = ONE / SUM( MoleConc * MW )  
      cvmixture = kilo * R * SUM( MoleConc * dUdT ) * ravgConc
    ELSE
      WRITE(*,*) '  The function needs either mass fraction or mole concentration! '
      STOP
    END IF
  END SUBROUTINE MassAveMixSpecHeat
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
    ! loop two times because +M could be on both sides
      idxM  = INDEX(String,'+M')
      idxM2 = INDEX(String,'(+M)')
      !
      IF (MAX(idxM,idxM2)>0) THEN
        IF (idxM2>0) THEN 
          ConstType='PRESSX'
          Factor='$(+M)'
          String(idxM2:idxM2+3)='    '
          nreakgtroe=nreakgtroe+1
        ELSE
          ConstType='TEMPX'
          Factor='$+M'
          String(idxM:idxM+1)='  '
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
    !print*, 'mockinput     string= ',string
    IF (SCAN(String(1:1),CKSetNumber)>0) THEN
      NextSpecies=SCAN(String,CKSetSpecies)            ! find the next letter
      READ(String(1:NextSpecies-1),*) Koeff
    ELSE
      NextSpecies=1
      Koeff=1.0d0
    END IF
  END SUBROUTINE ExtractKoeff
  !
  !
  SUBROUTINE CountDooku(nE,edIndex,edKoeff,edNames,nP,prIndex,prKoeff,prNames,inLine,EqSign)
    !IN
    CHARACTER(*), INTENT(IN) :: inLine
    CHARACTER(*), INTENT(IN) :: EqSign
    !OUT
    INTEGER        :: nE, nP
    INTEGER        :: edIndex(6),prIndex(6)
    REAL(RealKind) :: edKoeff(6),prKoeff(6)
    CHARACTER(20)  :: edNames(6),prNames(6)
    !TEMP
    CHARACTER(LEN(inLine)) :: Line
    INTEGER :: PosPlus,PosEq, lenEq
    REAL(RealKind) :: fac
    INTEGER :: NxtSpc
    !
    !
    nE = 0            ;      nP = 0
    edIndex = 0       ; prIndex = 0
    edKoeff = 0.0d0   ; prKoeff = 0.0d0
    edNames = 'none'  ; prNames = 'none'
    lenEq   = 2
    IF (EqSign=='<=') lenEq=3

    Line = ADJUSTL(inLine)

    ! educt side (left)
    DO
      Line    = ADJUSTL(Line)
      PosPlus = INDEX(Line,'+')
      PosEq   = INDEX(Line,EqSign)

      IF ( PosEq > 0 ) THEN
        nE = nE + 1
        IF ( PosPlus>0 .AND. PosPlus<PosEq ) THEN
          CALL ExtractKoeff(NxtSpc,fac,Line(1:PosPlus-1))
          edIndex(nE) = PositionSpeciesAll(Line(NxtSpc:PosPlus-1))
          edKoeff(nE) = fac
          edNames(nE) = Line(NxtSpc:PosPlus-1)
          Line = ADJUSTL(Line(PosPlus+1:))
        ELSE
          ! last educt
          CALL ExtractKoeff(NxtSpc,fac,Line(:PosEq-1))
          edIndex(nE) = PositionSpeciesAll(Line(NxtSpc:PosEq-1))
          edKoeff(nE) = fac
          edNames(nE) = Line(NxtSpc:PosEq-1)
          ! get rid of <=>, => or =
          Line = ADJUSTL(Line(PosEq+lenEq:))
        END IF
      ELSE
        EXIT
      END IF
    END DO
    
    ! product side (right)
    DO
      PosPlus = INDEX(Line,'+')
      nP = nP + 1
      IF ( PosPlus > 0 ) THEN
        CALL ExtractKoeff(NxtSpc,fac,Line(1:PosPlus-1))
        prIndex(nP) = PositionSpeciesAll(Line(NxtSpc:PosPlus-1))
        prKoeff(nP) = fac
        prNames(nP) = Line(NxtSpc:PosPlus-1)
        Line = ADJUSTL(Line(PosPlus+1:))
      ELSE
        CALL ExtractKoeff(NxtSpc,fac,Line)
        prIndex(nP) = PositionSpeciesAll(Line(NxtSpc:))
        prKoeff(nP) = fac
        prNames(nP) = Line(NxtSpc:)
        EXIT
      END IF
    END DO

  END SUBROUTINE CountDooku
  !
  !
  SUBROUTINE FindSection(UnitReac,SecName,dummy)
    CHARACTER(*) :: SecName
    INTEGER :: UnitReac
    !
    CHARACTER(80) :: dummy
    INTEGER :: io_err

    DO
      READ(UnitReac,'(A80)',IOSTAT=io_err) dummy
      dummy=ADJUSTL(dummy)
      IF (io_err/=0) WRITE(*,*) 'Did not find ',SecName
      IF (dummy(1:LEN(SecName))==SecName) EXIT
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
    !IF (Species(1:1)=='p') THEN
    !  CALL InsertHash(ListPartic,TRIM(ADJUSTL(Species)),ntPart)
    !  Type='Partic'
    !ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
    !  CALL InsertHash(ListAqua,TRIM(ADJUSTL(Species)),ntAqua)
    !  Type='Aqua'
    !ELSE IF (Species(1:1)=='s') THEN
    !  CALL InsertHash(ListSolid,TRIM(ADJUSTL(Species)),ntSolid)
    !  Type='Solid'
    !ELSE IF (Species(1:1)=='['.AND.LEN(TRIM(Species))<10.AND. &
    !  &      Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
    !  CALL InsertHash(ListNonReac,TRIM(ADJUSTL(Species)),ntkat)
    !  Type='Inert'
    !ELSE IF (Species(1:1)=='(') THEN
    !ELSE
      CALL InsertHash(ListGas,TRIM(ADJUSTL(Species)),ntGas)
      Type='Gas'
    !END IF
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
        IF (ALLOCATED(ReacSys(i)%TBidx)) THEN
          WRITE(*,'(A2,5X,I4,F6.2)') (' |', ReacSys(i)%TBidx(j),ReacSys(i)%TBalpha(j),j=1,SIZE(ReacSys(i)%TBidx))
        ELSE
          WRITE(*,'(A)') ' | no special 3rd bodys '
        END IF
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
        WRITE(*,'(A2,5X,I4,F6.2)') (' |', ReacSys(i)%TBidx(j),ReacSys(i)%TBalpha(j),j=1,SIZE(ReacSys(i)%TBidx))
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

  SUBROUTINE PrintError(io_err,reaction,String)
    INTEGER :: io_err, reaction
    CHARACTER(*) :: String

    WRITE(*,*) '   Error reading reaction line!'
    WRITE(*,*) '      Error code:      ',io_err
    WRITE(*,*) '      Number reaction: ',reaction
    WRITE(*,*) '      String failed:   ',String
    STOP '  Reading chemical system '
  END SUBROUTINE PrintError
  !
END MODULE mo_ckinput
