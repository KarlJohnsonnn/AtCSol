!=============================================================================
!===   Initialization of Multiphase Chemistry 
!===   (Concentrations, Emissions, Deposition Velocities, ...)
!=============================================================================
!
SUBROUTINE InitChem(nc)
   USE Kind_Mod
   USE Meteo_Mod
   USE mo_MPI
   USE mo_reac
   USE mo_met
   USE mo_mphys
   USE mo_control
   USE Chemsys_Mod
   !
   !=========================================================================
   !===  Read Parameters for the Chemistry          
   !=========================================================================
   !
   IMPLICIT NONE
   !
   INTEGER :: nc
   !
   !----------------------------------------------
   INTEGER :: ntap12,ntap13,ntap14  
   INTEGER :: i,  ie, jt, indspc
   !
   INTEGER :: ntImp1, ImpModDim
   INTEGER, ALLOCATABLE :: ImpInd1(:), ImpInd_I(:)
   REAL(RealKind), ALLOCATABLE :: FracImpac1(:),FracImpac_I(:,:)
   !
   !
   REAL(RealKind) :: pi43=4.0d0/3.0d0*PI
   INTEGER, ALLOCATABLE :: cutRO2(:)
   !
   REAL(RealKind) :: actLWC
   !----------------------------------------------
   !
   INTEGER :: ierr, ityp, pos, pos1, iPos
   REAL(RealKind) :: mm, alpha, nue_otemp, dg
   CHARACTER(240) :: ctline
   CHARACTER(1), PARAMETER :: dkreuz = '#'
   CHARACTER(3), PARAMETER :: ctnext = 'END'
   !
   CHARACTER(60) :: string = ''
   CHARACTER(20) :: streof = 'EOF'         
   !
   INTEGER :: ifind
   EXTERNAL   ifind
   !
   ! für pH Start
   REAL(RealKind) :: kappa
   REAL(RealKind), EXTERNAL :: pHValue
   !----------------------------------------------
   !
   INTEGER, PARAMETER ::  nbeg = 12, ntyp = 3   ! Strings for type search
                                                ! Gas, Aqua, Solid
   INTEGER, PARAMETER ::  nbeg_d = 7
 
   CHARACTER(20) ::  ctbeg(nbeg),ctend(nbeg)
   DATA  ctbeg(1) /'BEGIN_GAS'  /   ,ctend(1) /'END_GAS'/      &
   &    ,ctbeg(2) /'BEGIN_AQUA'/    ,ctend(2) /'END_AQUA'/     &
   &    ,ctbeg(3) /'BEGIN_SOLID'/   ,ctend(3) /'END_SOLID'/    &
   &    ,ctbeg(4) /'BEGIN_INITIAL'/ ,ctend(4) /'END_INITIAL'/  &
   &    ,ctbeg(5) /'BEGIN_AFRAC'/   ,ctend(5) /'END_AFRAC'/    &
   &    ,ctbeg(6) /'BEGIN_IMPAC'/   ,ctend(6) /'END_IMPAC'/    &
   &    ,ctbeg(7) /'BEGIN_MODE'/    ,ctend(7) /'END_MODE'/     &
   &    ,ctbeg(8) /'BEGIN_EMISS'/   ,ctend(8) /'END_EMISS'/    &
   &    ,ctbeg(9) /'BEGIN_DEPOS'/   ,ctend(9) /'END_DEPOS'/    &
   &    ,ctbeg(10) /'BEGIN_DIAG'/   ,ctend(10)/'END_DIAG'/     &
   &    ,ctbeg(11)/'BEGIN_OUT'/     ,ctend(11)/'END_OUT'/      &
   &    ,ctbeg(12)/'BEGIN_SUM' /    ,ctend(12)/'END_SUM'/
   CHARACTER(20) ::  ctbeg_d(nbeg_d),ctend_d(nbeg_d)
   DATA  ctbeg_d(1) /'BEGIN_DATAGAS'/  ,ctend_d(1) /'END_DATAGAS'/ &
   &    ,ctbeg_d(2) /'BEGIN_DATAQUA'/  ,ctend_d(2) /'END_DATAQUA'/ &
   &    ,ctbeg_d(3) /'BEGIN_DATASOL'/  ,ctend_d(3) /'END_DATASOL'/ &
   &    ,ctbeg_d(4) /'BEGIN_DATARO2'/  ,ctend_d(4) /'END_DATARO2'/ &
   &    ,ctbeg_d(5) /'BEGIN_DATARO2aq'/  ,ctend_d(5) /'END_DATARO2aq'/ &
   &    ,ctbeg_d(6) /'BEGIN_QTGAS'  /  ,ctend_d(6) /'END_QTGAS'  / &
   &    ,ctbeg_d(7) /'BEGIN_QTAQUA' /  ,ctend_d(7) /'END_QTAQUA' /  
   !
   !=========================================================================
   !===  Set Default Values
   !=========================================================================
   !
   dkmt = 1.0e+03          ! Standard Mass Trasfer Coefficient
   kBig = 1.0e+06          ! Constant for Equilibria
   !
   !=========================================================================
   !---  Allocation of Arrays and Setting of Standard Values
   !=========================================================================
   !
   !
   !---  ASCI output and diagnose
   ALLOCATE(FracOut(nspc,2))
   ALLOCATE(FracDiag(nspc,2))
   FracOut(:,:)  = 0
   FracDiag(:,:) = 0
   !
   !---  Initial cconcentrations, emissions and deposition
   ALLOCATE(y_iconc(nspc+ntkat))
   ALLOCATE(y_emi(nspc+ntkat))
   ALLOCATE(y_udepo(nspc+ntkat))
   ALLOCATE(y_c1(nspc+ntkat))
   ALLOCATE(y_c2(nspc+ntkat))
   !
   y_iconc = 1.d-20    
   y_emi   = 0.d0
   y_udepo = 0.d0
   y_c1    = 0.d0
   y_c2    = 0.d0
   !
   !--  Gas Phase Species ([H2O] ??)
   y_c1(1:ntGas)=5.d-6
   y_c2(1:ntGas)=5.d-5
   !
   !--  Read Species Names
   OPEN(UNIT=89,FILE=ADJUSTL(TRIM(ChemFile))//'.chem',STATUS='UNKNOWN')
   DO i=1,24
     READ(89,*)
   END DO

  !-------------------------------------------------------------------------
  !--- Species names      
   !-------------------------------------------------------------------------
   !
   !---  set indices for pH and water dissociation 
   Hp_ind   = 0
   OHm_ind  = 0
   aH2O_ind = 0
   !
   ALLOCATE(y_name(nspc+ntkat))
   !
   DO jt=1,ntGas
     READ(89,*)  y_name(jt)
   END DO
   DO jt=1,ntAqua
     READ(89,*)  y_name(ntGas+jt)
     IF (y_name(ntGas+jt) == 'Hp')   Hp_ind = jt + ntGas
     IF (y_name(ntGas+jt) == 'OHm') OHm_ind = jt + ntGas
   END DO
   DO jt=1,ntkat
     READ(89,*)  y_name(ntGas+ntAqua+jt)
     IF (y_name(ntGas+ntAqua+jt) == '[aH2O]') aH2O_ind = jt 
   END DO
   !
   IF (hp_ind == 0)  THEN
     IF (MPI_ID==0) WRITE(*,*) 'ReadChem...Warning: Hp  not in mechanism!' 
   END IF
   IF (ohm_ind == 0)  THEN
     IF (MPI_ID==0)  WRITE(*,*) 'ReadChem...Warning: OHm  not in mechanism!' 
   END IF
   IF (ah2o_ind == 0)  THEN
     IF (MPI_ID==0)  WRITE(*,*) 'ReadChem...Warning: aH2O not in mechanism!' 
   END IF
   CLOSE(89)
   !
   DO i=1,ntkat
      IF (y_name(nspc+i) /= '[aH2O]')  THEN
         y_c1(nspc+i) = 5.d-6
         y_c2(nspc+i) = 5.d-5
      END IF
   END DO
   !
   !=========================================================================
   !===  Split Species Names
   !=========================================================================
   !
   IF (ntGas >= 1)   ALLOCATE (GasName(ntGas))          ! gas phase species
   IF (ntAqua >= 1)  ALLOCATE (AquaName(ntAqua))        ! aqueous phase species
   IF (ntkat >= 1)   ALLOCATE (PassName(ntkat))         ! passive phase species

   DO jt=1,ntGas
      GasName(jt) = ADJUSTL(y_name(jt))
   END DO
   !
   DO jt=1,ntAqua
      AquaName(jt) = ADJUSTL(y_name(ntGas+jt))
   END DO
   !
   DO jt=1,ntkat
      PassName(jt) = ADJUSTL(y_name(ntGas+ntAqua+jt))
   END DO
   !
   !=========================================================================
   !===  Set  Chemical DATA  
   !===  (Molar Mass, Charges, Densities) 
   !=========================================================================
   !
   !---  Allocate arrays
   ALLOCATE (Charge(ntAqua))             ! charge of ions
   ALLOCATE (Ladung(ntAqua))             ! charge of ions
   ALLOCATE (SolubInd(ntAqua))           ! solubility
   ALLOCATE (MolMass(ntAqua))            ! molar mass of species
   ALLOCATE (SpcDens(ntAqua))            ! density of species
   ALLOCATE (OrgIndex(ntAqua))           ! carbon atoms
   ALLOCATE (CC(ntAqua))                 ! compound class
   ALLOCATE (ActIndex(ntAqua))           ! index for calculation of activity coefficient

   !---  Set default values
   Charge(:)   = 0.d0
   Ladung(:)   = 0.d0
   SolubInd(:) = 0.d0
   MolMass(:)  = 0.d0
   SpcDens(:)  = 1.d0
   OrgIndex(:) = 0.d0
   CC(:) = '  '
   ActIndex(:) = 0.d0

   !--------------------------------------------------------------
   !---  Determine Charges of Ions
   DO jt=1,ntAqua
      string = AquaName(jt)

   !--- cations
      pos1 = 1
      DO i=1,5
         pos = index(string(pos1:20),'p')
         IF (pos <= 0)  EXIT 
         Charge(jt) = Charge(jt) + 1
         pos1 = pos1 + pos
      END DO

   !--- anions
      pos1 = 1
      DO i=1,5
         pos = index(string(pos1:20),'m')
         IF (pos <= 0)  EXIT 
         Charge(jt) = Charge(jt) - 1
         pos1 = pos1 + pos
      END DO
   END DO
   Ladung=Charge
   !
   !=======================================================================
   !===  Read Species Properties Parameters  (Tape13 = InitFile)
   !=======================================================================
   !
   ntap13 = DataUnit
   OPEN(unit=ntap13,file=TRIM(DataFile),status='old')
   !
   !----------------------------------------------------------------------------
   !---Gas phase species
   !
   call PosTape(ntap13,ctbeg_d(1),ctend_d(1),ierr)
   IF (ierr == 0)  THEN
   !
 GAS: DO 
         READ(ntap13,'(A240)',END=77) ctline
         ctline = ADJUSTL(ctline)                          ! spaces
         IF (ctline(1:1) == dkreuz)  CYCLE GAS            ! comment
         IF (ctline(1:3) == ctnext)  EXIT                 ! end of par loop
         pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
         !
         string = ctline(1:pos-1)
         indspc = ifind(ntGas, string, GasName)
         IF (indspc > 0)  THEN
            READ(ctline(pos+1:),*)  mm,alpha,dg
            if (alpha == 0.d0 .and. dg == 0.d0) CYCLE GAS 

            y_c1(indspc) = 1.d-12 / (3.d0*dg)              ! 1.term k_t
            nue_oTemp    = SQRT(8.d3 * 8.313d0 / Pi / mm) ! constant nu
            y_c2(indspc) = 4.d-06 / (3.d0*alpha*nue_oTemp) ! 2.term k_t
         ELSE

!--- non-reactive gas phase  species
            indspc = ifind(ntKat, string, PassName)
            IF (indspc > 0)  THEN
               indspc = indspc + ntGas + ntAqua
               READ(ctline(pos+1:),*)  mm,alpha,dg
               if (alpha == 0.d0 .and. dg == 0.d0) CYCLE GAS 

               y_c1(indspc) = 1.d-12 / (3.d0*dg)              ! 1.term k_t
               nue_oTemp    = SQRT(8.d3 * 8.313d0 / Pi / mm) ! constant nu
               y_c2(indspc) = 4.d-06 / (3.d0*alpha*nue_oTemp) ! 2.term k_t
            ELSE
               !WRITE(*,*) 'InitChem...Warning (Constants): Species '    &
               !&                ,string(1:pos-1),' NOT in Reaction System!'
            END IF

         END IF
      END DO  GAS
   END IF
   !
   !----------------------------------------------------------------------------
   !---Gas phase species (RO2)
   !
   call PosTape(ntap13,ctbeg_d(4),ctend_d(4),ierr)
   IF (ierr == 0)  THEN
   !
   i = 0
 PEROXY: DO 
         READ(ntap13,'(A240)',END=77) ctline
         ctline = ADJUSTL(ctline)                         ! spaces
         IF (ctline(1:1) == dkreuz)  CYCLE PEROXY            ! comment
         IF (ctline(1:3) == ctnext)  EXIT                 ! end of par loop
         IF (ctline(1:5) == 'nRO2:')  THEN                ! number of RO2 in MCM-mechanism
            READ(ctline(INDEX(ctline,':')+1:),*) nRO2
            ALLOCATE (cutRO2(nRO2))
            cutRO2=0
            CYCLE
         ENDIF

         pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
         !
         string = ctline(1:pos-1)
         indspc = ifind(ntGas, string, GasName)
         IF (indspc > 0)  THEN  
            i = i + 1       
            cutRO2(i) = indspc
         ELSE
            !WRITE(*,*) 'InitChem...Warning (PEROXY ==> RO2_Gas): Species '   &
         !&             ,string(1:pos-1),' NOT in Reaction System!'
         END IF
      END DO  PEROXY
      nRO2 = i
   END IF
   ALLOCATE(RO2(nRO2))
   RO2=cutRO2(1:nRO2)
   IF (ALLOCATED(cutRO2)) DEALLOCATE(cutRO2)

   !
   !----------------------------------------------------------------------------
   !--  Aqueous Phase species 
   !
   call PosTape(ntap13,ctbeg_d(2),ctend_d(2),ierr)
   IF (ierr == 0)  THEN

AQUA: DO 
         READ(ntap13,'(A240)',END=77) ctline
         ctline = ADJUSTL(ctline)                        ! spaces
         IF (ctline(1:1) == dkreuz)  CYCLE AQUA          ! comment
         IF (ctline(1:3) == ctnext)  EXIT                ! end

         pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
         !
         string = ctline(1:pos-1)
         indspc = ifind(ntAqua, string, AquaName)
         IF (indspc > 0)  THEN
            IF (mAcoeff >= 1)  THEN
               READ(ctline(pos+1:),*)  MolMass(indspc), Charge(indspc),     &
         &                             SolubInd(indspc), OrgIndex(indspc),  &
         &                             ActIndex(indspc), CC(indspc)
            ELSEIF (mAcoeff == 0)  THEN
               READ(ctline(pos+1:),*)  MolMass(indspc), Charge(indspc),     &
         &                             SolubInd(indspc), OrgIndex(indspc)!,  &
         !&                             CC(indspc)
            END IF
         ELSE
            !WRITE(*,*) 'InitChem...Warning (AquaFrac): Species '   &
         !&             ,string(1:pos-1),' NOT in Reaction System!'
         END IF
      END DO  AQUA
!      DO i = 1, ntAqua
!        WRITE(*,'(A2,",")',ADVANCE='NO') CC(i)
!      ENDDO
   !
   !--  Check initialization of aqueous species
      DO jt=1,ntAqua
         IF (MolMass(jt) <= 1.E-3)  THEN
            !WRITE(*,*) 'InitChem...Warning (Aqua-Data): Species '   &
         !&             ,TRIM(ADJUSTL(AquaName(jt))),' NOT initialized!'
         !   WRITE(*,*) '   MolMass = ',MolMass(jt),'    Charge = ',Charge(jt),       &
         !&             '   SolubInd = ', SolubInd(jt), ' OrgIndex = ', OrgIndex(jt), &
         !&             ' Compound class = ', CC(jt), '   ActIndex = ', ActIndex(jt)
         END IF
      END DO 
   END IF
   !
   !----------------------------------------------------------------------------
   !---Aqueous phase species (RO2)
   !
   call PosTape(ntap13,ctbeg_d(5),ctend_d(5),ierr)
   IF (ierr == 0)  THEN
   !
   i = 0
 ROtwo: DO 
         READ(ntap13,'(A240)',END=77) ctline
         ctline = ADJUSTL(ctline)                         ! spaces
         IF (ctline(1:1) == dkreuz)  CYCLE ROtwo          ! comment
         IF (ctline(1:3) == ctnext)  EXIT                 ! end of par loop
         IF (ctline(1:7) == 'nRO2aq:')  THEN              ! number of RO2 in CAPRAM-mechanism
            READ(ctline(INDEX(ctline,':')+1:),*) nRO2aq
            ALLOCATE (cutRO2(nRO2aq))
            cutRO2=0
            CYCLE
         ENDIF

         pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
         !
         string = ctline(1:pos-1)
         indspc = ifind(ntAqua, string, AquaName)
         !print*, 'debug:: aq ind spc', indspc
         IF (indspc > 0)  THEN  
            i = i + 1       
            cutRO2(i) = indspc
         ELSE
            !WRITE(*,*) 'InitChem...Warning (ROtwo ==> RO2_Aqua): Species '   &
         !&             ,string(1:pos-1),' NOT in Reaction System!'
         END IF
      END DO  ROtwo
      nRO2aq = i
      ALLOCATE (RO2aq(nRO2aq))
      RO2aq=cutRO2(1:nRO2aq)
      DEALLOCATE(cutRO2)
   END IF
   !
77 close (ntap13)   
   !
   !
   !=======================================================================
   !===  Read Initial Scenario Data  (Tape12 = InitFile)
   !=======================================================================
   !
   ntap12 = InitUnit
   OPEN(unit=ntap12,file=TRIM(InitFile),status='old')
   DO i=1,10 
      READ(ntap12,*)    ! 10 line for the header of the inifile (begin blank)
   END DO
   !
   !---  Units for Chemistry
   READ(ntap12,*)  GasUnit 
   READ(ntap12,*)  AquaUnit
   !
   IF (GasUnit == GasRateUnit)  THEN
      GasFac = 1.e0
   ELSE IF (GasUnit > GasRateUnit)  THEN
      GasFac = mol2part
   ELSE IF (GasUnit < GasRateUnit)  THEN
      GasFac = 1.e0 / mol2part
   END IF
   !
   ConvGas = 1.e0 / mol2part
   IF (GasRateUnit == 1)  THEN
      ConvGas = 1.e0 
      ConvAir = 1.e0 / mol2part
   ELSE
      ConvGas = 1.e0 / mol2part
      ConvAir = 1.e0 
   END IF
   !
   !
   !WRITE(*,*) 'edit InitChem'
   !WRITE(*,*) 'gasunit','aquaunit','gasrateunit'
   !WRITE(*,*) gasunit,aquaunit,gasrateunit
   !STOP
   !
   !=========================================================================
   !===  READ  Chemical DATA  
   !===  (Initial Concentration, Emissions, Deposition Velocities, 
   !===   Species Constants)
   !=========================================================================
   !
TYP: DO ityp=1,ntyp
      string = ctbeg(ityp)
      call PosTape(ntap12,ctbeg(ityp),streof,ierr)        !begin_gas and so on...
   !
      IF (ierr /= 0)  THEN
         !WRITE(*,*) 'InitChem: No Changes of ',string(7:11), ' in Input!'
         CYCLE TYP
      END IF
   !
   !--------------------------------------------------------------
   !-- Initial Concentrations (Gas Phase or Mass Fraction)
      call PosTape(ntap12,ctbeg(ntyp+1),ctend(ityp),ierr)
      IF (ierr == 0)  THEN
   !
     INI: DO 
            READ(ntap12,'(A240)',END=99) ctline
            ctline = ADJUSTL(ctline)                         ! spaces
            IF (ctline(1:1) == dkreuz)  CYCLE INI            ! comment
            IF (ctline(1:3) == ctnext)  EXIT                 ! end
            pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
            !
            string = ctline(1:pos-1)
            indspc = ifind(nspc+ntkat, string, y_name)
            IF (indspc > 0)  THEN
               READ(ctline(pos+1:),*)  y_iconc(indspc)
            ELSE
               !WRITE(*,*) 'InitChem...Warning (Initial): Species '   &
            !&             ,string(1:pos-1),' NOT in Reaction System!'
            END IF
         END DO  INI
   !
      END IF
      call PosTape(ntap12,ctbeg(ityp),streof,ierr)
   !
   !--------------------------------------------------------------
   !-- Emissions
      call PosTape(ntap12,ctbeg(ntyp+5),ctend(ityp),ierr)
      IF (ierr == 0)  THEN
   !
     EMI: DO 
            READ(ntap12,'(A240)',END=99) ctline
            ctline = ADJUSTL(ctline)                         ! spaces
            IF (ctline(1:1) == dkreuz)  CYCLE EMI            ! comment
            IF (ctline(1:3) == ctnext)  EXIT                 ! end
            pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
            !
            string = ctline(1:pos-1)
            indspc = ifind(nspc+ntkat, string, y_name)
            IF (indspc > 0)  THEN
               READ(ctline(pos+1:),*)  y_emi(indspc)
            ELSE
               !WRITE(*,*) 'InitChem...Warning (Emission): Species '  &
            !&             ,string(1:pos-1),' NOT in Reaction System!'
            END IF
         END DO  EMI
   !
      END IF
      call PosTape(ntap12,ctbeg(ityp),streof,ierr)
   !
   !--------------------------------------------------------------
   !-- Deposition Velocities
      call PosTape(ntap12,ctbeg(ntyp+6),ctend(ityp),ierr)
      IF (ierr == 0)  THEN
   !
     DEP: DO 
            READ(ntap12,'(A240)',END=99) ctline
            ctline = ADJUSTL(ctline)                         ! spaces
            IF (ctline(1:1) == dkreuz)  CYCLE DEP            ! comment
            IF (ctline(1:3) == ctnext)  EXIT                 ! end
            pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
            !
            string = ctline(1:pos-1)
            indspc = ifind(nspc, string, y_name)
            IF (indspc > 0)  THEN
               READ(ctline(pos+1:),*)  y_udepo(indspc)
            ELSE
               !WRITE(*,*) 'InitChem...Warning (Deposition): Species '  &
            !&             ,string(1:pos-1),' NOT in Reaction System!'
            END IF
         END DO  DEP
   !
      END IF
      call PosTape(ntap12,ctbeg(ityp),streof,ierr)
   !
   !--------------------------------------------------------------
   !-- Diagnose Species
       call PosTape(ntap12,ctbeg(ntyp+7),ctend(ityp),ierr)
       IF (ierr == 0)  THEN
   !
     DIA: DO 
            READ(ntap12,'(A240)',END=99) ctline
            ctline = ADJUSTL(ctline)                         ! spaces
            IF (ctline(1:1) == dkreuz)  CYCLE DIA            ! comment
            IF (ctline(1:3) == ctnext)  EXIT                 ! end
            pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
            !
            string = ctline(1:pos-1)
            indspc = ifind(nspc, string, y_name)
            IF (indspc > 0)  THEN
               READ(ctline(pos+1:),*)  FracDiag(indspc,1:2)
            ELSE
               !WRITE(*,*) 'InitChem...Warning (Diagnose): Species '    &
            !&             ,string(1:pos-1),' NOT in Reaction System!'
            END IF
         END DO  DIA
   !
      END IF
      call PosTape(ntap12,ctbeg(ityp),streof,ierr)
   !
   !--------------------------------------------------------------
   !-- Output Species
        call PosTape(ntap12,ctbeg(ntyp+8),ctend(ityp),ierr)
        IF (ierr == 0)  THEN
   !
     OUT: DO 
            READ(ntap12,'(A240)',END=99) ctline
            ctline = ADJUSTL(ctline)                         ! spaces
            IF (ctline(1:1) == dkreuz)  CYCLE OUT            ! comment
            IF (ctline(1:3) == ctnext)  EXIT                 ! end
            pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
            !
            string = ctline(1:pos-1)
            indspc = ifind(nspc, string, y_name)
            IF (indspc > 0)  THEN
               READ(ctline(pos+1:),*)  FracOut(indspc,1:2)  
            ELSE
               WRITE(*,*) 'InitChem...Error (Output): Species '       &
            &             ,string(1:pos-1),' NOT in Reaction System!'
               STOP 'InitChem...Error: Wrong Species (Output)!!'
            END IF
         END DO  OUT
   !
      END IF
   END DO  TYP
   !
   !==============================================================
   !== Read Mass Fractions for the Whole Spectrum
   !==============================================================
   !-- 
   ntImp1 = 0
   call PosTape(ntap12,ctbeg(ntyp+2),streof,ierr)
   IF (ierr == 0)  THEN
      ALLOCATE(ImpInd1(ntAqua))
      ImpInd1(:) = 0
      ALLOCATE(FracImpac1(ntAqua))
      FracImpac1(:) = 0.d0

 AFRAC: DO
         READ(ntap12,'(A240)',END=99) ctline
         ctline = ADJUSTL(ctline)                         ! spaces
         IF (ctline(1:1) == dkreuz)  CYCLE AFRAC          ! comment
         IF (ctline(1:3) == ctnext)  EXIT                 ! end
         pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
         !
         string = ctline(1:pos-1)
         indspc = ifind(ntAqua, string, AquaName)

         IF (indspc > 0)  THEN
            ntImp1 = ntImp1 + 1
            ImpInd1(ntImp1) = indspc
            READ(ctline(pos+1:),*)  MolMass(ntImp1), Charge(ntImp1), &
         &                          SolubInd(ntImp1),FracImpac1(ntImp1)
         ELSE
            !WRITE(*,*) 'InitChem...Warning (AquaFrac): Species '   &
         !&             ,string(1:pos-1),' NOT in Reaction System!'
         END IF
      END DO  AFRAC
   END IF
   !==============================================================
   !== Read Mass Fractions from Impactor Measurements
   !==============================================================
   !
   !---  Input of Impactor Measurements
   nImpac = 0
   call PosTape(ntap12,ctbeg(ntyp+3),streof,ierr)
   IF (ierr == 0)  THEN
      READ(ntap12,*)  nImpac
      ALLOCATE(BdImpac(nImpac+1))
      READ(ntap12,*)  (BdImpac(i),i=1,nImpac+1)
      ALLOCATE(ImpInd_I(ntAqua))
      ImpInd_I(:) = 0
      ALLOCATE(FracImpac_I(ntAqua,nImpac))
      FracImpac_I(:,:) = 0.e0
   !
      ntImp = 0
IMPAC: DO 
         READ(ntap12,'(A240)',END=99) ctline
         ctline = ADJUSTL(ctline)                         ! spaces
         IF (ctline(1:1) == dkreuz)  CYCLE IMPAC          ! comment
         IF (ctline(1:3) == ctnext)  EXIT                 ! end
         pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
         !
         string = ctline(1:pos-1)
         indspc = ifind(ntAqua, string, AquaName)
         IF (indspc > 0)  THEN
            ntImp = ntImp + 1
            ImpInd_I(ntImp) = indspc
            READ(ctline(pos+1:),*)  MolMass(indspc), Charge(indspc), &
         &                          SolubInd(indspc),                &
         &                         (FracImpac_I(ntImp,i),i=1,nImpac)
         ELSE
            !WRITE(*,*) 'InitChem...Warning (Impactor): Species '   &
         !&             ,string(1:pos-1),' NOT in Reaction System!'
         END IF
      END DO  IMPAC

      ImpModDim = nImpac
   END IF
   !
   !==============================================================
   !== Read Mass Fractions from mode Measurements
   !==============================================================
   !
   !
   !---  Input of mode Measurements
   nMode = 0
   call PosTape(ntap12,ctbeg(ntyp+4),streof,ierr)
   IF (ierr == 0)  THEN
      IF (nImpac >= 1)  THEN
         WRITE(*,*) 'InitChem..Error: Mixing of Impactor and Mode Input!'
         WRITE(*,*) '                 Delete Either of Them!!'
         STOP       'InitChem..Error: Mixing of Impactor and Mode Input!'
      END IF
      READ(ntap12,*)  nMode
      ALLOCATE(ImpInd_I(ntAqua))
      ImpInd_I(:) = 0
      ALLOCATE(FracImpac_I(ntAqua,nMode))
      FracImpac_I(:,:) = 0.e0
   !
      ntImp = 0
 MODE:DO
         READ(ntap12,'(A240)',END=99) ctline
         ctline = ADJUSTL(ctline)                         ! spaces
         IF (ctline(1:1) == dkreuz)  CYCLE MODE           ! comment
         IF (ctline(1:3) == ctnext)  EXIT                 ! end
         pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
        !
         string = ctline(1:pos-1)
         indspc = ifind(ntAqua, string, AquaName)
         IF (indspc > 0)  THEN
            ntImp = ntImp + 1
            ImpInd_I(ntImp) = indspc
            READ(ctline(pos+1:),*)  MolMass(indspc), Charge(indspc), &
         &                          SolubInd(indspc),                &
         &                         (FracImpac_I(ntImp,i),i=1,nMode)
         ELSE
            !WRITE(*,*) 'InitChem...Warning (Mode): Species '   &
         !&             ,string(1:pos-1),' NOT in Reaction System!'
         END IF
      END DO  MODE

      ImpModDim = nMode
   END IF
   !
   !==============================================================
   !==  Output for Comparison with Impactor Measurements
   !==============================================================
   !
   nOutImp = 0
   !
99 close (ntap12)
   !
   !
   !=======================================================================
   !===  Set QtOutput Species 
   !=======================================================================
   !
   !== Set Default: All Species
   ALLOCATE (IndQtGas (ntGas))
   IndQtGas=0
   ALLOCATE (IndQtAqua(ntAqua))
   IndQtAqua=0
   IF (TRIM(ADJUSTL(QtListFile)) == 'WITHOUT')  THEN
      nQtGas  = ntGas 
      nQtAqua = ntAqua
 
      IndQtGas (:) = 1
      IndQtAqua(:) = 1
   ELSE
      nQtGas  = 0 
      nQtAqua = 0
 
      IndQtGas (:) = 0
      IndQtAqua(:) = 0

   !----------------------------------------------------------------------------
   !== Read Species List
      ntap14 = InitUnit
      OPEN(unit=ntap14,file=TRIM(ADJUSTL(QtListFile)),status='old')
   !
   !   WRITE(*,*) '==========================================================================='
   !   WRITE(*,*) 'InitCem:  QtList-File opend: ',TRIM(QtListFile),'   Unit =', ntap14
   !
   !---Gas phase species
   !
      call PosTape(ntap14,ctbeg_d(6),ctend_d(6),ierr)
      IF (ierr == 0)  THEN
   !
QtGAS: DO 
         READ(ntap14,'(A240)',END=78) ctline
         ctline = ADJUSTL(ctline)                          ! spaces
         IF (ctline(1:1) == dkreuz)  CYCLE QtGAS           ! comment
         IF (ctline(1:3) == ctnext)  EXIT                  ! end of par loop
         pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))    ! separator
         !
         string = ctline(1:pos-1)
         indspc = ifind(ntGas, string, GasName)
         IF (indspc > 0)  THEN
            nQtGas = nQtGas + 1
            IndQtGas(indspc) = 1
         ELSE
            !WRITE(*,*) 
            !WRITE(*,*) 'InitChem...Warning (QtGas): Species '   &
         !&             ,string(1:pos-1),' NOT in Reaction System!'
         END IF
       END DO  QtGAS
      END IF
   !
   !----------------------------------------------------------------------------
   !--  Aqueous Phase species 
   !
      call PosTape(ntap14,ctbeg_d(7),ctend_d(7),ierr)
      IF (ierr == 0)  THEN

QtAQUA: DO 
         READ(ntap14,'(A240)',END=78) ctline
         ctline = ADJUSTL(ctline)                        ! spaces
         IF (ctline(1:1) == dkreuz)  CYCLE QtAQUA          ! comment
         IF (ctline(1:3) == ctnext)  EXIT                ! end

         pos = MAX(INDEX(ctline,' '),INDEX(ctline,','))   ! separator
         !
         string = ctline(1:pos-1)
         indspc = ifind(ntAqua, string, AquaName)
         IF (indspc > 0)  THEN
            nQtAqua = nQtAqua + 1
            IndQtAqua(indspc) = 1
         ELSE
            !WRITE(*,*) 
            !WRITE(*,*) 'InitChem...Warning (QtAqua): Species '   &
         !&             ,string(1:pos-1),' NOT in Reaction System!'
         END IF
       END DO  QtAQUA
      END IF
   !
   END IF
   !
78 close (ntap14)   
   !WRITE(*,*) '   nQtGas = ',nQtGas, '    nQtAqua =', nQtAqua
   !WRITE(*,*) 'InitCem:  QtList-File closed: ',TRIM(QtListFile),'   Unit =', ntap14
   !WRITE(*,*) '==========================================================================='
   !
   !============================================================================================
   !============================================================================================
   !===  Set Initial Values
   !============================================================================================
   !============================================================================================
   !
   !--------------------------------------------------
   !---   Initial concentrations 
   !--------------------------------------------------
   !
   !--- set gas phase 
   !y(:,:) = 0.0d0
   !DO i=1,ntGas
   !  print*, i , y(:,i),'----',y_iconc(i) , GasFac
   !  y(1,i)=y_iconc(i)*GasFac
   !END DO

   !--- set initial aqueous phase concentrations
   ntAConc = 0
   DO i=1,ntAqua
      IF (y_iconc(ntGas+i) >= 1.d-30)  ntAConc = ntAConc + 1 
   END DO
   !
   ALLOCATE(AConcInd(ntAConc))
   AConcInd(:) = 0
   ALLOCATE(AquaConc(ntAConc))
   AquaConc(:) = 0.d0
   !
   jt = 0
   DO i=1,ntAqua
      IF (y_iconc(ntGas+i)>=1.d-30) THEN
         jt=jt+1
         AConcInd(jt)=i
         AquaConc(jt)=y_iconc(ntGas+i)
      END IF
   END DO
   IF (jt/=ntAConc) THEN
      WRITE(*,*)  'InitChem ... Error: Wrong ntAConc =',ntAConc,' /= ',jt
      STOP  'InitChem ... Error: Wrong NTACONC!'
   END IF
   !--------------------------------------------------
   !--   time-constant emissions
   !--------------------------------------------------
   !
   ALLOCATE(y_e(nc,nspc))
   y_e(:,:) = 0.0d0
   DO i=1,ntGas
      y_e(1,i) = y_emi(i) * GasFac
   END DO
   DO i=ntGas+1,ntGas+ntAqua
      y_e(1,i) = y_emi(i)
   END DO
   !
   !
   !--------------------------------------------------
   !--  passiv species                             
   !--------------------------------------------------
   ! Set maximum LWC level
   !LWClevelmax=3.0d-04 * 1.d-6 ![l_w/m3_air]
   !
   ! calculate actual LWC level
   actLWC=pseudoLWC(tAnf)
   !
   ! set aH2o value 
   aH2O=aH2OmolperL*actLWC*mol2part
   !
   ALLOCATE(ykat(nc,ntkat))
   ykat(:,:) = 0.0d0
   IF (ntkat >= 1)  THEN
      DO i=1,ntKat
        IF (y_name(nspc+i)=='[aH2O]') THEN
          ykat(:,ntkat) = y_iconc(nspc+i)
          aH2O=y_iconc(nspc+i)*actLWC*mol2part
        END IF
      END DO
   END IF
   !
   !--------------------------------------------------
   !-- Compute initial aqua concentrations
   !--------------------------------------------------
   !
   ! Read Microphysic
   CALL Read_SPEK(SPEK,InitFile)
   !
   !
   y_iconc(ntGas+1:)=1.0d-16
   y_iconc(1:ntGas)=y_iconc(1:ntGas)
   DO i=1,SIZE(ImpInd1)
     iPos=ImpInd1(i)
     IF (iPos>0) THEN
       !
       y_iconc(ntGas+iPos) = SPEK(1)%Number*1.0d+03      &  ! [#/m3]
       &                   * FracImpac1(i)               &  ! [g/g]
       &                   * (pi43*(SPEK(1)%Radius)**3)  &  ! [m3]
       &                   * SPEK(1)%Density             &  ! [kg/m3]
       &                   / MolMass(i)             ! 1/[kg/mol] 
       y_iconc(ntGas+iPos) = y_iconc(ntGas+iPos) * mol2part 
     END IF
   END DO
   !======================================================
   !---  Compute pH value and number of ions
   !---  Initial pH by charge balance 
   IF (pHSet >= 1.AND.ntAqua>0)  THEN
     Kappa = pHValue(y_iconc(ntGAs+1:ntGas+ntAqua))
     !print*, 'kappa=',kappa, hp_ind, ohm_ind
     IF (Kappa > 0.d0)  THEN
       y_iconc(Hp_ind)  = Kappa
     ELSE 
       y_iconc(OHm_ind) = y_iconc(OHm_ind) + y_iconc(Hp_ind) - Kappa
     END IF
   END IF
   !
END SUBROUTINE InitChem

!================================================================
!
!--------------------------------------------------
!---     Integer-Function  IFIND                ---        
!--------------------------------------------------
   INTEGER FUNCTION ifind(nspc, string, namen)
!
!    #### Suchen des Indizes fuer Spezi  ####
!
      IMPLICIT NONE
!
      INTEGER  nspc, indspc
      character(60) :: string, namen(nspc)
!
      ifind = -1

      DO indspc = 1,nspc
         IF (string.eq.namen(indspc))  THEN
            ifind = indspc
            exit
         END IF
      END DO
!
!--------------------------------------------------
   END FUNCTION ifind

!================================================================
!
    SUBROUTINE PosTape(unit, string, strend, ierr)
!
!-----------------------------------------------------------
!---  Search String within a File
!-----------------------------------------------------------

       IMPLICIT NONE
!
       INTEGER ::  unit           &  ! chanal number of tape
&                 ,ierr              ! error code 
!
       CHARACTER(20) :: string    &  ! searched string
&                      ,strend       ! end of search
!
!-- internal parameter
       INTEGER ::  Loop, pos
       CHARACTER(80) :: str80    
!
!-----------------------------------------------------------
!-- Initialization
!
       ierr   = 1
       str80  = ' '
       string = ADJUSTL(string)
       pos    = INDEX(string,' ') - 1

       IF ( pos <= 0 )  THEN
          ierr = -1                           ! string empty
          WRITE(*,*) 'PosTape: Empty String ',string,'!'
          RETURN
       ELSE
          string = string(1:pos)
       END IF

       strend = ADJUSTL(strend)
       IF (TRIM(strend) == 'EOF')  THEN
          Loop = 1
       ELSE
          Loop = 0
       END IF
!
!-----------------------------------------------------------
!-- Search Loop
!
10     DO  
          IF ( INDEX(str80,TRIM(string)) /= 0)  EXIT
          READ (unit,'(a80)',END=20,ERR=20)  str80
          IF ( INDEX(str80,TRIM(strend) ) >= 1)  RETURN
       END DO
       ierr = 0
       RETURN
!
!-----------------------------------------------------------
!-- Rewind
!
20     IF (Loop >= 1)  THEN
          REWIND (unit)
          Loop = Loop - 1
          GO TO 10
       END IF
       ierr = 2
!
       END SUBROUTINE PosTape
            
!=======================================================================

   INTEGER FUNCTION ReadLine(ntap,MetName,nDim,Values)
!
      USE Kind_Mod
!----------------------------------------------------------------
!---  Function for Reading and Analyzing one Input Line
!----------------------------------------------------------------
!
      IMPLICIT NONE

!---  External variables
      INTEGER :: ntap,nDim
      REAL(RealKind) :: Values(nDim)
      CHARACTER(20) :: MetName

!---  Internal variables
      INTEGER :: i, pos, ierr
      CHARACTER(120) :: ctLine, ctRest
      CHARACTER(1)  :: char
      CHARACTER(1), PARAMETER :: dkreuz = '#'
   
!----------------------------------------------------------------
!
!---  read line
      READ(ntap,'(A240)',ERR=999,END=999)  ctLine
      ctLine  = ADJUSTL(ctLine)
      READ(ctLine,'(A1)')  char

!---  check comment line
      IF (char == 'dkreuz')  THEN
         ierr = 0
      ELSE

!---  check END of record
         IF (ctLine(1:3) == 'END')  THEN
            ReadLine = 2
            RETURN
         END IF

!---  read name of values
         pos  = index(ctLine,' ') - 1
         IF (pos <= 0)  pos = index(ctLine,',') - 1
         IF (pos <= 0)  pos = index(ctLine,';') - 1
         IF (pos <= 0)  pos = index(ctLine,':') - 1
         MetName  = ctLine(1:pos)

!---  read name of values
         pos = pos + 1
         ctRest = ctLine(pos:)
         READ(ctRest,*,END=999,ERR=999)  (Values(i),i=1,nDim)
         ierr = 1
      END IF
      ReadLine = ierr
      RETURN

999   ReadLine = -1
!------------------------------------------------------------
   END FUNCTION ReadLine
       
