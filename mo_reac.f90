!==============================================================
!===  MODULE for the Description of Chemistry, Deposition
!===  and Emissions
!==============================================================
!
  MODULE mo_reac
    USE Kind_Mod
    USE mo_control

    CHARACTER(9) :: measure_gas_ph(2)=(/"molec/cm3", "mol/m3   "/)
    CHARACTER(9) :: measure_aqua_ph(1)=(/"mol/l"/)
    CHARACTER(12):: units(2)

!--------------------------------------------------------------
!---  Reaction Mechanism
!--------------------------------------------------------------
!
!---  reaction types
    TYPE def_para
      CHARACTER(8) :: name_type   ! name of the reaction type
      INTEGER      :: n_par       ! number of parameters for corresponding reaction type
      INTEGER      :: n_reac      ! number of corresponding reactions
      LOGICAL      :: act         ! true if n_reac>0
    END TYPE def_para

    TYPE val_para
      INTEGER,  ALLOCATABLE :: iR(:,:)
      REAL(dp), ALLOCATABLE :: vR(:,:)
    END TYPE val_para


    INTEGER, PARAMETER ::           &
    &            iPHOTABC   = 1  &
    &         ,  iPHOTMCM   = 2  &
    &         ,  iPHOTAB    = 3  &
    &         ,  iCONST     = 4  &
    &         ,  iTEMP      = 5  &
    &         ,  iTEMP1     = 6  &
    &         ,  iTEMP2     = 7  &
    &         ,  iTEMP3     = 8  &
    &         ,  iTEMP4     = 9  &
    &         ,  iASPEC1    = 10  &
    &         ,  iASPEC2    = 11  &
    &         ,  iASPEC3    = 12  &
    &         ,  iDCONST    = 13  &
    &         ,  iDTEMP     = 14  &
    &         ,  iDTEMP2    = 15  &
    &         ,  iDTEMP3    = 16  &
    &         ,  iDTEMP4    = 17  &
    &         ,  iDTEMP5    = 18  &
    &         ,  iMESKHIDZE = 19  &
    &         ,  iT1H2O     = 20  &
    &         ,  iS4H2O     = 21  &
    &         ,  iTROE      = 22  &
    &         ,  iTROEQ     = 23  &
    &         ,  iTROEF     = 24  &
    &         ,  iTROEQF    = 25  &
    &         ,  iTROEXP    = 26  &
    &         ,  iTROEMCM   = 27  &
    &         ,  iSPEC1     = 28  &
    &         ,  iSPEC2     = 29  &
    &         ,  iSPEC3     = 30  &
    &         ,  iSPEC4     = 31  &
    &         ,  iSPEC1MCM  = 32  &
    &         ,  iSPEC2MCM  = 33  &
    &         ,  iSPEC3MCM  = 34  &
    &         ,  iSPEC4MCM  = 35  &
    &         ,  iSPEC5MCM  = 36  &
    &         ,  iSPEC6MCM  = 37  &
    &         ,  iSPEC7MCM  = 38  &
    &         ,  iSPEC8MCM  = 39  &
    &         ,  iSPEC9MCM  = 40  &
    &         ,  iHOM1      = 41  &
    &         ,  iPHOTO     = 42  &
    &         ,  iPHOTO2    = 43  &
    &         ,  iPHOTO3    = 44  &
    &         ,  iEQUI      = 45  &
    &         ,  iSPECIAL   = 46
  
    INTEGER, PARAMETER :: nReacTypes = 46

    TYPE(val_para), DIMENSION(nReacTypes) :: reac_val

    TYPE(def_para), DIMENSION(nReacTypes) :: reac_par = &
                   (/def_para("PHOTABC",   3,  0, .FALSE.), &
                     def_para("PHOTMCM",   3,  0, .FALSE.), &
                     def_para("PHOTAB",    2,  0, .FALSE.), &
                     def_para("CONST",     1,  0, .FALSE.), &
                     def_para("TEMP",      3,  0, .FALSE.), &
                     def_para("TEMP1",     2,  0, .FALSE.), &
                     def_para("TEMP2",     2,  0, .FALSE.), &
                     def_para("TEMP3",     2,  0, .FALSE.), &
                     def_para("TEMP4",     2,  0, .FALSE.), &
                     def_para("ASPEC1",    2,  0, .FALSE.), &
                     def_para("ASPEC2",    3,  0, .FALSE.), &
                     def_para("ASPEC3",    2,  0, .FALSE.), &
                     def_para("DCONST",    2,  0, .FALSE.), &
                     def_para("DTEMP",     3,  0, .FALSE.), &
                     def_para("DTEMP2",    4,  0, .FALSE.), &
                     def_para("DTEMP3",    4,  0, .FALSE.), &
                     def_para("DTEMP4",    4,  0, .FALSE.), &
                     def_para("DTEMP5",    4,  0, .FALSE.), &
                     def_para("MESKHIDZE", 7,  0, .FALSE.), &
                     def_para("T1H2O",     2,  0, .FALSE.), &
                     def_para("S4H2O",     4,  0, .FALSE.), &
                     def_para("TROE",      4,  0, .FALSE.), &
                     def_para("TROEQ",     6,  0, .FALSE.), &
                     def_para("TROEF",     5,  0, .FALSE.), &
                     def_para("TROEQF",    7,  0, .FALSE.), &
                     def_para("TROEXP",    5,  0, .FALSE.), &
                     def_para("TROEMCM",  10,  0, .FALSE.), &
                     def_para("SPEC1",     2,  0, .FALSE.), &
                     def_para("SPEC2",     2,  0, .FALSE.), &
                     def_para("SPEC3",     6,  0, .FALSE.), &
                     def_para("SPEC4",     4,  0, .FALSE.), &
                     def_para("SPEC1MCM",  3,  0, .FALSE.), &
                     def_para("SPEC2MCM",  3,  0, .FALSE.), &
                     def_para("SPEC4MCM",  4,  0, .FALSE.), &
                     def_para("SPEC3MCM",  2,  0, .FALSE.), &
                     def_para("SPEC5MCM",  4,  0, .FALSE.), &
                     def_para("SPEC6MCM",  4,  0, .FALSE.), &
                     def_para("SPEC7MCM",  6,  0, .FALSE.), &
                     def_para("SPEC8MCM",  4,  0, .FALSE.), &
                     def_para("SPEC9MCM",  4,  0, .FALSE.), &
                     def_para("HOM1",      4,  0, .FALSE.), &
                     def_para("PHOTO",     4,  0, .FALSE.), &
                     def_para("PHOTO2",    4,  0, .FALSE.), &
                     def_para("PHOTO3",    4,  0, .FALSE.), &
                     def_para("EQUI",      1,  0, .FALSE.), &
                     def_para("SPECIAL",   0,  0, .FALSE.)  /)


!---  reaction structures
    TYPE reactant
       INTEGER :: i_spc
       REAL(dp) :: d_koef
    END TYPE reactant

    TYPE reaction
       CHARACTER(8):: str_class
       CHARACTER(12):: str_type
       CHARACTER(20), POINTER :: factor

       INTEGER :: type
       INTEGER :: anz_p
       INTEGER :: n_so
       INTEGER :: n_si
       INTEGER :: n_siso
       INTEGER :: n_so_a
       INTEGER :: n_si_a
       INTEGER :: n_block
       INTEGER :: diag_flag
       INTEGER :: nstr

       INTEGER, POINTER :: so(:)
       INTEGER, POINTER :: si(:)
       INTEGER, POINTER :: struct(:,:)
       INTEGER, POINTER :: reac_si_nr(:)
       INTEGER, POINTER :: reac_so_nr(:)

       REAL(dp), POINTER  :: last_rate(:,:)
       REAL(dp), POINTER  :: back_rate(:,:)
       REAL(dp), POINTER  :: v(:,:)
       REAL(dp), POINTER  :: w(:,:)
       REAL(dp), POINTER :: dparam(:)  

       LOGICAL :: odd = .FALSE.

       REAL(dp) :: fac_exp
       REAL(dp) :: fac_A             

       TYPE (reactant), POINTER :: reactant(:)

       TYPE (reaction), POINTER :: next
       TYPE (reaction), POINTER :: next_all
    END TYPE reaction

    TYPE reactype
      CHARACTER(80) :: name=''
      INTEGER       :: n=0
      LOGICAL       :: exist=.FALSE.
    END TYPE reactype

!--------------------------------------------------------------
!--   dimensions
    INTEGER :: nt=0           ! unused?
    INTEGER :: nFrac=0        ! number of aquatic droplett classes
		INTEGER :: nspc=0         ! Number of all species excluding katalytic/passive species
    INTEGER :: nr=0           ! number of all reactions, HENRY/DISS counts as two reactions
    INTEGER :: neq=0          ! number of all reactions, HENRY/DISS counts as two reactions
    INTEGER :: nPhases=0      ! number of involed species phases
    INTEGER :: nDIM=0         ! Dimension of ODE system
    INTEGER :: nDIMcl=0       ! Dimension of ODE system for classic linear algebra
    INTEGER :: nDIMex=0       ! Dimension of ODE system for extended linear algebra 
    REAL(dp) :: rNspc					! real value 1/nspc for error calcualtion

    ! number of species in each phase/state
    INTEGER :: ns = 0, ns_KAT = 0, ns_G_KAT = 0, ns_A_KAT = 0
    INTEGER :: ns_GAS = 0, ns_AQUA = 0, ns_SOLID = 0, ns_PARTI = 0

    ! number of each reaction class
    INTEGER :: nr_gas = 0, nr_aqua = 0, nr_henry = 0, nr_diss = 0, nr_solid = 0, nr_parti = 0, nr_micphys = 0
    INTEGER :: nr_special = 0 

    INTEGER :: nr_liquid = 0  ! = nr_diss+nr_aqua


    LOGICAL :: hasGasSpc=.FALSE. , hasAquaSpc=.FALSE. , hasSolidSpc=.FALSE. , hasPartiSpc=.FALSE.
    LOGICAL :: hasGasReac=.FALSE. ,   hasAquaReac=.FALSE. , hasSolidReac=.FALSE. , hasPartiReac=.FALSE. &
		&        , hasHenryReac=.FALSE. , hasDissReac=.FALSE. , hasMicroReac=.FALSE. , hasLiquidReac=.FALSE.

    LOGICAL :: hasPhotoReac=.FALSE. , hasFactorReac=.FALSE.

    ! number of each reaction type
    INTEGER :: nr_G_photo = 0, nr_G_const = 0, nr_G_temp = 0, nr_G_troe =0, nr_G_spec = 0, nr_G_lind = 0
    INTEGER :: nr_A_photo = 0, nr_A_const = 0, nr_A_temp = 0, nr_A_spec =0
    INTEGER :: nr_S_temp  = 0, nr_S_equi  = 0, nr_S_spec = 0
    INTEGER :: nr_SimpTB  = 0, nr_press   = 0
    INTEGER :: nr_G_special = 0, nr_A_special = 0,nr_S_special = 0,nr_P_special = 0,nr_M_special = 0
    INTEGER :: nr_H_special = 0, nr_D_special = 0, nr_D_Temp

    INTEGER :: nr_PHOTabc = 0, nr_PHOTab = 0, nr_PHOTmcm = 0, nr_CONST = 0
    INTEGER :: nr_TEMP = 0, nr_TEMP1 = 0, nr_TEMP2 = 0, nr_TEMP3 = 0, nr_TEMP4 = 0
    INTEGER :: nr_TROE  = 0, nr_TROEf = 0, nr_TROEq = 0, nr_TROEqf = 0, nr_TROExp = 0, nr_TROEmcm = 0
    INTEGER :: nr_SPEC1 = 0, nr_SPEC2 = 0, nr_SPEC3 = 0, nr_SPEC4 = 0
    INTEGER :: nr_SPEC1mcm = 0, nr_SPEC2mcm = 0, nr_SPEC3mcm = 0, nr_SPEC4mcm = 0
    INTEGER :: nr_SPEC5mcm = 0, nr_SPEC6mcm = 0, nr_SPEC7mcm = 0, nr_SPEC8mcm = 0, nr_SPEC9mcm = 0
    INTEGER :: nr_S4H2O  = 0, nr_T1H2O   = 0, nr_Meskhidze = 0
    INTEGER :: nr_ASPEC1 = 0, nr_ASPEC2  = 0, nr_ASPEC3   = 0, nr_ASPEC4 = 0
    INTEGER :: nr_DTEMP  = 0, nr_DTEMP2  = 0, nr_DTEMP3   = 0, nr_DTEMP4 = 0, nr_DTEMP5 = 0, nr_DCONST = 0
    INTEGER :: nr_HENRYga = 0, nr_HENRYag  = 0
    INTEGER :: nr_FACTOR = 0, nr_FAC_H2  = 0, nr_FAC_O2N2 = 0, nr_FAC_M  = 0, nr_FAC_O2 = 0, nr_FAC_N2 = 0
    INTEGER :: nr_FAC_H2O = 0, nr_FAC_RO2 = 0, nr_FAC_O2O2 = 0, nr_FAC_aH2O = 0, nr_FAC_RO2aq = 0
    INTEGER :: nr_HOaqua = 0      ! higher order aqueous reactions
    INTEGER :: nr_PHOTOkpp = 0, nr_PHOTO2kpp = 0, nr_PHOTO3kpp = 0
    INTEGER :: nr_HOM1 = 0

    LOGICAL :: PHOTO=.FALSE.

!    INTEGER :: nreakstemp,nreaksequi,nreaksspec
!--    define indices of special species
    INTEGER ::   Hp_ind,         &   ! Index Hp (aqua)
&               OHm_ind,         &   ! Index OHm  (aqua)
&               H2O_ind,         &   ! Index H2O (gas)
&              aH2O_ind,         &   ! Index aH2O (aqua)
&              SO4mm_ind,        &   ! Index SO4mm_ind (aqua)
&              HSO4m_ind,        &   ! Index HSO4m_ind (aqua)
&              Temp_ind              ! Index Temperatur

!--------------------------------------------------------------
!--    passive species, indices of Henry species
    REAL(dp), ALLOCATABLE :: ykat(:,:)
    INTEGER, ALLOCATABLE :: ind_henry(:,:)

!--    mass transfer coefficient, kBig
    REAL(dp) :: dkmt = 1.0d3      ! Standard Mass Trasfer Coefficient
    REAL(dp) :: kBig = 1.0d6    ! Constant for Equilibria
!
!--------------------------------------------------------------
!--    species names 
    CHARACTER(60), ALLOCATABLE :: y_name(:)
    CHARACTER(60), ALLOCATABLE :: GasName(:)    & ! gas phase species
&                                ,AquaName(:)   & ! aqueous phase species
&                                ,SolidName(:)  & ! solid species
&                                ,ParticName(:)  & ! solid species
&                                ,PassName(:)     ! passive species
!
!
!--------------------------------------------------------------
!-- Netcdf output and diagnostic variables 
    CHARACTER(60), ALLOCATABLE  :: Gas_Name_Netcdf(:)        &
&                                 ,Gas_LongName(:)           &
&                                 ,Noreac_Name_Netcdf(:)     &
&                                 ,Aqua_Name_Netcdf(:)       &
&                                 ,Aqua_LongName(:)          &
&                                 ,Sum_Aqua_Name_Netcdf(:,:)
!    CHARACTER(20), ALLOCATABLE  :: ReacName_Gas(:)  &
!&                                 ,ReacName_Hen(:)  &
!&                                 ,ReacName_Diss(:) &
!&                                 ,ReacName_Aqua(:) &
!&                                 ,ReacName_Sol(:)
!    INTEGER,  ALLOCATABLE       :: Gas_anz(:)  &
!&                                 ,Aqua_anz(:) &
!&                                 ,Hen_anz(:)  &
!&                                 ,Diss_anz(:) &
!&                                 ,Sol_anz(:)
!    INTEGER, PARAMETER          :: ndiag_gas = 11  &
!&                                 ,ndiag_aq = 16  
    CHARACTER(60), ALLOCATABLE  :: Diag_Name_Netcdf(:) &
&                                , Diag_LongName(:)    &
&                                , DiagERR_Name_Netcdf(:)

!    CHARACTER(60)               :: Diag_Name_Gas(ndiag_gas)  &
!&                                 ,Diag_Name_Aqua(ndiag_aq)  

    INTEGER, ALLOCATABLE  :: iDiag_Schwefel(:)  
!--------------------------------------------------------------
!--    Peroxyradicals
    LOGICAL              :: hasRO2
    INTEGER, ALLOCATABLE :: RO2(:)          ! Species-Index of peroxyradical
    INTEGER              :: nRO2            ! Number of Gasphase-Peroxyradicals in mechanism
!    REAL(dp), ALLOCATABLE :: SumRO2(:)       ! Summation of all gasphase peroxyradicals in everey cell
!--    aqueous phase
    LOGICAL              :: hasRO2aq
    INTEGER, ALLOCATABLE :: RO2aq(:)          ! Species-Index of peroxyradical
    INTEGER              :: nRO2aq            ! Number of aqueos phase peroxy radicals in mechanism
!    REAL(dp), ALLOCATABLE :: SumRO2aq(:)       ! Summation of all aqueopus phase peroxyradicals in everey cell
!
!--------------------------------------------------------------
!--    Aerosol species properties
    REAL(dp), ALLOCATABLE :: Charge(:)       ! charge of ions
    REAL(dp), ALLOCATABLE :: SolubInd(:)     ! solubility index
    REAL(dp), ALLOCATABLE :: MolMass(:)      ! molar mass of species
    REAL(dp), ALLOCATABLE :: SpcDens(:)      ! density of species
    REAL(dp), ALLOCATABLE :: OrgIndex(:)     ! carbon atoms
    CHARACTER(2), ALLOCATABLE :: CC(:)      ! compound class
    REAL(dp), ALLOCATABLE :: ActIndex(:)     ! Computation of Activity coefficients index 

!--    Input from impactor measurements
    INTEGER :: nImpac = 0                 & ! number of impactor stages
              ,nMode  = 0                 & ! number of input modes
              ,ntImp  = 0                   ! number of impactor species
    INTEGER, ALLOCATABLE :: ImpInd(:)       ! indices of impactor species
    REAL(dp), ALLOCATABLE :: BdImpac(:)      ! boundaries of stages
    REAL(dp), ALLOCATABLE :: FracImpac(:,:)  ! mass fractions per stage
    REAL(dp), ALLOCATABLE :: DryMod(:,:)     ! modal input masses per fractions

!!--    Input from impactor measurements
!    INTEGER :: ntAConc = 0                   ! number  of initialized concentrations
!    INTEGER, ALLOCATABLE :: AConcInd(:)      ! indices of initialized concentrations
!    REAL(dp), ALLOCATABLE :: AquaConc(:)      ! initial aqueous phase concentrations 
                                             ! after activation
!--------------------------------------------------------------
!--    input arrays
    !REAL(dp):: tcur
    REAL(dp), ALLOCATABLE :: InitValAct(:)
    REAL(dp), ALLOCATABLE :: InitValKat(:)

    REAL(dp), ALLOCATABLE :: henry_diff(:)
    REAL(dp), ALLOCATABLE :: henry_accom(:)

!--------------------------------------------------------------
!--    deposition and emissions
    REAL(dp), ALLOCATABLE :: vd(:)
    REAL(dp), ALLOCATABLE :: y_emi(:)
    REAL(dp), ALLOCATABLE :: y_depos(:)

!--------------------------------------------------------------
    
    CHARACTER(18), ALLOCATABLE :: ThSpecies(:)
    CHARACTER(6),  ALLOCATABLE :: ThRefDataCode(:)
    CHARACTER(2),  ALLOCATABLE :: ThAtoms(:)           ! (:,4)
    INTEGER,       ALLOCATABLE :: ThnAtoms(:,:)          ! (:,4)
    REAL(dp),      ALLOCATABLE :: ThTempRange(:,:)      ! (:,2)
    REAL(dp),      ALLOCATABLE :: ThMolMass(:)
    CHARACTER(1),  ALLOCATABLE :: ThPhase(:)
    REAL(dp),      ALLOCATABLE :: lowA(:),lowB(:),lowC(:),lowD(:),lowE(:),lowF(:),lowG(:)
    REAL(dp),      ALLOCATABLE :: highA(:),highB(:),highC(:),highD(:),highE(:),highF(:),highG(:)
    REAL(dp),      ALLOCATABLE :: H0_29815R(:)
    REAL(dp),      ALLOCATABLE :: SwitchTemp(:)
    !
    REAL(dp) :: y_total
    !
    !
    INTEGER :: loc_diagCnt, loc_rateCnt, loc_concCnt, glob_diagCnt
    INTEGER, ALLOCATABLE :: loc_ratePtr(:), loc_concPtr(:), glob_ratePtr(:)
    !
    INTEGER :: LowerRateLim, UpperRateLim
    !
    INTEGER, ALLOCATABLE :: first_ReacPtr(:)       ! array =(/ 1,2,3,4....,neq/)
    !
    LOGICAL :: combustion=.FALSE.                      ! flag for combustion mechanism
    !
    REAL(dp), ALLOCATABLE :: GFE(:), DGFEdT(:)
    REAL(dp), ALLOCATABLE :: DelGFE(:), DDelGFEdT(:)
    !
    REAL(dp), ALLOCATABLE :: MW(:)    ! molecular weight (combustion)
    REAL(dp), ALLOCATABLE :: rMW(:)    !1/ molecular weight (combustion)
    !
    ! more speedchem stuff
    REAL(dp) :: rho, rRho   ! rho = density, rRho=kilo/rho
    INTEGER, ALLOCATABLE :: SCperm(:)

    !
    ! indix arrays for different reaction types 
    ! Types for vectorised version
    TYPE ReacTypeIndex_CK
      INTEGER, ALLOCATABLE :: iArr(:),  iLind(:), iTroe(:)
      INTEGER, ALLOCATABLE :: iEqui(:), iXrev(:),  iTBody(:)
      INTEGER, ALLOCATABLE :: iTBodyExtra(:), iHigh(:), iLow(:)
      INTEGER :: nArr, nLind, nTroe
      INTEGER :: nEqui,nXrev, nTBody
      INTEGER :: nTBodyExtra, nHigh, nLow
    END TYPE ReacTypeIndex_CK

    TYPE(ReacTypeIndex_CK) :: RTind

    !
    ! indix arrays for different reaction types 
    ! Types for vectorised version
    TYPE ReacTypeIndex_TR
      INTEGER, ALLOCATABLE :: iPHOTabc(:),  iPHOTab(:),   iPHOTmcm(:), iCONST(:)  &
      &                     , iTEMP(:),    iTEMP1(:),    iTEMP2(:),    iTEMP3(:), iTEMP4(:)  &
      &                     , iTROE(:),     iTROEf(:),    iTROEq(:),   iTROEqf(:), iTROExp(:) &
      &                     , iTROEmcm(:),  iSPEC1(:),    iSPEC2(:),   iSPEC3(:),  iSPEC4(:)  &
      &                     , iSPEC1mcm(:), iSPEC2mcm(:), iSPEC3mcm(:), iSPEC4mcm(:)  &
      &                     , iSPEC5mcm(:), iSPEC6mcm(:), iSPEC7mcm(:), iSPEC8mcm(:), iSPEC9mcm(:)  &
      &                     , iS4H2O(:),    iT1H2O(:),    iASPEC1(:),   iASPEC2(:)    &
      &                     , iASPEC3(:),   iASPEC4(:),   iDCONST(:,:),   iDTEMP(:,:) &
      &                     , iDTEMP2(:,:), iDTEMP3(:,:), iDTEMP4(:,:), iDTEMP5(:,:), iMeskhidze(:,:)&
      &                     , iHENRY(:,:) , iHOM1(:) &
      ! nFACTOR x 2, with [H2,O2N2,M,O2,N2,H2O,RO2,O2O2,aH2O,RO2aq,+M/(+M)]
      &                     , iFAC_H2(:),  iFAC_O2N2(:), iFAC_M(:),    iFAC_O2(:),   iFAC_N2(:) &
      &                     , iFAC_H2O(:), iFAC_RO2(:),  iFAC_O2O2(:), iFAC_aH2O(:), iFAC_RO2aq(:) &
      &                     , iHOaqua(:) &
      &                     , iPHOTOkpp(:), iPHOTO2kpp(:), iPHOTO3kpp(:) &
      &                     , iSPECIAL(:)
      REAL(dp), ALLOCATABLE :: PHOTabc(:,:)   &  ! nPHOTabc x 3
      &                     ,  PHOTab(:,:)    &  ! nPHOTab x 2
      &                     ,  PHOTmcm(:,:)   &  ! nPHOTmcm x 3
      &                     ,  CONST(:)       &  ! nCONST x 1
      &                     ,  TEMP(:,:), TEMP1(:,:), TEMP2(:,:), TEMP3(:,:), TEMP4(:,:) & ! nTEMPi x 2
      &                     ,  TROE(:,:)      &  ! nTROE x 4
      &                     ,  TROEf(:,:)     &  ! nTROEf x 5
      &                     ,  TROEq(:,:)     &  ! nTROEq x 6
      &                     ,  TROEqf(:,:)    &  ! nTROEqf x 7
      &                     ,  TROExp(:,:)    &  ! nTROExp x 5
      &                     ,  TROEmcm(:,:)   &  ! nTROExp x 10
      &                     ,  SPEC1(:,:), SPEC2(:,:)  &! nSPEC1,2 x 2
      &                     ,  SPEC3(:,:)     &  ! nSPEC3 x 6
      &                     ,  SPEC4(:,:)     &  ! nSPEC4 x 4
      &                     ,  SPEC1mcm(:,:), SPEC2mcm(:,:)  &! nSPEC1mcm,2 x 3
      &                     ,  SPEC3mcm(:,:)   & ! nSPEC3mcm x 2
      &                     ,  SPEC4mcm(:,:), SPEC5mcm(:,:), SPEC6mcm(:,:), SPEC8mcm(:,:) & ! nSPEC4,5,6mcm x 4
      &                     ,  SPEC7mcm(:,:)  &  ! nSPEC7mcm x 6
      &                     ,  SPEC9mcm(:,:)  &  ! nSPEC9mcm x 10
      &                     ,  S4H2O(:,:)     &  ! nS4H2o x 4
      &                     ,  T1H2O(:,:)     &  ! nT1H2o x 2
      &                     ,  Meskhidze(:,:) &  ! nMeskhidze x 7
      &                     ,  ASPEC1(:,:), ASPEC3(:,:) &! nASPEC1,3 x 2
      &                     ,  ASPEC2(:,:), ASPEC4(:,:) &! nAPSEC2,4 x 3
      &                     ,  DTEMP(:,:), DTEMP5(:,:)  &! nDTEMP,5 x 3
      &                     ,  DTEMP2(:,:), DTEMP3(:,:), DTEMP4(:,:) & ! nDTEMP2,3,4 x 4
      &                     ,  DCONST(:,:)    &    ! nDCONST x 2
      &                     ,  HENRY(:,:)     &    ! nHENRY x 2
      &                     ,  HOaqua(:)      &
      &                     ,  PHOTOkpp(:), PHOTO2kpp(:), PHOTO3kpp(:) &! nPHOTOi x 1
      &                     ,  HOM1(:,:)
    END TYPE ReacTypeIndex_TR

    TYPE(ReacTypeIndex_TR) :: iR


    TYPE ReacTypeParameter_CK
      ! different kinds of arrhenius parameter
      REAL(dp), ALLOCATABLE :: A(:),    b(:),    E(:)
      REAL(dp), ALLOCATABLE :: A0(:),   b0(:),   E0(:)
      REAL(dp), ALLOCATABLE :: AX(:),   bX(:),   EX(:)
      REAL(dp), ALLOCATABLE :: Ainf(:), binf(:), Einf(:)
      ! Troe parameter
      REAL(dp), ALLOCATABLE :: T1(:), T2(:), T3(:), T4(:)
    END TYPE ReacTypeParameter_CK

    TYPE(ReacTypeParameter_CK) :: RTpar

    INTEGER, ALLOCATABLE  :: AtomicMatrix(:,:)  ! dim = (nspc, natoms), where natoms = number of different elements in the system

    INTEGER, ALLOCATABLE  :: iFO_kat(:,:) ! reaction number where stoech coef == ONE
    INTEGER, ALLOCATABLE  :: iFO(:,:) ! reaction number where stoech coef == ONE
    INTEGER, ALLOCATABLE  :: iSO(:,:) ! reactions where species have second order 
    INTEGER, ALLOCATABLE  :: iHO(:,:) ! higher order reactions
    REAL(dp), ALLOCATABLE :: aHO(:) ! higher order reactions contains also fractions and noninteger values
    INTEGER :: nFirst_order=0, nSecond_order=0, nHigher_order=0, nFirst_orderKAT=0

    ! array for Troe reactions
    REAL(dp), ALLOCATABLE :: vlog10_Pr(:),      &
                           & vlog10_Fcent(:),   &
                           & vPr(:),            &
                           & vcTroe(:), vn1Troe(:)
    !
    ! stuff for the mechanism reduction
    CHARACTER(80), ALLOCATABLE :: Target_Names(:)
    INTEGER,       ALLOCATABLE :: Target_Index(:)


    ! index set of species kind
    TYPE SPC_Index_T
      CHARACTER(100), ALLOCATABLE :: Name(:,:)         ! name of species (iSpc,iFrac)
      INTEGER,  ALLOCATABLE :: Ind(:,:)   ! species index   (iSpc,iFrac)
      REAL(dp), ALLOCATABLE :: Conc(:,:)  ! in molec/cm3    (iSpc,iFrac)
    END TYPE SPC_Index_T

    TYPE (SPC_Index_T) :: iGas, iAqua, iKat


    INTEGER, SAVE        :: bGs(2),bAs(2),bSs(2),bPs(2)          ! phase boundaries
		INTEGER, ALLOCATABLE :: bGr(:),bHr(:),bAr(:),bSr(:),bPr(:)   ! phase boundaries

		INTEGER, ALLOCATABLE :: iGs(:),iAs(:),iSs(:),iPs(:)          ! indices phases
    INTEGER, ALLOCATABLE :: iGr(:),iHr(:),iAr(:),iSr(:),iPr(:)   ! indices phases



  TYPE AFrac_T
    CHARACTER(LenName), ALLOCATABLE :: Species(:)
    REAL(dp),           ALLOCATABLE :: MolMass(:)   ! [g/mol]
    INTEGER,            ALLOCATABLE :: Charge(:)    ! ladung (+,-,++,--,...)
    REAL(dp),           ALLOCATABLE :: SolubInd(:)  ! Löslichkeitsindex
    REAL(dp),           ALLOCATABLE :: Frac1(:)     ! [g/g]
  END TYPE AFrac_T

  TYPE(AFRAC_T) :: AFrac


  TYPE Modes_T
    REAL(dp), ALLOCATABLE :: Radius(:)    ! [m] radius particle
    REAL(dp), ALLOCATABLE :: wetRadius(:) ! [m] radius droplett
    REAL(dp), ALLOCATABLE :: Number(:)    ! [#/cm3]
    REAL(dp), ALLOCATABLE :: Density(:)   ! [kg/m3]
  END TYPE Modes_T

  TYPE(Modes_T) :: Mode

END MODULE mo_reac
