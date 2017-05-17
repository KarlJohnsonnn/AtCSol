!==============================================================
!===  MODULE for the Description of Chemistry, Deposition
!===  and Emissions
!==============================================================
!
  MODULE mo_reac
    USE Kind_Mod

!---  Files and Units
    CHARACTER(30) :: file_name

    CHARACTER(9) :: measure_gas_ph(2)=(/"molec/cm3", "mol/m3   "/)
    CHARACTER(9) :: measure_aqua_ph(1)=(/"mol/l"/)
    CHARACTER(12):: units(2)

    INTEGER, PARAMETER     :: a_unit  = 2   &
                             ,a_gas_u = 2
    REAL(dp) :: dcon1,dcon2

!--------------------------------------------------------------
!---  Reaction Mechanism
!--------------------------------------------------------------
!
!---  reaction types
    TYPE def_para
      CHARACTER(8) :: str_type
      INTEGER :: anzp
      INTEGER :: type
    END TYPE def_para
  
    INTEGER, PARAMETER ::        &
&              nr_reac  = 37     &  ! number of reaction types
&             ,nr_cons  =  4     &  ! constant rate
&             ,nr_sgasA = 16     &  ! first special gas reaction 
&             ,nr_sgasE = 35     &  ! last  special gas reaction 
&             ,nr_equi  = 36     &  ! equilibrium (for solubility)
&             ,nr_dtemp3  = 13        ! equilibrium (for solubility)


    TYPE(def_para), DIMENSION(nr_reac) :: var_par = &
                   (/def_para("PHOTABC", 3,  1), &
                     def_para("PHOTMCM", 3,  2), &
                     def_para("PHOTAB",  2,  3), &
                     def_para("CONST",   1,  4), &
                     def_para("TEMP1",   2,  5), &
                     def_para("TEMP2",   2,  6), &
                     def_para("TEMP3",   2,  7), &
                     def_para("TEMP4",   2,  8), &
                     def_para("TEMPX",   3,  9), &
                     def_para("ASPEC1",  2, 10), &
                     def_para("ASPEC2",  3, 11), &
                     def_para("ASPEC3",  2, 12), &
                     def_para("DCONST",  2, 13), &
                     def_para("DTEMP",   3, 14), &
                     def_para("DTEMP2",  4, 15), &
                     def_para("DTEMP3",  4, 16), &
                     def_para("TROE",    4, 17), &
                     def_para("TROEQ",   6, 18), &
                     def_para("TROEF",   5, 19), &
                     def_para("TROEQF",  7, 20), &
                     def_para("TROEXP",  5, 21), &
                     def_para("TROEMCM",10, 22), &
                     def_para("SPEC1",   2, 23), &
                     def_para("SPEC2",   2, 24), &
                     def_para("SPEC3",   6, 25), &
                     def_para("SPEC4",   4, 26), &
                     def_para("SPEC1MCM",3, 27), &
                     def_para("SPEC2MCM",3, 28), &
                     def_para("SPEC4MCM",4, 29), &
                     def_para("SPEC3MCM",2, 30), &
                     def_para("SPEC5MCM",4, 31), &
                     def_para("SPEC6MCM",4, 32), &
                     def_para("SPEC7MCM",6, 33), &
                     def_para("SPEC8MCM",4, 34), &
                     def_para("T1H2O",   2, 35), &
                     def_para("S4H2O",   4, 36), &
                     def_para("EQUI",    1, 37) /)
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

!--------------------------------------------------------------
!--   dimensions
    INTEGER :: nt=0, ntsolid=0, ntpart=0
    INTEGER :: ntGas=0                 ! number of gaseous species
    INTEGER :: ntAqua=0               ! number of aqueeous species
    INTEGER :: ntKat=0                ! number of katalysator species
    INTEGER :: ntKatGas=0,ntKatAqua=0   ! number of katalysator species gaseous ,aqueous
    INTEGER :: neq=0,nspc=0,nreak=0
    INTEGER :: nDIM=0
    INTEGER :: nDIMcl=0
    INTEGER :: nDIMex=0
    REAL(dp) :: rNspc

    INTEGER :: nreakgas=0,nreakhenry=0,nreakdissoc=0,nreakaqua=0, nreaksolid=0
    INTEGER :: nreakgphoto=0,nreakgconst=0,nreakgtemp=0,nreakgtroe=0,nreakgspec=0
    INTEGER :: nreakaphoto=0,nreakaconst=0,nreakatemp=0,nreakaspec=0,nreakglind=0
    INTEGER :: nreaksolidtemp=0,nreaksolidequi=0,nreaksolidspec=0,nreakSimpTB=0
    INTEGER :: nreakpress=0
!    INTEGER :: nreakstemp,nreaksequi,nreaksspec
!--    define indices of special species
    INTEGER ::   Hp_ind,         &   ! Index Hp
&               OHm_ind,         &   ! Index OHm
&              aH2O_ind,         &   ! Index aH2O
&              Temp_ind              ! Index Temperatur

!--------------------------------------------------------------
!--    passive species, indices of Henry species
    REAL(dp), ALLOCATABLE :: ykat(:,:)
    INTEGER, ALLOCATABLE :: ind_henry(:,:)

!--    mass transfer coefficient, kBig
    REAL(dp):: dkmt, kBig
!
!--------------------------------------------------------------
!--    species names 
    CHARACTER(60), ALLOCATABLE :: y_name(:)
    CHARACTER(60), ALLOCATABLE :: GasName(:)    & ! gas phase species
&                                ,AquaName(:)   & ! aqueous phase species
&                                ,SolidName(:)  & ! solid species
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
    CHARACTER(20), ALLOCATABLE  :: ReacName_Gas(:)  &
&                                 ,ReacName_Hen(:)  &
&                                 ,ReacName_Diss(:) &
&                                 ,ReacName_Aqua(:) &
&                                 ,ReacName_Sol(:)
    INTEGER,  ALLOCATABLE       :: Gas_anz(:)  &
&                                 ,Aqua_anz(:) &
&                                 ,Hen_anz(:)  &
&                                 ,Diss_anz(:) &
&                                 ,Sol_anz(:)
    INTEGER, PARAMETER          :: ndiag_gas = 11  &
&                                 ,ndiag_aq = 16  
    CHARACTER(60), ALLOCATABLE  :: Diag_Name_Netcdf(:) &
&                                , Diag_LongName(:)    &
&                                , DiagERR_Name_Netcdf(:)

    CHARACTER(60)               :: Diag_Name_Gas(ndiag_gas)  &
&                                 ,Diag_Name_Aqua(ndiag_aq)  

    INTEGER, ALLOCATABLE  :: iDiag_Schwefel(:)  
!--------------------------------------------------------------
!--    Peroxyradicals
    LOGICAL              :: hasRO2
    INTEGER, ALLOCATABLE :: RO2(:)          ! Species-Index of peroxyradical
    INTEGER              :: nRO2            ! Number of Gasphase-Peroxyradicals in mechanism
    REAL(dp), ALLOCATABLE :: SumRO2(:)       ! Summation of all gasphase peroxyradicals in everey cell
!--    aqueous phase
    LOGICAL              :: hasRO2aq
    INTEGER, ALLOCATABLE :: RO2aq(:)          ! Species-Index of peroxyradical
    INTEGER              :: nRO2aq            ! Number of aqueos phase peroxy radicals in mechanism
    REAL(dp), ALLOCATABLE :: SumRO2aq(:)       ! Summation of all aqueopus phase peroxyradicals in everey cell
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

!--    Input from impactor measurements
    INTEGER :: ntAConc = 0                   ! number  of initialized concentrations
    INTEGER, ALLOCATABLE :: AConcInd(:)      ! indices of initialized concentrations
    REAL(dp), ALLOCATABLE :: AquaConc(:)      ! initial aqueous phase concentrations 
                                             ! after activation
!--------------------------------------------------------------
!--    input arrays
    !REAL(dp):: tcur
    !REAL(dp), ALLOCATABLE :: y_iconc(:)
    REAL(dp), ALLOCATABLE :: InitValAct(:)
    REAL(dp), ALLOCATABLE :: InitValKat(:)

    !REAL(dp), ALLOCATABLE :: y_emi(:)
    REAL(dp), ALLOCATABLE :: henry_diff(:)
    REAL(dp), ALLOCATABLE :: henry_accom(:)

!--------------------------------------------------------------
!--    deposition and emissions
    REAL(dp), ALLOCATABLE :: vd(:)
    REAL(dp), ALLOCATABLE :: y_e(:)

!--------------------------------------------------------------
    
    CHARACTER(18),ALLOCATABLE :: ThSpecies(:)
    CHARACTER(6),ALLOCATABLE  :: ThRefDataCode(:)
    CHARACTER(2),ALLOCATABLE  :: ThAtoms(:,:)           ! (:,4)
    INTEGER,ALLOCATABLE       :: ThnAtoms(:,:)          ! (:,4)
    REAL(dp),ALLOCATABLE :: ThTempRange(:,:)      ! (:,2)
    REAL(dp),ALLOCATABLE :: ThMolMass(:)
    CHARACTER(1),ALLOCATABLE  :: ThPhase(:)
    REAL(dp),ALLOCATABLE :: lowA(:),lowB(:),lowC(:),lowD(:),lowE(:),lowF(:),lowG(:)
    REAL(dp),ALLOCATABLE :: highA(:),highB(:),highC(:),highD(:),highE(:),highF(:),highG(:)
    REAL(dp),ALLOCATABLE :: H0_29815R(:)
    REAL(dp),ALLOCATABLE :: SwitchTemp(:)
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
    REAL(dp), ALLOCATABLE :: SCperm(:)

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
      INTEGER, ALLOCATABLE :: iPHOTabc(:),  iPHOTab(:),   iPHOTmcm(:), iCONST(:)
      INTEGER, ALLOCATABLE :: iTEMP1(:),    iTEMP2(:),    iTEMP3(:),   iTEMP4(:)
      INTEGER, ALLOCATABLE :: iTROE(:),     iTROEf(:),    iTROEq(:),   iTROEqf(:), iTROExp(:)
      INTEGER, ALLOCATABLE :: iTROEmcm(:),  iSPEC1(:),    iSPEC2(:),   iSPEC3(:),  iSPEC4(:)
      INTEGER, ALLOCATABLE :: iSPEC1mcm(:), iSPEC2mcm(:), iSPEC3mcm(:), iSPEC4mcm(:)
      INTEGER, ALLOCATABLE :: iSPEC5mcm(:), iSPEC6mcm(:), iSPEC7mcm(:), iSPEC8mcm(:)
      INTEGER, ALLOCATABLE :: iS4H2O(:),    iT1H2O(:),    iASPEC1(:),   iASPEC2(:)
      INTEGER, ALLOCATABLE :: iASPEC3(:),   iASPEC4(:),   iDCONST(:,:),   iDTEMP(:,:)
      INTEGER, ALLOCATABLE :: iDTEMP2(:,:), iDTEMP3(:,:), iDTEMP4(:,:), iDTEMP5(:,:), iMeskhidze(:,:)
      INTEGER, ALLOCATABLE :: iHENRY(:,:)
      ! nFACTOR x 2, with [H2,O2N2,M,O2,N2,H2O,RO2,O2O2,aH2O,RO2aq,+M/(+M)]
      INTEGER, ALLOCATABLE :: iFAC_H2(:),  iFAC_O2N2(:), iFAC_M(:),    iFAC_O2(:),   iFAC_N2(:)
      INTEGER, ALLOCATABLE :: iFAC_H2O(:), iFAC_RO2(:),  iFAC_O2O2(:), iFAC_aH2O(:), iFAC_RO2aq(:)
      INTEGER, ALLOCATABLE :: iHOaqua(:)
      INTEGER, ALLOCATABLE :: iPHOTOkpp(:), iPHOTO2kpp(:), iPHOTO3kpp(:)
    END TYPE ReacTypeIndex_TR

    INTEGER :: nPHOTabc=0, nPHOTab=0, nPHOTmcm=0, nCONST=0, nTEMP1=0, nTEMP2=0, nTEMP3=0, nTEMP4=0
    INTEGER :: nTROE=0, nTROEf=0, nTROEq=0, nTROEqf=0, nTROExp=0, nTROEmcm=0, nSPEC1=0, nSPEC2=0
    INTEGER :: nSPEC3=0, nSPEC4=0, nSPEC1mcm=0, nSPEC2mcm=0, nSPEC3mcm=0, nSPEC4mcm=0, nSPEC5mcm=0
    INTEGER :: nSPEC6mcm=0, nSPEC7mcm=0, nSPEC8mcm=0, nS4H2O=0, nT1H2O=0, nASPEC1=0, nASPEC2=0
    INTEGER :: nASPEC3=0, nASPEC4=0
    INTEGER :: nDTEMP=0, nDTEMP2=0, nDTEMP3=0, nDTEMP4=0, nDTEMP5=0, nDCONST=0
    INTEGER :: nMeskhidze=0, nHENRY=0, nHENRYga=0, nHENRYag=0
    INTEGER :: nFACTOR=0, nFAC_H2=0, nFAC_O2N2=0, nFAC_M=0, nFAC_O2=0, nFAC_N2=0
    INTEGER :: nFAC_H2O=0, nFAC_RO2=0, nFAC_O2O2=0, nFAC_aH2O=0, nFAC_RO2aq=0
    INTEGER :: nHOaqua=0
    INTEGER :: nPHOTOkpp=0, nPHOTO2kpp=0, nPHOTO3kpp=0

    TYPE(ReacTypeIndex_TR) :: RTind2


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

    TYPE ReacTypeParameter_TR
      ! different kinds of arrhenius parameter
      REAL(dp), ALLOCATABLE :: PHOTabc(:,:)     ! nPHOTabc x 3
      REAL(dp), ALLOCATABLE :: PHOTab(:,:)      ! nPHOTab x 2
      REAL(dp), ALLOCATABLE :: PHOTmcm(:,:)     ! nPHOTmcm x 3
      REAL(dp), ALLOCATABLE :: CONST(:)         ! nCONST x 1
      REAL(dp), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:), TEMP3(:,:), TEMP4(:,:) ! nTEMPi x 2
      REAL(dp), ALLOCATABLE :: TROE(:,:)        ! nTROE x 4
      REAL(dp), ALLOCATABLE :: TROEf(:,:)       ! nTROEf x 5
      REAL(dp), ALLOCATABLE :: TROEq(:,:)       ! nTROEq x 6
      REAL(dp), ALLOCATABLE :: TROEqf(:,:)      ! nTROEqf x 7
      REAL(dp), ALLOCATABLE :: TROExp(:,:)      ! nTROExp x 5
      REAL(dp), ALLOCATABLE :: TROEmcm(:,:)     ! nTROExp x 10
      REAL(dp), ALLOCATABLE :: SPEC1(:,:), SPEC2(:,:)  ! nSPEC1,2 x 2
      REAL(dp), ALLOCATABLE :: SPEC3(:,:)       ! nSPEC3 x 6
      REAL(dp), ALLOCATABLE :: SPEC4(:,:)       ! nSPEC4 x 4
      REAL(dp), ALLOCATABLE :: SPEC1mcm(:,:), SPEC2mcm(:,:)  ! nSPEC1mcm,2 x 3
      REAL(dp), ALLOCATABLE :: SPEC3mcm(:,:)    ! nSPEC3mcm x 2
      REAL(dp), ALLOCATABLE :: SPEC4mcm(:,:), SPEC5mcm(:,:), SPEC6mcm(:,:), SPEC8mcm(:,:)  ! nSPEC4,5,6mcm x 4
      REAL(dp), ALLOCATABLE :: SPEC7mcm(:,:)    ! nSPEC7mcm x 6
      REAL(dp), ALLOCATABLE :: S4H2O(:,:)       ! nS4H2o x 4
      REAL(dp), ALLOCATABLE :: T1H2O(:,:)       ! nT1H2o x 2
      REAL(dp), ALLOCATABLE :: Meskhidze(:,:)   ! nMeskhidze x 7
      REAL(dp), ALLOCATABLE :: ASPEC1(:,:), ASPEC3(:,:) ! nASPEC1,3 x 2
      REAL(dp), ALLOCATABLE :: ASPEC2(:,:), ASPEC4(:,:) ! nAPSEC2,4 x 3
      REAL(dp), ALLOCATABLE :: DTEMP(:,:), DTEMP5(:,:)  ! nDTEMP,5 x 3
      REAL(dp), ALLOCATABLE :: DTEMP2(:,:), DTEMP3(:,:), DTEMP4(:,:)  ! nDTEMP2,3,4 x 4
      REAL(dp), ALLOCATABLE :: DCONST(:,:)        ! nDCONST x 2
      REAL(dp), ALLOCATABLE :: HENRY(:,:)         ! nHENRY x 2
      REAL(dp), ALLOCATABLE :: HOaqua(:)
      REAL(dp), ALLOCATABLE :: PHOTOkpp(:), PHOTO2kpp(:), PHOTO3kpp(:) ! nPHOTOi x 1
    END TYPE ReacTypeParameter_TR

    TYPE(ReacTypeParameter_TR) :: RTpar2


    INTEGER, ALLOCATABLE  :: first_orderKAT(:,:) ! reaction number where stoech coef == ONE
    INTEGER, ALLOCATABLE  :: first_order(:,:) ! reaction number where stoech coef == ONE
    INTEGER, ALLOCATABLE  :: second_order(:,:) ! reactions where species have second order 
    INTEGER, ALLOCATABLE  :: higher_order(:,:) ! higher order reactions
    REAL(dp), ALLOCATABLE :: ahigher_order(:) ! higher order reactions contains also fractions and noninteger values
    INTEGER :: nFirst_order=0, nSecond_order=0, nHigher_order=0, nFirst_orderKAT=0

    ! array for Troe reactions
    REAL(dp), ALLOCATABLE :: vlog10_Pr(:),      &
                           & vlog10_Fcent(:),   &
                           & vPr(:),            &
                           & vcTroe(:), vn1Troe(:)
    !
END MODULE mo_reac
