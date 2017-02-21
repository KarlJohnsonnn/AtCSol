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
    REAL(RealKind) :: dcon1,dcon2

!--------------------------------------------------------------
!---  Reaction Mechanism
!--------------------------------------------------------------
!
!---  reaction types
    TYPE def_para
      CHARACTER(RealKind) :: str_type
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
       REAL(RealKind) :: d_koef
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

       REAL(RealKind), POINTER  :: last_rate(:,:)
       REAL(RealKind), POINTER  :: back_rate(:,:)
       REAL(RealKind), POINTER  :: v(:,:)
       REAL(RealKind), POINTER  :: w(:,:)
       REAL(RealKind), POINTER :: dparam(:)  

       LOGICAL :: odd = .FALSE.

       REAL(RealKind) :: fac_exp
       REAL(RealKind) :: fac_A             

       TYPE (reactant), POINTER :: reactant(:)

       TYPE (reaction), POINTER :: next
       TYPE (reaction), POINTER :: next_all
    END TYPE reaction

!--    reaction pointers
    TYPE (reaction), POINTER :: all_first
    TYPE (reaction), POINTER :: gas_first
    TYPE (reaction), POINTER :: henry_first
    TYPE (reaction), POINTER :: dissoc_first
    TYPE (reaction), POINTER :: solid_first
    TYPE (reaction), POINTER :: aqua_first
    TYPE (reaction), POINTER :: gphoto_first
    TYPE (reaction), POINTER :: aphoto_first
    TYPE (reaction), POINTER :: gconst_first
    TYPE (reaction), POINTER :: aconst_first
    TYPE (reaction), POINTER :: buf_first

!--------------------------------------------------------------
!--   dimensions
    INTEGER :: nt, ntsolid, ntpart
    INTEGER :: ntGas          ! number of gaseous species
    INTEGER :: ntAqua         ! number of aqueeous species
    INTEGER :: ntKat          ! number of katalysator species
    INTEGER :: neq,nspc,nkat,nreak,nHenry
    INTEGER :: nDIM
    INTEGER :: nDIMcl
    INTEGER :: nDIMex

    INTEGER :: nreakgas,nreakhenry,nreakdissoc,nreakaqua, nreaksolid
    INTEGER :: nreakgphoto,nreakgconst,nreakgtemp,nreakgtroe,nreakgspec
    INTEGER :: nreakaphoto,nreakaconst,nreakatemp,nreakaspec
    INTEGER :: nreaksolidtemp,nreaksolidequi,nreaksolidspec
!    INTEGER :: nreakstemp,nreaksequi,nreaksspec
!--    define indices of special species
    INTEGER ::   hp_ind,         &   ! Index Hp
&               ohm_ind,         &   ! Index OHm
&              ah2o_ind,         &   ! Index aH2O
&              Temp_ind              ! Index Temperatur

!--------------------------------------------------------------
!--    passive species, indices of Henry species
    REAL(RealKind), ALLOCATABLE :: ykat(:,:)
    INTEGER, ALLOCATABLE :: ind_henry(:,:)
    CHARACTER(20), ALLOCATABLE :: ReactionTypes(:)

!--    mass transfer coefficient, kBig
    REAL(RealKind):: dkmt, kBig
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
    REAL(RealKind), ALLOCATABLE :: SumRO2(:)       ! Summation of all gasphase peroxyradicals in everey cell
!--    aqueous phase
    LOGICAL              :: hasRO2aq
    INTEGER, ALLOCATABLE :: RO2aq(:)          ! Species-Index of peroxyradical
    INTEGER              :: nRO2aq            ! Number of aqueos phase peroxy radicals in mechanism
    REAL(RealKind), ALLOCATABLE :: SumRO2aq(:)       ! Summation of all aqueopus phase peroxyradicals in everey cell
!
!--------------------------------------------------------------
!--    Aerosol species properties
    REAL(RealKind), ALLOCATABLE :: Charge(:)       ! charge of ions
    REAL(RealKind), ALLOCATABLE :: SolubInd(:)     ! solubility index
    REAL(RealKind), ALLOCATABLE :: MolMass(:)      ! molar mass of species
    REAL(RealKind), ALLOCATABLE :: SpcDens(:)      ! density of species
    REAL(RealKind), ALLOCATABLE :: OrgIndex(:)     ! carbon atoms
    CHARACTER(2), ALLOCATABLE :: CC(:)      ! compound class
    REAL(RealKind), ALLOCATABLE :: ActIndex(:)     ! Computation of Activity coefficients index 

!--    Input from impactor measurements
    INTEGER :: nImpac = 0                 & ! number of impactor stages
              ,nMode  = 0                 & ! number of input modes
              ,ntImp  = 0                   ! number of impactor species
    INTEGER, ALLOCATABLE :: ImpInd(:)       ! indices of impactor species
    REAL(RealKind), ALLOCATABLE :: BdImpac(:)      ! boundaries of stages
    REAL(RealKind), ALLOCATABLE :: FracImpac(:,:)  ! mass fractions per stage
    REAL(RealKind), ALLOCATABLE :: DryMod(:,:)     ! modal input masses per fractions

!--    Input from impactor measurements
    INTEGER :: ntAConc = 0                   ! number  of initialized concentrations
    INTEGER, ALLOCATABLE :: AConcInd(:)      ! indices of initialized concentrations
    REAL(RealKind), ALLOCATABLE :: AquaConc(:)      ! initial aqueous phase concentrations 
                                             ! after activation
!--------------------------------------------------------------
!--    input arrays
    !REAL(RealKind):: tcur
    !REAL(RealKind), ALLOCATABLE :: y_iconc(:)
    REAL(RealKind), ALLOCATABLE :: InitValAct(:)
    REAL(RealKind), ALLOCATABLE :: InitValKat(:)

    !REAL(RealKind), ALLOCATABLE :: y_emi(:)
    REAL(RealKind), ALLOCATABLE :: y_udepo(:)
    REAL(RealKind), ALLOCATABLE :: y_c1(:)
    REAL(RealKind), ALLOCATABLE :: y_c2(:)

!--------------------------------------------------------------
!--    deposition and emissions
    REAL(RealKind), ALLOCATABLE :: vd(:)
    REAL(RealKind), ALLOCATABLE :: y_e(:)

!--------------------------------------------------------------
    
    CHARACTER(18),ALLOCATABLE :: ThSpecies(:)
    CHARACTER(6),ALLOCATABLE  :: ThRefDataCode(:)
    CHARACTER(2),ALLOCATABLE  :: ThAtoms(:,:)           ! (:,4)
    INTEGER,ALLOCATABLE       :: ThnAtoms(:,:)          ! (:,4)
    REAL(RealKind),ALLOCATABLE :: ThTempRange(:,:)      ! (:,2)
    REAL(RealKind),ALLOCATABLE :: ThMolMass(:)
    CHARACTER(1),ALLOCATABLE  :: ThPhase(:)
    REAL(RealKind),ALLOCATABLE :: lowA(:),lowB(:),lowC(:),lowD(:),lowE(:),lowF(:),lowG(:)
    REAL(RealKind),ALLOCATABLE :: highA(:),highB(:),highC(:),highD(:),highE(:),highF(:),highG(:)
    REAL(RealKind),ALLOCATABLE :: H0_29815R(:)
    REAL(RealKind),ALLOCATABLE :: SwitchTemp(:)
    !
    REAL(RealKind) :: y_total
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
    REAL(RealKind), ALLOCATABLE :: GFE(:), DGFEdT(:)
    REAL(RealKind), ALLOCATABLE :: DelGFE(:), DDelGFEdT(:)
    !
END MODULE mo_reac
