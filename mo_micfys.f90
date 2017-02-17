!==============================================================
!===  MODULE for Microphysics
!==============================================================
!
   MODULE mo_micfys

!------------------------------------------------------------------
!---  dimensions
      INTEGER :: JMAX=66, NMAX=2
!===      JMAX=132 ==> NMAX=3, JMAX=264 ==> NMAX=5
!
      INTEGER :: JHMAX,NHMAX,JHFAK,JHFAK2
!===      JHMAX=JMAX, NHMAX=NMAX, JHFAK=JHMAX/JMAX, JHFAK2=JHFAK**2

      INTEGER :: dimS = 500  & ! Dimension of Kernel Table (Solid)
&               ,dimW = 500    ! Dimension of Kernel Table (Water)

      INTEGER :: nComp = 7     ! Number of Compounds
!
!------------------------------------------------------------------
!--- Microphysics control variables
      INTEGER :: iruf, ifirst, icall, irhum, iclose, iwahl, istop  &
     &   ,iextern, icoarse, icond, ikoll, ibrea, idepo, imoving    &
     &   ,iideal, inacl, ifeedback, isolutens   
!
!--- Counters and flags
      INTEGER :: ISTEP, ITIME, ITIM1, ITIM2
      INTEGER :: IterMax, i_akt

      INTEGER :: J12, J2l, J2u, J5l, J5u
      INTEGER :: INCOUT,INCOUT1,INCOUT2,IOUT2

!--- INTEGER arrays
      INTEGER, ALLOCATABLE :: IAKT(:)
      INTEGER, ALLOCATABLE :: KZIEL(:,:)

!------------------------------------------------------------------
!---  REAL(8) variables
      REAL(8) :: TABS,DELT,QV,DQV,PTOT,DP,RHOTOT,DRHO,RHOA  &
&        ,TIM,SATT,SATTini,ENTR_FAC,NWAKT, ES,SATT_Start, Rho0

!---  Sums of spectral variables
      REAL(8) :: QSSUM,QWSUM,NWSUM,SumNW, SumQW, SumQS

!---  Prescribed Supersaturation (only if mCase = 5)
      REAL(8) :: SSAT = 0.E0   

!---  Microphysics constants 
      REAL(8) :: CP,GRAV,KTH,LS,LV,mol_w,PI,RG_DRY,RW,RHOW,SIGMA0  &
     &   ,mol_AS,mol_AI,rho_AS,rho_AI,vantS,E0,C1,C2               &
     &   ,CAPPA,FAC34,FACT,fact_AS,fact_AI,PC1,PC2,QC1,QC2       &
     &   ,SMALL1,RHOTOT0,QU1D3,QU4D3,QU5D3                       

!---  Microphysics help variables 
      REAL(8) ::    &
     &    DELTAT,TIMMAX,TIM1,TIM2,size_fact,accu_fact,step_fact  &
     &   ,R0,RAmin,RAmax,PR,TC,mu,bb1,b2,bb,ALPHA_C,ALPHA_T      &
     &   ,DEL_V,DEL_T,WEI_Q,WEI_N,SC_cond,SC_ggw,SATTabs         &
     &   ,TIMEall,dpdh,p_akt,v_akt,l_akt,h_akt
!
!------------------------------------------------------------------
!---  Allocatable arrays
!------------------------------------------------------------------
!
!--- Microphysical variables
      REAL(8), ALLOCATABLE ::   &
     &         QW(:),NW(:),QS(:) ,DQW(:),DNW(:),DQS(:)        &
     &        ,MQUER(:),CQUER(:),AQUER(:), epsilon(:),RHOT(:) &
     &        ,QWteil(:,:),QSteil1(:,:),QSteil2(:,:)    
!
!--- Microphysical variables: Condensation
      REAL(8), ALLOCATABLE :: r_dry(:),r_wet(:),m_dry(:),m_sol(:)  &
     &   ,r_sol(:),SATTeq(:),DMFAK(:),AAA(:),BBB(:)                &
     &   ,S_krit(:),r_krit(:),m_krit(:)

!--- Mean particle density
      REAL(8), ALLOCATABLE :: DensPart(:)

!--- Feedback variables
      REAL(8), ALLOCATABLE :: Raou(:), RaouBBB(:), RaoultTerm(:)
      REAL(8), ALLOCATABLE :: SurfTens(:)
 
!--- Spectrum and Cernel
      REAL(8), ALLOCATABLE ::   &
     &    DIFF21(:),MGRENZ(:),MMITTE(:),RGRENZ(:),RMITTE(:)  &
     &   ,KOLK(:,:),VMITTE(:),KOLK1D(:)
!
!---  Konstanten
      REAL(8), ALLOCATABLE :: ALPHA(:),CFACT(:),QQQ(:)
      REAL(8), ALLOCATABLE :: QS_CTM(:), QW_CTM(:)
!
!---  Coumpounds characterizing constanms
      REAL(8), ALLOCATABLE :: mol_g(:),phi_g(:),a_phi(:),b_phi(:)

!---  Character
      CHARACTER(13) :: ch_path

!------------------------------------------------------------------
!-- Trajectory values
      INTEGER ::  nTraj, ModQV
      REAL(8), ALLOCATABLE :: location(:), altitude(:), pressure(:)
      REAL(8), ALLOCATABLE :: slope(:)
      REAL(8), ALLOCATABLE :: ttime(:), twater(:), temperature(:)
      REAL(8), ALLOCATABLE :: GradQV(:)

      REAL(8) :: hwind  =   2.e0                  !  horizontal wind [m/s]

!------------------------------------------------------------------
!-- Initial aerosol distribution
      INTEGER ::  nModes
      REAL(8), ALLOCATABLE :: ApInit(:,:)
!
!-- Initial aerosol spectrum
      INTEGER ::  nSpek
      REAL(8), ALLOCATABLE :: SpekInit(:,:)
!
!------------------------------------------------------------------

!---  LOGICAL
      LOGICAL L_Feed
!
   END MODULE mo_micfys
!
!==============================================================
!===  MODULE for FEBUKO scenarios
!==============================================================
!
   MODULE mo_febuko
!
      INTEGER :: ifebuko,ihour,iquart,idynamik,iinvers,istation  &
     &          ,iwind_gl,iwind_sm,iwind_gb
!
      INTEGER :: ProfMax = 23      & ! Dimension of orography profile
     &          ,mMax    = 5
!
      REAL(8) ::   &                                          
     &    TABS0gl,PTOT0gl,RHOTOT0gl,SATT0gl,WIND0gl,QVgl   &
     &   ,TABS0sm,PTOT0sm,RHOTOT0sm,SATT0sm,WIND0sm,QVsm   &
     &   ,TABS0gb,PTOT0gb,RHOTOT0gb,SATT0gb,WIND0gb,QVgb   &
     &   ,LWCAsm,LWCBsm,LWCCsm,LWCFsm,PSADsm               &
     &   ,REFFsm,REFFAsm,REFFFsm,DROPsm,NAPBsm
      REAL(8) :: APm_gl,APm_sm,APm_gb
!
      REAL(8), ALLOCATABLE ::   &                                          
     &    APn_gl(:),APs_gl(:),APd_gl(:)     &
     &   ,APn_sm(:),APs_sm(:),APd_sm(:)     &
     &   ,APn_gb(:),APs_gb(:),APd_gb(:)     &
     &   ,l_pro(:),h_pro(:),steep(:)
!
   END MODULE mo_febuko
