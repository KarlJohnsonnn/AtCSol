  !=========================================================================!
  !                                                                         !         
  !                                                                         !
  !               Calculating the Rates of the chemical system              ! 
  !                                                                         ! 
  !                                                                         !
  !=========================================================================!
  !
  MODULE Rates_Mod
    !
    USE Kind_Mod
    USE mo_reac
    USE mo_control
    USE mo_MPI
    USE Sparse_Mod
    USE Chemsys_Mod
    USE Meteo_Mod
    USE mo_ckinput, ONLY: lowA,lowB,lowC,lowD,lowE,lowF,lowG       &
    &                  , highA,highB,highC,highD,highE,highF,highG
    !
    IMPLICIT NONE
    !
    ! 
    REAL(dp) :: LAT=45.0_dp
    REAL(dp) :: LONG=0.0_dp
    REAL(dp) :: fac_exp=1.0_dp, fac_A=1.0_dp
    INTEGER  :: IDAT=010619
    !
    ! some factors for calculating Troe press dep. reactions
    REAL(dp) :: rFacEq      ! factor nessesary for equilibrium constant
    REAL(dp) :: cTroe, n1Troe
    REAL(dp), PARAMETER :: dTroe =  0.14_dp

    REAL(dp) :: log10_Pr, log10_Fcent, Pr
    INTEGER :: globi
    !
    CONTAINS
    !
    !
    !======================================================================!
    !      Calculate the Rates for current concentraions and time
    !======================================================================!
    SUBROUTINE ReactionRatesAndDerivative_ChemKin(Time,Y_in,Rate,DRatedT)
      !--------------------------------------------------------------------!
      ! Input: 
      !   - Time
      !   - y_conc
!     REAL(dp), INTENT(IN) :: Time
      REAL(dp) :: Time
      REAL(dp), INTENT(IN) :: Y_in(nDIM)
      !--------------------------------------------------------------------!
      ! Output:
      !   - Rate vector
      REAL(dp) :: Rate(neq)
      REAL(dp) :: DRatedT(neq)
      !--------------------------------------------------------------------!
      ! Temporary variables:
      REAL(dp) :: Conc(nspc)
      REAL(dp) :: T(10)
      REAL(dp) :: k,       vK(neq),   vKaqua(nAreak,ntFrac)
      REAL(dp) :: Prod,    vProd(neq) 
      REAL(dp) :: DkdT,    dvK_dT(neq)   
      REAL(dp) :: Meff,    vMeff(neq)
      INTEGER :: iReac, ii, j
      ! 
      REAL(dp) :: tmpK(neq), dtmpK(neq)
      
      ! temp arrays for chemkin input (temp depended)
      REAL(dp) :: vMeffX(RTind%nTBodyExtra)
      REAL(dp) :: vFTL(neq),      DF_PDdT(neq)
      REAL(dp) :: Dk0dT(neq),     DkinfdT(neq),      vdFTL_dT(neq)
      REAL(dp) :: k0(RTind%nLow), kinf(RTind%nHigh), k0M(RTind%nLow)
      REAL(dp) :: vrkinfpk0M(RTind%nLow)
      REAL(dp) :: vrKeq(RTind%nEqui), DeRdT(RTind%nEqui)
      REAL(dp) :: rRcT

      !==================================================================!
      !===============calc rates for ReactionSystem======================!
      !==================================================================!
     
      TimeRateA = MPI_WTIME()
      
      Rate    = ZERO
      DRatedT = ZERO

      !--- Concentration
      Conc = Y_in(1:nspc)
      
      T = UpdateTempArray ( Y_in(nDIM) )
      CALL GibbsFreeEnergie( GFE , T )
      CALL CalcDeltaGibbs  ( DelGFE )
      rFacEq = mega * R * T(1) * rPatm  ! in [cm3/mol]
      !
      CALL CalcDiffGibbsFreeEnergie( DGFEdT , T )
      CALL CalcDiffDeltaGibbs( DDelGFEdT )

      !*************************************************************
      ! --- compute rate of all reactions ---
      !*************************************************************
      !
      ! initializing vector components
      vMeff = ONE; DF_PDdT = Zero; 

      ! calculate effective molecularity of the reactions
      IF ( RTind%nTBody>0 ) THEN
        vMeff( RTind%iTBody ) = SUM( Conc )

        IF ( RTind%nTBodyExtra>0) THEN
          CALL DAX_sparse( vMeffX, TB_sparse, Conc )
          vMeff(RTind%iTBodyExtra) = vMeff(RTind%iTBodyExtra) - vMeffX
        END IF

      END IF

      ! ====== Compute the rate constant for specific reaction type
      !
      rRcT = rRcal*T(6)

      ! normal Arrhenius
      IF ( RTind%nArr>0 )  THEN
        tmpK(RTind%iArr)  = RTpar%A*EXP(RTpar%b*T(8) - RTpar%E*rRcT)
        dtmpK(RTind%iArr) = (RTpar%b + RTpar%E*rRcT)*T(6)
      END IF
     
      ! Backward reaction with explicitly given Arrhenius parameters
      IF ( RTind%nXrev>0 ) THEN
        tmpK(RTind%iXrev)  = RTpar%AX*EXP(RTpar%bX*T(8) - RTpar%EX*rRcT)
        dtmpK(RTind%iXrev) = (RTpar%bX + RTpar%EX*rRcT)*T(6)
      END IF

      ! reduced pressure value for falloff reactions
      vFTL = ONE;        vdFTL_dT = ONE
      Dk0dT = ZERO;      DkinfdT  = ZERO

      IF ( RTind%nLow>0 ) THEN
        ! Low pressure Arrhenius
        k0    = RTpar%A0*EXP(RTpar%b0*T(8) - RTpar%E0*rRcT)
        Dk0dT(RTind%iLow)    = RTpar%b0 + RTpar%E0*rRcT
        ! High pressure Arrhenius
        kinf    = RTpar%Ainf*EXP(RTpar%binf*T(8) - RTpar%Einf*rRcT)
        DkinfdT(RTind%iHigh) = RTpar%binf + RTpar%Einf*rRcT

        vPr = ONE

        k0M = k0 * vMeff(RTind%iLow)
        vrkinfpk0M = ONE / (kinf + k0M)

        vdFTL_dT(RTind%iLow) = vMeff(RTind%iTroe)*k0*kinf*(vrkinfpk0M*vrkinfpk0M) &
        &                    * (Dk0dT(RTind%iLow) - DkinfdT(RTind%iHigh))*T(6)

        vPr(RTind%iLow)    = k0M / (kinf+k0M)

        tmpK(RTind%iLow)   = kinf
        dtmpK(RTind%iHigh) = DkinfdT(RTind%iHigh)
      END IF
      
      ! reverse reaction, calculating the equi constant 
      IF (RTind%nEqui>0) THEN
        vrKeq = EXP(DelGFE(RTind%iEqui)) * rFacEq**(-sumBAT(RTind%iEqui))
        DeRdT = sumBAT(RTind%iEqui)*T(6) + DDelGFEdT(RTind%iEqui) 

        tmpK(RTind%iEqui)  = tmpK(RTind%iEqui) * vrKeq
        dtmpK(RTind%iEqui) = dtmpK(RTind%iEqui) + DeRdT 
      END IF

      ! rate constants according to Lindemann's form
      IF ( RTind%nLind>0 ) THEN
        vFTL(RTind%iLind)    = vPr(RTind%iLind)
        DF_PDdT(RTind%iLind) = vdFTL_dT(RTind%iLind)
        vMeff(RTind%iLind)   = ONE
      END IF

      ! rate constants according to Troe's form
      IF ( RTind%nTroe>0 ) THEN
        vFTL(RTind%iTroe)    = TroeFactorVec(T)
        DF_PDdT(RTind%iTroe) = DTroeFactorVec(Dk0dT,DkinfdT,vdFTL_dT,T)
        vMeff(RTind%iTroe)   = ONE
      END IF

      vK      = vFTL * tmpK
      dvK_dT  = (DF_PDdT + dtmpK)*T(6)
      DRatedT = dvK_dT


      ! ==== Law of mass action productories
      vProd = MassActionProducts( Conc )
      
      Rate  = vMeff * vk * vProd
      
      !WRITE(*,*) ''
      !DO j=1,neq
      !  WRITE(*,*) 'j,Time,Meff,k,prd,Rate =', &
      !  &           j,Time,vMeff(j),vK(j),vProd(j),Rate(j)
      !END DO

      TimeRates = TimeRates + MPI_WTIME() - TimeRateA

      !stop 'ratesmod'
    END SUBROUTINE ReactionRatesAndDerivative_ChemKin

    !======================================================================!
    !      Calculate the Rates for current concentraions and time
    !======================================================================!
    SUBROUTINE ReactionRates_Tropos(Time,Y_in,Rate)
      !--------------------------------------------------------------------!
      ! Input: 
      !   - Time
      !   - y_conc
!     REAL(dp), INTENT(IN) :: Time
      REAL(dp) :: Time
      REAL(dp), INTENT(IN) :: Y_in(nDIM)
      !--------------------------------------------------------------------!
      ! Output:
      !   - Rate vector
      REAL(dp) :: Rate(neq)
      !--------------------------------------------------------------------!
      ! Temporary variables:
      REAL(dp) :: Conc(nspc)
      REAL(dp) :: chi(3), LWC
      REAL(dp) :: T(10)
      REAL(dp) :: Meff(neq) , k(neq) , Prod(neq)
      REAL(dp) :: AquaFac(nHOAqua)
      
      !REAL(dp) :: tHenry(nHENRY,2)
      REAL(dp) :: tHenry(nHENRY,2,ntFrac)
      !==================================================================!
      !===============calc rates for ReactionSystem======================!
      !==================================================================!
     
      TimeRateA = MPI_WTIME()
      
      Rate = ZERO

      ! --- Compute zenith for photo reactions
      IF (nreakgphoto>0) THEN
        ! tropos syntax calculation of zenith
        chi(1) = Zenith(Time)
        chi(2) = ONE/COS(Chi(1))
        ! kkp syntax calculation of sun for photo reactions
        chi(3) = UpdateSun(Time)
      END IF
     
      
      Conc = Y_in
      T = UpdateTempArray ( 280.0_dp )
      
      !*************************************************************
      ! --- compute rate of all reactions (gas,henry,aqua,diss) ---
      !*************************************************************

      ! ====== Computing effective molecularity 
      Meff = ONE
      IF ( nFACTOR > 0 ) Meff = EffectiveMolecularity( Conc )
      
      ! ====== Compute the rate constant for specific reaction type
      CALL ComputeRateConstant( k, T, Time, chi, mAir, Conc, Meff )

      ! ===== correct unit of concentrations for higher order aqueous reactions
      IF ( ntAqua > 0 ) THEN
        LWC = pseudoLWC(Time)
        InitValKat(aH2O_ind) = aH2O*LWC
        AquaFac = ONE / (LWC*mol2part)**RTpar2%HOaqua
        k(RTind2%iHOaqua) = k(RTind2%iHOaqua) * AquaFac
      END IF

      !=== Compute mass transfer coefficient 
      IF ( nHENRY > 0 ) THEN
        thenry = MassTransfer( k(RTind2%iHENRY(:,1)), T, LWC )
        k(RTind2%iHENRY(:,1)) = thenry(:,1,1)
        k(RTind2%iHENRY(:,3)) = thenry(:,2,1)
      END IF

      ! ==== Law of mass action productories
      Prod = MassActionProducts( Conc )
      
      Rate  = Meff * k * Prod
      
      !WRITE(*,*) ''
      !DO j=1,neq
      !  WRITE(*,*) 'j,Time,Meff,k,prd,Rate =', &
      !  &           j,Time,vMeff(j),vK(j),vProd(j),Rate(j)
      !END DO

      TimeRates = TimeRates + MPI_WTIME() - TimeRateA
      
    END SUBROUTINE ReactionRates_Tropos


    FUNCTION MassActionProducts(Conc) RESULT(Prod)
      REAL(dp) :: Prod(neq)
      REAL(dp) :: Conc(nspc)
      INTEGER :: i

      Prod = ONE
      !
      ! stoechometric coefficients equal 1
      DO i=1,nFirst_order
        Prod(first_order(i,1)) = Prod(first_order(i,1)) &
         &                     * Conc(first_order(i,2))
      END DO
      !
      ! stoechometric coefficients equal 2
      DO i=1,nSecond_order
        Prod(second_order(i,1)) = Prod(second_order(i,1))    &
        &                       * Conc(second_order(i,2))    &
        &                       * ABS(Conc(second_order(i,2)))
      END DO
      !
      ! stoechometric coefficients not equal 1 or 2
      DO i=1,nHigher_order
        Prod(higher_order(i,1)) =  Prod(higher_order(i,1)) &
        &                       *  Conc(higher_order(i,2)) &
        &                       ** ahigher_order(i)
      END DO
      !
      ! if there are passive (katalytic) species e.g. [N2], [O2] or [aH2O]
      DO i=1,nfirst_orderKAT
        Prod(first_orderKAT(i,1)) = Prod(first_orderKAT(i,1))     & 
        &                         * InitValKat(first_orderKAT(i,2)) 
      END DO

    END FUNCTION MassActionProducts

    FUNCTION MassTransfer(kin,Temp,LWC) RESULT(k)
      REAL(dp) :: k(nHENRY,2,ntFrac), kin(nHENRY)
      REAL(dp) :: Temp(:), LWC
      ! TEMO
      REAL(dp) :: kmt(nHENRY,ntFrac)
      REAL(dp) :: term_diff(nHENRY), term_accom(nHENRY)
      REAL(dp) :: wetRadius(ntFrac)
      INTEGER :: i
      !
      !
      !---------------------------------------------------------------------------
      term_diff  = henry_diff(  RTind2%iHENRY(:,2) )               ! diffusion term
      term_accom = henry_accom( RTind2%iHENRY(:,2) ) * Temp(10)  ! accom term
      !--------------------------------------------------------------------------!
      !
      ! Compute new wet radius for droplett class iFrac
      wetRadius(:) = (Pi34*LWC/Frac%Number(:))**(rTHREE)
      !
      !--  mass transfer coefficient
      kmt = dkmt  ! set minimal transfer coefficient
      FORALL ( i = 1:nHENRY , term_diff(i) /= ZERO )
        kmt(i,:) = ONE / ( term_diff(i)*wetRadius(:)*wetRadius(:) + term_accom(i)*wetRadius(:) )
      END FORALL

      ! direaction GasSpecies-->AquaSpecies
      k(:,1,:) = milli * kmt(:,:) * LWC

      ! direaction AquaSpecies-->GasSpecies  
      DO i=1,ntFrac
        k(:,2,i) = kmt(:,i) / ( kin(:) * GasConst_R * Temp(1))  ! (...) = HenryConst*GasConstant*Temperatur
      END DO

    END FUNCTION MassTransfer
    
    FUNCTION EffectiveMolecularity(Conc) RESULT(M)
      !OUT
      REAL(dp) :: M(neq)
      !IN
      REAL(dp) :: Conc(:)
      !
      M = ONE

      IF(nFAC_H2>0)    M(RTind2%iFAC_H2)   = ((mH2*mair)**fac_exp)*fac_A
      IF(nFAC_O2N2>0)  M(RTind2%iFAC_O2N2) = (((mO2*mair)*(mN2*mair))**fac_exp)*fac_A
      IF(nFAC_M>0)     M(RTind2%iFAC_M)    = (mair**fac_exp)*fac_A
      IF(nFAC_O2>0)    M(RTind2%iFAC_O2)   = ((mO2*mair)**fac_exp)*fac_A
      IF(nFAC_N2>0)    M(RTind2%iFAC_N2)   = ((mN2*mair)**fac_exp)*fac_A
      IF(nFAC_H2O>0)   M(RTind2%iFAC_H2O)  = (mH2O**fac_exp)*fac_A
      IF(nFAC_RO2>0)   M(RTind2%iFAC_RO2)  = SUM(Conc(RO2))
      IF(nFAC_O2O2>0)  M(RTind2%iFAC_O2O2) = (((mO2*mair)*(mO2*mair))**fac_exp)*fac_A
      !IF(nFAC_aH2O>0) M(RTind2%iFAC_aH2O) = aH2OmolperL*LWC*mol2part
      IF(nFAC_RO2aq>0) M(RTind2%iFAC_RO2aq) = SUM(Conc(RO2aq))

    END FUNCTION EffectiveMolecularity

    SUBROUTINE ComputeRateConstant(k,T,Time,chi,mAir,Conc,Meff)
      REAL(dp), INTENT(OUT) :: k(neq)
      REAL(dp), INTENT(IN) :: Time, mAir, chi(:)
      REAL(dp), INTENT(IN) :: T(:)
      REAL(dp), INTENT(IN) :: Conc(:)
      REAL(dp), INTENT(INOUT) :: Meff(neq)

      REAL(dp) :: k_DC(nDCONST,2), k_T1(nDTEMP,2), k_T2(nDTEMP2,2), k_T3(nDTEMP3,2)
      REAL(dp) :: k_T4(nDTEMP4,2), k_T5(nDTEMP5,2), mesk(nMeskhidze,2)

      k = ZERO

      ! *** Photolytic reactions
      IF (nPHOTAB>0)  k(RTind2%iPHOTab)  = vPhoABCompute ( Time, chi )
      IF (nPHOTabc>0) k(RTind2%iPHOTabc) = vPhoABCCompute( Time, chi )
      IF (nPHOTMCM>0) k(RTind2%iPHOTmcm) = vPhoMCMCompute( Time, chi )

      ! *** Constant reactions
      IF (nCONST>0) k(RTind2%iCONST) = vConstCompute( )

      ! *** Temperature dependend reaction
      IF (nTEMP1>0) k(RTind2%iTEMP1) = vTemp1Compute( T )
      IF (nTEMP2>0) k(RTind2%iTEMP2) = vTemp2Compute( T )
      IF (nTEMP3>0) k(RTind2%iTEMP3) = vTemp3Compute( T )
      IF (nTEMP4>0) k(RTind2%iTEMP4) = vTemp4Compute( T )

      ! *** specieal aqua reactions
      IF (nASPEC1>0) k(RTind2%iASPEC1) = vAspec1Compute( T, Conc(Hp_ind) )
      IF (nASPEC2>0) k(RTind2%iASPEC2) = vAspec2Compute( T, Conc(Hp_ind) )
      IF (nASPEC3>0) k(RTind2%iASPEC3) = vAspec3Compute( T, Conc(Hp_ind) )
      
      ! *** dissociation reactions
      IF (nDCONST>0) THEN
        k_DC = vDConstCompute( )
        k(RTind2%iDCONST(:,1)) = k_DC(:,1) ! forward reactions
        k(RTind2%iDCONST(:,2)) = k_DC(:,2) ! backward reactions
      END IF
      IF (nDTEMP>0)  THEN 
        k_T1 = vDTempCompute( T )
        k(RTind2%iDTEMP(:,1)) = k_T1(:,1) ! forward reactions
        k(RTind2%iDTEMP(:,2)) = k_T1(:,2) ! backward reactions
      END IF
      IF (nDTEMP2>0) THEN 
        k_T2 = vDTemp2Compute( T )
        k(RTind2%iDTEMP2(:,1)) = k_T2(:,1) ! forward reactions
        k(RTind2%iDTEMP2(:,2)) = k_T2(:,2) ! backward reactions
      END IF
      IF (nDTEMP3>0) THEN 
        k_T3 = vDTemp3Compute( T )
        k(RTind2%iDTEMP3(:,1)) = k_T3(:,1) ! forward reactions
        k(RTind2%iDTEMP3(:,2)) = k_T3(:,2) ! backward reactions
      END IF
      IF (nDTEMP4>0) THEN 
        k_T4 = vDTemp4Compute( T )
        k(RTind2%iDTEMP4(:,1)) = k_T4(:,1) ! forward reactions
        k(RTind2%iDTEMP4(:,2)) = k_T4(:,2) ! backward reactions
      END IF
      IF (nDTEMP5>0) THEN 
        k_T5 = vDTemp5Compute( T )
        k(RTind2%iDTEMP5(:,1)) = k_T5(:,1) ! forward reactions
        k(RTind2%iDTEMP5(:,2)) = k_T5(:,2) ! backward reactions
      END IF
      IF (nMeskhidze>0) THEN
        mesk = vMeskhidzeCompute( T )
        k(RTind2%iMeskhidze(:,1)) = mesk(:,1) ! forward reactions
        k(RTind2%iMeskhidze(:,2)) = mesk(:,2) ! backward reactions
      END IF
     
      ! *** Troe reactions
      IF (nTROE>0)    k(RTind2%iTROE)    = vTroeCompute   ( T, mAir )
      IF (nTROEQ>0)   k(RTind2%iTROEq)   = vTroeEqCompute ( T, mAir )
      IF (nTROEF>0)   k(RTind2%iTROEf)   = vTroeFCompute  ( T, mAir )
      IF (nTROEQF>0)  k(RTind2%iTROEqf)  = vTroeEqfCompute( T, mAir )
      IF (nTROEXP>0)  k(RTind2%iTROExp)  = vTroeXPCompute ( T, mAir )
      IF (nTROEMCM>0) k(RTind2%iTROEmcm) = vTroeMCMCompute( T, mAir )

      IF (nSPEC1>0) k(RTind2%iSPEC1) = vSpec1Compute(    mAir )
      IF (nSPEC2>0) k(RTind2%iSPEC2) = vSpec2Compute( T, mAir )
      IF (nSPEC3>0) k(RTind2%iSPEC3) = vSpec3Compute( T, mAir )
      IF (nSPEC4>0) k(RTind2%iSPEC4) = vSpec4Compute( T, mAir )

      IF (nSPEC1MCM>0) k(RTind2%iSPEC1mcm) = vSPEC1MCMCompute( T, mAir )
      IF (nSPEC2MCM>0) k(RTind2%iSPEC2mcm) = vSPEC2MCMCompute( T       )
      IF (nSPEC3MCM>0) k(RTind2%iSPEC3mcm) = vSPEC3MCMCompute(    mAir )
      IF (nSPEC4MCM>0) k(RTind2%iSPEC4mcm) = vSpec4MCMCompute( T       )
      IF (nSPEC5MCM>0) k(RTind2%iSPEC5mcm) = vSpec5MCMCompute( T, mAir )
      IF (nSPEC6MCM>0) k(RTind2%iSPEC6mcm) = vSpec6MCMCompute( T       )
      IF (nSPEC7MCM>0) k(RTind2%iSPEC7mcm) = vSpec7MCMCompute( T       )
      IF (nSPEC8MCM>0) k(RTind2%iSPEC8mcm) = vSpec8MCMCompute( T, mAir )

      IF (nT1H2O>0) k(RTind2%iT1H2O) = vT1H2OCompute( T )
      IF (nS4H2O>0) k(RTind2%iS4H2O) = vS4H2OCompute( T, mAir)
      
      ! *** KKP photolytic reactions
      IF (nPHOTOkpp>0)  k(RTind2%iPHOTOkpp)  = vPhotokppCompute ( Time, Chi(3) )
      IF (nPHOTO2kpp>0) k(RTind2%iPHOTO2kpp) = vPhoto2kppCompute( Time, Chi(3) )
      IF (nPHOTO3kpp>0) k(RTind2%iPHOTO3kpp) = vPhoto3kppCompute( Time, Chi(3) )

    END SUBROUTINE ComputeRateConstant
    

    !=========================================================================!
    !                  calculate sun 
    !=========================================================================!
    FUNCTION Zenith(Time) RESULT(Chi)
      !-----------------------------------------------------------------------!
      ! Input:
      !   - Time
      REAL(dp) :: Time
      !-----------------------------------------------------------------------!
      ! Output:
      !   - sun angle chi
      REAL(dp) :: Chi
      !-----------------------------------------------------------------------!
      ! Temporary variables:
      !INTEGER :: IDAT
      REAL(dp) :: LBGMT, LZGMT
      REAL(dp) :: ML
      ! 
      REAL(dp) :: GMT
      REAL(dp) :: RLT, RPHI
      !    
      INTEGER        :: IIYEAR, IYEAR, IMTH, IDAY, IIY, NYEARS, LEAP, NOLEAP
      REAL(dp) :: YREF,YR
      !   
      INTEGER        :: I, IJ, JD, IJD, IN
      REAL(dp) :: D, RML, W, WR, EC, EPSI, PEPSI, YT, CW, SW, SSW  & 
      &                  , EYT, FEQT1, FEQT2, FEQT3, FEQT4, FEQT5, FEQT6 &
      &                  , FEQT7, FEQT, EQT
      !         
      REAL(dp) :: REQT, RA, RRA, TAB, RDECL, DECL, ZPT, CSZ, ZR    &
      &                     , CAZ, RAZ, AZIMUTH
      !           
      INTEGER :: IMN(12)
      DATA IMN/31,28,31,30,31,30,31,31,30,31,30,31/
      !
      !----------------------------------------------------------------------!
      !
      ! set GMT
      GMT = Time / HOUR
      !
      !  convert to radians
      RLT = LAT*DR
      RPHI = LONG*DR
      !
      !  parse date
      IIYEAR = IDAT/10000
      IYEAR = 19*100 + IIYEAR
      IF (IIYEAR <= 50) IYEAR = IYEAR + 100 
      IMTH = (IDAT - IIYEAR*10000)/100
      IDAY = IDAT - IIYEAR*10000 - IMTH*100
      !
      !  identify and correct leap years
      IIY = (IIYEAR/4)*4
      IF(IIY.EQ.IIYEAR) IMN(2) = 29
      !
      !  count days from Dec.31,1973 to Jan 1, YEAR, then add to 2,442,047.5
      YREF =  2442047.5_dp
      NYEARS = IYEAR - 1974
      LEAP = (NYEARS+1)/4
      IF(NYEARS.LE.-1) LEAP = (NYEARS-2)/4
      NOLEAP = NYEARS - LEAP
      YR = YREF + 365.0_dp*NOLEAP + 366.0_dp*LEAP
      !
      IJD = 0
      IN = IMTH - 1
      IF(IN.EQ.0) GO TO 40
      DO 30 I=1,IN
      IJD = IJD + IMN(I)
    30   CONTINUE
      IJD = IJD + IDAY
      GO TO 50
    40   IJD = IDAY
    50   IJ = IYEAR - 1973
      !
      !      print julian days current "ijd"
      JD = IJD + (YR - YREF)
      D = JD + GMT/24.0_dp
      !
      !      calc geom mean longitude
      ML = 279.2801988_dp + .9856473354_dp*D + 2.267e-13_dp*D*D
      RML = ML*DR
      !
      !      calc equation of time in sec
      !      w = mean long of perigee
      !      e = eccentricity
      !      epsi = mean obliquity of ecliptic
      W = 282.4932328_dp + 4.70684e-5_dp*D + 3.39e-13_dp*D*D
      WR = W*DR
      EC = 1.6720041e-2_dp - 1.1444e-9_dp*D - 9.4e-17_dp*D*D
      EPSI = 23.44266511_dp - 3.5626e-7_dp*D - 1.23e-15_dp*D*D
      PEPSI = EPSI*DR
      YT = (TAN(PEPSI*rTWO))**2
      CW = COS(WR)
      SW = SIN(WR)
      SSW = SIN(TWO*WR)
      EYT = TWO*EC*YT
      FEQT1 = SIN(RML)*(-EYT*CW - TWO*EC*CW)
      FEQT2 = COS(RML)*(TWO*EC*SW - EYT*SW)
      FEQT3 = SIN(TWO*RML)*(YT - (FIVE*EC*EC*rFOUR)*(CW*CW-SW*SW))
      FEQT4 = COS(TWO*RML)*(FIVE*EC**2*SSW*rFOUR)
      FEQT5 = SIN(THREE*RML)*(EYT*CW)
      FEQT6 = COS(THREE*RML)*(-EYT*SW)
      FEQT7 = -SIN(FOUR*RML)*(rTWO*YT*YT)
      FEQT = FEQT1 + FEQT2 + FEQT3 + FEQT4 + FEQT5 + FEQT6 + FEQT7
      EQT = FEQT*13751.0_dp
      !
      !   convert eq of time from sec to deg
      REQT = EQT/240.0_dp
      !
      !   calc right ascension in rads
      RA = ML - REQT
      RRA = RA*DR
      !
      !   calc declination in rads, deg
      TAB = 0.43360_dp*SIN(RRA)
      RDECL = ATAN(TAB)
      DECL = RDECL/DR
      !
      !   calc local hour angle
      LBGMT = 12.0_dp - EQT/HOUR + LONG*24.0_dp/360.0_dp
      LZGMT = 15.0_dp*(GMT - LBGMT)
      ZPT = LZGMT*DR
      CSZ = SIN(RLT)*SIN(RDECL) + COS(RLT)*COS(RDECL)*COS(ZPT)
      ZR = ACOS(CSZ)
      ! 
      !   calc local solar azimuth
      CAZ = (SIN(RDECL) - SIN(RLT)*COS(ZR))/(COS(RLT)*SIN(ZR))
      RAZ = ACOS(CAZ)
      AZIMUTH = RAZ/DR
      !
      !--- set Zenith Angle
      Chi =  1.745329252e-02_dp * ZR/DR
    END FUNCTION Zenith

    FUNCTION UpdateSun(Time) RESULT(Sun)
      !--------------------------------------------------------------------
      ! Input:
      !   - Time
      REAL(dp) :: Time
      !--------------------------------------------------------------------!
      ! Output:
      !   - Sun
      REAL(dp)  :: Sun
      !--------------------------------------------------------------------!
      ! Temporary variables:
      REAL(dp), PARAMETER :: SunRise=4.50_dp, SunSet=19.50_dp
      REAL(dp) :: Thour, Tlocal, Ttmp
      !
      Thour  = Time / HOUR
      Tlocal = Thour - FLOOR(Thour/hourday)*hourday
      !
      IF( (Tlocal>=SunRise) .AND. (Tlocal<=SunSet) ) THEN
        Ttmp = (TWO*Tlocal-SunRise-SunSet) / (SunSet-SunRise);
        IF ( Ttmp >ZERO ) THEN
          Ttmp =  Ttmp * Ttmp
        ELSE
          Ttmp = -Ttmp * Ttmp
        END IF
        Sun = (ONE+COS(Pi*Ttmp)) * rTWO
      ELSE
        Sun = ZERO
      END IF
    END FUNCTION UpdateSun
   
  !========================================================================!
  !                          GASEOUS REACTIONS                             !
  !                          -----------------                             !
  !========================================================================!
  ! ===  Photolysis reactions                                              !
  !========================================================================!
  ! Input:                                                                 !
  !   - Contants                                                           !
  !   - Time                                                               !
  !------------------------------------------------------------------------!
  ! Output:                                                                !
  !   - Reaction Constant                                                  !
  !------------------------------------------------------------------------!

    FUNCTION vPhoABCCompute( Time, chi ) RESULT(kPHOTabc)
      REAL(dp)  :: kPHOTabc(nPHOTabc)
      REAL(dp), INTENT(IN)  :: Time, Chi(:)
      REAL(dp), DIMENSION(nPHOTabc) :: ChiZ, yChiZ, EyChiZ
      INTEGER :: i, j
      
      IF ( Chi(1) < PiHalf ) THEN
        ChiZ = Chi(1) * RTpar2%PHOTabc(:,3) 
        DO i = 1,nPHOTabc
          IF (ChiZ(i) < PiHalf) THEN
            yChiZ(i) = RTpar2%PHOTabc(i,2) * (One - One/COS(ChiZ(i)))
            IF ( yChiZ(i) > mTHIRTY ) THEN
              EyChiZ(i) = EXP(yChiZ(i))
            ELSE
              EyChiZ(i) = EyChiZmin   ! = 9.357d-14  
            END IF
          ELSE
            EyChiZ(i) = EyChiZmin   ! = 9.357d-14 
          END IF
        END DO
        kPHOTabc = Dust * RTpar2%PHOTabc(:,1) * EyChiz
      ELSE
        kPHOTabc = ZERO
      END IF
    END FUNCTION vPhoABCCompute

    FUNCTION vPhoABCompute(Time,Chi) RESULT(kPHOTab)
      REAL(dp)  :: kPHOTab(nPHOTab)
      REAL(dp), INTENT(IN)  :: Time, Chi(:)
     
      IF ( Chi(1) < PiHalf ) THEN
        kPHOTab = Dust * RTpar2%PHOTab(:,1)*EXP(-RTpar2%PHOTab(:,2)*Chi(2))
      ELSE
        kPHOTab = ZERO
      END IF
    END FUNCTION vPhoABCompute
    
    FUNCTION vPhoMCMCompute(Time,Chi) RESULT(kPHOTmcm)
      REAL(dp)  :: kPHOTmcm(nPHOTmcm)
      REAL(dp), INTENT(IN)  :: Time, Chi(:)
      REAL(dp), DIMENSION(nPHOTmcm) :: ChiZ, yChiZ
      
      !---  MCM version
      IF ( Chi(1) < PiHalf ) THEN
        ChiZ  = EXP( -RTpar2%PHOTmcm(:,3) * Chi(2) )
        yChiZ = Chi(2) ** RTpar2%PHOTmcm(:,2)
        kPHOTmcm = Dust * RTpar2%PHOTmcm(:,1) * yChiZ * ChiZ
      ELSE
        kPHOTmcm = ZERO
      END IF
    END FUNCTION vPhoMCMCompute
    
    !
    !==========================================================================!
    ! ===  Constant reactions
    !==========================================================================!
    FUNCTION vConstCompute() RESULT(kCONST)
      REAL(dp)  :: kCONST(nCONST)
      kCONST = RTpar2%CONST
    END FUNCTION vConstCompute
    !
    !
    !==========================================================================!
    ! ===  Temperature-Dependent  (Arrhenius)
    !==========================================================================!
    ! Input: 
    !   - Contants
    !   - Temperature
    !--------------------------------------------------------------------------!
    ! Output:
    !   - Reaction constant
    !--------------------------------------------------------------------------!
    FUNCTION vTemp1Compute(Temp) RESULT(kTEMP1)
      REAL(dp)  :: kTEMP1(nTEMP1)
      REAL(dp), INTENT(IN)  :: Temp(:)
      
      kTEMP1 = RTpar2%TEMP1(:,1) * EXP(-RTpar2%TEMP1(:,2)*Temp(6))
    END FUNCTION vTemp1Compute
    
    FUNCTION vTemp2Compute(Temp) RESULT(kTEMP2)
      REAL(dp)  :: kTEMP2(nTEMP2)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      kTEMP2 = RTpar2%TEMP2(:,1) * Temp(2) * EXP( -RTpar2%TEMP2(:,2)*Temp(6) )
    END FUNCTION vTemp2Compute
    
    FUNCTION vTemp3Compute(Temp) RESULT(kTEMP3)
      REAL(dp)  :: kTEMP3(nTEMP3)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      kTEMP3 = RTpar2%TEMP3(:,1) * EXP( RTpar2%TEMP3(:,2)*(Temp(6) - InvRefTemp) )
    END FUNCTION vTemp3Compute
    
    FUNCTION vTemp4Compute(Temp) RESULT(kTEMP4)
      REAL(dp)  :: kTEMP4(nTEMP4)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      kTEMP4 = RTpar2%TEMP4(:,1) * Temp(1) * EXP( -RTpar2%TEMP4(:,2)*Temp(6) )
    END FUNCTION vTemp4Compute
    
   
    !   ***************************************************************
    !   ** Species nondimensional gibbs potentials                   **
    !   ***************************************************************
    SUBROUTINE GibbsFreeEnergie(Gibbs,T)
      REAL(dp) :: Gibbs(:)
      REAL(dp) :: T(:)
      !
      Gibbs(:)=ZERO
      ! WHERE wird bald abgeschaft, -> vektorisieren
      WHERE (SwitchTemp>T(1))
        Gibbs = lowA*(ONE-T(8)) - rTWO*lowB*T(1) - rSIX*lowC*T(2)        &
        &        - rTWELV*lowD*T(3) - rTWENTY*lowE*T(4) + lowF*T(6) - lowG 
      ELSEWHERE
        Gibbs = highA*(ONE-T(8)) - rTWO*highB*T(1) - rSIX*highC*T(2)         &
        &        - rTWELV*highD*T(3) - rTWENTY*highE*T(4) + highF*T(6) - highG 
      END WHERE
    END SUBROUTINE GibbsFreeEnergie
    !
    !
    SUBROUTINE CalcDiffGibbsFreeEnergie(DGibbsdT,T)
      REAL(dp) :: DGibbsdT(:)
      REAL(dp) :: T(:)
      !
      DGibbsdT(:)=ZERO
      ! WHERE wird bald abgeschaft, -> vektorisieren
      WHERE (SwitchTemp>T(1))
        DGibbsdT=-(lowA*T(6) + 0.5d0*lowB + lowC*T(1)/3.0d0 +            &
        &              0.25d0*lowD*T(2) + 0.2d0*lowE*T(3) + lowF*T(7) )
      ELSEWHERE
        DGibbsdT=-(highA*T(6) + 0.5d0*highB + highC*T(1)/3.0d0 +         &
        &              0.25d0*highD*T(2) + 0.2d0*highE*T(3) + highF*T(7) )
      END WHERE  
    END SUBROUTINE CalcDiffGibbsFreeEnergie
    !
    !
    SUBROUTINE CalcDeltaGibbs(DelGibbs)
      REAL(dp) :: DelGibbs(:)
      !
      INTEGER :: iR
      INTEGER :: from, to
      !
      DelGibbs = ZERO         
      !
      DO iR=1,neq
        from = BA%RowPtr(iR);   to = BA%RowPtr(iR+1)-1
        DelGibbs(iR) = DelGibbs(iR)   &
        &            - SUM( BA%Val(from:to) * GFE(BA%ColInd(from:to)) )
      END DO

    END SUBROUTINE CalcDeltaGibbs
    !
    SUBROUTINE CalcDiffDeltaGibbs(DiffDelGibbs)
      REAL(dp) :: DiffDelGibbs(:)
      !
      INTEGER :: iR, jS, jj
      !
      DiffDelGibbs(:)=ZERO         
      DO iR=1,BA%m
        DO jj=BA%RowPtr(iR),BA%RowPtr(iR+1)-1
          jS=BA%ColInd(jj)
          DiffDelGibbs(iR)=DiffDelGibbs(iR)+BA%Val(jj)*DGFEdT(jS)      
        END DO
      END DO
    END SUBROUTINE CalcDiffDeltaGibbs
    !
    !
    SUBROUTINE scTHERMO(C,H,S,T)
      !
      ! IN
      REAL(dp) :: T(8)
      !
      ! OUT
      REAL(dp) :: C(nspc)       ! molar heat capacities at constant pressure
      REAL(dp) :: H(nspc)       ! the standardstate molar enthalpy
      REAL(dp) :: S(nspc)       ! standard-state entropy at 298 K
      !
      REAL(dp) :: dHdT(nspc)    ! Enthaply derivative in dT [J/mol/K^2]
      REAL(dp) :: dGdT(nspc)    ! Gibbs potential derivative in dT [J/mol/K^2]
      REAL(dp) :: dCvdT(nspc)   ! Constant volume specific heat derivative in dT [J/mol/K]
      
      WHERE (SwitchTemp>T(1))
        C = lowA + lowB*T(1) + lowC*T(2) + lowD*T(3) + lowE*T(4)
        H = lowA + rTWO*lowB*T(1) + rTHREE*lowC*T(2) + rFOUR*lowD*T(3)    &
        &        + rFIVE*lowE*T(4) + lowF*T(6)
        S = lowA*T(8) + lowB*T(1) + rTWO*lowC*T(2) + rTHREE*lowD*T(3)     &
        &        + rFIVE*lowE*T(4) + lowG
        !
        dCvdT = R * (highB + TWO*highC*T(1) + THREE*highD*T(2) + FOUR*highE*T(3))
        dHdT = R * (highA + highB*T(1) + highC*T(2) + highD*T(3) + highE*T(4))
        dGdT = - R * (highG + highA*T(8) + highB*T(1) + rTWO*highC*T(2)      &
        &                   + rTHREE*highD*T(3) +rFOUR*highE*T(4))
      ELSEWHERE
        C = highA + highB*T(1) + highC*T(2) + highD*T(3) + highE*T(4)
        H = highA + rTWO*highB*T(1) + rTHREE*highC*T(2) + rFOUR*highD*T(3)    &
        &         + rFIVE*highE*T(4) + highF*T(6)
        S = highA*T(8) + highB*T(1) + rTWO*highC*T(2) + rTHREE*highD*T(3)     &
        &         + rFIVE*highE*T(4) + highG
        !
        dCvdT = R * (highB + TWO*highC*T(1) + THREE*highD*T(2) + FOUR*highE*T(3))
        dHdT = R * (highA + highB*T(1) + highC*T(2) + highD*T(3) + highE*T(4))
        dGdT = - R * (highG + highA*T(8) + highB*T(1) + rTWO*highC*T(2)      &
        &                   + rTHREE*highD*T(3) +rFOUR*highE*T(4))
      END WHERE  

      !print*, 'debug::rates    H=',H
      !print*, 'debug::rates    C=',C
      !print*, 'debug::rates    S=',S

      !stop
    END SUBROUTINE scTHERMO


    FUNCTION TroeFactorVec(T)
      !OUT
      REAL(dp) :: TroeFactorVec(RTind%nTroe)
      !IN
      REAL(dp) :: T(:)
      !TEMP
      REAL(dp) :: Fcent(RTind%nTroe),     &
                      & logF_Troe(RTind%nTroe), &
                      & FACtroe(RTind%nTroe)


    
      Fcent = (ONE - RTpar%T1) * EXP(-T(1)/RTpar%T2) &
            &      + RTpar%T1  * EXP(-T(1)/RTpar%T3) &
            &      +             EXP(-T(6)*RTpar%T4)
      
      vlog10_Fcent = LOG10(Fcent)
      vlog10_Pr    = LOG10(vPr(RTind%iTroe)/(ONE-vPr(RTind%iTroe)))

      vcTroe  = -0.40e0_dp - 0.67e0_dp*vlog10_Fcent    ! will be used for deriv too 
      vn1Troe =  0.75e0_dp - 1.27e0_dp*vlog10_Fcent
      
      FACtroe = (vlog10_Pr+vcTroe) / (vn1Troe-dTroe*(vlog10_Pr+vcTroe))

      logF_Troe  = vlog10_Fcent / (ONE+FACtroe*FACtroe)
      TroeFactorVec = vPr(RTind%iTroe) * TEN**logF_Troe

      !print*,' vlog10_Fcent(1),vlog10_Pr(1),FACtroe(1),logF_Troe(1)'
      !print*, vlog10_Fcent(1),vlog10_Pr(1),FACtroe(1),logF_Troe(1)
      !print*, vcTroe(1), vn1Troe(1), dTroe
      !stop
    
    END FUNCTION TroeFactorVec

    FUNCTION DTroeFactorVec(Dk0dT,DkinfdT,dFTL_dT,T)
      REAL(dp) :: DTroeFactorVec(RTind%nTroe)
      !IN
      REAL(dp) :: Dk0dT(:), DkinfdT(:), dFTL_dT(:)
      REAL(dp) :: T(:)
      !TEMP
      REAL(dp), DIMENSION(RTind%nTroe) :: dlog10_PrdT, log10_Prc, tmplog10
      REAL(dp), DIMENSION(RTind%nTroe) :: DlogF_Troedlog_Pr, Dlog_F_TroedT
      REAL(dp), DIMENSION(RTind%nTroe) :: log10_FTroe

      ! Troe mechanism

      dlog10_PrdT = rln10*(Dk0dT(RTind%iTroe) - DkinfdT(RTind%iTroe))*T(6)

      log10_Prc   = vlog10_Pr + cTroe

      DlogF_Troedlog_Pr = - TWO*n1Troe*vlog10_Fcent*log10_Prc             &
                        &       * ( n1Troe - dTroe*log10_Prc )           &
                        &  / (  n1Troe**2 - TWO*n1Troe*dTroe*log10_Prc + &
                        &      (dTroe**2 + ONE)*log10_Prc*log10_Prc  )**2

      Dlog_F_TroedT  = DlogF_Troedlog_Pr * dlog10_PrdT
      tmplog10       = log10_Prc   / (n1Troe - dTroe*log10_Prc)
      log10_FTroe    = vlog10_Fcent / (ONE   + (tmplog10*tmplog10))

      DTroeFactorVec = TEN**log10_FTroe * ( dFTL_dT(RTind%iTroe) &
      &              + vPr(RTind%iTroe) * ln10*Dlog_F_TroedT )

    END FUNCTION DTroeFactorVec

    FUNCTION UpdateTempArray(Temperature) RESULT(TempArr)
      REAL(dp) :: Temperature 
      REAL(dp) :: TempArr(10)
      !
      INTEGER :: i
      !
      TempArr(1) = Temperature                ! T
      DO i=2,5
        TempArr(i) = TempArr(i-1)*Temperature   ! T^2 ... T^5
      END DO
      TempArr(6)  = ONE / Temperature          ! 1/T
      TempArr(7)  = ONE / TempArr(2)           ! 1/T^2
      TempArr(8)  = LOG(Temperature)           ! ln(T)
      TempArr(9)  = SQRT(Temperature)           ! ln(T)
      TempArr(10) = ONE / TempArr(9)
    END FUNCTION UpdateTempArray
    !
    !
    !-------------------------------------------------------------------------
    !---  Species internal energies in moles [J/mol]  
    !-------------------------------------------------------------------------
    SUBROUTINE InternalEnergy(U,T)
      !OUT
      REAL(dp) :: U(nspc)   
      !IN
      REAL(dp) :: T(:)
      
      WHERE (SwitchTemp>T(1))
        U  = (  (lowA-ONE)*T(1) + rTWO*lowB*T(2)  + rTHREE*lowC*T(3)    & 
        &     + rFOUR*lowD*T(4) + rFIVE*lowE*T(5) + lowF )
      ELSEWHERE
        U  = (  (highA-ONE)*T(1) + rTWO*highB*T(2)  + rTHREE*highC*T(3) & 
        &     + rFOUR*highD*T(4) + rFIVE*highE*T(5) + highF )
      END WHERE
      U    = U * R
    END SUBROUTINE InternalEnergy
    !
    !-------------------------------------------------------------------------------
    !--- Nondimensionl Derivatives of specific heats at constant volume in moles [-]
    !-------------------------------------------------------------------------------
    SUBROUTINE DiffInternalEnergy(dUdT,T)
      !OUT
      REAL(dp) :: dUdT(nspc)   
      !IN
      REAL(dp) :: T(:)                      
      !
      WHERE (SwitchTemp>T(1))
        dUdT  = (lowA + lowB*T(1) + lowC*T(2)     & 
        &             + lowD*T(3) + lowE*T(4)) 
      ELSEWHERE
        dUdT  = (highA + highB*T(1) + highC*T(2)  & 
        &              + highD*T(3) + highE*T(4))  
      END WHERE
      dUdT = (dUdT - ONE)        ! speedchem SCmodule.f ~2618
    END SUBROUTINE DiffInternalEnergy
    !
    !
    SUBROUTINE Diff2InternalEnergy(d2UdT2,T)
      REAL(dp) :: d2UdT2(nspc)     !Constant volume specific heat’s derivative [J/mol/K2] 
      REAL(dp) :: T(:)                      
      !
      WHERE (SwitchTemp>T(1))
        d2UdT2 = lowB + TWO*lowC*T(1) + THREE*lowD*T(2) + FOUR*lowE*T(3) 
      ELSEWHERE    
        d2UdT2 = highB + TWO*highC*T(1) + THREE*highD*T(2) + FOUR*highE*T(3)
      END WHERE
      d2UdT2 = d2UdT2 * R
    END SUBROUTINE Diff2InternalEnergy
    
    !--------------------------------------------------------------------------!
    
    FUNCTION vTroeCompute(Temp,mAir) RESULT(kTROE)
      REAL(dp)  :: kTROE(nTROE)
      REAL(dp), INTENT(IN)  :: Temp(:) , mAir
      REAL(dp), DIMENSION(nTROE) :: k1, k2, log10_k1k2
      REAL(dp), PARAMETER   :: F = 0.6_dp
      
      k1 = RTpar2%TROE(:,1) * (Temp(1)*r300)**(-RTpar2%TROE(:,2)) * mAir
      k2 = RTpar2%TROE(:,3) * (Temp(1)*r300)**(-RTpar2%TROE(:,4))
      log10_k1k2 = LOG10(k1/k2)
      kTROE = k1 / (One + k1/k2) * F**(One / (One + log10_k1k2*log10_k1k2) )
    END FUNCTION vTroeCompute
      
    FUNCTION vTroeMCMCompute(Temp,mAir) RESULT(kTROEmcm)
      REAL(dp)  :: kTROEmcm(nTROEmcm)
      REAL(dp), INTENT(IN)  :: Temp(:) , mAir
      REAL(dp), DIMENSION(nTROEmcm) :: k1, k2, Fc, tmpTROE
      REAL(dp), PARAMETER   :: n=0.75_dp, d=1.27_dp
      !
      k1 = RTpar2%TROEmcm(:,1)*(Temp(1)*InvRefTemp)**RTpar2%TROEmcm(:,2)  &
      &  * EXP(RTpar2%TROEmcm(:,3)*Temp(6))*mAir
      k2 = RTpar2%TROEmcm(:,4)*(Temp(1)*InvRefTemp)**RTpar2%TROEmcm(:,5)  &
      &  * EXP(RTpar2%TROEmcm(:,6)*Temp(6))
      Fc = RTpar2%TROEmcm(:,7)*EXP(RTpar2%TROEmcm(:,8) &
      &  * Temp(6))+RTpar2%TROEmcm(:,9)*EXP(Temp(1)/RTpar2%TROEmcm(:,10))
      tmpTROE  = LOG10(k1/k2)/(n - d*LOG10(fc))
      kTROEmcm = k1 / (One + k1/k2) * Fc**(One / (One + tmpTROE*tmpTROE))
    END FUNCTION vTroeMCMCompute

    !---  Troe with variable F factor
    FUNCTION vTroeFCompute(Temp,mAir) RESULT(kTROEf)  ! Barthel
      REAL(dp)  :: kTROEf(nTROEf)
      REAL(dp), INTENT(IN)  :: Temp(:) , mAir
      REAL(dp), DIMENSION(nTROEf) :: k1, k2, log10_k1k2
      !
      k1 = RTpar2%TROEf(:,1)*(Temp(1)*r300)**(-RTpar2%TROEf(:,2)) * mAir
      k2 = RTpar2%TROEf(:,3)*(Temp(1)*r300)**(-RTpar2%TROEf(:,4))
      log10_k1k2 = LOG10(k1/k2)
      kTROEf = k1 / (One+k1/k2) * RTpar2%TROEf(:,5)**(One / (One + log10_k1k2*log10_k1k2))
    END FUNCTION vTroeFCompute
    
    !---  Troe equilibrium
    FUNCTION vTroeEqCompute(Temp,mAir) RESULT(kTROEq)
      REAL(dp)  :: kTROEq(nTROEq)
      REAL(dp), INTENT(IN)  :: Temp(:) , mAir
      REAL(dp), DIMENSION(nTROEq) :: k1, k2, k3, log10_k1k2
      REAL(dp), PARAMETER   :: F = 0.6_dp
      !
      k1 = RTpar2%TROEq(:,1)*(Temp(1)*r300)**(-RTpar2%TROEq(:,2)) * mAir
      k2 = RTpar2%TROEq(:,3)*(Temp(1)*r300)**(-RTpar2%TROEq(:,4))
      log10_k1k2 = LOG10(k1/k2)
      k3 = k1 / (One + k1/k2) * F**(One / (One + log10_k1k2*log10_k1k2))
      kTROEq = k3 / ( RTpar2%TROEq(:,5) * EXP(RTpar2%TROEq(:,6)*Temp(6)) )          
    END FUNCTION vTroeEqCompute
    
    FUNCTION vTroeEqfCompute(Temp,mAir) RESULT(kTROEqf) ! Barthel
      REAL(dp)  :: kTROEqf(nTROEqf)
      REAL(dp), INTENT(IN)  :: Temp(:) , mAir
      REAL(dp), DIMENSION(nTROEqf) :: k1, k2, k3, log10_k1k2
      ! 
      k1 = RTpar2%TROEqf(:,1) * (Temp(1)*r300)**(-RTpar2%TROEqf(:,2)) * mAir
      k2 = RTpar2%TROEqf(:,3) * (Temp(1)*r300)**(-RTpar2%TROEqf(:,4))
      log10_k1k2 = LOG10(k1/k2)
      k3 = k1 / (One + k1/k2) * RTpar2%TROEqf(:,7)**(One / (One + log10_k1k2*log10_k1k2))
      kTROEqf = k3 / (RTpar2%TROEqf(:,5) * EXP(RTpar2%TROEqf(:,6)*Temp(6)))          
    END FUNCTION vTroeEqfCompute

    !---  modified Troe with variable F factor
    FUNCTION vTroeXPCompute(Temp,mAir) RESULT(kTROExp)  ! Barthel
      REAL(dp)  :: kTROExp(nTROExp)
      REAL(dp), INTENT(IN)  :: Temp(:) , mAir
      REAL(dp), DIMENSION(nTROExp) :: k1, k2, log10_k1k2
      !
      k1 = RTpar2%TROExp(:,1)*EXP(-RTpar2%TROExp(:,2)*Temp(6)) * mAir
      k2 = RTpar2%TROExp(:,3)*EXP(-RTpar2%TROExp(:,4)*Temp(6))
      log10_k1k2 = LOG10(k1/k2)
      kTROExp = k1 / (One + k1/k2) * RTpar2%TROExp(:,5)**(One / (One + log10_k1k2*log10_k1k2)) 
    END FUNCTION vTroeXPCompute
    
    FUNCTION vSpec1Compute(mAir) RESULT(kSPEC1)
      REAL(dp)  :: kSPEC1(nSPEC1)
      REAL(dp), INTENT(IN)  :: mAir
      !
      kSPEC1 = RTpar2%SPEC1(:,1) * (ONE + mAir * RTpar2%SPEC1(:,2))
    END FUNCTION vSpec1Compute
    
    FUNCTION vSpec2Compute(Temp,mAir) RESULT(kSPEC2)
      REAL(dp)  :: kSPEC2(nSPEC2)
      REAL(dp), INTENT(IN)  :: Temp(:), mAir
      ! 
      kSPEC2 = mAir * RTpar2%SPEC2(:,1) * (Temp(1)*r300)**RTpar2%SPEC2(:,2)
    END FUNCTION vSpec2Compute
    
    FUNCTION vSpec3Compute(Temp,mAir) RESULT(kSPEC3)
      REAL(dp)  :: kSPEC3(nSPEC3)
      REAL(dp), INTENT(IN)  :: Temp(:), mAir
      REAL(dp), DIMENSION(nSPEC3) :: k1, k2 ,k3
      !
      k1 = RTpar2%SPEC3(:,1)*EXP(RTpar2%SPEC3(:,2)*Temp(6))
      k2 = RTpar2%SPEC3(:,3)*EXP(RTpar2%SPEC3(:,4)*Temp(6))
      k3 = RTpar2%SPEC3(:,5)*EXP(RTpar2%SPEC3(:,6)*Temp(6)) * mAir
      kSPEC3 = k1 + k3 / (One + k3/k2)
    END FUNCTION vSpec3Compute
    
    FUNCTION vSpec4Compute(Temp,mAir) RESULT(kSPEC4)
      REAL(dp)  :: kSPEC4(nSPEC4)
      REAL(dp), INTENT(IN)  :: Temp(:), mAir
      REAL(dp), DIMENSION(nSPEC4) :: k1, k2 ,k3
      !
      kSPEC4 = RTpar2%SPEC4(:,1) * EXP(RTpar2%SPEC4(:,2)*Temp(6)) +   &
      &        RTpar2%SPEC4(:,3) * EXP(RTpar2%SPEC4(:,4)*Temp(6)) * mAir 
    END FUNCTION vSpec4Compute
    
    FUNCTION vSpec1MCMCompute(Temp,mAir) RESULT(kSPEC1mcm)
      REAL(dp)  :: kSPEC1mcm(nSPEC1mcm)
      REAL(dp), INTENT(IN)  :: Temp(:), mAir
      
      kSPEC1mcm = RTpar2%SPEC1mcm(:,1)*(ONE + mAir*RTpar2%SPEC1mcm(:,2)*r300*Temp(1)/RTpar2%SPEC1mcm(:,3))
    END FUNCTION vSpec1MCMCompute
    
    FUNCTION vSpec2MCMCompute(Temp) RESULT(kSPEC2mcm)
      REAL(dp)  :: kSPEC2mcm(nSPEC2mcm)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      kSPEC2MCM = RTpar2%SPEC2mcm(:,1)*(Temp(1)*r300)**RTpar2%SPEC2mcm(:,2) &
                * EXP(RTpar2%SPEC2mcm(:,3)*Temp(6))
    END FUNCTION vSpec2MCMCompute
    
    FUNCTION vSpec3MCMCompute(mAir) RESULT(kSPEC3mcm)
      REAL(dp)  :: kSPEC3mcm(nSPEC3mcm)
      REAL(dp), INTENT(IN)  :: mAir
      !
      kSPEC3mcm = RTpar2%SPEC3mcm(:,1)*(ONE + mAir/RTpar2%SPEC3mcm(:,2))
    END FUNCTION vSpec3MCMCompute
    
    FUNCTION vSpec4MCMCompute(Temp) RESULT(kSPEC4mcm)
      REAL(dp)  :: kSPEC4mcm(nSPEC4mcm)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      kSPEC4mcm = RTpar2%SPEC4mcm(:,1)*(ONE + RTpar2%SPEC4mcm(:,2) &
      &         * EXP(RTpar2%SPEC4mcm(:,3)*Temp(6))*H2O)           &
      &         * EXP(RTpar2%SPEC4mcm(:,4)*Temp(6))
    END FUNCTION vSpec4MCMCompute
    
    FUNCTION vSpec5MCMCompute(Temp,mAir) RESULT(kSPEC5mcm)
      REAL(dp)  :: kSPEC5mcm(nSPEC5mcm)
      REAL(dp), INTENT(IN)  :: Temp(:), mAir
      REAL(dp), DIMENSION(nSPEC5mcm) :: k1, k2
      REAL(dp), PARAMETER   :: F = 0.21_dp
      !
      k1 = RTpar2%SPEC5mcm(:,1) * mAir * F * EXP(RTpar2%SPEC5mcm(:,2)*Temp(6))
      k2 = RTpar2%SPEC5mcm(:,3) * mAir * F * EXP(RTpar2%SPEC5mcm(:,4)*Temp(6))
      kSPEC5mcm = k1 * (One - k2)
    END FUNCTION vSpec5MCMCompute
    
    FUNCTION vSpec6MCMCompute(Temp) RESULT(kSPEC6mcm)
      REAL(dp)  :: kSPEC6mcm(nSPEC6mcm)
      REAL(dp), INTENT(IN)  :: Temp(:)
      REAL(dp), DIMENSION(nSPEC6mcm) :: k1, k2
      !
      k1 = RTpar2%SPEC6mcm(:,1)*EXP(RTpar2%SPEC6mcm(:,2)*Temp(6))
      k2 = RTpar2%SPEC6mcm(:,3)*EXP(RTpar2%SPEC6mcm(:,4)*Temp(6))
      kSPEC6mcm = k1 * (ONE - k2)
    END FUNCTION vSpec6MCMCompute
    
    FUNCTION vSpec7MCMCompute(Temp) RESULT(kSPEC7mcm)
      REAL(dp)  :: kSPEC7mcm(nSPEC7mcm)
      REAL(dp), INTENT(IN)  :: Temp(:)
      REAL(dp), DIMENSION(nSPEC7mcm) :: k1, k2
      !
      k1 = RTpar2%SPEC7mcm(:,1)*EXP(RTpar2%SPEC7mcm(:,2)*Temp(6))
      k2 = RTpar2%SPEC7mcm(:,3)*EXP(RTpar2%SPEC7mcm(:,4)*Temp(6))
      kSPEC7mcm = k1 * (RTpar2%SPEC7mcm(:,5) - RTpar2%SPEC7mcm(:,6) / (One + k2))
    END FUNCTION vSpec7MCMCompute
    
    FUNCTION vSpec8MCMCompute(Temp,mAir) RESULT(kSPEC8mcm)
      REAL(dp)  :: kSPEC8mcm(nSPEC8mcm)
      REAL(dp), INTENT(IN)  :: Temp(:), mAir
      REAL(dp), DIMENSION(nSPEC8mcm) :: k1, k2
      REAL(dp), PARAMETER   :: F = 0.21_dp
      !
      k1 = RTpar2%SPEC8mcm(:,1) * mAir * F * EXP(RTpar2%SPEC8mcm(:,2)*Temp(6))
      k2 = RTpar2%SPEC8mcm(:,3) * mAir * F * EXP(RTpar2%SPEC8mcm(:,4)*Temp(6))
      kSPEC8mcm = k1 / (One + k2) * Temp(6)
    END FUNCTION vSpec8MCMCompute
    
    FUNCTION vS4H2OCompute(Temp,mAir) RESULT(kS4H2O)
      REAL(dp)  :: kS4H2O(nS4H2O)
      REAL(dp), INTENT(IN)  :: Temp(:), mAir
      
      kS4H2O = RTpar2%S4H2O(:,1) * EXP(RTpar2%S4H2O(:,2)*Temp(6))      & 
      &      + RTpar2%S4H2O(:,3) * EXP(RTpar2%S4H2O(:,4)*Temp(6)) * mAir
    END FUNCTION vS4H2OCompute
    
    FUNCTION vPhotokppCompute(Time,Sun) RESULT(kPHOTO)
      REAL(dp) :: kPHOTO(nPHOTOkpp)
      REAL(dp) :: Time, Sun
      !
      kPHOTO = RTpar2%PHOTOkpp(:) * Sun
    END FUNCTION vPhotokppCompute
    
    FUNCTION vPhoto2kppCompute(Time,Sun) RESULT(kPHOTO2)
      REAL(dp) :: kPHOTO2(nPHOTO2kpp)
      REAL(dp) :: Time, Sun
      !
      kPHOTO2 = RTpar2%PHOTO2kpp(:) * Sun*Sun
    END FUNCTION vPhoto2kppCompute
    
    FUNCTION vPhoto3kppCompute(Time,Sun) RESULT(kPHOTO3)
      REAL(dp) :: kPHOTO3(nPHOTO3kpp)
      REAL(dp) :: Time, Sun
      !
      kPHOTO3 = RTpar2%PHOTO3kpp(:) * Sun*Sun*Sun
    END FUNCTION vPhoto3kppCompute
    !======================================================================!
    !                          AQUEOUS REACTIONS                           !
    !                          -----------------                           !
    !======================================================================!
    ! ===  Photolysis reactions 
    !======================================================================!
    ! Input:
    !   - Contants
    !   - Time
    !----------------------------------------------------------------------!
    ! Output:
    !   - Reaction Constant
    !----------------------------------------------------------------------!
    
    FUNCTION vDConstCompute() RESULT(k)
      REAL(dp)  :: k(nDCONST,2)
      !
      k(:,2) = RTpar2%DCONST(:,2)
      k(:,1) = RTpar2%DCONST(:,1) * k(:,2)
    END FUNCTION vDConstCompute
    
    FUNCTION vDTempCompute(Temp) RESULT(k)
      REAL(dp)  :: k(nDTEMP,2)
      REAL(dp), INTENT(IN) :: Temp(:)
      !
      k(:,2) = RTpar2%DTEMP(:,3)
      k(:,1) = RTpar2%DTEMP(:,1)*EXP(RTpar2%DTEMP(:,2)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vDTempCompute
    
    FUNCTION vDTemp2Compute(Temp)   RESULT(k)
      REAL(dp)  :: k(nDTEMP2,2)
      REAL(dp), INTENT(IN) :: Temp(:)
      !
      k(:,2) = RTpar2%DTEMP2(:,3)*EXP(RTpar2%DTEMP2(:,4)*(Temp(6)-InvRefTemp))
      k(:,1) = RTpar2%DTEMP2(:,1)*EXP(RTpar2%DTEMP2(:,2)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vDTemp2Compute
    
    FUNCTION vDTemp3Compute(Temp) RESULT(k)
      REAL(dp)  :: k(nDTEMP3,2)
      REAL(dp), INTENT(IN)  :: Temp(:)
      ! 
      k(:,2) = RTpar2%DTEMP3(:,4)
      k(:,1) = RTpar2%DTEMP3(:,1) * EXP(RTpar2%DTEMP3(:,2)*(RefTemp*Temp(6)-ONE)         &
      &      + RTpar2%DTEMP3(:,3)*(One-RefTemp*Temp(6) + LOG10(RefTemp*Temp(6)))) * k(:,2)
    END FUNCTION vDTemp3Compute
    
    FUNCTION vDTemp4Compute(Temp) RESULT(k)
      REAL(dp)  :: k(nDTEMP4,2)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      k(:,2) = ONE  ! OSSI
      k(:,1) = RTpar2%DTEMP4(:,1)*EXP(RTpar2%DTEMP4(:,2)                &
      &      * (Temp(1)*InvRefTemp-One) + RTpar2%DTEMP4(:,3)            &
      &      * (ONE+LOG(Temp(1)*InvRefTemp)-Temp(1)*InvRefTemp)) * k(:,2)
    END FUNCTION vDTemp4Compute
    
    FUNCTION vDTemp5Compute(Temp) RESULT(k)
      REAL(dp)  :: k(nDTEMP5,2)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      k(:,2) = ONE  ! OSSI
      k(:,1) = RTpar2%DTEMP5(:,1)*(Temp*InvRefTemp)**RTpar2%DTEMP4(:,2)  &
      &      * EXP(RTpar2%DTEMP5(:,3)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vDTemp5Compute      
    
    FUNCTION vMeskhidzeCompute(Temp) RESULT(k)
      REAL(dp)  :: k(nMeskhidze,2)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      k(:,2) = RTpar2%Meskhidze(:,4) * EXP(RTpar2%Meskhidze(:,5)  &
      &      * (Temp(6)-InvRefTemp)) * RTpar2%Meskhidze(:,7)!*(ActivityHp)**m
      k(:,1) = RTpar2%Meskhidze(:,1) * (Temp*InvRefTemp)**RTpar2%Meskhidze(:,2)  &
      &      * EXP(RTpar2%Meskhidze(:,3)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vMeskhidzeCompute
    
    FUNCTION vAspec1Compute(Temp,Hp) RESULT(kASPEC1)
      REAL(dp)  :: kASPEC1(nASPEC1)
      REAL(dp), INTENT(IN)  :: Temp(:), Hp
      REAL(dp), PARAMETER   :: x = 13.0d0
      !
      kASPEC1 = Hp*(RTpar2%ASPEC1(:,1)*EXP(RTpar2%ASPEC1(:,2) & 
      &       * (Temp(6)-InvRefTemp)))/(ONE + x*Hp)
    END FUNCTION vAspec1Compute
    
    FUNCTION vAspec2Compute(Temp,Hp) RESULT(kASPEC2)
      REAL(dp)  :: kASPEC2(nASPEC2)
      REAL(dp), INTENT(IN)  :: Temp(:), Hp
      !
      kASPEC2 = Hp**RTpar2%ASPEC2(:,2) *                       &
      &         (RTpar2%ASPEC2(:,1)*EXP(RTpar2%ASPEC2(:,3)*(Temp(6)-InvRefTemp)))
    END FUNCTION vAspec2Compute
    
    FUNCTION vAspec3Compute(Temp,Hp) RESULT(kASPEC3)
      REAL(dp)  :: kASPEC3(nASPEC3)
      REAL(dp), INTENT(IN)  :: Temp(:), Hp

      kASPEC3 = RTpar2%ASPEC3(:,1)*EXP(RTpar2%ASPEC3(:,2)*(-LOG10(Hp)))
    END FUNCTION vAspec3Compute
    
    FUNCTION vT1H2OCompute(Temp) RESULT(kT1H2O)
      REAL(dp)  :: kT1H2O(nT1H2O)
      REAL(dp), INTENT(IN)  :: Temp(:)
      
      kT1H2O = RTpar2%T1H2O(:,1)*EXP(-RTpar2%T1H2O(:,2)*Temp(6))
    END FUNCTION vT1H2OCompute


  END MODULE Rates_Mod
