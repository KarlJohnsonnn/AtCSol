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
    REAL(dp) :: LAT  = 45.0_dp
    REAL(dp) :: LONG =  0.0_dp
    REAL(dp) :: fac_exp = 1.0_dp, fac_A = 1.0_dp
    INTEGER  :: IDAT = 010619
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
      REAL(dp) :: k,       vK(neq)
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
          vMeffX = DAX_sparse( TB_sparse, Conc )
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
      !  WRITE(*,101) '  iR = ',j ,'  t = ',Time,'  Meff = ',vMeff(j), &
      !  &            '  k = ',vK(j),'  prd = ',vProd(j),'  Rate = ',Rate(j)
      !END DO
      !101 FORMAT(A,I0,A,F6.2,4(A,Es10.2))

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
      REAL(dp) :: AquaFac(nr_HOAqua)
      
      !REAL(dp) :: tHenry(nr_HENRY,2)
      REAL(dp) :: tHenry(nr_HENRY,2,ntFrac)
      INTEGER  :: j,i
      !==================================================================!
      !===============calc rates for ReactionSystem======================!
      !==================================================================!
     
      TimeRateA = MPI_WTIME()
      
      Rate = ZERO

      ! --- Compute zenith for photo reactions
      IF ( PHOTO ) THEN
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
      IF ( nr_FACTOR > 0 ) Meff = EffectiveMolecularity( Conc )
      
      ! ====== Compute the rate constant for specific reaction type
      CALL ComputeRateConstant( k, T, Time, chi, mAir, Conc, Meff )

      ! ===== correct unit of concentrations for higher order aqueous reactions
      IF ( ns_AQUA > 0 ) THEN
        LWC = pseudoLWC(Time)
        InitValKat(aH2O_ind) = aH2O*LWC
        k(iR%iHOaqua) = k(iR%iHOaqua) / (LWC*mol2part)**iR%HOaqua
      END IF

      !=== Compute mass transfer coefficient 
      IF ( nr_HENRY > 0 ) THEN
        thenry = MassTransfer( k(iR%iHENRY(:,1)), T, LWC )
        k(iR%iHENRY(:,1)) = thenry(:,1,1)
        k(iR%iHENRY(:,3)) = thenry(:,2,1)
      END IF

      ! ==== Law of mass action productories
      Prod = MassActionProducts( Conc )

      Rate = Meff * k * Prod

      !i = 0
      !DO j=1,neq
      !  WRITE(*,101) '  iR = ',j ,'  t = ',Time,'  Meff = ',Meff(j), &
      !  &            '  k = ',K(j),'  prd = ',Prod(j),'  Rate = ',Rate(j), &
      !  &            '   '//TRIM(ReactionSystem(j)%Line1)
      !END DO
      !101 FORMAT(A,I0,A,F6.2,4(A,Es10.2),A)
      !STOP 'Rates_Mod'
      
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
      REAL(dp) :: k(nr_HENRY,2,ntFrac), kin(nr_HENRY)
      REAL(dp) :: Temp(:), LWC
      ! TEMO
      REAL(dp) :: kmt(nr_HENRY,ntFrac)
      REAL(dp) :: term_diff(nr_HENRY), term_accom(nr_HENRY)
      REAL(dp) :: wetRadius(ntFrac)
      INTEGER :: i
      !
      REAL(dp) :: t1
      !
      !---------------------------------------------------------------------------
      term_diff  = henry_diff(  iR%iHENRY(:,2) )             ! diffusion term
      term_accom = henry_accom( iR%iHENRY(:,2) ) * Temp(10)  ! accom term
      !--------------------------------------------------------------------------!
      !
      ! Compute new wet radius for droplett class iFrac
      wetRadius(:) = (Pi34*LWC/Frac%Number(:))**(rTHREE)
      !
      !--  mass transfer coefficient
      kmt = dkmt  ! set minimal transfer coefficient
      DO  i = 1,nr_HENRY 
        IF (term_diff(i) /= ZERO) THEN
          kmt(i,:) = term_diff(i)*wetRadius(:)*wetRadius(:) + term_accom(i)*wetRadius(:)
        END IF
      END DO 
      kmt = ONE / kmt

      ! direaction GasSpecies-->AquaSpecies
      k(:,1,:) = milli * kmt(:,:) * LWC

      ! direaction AquaSpecies-->GasSpecies  
      DO i=1,ntFrac
        k(:,2,i) = kmt(:,i) / (kin(:) * GasConst_R * Temp(1))  ! (...) = HenryConst*GasConstant*Temperatur
      END DO

    END FUNCTION MassTransfer
    
    FUNCTION EffectiveMolecularity(Conc) RESULT(M)
      !OUT
      REAL(dp) :: M(neq)
      !IN
      REAL(dp) :: Conc(:)
      !
      M = ONE

      IF(nr_FAC_H2>0)    M(iR%iFAC_H2)   = ((mH2*mair)**fac_exp)*fac_A
      IF(nr_FAC_O2N2>0)  M(iR%iFAC_O2N2) = (((mO2*mair)*(mN2*mair))**fac_exp)*fac_A
      IF(nr_FAC_M>0)     M(iR%iFAC_M)    = (mair**fac_exp)*fac_A
      IF(nr_FAC_O2>0)    M(iR%iFAC_O2)   = ((mO2*mair)**fac_exp)*fac_A
      IF(nr_FAC_N2>0)    M(iR%iFAC_N2)   = ((mN2*mair)**fac_exp)*fac_A
      IF(nr_FAC_H2O>0)   M(iR%iFAC_H2O)  = (mH2O**fac_exp)*fac_A
      IF(nr_FAC_RO2>0)   M(iR%iFAC_RO2)  = SUM(Conc(RO2))
      IF(nr_FAC_O2O2>0)  M(iR%iFAC_O2O2) = (((mO2*mair)*(mO2*mair))**fac_exp)*fac_A
      !IF(nr_FAC_aH2O>0) M(iR%iFAC_aH2O) = aH2OmolperL*LWC*mol2part
      IF(nr_FAC_RO2aq>0) M(iR%iFAC_RO2aq) = SUM(Conc(RO2aq))

    END FUNCTION EffectiveMolecularity

    SUBROUTINE ComputeRateConstant(k,T,Time,chi,mAir,Conc,Meff)
      USE fparser

      REAL(dp), INTENT(OUT)   :: k(neq)
      REAL(dp), INTENT(IN)    :: Time, mAir, chi(:)
      REAL(dp), INTENT(IN)    :: T(:)
      REAL(dp), INTENT(IN)    :: Conc(:)
      REAL(dp), INTENT(INOUT) :: Meff(neq)

      REAL(dp) :: k_DC(nr_DCONST,2), k_T1(nr_DTEMP,2), k_T2(nr_DTEMP2,2), k_T3(nr_DTEMP3,2)
      REAL(dp) :: k_T4(nr_DTEMP4,2), k_T5(nr_DTEMP5,2), mesk(nr_Meskhidze,2)

      ! Photoabc tempo parameter
      REAL(dp), DIMENSION(nr_PHOTabc) :: ChiZabc, yChiZabc, EyChiZabc
      REAL(dp), DIMENSION(nr_PHOTmcm) :: ChiZmcm, yChiZmcm
      INTEGER :: i, j
      ! Aspec tempo parameter
      REAL(dp), PARAMETER   :: x = 13.0d0
      ! Troe tempo parameters
      REAL(dp), DIMENSION(nr_TROE)    :: k1, k2, log10_k1k2
      REAL(dp), DIMENSION(nr_TROEq)   :: k1q, k2q, k3q, log10_k1k2q
      REAL(dp), DIMENSION(nr_TROEf)   :: k1f, k2f, log10_k1k2f
      REAL(dp), DIMENSION(nr_TROEqf)  :: k1qf, k2qf, k3qf, log10_k1k2qf
      REAL(dp), DIMENSION(nr_TROExp)  :: k1xp, k2xp, log10_k1k2xp
      REAL(dp), DIMENSION(nr_TROEmcm) :: k1mcm, k2mcm, Fc, tmpTROE
      REAL(dp), PARAMETER   :: n = 0.75_dp, d = 1.27_dp
      REAL(dp) :: F, Tr300
      ! Spec tempo parameters
      REAL(dp), DIMENSION(nr_SPEC3)    :: k1s, k2s ,k3s
      REAL(dp), DIMENSION(nr_SPEC5mcm) :: k1smcm5, k2smcm5
      REAL(dp), DIMENSION(nr_SPEC6mcm) :: k1smcm6, k2smcm6
      REAL(dp), DIMENSION(nr_SPEC7mcm) :: k1smcm7, k2smcm7
      REAL(dp), DIMENSION(nr_SPEC8mcm) :: k1smcm8, k2smcm8

      k = ZERO
      
      ! *************************************************************************
      ! *** Photolytic reactions
      !
      IF (nr_PHOTAB>0) THEN
        IF ( chi(1) < PiHalf ) THEN
          k(iR%iPHOTab) = Dust * iR%PHOTab(:,1)*EXP(-iR%PHOTab(:,2)*chi(2))
        END IF
      END IF

      IF (nr_PHOTabc>0) THEN
        IF ( chi(1) < PiHalf ) THEN
          ChiZabc = chi(1) * iR%PHOTabc(:,3) 
          DO i = 1,nr_PHOTabc
            IF (ChiZabc(i) < PiHalf) THEN
              yChiZabc(i) = iR%PHOTabc(i,2) * (One - One/COS(ChiZabc(i)))
              IF ( yChiZabc(i) > mTHIRTY ) THEN
                EyChiZabc(i) = EXP(yChiZabc(i))
              ELSE
                EyChiZabc(i) = EyChiZmin   ! = 9.357d-14  
              END IF
            ELSE
              EyChiZabc(i) = EyChiZmin   ! = 9.357d-14 
            END IF
          END DO
          k(iR%iPHOTabc) = Dust * iR%PHOTabc(:,1) * EyChizabc
        END IF
      END IF

      IF (nr_PHOTMCM>0) THEN
        IF ( chi(1) < PiHalf ) THEN
          ChiZmcm  = EXP( -iR%PHOTmcm(:,3) * chi(2) )
          yChiZmcm = chi(2) ** iR%PHOTmcm(:,2)
          k(iR%iPHOTmcm) = Dust * iR%PHOTmcm(:,1) * yChiZmcm * ChiZmcm
        END IF
      END IF
      ! ************************************************************************

      ! ************************************************************************
      ! *** Constant reactions
      !
      IF (nr_CONST>0) k(iR%iCONST) = iR%CONST
      ! ************************************************************************

      ! ************************************************************************
      ! *** Temperature dependend reaction
      !
      IF (nr_TEMP1>0) THEN
        k(iR%iTEMP1) = iR%TEMP1(:,1)*EXP(-iR%TEMP1(:,2)*T(6))
      END IF
      IF (nr_TEMP2>0) THEN
        k(iR%iTEMP2) = iR%TEMP2(:,1)*T(2)*EXP(-iR%TEMP2(:,2)*T(6))
      END IF
      IF (nr_TEMP3>0) THEN
        k(iR%iTEMP3) = iR%TEMP3(:,1)*EXP(iR%TEMP3(:,2)*(T(6)-InvRefTemp))
      END IF
      IF (nr_TEMP4>0) THEN
        k(iR%iTEMP4) = iR%TEMP4(:,1)*T(1)*EXP(-iR%TEMP4(:,2)*T(6))
      END IF
      ! ************************************************************************

      ! ************************************************************************
      ! *** specieal aqua reactions
      !
      IF (nr_ASPEC1>0) THEN
        k(iR%iASPEC1) = Conc(Hp_ind)*(iR%ASPEC1(:,1)*EXP(iR%ASPEC1(:,2)       & 
        &             * (T(6)-InvRefTemp)))/(ONE + x*Conc(Hp_ind))
      END IF
      IF (nr_ASPEC2>0) THEN 
        k(iR%iASPEC2) = Conc(Hp_ind)**iR%ASPEC2(:,2)                          &
        &             * (iR%ASPEC2(:,1)*EXP(iR%ASPEC2(:,3)*(T(6)-InvRefTemp)))
      END IF
      IF (nr_ASPEC3>0) THEN 
        k(iR%iASPEC3) = iR%ASPEC3(:,1)*EXP(iR%ASPEC3(:,2)*(-LOG10(Conc(Hp_ind))))
      END IF
      ! ************************************************************************
      
      ! ************************************************************************
      ! *** dissociation reactions
      !
      IF (nr_DCONST>0) THEN
        k_DC = vDConstCompute( )
        k(iR%iDCONST(:,1)) = k_DC(:,1) ! forward reactions
        k(iR%iDCONST(:,2)) = k_DC(:,2) ! backward reactions
      END IF
      IF (nr_DTEMP>0)  THEN 
        k_T1 = vDTempCompute( T )
        k(iR%iDTEMP(:,1)) = k_T1(:,1) ! forward reactions
        k(iR%iDTEMP(:,2)) = k_T1(:,2) ! backward reactions
      END IF
      IF (nr_DTEMP2>0) THEN 
        k_T2 = vDTemp2Compute( T )
        k(iR%iDTEMP2(:,1)) = k_T2(:,1) ! forward reactions
        k(iR%iDTEMP2(:,2)) = k_T2(:,2) ! backward reactions
      END IF
      IF (nr_DTEMP3>0) THEN 
        k_T3 = vDTemp3Compute( T )
        k(iR%iDTEMP3(:,1)) = k_T3(:,1) ! forward reactions
        k(iR%iDTEMP3(:,2)) = k_T3(:,2) ! backward reactions
      END IF
      IF (nr_DTEMP4>0) THEN 
        k_T4 = vDTemp4Compute( T )
        k(iR%iDTEMP4(:,1)) = k_T4(:,1) ! forward reactions
        k(iR%iDTEMP4(:,2)) = k_T4(:,2) ! backward reactions
      END IF
      IF (nr_DTEMP5>0) THEN 
        k_T5 = vDTemp5Compute( T )
        k(iR%iDTEMP5(:,1)) = k_T5(:,1) ! forward reactions
        k(iR%iDTEMP5(:,2)) = k_T5(:,2) ! backward reactions
      END IF
      IF (nr_Meskhidze>0) THEN
        mesk = vMeskhidzeCompute( T )
        k(iR%iMeskhidze(:,1)) = mesk(:,1) ! forward reactions
        k(iR%iMeskhidze(:,2)) = mesk(:,2) ! backward reactions
      END IF
      ! ************************************************************************
     
      ! ************************************************************************
      ! *** Troe reactions
      !
      Tr300 = T(1)*r300
      IF (nr_TROE>0) THEN   
        F = 0.6_dp
        k1 = iR%TROE(:,1) * (Tr300)**(-iR%TROE(:,2))*mAir
        k2 = iR%TROE(:,3) * (Tr300)**(-iR%TROE(:,4))
        log10_k1k2  = LOG10(k1/k2)
        k(iR%iTROE) = k1/(One+k1/k2)*F**(One/(One+log10_k1k2*log10_k1k2))
      END IF

      IF (nr_TROEQ>0) THEN  
        F = 0.6_dp
        k1q = iR%TROEq(:,1)*(Tr300)**(-iR%TROEq(:,2))*mAir
        k2q = iR%TROEq(:,3)*(Tr300)**(-iR%TROEq(:,4))
        log10_k1k2q = LOG10(k1q/k2q)
        k3q = k1q/(One+k1q/k2q)*F**(One/(One+log10_k1k2q*log10_k1k2q))
        k(iR%iTROEq) = k3q/(iR%TROEq(:,5)*EXP(iR%TROEq(:,6)*T(6)))
      END IF
      
      IF (nr_TROEF>0) THEN  
        k1f = iR%TROEf(:,1)*(Tr300)**(-iR%TROEf(:,2))*mAir
        k2f = iR%TROEf(:,3)*(Tr300)**(-iR%TROEf(:,4))
        log10_k1k2f = LOG10(k1f/k2f)
        k(iR%iTROEf) = k1f/(One+k1f/k2f)*iR%TROEf(:,5)**(One/(One+log10_k1k2f*log10_k1k2f))
      END IF
      
      IF (nr_TROEQF>0) THEN 
        k1qf = iR%TROEqf(:,1) * (Tr300)**(-iR%TROEqf(:,2))*mAir
        k2qf = iR%TROEqf(:,3) * (Tr300)**(-iR%TROEqf(:,4))
        log10_k1k2qf = LOG10(k1qf/k2qf)
        k3qf = k1qf/(One+k1qf/k2qf)*iR%TROEqf(:,7)**(One/(One+log10_k1k2qf*log10_k1k2qf))
        k(iR%iTROEqf) = k3qf/(iR%TROEqf(:,5)*EXP(iR%TROEqf(:,6)*T(6)))          
      END IF
      
      IF (nr_TROEXP>0) THEN 
        k1xp = iR%TROExp(:,1)*EXP(-iR%TROExp(:,2)*T(6))*mAir
        k2xp = iR%TROExp(:,3)*EXP(-iR%TROExp(:,4)*T(6))
        log10_k1k2xp  = LOG10(k1xp/k2xp)
        k(iR%iTROExp) = k1xp/(One+k1xp/k2xp)*iR%TROExp(:,5)**(One/(One+log10_k1k2xp*log10_k1k2xp)) 
      END IF
      
      IF (nr_TROEMCM>0) THEN 
        k1mcm = iR%TROEmcm(:,1)*(T(1)*InvRefTemp)**iR%TROEmcm(:,2)*EXP(iR%TROEmcm(:,3)*T(6))*mAir
        k2mcm = iR%TROEmcm(:,4)*(T(1)*InvRefTemp)**iR%TROEmcm(:,5)*EXP(iR%TROEmcm(:,6)*T(6))
        Fc    = iR%TROEmcm(:,7)*EXP(iR%TROEmcm(:,8)*T(6))+iR%TROEmcm(:,9)*EXP(T(1)/iR%TROEmcm(:,10))
        tmpTROE = LOG10(k1mcm/k2mcm)/(n-d*LOG10(Fc))
        k(iR%iTROEmcm) = k1mcm/(One+k1mcm/k2mcm)*Fc**(One/(One+tmpTROE*tmpTROE))
      END IF
      ! ************************************************************************

      ! ************************************************************************
      ! *** Spec reactions
      !
      IF (nr_SPEC1>0) k(iR%iSPEC1) = iR%SPEC1(:,1)*(ONE+mAir*iR%SPEC1(:,2))
      IF (nr_SPEC2>0) k(iR%iSPEC2) = iR%SPEC2(:,1)*(Tr300)**iR%SPEC2(:,2)*mAir
      IF (nr_SPEC3>0) THEN 
        k1s = iR%SPEC3(:,1)*EXP(iR%SPEC3(:,2)*T(6))
        k2s = iR%SPEC3(:,3)*EXP(iR%SPEC3(:,4)*T(6))
        k3s = iR%SPEC3(:,5)*EXP(iR%SPEC3(:,6)*T(6))*mAir
        k(iR%iSPEC3) = k1s+k3s/(One+k3s/k2s)
      END IF
      IF (nr_SPEC4>0) THEN 
        k(iR%iSPEC4) = iR%SPEC4(:,1)*EXP(iR%SPEC4(:,2)*T(6))+iR%SPEC4(:,3)*EXP(iR%SPEC4(:,4)*T(6))*mAir 
      END IF
      IF (nr_SPEC1MCM>0) THEN 
        k(iR%iSPEC1mcm) = iR%SPEC1mcm(:,1)*(ONE+mAir*iR%SPEC1mcm(:,2)*Tr300/iR%SPEC1mcm(:,3))
      END IF
      IF (nr_SPEC2MCM>0) THEN 
        k(iR%iSPEC2mcm) = iR%SPEC2mcm(:,1)*(Tr300)**iR%SPEC2mcm(:,2)*EXP(iR%SPEC2mcm(:,3)*T(6))
      END IF
      IF (nr_SPEC3MCM>0) THEN 
        k(iR%iSPEC3mcm) = iR%SPEC3mcm(:,1)*(ONE+mAir/iR%SPEC3mcm(:,2))
      END IF
      IF (nr_SPEC4MCM>0) THEN 
        k(iR%iSPEC4mcm) = iR%SPEC4mcm(:,1)*(ONE+iR%SPEC4mcm(:,2) &
        &               * EXP(iR%SPEC4mcm(:,3)*T(6))*H2O)*EXP(iR%SPEC4mcm(:,4)*T(6))
      END IF
      IF (nr_SPEC5MCM>0) THEN 
        F = 0.21_dp
        k1smcm5 = iR%SPEC5mcm(:,1)*mAir*F*EXP(iR%SPEC5mcm(:,2)*T(6))
        k2smcm5 = iR%SPEC5mcm(:,3)*mAir*F*EXP(iR%SPEC5mcm(:,4)*T(6))
        k(iR%iSPEC5mcm) = k1smcm5*(One-k2smcm5)
      END IF
      IF (nr_SPEC6MCM>0) THEN 
        k1smcm6 = iR%SPEC6mcm(:,1)*EXP(iR%SPEC6mcm(:,2)*T(6))
        k2smcm6 = iR%SPEC6mcm(:,3)*EXP(iR%SPEC6mcm(:,4)*T(6))
        k(iR%iSPEC6mcm) = k1smcm6*(ONE-k2smcm6)
      END IF
      IF (nr_SPEC7MCM>0) THEN 
        k1smcm7 = iR%SPEC7mcm(:,1)*EXP(iR%SPEC7mcm(:,2)*T(6))
        k2smcm7 = iR%SPEC7mcm(:,3)*EXP(iR%SPEC7mcm(:,4)*T(6))
        k(iR%iSPEC7mcm) = k1smcm7*(iR%SPEC7mcm(:,5)-iR%SPEC7mcm(:,6)/(One+k2smcm7))
      END IF
      IF (nr_SPEC8MCM>0) THEN 
        F = 0.21_dp
        k1smcm8 = iR%SPEC8mcm(:,1)*mAir*F*EXP(iR%SPEC8mcm(:,2)*T(6))
        k2smcm8 = iR%SPEC8mcm(:,3)*mAir*F*EXP(iR%SPEC8mcm(:,4)*T(6))
        k(iR%iSPEC8mcm) = k1smcm8/(One+k2smcm8)*T(6)
      END IF
      ! ************************************************************************

      IF (nr_T1H2O>0) k(iR%iT1H2O) = iR%T1H2O(:,1)*EXP(-iR%T1H2O(:,2)*T(6))
      IF (nr_S4H2O>0) THEN 
        k(iR%iS4H2O) = iR%S4H2O(:,1)*EXP(iR%S4H2O(:,2)*T(6))      & 
        &            + iR%S4H2O(:,3)*EXP(iR%S4H2O(:,4)*T(6))*mAir
      END IF
      
      ! *** KKP photolytic reactions
      IF (nr_PHOTOkpp>0)  THEN 
        k(iR%iPHOTOkpp)  = vPhotokppCompute ( Time, Chi(3) )
      END IF
      IF (nr_PHOTO2kpp>0) THEN 
        k(iR%iPHOTO2kpp) = vPhoto2kppCompute( Time, Chi(3) )
      END IF
      IF (nr_PHOTO3kpp>0) THEN 
        k(iR%iPHOTO3kpp) = vPhoto3kppCompute( Time, Chi(3) )
      END IF
      
      ! special reactions
      IF (nr_special>0) THEN 
        DO i = 1,nr_SPECIAL
          j = iR%iSPECIAL(i)
          IF (ReactionSystem(j)%Special%Temp) THEN
            k(j) = evalf(i,[ Conc(ReactionSystem(j)%Special%iVariables),T(1)])
          ELSE
            k(j) = evalf(i,Conc(ReactionSystem(j)%Special%iVariables))
          END IF
        END DO
      END IF

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
      INTEGER  :: IIYEAR, IYEAR, IMTH, IDAY, IIY, NYEARS, LEAP, NOLEAP
      REAL(dp) :: YREF,YR
      !   
      INTEGER  :: I, IJ, JD, IJD, IN
      REAL(dp) :: D, RML, W, WR, EC, EPSI, PEPSI, YT, CW, SW, SSW  & 
      &         , EYT, FEQT1, FEQT2, FEQT3, FEQT4, FEQT5, FEQT6 &
      &         , FEQT7, FEQT, EQT
      !         
      REAL(dp) :: REQT, RA, RRA, TAB, RDECL, DECL, ZPT, CSZ, ZR    &
      &         , CAZ, RAZ, AZIMUTH
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
      TempArr(1) = Temperature                 ! T
      DO i=2,5
        TempArr(i) = TempArr(i-1)*Temperature  ! T^2 ... T^5
      END DO
      TempArr(6)  = ONE / Temperature          ! 1/T
      TempArr(7)  = ONE / TempArr(2)           ! 1/T^2
      TempArr(8)  = LOG(Temperature)           ! ln(T)
      TempArr(9)  = SQRT(Temperature)          ! sqrt(T)
      TempArr(10) = ONE / TempArr(9)           ! 1/sqrt(T)
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
    
    FUNCTION vPhotokppCompute(Time,Sun) RESULT(kPHOTO)
      REAL(dp) :: kPHOTO(nr_PHOTOkpp)
      REAL(dp) :: Time, Sun
      !
      kPHOTO = iR%PHOTOkpp(:) * Sun
    END FUNCTION vPhotokppCompute
    
    FUNCTION vPhoto2kppCompute(Time,Sun) RESULT(kPHOTO2)
      REAL(dp) :: kPHOTO2(nr_PHOTO2kpp)
      REAL(dp) :: Time, Sun
      !
      kPHOTO2 = iR%PHOTO2kpp(:) * Sun*Sun
    END FUNCTION vPhoto2kppCompute
    
    FUNCTION vPhoto3kppCompute(Time,Sun) RESULT(kPHOTO3)
      REAL(dp) :: kPHOTO3(nr_PHOTO3kpp)
      REAL(dp) :: Time, Sun
      !
      kPHOTO3 = iR%PHOTO3kpp(:) * Sun*Sun*Sun
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
      REAL(dp)  :: k(nr_DCONST,2)
      !
      k(:,2) = iR%DCONST(:,2)
      k(:,1) = iR%DCONST(:,1) * k(:,2)
    END FUNCTION vDConstCompute
    
    FUNCTION vDTempCompute(Temp) RESULT(k)
      REAL(dp)  :: k(nr_DTEMP,2)
      REAL(dp), INTENT(IN) :: Temp(:)
      !
      k(:,2) = iR%DTEMP(:,3)
      k(:,1) = iR%DTEMP(:,1)*EXP(iR%DTEMP(:,2)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vDTempCompute
    
    FUNCTION vDTemp2Compute(Temp)   RESULT(k)
      REAL(dp)  :: k(nr_DTEMP2,2)
      REAL(dp), INTENT(IN) :: Temp(:)
      !
      k(:,2) = iR%DTEMP2(:,3)*EXP(iR%DTEMP2(:,4)*(Temp(6)-InvRefTemp))
      k(:,1) = iR%DTEMP2(:,1)*EXP(iR%DTEMP2(:,2)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vDTemp2Compute
    
    FUNCTION vDTemp3Compute(Temp) RESULT(k)
      REAL(dp)  :: k(nr_DTEMP3,2)
      REAL(dp), INTENT(IN)  :: Temp(:)
      ! 
      k(:,2) = iR%DTEMP3(:,4)
      k(:,1) = iR%DTEMP3(:,1) * EXP(iR%DTEMP3(:,2)*(RefTemp*Temp(6)-ONE)         &
      &      + iR%DTEMP3(:,3)*(One-RefTemp*Temp(6) + LOG10(RefTemp*Temp(6)))) * k(:,2)
    END FUNCTION vDTemp3Compute
    
    FUNCTION vDTemp4Compute(Temp) RESULT(k)
      REAL(dp)  :: k(nr_DTEMP4,2)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      k(:,2) = ONE  ! OSSI
      k(:,1) = iR%DTEMP4(:,1)*EXP(iR%DTEMP4(:,2)                &
      &      * (Temp(1)*InvRefTemp-One) + iR%DTEMP4(:,3)            &
      &      * (ONE+LOG(Temp(1)*InvRefTemp)-Temp(1)*InvRefTemp)) * k(:,2)
    END FUNCTION vDTemp4Compute
    
    FUNCTION vDTemp5Compute(Temp) RESULT(k)
      REAL(dp)  :: k(nr_DTEMP5,2)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      k(:,2) = ONE  ! OSSI
      k(:,1) = iR%DTEMP5(:,1)*(Temp*InvRefTemp)**iR%DTEMP4(:,2)  &
      &      * EXP(iR%DTEMP5(:,3)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vDTemp5Compute      
    
    FUNCTION vMeskhidzeCompute(Temp) RESULT(k)
      REAL(dp)  :: k(nr_Meskhidze,2)
      REAL(dp), INTENT(IN)  :: Temp(:)
      !
      k(:,2) = iR%Meskhidze(:,4) * EXP(iR%Meskhidze(:,5)  &
      &      * (Temp(6)-InvRefTemp)) * iR%Meskhidze(:,7)!*(ActivityHp)**m
      k(:,1) = iR%Meskhidze(:,1) * (Temp*InvRefTemp)**iR%Meskhidze(:,2)  &
      &      * EXP(iR%Meskhidze(:,3)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vMeskhidzeCompute
    
  END MODULE Rates_Mod
