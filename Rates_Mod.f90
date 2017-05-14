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
    REAL(RealKind) :: LAT=45.0D0
    REAL(RealKind) :: LONG=0.0D0
    REAL(RealKind) :: fac_exp=1.0d0, fac_A=1.0d0
    Integer :: iDat=010619
    !
    ! some factors for calculating Troe press dep. reactions
    REAL(RealKind) :: rFacEq      ! factor nessesary for equilibrium constant
    REAL(RealKind) :: cTroe, n1Troe
    REAL(RealKind), PARAMETER :: dTroe =  0.14d0

    REAL(RealKind) :: log10_Pr, log10_Fcent, Pr
    INTEGER :: globi
    !
    CONTAINS
    !
    !
    !======================================================================!
    !      Calculate the Rates for current concentraions and time
    !======================================================================!
    SUBROUTINE ReactionRatesAndDerivative(Time,Y_in,Rate,DRatedT)
      !--------------------------------------------------------------------!
      ! Input: 
      !   - Time
      !   - y_conc
!     REAL(RealKind), INTENT(IN) :: Time
      REAL(RealKind) :: Time
      REAL(RealKind), INTENT(IN) :: Y_in(nDIM)
      !--------------------------------------------------------------------!
      ! Output:
      !   - Rate vector
      REAL(RealKind) :: Rate(neq)
      REAL(RealKind) :: DRatedT(neq)
      !--------------------------------------------------------------------!
      ! Temporary variables:
      REAL(RealKind) :: y_conc(nDIM)
      REAL(RealKind) :: chi(2), LWC
      REAL(RealKind) :: T(10)
      REAL(RealKind) :: k,       vK(neq)
      REAL(RealKind) :: Prod,    vProd(neq) 
      REAL(RealKind) :: AquaFac, vAquaFac(nHOAqua)
      REAL(RealKind) :: DkdT,    dvK_dT(neq)   
      REAL(RealKind) :: Meff,    vMeff(neq)
      INTEGER :: iReac, ii, j
      ! 
      ! TempEq
      REAL(RealKind) :: Temp_in
      !
      REAL(RealKind) :: tmpK(neq), dtmpK(neq)
      
      ! temp arrays for chemkin input (temp depended)
      REAL(RealKind) :: vMeffX(RTind%nTBodyExtra)
      REAL(RealKind) :: vFTL(neq),      DF_PDdT(neq)
      REAL(RealKind) :: Dk0dT(neq),     DkinfdT(neq),      vdFTL_dT(neq)
      REAL(RealKind) :: k0(RTind%nLow), kinf(RTind%nHigh), k0M(RTind%nLow)
      REAL(RealKind) :: vrkinfpk0M(RTind%nLow)
      REAL(RealKind) :: vrKeq(RTind%nEqui), DeRdT(RTind%nEqui)
      REAL(RealKind) :: rRcT

      REAL(RealKind) :: tHenry(nHENRY,2)
      !==================================================================!
      !===============calc rates for ReactionSystem======================!
      !==================================================================!
     
      TimeRateA = MPI_WTIME()
      
      Rate    = ZERO
      DRatedT = ZERO

      !print*, 'neq=',neq
      ! --- Compute zenith for photo reactions
      IF (nreakgphoto>0) chi(1) = Zenith(Time); chi(2) = ONE/COS(Chi(1))
     
      ! --- Compute the liquid water content 
      IF (nreakaqua>0)   LWC = pseudoLWC(Time)
      
      IF ( TempEq ) THEN
        
        Temp_in  = Y_in(nDIM)
        !IF (.NOT. Temp_in > ZERO) STOP ' Temperatur not > 0 '

        !--- Concentration
        y_conc(1:nspc) = Y_in(1:nspc)
        
        !--- Temperature
        y_conc(nDIM)   = Temp_in 


        CALL UpdateTempArray ( T   , Temp_in )
        CALL GibbsFreeEnergie( GFE , T )
        CALL CalcDeltaGibbs  ( DelGFE )
        rFacEq = mega * R * T(1) * rPatm  ! in [cm3/mol]
        !
        CALL CalcDiffGibbsFreeEnergie( DGFEdT , T )
        CALL CalcDiffDeltaGibbs( DDelGFEdT )

        !do ii=1,nspc
        !  write(985,*) 'DBg ::conc = ',y_conc(SCperm(ii))
        !  write(986,*) 'DBg ::GFE = ',GFE(SCperm(ii))
        !end do
        !do ii=1,neq
        !  write(984,*) 'DBg ::DelGFE = ',DelGFE(ii)
        !end do

      ELSE

        T(1)   = 280.0d0    ! = 280 [K]
        CALL UpdateTempArray ( T   , T(1) )
        y_conc = Y_in

      END IF
      
      !*************************************************************
      ! --- compute rate of all reactions (gas,henry,aqua,diss) ---
      !*************************************************************
      IF ( Vectorized ) THEN
        IF ( Chemkin ) THEN 
          ! ====== Computing effective molecularity 
          !
          ! initializing vector components
          vMeff = ONE; DF_PDdT = Zero; 

          ! calculate effective molecularity of the reactions
          IF ( RTind%nTBody>0 ) THEN
            vMeff( RTind%iTBody ) = SUM( y_conc(1:nspc) )

            IF ( RTind%nTBodyExtra>0) THEN
              CALL DAX_sparse(vMeffX,TB_sparse,y_conc(1:nspc))
              vMeff(RTind%iTBodyExtra) = vMeff(RTind%iTBodyExtra) - vMeffX
            END IF

          END IF

          ! ====== Compute the rate constant for specific reaction type
          !
          !print*, ''
          !print*, 'RTind%nLind, RTind%nHigh,RTind%nLow,RTind%nxrev,RTind%nArr, RTind%nEqui'
          !print*,RTind%nLind, RTind%nHigh,RTind%nLow,RTind%nxrev,RTind%nArr, RTind%nEqui
          !print*, 'SIZE(RTpar%A), SIZE(RTpar%A0),SIZE(RTpar%Ainf),SIZE(RTpar%T1)'
          !print*,SIZE(RTpar%A), SIZE(RTpar%A0),SIZE(RTpar%Ainf),SIZE(RTpar%T1)
          !stop

          ! normal Arrhenius
          rRcT = rRcal*T(6)
          !print*, ' ind arr =',RTind%iArr
          !print*, ' par arr =',RTpar%A(1),RTpar%b(1),RTpar%E(1)
          IF ( RTind%nArr>0 )  THEN
            tmpK(RTind%iArr)  = RTpar%A*EXP(RTpar%b*T(8) - RTpar%E*rRcT)
            dtmpK(RTind%iArr) = (RTpar%b + RTpar%E*rRcT)*T(6)
            !print*, ' k 59 arr =',tmpk(67)
          END IF
          !print *, 'k arr = ',tmpK(RTind%iArr)
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

            !print*, ' k par::  = ',RTpar%A0(1),RTpar%b0(1),RTpar%E0(1)
            !print*, ' k par::  = ',RTpar%Ainf(1),RTpar%binf(1),RTpar%Einf(1)
            !print*, ' k par::  = ',kinf(1), k0(1)
            !stop
            
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
            !print*, ' troefac = ',vFTL(17)
          END IF

          ! Backward reaction konstant for equilibrium reactions
          !print*, ' ind arr =',RTind%iEqui
          !print*, ' par arr =',DelGFE(RTind%iEqui), sumBAT(RTind%iEqui)
          !print*,' R67  = ',ReactionSystem(67)%Line1
          !print*, ' ind arr =',RTind%iArr
          !print*, ' par arr =',RTpar%A(67),RTpar%b(67),RTpar%E(67)
          !print*, ' k 59 arr =',tmpk(67),tmpk(68)


          vK     = vFTL * tmpK
          dvK_dT = (DF_PDdT + dtmpK)*T(6)
          
          !Rate = vMeff * vk
          DRatedT = dvK_dT


        ELSE !if tropos data syntax

          ! ====== Computing effective molecularity 
          vMeff = ONE
          IF ( nFACTOR > 0 ) CALL vEffectiveMolecularity( vMeff, y_conc(1:nspc), LWC )
          
          ! ====== Compute the rate constant for specific reaction type ===
          CALL vComputeRateConstant( vk, dvK_dT, T, Time, chi, mAir, y_conc, vMeff )

          ! ===== correct unit of concentrations higher order aqueous reactions
          IF ( ntAqua > 0 ) THEN
            InitValKat(aH2O_ind) = aH2O*LWC
            vAquaFac = ONE / (LWC*mol2part)**RTpar2%HOaqua
            vk(RTind2%iHOaqua) = vk(RTind2%iHOaqua) * vAquaFac
          END IF

          !=== Compute mass transfer coefficient 
          IF ( nHENRY > 0 ) THEN
            CALL vMassTransfer( thenry, vk(RTind2%iHENRY(:,1)), T, LWC )
            vk(RTind2%iHENRY(:,1)) = thenry(:,1)
            vk(RTind2%iHENRY(:,3)) = thenry(:,2)
          END IF

          IF (TempEq) DRatedT = dvK_dT

        END IF

        ! ==== Law of mass action productories
        !
        vProd = MassActionProducts(y_conc(1:nspc))
        !
        Rate  = vMeff * vk * vProd
        
        !do j=1,neq; WRITE(987,*) 'DBg:: i, k , prd, Rate =', j, vK(j), vProd(j), Rate(j); enddo
        !stop ' rates mod '

      ELSE

        LOOP_OVER_ALL_REACTIONS: DO ii = 1 , loc_rateCnt
          
          ! if more than one processor is used split rate calculation
          !IF ( MPI_np>1 .AND. ParOrdering>=0 )  THEN
          !  iReac = loc_RatePtr(ii)
          !ELSE
          iReac = ii
          !END IF
         
          ! ====== Computing effective molecularity 
          CALL EffectiveMolecularity( Meff , y_conc(1:nspc) , iReac , LWC )
          
          ! ====== Compute the rate constant for specific reaction type ===
          CALL ComputeRateConstant( k , DkdT , T , Time , chi , mAir , iReac , y_conc , Meff )
  
          AQUATIC_REACTIONS: IF ( ntAqua > 0 ) THEN
  
            IF ( ReactionSystem(iReac)%Type == 'DISS' .OR. &
            &    ReactionSystem(iReac)%Type == 'AQUA' ) THEN
              ! ===== correct unit of concentrations
              AquaFac = ONE / ( (LWC*mol2part)**ReactionSystem(iReac)%SumAqCoef )
              k       = k * AquaFac
            END IF
  
            IF ( ReactionSystem(iReac)%Type == 'HENRY' ) THEN
              !=== Compute Henry mass transfer coefficient
              CALL MassTransfer( k , T(1) , Time , iReac , LWC )
            END IF
  
          END IF AQUATIC_REACTIONS

          ! === Calculate the product of concentrations 
          Prod = ONE
          IF (ReactionSystem(iReac)%nInActEd/=0) THEN
            SELECT CASE (ReactionSystem(iReac)%InActEductSpc(1))
              CASE ('[H2O]')
                Prod = Prod*H2O
              CASE ('[N2]')
                Prod = Prod*N2
              CASE ('[O2]')
                Prod = Prod*O2
              CASE ('[aH2O]')
                Prod = Prod*aH2O*LWC
            END SELECT
          END IF
          !
          ! ====================== Mass action products ===================
          DO j = A%RowPtr(iReac) , A%RowPtr(iReac+1)-1
            IF ( A%Val(j) == ONE ) THEN
              Prod = Prod * y_conc(A%ColInd(j))
            ELSE IF ( A%Val(j) == TWO ) THEN
              Prod = Prod * y_conc(A%ColInd(j))*ABS(y_conc(A%ColInd(j)))
            ELSE
              Prod = Prod * y_conc(A%ColInd(j))**A%Val(j)
            END IF
          END DO
      
          Rate(iReac) = Meff * k * Prod
          IF (TempEq) DRatedT(iReac) = DkdT
  
          !WRITE(987,*) 'DBg:: i, k , prd, Rate =', iReac, k, Prod, Rate(iReac), Meff
        END DO LOOP_OVER_ALL_REACTIONS
      END IF

      TimeRates = TimeRates + MPI_WTIME() - TimeRateA

      !stop 'ratesmod'
      !
      
      ! gather the values of the other processes
      !CALL GatherAllPartitions(Rate,MyParties)
      !IF (TempEq) CALL GatherAllPartitions(DRatedT,MyParties)
    END SUBROUTINE ReactionRatesAndDerivative


    FUNCTION MassActionProducts(Conc) RESULT(vProd)
      REAL(RealKind) :: vProd(neq)
      REAL(RealKind), INTENT(IN)  :: Conc(nspc)
      INTEGER :: i

      vProd = ONE
      !
      ! stoechometric coefficients equal 1
      DO i=1,nFirst_order
        vProd(first_order(i,1)) = vProd(first_order(i,1)) &
         &                      * Conc(first_order(i,2))
      END DO
      !
      ! stoechometric coefficients equal 2
      DO i=1,nSecond_order
        vProd(second_order(i,1)) = vProd(second_order(i,1))   &
        &                        * Conc(second_order(i,2))    &
        &                        * ABS(Conc(second_order(i,2)))
      END DO
      !
      ! stoechometric coefficients not equal 1 or 2
      DO i=1,nHigher_order
        vProd(higher_order(i,1)) =  vProd(higher_order(i,1)) &
        &                        *  Conc(higher_order(i,2))  &
        &                        ** ahigher_order(i)
      END DO
      !
      ! if there are passive (katalytic) species e.g. [N2], [O2] or [aH2O]
      DO i=1,nfirst_orderKAT
        vProd(first_orderKAT(i,1)) = vProd(first_orderKAT(i,1))    & 
        &                          * InitValKat(first_orderKAT(i,2)) 
      END DO
      

      !stop ' massen prod'
    END FUNCTION MassActionProducts


    !=====================================================================!
    ! === Converts the mass for Henry reactions Gas->Aqua , Aqua->Gas
    !=====================================================================!
    SUBROUTINE MassTransfer(k,Temp,Time,iReac,LWC)
      REAL(RealKind) :: k
      REAL(RealKind) :: Temp, Time, LWC
      INTEGER :: iReac
      !
      REAL(RealKind) :: term_diff
      REAL(RealKind) :: term_accom
      REAL(RealKind) :: kmt
      !
      !
      !---------------------------------------------------------------------------
      term_diff   = henry_diff(ReactionSystem(iReac)%HenrySpc)               ! diffusion term
      term_accom  = henry_accom(ReactionSystem(iReac)%HenrySpc) / SQRT(Temp)  ! accom term
      !--------------------------------------------------------------------------!
      !
      ! Compute new wet radius for droplett class iFrac
      SPEK(1)%wetRadius = (Pi34*LWC/SPEK(1)%Number)**(rTHREE)
      !
      !--  mass transfer coefficient
      IF ( term_diff /= ZERO )  THEN   
        kmt = ONE/(term_diff*SPEK(1)%wetRadius*SPEK(1)%wetRadius &
        &                         + term_accom*SPEK(1)%wetRadius  )
      ELSE
        kmt = dkmt
      END IF
      !
      !print*, 'Debug::   vor  k=',k
      !print*, 'Debug::      kmt=',kmt
      ! direaction GasSpecies-->AquaSpecies
      IF (ReactionSystem(iReac)%direction=='GA') THEN  
        k = milli * kmt * LWC ! orginal
        !
      ! direaction AquaSpecies-->GasSpecies  
      ELSE 
        k = kmt / ( k * GasConst_R * Temp)   !()=HenryConst*GasConstant*Temperatur
      END IF
      !print*, 'Debug:: gasconst=',GasConst_R
      !print*, 'Debug::     henry_diff=',henry_diff(ReactionSystem(iReac)%HenrySpc)    
      !print*, 'Debug::     henry_accom=',henry_accom(ReactionSystem(iReac)%HenrySpc)
    END SUBROUTINE MassTransfer
    SUBROUTINE vMassTransfer(k,kin,Temp,LWC)
      REAL(RealKind) :: k(nHENRY,2), kin(nHENRY)
      REAL(RealKind) :: Temp(:), LWC
      ! TEMO
      REAL(RealKind), DIMENSION(nHENRY) :: term_diff, term_accom, kmt, kTemp
      REAL(RealKind) :: wetRadius
      INTEGER :: i
      !
      !
      !---------------------------------------------------------------------------
      term_diff  = henry_diff(  RTind2%iHENRY(:,2) )               ! diffusion term
      term_accom = henry_accom( RTind2%iHENRY(:,2) ) * Temp(10)  ! accom term
      !--------------------------------------------------------------------------!
      !
      ! Compute new wet radius for droplett class iFrac
      wetRadius = (Pi34*LWC/SPEK(1)%Number)**(rTHREE)
      !
      !--  mass transfer coefficient
      kmt = dkmt  ! set minimal transfer coefficient
      FORALL ( i = 1:nHENRY , term_diff(i) /= ZERO )
        kmt(i) = ONE / ( term_diff(i)*wetRadius*wetRadius + term_accom(i)*wetRadius )
      END FORALL

      ! direaction GasSpecies-->AquaSpecies
      k(:,1) = milli * kmt * LWC

      ! direaction AquaSpecies-->GasSpecies  
      k(:,2) = kmt / ( kin * GasConst_R * Temp(1))  ! (...) = HenryConst*GasConstant*Temperatur

    END SUBROUTINE vMassTransfer
    !
    !=======================================================================!
    ! ===  Select the Type of the Constant and calculate the value
    !=======================================================================!
    SUBROUTINE ComputeRateConstant(k,DkdT,T,Time,chi,mAir,iReac,y_conc,Meff)
      REAL(RealKind), INTENT(INOUT) :: k, DkdT
      REAL(RealKind), INTENT(IN) :: Time, mAir, chi(:)
      REAL(RealKind), INTENT(IN) :: T(10)
      REAL(RealKind), INTENT(IN) :: y_conc(:)
      INTEGER,        INTENT(IN) :: iReac
      REAL(RealKind), INTENT(INOUT), OPTIONAL :: Meff
      REAL(RealKind) :: EqRate,BaRate,FoRate
      !
      REAL(RealKind) :: k0, kinf
      ! calc reaction constant
      ! Skip photochemical reactions at night
      BaRate = ZERO
      DkdT   = ZERO
      !
      SELECT CASE (ReactionSystem(iReac)%TypeConstant)
        CASE ('PHOTABC')
          CALL PhoABCCompute(k,ReactionSystem(iReac)%Constants,Time,chi)  
        CASE ('PHOTMCM')
          CALL PhoMCMCompute(k,ReactionSystem(iReac)%Constants,Time,chi)  
        CASE ('PHOTAB')
          CALL PhoABCompute(k,ReactionSystem(iReac)%Constants,Time,chi)  
        CASE ('CONST')
          CALL ConstCompute(k,ReactionSystem(iReac)%Constants)
        CASE ('TEMP1')
          CALL Temp1Compute(k,ReactionSystem(iReac)%Constants,T(1))
        CASE ('TEMP2')
          CALL Temp2Compute(k,ReactionSystem(iReac)%Constants,T(1))
        CASE ('TEMP3')                    ! Henry Rate calculation
          CALL Temp3Compute(k,ReactionSystem(iReac)%Constants,T(1)) ! compute HenryConst
        CASE ('TEMP4')
          CALL Temp4Compute(k,ReactionSystem(iReac)%Constants,T(1))
        CASE ('TEMPX')
          ! new chemkin based routines
          CALL TempXCompute(k,EqRate,FoRate,ReactionSystem(iReac)%Constants,T,iReac)
          CALL DiffTempXCompute(DkdT,ReactionSystem(iReac)%Constants,T,iReac)
        CASE ('PRESSX')
          ! new chemkin based routines
          CALL PressXCompute(k,k0,kinf,iReac,T,Meff)
          CALL DiffPressXCompute(DkdT,k0,kinf,iReac,T,Meff)
          Meff = ONE
        CASE ('ASPEC1')
          CALL Aspec1Compute(k,ReactionSystem(iReac)%Constants,T(1),y_conc)
        CASE ('ASPEC2')
          CALL Aspec2Compute(k,ReactionSystem(iReac)%Constants,T(1),y_conc)
        CASE ('ASPEC3')
          CALL Aspec3Compute(k,ReactionSystem(iReac)%Constants,T(1),y_conc)
        CASE ('DCONST')
          CALL DConstCompute(EqRate,BaRate,ReactionSystem(iReac)%Constants)  
          CALL ReactionDirection(k,EqRate,BaRate,iReac)
        CASE ('DTEMP')
          CALL DTempCompute(EqRate,BaRate,ReactionSystem(iReac)%Constants,T(1))     
          CALL ReactionDirection(k,EqRate,BaRate,iReac)
        CASE ('DTEMP2')
          CALL DTemp2Compute(EqRate,BaRate,ReactionSystem(iReac)%Constants,T(1))     
          CALL ReactionDirection(k,EqRate,BaRate,iReac)
        CASE ('DTEMP3')
          CALL DTemp3Compute(EqRate,BaRate,ReactionSystem(iReac)%Constants,T(1))     
          CALL ReactionDirection(k,EqRate,BaRate,iReac)
        CASE ('DTEMP4')
          CALL DTemp4Compute(EqRate,BaRate,ReactionSystem(iReac)%Constants,T(1))     
          CALL ReactionDirection(k,EqRate,BaRate,iReac)
        CASE ('DTEMP5')
          CALL DTemp5Compute(EqRate,BaRate,ReactionSystem(iReac)%Constants,T(1))     
          CALL ReactionDirection(k,EqRate,BaRate,iReac)
        CASE ('TROE')
          CALL TroeCompute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('TROEQ')
          CALL TroeEqCompute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('TROEF')
          CALL TroeFCompute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('TROEQF')
          CALL TroeEqfCompute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('TROEXP')
          CALL TroeXPCompute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('TROEMCM')
          CALL TroeMCMCompute(k,ReactionSystem(iReac)%Constants,T(1),mAir)!
        CASE ('SPEC1')
          CALL Spec1Compute(k,ReactionSystem(iReac)%Constants,mAir)
        CASE ('SPEC2')
          CALL Spec2Compute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('SPEC3')
          CALL Spec3Compute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('SPEC4')
          CALL Spec4Compute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('SPEC1MCM')
          CALL SPEC1MCMCompute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('SPEC2MCM')
          CALL SPEC2MCMCompute(k,ReactionSystem(iReac)%Constants,T(1))
        CASE ('SPEC3MCM')
          CALL SPEC3MCMCompute(k,ReactionSystem(iReac)%Constants,mAir)
        CASE ('SPEC4MCM')
          CALL Spec4MCMCompute(k,ReactionSystem(iReac)%Constants,T(1))
        CASE ('SPEC5MCM')
          CALL Spec5MCMCompute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('SPEC6MCM')
          CALL Spec6MCMCompute(k,ReactionSystem(iReac)%Constants,T(1))
        CASE ('SPEC7MCM')
          CALL Spec7MCMCompute(k,ReactionSystem(iReac)%Constants,T(1))
        CASE ('SPEC8MCM')
          CALL Spec8MCMCompute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('T1H2O')
          CALL T1H2OCompute(k,ReactionSystem(iReac)%Constants,T(1))
        CASE ('S4H2O')
          CALL S4H2OCompute(k,ReactionSystem(iReac)%Constants,T(1),mAir)
        CASE ('EQUI')
          CALL EquiCompute(k,ReactionSystem(iReac)%Constants)
        CASE ('PHOTO')    ! kpp
          CALL PHOTO(k,ReactionSystem(iReac)%Constants,Time)
        CASE ('PHOTO2')   ! kpp
          CALL PHOTO2(k,ReactionSystem(iReac)%Constants,Time)
        CASE ('PHOTO3')   ! kpp
          CALL PHOTO3(k,ReactionSystem(iReac)%Constants,Time)
        CASE DEFAULT
          WRITE(*,*) 'Reaction  ',iReac,':   ',TRIM(ReactionSystem(iReac)%Line1)
          WRITE(*,*) 'Unknown ReactionType: ',ADJUSTL(ReactionSystem(iReac)%TypeConstant)
          CALL FinishMPI()
          STOP 'STOP Rates_Mod'
      END SELECT
    END SUBROUTINE ComputeRateConstant
    SUBROUTINE vComputeRateConstant(k,DkdT,T,Time,chi,mAir,y_conc,Meff)
      REAL(RealKind), INTENT(OUT) :: k(neq), DkdT(neq)
      REAL(RealKind), INTENT(IN) :: Time, mAir, chi(2)
      REAL(RealKind), INTENT(IN) :: T(10)
      REAL(RealKind), INTENT(IN) :: y_conc(:)
      REAL(RealKind), INTENT(INOUT) :: Meff(neq)

      REAL(RealKind) :: k_DC(nDCONST,2), k_T1(nDTEMP,2), k_T2(nDTEMP2,2), k_T3(nDTEMP3,2)
      REAL(RealKind) :: k_T4(nDTEMP4,2), k_T5(nDTEMP5,2), mesk(nMeskhidze,2)

      
      DkdT   = ZERO

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
      IF (nASPEC1>0) k(RTind2%iASPEC1) = vAspec1Compute( T, y_conc(Hp_ind) )
      IF (nASPEC2>0) k(RTind2%iASPEC2) = vAspec2Compute( T, y_conc(Hp_ind) )
      IF (nASPEC3>0) k(RTind2%iASPEC3) = vAspec3Compute( T, y_conc(Hp_ind) )
      
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

      IF (nT1H2O>0) k(RTind2%iT1H2O) =  vT1H2OCompute( T )
      IF (nS4H2O>0) k(RTind2%iS4H2O) = vS4H2OCompute ( T, mAir)
      
      ! *** KKP photolytic reactions
      !IF (nEQUI>0)    CALL vEquiCompute(k,ReactionSystem(iReac)%Constants)
      !IF (nPHOTO>0)   CALL vPHOTO(k,ReactionSystem(iReac)%Constants,Time)
      !IF (nPHOTO2>0)  CALL vPHOTO2(k,ReactionSystem(iReac)%Constants,Time)
      !IF (nPHOTO3>0)  CALL vPHOTO3(k,ReactionSystem(iReac)%Constants,Time)

      ! *** chemkin stuff
      !IF (nTEMPX>0) THEN
      !  ! new chemkin based routines
      !  CALL TempXCompute(k,EqRate,FoRate,ReactionSystem(iReac)%Constants,T,iReac)
      !  CALL DiffTempXCompute(DkdT,ReactionSystem(iReac)%Constants,T,iReac)
      !END IF
      !IF (nPRESSX>0) THEN
      !  ! new chemkin based routines
      !  CALL PressXCompute(k,k0,kinf,iReac,T,Meff)
      !  CALL DiffPressXCompute(DkdT,k0,kinf,iReac,T,Meff)
      !  Meff = ONE
      !END IF
      !
    END SUBROUTINE vComputeRateConstant
    !
    !
    SUBROUTINE ReactionDirection(ReacConst,EquiRate,BackRate,iReac)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: EquiRate, BackRate
      INTEGER :: iReac
      !iReac
      IF ( ReactionSystem(iReac)%bR ) THEN
        ReacConst = BackRate
      ELSE
        ReacConst = EquiRate * BackRate
      END IF
    END SUBROUTINE ReactionDirection   
    ! not nessessary anymore
    !SUBROUTINE vReactionDirection(k,kEq,kBr)
    !  REAL(RealKind)  :: k(:)
    !  REAL(RealKind), INTENT(IN)    :: kEq(:), kBr(:)
    !  
    !  k(RTind2%iDback) = kBr
    !  k(RTind2%iDforw) = kEq * kBr
    !END SUBROUTINE vReactionDirection   
    !
    !=====================================================================!
    ! ===  Multiplication with FACTOR
    !=====================================================================!
    !*****************************************************************
    SUBROUTINE EffectiveMolecularity(M,Conc,iReac,LWC)
      !OUT
      REAL(RealKind), INTENT(OUT) :: M
      !IN
      REAL(RealKind), INTENT(IN)  :: Conc(:)
      REAL(RealKind), INTENT(IN)  :: LWC
      INTEGER, INTENT(IN) :: iReac
      INTEGER :: i
      !TEMP
      !
      !print*, 'debugg ro2=',SUM(y_conc(RO2))
      !END DO
      !stop
      !
      !
      M = ONE
      IF (ReactionSystem(iReac)%Factor(1:1)=='$') THEN

        SELECT CASE (ReactionSystem(iReac)%Factor)
          CASE ('$H2')
            M=((mH2*mair)**fac_exp)*fac_A
          CASE ('$O2N2')
            M=(((mO2*mair)*(mN2*mair))**fac_exp)*fac_A
          CASE ('$M')
            M=(mair**fac_exp)*fac_A
          CASE ('$O2')
            M=((mO2*mair)**fac_exp)*fac_A
          CASE ('$N2')
            M=((mN2*mair)**fac_exp)*fac_A
          CASE ('$H2O')
            M=(mH2O**fac_exp)*fac_A
          CASE ('$RO2')
            M=SUM(Conc(RO2))
          CASE ('$O2O2')
            M=(((mO2*mair)**TWO)**fac_exp)*fac_A
          CASE ('$aH2O')
            !k=k*aH2OmolperL*LWC*mol2part
          CASE ('$RO2aq')
            M=SUM(Conc(RO2aq))
          CASE ('$+M','$(+M)')
            M=SUM(Conc)
            IF (ALLOCATED(ReactionSystem(iReac)%TBidx)) THEN
              M = M - SUM(( ONE - ReactionSystem(iReac)%TBalpha) &
                &             * Conc(ReactionSystem(iReac)%TBidx))
            END IF
            !if (ireac==590.or.ireac==591) print*, '590,591 M = ',M
          CASE DEFAULT
            WRITE(*,*) 'Reaction: ',iReac
            CALL FinishMPI()
            STOP 'Unknown FACTOR (error at Rate calc)'
        END SELECT

      END IF

    END SUBROUTINE EffectiveMolecularity
    SUBROUTINE vEffectiveMolecularity(M,Conc,LWC)
      !OUT
      REAL(RealKind), INTENT(OUT) :: M(neq)
      !IN
      REAL(RealKind), INTENT(IN)  :: Conc(:)
      REAL(RealKind), INTENT(IN)  :: LWC
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

    END SUBROUTINE vEffectiveMolecularity
    !*****************************************************************
    !
    !=========================================================================!
    !                  calculate sun 
    !=========================================================================!
    FUNCTION Zenith(Time) RESULT(Chi)
      !-----------------------------------------------------------------------!
      ! Input:
      !   - Time
      REAL(RealKind) :: Time
      !-----------------------------------------------------------------------!
      ! Output:
      !   - sun angle chi
      REAL(RealKind) :: Chi
      !-----------------------------------------------------------------------!
      ! Temporary variables:
      !INTEGER :: IDAT
      REAL(RealKind) :: LBGMT, LZGMT
      REAL(RealKind) :: ML
      ! 
      REAL(RealKind) :: GMT
      REAL(RealKind) :: RLT, RPHI
      !    
      INTEGER        :: IIYEAR, IYEAR, IMTH, IDAY, IIY, NYEARS, LEAP, NOLEAP
      REAL(RealKind) :: YREF,YR
      !   
      INTEGER        :: I, IJ, JD, IJD, IN
      REAL(RealKind) :: D, RML, W, WR, EC, EPSI, PEPSI, YT, CW, SW, SSW  & 
      &                  , EYT, FEQT1, FEQT2, FEQT3, FEQT4, FEQT5, FEQT6 &
      &                  , FEQT7, FEQT, EQT
      !         
      REAL(RealKind) :: REQT, RA, RRA, TAB, RDECL, DECL, ZPT, CSZ, ZR    &
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
      YREF =  2442047.5
      NYEARS = IYEAR - 1974
      LEAP = (NYEARS+1)/4
      IF(NYEARS.LE.-1) LEAP = (NYEARS-2)/4
      NOLEAP = NYEARS - LEAP
      YR = YREF + 365.*NOLEAP + 366.*LEAP
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
      D = JD + GMT/24.0
      !
      !      calc geom mean longitude
      ML = 279.2801988 + .9856473354*D + 2.267d-13*D*D
      RML = ML*DR
      !
      !      calc equation of time in sec
      !      w = mean long of perigee
      !      e = eccentricity
      !      epsi = mean obliquity of ecliptic
      W = 282.4932328 + 4.70684d-5*D + 3.39d-13*D*D
      WR = W*DR
      EC = 1.6720041d-2 - 1.1444d-9*D - 9.4d-17*D*D
      EPSI = 23.44266511 - 3.5626d-7*D - 1.23d-15*D*D
      PEPSI = EPSI*DR
      YT = (TAN(PEPSI/2.0))**2
      CW = COS(WR)
      SW = SIN(WR)
      SSW = SIN(2.0*WR)
      EYT = 2.*EC*YT
      FEQT1 = SIN(RML)*(-EYT*CW - 2.*EC*CW)
      FEQT2 = COS(RML)*(2.*EC*SW - EYT*SW)
      FEQT3 = SIN(2.*RML)*(YT - (5.*EC*EC/4.)*(CW*CW-SW*SW))
      FEQT4 = COS(2.*RML)*(5.*EC**2*SSW/4.)
      FEQT5 = SIN(3.*RML)*(EYT*CW)
      FEQT6 = COS(3.*RML)*(-EYT*SW)
      FEQT7 = -SIN(4.*RML)*(.5*YT*YT)
      FEQT = FEQT1 + FEQT2 + FEQT3 + FEQT4 + FEQT5 + FEQT6 + FEQT7
      EQT = FEQT*13751.0
      !
      !   convert eq of time from sec to deg
      REQT = EQT/240.
      !
      !   calc right ascension in rads
      RA = ML - REQT
      RRA = RA*DR
      !
      !   calc declination in rads, deg
      TAB = 0.43360*SIN(RRA)
      RDECL = ATAN(TAB)
      DECL = RDECL/DR
      !
      !   calc local hour angle
      LBGMT = 12.0 - EQT/3600. + LONG*24./360.
      LZGMT = 15.0*(GMT - LBGMT)
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
      Chi =  1.745329252D-02 * ZR/DR
    END FUNCTION Zenith
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

    !*****************************************************************
    SUBROUTINE PhoABCCompute(ReacConst,Constants,Time,chi)
      REAL(RealKind) :: ReacConst,Time
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Chi(:)
      !
      REAL(RealKind) :: ChiZ,yChiZ,EyChiZ
      !
      !CALL Zenith(Time,chi)
      IF (Chi(1)<PiHalf) THEN
        ChiZ=Chi(1)*Constants(3)
        IF (ChiZ<PiHalf) THEN
          yChiZ=Constants(2)*(One-One/COS(ChiZ))
          IF (yChiZ>-30.0d0) THEN
            EyChiZ=EXP(yChiZ)
          ELSE
            EyChiZ=9.357d-14
          END IF
        ELSE
          EyChiZ=9.357d-14
        END IF
        ReacConst=Dust*Constants(1)*EyChiz
      ELSE
        ReacConst=ZERO
      END IF
    END SUBROUTINE PhoABCCompute
    !
    FUNCTION vPhoABCCompute( Time, chi ) RESULT(kPHOTabc)
      REAL(RealKind)  :: kPHOTabc(nPHOTabc)
      REAL(RealKind), INTENT(IN)  :: Time, Chi(:)
      REAL(RealKind), DIMENSION(nPHOTabc) :: ChiZ, yChiZ, EyChiZ
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
    !*****************************************************************
   
    
    !*****************************************************************
    SUBROUTINE PhoABCompute(ReacConst,Constants,Time,Chi)
      REAL(RealKind) :: ReacConst,Time,Chi(:)
      REAL(RealKind) :: Constants(:)
      !
      !CALL Zenith(Time,chi)
      IF (Chi(1)<PiHalf) THEN
        ReacConst=Dust*Constants(1)*EXP(-Constants(2)/COS(Chi(1)))
      ELSE
        ReacConst=ZERO
      END IF
    END SUBROUTINE PhoABCompute
    FUNCTION vPhoABCompute(Time,Chi) RESULT(kPHOTab)
      REAL(RealKind)  :: kPHOTab(nPHOTab)
      REAL(RealKind), INTENT(IN)  :: Time, Chi(:)
     
      IF ( Chi(1) < PiHalf ) THEN
        kPHOTab = Dust * RTpar2%PHOTab(:,1)*EXP(-RTpar2%PHOTab(:,2)*Chi(2))
      ELSE
        kPHOTab = ZERO
      END IF
    END FUNCTION vPhoABCompute
    !*****************************************************************
   
    
    !*****************************************************************
    SUBROUTINE PhoMCMCompute(ReacConst,Constants,Time,chi)
      REAL(RealKind) :: ReacConst,Time
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Chi(:),chiz,ychiz
      !
      !CALL Zenith(Time,chi)
      !---  MCM version
      IF (Chi(1)<PiHalf) THEN
        chiz=EXP(-Constants(3)*(One/COS(chi(1))))
        ychiz=(COS(chi(1)))**(Constants(2))
        ReacConst=Dust*Constants(1)*ychiz*chiz
      ELSE
        ReacConst=ZERO
      END IF
    END SUBROUTINE PhoMCMCompute
    FUNCTION vPhoMCMCompute(Time,Chi) RESULT(kPHOTmcm)
      REAL(RealKind)  :: kPHOTmcm(nPHOTmcm)
      REAL(RealKind), INTENT(IN)  :: Time, Chi(:)
      REAL(RealKind), DIMENSION(nPHOTmcm) :: ChiZ, yChiZ
      
      !---  MCM version
      IF ( Chi(1) < PiHalf ) THEN
        ChiZ  = EXP( -RTpar2%PHOTmcm(:,3) * (Chi(2)) )
        yChiZ = Chi(2) ** RTpar2%PHOTmcm(:,2)
        kPHOTmcm = Dust * RTpar2%PHOTmcm(:,1) * yChiZ * ChiZ
      ELSE
        kPHOTmcm = ZERO
      END IF
    END FUNCTION vPhoMCMCompute
    !*****************************************************************
    !
    !
    !==========================================================================!
    ! ===  Constant reactions
    !==========================================================================!
    !*****************************************************************
    SUBROUTINE ConstCompute(ReacConst,Constants)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      ReacConst=Constants(1)
    END SUBROUTINE ConstCompute
    FUNCTION vConstCompute() RESULT(kCONST)
      REAL(RealKind)  :: kCONST(nCONST)
      kCONST = RTpar2%CONST
    END FUNCTION vConstCompute
    !*****************************************************************
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
    !*****************************************************************
    SUBROUTINE Temp1Compute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*EXP(-Constants(2)/Temp)
    END SUBROUTINE Temp1Compute
    FUNCTION vTemp1Compute(Temp) RESULT(kTEMP1)
      REAL(RealKind)  :: kTEMP1(nTEMP1)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      
      !print *, ' erste param = ',RTpar2%TEMP1(1,1) ,RTpar2%TEMP1(1,2)
      !print *, ' erste temp1 = ',RTpar2%TEMP1(1,1) * EXP(-RTpar2%TEMP1(1,2)*Temp(6))
      kTEMP1 = RTpar2%TEMP1(:,1) * EXP(-RTpar2%TEMP1(:,2)*Temp(6))
    END FUNCTION vTemp1Compute
    !*****************************************************************
    !
    !
    !
    !*****************************************************************
    SUBROUTINE Temp2Compute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*Temp*Temp*EXP(-Constants(2)/Temp)
    END SUBROUTINE Temp2Compute
    FUNCTION vTemp2Compute(Temp) RESULT(kTEMP2)
      REAL(RealKind)  :: kTEMP2(nTEMP2)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      !
      kTEMP2 = RTpar2%TEMP2(:,1) * Temp(2) * EXP( -RTpar2%TEMP2(:,2)*Temp(6) )
    END FUNCTION vTemp2Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Temp3Compute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*EXP(Constants(2)*(One/Temp-InvRefTemp))
    END SUBROUTINE Temp3Compute
    FUNCTION vTemp3Compute(Temp) RESULT(kTEMP3)
      REAL(RealKind)  :: kTEMP3(nTEMP3)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      !
      kTEMP3 = RTpar2%TEMP3(:,1) * EXP( RTpar2%TEMP3(:,2)*(Temp(6) - InvRefTemp) )
    END FUNCTION vTemp3Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Temp4Compute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*Temp*EXP(-Constants(2)/Temp)
    END SUBROUTINE Temp4Compute
    FUNCTION vTemp4Compute(Temp) RESULT(kTEMP4)
      REAL(RealKind)  :: kTEMP4(nTEMP4)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      !
      kTEMP4 = RTpar2%TEMP4(:,1) * Temp(1) * EXP( -RTpar2%TEMP4(:,2)*Temp(6) )
    END FUNCTION vTemp4Compute
    !*****************************************************************
    !
    !
    !   ***************************************************************
    !   ** Species nondimensional gibbs potentials                   **
    !   ***************************************************************
    SUBROUTINE GibbsFreeEnergie(Gibbs,T)
      REAL(RealKind) :: Gibbs(:)
      REAL(RealKind) :: T(:)
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
      REAL(RealKind) :: DGibbsdT(:)
      REAL(RealKind) :: T(:)
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
      REAL(RealKind) :: DelGibbs(:)
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
      REAL(RealKind) :: DiffDelGibbs(:)
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
      REAL(RealKind) :: T(8)
      !
      ! OUT
      REAL(RealKind) :: C(nspc)       ! molar heat capacities at constant pressure
      REAL(RealKind) :: H(nspc)       ! the standardstate molar enthalpy
      REAL(RealKind) :: S(nspc)       ! standard-state entropy at 298 K
      !
      REAL(RealKind) :: dHdT(nspc)    ! Enthaply derivative in dT [J/mol/K^2]
      REAL(RealKind) :: dGdT(nspc)    ! Gibbs potential derivative in dT [J/mol/K^2]
      REAL(RealKind) :: dCvdT(nspc)   ! Constant volume specific heat derivative in dT [J/mol/K]
      
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
    !
    !
    SUBROUTINE TempXCompute(k,rK_eq,kf,Const,T,iReac)
      !OUT:
      REAL(RealKind) :: k               ! Reaction rate constant
      REAL(RealKind) :: kf              ! Forward rate constant
      REAL(RealKind) :: rK_eq           ! reciprocal Equilibrium rate constant
      !IN:
      REAL(RealKind) :: Const(:)    ! Arrhenius parameter A,b,E_a
      REAL(RealKind) :: T(:)            ! Temperature array SIZE=8
      INTEGER :: iReac                  
      !temp
      REAL(RealKind) :: rRcT

      rRcT = rRcal*T(6)
      kf   = Const(1) * EXP(Const(2)*T(8) - Const(3)*rRcT)

      IF ( ReactionSystem(iReac)%bR ) THEN
        ! backward without extra coef, equiv constant nessesarry
        ! compute inverse equiv. constant
        ! explicite reverse reaction will be computed like 
        ! standart forward reaction
        rK_eq = EXP(+DelGFE(iReac)) * rFacEq**(-sumBAT(iReac))
        k     = kf * rK_eq
        !print*, ' deb:: rK_eq,Delgfe,rFacEq = ',rK_eq,rK_eq,rFacEq
      ELSE
        k     = Const(1) * EXP(Const(2)*T(8) - Const(3)*rRcT)
      END IF
    END SUBROUTINE TempXCompute
    !
    !
    SUBROUTINE DiffTempXCompute(Dkcoef,Const,T,iReac)
      REAL(RealKind) :: Dkcoef
      REAL(RealKind) :: DfRdT, DbRdT, DeRdT
      REAL(RealKind) :: Const(:)
      REAL(RealKind) :: T(:)              ! temperatur array
      INTEGER :: iReac
      !
      !TEMP
      REAL(RealKind) :: rRcT

      rRcT  = rRcal*T(6)
      DfRdT = Const(2) + Const(3)*rRcT     ! (21) Perini
     
      IF ( ReactionSystem(iReac)%bR ) THEN
        ! reverse reaction using equiv constant
        DeRdT  = sumBAT(iReac)*T(6) + DDelGFEdT(iReac)   ! (23) Perini
        Dkcoef = (DfRdT + DeRdT)*T(6)
      ELSE
        Dkcoef = DfRdT*T(6)
      END IF
  
    END SUBROUTINE DiffTempXCompute
    
   
    SUBROUTINE PressXCompute(k,k0,kinf,iR,T,Meff)
      REAL(RealKind) :: k
      !
      REAL(RealKind) :: T(:)
      REAL(RealKind) :: Meff
      INTEGER :: iR
      REAL(RealKind) :: FACtroe
      REAL(RealKind) :: logF_Troe, log10_Fcent, log10_Pr
      REAL(RealKind) :: cnd(3)
      !
      REAL(RealKind) :: k0 , k0M , kinf, rK_eq, Fcent, rRcT, FTL
      REAL(RealKind) :: High(3), Low(3)
    
      IF (ReactionSystem(iR)%Line3=='LOW') THEN
        High = ReactionSystem(iR)%Constants
        Low  = ReactionSystem(iR)%LowConst
      ELSE
        High = ReactionSystem(iR)%HighConst
        Low  = ReactionSystem(iR)%Constants
      END IF

      rRcT = rRcal*T(6)
      k0   = Low(1)  * EXP(  Low(2)*T(8) -  Low(3)*rRcT )
      kinf = High(1) * EXP( High(2)*T(8) - High(3)*rRcT )
      
      k0M  = k0*Meff
      Pr   = k0M / (kinf+k0M)
   
      IF (ALLOCATED(ReactionSystem(iR)%TroeConst)) THEN
        ! Troe form
        FTL = TroeFactor(iR,T)
      ELSE
        ! Lind form
        FTL = Pr
      END IF
     
      ! direction (forward - backward)
      IF ( ReactionSystem(iR)%bR ) THEN
        rK_eq = EXP(+DelGFE(iR)) * rFacEq**(-sumBAT(iR))
        k = kinf * FTL * rK_eq
      ELSE
        k = kinf * FTL
      END IF

    END SUBROUTINE PressXCompute


    SUBROUTINE NumDiffPressXCompute(DiffFactor_PD,k,iR,T,Meff)
      REAL(RealKind), INTENT(OUT) :: DiffFactor_PD

      REAL(RealKind), INTENT(IN)  :: k, Meff, T(:)

      REAL(RealKind) :: kdel, k0, kinf
      INTEGER        :: iR
      REAL(RealKind) :: Tdel(10)
      REAL(RealKind) :: delTemp

      delTemp = 1.0d-08

      CALL UpdateTempArray(Tdel , T(1)+delTemp)
      CALL PressXCompute(kdel, k0, kinf, iR, Tdel, Meff)

      DiffFactor_PD = (kdel - k) / delTemp
      DiffFactor_PD = DiffFactor_PD/k

    END SUBROUTINE NumDiffPressXCompute


    SUBROUTINE DiffPressXCompute(DiffFactor_PD,k0,kinf,iR,T,Meff)
      !OUT
      REAL(RealKind) :: DiffFactor_PD
      !IN
      INTEGER :: iR 
      REAL(RealKind) :: k0, kinf, Meff
      REAL(RealKind) :: cnd(3)
      REAL(RealKind) :: T(10)
      !Temp
      REAL(RealKind) :: Dk0dT, DeRdT, DkinfdT, rRcT, rkinfpk0M
      REAL(RealKind) :: dFTL_dT
      REAL(RealKind) :: log10_Prc, log10_FTroe, dlog10_PrdT

      REAL(RealKind) :: TempDiffFactor    ! forward or backward
      REAL(RealKind) :: tmplog10
      REAL(RealKind) :: DlogF_Troedlog_Pr, Dlog_F_TroedT
      REAL(RealKind) :: DF_PDdT

      REAL(RealKind) :: DF_PDdT_Troe
      !
      REAL(RealKind) :: High(3), Low(3)
      !
      IF (ReactionSystem(iR)%Line3=='LOW') THEN
        High = ReactionSystem(iR)%Constants
        Low  = ReactionSystem(iR)%LowConst
      ELSE
        High = ReactionSystem(iR)%HighConst
        Low  = ReactionSystem(iR)%Constants
      END IF
     
      ! Lind mechnism also nessesarry for Troe 
      rRcT    = rRcal*T(6)
      Dk0dT   =  Low(2) +  Low(3)*rRcT
      DkinfdT = High(2) + High(3)*rRcT

      rkinfpk0M = ONE / (kinf + k0*Meff)
      dFTL_dT   = Meff*k0*kinf*(rkinfpk0M*rkinfpk0M) * (Dk0dT - DkinfdT)*T(6)
      
      IF (ALLOCATED(ReactionSystem(iR)%TroeConst)) THEN
        ! Troe mechanism
        DF_PDdT = DiffTroeFactor(Dk0dT,DkinfdT,dFTL_dT,T)
      ELSE
        ! Lind mechnism
        DF_PDdT = dFTL_dT
      END IF
      
      IF (ReactionSystem(iR)%bR) THEN
        DeRdT  = - sumBAT(iR)*T(6) - DDelGFEdT(iR)
        TempDiffFactor = DkinfdT - DeRdT
      ELSE
        TempDiffFactor = DkinfdT
      END IF
      !
      DiffFactor_PD  = ( DF_PDdT + TempDiffFactor )*T(6)
    
    END SUBROUTINE DiffPressXCompute


    FUNCTION TroeFactor(iR,T)
      !OUT
      REAL(RealKind) :: TroeFactor
      !IN
      INTEGER        :: iR
      REAL(RealKind) :: T(:)
      !TEMP
      REAL(RealKind) :: Fcent, logF_Troe, FACtroe
      REAL(RealKind) :: Troe(4)


      Troe = ReactionSystem(iR)%TroeConst
    
      Fcent = (ONE - Troe(1)) * EXP(-T(1)/Troe(2)) &
            &      + Troe(1)  * EXP(-T(1)/Troe(3)) &
            &      +            EXP(-T(6)*Troe(4))
      
      log10_Fcent = LOG10(Fcent)
      log10_Pr    = LOG10(Pr/(ONE-Pr))

      cTroe = -0.40d0 - 0.67d0*log10_Fcent    ! will be used for deriv too 
      n1Troe =  0.75d0 - 1.27d0*log10_Fcent
      
      FACtroe = (log10_Pr + cTroe) / (n1Troe - dTroe*(log10_Pr + cTroe))

      logF_Troe  = log10_Fcent / (ONE + FACtroe*FACtroe)
      TroeFactor = Pr * TEN**logF_Troe
    
    END FUNCTION TroeFactor


    FUNCTION TroeFactorVec(T)
      !OUT
      REAL(RealKind) :: TroeFactorVec(RTind%nTroe)
      !IN
      REAL(RealKind) :: T(:)
      !TEMP
      REAL(RealKind) :: Fcent(RTind%nTroe),     &
                      & logF_Troe(RTind%nTroe), &
                      & FACtroe(RTind%nTroe)


    
      Fcent = (ONE - RTpar%T1) * EXP(-T(1)/RTpar%T2) &
            &      + RTpar%T1  * EXP(-T(1)/RTpar%T3) &
            &      +             EXP(-T(6)*RTpar%T4)
      
      vlog10_Fcent = LOG10(Fcent)
      vlog10_Pr    = LOG10(vPr(RTind%iTroe)/(ONE-vPr(RTind%iTroe)))

      vcTroe  = -0.40d0 - 0.67d0*vlog10_Fcent    ! will be used for deriv too 
      vn1Troe =  0.75d0 - 1.27d0*vlog10_Fcent
      
      FACtroe = (vlog10_Pr+vcTroe) / (vn1Troe-dTroe*(vlog10_Pr+vcTroe))

      logF_Troe  = vlog10_Fcent / (ONE+FACtroe*FACtroe)
      TroeFactorVec = vPr(RTind%iTroe) * TEN**logF_Troe

      !print*,' vlog10_Fcent(1),vlog10_Pr(1),FACtroe(1),logF_Troe(1)'
      !print*, vlog10_Fcent(1),vlog10_Pr(1),FACtroe(1),logF_Troe(1)
      !print*, vcTroe(1), vn1Troe(1), dTroe
      !stop
    
    END FUNCTION TroeFactorVec

    FUNCTION DiffTroeFactor(Dk0dT,DkinfdT,dFTL_dT,T)
      REAL(RealKind) :: DiffTroeFactor
      !IN
      REAL(RealKind) :: Dk0dT, DkinfdT, dFTL_dT, Pr
      REAL(RealKind) :: T(:)
      !TEMP
      REAL(RealKind) :: dlog10_PrdT, log10_Prc, tmplog10
      REAL(RealKind) :: DlogF_Troedlog_Pr, Dlog_F_TroedT
      REAL(RealKind) :: log10_FTroe

      ! Troe mechanism

      dlog10_PrdT = rln10 * (Dk0dT - DkinfdT) * T(6)

      log10_Prc   = log10_Pr + cTroe

      DlogF_Troedlog_Pr = - TWO*n1Troe*log10_Fcent*log10_Prc             &
                        &       * ( n1Troe - dTroe*log10_Prc )           &
                        &  / (  n1Troe**2 - TWO*n1Troe*dTroe*log10_Prc +  &
                        &      (dTroe**2 + ONE)*log10_Prc*log10_Prc  )**2

      Dlog_F_TroedT     = DlogF_Troedlog_Pr * dlog10_PrdT
      tmplog10          = log10_Prc   / (n1Troe - dTroe*log10_Prc)
      log10_FTroe       = log10_Fcent / (ONE   + (tmplog10*tmplog10))

      DiffTroeFactor    = TEN**log10_FTroe * (dFTL_dT + Pr*ln10*Dlog_F_TroedT)

    END FUNCTION DiffTroeFactor

    FUNCTION DTroeFactorVec(Dk0dT,DkinfdT,dFTL_dT,T)
      REAL(RealKind) :: DTroeFactorVec(RTind%nTroe)
      !IN
      REAL(RealKind) :: Dk0dT(:), DkinfdT(:), dFTL_dT(:)
      REAL(RealKind) :: T(:)
      !TEMP
      REAL(RealKind), DIMENSION(RTind%nTroe) :: dlog10_PrdT, log10_Prc, tmplog10
      REAL(RealKind), DIMENSION(RTind%nTroe) :: DlogF_Troedlog_Pr, Dlog_F_TroedT
      REAL(RealKind), DIMENSION(RTind%nTroe) :: log10_FTroe

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

    SUBROUTINE UpdateTempArray(TempArr,Temperature)
      REAL(RealKind) :: Temperature 
      REAL(RealKind) :: TempArr(10)
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
    END SUBROUTINE UpdateTempArray
    !
    !
    !-------------------------------------------------------------------------
    !---  Species internal energies in moles [J/mol]  
    !-------------------------------------------------------------------------
    SUBROUTINE InternalEnergy(U,T)
      !OUT
      REAL(RealKind) :: U(nspc)   
      !IN
      REAL(RealKind) :: T(:)
      
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
      REAL(RealKind) :: dUdT(nspc)   
      !IN
      REAL(RealKind) :: T(:)                      
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
      REAL(RealKind) :: d2UdT2(nspc)     !Constant volume specific heat’s derivative [J/mol/K2] 
      REAL(RealKind) :: T(:)                      
      !
      WHERE (SwitchTemp>T(1))
        d2UdT2 = lowB + TWO*lowC*T(1) + THREE*lowD*T(2) + FOUR*lowE*T(3) 
      ELSEWHERE    
        d2UdT2 = highB + TWO*highC*T(1) + THREE*highD*T(2) + FOUR*highE*T(3)
      END WHERE
      d2UdT2 = d2UdT2 * R
    END SUBROUTINE Diff2InternalEnergy
    !
    !===  Special reactions for Gas Phase Chemistry 
    !===  (Check: current%str_type = 'GAS' in ReadChem, nfra=1)
    !==========================================================================!
    ! Input:
    !   - Constants
    !   - Temperatur
    !   - mAir
    !--------------------------------------------------------------------------!
    ! Output:
    !   - Reaction constant
    !--------------------------------------------------------------------------!
    !
    !
    !
    !
    !
    !*****************************************************************
    SUBROUTINE TroeCompute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp,k1,k2
      REAL(RealKind) :: mAir
      !
      k1=mAir*Constants(1)*(Temp/300.0d0)**(-Constants(2))
      k2=Constants(3)*(Temp/300.0d0)**(-Constants(4))
      ReacConst=k1/(One+k1/k2)*0.6d0**(One/(One+(LOG10(k1/k2))*(LOG10(k1/k2))))
    END SUBROUTINE TroeCompute
    FUNCTION vTroeCompute(Temp,mAir) RESULT(kTROE)
      REAL(RealKind)  :: kTROE(nTROE)
      REAL(RealKind), INTENT(IN)  :: Temp(:) , mAir
      REAL(RealKind), DIMENSION(nTROE) :: k1, k2, log10_k1k2
      REAL(RealKind), PARAMETER   :: F = 0.6d0 
      
      k1 = RTpar2%TROE(:,1) * (Temp(1)*r300)**(-RTpar2%TROE(:,2)) * mAir
      k2 = RTpar2%TROE(:,3) * (Temp(1)*r300)**(-RTpar2%TROE(:,4))
      log10_k1k2 = LOG10(k1/k2)
      kTROE = k1 / (One + k1/k2) * F**(One / (One + log10_k1k2*log10_k1k2) )
    END FUNCTION vTroeCompute
    !*****************************************************************
   
    
    !*****************************************************************
    SUBROUTINE TroeMCMCompute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp,k1,k2,fc
      REAL(RealKind) :: mAir
      !
      k1=mAir*Constants(1)*(Temp/298.d0)**(Constants(2))*EXP(Constants(3)/Temp)
      k2=Constants(4)*(Temp/298.d0)**(Constants(5))*EXP(Constants(6)/Temp)
      fc=Constants(7)*EXP(Constants(8)/Temp)+Constants(9)*EXP(Temp/Constants(10))
      ReacConst=k1/(One+k1/k2)*                                                  &
      &         fc**(One/(One+(LOG10(k1/k2)/(0.75d0-1.27d0*LOG10(fc)))**2.0d0))
    END SUBROUTINE TroeMCMCompute
    FUNCTION vTroeMCMCompute(Temp,mAir) RESULT(kTROEmcm)
      REAL(RealKind)  :: kTROEmcm(nTROEmcm)
      REAL(RealKind), INTENT(IN)  :: Temp(:) , mAir
      REAL(RealKind), DIMENSION(nTROEmcm) :: k1, k2, Fc, tmpTROE
      REAL(RealKind), PARAMETER   :: n=0.75d0, d=1.27d0
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
    !*****************************************************************
  
 
    !*****************************************************************
    !---  Troe with variable F factor
    SUBROUTINE TroeFCompute(ReacConst,Constants,Temp,mAir) ! Barthel
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp,k1,k2
      REAL(RealKind) :: mAir
      !
      k1=mAir*Constants(1)*(Temp/300.0d0)**(-Constants(2))
      k2=Constants(3)*(Temp/300.0d0)**(-Constants(4))
      ReacConst=k1/(One+k1/k2)*Constants(5)**(One/(One+(LOG10(k1/k2))*(LOG10(k1/k2))))
    END SUBROUTINE TroeFCompute
    FUNCTION vTroeFCompute(Temp,mAir) RESULT(kTROEf)  ! Barthel
      REAL(RealKind)  :: kTROEf(nTROEf)
      REAL(RealKind), INTENT(IN)  :: Temp(:) , mAir
      REAL(RealKind), DIMENSION(nTROEf) :: k1, k2, log10_k1k2
      !
      k1 = RTpar2%TROEf(:,1)*(Temp(1)*r300)**(-RTpar2%TROEf(:,2)) * mAir
      k2 = RTpar2%TROEf(:,3)*(Temp(1)*r300)**(-RTpar2%TROEf(:,4))
      log10_k1k2 = LOG10(k1/k2)
      kTROEf = k1 / (One+k1/k2) * RTpar2%TROEf(:,5)**(One / (One + log10_k1k2*log10_k1k2))
    END FUNCTION vTroeFCompute
    !*****************************************************************
    
   
    !*****************************************************************
    !---  Troe equilibrium
    SUBROUTINE TroeEqCompute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp,k1,k2,k3
      REAL(RealKind) :: mAir
      !
      k1=mAir*Constants(1)*(Temp/300.d0)**(-Constants(2))
      k2=Constants(3)*(Temp/300.d0)**(-Constants(4))
      k3=k1/(One+k1/k2)*0.6d0**(One/(One+(LOG10(k1/k2))*(LOG10(k1/k2))))
      ReacConst=k3/(Constants(5)*EXP(Constants(6)/Temp))          
    END SUBROUTINE TroeEqCompute
    FUNCTION vTroeEqCompute(Temp,mAir) RESULT(kTROEq)
      REAL(RealKind)  :: kTROEq(nTROEq)
      REAL(RealKind), INTENT(IN)  :: Temp(:) , mAir
      REAL(RealKind), DIMENSION(nTROEq) :: k1, k2, k3, log10_k1k2
      REAL(RealKind), PARAMETER   :: F = 0.6d0 
      !
      k1 = RTpar2%TROEq(:,1)*(Temp(1)*r300)**(-RTpar2%TROEq(:,2)) * mAir
      k2 = RTpar2%TROEq(:,3)*(Temp(1)*r300)**(-RTpar2%TROEq(:,4))
      log10_k1k2 = LOG10(k1/k2)
      k3 = k1 / (One + k1/k2) * F**(One / (One + log10_k1k2*log10_k1k2))
      kTROEq = k3 / ( RTpar2%TROEq(:,5) * EXP(RTpar2%TROEq(:,6)*Temp(6)) )          
    END FUNCTION vTroeEqCompute
    !*****************************************************************
    
    !
    !*****************************************************************
    SUBROUTINE TroeEqfCompute(ReacConst,Constants,Temp,mAir) ! Barthel
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp,k1,k2,k3
      REAL(RealKind) :: mAir
      ! 
      k1=mAir*Constants(1)*(Temp/300.d0)**(-Constants(2))
      k2=Constants(3)*(Temp/300.0d0)**(-Constants(4))
      k3=k1/(One+k1/k2)*Constants(7)**(One/(One+(LOG10(k1/k2))*(LOG10(k1/k2))))
      ReacConst=k3/(Constants(5)*EXP(Constants(6)/Temp))          
    END SUBROUTINE TroeEqfCompute
    FUNCTION vTroeEqfCompute(Temp,mAir) RESULT(kTROEqf) ! Barthel
      REAL(RealKind)  :: kTROEqf(nTROEqf)
      REAL(RealKind), INTENT(IN)  :: Temp(:) , mAir
      REAL(RealKind), DIMENSION(nTROEqf) :: k1, k2, k3, log10_k1k2
      ! 
      k1 = RTpar2%TROEqf(:,1) * (Temp(1)*r300)**(-RTpar2%TROEqf(:,2)) * mAir
      k2 = RTpar2%TROEqf(:,3) * (Temp(1)*r300)**(-RTpar2%TROEqf(:,4))
      log10_k1k2 = LOG10(k1/k2)
      k3 = k1 / (One + k1/k2) * RTpar2%TROEqf(:,7)**(One / (One + log10_k1k2*log10_k1k2))
      kTROEqf = k3 / (RTpar2%TROEqf(:,5) * EXP(RTpar2%TROEqf(:,6)*Temp(6)))          
    END FUNCTION vTroeEqfCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    !---  modified Troe with variable F factor
    SUBROUTINE TroeXPCompute(ReacConst,Constants,Temp,mAir) ! Barthel
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp,k1,k2
      REAL(RealKind) :: mAir
      !
      k1=mAir*Constants(1)*EXP(-Constants(2)/Temp)
      k2=Constants(3)*EXP(-Constants(4)/Temp)
      ReacConst=k1/(One+k1/k2)*Constants(5)**(One/(One+          &
      &         (LOG10(k1/k2))*(LOG10(k1/k2)))) 
    END SUBROUTINE TroeXPCompute
    FUNCTION vTroeXPCompute(Temp,mAir) RESULT(kTROExp)  ! Barthel
      REAL(RealKind)  :: kTROExp(nTROExp)
      REAL(RealKind), INTENT(IN)  :: Temp(:) , mAir
      REAL(RealKind), DIMENSION(nTROExp) :: k1, k2, log10_k1k2
      !
      k1 = RTpar2%TROExp(:,1)*EXP(-RTpar2%TROExp(:,2)*Temp(6)) * mAir
      k2 = RTpar2%TROExp(:,3)*EXP(-RTpar2%TROExp(:,4)*Temp(6))
      log10_k1k2 = LOG10(k1/k2)
      kTROExp = k1 / (One + k1/k2) * RTpar2%TROExp(:,5)**(One / (One + log10_k1k2*log10_k1k2)) 
    END FUNCTION vTroeXPCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Spec1Compute(ReacConst,Constants,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: mAir
      !
      ReacConst=Constants(1)*(1.0d0+mAir*Constants(2))
    END SUBROUTINE Spec1Compute
    FUNCTION vSpec1Compute(mAir) RESULT(kSPEC1)
      REAL(RealKind)  :: kSPEC1(nSPEC1)
      REAL(RealKind), INTENT(IN)  :: mAir
      !
      kSPEC1 = RTpar2%SPEC1(:,1) * (ONE + mAir * RTpar2%SPEC1(:,2))
      !do globi=1,nSPEC1; print*, ' spec1 = ',RTind2%iSPEC1(globi), kSPEC1(globi); enddo
    END FUNCTION vSpec1Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Spec2Compute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: mAir
      ! 
      ReacConst=mAir*Constants(1)*(Temp/300.0d0)**Constants(2)
    END SUBROUTINE Spec2Compute
    !
    FUNCTION vSpec2Compute(Temp,mAir) RESULT(kSPEC2)
      REAL(RealKind)  :: kSPEC2(nSPEC2)
      REAL(RealKind), INTENT(IN)  :: Temp(:), mAir
      ! 
      !print*, 'mAir ,temp(1), pars :: ',mAir,temp(1),RTpar2%SPEC2(:,:)
      kSPEC2 = mAir * RTpar2%SPEC2(:,1) * (Temp(1)*r300)**RTpar2%SPEC2(:,2)
      !do globi=1,nSPEC2; print*, ' spec2 = ',RTind2%iSPEC2(globi), kSPEC2(globi); enddo
    END FUNCTION vSpec2Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Spec3Compute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp,k1,k2,k3
      REAL(RealKind) :: mAir
      !
      k1=Constants(1)*EXP(Constants(2)/Temp)
      k2=Constants(3)*EXP(Constants(4)/Temp)
      k3=mAir*Constants(5)*EXP(Constants(6)/Temp)
      ReacConst=k1+k3/(One+k3/k2)
    
    END SUBROUTINE Spec3Compute
    FUNCTION vSpec3Compute(Temp,mAir) RESULT(kSPEC3)
      REAL(RealKind)  :: kSPEC3(nSPEC3)
      REAL(RealKind), INTENT(IN)  :: Temp(:), mAir
      REAL(RealKind), DIMENSION(nSPEC3) :: k1, k2 ,k3
      !
      k1 = RTpar2%SPEC3(:,1)*EXP(RTpar2%SPEC3(:,2)*Temp(6))
      k2 = RTpar2%SPEC3(:,3)*EXP(RTpar2%SPEC3(:,4)*Temp(6))
      k3 = RTpar2%SPEC3(:,5)*EXP(RTpar2%SPEC3(:,6)*Temp(6)) * mAir
      kSPEC3 = k1 + k3 / (One + k3/k2)
      !do globi=1,nSPEC3; print*, ' spec3 = ',RTind2%iSPEC3(globi), kSPEC3(globi); enddo
    END FUNCTION vSpec3Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Spec4Compute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: mAir
      !
      ReacConst=Constants(1)*EXP(Constants(2)/Temp)+       &
      &         mAir*Constants(3)*EXP(Constants(4)/Temp)
    END SUBROUTINE Spec4Compute
    FUNCTION vSpec4Compute(Temp,mAir) RESULT(kSPEC4)
      REAL(RealKind)  :: kSPEC4(nSPEC4)
      REAL(RealKind), INTENT(IN)  :: Temp(:), mAir
      REAL(RealKind), DIMENSION(nSPEC4) :: k1, k2 ,k3
      !
      kSPEC4 = RTpar2%SPEC4(:,1) * EXP(RTpar2%SPEC4(:,2)*Temp(6)) +   &
      &        RTpar2%SPEC4(:,3) * EXP(RTpar2%SPEC4(:,4)*Temp(6)) * mAir 
    END FUNCTION vSpec4Compute
    !
    !*****************************************************************
    !
    !*****************************************************************
    SUBROUTINE Spec1MCMCompute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: mAir
      !
      ReacConst=Constants(1)*(1.0d0+mAir*Constants(2)/     &
      &         (Constants(3)*300.0d0/Temp))
    END SUBROUTINE Spec1MCMCompute
    FUNCTION vSpec1MCMCompute(Temp,mAir) RESULT(kSPEC1mcm)
      REAL(RealKind)  :: kSPEC1mcm(nSPEC1mcm)
      REAL(RealKind), INTENT(IN)  :: Temp(:), mAir
      
      kSPEC1mcm = RTpar2%SPEC1mcm(:,1)*(ONE + mAir*RTpar2%SPEC1mcm(:,2)*r300*Temp(1)/RTpar2%SPEC1mcm(:,3))
    END FUNCTION vSpec1MCMCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Spec2MCMCompute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*(Temp/300.0d0)**              &
      &         Constants(2)*EXP(Constants(3)/Temp)
    END SUBROUTINE Spec2MCMCompute
    FUNCTION vSpec2MCMCompute(Temp) RESULT(kSPEC2mcm)
      REAL(RealKind)  :: kSPEC2mcm(nSPEC2mcm)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      !
      kSPEC2MCM = RTpar2%SPEC2mcm(:,1)*(Temp(1)*r300)**RTpar2%SPEC2mcm(:,2) &
                * EXP(RTpar2%SPEC2mcm(:,3)*Temp(6))
    END FUNCTION vSpec2MCMCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Spec3MCMCompute(ReacConst,Constants,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: mAir
      !
      ReacConst=Constants(1)*(1.0d0+mAir/Constants(2))
    END SUBROUTINE Spec3MCMCompute
    FUNCTION vSpec3MCMCompute(mAir) RESULT(kSPEC3mcm)
      REAL(RealKind)  :: kSPEC3mcm(nSPEC3mcm)
      REAL(RealKind), INTENT(IN)  :: mAir
      !
      kSPEC3mcm = RTpar2%SPEC3mcm(:,1)*(ONE + mAir/RTpar2%SPEC3mcm(:,2))
    END FUNCTION vSpec3MCMCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Spec4MCMCompute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*(1.0d0+Constants(2)*                  &
      &         EXP(Constants(3)/Temp)*H2O)*EXP(Constants(4)/Temp)
    END SUBROUTINE Spec4MCMCompute
    FUNCTION vSpec4MCMCompute(Temp) RESULT(kSPEC4mcm)
      REAL(RealKind)  :: kSPEC4mcm(nSPEC4mcm)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      !
      kSPEC4mcm = RTpar2%SPEC4mcm(:,1)*(ONE + RTpar2%SPEC4mcm(:,2) &
      &         * EXP(RTpar2%SPEC4mcm(:,3)*Temp(6))*H2O)           &
      &         * EXP(RTpar2%SPEC4mcm(:,4)*Temp(6))
    END FUNCTION vSpec4MCMCompute
    !*****************************************************************
    !
    !*****************************************************************
    SUBROUTINE Spec5MCMCompute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp,k1,k2
      REAL(RealKind) :: mAir
      !
      k1=Constants(1)*mAir*0.21d0*EXP(Constants(2)/Temp)
      k2=Constants(3)*mAir*0.21d0*EXP(Constants(4)/Temp)
      ReacConst=(k1*(One-k2))
    END SUBROUTINE Spec5MCMCompute
    FUNCTION vSpec5MCMCompute(Temp,mAir) RESULT(kSPEC5mcm)
      REAL(RealKind)  :: kSPEC5mcm(nSPEC5mcm)
      REAL(RealKind), INTENT(IN)  :: Temp(:), mAir
      REAL(RealKind), DIMENSION(nSPEC5mcm) :: k1, k2
      REAL(RealKind), PARAMETER   :: F = 0.21d0
      !
      k1 = RTpar2%SPEC5mcm(:,1) * mAir * F * EXP(RTpar2%SPEC5mcm(:,2)*Temp(6))
      k2 = RTpar2%SPEC5mcm(:,3) * mAir * F * EXP(RTpar2%SPEC5mcm(:,4)*Temp(6))
      kSPEC5mcm = k1 * (One - k2)
    END FUNCTION vSpec5MCMCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Spec6MCMCompute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: K_1,K_2
      !
      K_1=Constants(1)*EXP(Constants(2)/Temp)
      K_2=Constants(3)*EXP(Constants(4)/Temp)
      ReacConst=K_1*(1.0d0-K_2)
    END SUBROUTINE Spec6MCMCompute
    FUNCTION vSpec6MCMCompute(Temp) RESULT(kSPEC6mcm)
      REAL(RealKind)  :: kSPEC6mcm(nSPEC6mcm)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      REAL(RealKind), DIMENSION(nSPEC6mcm) :: k1, k2
      !
      k1 = RTpar2%SPEC6mcm(:,1)*EXP(RTpar2%SPEC6mcm(:,2)*Temp(6))
      k2 = RTpar2%SPEC6mcm(:,3)*EXP(RTpar2%SPEC6mcm(:,4)*Temp(6))
      kSPEC6mcm = k1 * (ONE - k2)
    END FUNCTION vSpec6MCMCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Spec7MCMCompute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: K_1,K_2
      !
      K_1=Constants(1)*EXP(Constants(2)/Temp)
      K_2=Constants(3)*EXP(Constants(4)/Temp)
      ReacConst=K_1*(Constants(5)-Constants(6)/(One+K_2))
    END SUBROUTINE Spec7MCMCompute
    FUNCTION vSpec7MCMCompute(Temp) RESULT(kSPEC7mcm)
      REAL(RealKind)  :: kSPEC7mcm(nSPEC7mcm)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      REAL(RealKind), DIMENSION(nSPEC7mcm) :: k1, k2
      !
      k1 = RTpar2%SPEC7mcm(:,1)*EXP(RTpar2%SPEC7mcm(:,2)*Temp(6))
      k2 = RTpar2%SPEC7mcm(:,3)*EXP(RTpar2%SPEC7mcm(:,4)*Temp(6))
      kSPEC7mcm = k1 * (RTpar2%SPEC7mcm(:,5) - RTpar2%SPEC7mcm(:,6) / (One + k2))
    END FUNCTION vSpec7MCMCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Spec8MCMCompute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: mAir
      REAL(RealKind) :: K_1,K_2
      !
      K_1=Constants(1)*mAir*0.21d0*EXP(Constants(2)/Temp)
      K_2=Constants(3)*mAir*0.21d0*EXP(Constants(4)/Temp)
      ReacConst=K_1/((One+K_2)*Temp)
    END SUBROUTINE Spec8MCMCompute
    FUNCTION vSpec8MCMCompute(Temp,mAir) RESULT(kSPEC8mcm)
      REAL(RealKind)  :: kSPEC8mcm(nSPEC8mcm)
      REAL(RealKind), INTENT(IN)  :: Temp(:), mAir
      REAL(RealKind), DIMENSION(nSPEC8mcm) :: k1, k2
      REAL(RealKind), PARAMETER   :: F = 0.21d0
      !
      k1 = RTpar2%SPEC8mcm(:,1) * mAir * F * EXP(RTpar2%SPEC8mcm(:,2)*Temp(6))
      k2 = RTpar2%SPEC8mcm(:,3) * mAir * F * EXP(RTpar2%SPEC8mcm(:,4)*Temp(6))
      kSPEC8mcm = k1 / (One + k2) * Temp(6)
    END FUNCTION vSpec8MCMCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE S4H2OCompute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: mAir
      !
      ReacConst=Constants(1)*EXP(Constants(2)/Temp)+      & 
      &         mAir*Constants(3)*EXP(Constants(4)/Temp)
    END SUBROUTINE S4H2OCompute
    FUNCTION vS4H2OCompute(Temp,mAir) RESULT(kS4H2O)
      REAL(RealKind)  :: kS4H2O(nS4H2O)
      REAL(RealKind), INTENT(IN)  :: Temp(:), mAir
      
      kS4H2O = RTpar2%S4H2O(:,1) * EXP(RTpar2%S4H2O(:,2)*Temp(6))      & 
      &      + RTpar2%S4H2O(:,3) * EXP(RTpar2%S4H2O(:,4)*Temp(6)) * mAir
    END FUNCTION vS4H2OCompute
    !*****************************************************************
    !
    !
    !====================================================================!
    ! === Photolysis reactions for KPP systems
    !====================================================================!
    SUBROUTINE UpdateSun(Time,Sun)
      !--------------------------------------------------------------------
      ! Input:
      !   - Time
      REAL(Realkind) :: Time
      !--------------------------------------------------------------------!
      ! Output:
      !   - Sun
      REAL(Realkind)  :: Sun
      !--------------------------------------------------------------------!
      ! Temporary variables:
      REAL(Realkind) :: SunRise, SunSet
      REAL(Realkind) :: Thour, Tlocal, Ttmp
      !
      SunRise=4.50d0
      SunSet=19.50d0
      Thour=Time/3600.0d0
      Tlocal=Thour-FLOOR(Thour/24.0d0)*24.0d0
      !
      IF((Tlocal>=SunRise).AND.(Tlocal<=SunSet)) THEN
        Ttmp=(2.0d0*Tlocal-SunRise-SunSet)/(SunSet-SunRise);
        IF (Ttmp>0.0d0) THEN
          Ttmp=Ttmp*Ttmp;
        ELSE
          Ttmp=-Ttmp*Ttmp;
        END IF
        Sun=(1.0d0+COS(Pi*Ttmp))/2.0d0
      ELSE
        Sun=ZERO;
      END IF
    END SUBROUTINE UpdateSun
    !
    !
    SUBROUTINE PHOTO(ReacConst,Constants,Time)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(Realkind) :: Time
      REAL(Realkind) :: Sun
      !
      CALL UpdateSun(Time,Sun)
      ReacConst=Constants(1)*Sun
    END SUBROUTINE PHOTO
    !
    !
    SUBROUTINE PHOTO2(ReacConst,Constants,Time)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(Realkind) :: Time
      REAL(Realkind) :: Sun
      !
      CALL UpdateSun(Time,Sun)
      ReacConst=Constants(1)*Sun*Sun
    END SUBROUTINE PHOTO2
    !
    !
    SUBROUTINE PHOTO3(ReacConst,Constants,Time)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(Realkind) :: Time
      REAL(Realkind) :: Sun
      !
      CALL UpdateSun(Time,Sun)
      ReacConst=Constants(1)*Sun*Sun*Sun
    END SUBROUTINE PHOTO3
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
    SUBROUTINE T2H2OCompute(ReacConst,Constants,Temp)!,rHum)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
     ! REAL(RealKind) :: rHum
      !
     ! ReacConst=Rhum*Constants(1)*EXP(-Constants(2)/Temp)
    END SUBROUTINE T2H2OCompute
    !
    !
    !*****************************************************************
    SUBROUTINE DConstCompute(EquiRate,BackRate,Constants)
      REAL(RealKind) :: EquiRate
      REAL(RealKind) :: BackRate
      REAL(RealKind) :: Constants(:)
      !
      EquiRate=Constants(1)
      BackRate=Constants(2)
    END SUBROUTINE DConstCompute
    !SUBROUTINE vDConstCompute(kf,kb)
    FUNCTION vDConstCompute() RESULT(k)
      !REAL(RealKind)  :: kf(nDCONST_f), kb(nDCONST_b)
      REAL(RealKind)  :: k(nDCONST,2)
      !print*, ' index const f =', RTind2%iDCONST_f
      !print*, ' index const b =', RTind2%iDCONST_b
      !
      !kb = RTpar2%DCONST(RTind2%iDCONST_b,2)
      !kf = RTpar2%DCONST(RTind2%iDCONST_f,1) * kb
      k(:,2) = RTpar2%DCONST(:,2)
      k(:,1) = RTpar2%DCONST(:,1) * k(:,2)
    END FUNCTION vDConstCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE DTempCompute(EquiRate,BackRate,Constants,Temp)
      REAL(RealKind) :: EquiRate
      REAL(RealKind) :: BackRate
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      EquiRate=Constants(1)*EXP(Constants(2)*(ONE/Temp-InvRefTemp))
      BackRate=Constants(3)
    END SUBROUTINE DTempCompute
    FUNCTION vDTempCompute(Temp) RESULT(k)
      REAL(RealKind)  :: k(nDTEMP,2)
      REAL(RealKind), INTENT(IN) :: Temp(:)
      !
      k(:,2) = RTpar2%DTEMP(:,3)
      k(:,1) = RTpar2%DTEMP(:,1)*EXP(RTpar2%DTEMP(:,2)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vDTempCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE DTemp2Compute(EquiRate,BackRate,Constants,Temp)  
      REAL(RealKind) :: EquiRate
      REAL(RealKind) :: BackRate
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      EquiRate=Constants(1)*EXP(Constants(2)*(ONE/Temp-InvRefTemp))
      BackRate=Constants(3)*EXP(Constants(4)*(ONE/Temp-InvRefTemp))
    END SUBROUTINE DTemp2Compute
    FUNCTION vDTemp2Compute(Temp)   RESULT(k)
      REAL(RealKind)  :: k(nDTEMP2,2)
      REAL(RealKind), INTENT(IN) :: Temp(:)
      !
      k(:,2) = RTpar2%DTEMP2(:,3)*EXP(RTpar2%DTEMP2(:,4)*(Temp(6)-InvRefTemp))
      k(:,1) = RTpar2%DTEMP2(:,1)*EXP(RTpar2%DTEMP2(:,2)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vDTemp2Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE DTemp3Compute(EquiRate,BackRate,Constants,Temp)
      REAL(RealKind) :: EquiRate
      REAL(RealKind) :: BackRate
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      ! 
      EquiRate=Constants(1)*EXP(Constants(2)*(RefTemp/Temp-ONE)              &
      &                         +Constants(3)*(One-RefTemp/Temp              &
      &                                           +LOG10(RefTemp/Temp)))
      BackRate=Constants(4)
    END SUBROUTINE DTemp3Compute
    FUNCTION vDTemp3Compute(Temp) RESULT(k)
      REAL(RealKind)  :: k(nDTEMP3,2)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      ! 
      k(:,2) = RTpar2%DTEMP3(:,4)
      k(:,1) = RTpar2%DTEMP3(:,1) * EXP(RTpar2%DTEMP3(:,2)*(RefTemp*Temp(6)-ONE)         &
      &      + RTpar2%DTEMP3(:,3)*(One-RefTemp*Temp(6) + LOG10(RefTemp*Temp(6)))) * k(:,2)
    END FUNCTION vDTemp3Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE MeskhidzeCompute(EquiRate,BackRate,Constants,Temp)
      REAL(RealKind) :: EquiRate
      REAL(RealKind) :: BackRate
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      EquiRate=Constants(1)*(Temp/RefTemp)**Constants(2)        &
      &        *EXP(Constants(3)*(One/Temp-One/RefTemp)) 
      BackRate=Constants(4)*EXP(Constants(5)                    &
      &        *(One/Temp-One/RefTemp))*Constants(7)!*(ActivityHp)**m
    END SUBROUTINE MeskhidzeCompute
    FUNCTION vMeskhidzeCompute(Temp) RESULT(k)
      REAL(RealKind)  :: k(nMeskhidze,2)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      !
      k(:,2) = RTpar2%Meskhidze(:,4) * EXP(RTpar2%Meskhidze(:,5)  &
      &      * (Temp(6)-InvRefTemp)) * RTpar2%Meskhidze(:,7)!*(ActivityHp)**m
      k(:,1) = RTpar2%Meskhidze(:,1) * (Temp*InvRefTemp)**RTpar2%Meskhidze(:,2)  &
      &      * EXP(RTpar2%Meskhidze(:,3)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vMeskhidzeCompute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE DTemp4Compute(EquiRate,BackRate,Constants,Temp)
      REAL(RealKind) :: EquiRate
      REAL(RealKind) :: BackRate
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      EquiRate=Constants(1)*EXP(Constants(2)*(Temp*InvRefTemp-One)     &
      &        +Constants(3)*(ONE+LOG(Temp*InvRefTemp)-Temp*InvRefTemp))
      BackRate=ONE  ! OSSI
    END SUBROUTINE DTemp4Compute
    FUNCTION vDTemp4Compute(Temp) RESULT(k)
      REAL(RealKind)  :: k(nDTEMP4,2)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      !
      k(:,2) = ONE  ! OSSI
      k(:,1) = RTpar2%DTEMP4(:,1)*EXP(RTpar2%DTEMP4(:,2)                &
      &      * (Temp(1)*InvRefTemp-One) + RTpar2%DTEMP4(:,3)            &
      &      * (ONE+LOG(Temp(1)*InvRefTemp)-Temp(1)*InvRefTemp)) * k(:,2)
    END FUNCTION vDTemp4Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE DTemp5Compute(EquiRate,BackRate,Constants,Temp)
      REAL(RealKind) :: EquiRate
      REAL(RealKind) :: BackRate
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      EquiRate=Constants(1)*(Temp*InvRefTemp)**Constants(2)  &
      &        *EXP(Constants(3)*(ONE/Temp-InvRefTemp))
      BackRate=ONE  ! OSSI
    END SUBROUTINE DTemp5Compute      
    FUNCTION vDTemp5Compute(Temp) RESULT(k)
      REAL(RealKind)  :: k(nDTEMP5,2)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      !
      k(:,2) = ONE  ! OSSI
      k(:,1) = RTpar2%DTEMP5(:,1)*(Temp*InvRefTemp)**RTpar2%DTEMP4(:,2)  &
      &      * EXP(RTpar2%DTEMP5(:,3)*(Temp(6)-InvRefTemp)) * k(:,2)
    END FUNCTION vDTemp5Compute      
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Aspec1Compute(ReacConst,Constants,Temp,y_conc)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: y_conc(:)
      !
      IF (y_conc(Hp_ind)>ZERO) THEN
          ReacConst=y_conc(Hp_ind)*(Constants(1)*               &
            &              EXP(Constants(2)*(ONE/Temp-InvRefTemp)))/  &
            &              (ONE+13.0d0*y_conc(Hp_ind))
     END IF
    END SUBROUTINE Aspec1Compute
    FUNCTION vAspec1Compute(Temp,Hp) RESULT(kASPEC1)
      REAL(RealKind)  :: kASPEC1(nASPEC1)
      REAL(RealKind), INTENT(IN)  :: Temp(:), Hp
      REAL(RealKind), PARAMETER   :: x = 13.0d0
      !
      kASPEC1 = Hp*(RTpar2%ASPEC1(:,1)*EXP(RTpar2%ASPEC1(:,2) & 
      &       * (Temp(6)-InvRefTemp)))/(ONE + x*Hp)
    END FUNCTION vAspec1Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Aspec2Compute(ReacConst,Constants,Temp,y_conc)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: y_conc(:)
      !
      IF (y_conc(Hp_ind)>ZERO) THEN
          ReacConst=y_conc(Hp_ind)**Constants(2) *                       &
            &         (Constants(1)*EXP(Constants(3)*(ONE/Temp-InvRefTemp)))
     END IF
    END SUBROUTINE Aspec2Compute
    FUNCTION vAspec2Compute(Temp,Hp) RESULT(kASPEC2)
      REAL(RealKind)  :: kASPEC2(nASPEC2)
      REAL(RealKind), INTENT(IN)  :: Temp(:), Hp
      !
      kASPEC2 = Hp**RTpar2%ASPEC2(:,2) *                       &
      &         (RTpar2%ASPEC2(:,1)*EXP(RTpar2%ASPEC2(:,3)*(Temp(6)-InvRefTemp)))
    END FUNCTION vAspec2Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE Aspec3Compute(ReacConst,Constants,Temp,y_conc)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: y_conc(:)
      !
      IF (y_conc(Hp_ind)>ZERO) THEN
        ReacConst=Constants(1)*EXP(Constants(2)*(-LOG10(y_conc(Hp_ind))))
      END IF
    END SUBROUTINE Aspec3Compute
    FUNCTION vAspec3Compute(Temp,Hp) RESULT(kASPEC3)
      REAL(RealKind)  :: kASPEC3(nASPEC3)
      REAL(RealKind), INTENT(IN)  :: Temp(:), Hp

      kASPEC3 = RTpar2%ASPEC3(:,1)*EXP(RTpar2%ASPEC3(:,2)*(-LOG10(Hp)))
    END FUNCTION vAspec3Compute
    !*****************************************************************
    !
    !
    !*****************************************************************
    SUBROUTINE T1H2OCompute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*EXP(-Constants(2)/Temp)
    END SUBROUTINE T1H2OCompute
    FUNCTION vT1H2OCompute(Temp) RESULT(kT1H2O)
      REAL(RealKind)  :: kT1H2O(nT1H2O)
      REAL(RealKind), INTENT(IN)  :: Temp(:)
      
      kT1H2O = RTpar2%T1H2O(:,1)*EXP(-RTpar2%T1H2O(:,2)*Temp(6))
    END FUNCTION vT1H2OCompute
    !*****************************************************************
    !
    !
    ! ===============================================================
    ! === Solid reactions
    ! ===============================================================
    SUBROUTINE EquiCompute(ReacConst,Constants)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      !
      ReacConst=Constants(1)
    END SUBROUTINE EquiCompute

    SUBROUTINE PrintReaction(iR)
      INTEGER :: iR

      WRITE(*,*) ' ********************************************************************************************'
      WRITE(*,*) '  ReactionNumber= ', iR
      WRITE(*,*) '  Reaction      = ', TRIM(ReactionSystem(iR)%Line1)
      WRITE(*,*) '  FACTOR        = ', TRIM(ReactionSystem(iR)%Line2)
      WRITE(*,*) '  TYPE          = ', TRIM(ReactionSystem(iR)%Type)
      WRITE(*,*) '  TYPE Constant = ', TRIM(ReactionSystem(iR)%TypeConstant)
      WRITE(*,*) '  Constants     = ', ReactionSystem(iR)%Constants
      WRITE(*,*) ' ********************************************************************************************'
      WRITE(*,*) 
    END SUBROUTINE PrintReaction

  END MODULE Rates_Mod
