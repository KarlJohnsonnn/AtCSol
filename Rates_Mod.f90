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
    REAL(RealKind) :: cTroe, nTroe, dTroe
    REAL(RealKind) :: log10_Pr, log10_Fcent, Pr
    !
    CONTAINS
    !
    !
    !======================================================================!
    !      Calculate the Rates for current concentraions and time
    !======================================================================!
    SUBROUTINE Rates(Time,Y_in,Rate,DRatedT)
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
      REAL(RealKind) :: chi, actLWC
      REAL(RealKind) :: T(8)
      REAL(RealKind) :: k    ! tmp buffer for reaction const
      REAL(RealKind) :: Prod, AquaFac
      REAL(RealKind) :: DkdT    ! tmp buffer for reaction const
      REAL(RealKind) :: Meff
      INTEGER :: iReac, ii, j
      ! 
      ! combustion
      REAL(RealKind) :: Temp_in

      !==================================================================!
      !===============calc rates for ReactionSystem======================!
      !==================================================================!
     
      TimeRateA = MPI_WTIME()
      
      Rate    = ZERO
      DRatedT = ZERO

      ! --- Compute zenith for photo reactions
      chi     = Zenith(Time)
     
      ! --- Compute the liquid water content 
      actLWC  = pseudoLWC(Time)
      
      IF ( combustion ) THEN
        
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

      ELSE

        T(1)   = Temp    ! = 280 [K]
        y_conc = Y_in

      END IF
      !print*, '    C = ',y_conc([2831,945,2629])
      !print*, ' ctot = ', SUM(y_conc(1:nspc))

      !print*,
      ! --- compute rate of all reactions (gas,henry,aqua,diss)
      LOOP_OVER_ALL_REACTIONS: DO ii = 1 , loc_rateCnt
        
        ! if more than one processor is used split rate calculation
        !IF ( MPI_np>1 .AND. ParOrdering>=0 )  THEN
        !  iReac = loc_RatePtr(ii)
        !ELSE
        iReac = ii
        !END IF
       
        ! ====== Computing effective molecularity 
        CALL EffectiveMolecularity( Meff , y_conc(1:nspc) , iReac , actLWC )
        
        ! ====== Compute the rate constant for specific reaction type ===
        CALL ComputeRateConstant( k , DkdT , T , Time , chi , mAir , iReac , y_conc , Meff )

        AQUATIC_REACTIONS: IF ( ntAqua > 0 ) THEN

          IF      ( ReactionSystem(iReac)%Type == 'DISS' .OR. &
          &         ReactionSystem(iReac)%Type == 'AQUA' ) THEN
            ! ===== correct unit of concentrations
            AquaFac = ONE / ( (actLWC*mol2part)**ReactionSystem(iReac)%SumAqCoef )
            k       = k * AquaFac
          END IF

          IF ( ReactionSystem(iReac)%Type == 'HENRY' ) THEN
            !=== Compute Henry mass transfer coefficient
            CALL MassTransfer( k , T(1) , Time , iReac , actLWC )
          END IF

        END IF AQUATIC_REACTIONS
  
        ! === Calculate the product of concentrations 
        Prod = ONE
        DO j=A%RowPtr(iReac),A%RowPtr(iReac+1)-1
          IF      ( A%Val(j) == ONE ) THEN
            Prod = Prod * y_conc(A%ColInd(j))
          ELSE IF ( A%Val(j) == TWO ) THEN
            Prod = Prod * y_conc(A%ColInd(j))*ABS(y_conc(A%ColInd(j)))
          ELSE
            Prod = Prod * y_conc(A%ColInd(j))**A%Val(j)
          END IF  
        END DO
       
        Rate(iReac) = Meff * k * Prod
        IF (combustion) DRatedT(iReac) = DkdT

        !WRITE(987,*) 'DBg:: i, k , prd, Rate =', iReac, k, Prod, Rate(iReac)

      END DO LOOP_OVER_ALL_REACTIONS
      !stop 'ratesmod'

      TimeRateE=MPI_WTIME()
      TimeRates=TimeRates+(TimeRateE-TimeRateA)
      !
      
      ! gather the values of the other processes
      !CALL GatherAllPartitions(Rate,MyParties)
      !IF (combustion) CALL GatherAllPartitions(DRatedT,MyParties)
    END SUBROUTINE Rates
   !
    !
    !=====================================================================!
    ! === Converts the mass for Henry reactions Gas->Aqua , Aqua->Gas
    !=====================================================================!
    SUBROUTINE MassTransfer(k,Temp,Time,iReac,actLWC)
      REAL(RealKind) :: k
      REAL(RealKind) :: Temp, Time, actLWC
      INTEGER :: iReac
      !
      REAL(RealKind) :: term_diff
      REAL(RealKind) :: term_accom
      REAL(RealKind) :: kmt
      !
      !
      !---------------------------------------------------------------------------
      term_diff   = y_c1(ReactionSystem(iReac)%HenrySpc)               ! diffusion term
      term_accom  = y_c2(ReactionSystem(iReac)%HenrySpc) / SQRT(Temp)  ! accom term
      !--------------------------------------------------------------------------!
      !
      ! Compute new wet radius for droplett class iFrac
      SPEK(1)%wetRadius = (Pi34*actLWC/SPEK(1)%Number)**(rTHREE)
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
        k = milli * kmt * actLWC ! orginal
        !
      ! direaction AquaSpecies-->GasSpecies  
      ELSE 
        k = kmt / ( k * GasConst_R * Temp)   !()=HenryConst*GasConstant*Temperatur
      END IF
      !print*, 'Debug:: gasconst=',GasConst_R
      !print*, 'Debug::     y_c1=',y_c1(ReactionSystem(iReac)%HenrySpc)    
      !print*, 'Debug::     y_c2=',y_c2(ReactionSystem(iReac)%HenrySpc)
    END SUBROUTINE MassTransfer
    !
    !=======================================================================!
    ! ===  Select the Type of the Constant and calculate the value
    !=======================================================================!
    SUBROUTINE ComputeRateConstant(k,DkdT,T,Time,chi,mAir,iReac,y_conc,Meff)
      REAL(RealKind), INTENT(INOUT) :: k, DkdT
      REAL(RealKind), INTENT(IN) :: Time, mAir, chi
      REAL(RealKind), INTENT(IN) :: T(8)
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
          CALL NumDiffPressXCompute(DkdT,k,iReac,T,Meff)
          !CALL DiffPressXCompute(DkdT,k0,kinf,iReac,T,Meff)
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
    !
    !
    SUBROUTINE ReactionDirection(ReacConst,EquiRate,BackRate,iReac)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: EquiRate, Backrate
      INTEGER :: iReac
      !
      IF ( ReactionSystem(iReac)%bR ) THEN
        ReacConst = BackRate
      ELSE
        ReacConst = EquiRate * BackRate
      END IF
    END SUBROUTINE ReactionDirection   
    !
    !=====================================================================!
    ! ===  Multiplication with FACTOR
    !=====================================================================!
    SUBROUTINE EffectiveMolecularity(M,Conc,iReac,LWC)
      !OUT
      REAL(RealKind), INTENT(OUT) :: M
      !IN
      REAL(RealKind), INTENT(IN)  :: Conc(:)
      REAL(RealKind), INTENT(IN)  :: LWC
      INTEGER, INTENT(IN) :: iReac
      INTEGER :: j
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
            !k=k*aH2OmolperL*actLWC*mol2part
          CASE ('$RO2aq')
            M=SUM(Conc(RO2aq))
          CASE ('$+M','$(+M)')
            M=SUM(Conc)
            IF (ALLOCATED(ReactionSystem(iReac)%TBidx)) THEN
              M = M - SUM(( ONE - ReactionSystem(iReac)%TBalpha) &
                &             * Conc(ReactionSystem(iReac)%TBidx))
              !do j=1,SIZE(ReactionSystem(iReac)%TBalpha)
               ! print*, 'tb = ',iReac,ReactionSystem(iReac)%TBidx(j),&
              !    &             ONE - ReactionSystem(iReac)%TBalpha(j)
              !end do
            END IF
            !IF( mod(iReac,2)==0.) print*, 'DEBUG::SCmodule    M=',iReac,M
          CASE DEFAULT
            WRITE(*,*) 'Reaction: ',iReac
            CALL FinishMPI()
            STOP 'Unknown FACTOR (error at Rate calc)'
        END SELECT

      END IF

      IF (ReactionSystem(iReac)%nInActEd/=0) THEN
        SELECT CASE (ReactionSystem(iReac)%InActEductSpc(1))
          CASE ('[H2O]')
            M=M*H2O
          CASE ('[N2]')
            M=M*N2
          CASE ('[O2]')
            M=M*O2
          CASE ('[aH2O]')
            M=M*aH2OmolperL*LWC*mol2part
          CASE DEFAULT
            !STOP
        END SELECT
      END IF
      !
    END SUBROUTINE EffectiveMolecularity
    !
    !=========================================================================!
    !                  calculate sun 
    !=========================================================================!
    REAL(RealKind) FUNCTION Zenith(Time)
      !-----------------------------------------------------------------------!
      ! Input:
      !   - Time
      REAL(RealKind) :: Time
      !-----------------------------------------------------------------------!
      ! Output:
      !   - sun angle chi
      !REAL(RealKind) :: Zenith 
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
      GMT = Time / 3600.0D0
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
      !chi =  1.745329252D-02 * ZR/DR
      Zenith =  1.745329252D-02 * ZR/DR
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
    SUBROUTINE PhoABCCompute(ReacConst,Constants,Time,chi)
      REAL(RealKind) :: ReacConst,Time
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Chi
      !
      REAL(RealKind) :: ChiZ,yChiZ,EyChiZ
      !
      !CALL Zenith(Time,chi)
      IF (Chi<PiHalf) THEN
        ChiZ=Chi*Constants(3)
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
    !
    SUBROUTINE PhoABCompute(ReacConst,Constants,Time,chi)
      REAL(RealKind) :: ReacConst,Time,Chi
      REAL(RealKind) :: Constants(:)
      !
      !CALL Zenith(Time,chi)
      IF (Chi<PiHalf) THEN
        ReacConst=Dust*Constants(1)*EXP(-Constants(2)/COS(Chi))
      ELSE
        ReacConst=ZERO
      END IF
    END SUBROUTINE PhoABCompute
    !
    !
    SUBROUTINE PhoMCMCompute(ReacConst,Constants,Time,chi)
      REAL(RealKind) :: ReacConst,Time
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Chi,chiz,ychiz
      !
      !CALL Zenith(Time,chi)
      !---  MCM version
      IF (Chi<PiHalf) THEN
        chiz=EXP(-Constants(3)*(One/COS(chi)))
        ychiz=(COS(chi))**(Constants(2))
        ReacConst=Dust*Constants(1)*ychiz*chiz
      ELSE
        ReacConst=ZERO
      END IF
    END SUBROUTINE PhoMCMCompute
    !
    !
    !==========================================================================!
    ! ===  Constant reactions
    !==========================================================================!
    SUBROUTINE ConstCompute(ReacConst,Constants)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      ReacConst=Constants(1)
    END SUBROUTINE ConstCompute
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
    SUBROUTINE Temp1Compute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*EXP(-Constants(2)/Temp)
    END SUBROUTINE Temp1Compute
    !
    !
    SUBROUTINE Temp2Compute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*Temp*Temp*EXP(-Constants(2)/Temp)
    END SUBROUTINE Temp2Compute
    !
    !
    SUBROUTINE Temp3Compute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*EXP(Constants(2)*(One/Temp-InvRefTemp))
    END SUBROUTINE Temp3Compute
    !
    !
    SUBROUTINE Temp4Compute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*Temp*EXP(-Constants(2)/Temp)
    END SUBROUTINE Temp4Compute
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
        &               - SUM( BA%Val(from:to) * GFE(BA%ColInd(from:to)) )
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
      !TEMP:
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
      ELSE 
        k     = kf
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
        DeRdT  = - sumBAT(iReac)*T(6) - DDelGFEdT(iReac)   ! (23) Perini
        Dkcoef = (DfRdT - DeRdT)*T(6)
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
      REAL(RealKind) :: Tdel(8)
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
      REAL(RealKind) :: T(8)
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
      nTroe =  0.75d0 - 1.27d0*log10_Fcent
      dTroe =  0.14d0
      
      FACtroe = (log10_Pr + cTroe) / (nTroe - dTroe*(log10_Pr + cTroe))

      logF_Troe  = log10_Fcent / (ONE + FACtroe*FACtroe)
      TroeFactor = Pr * TEN**logF_Troe
    
    END FUNCTION TroeFactor


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

      DlogF_Troedlog_Pr = - TWO*nTroe*log10_Fcent*log10_Prc             &
                        &       * ( nTroe - dTroe*log10_Prc )           &
                        &  / (  nTroe**2 - TWO*nTroe*dTroe*log10_Prc +  &
                        &      (dTroe**2 + ONE)*log10_Prc*log10_Prc  )**2

      Dlog_F_TroedT     = DlogF_Troedlog_Pr * dlog10_PrdT
      tmplog10          = log10_Prc   / (nTroe - dTroe*log10_Prc)
      log10_FTroe       = log10_Fcent / (ONE   + (tmplog10*tmplog10))

      DiffTroeFactor    = TEN**log10_FTroe * (dFTL_dT + Pr*ln10*Dlog_F_TroedT)

    END FUNCTION DiffTroeFactor


    SUBROUTINE UpdateTempArray(TempArr,Temperature)
      REAL(RealKind) :: Temperature 
      REAL(RealKind) :: TempArr(8)
      !
      INTEGER :: i
      !
      TempArr(1) = Temperature                ! T
      DO i=2,5
        TempArr(i) = TempArr(i-1)*Temperature   ! T^2 ... T^5
      END DO
      TempArr(6) = ONE / Temperature          ! 1/T
      TempArr(7) = ONE / TempArr(2)           ! 1/T^2
      TempArr(8) = LOG(Temperature)           ! ln(T)
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
    !
    !
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
    !
    !
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
    !
    !
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
    !
    !
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
    !
    !
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
    !
    !
    SUBROUTINE Spec1Compute(ReacConst,Constants,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: mAir
      !
      ReacConst=Constants(1)*(1.0d0+mAir*Constants(2))
    END SUBROUTINE Spec1Compute
    !
    !
    SUBROUTINE Spec2Compute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: mAir
      ! 
      ReacConst=mAir*Constants(1)*(Temp/300.0d0)**Constants(2)
    END SUBROUTINE Spec2Compute
    !
    !
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
    !
    !
    SUBROUTINE Spec4Compute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: mAir
      !
      ReacConst=Constants(1)*EXP(Constants(2)/Temp)+       &
      &         mAir*Constants(3)*EXP(Constants(4)/Temp)
    END SUBROUTINE Spec4Compute
    !
    !
    SUBROUTINE Spec1MCMCompute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: mAir
      !
      ReacConst=Constants(1)*(1.0d0+mAir*Constants(2)/     &
      &         (Constants(3)*300.0d0/Temp))
    END SUBROUTINE Spec1MCMCompute
    !
    !
    SUBROUTINE Spec2MCMCompute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*(Temp/300.0d0)**              &
      &         Constants(2)*EXP(Constants(3)/Temp)
    END SUBROUTINE Spec2MCMCompute
    !
    !
    SUBROUTINE Spec3MCMCompute(ReacConst,Constants,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: mAir
      !
      ReacConst=Constants(1)*(1.0d0+mAir/Constants(2))
    END SUBROUTINE Spec3MCMCompute
    !
    !
    SUBROUTINE Spec4MCMCompute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*(1.0d0+Constants(2)*                  &
      &         EXP(Constants(3)/Temp)*H2O)*EXP(Constants(4)/Temp)
    END SUBROUTINE Spec4MCMCompute
    !
    !
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
    !
    !
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
    !
    !
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
    !
    !
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
    !
    !
    SUBROUTINE S4H2OCompute(ReacConst,Constants,Temp,mAir)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: mAir
      !
      ReacConst=Constants(1)*EXP(Constants(2)/Temp)+      & 
      &         mAir*Constants(3)*EXP(Constants(4)/Temp)
    END SUBROUTINE S4H2OCompute
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
      REAL(Realkind), INTENT(OUT) :: Sun
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
    SUBROUTINE DConstCompute(EquiRate,BackRate,Constants)
      REAL(RealKind) :: EquiRate
      REAL(RealKind) :: BackRate
      REAL(RealKind) :: Constants(:)
      !
      EquiRate=Constants(1)
      BackRate=Constants(2)
    END SUBROUTINE DConstCompute
    !
    !
    SUBROUTINE DTempCompute(EquiRate,BackRate,Constants,Temp)
      REAL(RealKind) :: EquiRate
      REAL(RealKind) :: BackRate
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      EquiRate=Constants(1)*EXP(Constants(2)*(ONE/Temp-InvRefTemp))
      BackRate=Constants(3)
    END SUBROUTINE DTempCompute
    !
    !
    SUBROUTINE DTemp2Compute(EquiRate,BackRate,Constants,Temp)  
      REAL(RealKind) :: EquiRate
      REAL(RealKind) :: BackRate
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      EquiRate=Constants(1)*EXP(Constants(2)*(ONE/Temp-InvRefTemp))
      BackRate=Constants(3)*EXP(Constants(4)*(ONE/Temp-InvRefTemp))
    END SUBROUTINE DTemp2Compute
    !
    !
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
    !
    !
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
    !
    !
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
    !
    !
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
    !
    !
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
    !
    !
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
    !
    !
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
    !
    !
    SUBROUTINE T1H2OCompute(ReacConst,Constants,Temp)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: Temp
      !
      ReacConst=Constants(1)*EXP(-Constants(2)/Temp)
    END SUBROUTINE T1H2OCompute
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
  END MODULE Rates_Mod
