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
    CONTAINS
    !
    !
    !======================================================================!
    !      Calculate the Rates for current concentraions and time
    !======================================================================!
    SUBROUTINE Rates(Time,y_conc,Rate,DRatedT)
      !--------------------------------------------------------------------!
      ! Input: 
      !   - Time
      !   - y_conc
!     REAL(RealKind), INTENT(IN) :: Time
      REAL(RealKind) :: Time
      REAL(RealKind), INTENT(IN) :: y_conc(nDIM)
      !--------------------------------------------------------------------!
      ! Output:
      !   - Rate vector
      REAL(RealKind) :: Rate(neq)
      REAL(RealKind) :: DRatedT(neq)
      !--------------------------------------------------------------------!
      ! Temporary variables:
      REAL(RealKind) :: chi, actLWC
      REAL(RealKind) :: T(7)
      REAL(RealKind) :: k    ! tmp buffer for reaction const
      REAL(RealKind) :: DkdT    ! tmp buffer for reaction const
      REAL(RealKind) :: Meff
      INTEGER :: iReac, ii, j
      ! 
      !
      !==================================================================!
      !===============calc rate for ReactionSystem=======================!
      !
      TimeRateA=MPI_WTIME()
      !
      Rate=ZERO
      DRatedT=ZERO
      ! --- Compute zenith for photo reactions
      chi=Zenith(Time)
      ! 
      ! --- Compute the liquid water content 
      actLWC=pseudoLWC(Time)
      !
      !print*, 'debug:: chi,lwc', chi, actLWC
      ! --- Update temperature array
      IF ( combustion ) THEN
        CALL UpdateTempArray(T,y_conc(nDIM))
      ELSE
        T(1)=Temp    ! = 280 [K]
      END IF
      !
      ! --- compute rate of all reactions (gas,henry,aqua,diss)
      !print*, MPI_ID, '--loc_rateCnt= ',loc_rateCnt
      !print*, MPI_ID, '--loc_ratePtr= ',loc_ratePtr(:)
      DO ii=1,loc_rateCnt
        IF (MPI_np>1.AND.ParOrdering>=0) THEN
          iReac=loc_RatePtr(ii)
        ELSE
          iReac=ii
        END IF
        !
        ! ====== Compute the rate constant for specific reaction type ===
        CALL ComputeRateConstant(k,DkdT,T,Time,chi,mAir,iReac,y_conc,Meff)
        !
        ! ================ Multiplication with FACTOR====================
        CALL CheckThirdBodys(Meff,y_conc,iReac)
        !
        ! ========= Multiplication with Inactiv (passiv) Educts ===========
        CALL CheckInactEducts(k,iReac,actLWC)
        !
        ! ===== mult with (LWC^SUM(EductsCoefs)^-1) for correct dim units =============
        IF (ReactionSystem(iReac)%Type=='DISS'.OR.ReactionSystem(iReac)%Type=='AQUA') THEN
          CALL AquaDimLWC(k,iReac,actLWC)
        END IF
        !
        ! === Calc the coefficient for Gas->Aqua or Aqua->Gas Reactions ===
        IF (ReactionSystem(iReac)%Type=='HENRY') THEN
          CALL HenryMassTransferCoef(k,T(1),Time,iReac,actLWC)
        END IF
        !
        ! ======================== Build rate =============================
        DO j=A%RowPtr(iReac),A%RowPtr(iReac+1)-1
          IF ( A%Val(j)==ONE ) THEN
            k=k*y_conc(A%ColInd(j))
          ELSE IF ( A%Val(j)==TWO ) THEN
            k=k*y_conc(A%ColInd(j))*ABS(y_conc(A%ColInd(j)))
          ELSE
            k=k*y_conc(A%ColInd(j))**A%Val(j)
          END IF  
        END DO
        ! Last step
        Rate(iReac) = Meff * k
       ! print*, 'debugg:: i,rate(i) = ', ireac, Rate(ireac)
        IF (combustion) DRatedT(iReac) = DkdT
      END DO
      !stop
      TimeRateE=MPI_WTIME()
      TimeRates=TimeRates+(TimeRateE-TimeRateA)
      !
      
      ! gather the values of the other processes
      !CALL GatherAllPartitions(Rate,MyParties)
       !print*, 'MPI_ID=',MPI_ID, 'time=',Time
      !IF (combustion) CALL GatherAllPartitions(DRatedT,MyParties)
    END SUBROUTINE Rates
   !
    !
    !=====================================================================!
    ! === Converts the mass for Henry reactions Gas->Aqua , Aqua->Gas
    !=====================================================================!
    SUBROUTINE HenryMassTransferCoef(k,Temp,Time,iReac,actLWC)
      REAL(RealKind) :: k
      REAL(RealKind) :: Temp, Time, actLWC
      INTEGER :: iReac
      !
      REAL(RealKind) :: term_diff
      REAL(RealKind) :: term_accom
      REAL(RealKind) :: kmt
      !
      REAL(RealKind), PARAMETER :: drittel=1.0d0/3.0d0
      !
      !---------------------------------------------------------------------------
      term_diff =y_c1(ReactionSystem(iReac)%HenrySpc)                       ! diffusion term
      term_accom=y_c2(ReactionSystem(iReac)%HenrySpc)/SQRT(Temp)            ! accom term
      !--------------------------------------------------------------------------!
      !
      ! Compute new wet radius for droplett class iFrac
      SPEK(1)%wetRadius=(3.0d0/4.0d0/PI*actLWC/SPEK(1)%Number)**(drittel)
      !
      !--  mass transfer coefficient
      IF (term_diff.ne.0.0d0)  THEN   
        kmt = 1.0d0/(term_diff*SPEK(1)%wetRadius*SPEK(1)%wetRadius &
        &                           + term_accom*SPEK(1)%wetRadius  )
      ELSE
        kmt = dkmt
      END IF
      !
      ! direaction GasSpecies-->AquaSpecies
      IF (ReactionSystem(iReac)%direction=='GA') THEN  
        k=kmt*actLWC*1.d-3 ! orginal
        !
      ! direaction AquaSpecies-->GasSpecies  
      ELSE 
        k=kmt/(k*GasConst_R*Temp)   !()=HenryConst*GasConstant*Temperatur
      END IF
    END SUBROUTINE HenryMassTransferCoef
    !
    !=======================================================================!
    ! Compute the correct unit dimension for higher order aqueous reactions 
    !=======================================================================!
    SUBROUTINE AquaDimLWC(k,iReac,actLWC)
      INTEGER, INTENT(IN) :: iReac
      REAL(RealKind), INTENT(IN) :: actLWC
      REAL(RealKind), INTENT(INOUT) :: k
      !
      k=k/((actLWC*mol2part)**ReactionSystem(iReac)%SumAqCoef)  ! SumAqCoef=SUM(educt%koeff)-1
      !
    END SUBROUTINE AquaDimLWC
    !
    !=======================================================================!
    ! ===  Select the Type of the Constant and calculate the value
    !=======================================================================!
    SUBROUTINE ComputeRateConstant(k,DkdT,T,Time,chi,mAir,iReac,y_conc,Meff)
      REAL(RealKind), INTENT(INOUT) :: k, DkdT
      REAL(RealKind), INTENT(IN) :: Time, mAir, chi
      REAL(RealKind), INTENT(IN) :: T(7)
      REAL(RealKind), INTENT(IN) :: y_conc(:)
      INTEGER, INTENT(IN) :: iReac
      REAL(RealKind), INTENT(IN), OPTIONAL :: Meff
      REAL(RealKind) :: EqRate,BaRate,FoRate
      !
      REAL(RealKind) :: k0, kinf, Fcent, logF
      ! calc reaction constant
      ! Skip photochemical reactions at night
      BaRate=ZERO
      DkdT=ZERO
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
          CALL TempXComputeCK(k,EqRate,FoRate,ReactionSystem(iReac)%Constants,T,iReac)
          CALL DiffTempXComputeCK(DkdT,EqRate,FoRate,ReactionSystem(iReac)%Constants,T,iReac)
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
        CASE ('PRESSX')
          CALL PressXCompute(k,k0,kinf,EqRate,Fcent,logF,iReac,T,mAir,Meff)
          CALL DiffPressXCompute(DkdT,k0,kinf,iReac,T,Meff,logF,Fcent)
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
    SUBROUTINE ReactionDirection(ReacConst,EquiRate,BackRate,iReac)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: EquiRate, Backrate
      INTEGER :: iReac
      !
      ! backward reaction: k_b = k_b
      IF (ReactionSystem(iReac)%Line2=='BackReaction') THEN
        ReacConst=BackRate
      !
      ! forwar reaction: k_f = k_eq * k_b
      ELSE
        ReacConst=EquiRate*BackRate
      END IF
    END SUBROUTINE ReactionDirection   
    !
    !=====================================================================!
    ! ==  Multiplication with passiv species
    !========================================================passiv species
    SUBROUTINE CheckInactEducts(k,iReac,actLWC)
      REAL(RealKind), INTENT(INOUT) :: k
      INTEGER, INTENT(IN) :: iReac
      REAL(RealKind), INTENT(IN) :: actLWC
      !
      IF (ReactionSystem(iReac)%nInActEd/=0) THEN
        SELECT CASE (ReactionSystem(iReac)%InActEductSpc(1))
          CASE ('[H2O]')
            k=k*H2O
          CASE ('[N2]')
            k=k*N2
          CASE ('[O2]')
            k=k*O2
          CASE ('[aH2O]')
            k=k*aH2OmolperL*actLWC*mol2part
          CASE DEFAULT
            !STOP
        END SELECT
      END IF
    END SUBROUTINE CheckInactEducts
    !
    !=====================================================================!
    ! ===  Multiplication with FACTOR
    !=====================================================================!
    SUBROUTINE CheckThirdBodys(M,y_conc,iReac)
      REAL(RealKind), INTENT(OUT) :: M
      REAL(RealKind), INTENT(IN)    :: y_conc(:)
      INTEGER, INTENT(IN) :: iReac
      INTEGER :: i
      !
      !print*, 'debugg ro2=',SUM(y_conc(RO2))
      !DO i=1,1201
      !  WRITE(334,*) RO2(i)
      !END DO
      !stop
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
            M=SUM(y_conc(RO2))
          CASE ('$O2O2')
            M=(((mO2*mair)**2)**fac_exp)*fac_A
          CASE ('$aH2O')
            !k=k*aH2OmolperL*actLWC*mol2part
          CASE ('$RO2aq')
            M=SUM(y_conc(RO2aq))
          CASE ('$+M')
            CALL ThirdBodyCompute(M,iReac,y_conc)
          CASE ('$(+M)')
            CALL ThirdBodyCompute(M,iReac,y_conc)
          CASE DEFAULT
            WRITE(*,*) 'Reaction: ',iReac
            CALL FinishMPI()
            STOP 'Unknown FACTOR (error at Rate calc)'
        END SELECT
      ELSE 
        M=1.0d0
      END IF
      !
    END SUBROUTINE CheckThirdBodys
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
    SUBROUTINE TempXComputeCK(ReacConst,EquiRate,ForwRate,Constants,T,iReac)
      REAL(RealKind) :: ReacConst
      REAL(RealKind) :: ForwRate, EquiRate!, BackRate
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: T(7)
      INTEGER :: iReac
      !
      REAL(RealKind) :: g0(nspc)   
      REAL(RealKind) :: delg0
      INTEGER :: jj,j
      !
      IF (ReactionSystem(iReac)%Factor/='rev') THEN
        !
        IF (ReactionSystem(iReac)%Line3=='rev'.OR.ReactionSystem(iReac)%Line3=='REV') THEN
          ! backward with backwardcoeffs reaction
          ReacConst=Constants(1)*T(1)**Constants(2)*EXP(-Constants(3)*T(6))
        ELSE
          ! backward calculated with eq. constant
          ! Standard non-dimensional Gibbs free energy
          !print*, 'debug:: switchtemp = ', SwitchTemp
          WHERE (SwitchTemp>T(1))
            g0=-(lowA*(LOG(T(1))-1) + 0.5d0*lowB*T(1) + lowC*T(2)/6.0d0 + &
            &        lowD*T(3)/12.0d0 + 0.05d0*lowE*T(4) + lowF*T(6) + lowG )           ! [-] dim.less
          ELSEWHERE
            g0=-(highA*(LOG(T(1))-1) + 0.5d0*highB*T(1) + highC*T(2)/6.0d0 + &
            &        highD*T(3)/12.0d0 + 0.05d0*highE*T(4) + highF*T(6) + highG )           ! [-] dim.less
          END WHERE
          !
          delg0=ZERO
          DO jj=BA%RowPtr(iReac),BA%RowPtr(iReac+1)-1
            j=BA%ColInd(jj)
            delg0=delg0+BA%Val(jj)*g0(j)      
          END DO
          !
          EquiRate=EXP(-delg0)*(PressR*T(6))**sumBAT(iReac)
          ! backward reaction
          ReacConst=Constants(1)*T(1)**Constants(2)*EXP(-Constants(3)*T(6))/EquiRate
        END IF
      ELSE
        ! forward reaction
        ReacConst=Constants(1)*T(1)**Constants(2)*EXP(-Constants(3)*T(6))
      END IF
    END SUBROUTINE TempXComputeCK
    !
    !
    SUBROUTINE DiffTempXComputeCK(DkdT,EquiRate,ForwRate,Constants,T,iReac)
      REAL(RealKind) :: DkDT
      REAL(RealKind) :: DfRdT, DbRdT, DeRdT
      REAL(RealKind) :: ForwRate, EquiRAte
      REAL(RealKind) :: Constants(:)
      REAL(RealKind) :: T(7)              ! temperatur array
      INTEGER :: iReac
      !
      REAL(RealKind) :: Dg0dT(nspc)       ! in [1/K]
      REAL(RealKind) :: Ddelg0dT
      INTEGER :: jj, i
      REAL(RealKind) :: Tmp
      !
      !hier ?????????????
      !???????
      !
      WHERE (SwitchTemp>T(1))
        Dg0dT=-(lowA*T(6) + 0.5d0*lowB + lowC*T(1)/3.0d0 +     &
        &              0.25d0*lowD*T(2) + 0.2d0*lowE*T(3) + lowF*T(7) )
      ELSEWHERE
        Dg0dT=-(highA*T(6) + 0.5d0*highB + highC*T(1)/3.0d0 +     &
        &              0.25d0*highD*T(2) + 0.2d0*highE*T(3) + highF*T(7) )
      END WHERE  
      !
      Tmp=0.0d0
      DO jj=A%RowPtr(iReac),A%RowPtr(iReac+1)-1
        Tmp=Tmp-Dg0dT(A%ColInd(jj))*A%Val(jj)
      END DO
      Ddelg0dT=Tmp
      Tmp=0.0d0
      DO jj=B%RowPtr(iReac),B%RowPtr(iReac+1)-1
        Tmp=Tmp+Dg0dT(B%ColInd(jj))*B%Val(jj)
      END DO
      Ddelg0dT=Ddelg0dT+Tmp
      !
      DfRdT=ForwRate*T(6)*(Constants(2)+Constants(3)*T(6))
      !
      DeRdT=-EquiRate*(sumBAT(iReac)*T(6)+Ddelg0dT)
      !
      DbRdT=(DfRdT-ForwRate*DeRdT/EquiRate)/EquiRate
      !
      IF (ReactionSystem(iReac)%Line3=='rev') THEN
        DkDT=DbRdT
      ELSE   
        DkDT=DfRdT
      END IF
      !
    END SUBROUTINE DiffTempXComputeCK
    ! 
    !
    SUBROUTINE ThirdBodyCompute(M,i,y_conc)
      REAL(RealKind) :: M
      REAL(RealKind) :: y_conc(nspc)
      INTEGER :: i        !iReac
      !
      INTEGER :: j, jj, alphI
      !
      M=0
      alphI=1
      DO jj=1,SIZE(ReactionSystem(i)%TB)
        j=ReactionSystem(i)%TB(jj)
        M=M+y_conc(j)*(1.0d0-ReactionSystem(i)%TBalpha(jj) )
      END DO
    END SUBROUTINE ThirdBodyCompute
    !
    !
    SUBROUTINE UpdateTempArray(TempArr,Temp)
      REAL(RealKind) :: Temp
      REAL(RealKind) :: TempArr(7)
      !
      INTEGER :: i
      !
      TempArr(1)=Temp                     ! T
      DO i=2,5
        TempArr(i)=TempArr(i-1)*Temp      ! T^2 ... T^5
      END DO
      TempArr(6)=1.0d0/Temp               ! 1/T
      TempArr(7)=1.0d0/TempArr(2)         ! 1/T^2
    END SUBROUTINE UpdateTempArray
    !
    !
    SUBROUTINE SpcInternalEnergy(Umol,T)
      REAL(RealKind) :: Umol(nspc)        ! Species internal energy in [J/mol]
      REAL(RealKind) :: T(7)
      !  
      !
      WHERE (SwitchTemp>T(1))
        Umol=((lowA-1.0d0)*T(1) + 0.5d0*lowB*T(2) + lowC*T(3)/3.0d0 & 
        &                + 0.25d0*lowD*T(4) + 0.2d0*lowE*T(5) + lowF)
      ELSEWHERE
        Umol=((highA-1.0d0)*T(1) + 0.5d0*highB*T(2) + highC*T(3)/3.0d0 & 
        &                + 0.25d0*highD*T(4) + 0.2d0*highE*T(5) + highF)
      END WHERE
      Umol=Umol*R_const
    END SUBROUTINE SpcInternalEnergy
    !
    !
    SUBROUTINE DiffSpcInternalEnergy(DUmoldT,T)
      REAL(RealKind) :: DUmoldT(nspc)     ! derivative of Species internal energyin [J/mol/K], (constant volume specific heat)
      REAL(RealKind) :: T(7)                      
      !
      WHERE (SwitchTemp>T(1))
        DUmoldT=((lowA-1.0d0) + lowB*T(1) + lowC*T(2) & 
        &                + lowD*T(3) + lowE*T(4) )
      ELSEWHERE
        DUmoldT=((highA-1.0d0) + highB*T(1) + highC*T(2) & 
        &                + highD*T(3) + highE*T(4) )
      END WHERE
      DUmoldT=DUmoldT*R_const
    END SUBROUTINE DiffSpcInternalEnergy
    !
    !
    SUBROUTINE Diff2SpcInternalEnergy(D2UmoldT2,T)
      REAL(RealKind) :: D2UmoldT2(nspc)     !Constant volume specific heat’s derivative [J/mol/K2] 
      REAL(RealKind) :: T(7)                      
      !
      WHERE (SwitchTemp>T(1))
        D2UmoldT2=(lowB + 2.0d0*lowC*T(1) & 
        &                + 3.0d0*lowD*T(2) + 4.0d0*lowE*T(3) )
      ELSEWHERE    
        D2UmoldT2=(highB + 2.0d0*highC*T(1) & 
        &                + 3.0d0*highD*T(2) + 4.0d0*highE*T(3) )
      END WHERE
      D2UmoldT2=D2UmoldT2*R_const
    END SUBROUTINE Diff2SpcInternalEnergy
    !
    !is the mass average mixture specific  heat at constant volume,
    SUBROUTINE MassAveMixSpecHeat(c_v,T,y,DUmoldT)
      REAL(RealKind) :: c_v
      !
      REAL(RealKind) :: T(:)
      REAL(RealKind) :: y(:)
      REAL(RealKind) :: DUmoldT(:)
      !
      INTEGER :: i
      !
      c_v=0.0d0
      DO i=1,nspc
        c_v=c_v+DUmoldT(i)*y(i) !/W(i) für andere einheit
      END DO
    END SUBROUTINE MassAveMixSpecHeat
    !
    !==========================================================================!
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
    SUBROUTINE PressXCompute(ReacConst,k0,kinf,EquiRate,Fcent,logF,i,T,mAir,Meff)
      REAL(RealKind) :: ReacConst, EquiRate
      !
      REAL(RealKind) :: HighConst(3), LowConst(3)
      REAL(RealKind) :: T(:)
      REAL(RealKind) :: mAir, Meff
      INTEGER :: i, jj, j
      !
      REAL(RealKind) :: g0(nspc)
      REAL(RealKind) :: delg0,  k
      REAL(RealKind) :: k0 ,kinf, Pr, logF, c, n, Fcent
      REAL(RealKind) :: TroeConst(4)
      !
      IF (ReactionSystem(i)%Line3=='low') THEN
        HighConst(:)=ReactionSystem(i)%Constants(:)
        LowConst(:)=ReactionSystem(i)%LowConst(:)
      ELSE
        HighConst(:)=ReactionSystem(i)%HighConst(:)
        LowConst(:)=ReactionSystem(i)%Constants(:)
      END IF
      !
      TroeConst(:)=ReactionSystem(i)%TroeConst(:)
      !
      k0=LowConst(1)*T(1)**LowConst(2)*EXP(-LowConst(3)/R_Const*T(6))
      kinf=HighConst(1)*T(1)**HighConst(2)*EXP(-HighConst(3)/R_Const*T(6))
      !
      Pr=k0/kinf*Meff
      !
      ! hier nochmal nach LIND oder TROE abfragen für H2 O2 nur troe
      ! IF (TROE) THEN
      Fcent=(1.0d0-TroeConst(1))*EXP(-T(1)/TroeConst(2))+                   &
      &      TroeConst(1)*EXP(-T(1)/TroeConst(3))+EXP(-TroeConst(4)*T(6))
      c=-0.4d0-0.67d0*LOG10(Fcent)
      n=0.75d0-1.27d0*LOG10(Fcent)
      logF=LOG10(Fcent)/(1.0d0 + ((LOG10(Pr)+c)/(n-0.14d0*(LOG10(Pr)+c)))**2)
      k=kinf*(Pr/(1.0d0+Pr))*10.0d0**logF
      !
      IF (ReactionSystem(i)%Line2=='BackReaction') THEN
        ! backward calculated with eq. constant
        ! Standard non-dimensional Gibbs free energy
        WHERE (SwitchTemp>T(1))
          g0=-(lowA*(LOG(T(1))-1) + 0.5d0*lowB*T(1) + lowC*T(2)/6.0d0 + &
          &        lowD*T(3)/12.0d0 + 0.05d0*lowE*T(4) + lowF*T(6) + lowG )           ! [-] dim.less
        ELSEWHERE
          g0=-(highA*(LOG(T(1))-1) + 0.5d0*highB*T(1) + highC*T(2)/6.0d0 + &
          &        highD*T(3)/12.0d0 + 0.05d0*highE*T(4) + highF*T(6) + highG )           ! [-] dim.less
        END WHERE
        !
        delg0=ZERO
        DO jj=BA%RowPtr(i),BA%RowPtr(i+1)-1
          j=BA%ColInd(jj)
          delg0=delg0+BA%Val(jj)*g0(j)      
        END DO
        !
        EquiRate=EXP(-delg0)*(PressR*T(1))**sumBAT(i)
        ! backward reaction
        ReacConst=k/EquiRate
      ELSE
        EquiRate=ONE
        ReacConst=k
      END IF
    END SUBROUTINE PressXCompute
    !
    !
    SUBROUTINE DiffPressXCompute(DiffCoef,kf0,kfoo,i,T,Meff,logF_troe,F_cent)
      REAL(RealKind) :: DiffCoef
      REAL(RealKind) :: logF_troe, kf0, kfoo, Meff, F_cent
      REAL(RealKind) :: Dkf0dT, DkfoodT
      REAL(RealKind) :: T(7)
      INTEGER :: i 
      REAL(RealKind) :: DF_troedT, DF_linddT
      REAL(RealKind) :: Pr, DlogFTdPr, DlogFTdT, DlogPrdT
      REAL(RealKind) :: c, n ,d, l10Prc
      !
      REAL(RealKind) :: HighConst(3), LowConst(3)
      !
      IF (ReactionSystem(i)%Line3=='low') THEN
        HighConst(:)=ReactionSystem(i)%Constants(:)
        LowConst(:)=ReactionSystem(i)%LowConst(:)
      ELSE
        HighConst(:)=ReactionSystem(i)%HighConst(:)
        LowConst(:)=ReactionSystem(i)%Constants(:)
      END IF
      !
      Pr=kf0/kfoo*Meff
      c=-0.4d0-0.67d0*LOG10(F_cent)
      n=0.75d0-1.27d0*LOG10(F_cent)
      d=0.14d0
      l10Prc=LOG10(Pr)+c
      !
      Dkf0dT=kf0*(LowConst(2)+LowConst(3)/R_Const*T(6))*T(6)
      DkfoodT=kfoo*(HighConst(2)+HighConst(3)/R_Const*T(6))*T(6)
      !
      DlogFTdPr= -2.0d0*LOG10(F_cent)*l10Prc                  &
      &          *(n-d*l10Prc)*(n*n-2.0d0*n*d*l10Prc          &
      &                         +(d+1)*l10Prc*l10Prc)**(-2.0d0)
      !
      DF_linddT= Meff/(kfoo+kf0*Meff)**2 * ( kfoo*Dkf0dT - kf0*DkfoodT )
      !
      DlogPrdT=1.0d0/LOG10(10.0d0)*(1.0d0/kf0*Dkf0dT+1.0d0/kfoo*DkfoodT)
      DlogFTdT=DlogFTdPr*DlogPrdT
      !
      DF_troedT=10.0d0**logF_troe * (DF_linddT+LOG10(10.0d0)*DlogFTdT)
      !
      DiffCoef=DF_troeDT+10.0d0**logF_troe*(HighConst(2)+HighConst(3)/R_Const*T(6))*T(6)
    END SUBROUTINE
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
