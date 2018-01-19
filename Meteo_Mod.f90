MODULE Meteo_Mod

  USE Kind_Mod

  IMPLICIT NONE
  
  REAL(dp), PARAMETER :: Temp = 280.0d0 !298.15d0
  REAL(dp), PARAMETER :: Cp   = 1004.0D0
  REAL(dp), PARAMETER :: Cv   = 717.0D0
  REAL(dp), PARAMETER :: RhoC = 1.14D0
  REAL(dp), PARAMETER :: Pres = 850.d0 ! hPa
  REAL(dp), PARAMETER :: p0   = 1013.25d0 ! hPa Normaldruck
  REAL(dp), PARAMETER :: Rd   = Cp-Cv
  !REAL(dp) :: mAir = 2.46d19/Temp*298.15d0*Pres/p0*1.0d0

  REAL(dp), PARAMETER :: mol2part   = 6.02295d17            &
  &                    , GasConst_R = 0.082056d0            &   ! [in l*atm/mol/K]
  &                    , InvRefTemp = 1.0D0/298.15D0        &
  &                    , RefTemp    = 298.15D0              &
  &                    , N2         = 1.960D19              &
  &                    , O2         = 5.100D18              &
  &                    , H2O        = 5.100D17              &
  &                    , M          = N2 + O2               &
  &                    , aH2OmolperL= 5.55D01               &
  &                    , aH2O       = aH2OmolperL*mol2part  &
  &                    , rhum       = 5.0920016900153622d-3 &
  &                    , mH2        = 0.0d0                 &
  &                    , mN2        = 0.7809d0              &
  &                    , mO2        = 0.2095d0              &
  &                    , MM_h2o     = 18.01534d0            &
  &                    , h2o_push   = 6.023d20/MM_h2o       &
  &                    , densi      = 1.0545184035426323d0 
  REAL(dp) :: mH2O

  REAL(dp), PARAMETER :: LWCconst=2.0d-8       ! [l/m3] no cloud 
  REAL(dp), PARAMETER :: NCC=1000d0

  !REAL(dp), PARAMETER :: R_const=8.344598d0 ! [ J / mol / K ] = [ kg * m2 / s2 / mol /K ]
  REAL(dp), PARAMETER :: R_const=8.3144621d0 ! [ J / mol / K ] = [ kg * m2 / s2 / mol /K ]
  REAL(dp) :: PressR=Pres/R_const

!-- more LWC stuff for pseudo function 
  REAL(dp) :: LWCb(6) ! boundaries for linear pseudo lwc function

  REAL(dp) :: LAT  = 45.0_dp
  REAL(dp) :: LONG =  0.0_dp
  INTEGER  :: IDAT = 010619

CONTAINS 

  !=========================================================================!
  !                  calculate sun 
  !=========================================================================!
  FUNCTION Zenith(Time) RESULT(Chi)
    USE mo_control
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

    
  FUNCTION Set_pseudoLWCbounds() RESULT(bounds)
    USE mo_control, ONLY: tBegin, HOUR
    REAL(dp) :: bounds(6)
    ! --- set cloud intervall
    bounds(1) = tBegin * HOUR 
    bounds(2) = bounds(1) + 1.00_dp  * HOUR
    bounds(3) = bounds(2) + 0.25_dp  * HOUR
    bounds(4) = bounds(3) + 9.50_dp  * HOUR
    bounds(5) = bounds(4) + 0.25_dp  * HOUR
    bounds(6) = bounds(5) + 1.00_dp  * HOUR
  END FUNCTION Set_pseudoLWCbounds

  SUBROUTINE ThirdPartyKoeff(Time)
    REAL(dp) :: Time

    mH2O=rhum*densi*h2o_push
    
  END SUBROUTINE ThirdPartyKoeff
    
  FUNCTION pseudoLWC(RealTime)  RESULT(LWC)
    USE Kind_Mod
    USE mo_reac
    USE mo_control
    !
    IMPLICIT NONE
    !
    REAL(dp) :: LWC
    REAL(dp) :: RealTime
    REAL(dp) :: Time
    !
    !         LWC
    !           /|\       fake LWC function (periodicly)
    !            |
    !            |
    ! maxlwclvl _|_ _ _ ____             ____
    !            |     /|  |\           /    \
    !            |    /      \         /      \
    !            |   /  |  |  \       /        \
    !            |  /          \     /          \
    ! minlwclvl _|_/    |  |    \___/            \__...... 
    !------------+----------------------------------------------->
    !            | |    |  |    | |                               Time
    !  LWCb 1 2    3  4    5 6
    !
    IF ( constLWC ) THEN
      ! Constant LWC level 
      LWC = LWCconst
    ELSE
      ! calculate LWC level
      Time = MODULO(RealTime,LWCb(6))
      !
      LWC = LwcLevelmin

      ! --- Constant minimum LWC level
      IF     (LWCb(1)<=Time.AND.Time<LWCb(2)) THEN
        LWC = LwcLevelmin

      ! --- linear increase of LWC level till maximum
      ELSEIF (LWCb(2)<=Time.AND.Time<LWCb(3)) THEN
        LWC = (Time-LWCb(2))/(LWCb(3)-LWCb(2))*LwcLevelmax  +  &
        &     (LWCb(3)-Time)/(LWCb(3)-LWCb(2))*LwcLevelmin

      ! --- Constant maximum LWC level
      ELSEIF (LWCb(3)<=Time.AND.Time<LWCb(4)) THEN
        LWC = LwcLevelmax

      ! --- linear decrease of LWC level till minimum
      ELSEIF (LWCb(4)<=Time.AND.Time<LWCb(5)) THEN
        LWC = (Time-LWCb(4))/(LWCb(5)-LWCb(4))*LwcLevelmin  +  &
        &     (LWCb(5)-Time)/(LWCb(5)-LWCb(4))*LwcLevelmax

      ! --- Constant minium LWC level
      ELSEIF (LWCb(5)<=Time.AND.Time<=LWCb(6)) THEN
        LWC = LwcLevelmin

      END IF
    END IF
  END FUNCTION

  !**************************************************************************!
  !***
  !***  Computation of pH-Value (Hp concentration)
  !***
  !**************************************************************************!
  !
  FUNCTION pHValue(ConcAqua) RESULT(pH)
    USE mo_reac,    ONLY: ns_AQUA, Hp_ind, Charge
    USE mo_control, ONLY: ZERO
    USE Kind_Mod,   ONLY: dp

    IMPLICIT NONE

    REAL(dp) :: pH
    REAL(dp) :: ConcAqua(ns_AQUA)   ! molar density
    INTEGER :: j

    pH = ZERO
    DO j=1,ns_AQUA!-1 
      IF (j/=Hp_ind) pH = pH + ConcAqua(j) * Charge(j)
    END DO
    pH = -pH
  END FUNCTION pHValue


END MODULE Meteo_Mod
