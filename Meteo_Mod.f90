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
  REAL(dp), PARAMETER :: mAir = 2.46d19/Temp*298.15d0*Pres/p0*1.0d0

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

CONTAINS 

    
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
