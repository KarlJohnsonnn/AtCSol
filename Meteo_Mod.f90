MODULE Meteo_Mod

  USE Kind_Mod

  IMPLICIT NONE
  
  REAL(RealKind), PARAMETER :: Temp=280.0d0 !298.15d0
  REAL(RealKind), PARAMETER :: Cp=1004.0D0
  REAL(RealKind), PARAMETER :: Cv=717.0D0
  REAL(RealKind), PARAMETER :: RhoC=1.14D0
  REAL(RealKind), PARAMETER :: Pres=850.d0 ! hPa
  REAL(RealKind), PARAMETER :: p0=1013.25d0 ! hPa Normaldruck
  REAL(RealKind), PARAMETER :: Rd=Cp-Cv
  REAL(RealKind), PARAMETER :: mAir=2.46d19/Temp*298.15d0*Pres/p0*1.0d0

  REAL(RealKind), PARAMETER :: InvRefTemp=1.0D0/298.15D0  &
  &                          , RefTemp=298.15D0           &
  &                          , N2=1.960D19                &
  &                          , O2=5.100D18                &
  &                          , H2O=5.100D17               &
  &                          , M=N2+O2                    &
  &                          , aH2OmolperL=5.55D01        &
  &                          , rhum=5.0920016900153622d-003 &
  &                          , mH2=0.0d0                  &
  &                          , mN2=0.7809d0               &
  &                          , mO2=0.2095d0               &
  &                          , MM_h2o=18.01534d0          &
  &                          , h2o_push=6.023d20/MM_h2o   &
  &                          , densi=1.0545184035426323d0 
   REAL(RealKind) :: mH2O

  REAL(RealKind), PARAMETER :: LWCconst=2.0d-8       ! [l/m3] no cloud 
  REAL(RealKind), PARAMETER :: NCC=1000d0

  !REAL(RealKind), PARAMETER :: R_const=8.344598d0 ! [ J / mol / K ] = [ kg * m2 / s2 / mol /K ]
  REAL(RealKind), PARAMETER :: R_const=8.3144621d0 ! [ J / mol / K ] = [ kg * m2 / s2 / mol /K ]
  REAL(RealKind) :: PressR=Pres/R_const

!-- more LWC stuff for pseudo function 
      REAL(RealKind) :: LWCbounds(6) ! boundaries for linear pseudo lwc function

CONTAINS 

SUBROUTINE ThirdPartyKoeff(Time)
  REAL(RealKind) :: Time

  mH2O=rhum*densi*h2o_push
  
END SUBROUTINE ThirdPartyKoeff
  
FUNCTION pseudoLWC(RealTime)  
  USE Kind_Mod
  USE mo_reac
  USE mo_control
  !
  IMPLICIT NONE
  !
  REAL(RealKind) :: pseudoLWC
  REAL(RealKind) :: RealTime
  REAL(RealKind) :: Time
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
  !  LWCbounds 1 2    3  4    5 6
  !
  IF (constLWC==0) THEN
    ! Constant LWC level 
    pseudoLWC=LWCconst
  ELSE
    ! calculate LWC level
    Time=MODULO(RealTime,LWCbounds(6))
    !
    pseudoLWC=LwcLevelmin
    ! --- Constant minimum LWC level
    IF     (LWCbounds(1)<=Time.AND.Time<LWCbounds(2)) THEN
      pseudoLWC=LwcLevelmin
    ! --- linear increase of LWC level till maximum
    ELSEIF (LWCbounds(2)<=Time.AND.Time<LWCbounds(3)) THEN
      pseudoLWC=(Time-LWCbounds(2))/(LWCbounds(3)-LWCbounds(2))*LwcLevelmax   +  &
      &         (LWCbounds(3)-Time)/(LWCbounds(3)-LWCbounds(2))*LwcLevelmin
    ! --- Constant maximum LWC level
    ELSEIF (LWCbounds(3)<=Time.AND.Time<LWCbounds(4)) THEN
      pseudoLWC=LwcLevelmax
    ! --- linear decrease of LWC level till minimum
    ELSEIF (LWCbounds(4)<=Time.AND.Time<LWCbounds(5)) THEN
      pseudoLWC=(Time-LWCbounds(4))/(LWCbounds(5)-LWCbounds(4))*LwcLevelmin   +  &
      &         (LWCbounds(5)-Time)/(LWCbounds(5)-LWCbounds(4))*LwcLevelmax
    ! --- Constant minium LWC level
    ELSEIF (LWCbounds(5)<=Time.AND.Time<=LWCbounds(6)) THEN
      pseudoLWC=LwcLevelmin
    END IF
  END IF
END FUNCTION

END MODULE Meteo_Mod
