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
