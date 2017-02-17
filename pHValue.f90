!**************************************************************************!
!***
!***  Computation of pH-Value (Hp concentration) in (iC,iFrac)
!***
!**************************************************************************!
!
   FUNCTION pHValue(ya)
     USE mo_reac
     USE mo_control
     USE Kind_Mod

     IMPLICIT NONE

     REAL(RealKind) :: pHValue
     REAL(RealKind) :: ya(ntAqua)   ! molar density
     INTEGER :: jt

!--------------------------------------------------------------------------!
!
     pHValue = 0.d0
     DO jt=1,ntAqua!-1 
       IF (jt==Hp_ind) CYCLE
       pHValue = pHValue + ya(jt) * Charge(jt)
     END DO
     pHValue = -pHValue

!--------------------------------------------------------------------------!
   END FUNCTION pHValue

