!**************************************************************************!
!***
!***  Computation of pH-Value (Hp concentration) in (iC,iFrac)
!***
!**************************************************************************!
!
   REAL(8) FUNCTION pHValue(ya)
      USE mo_reac
      USE mo_control

      IMPLICIT NONE

      REAL(8) :: ya(ntAqua)   ! molar density

      INTEGER :: jt

!--------------------------------------------------------------------------!
!
      pHValue = 0.d0
      DO jt=1,ntAqua!-1 
        !print*, 'jt,name,charge', jt,y_name(ntGas+jt), Ladung(jt)
        IF (jt==Hp_ind) CYCLE
        pHValue = pHValue + ya(jt) * Ladung(jt)
      END DO
      pHValue = -pHValue

!--------------------------------------------------------------------------!
   END FUNCTION pHValue

