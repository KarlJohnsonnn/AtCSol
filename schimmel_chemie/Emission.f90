!**************************************************************************!
!***   
!***  External Emissions
!***   
!**************************************************************************!
!
SUBROUTINE Emission(nc,y,f,t)
  USE mo_reac
  USE mo_met
  USE mo_control

  IMPLICIT NONE 

  INTEGER :: nc
  REAL(8) :: t
  REAL(8) :: y(nc,nt) 
  REAL(8) :: f(nc,nt) 

  INTEGER :: ic, ia, i
  REAL(8) :: efa
  EXTERNAL   efa

!--------------------------------------------------------------------------!
!
  DO ic=1,nc
     ia = ntFrac*ntAqua
     DO i=1,ntGas
        f(ic,ia+i) = f(ic,ia+i) + efa(t) * y_e(ic,i)
     END DO

     ia = ntGas
     DO i=1,ntFrac*ntAqua
        f(ic,i) = f(ic,i) + efa(t) * y_e(ic,ia+i)
     END DO
  END DO

END SUBROUTINE Emission

