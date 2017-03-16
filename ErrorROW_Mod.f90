MODULE ErrorROW_Mod
  !
  USE Kind_Mod,   ONLY: RealKind
  USE mo_reac,    ONLY: nspc, nDIM
  USE mo_control, ONLY: Error_Est
  !
  IMPLICIT NONE
  !
  CONTAINS
  !
  FUNCTION MAXNORM(Vec)
    REAL(RealKind) :: MAXNORM
    REAL(RealKind) :: Vec(:)
    MAXNORM=MAXVAL(ABS(Vec))  
  END FUNCTION
  !
  FUNCTION NORM2(Vec)
    REAL(RealKind) :: NORM2
    REAL(RealKind) :: Vec(:)
    NORM2=SUM(Vec(:)*Vec(:))
  END FUNCTION
  
 SUBROUTINE ERROR(err,En_Index,ynew,yhat,ATol,RTol,t)
    !
    REAL(RealKind) :: err
    REAL(RealKind), DIMENSION(:) :: ynew, yhat, ATol
    REAL(RealKind) :: RTol, t
    !
    REAL(RealKind) :: scalTol(nDIM), En_Values(nDIM)
    INTEGER :: En_Index(1,1)
    !
    scalTol       = ATol + MAX( ABS(yhat) , ABS(ynew) ) * RTol ! scaling strategie
    En_Values     = ABS( ynew - yhat ) / scalTol               ! local error est.
    En_Index(1,1) = MAXLOC( En_Values , 1 )                    ! max error component
    !
    IF ( Error_Est == 2 ) THEN
      err = NORM2(En_Values) / nspc   ! euclikd norm
    ELSE
      err = MAXNORM(En_Values)      ! maximum norm
    END IF
    !
  END SUBROUTINE ERROR
END MODULE ErrorROW_Mod



