MODULE Functions_Mod

  USE Kind_Mod
  
  IMPLICIT NONE
  
!  INTERFACE ERROR
!    REAL FUNCTION ERROR_MAX(ynew,yhat,y0,ATol,RTol)
!      REAL(8), DIMENSION(:), INTENT(IN) :: ynew, yhat, y0, ATol
!      REAL(8), INTENT(IN) :: RTol
!    END FUNCTION
!    REAL FUNCTION ERROR_NORM2(ynew,yhat,y0,ATol,RTol,Threshold,h)
!      REAL(8), DIMENSION(:), INTENT(IN) :: ynew, yhat, y0, ATol, Threshold
!      REAL(8), INTENT(IN) :: RTol, h
!    END FUNCTION
!  END INTERFACE ERROR
  
  CONTAINS
  
  REAL FUNCTION MAXNORM(Vec)
    REAL(RealKind) :: Vec(:)
    MAXNORM=MAXVAL(ABS(Vec))  
  END FUNCTION
  
  REAL FUNCTION ERROR_MAXNORM(ynew,yhat,y0,ATol,RTol)
    REAL(RealKind), DIMENSION(:) :: ynew, yhat, y0, ATol
    REAL(RealKind) :: RTol
    ERROR_MAXNORM=MAXNORM((ynew-yhat)/(ATol+MAX(ABS(y0),ABS(ynew))*RTol))  
  END FUNCTION
  
  REAL FUNCTION ERROR_NORM2(ynew,yhat,y0,ATol,RTol,Threshold,h)
    REAL(RealKind), DIMENSION(:) :: ynew, yhat, y0, ATol, Threshold
    REAL(RealKind) :: RTol, h
    ERROR_NORM2=(h/6.0d0)*NORM2((ynew-yhat)/MAX(MAX(ABS(y0),ABS(ynew)),(Threshold)))
  END FUNCTION  
  
END MODULE Functions_Mod