MODULE ErrorROW_Mod

  USE Kind_Mod
  USE mo_reac
  USE mo_control
  USE mo_MPI
  
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
  
 SUBROUTINE ERROR(err,En_Index,ynew,yhat,y0,ATol,RTol,t,errVals)
    !
    REAL(RealKind) :: err
    REAL(RealKind), DIMENSION(:) :: ynew, yhat, y0, ATol
    REAL(RealKind) :: RTol, t
    CHARACTER(60) :: errSpc
    REAL(RealKind), OPTIONAL :: errVals(nspc)
    !
    REAL(RealKind) :: scalTol(nspc), En_Values(nspc)
    INTEGER :: En_Index(1,1), yn,i
    !
    !do i=1,nspc
    !  if (i==24.OR.i==74) THEN
    !    scalTol(i)=1.d0-8+MAX(ABS(yhat(i)),ABS(ynew(i)))*RTol
    !  else
    !    scalTol(i)=ATol(i)+MAX(ABS(yhat(i)),ABS(ynew(i)))*RTol
    !  end if
    !end do
    scalTol(:)=ATol(:)+MAX(ABS(yhat(:)),ABS(ynew(:)))*RTol
    En_Values(:)=ABS(ynew(:)-yhat(:))/scalTol(:)
    En_Index(1,1)=MAXLOC(En_Values(:),1)
    !
    IF (Error_Est==2) THEN
      !err=NORM2(En_Values)/sqrtNSPC
      err=NORM2(En_Values)/nspc
    ELSE
      err=MAXNORM(En_Values)  
    END IF
    !
  END SUBROUTINE ERROR
END MODULE ErrorROW_Mod



