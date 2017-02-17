MODULE String_Mod

  IMPLICIT NONE
  INTERFACE OPERATOR(.NES.)
    MODULE PROCEDURE CompStringNE
  END INTERFACE
  INTERFACE OPERATOR(.EQS.)
    MODULE PROCEDURE CompStringEQ
  END INTERFACE


CONTAINS

FUNCTION CompStringNE(S1,S2)
  LOGICAL :: CompStringNE
  CHARACTER, INTENT(IN), TARGET :: S1(:),S2(:)

  INTEGER :: i
  CHARACTER, POINTER :: SLong(:),SShort(:)

  IF (SIZE(S1)>=SIZE(S2)) THEN
    SLong=>S1
    SShort=>S2
  ELSE   
    SLong=>S1
    SShort=>S2
  END IF  
  CompStringNE=.FALSE.
  DO i=1,SIZE(SLong)
    IF (i<=SIZE(SShort)) THEN
      IF (SLong(i)/=SShort(i)) THEN
        CompStringNE=.TRUE.
        EXIT
      END IF
    ELSE
      IF (Slong(i)/=' ') THEN
        CompStringNE=.TRUE.
        EXIT
      END IF
    END IF
  END DO

END FUNCTION CompStringNE  

FUNCTION CompStringEQ(S1,S2)
  LOGICAL :: CompStringEQ
  CHARACTER, INTENT(IN), TARGET :: S1(:),S2(:)

  INTEGER :: i
  CHARACTER, POINTER :: SLong(:),SShort(:)

  IF (SIZE(S1)>=SIZE(S2)) THEN
    SLong=>S1
    SShort=>S2
  ELSE   
    SLong=>S1
    SShort=>S2
  END IF  
  CompStringEQ=.TRUE.
  DO i=1,SIZE(SLong)
    IF (i<=SIZE(SShort)) THEN
      IF (SLong(i)/=SShort(i)) THEN
        CompStringEQ=.FALSE.
        EXIT
      END IF
    ELSE
      IF (Slong(i)/=' ') THEN
        CompStringEQ=.FALSE.
        EXIT
      END IF
    END IF
  END DO

END FUNCTION CompStringEQ  

END MODULE String_Mod

