MODULE mo_unirnk
INTEGER, PARAMETER :: kdp = SELECTed_REAL_KIND(15)
PUBLIC :: unirnk
PRIVATE :: kdp
PRIVATE :: R_unirnk, I_unirnk, D_unirnk
PRIVATE :: R_nearless, I_nearless, D_nearless, nearless
INTERFACE unirnk
  !MODULE PROCEDURE D_unirnk,  I_unirnk
  MODULE PROCEDURE D_unirnk, R_unirnk, I_unirnk
END INTERFACE unirnk
INTERFACE nearless
  MODULE PROCEDURE D_nearless, R_nearless, I_nearless
  !MODULE PROCEDURE D_nearless,  I_nearless
END INTERFACE nearless

contains

SUBROUTINE D_unirnk (XVALT, IRNGT, NUNI)
! __________________________________________________________
!   UNIRNK = Merge-sort ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it discards indices that correspond to
!   duplicate entries.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      REAL (KIND=8), DIMENSION (:), INTENT (In) :: XVALT  ! zu sortierENDer Vektor
      INTEGER, DIMENSION (:), INTENT (Out) :: IRNGT         ! Permutationsvektor
      INTEGER, INTENT (Out) :: NUNI                         ! laenge des sortierten Vektor
! __________________________________________________________
      INTEGER, DIMENSION (SIZE(IRNGT)) :: JWRKT
      INTEGER :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
      INTEGER :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
      REAL (KIND=8) :: XTST, XVALA, XVALB
!
!
      NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
      NUNI = NVAL
!
      SELECT CASE (NVAL)
      CASE (:0)
         RETURN
      CASE (1)
         IRNGT (1) = 1
         RETURN
      CASE DEFAULT
         CONTINUE
      END SELECT
!
!  Fill-in the index array, creating ordered couples
!
      DO IIND = 2, NVAL, 2
         IF (XVALT(IIND-1) < XVALT(IIND)) THEN
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         ELSE
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         END IF
      END DO
      IF (MODULO(NVAL, 2) /= 0) THEN
         IRNGT (NVAL) = NVAL
      END IF
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      DO
         IF (NVAL <= 4) EXIT
!
!   Loop on merges of A and B into C
!
         DO IWRKD = 0, NVAL - 1, 4
            IF ((IWRKD+4) > NVAL) THEN
               IF ((IWRKD+2) >= NVAL) EXIT
!
!   1 2 3
!
               IF (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) EXIT
!
!   1 3 2
!
               IF (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) THEN
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               ELSE
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               END IF
               EXIT
            END IF
!
!   1 2 3 4
!
            IF (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) CYCLE
!
!   1 3 x x
!
            IF (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) THEN
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               IF (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) THEN
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               ELSE
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               END IF
!
!   3 x x x
!
            ELSE
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               IF (XVALT(IRNG1) <= XVALT(IRNGT(IWRKD+4))) THEN
                  IRNGT (IWRKD+2) = IRNG1
                  IF (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) THEN
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  ELSE
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  END IF
               ELSE
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               END IF
            END IF
         END DO
!
!  The Cs become As and Bs
!
         LMTNA = 4
         EXIT
      END DO
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is DOubled.
!
      DO
         IF (2*LMTNA >= NVAL) EXIT
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         DO
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            IF (IWRKF >= NVAL) THEN
               IF (JINDA >= NVAL) EXIT
               IWRKF = NVAL
            END IF
            IINDA = 1
            IINDB = JINDA + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
            XVALA = XVALT (JWRKT(IINDA))
            XVALB = XVALT (IRNGT(IINDB))
!
            DO
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               IF (XVALA > XVALB) THEN
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  IF (IINDB > IWRKF) THEN
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     EXIT
                  END IF
                  XVALB = XVALT (IRNGT(IINDB))
               ELSE
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  IF (IINDA > LMTNA) EXIT! Only B still with unprocessed values
                  XVALA = XVALT (JWRKT(IINDA))
               END IF
!
            END DO
         END DO
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      END DO
!
!   Last merge of A and B into C, with removal of duplicates.
!
      IINDA = 1
      IINDB = LMTNA + 1
      NUNI = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
      IF (IINDB <= NVAL) THEN
        XTST = NEARLESS (Min(XVALT(JWRKT(1)), XVALT(IRNGT(IINDB))))
      ELSE
        XTST = NEARLESS (XVALT(JWRKT(1)))
      ENDIF
      DO IWRK = 1, NVAL
!
!  We still have unprocessed values in both A and B
!
         IF (IINDA <= LMTNA) THEN
            IF (IINDB <= NVAL) THEN
               IF (XVALT(JWRKT(IINDA)) > XVALT(IRNGT(IINDB))) THEN
                  IRNG = IRNGT (IINDB)
                  IINDB = IINDB + 1
               ELSE
                  IRNG = JWRKT (IINDA)
                  IINDA = IINDA + 1
               END IF
            ELSE
!
!  Only A still with unprocessed values
!
               IRNG = JWRKT (IINDA)
               IINDA = IINDA + 1
            END IF
         ELSE
!
!  Only B still with unprocessed values
!
            IRNG = IRNGT (IWRK)
         END IF
         IF (XVALT(IRNG) > XTST) THEN
            XTST = XVALT (IRNG)
            NUNI = NUNI + 1
            IRNGT (NUNI) = IRNG
         END IF
!
      END DO
!
      RETURN
!
END SUBROUTINE D_unirnk

SUBROUTINE R_unirnk (XVALT, IRNGT, NUNI)
! __________________________________________________________
!   UNIRNK = Merge-sort ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it discards indices that correspond to
!   duplicate entries.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      REAL(KIND=4), DIMENSION (:), INTENT (In) :: XVALT
      INTEGER, DIMENSION (:), INTENT (Out) :: IRNGT
      INTEGER, INTENT (Out) :: NUNI
! __________________________________________________________
      INTEGER, DIMENSION (SIZE(IRNGT)) :: JWRKT
      INTEGER :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
      INTEGER :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
      REAL(KIND=4) :: XTST, XVALA, XVALB
!
!
      NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
      NUNI = NVAL
!
      SELECT CASE (NVAL)
      CASE (:0)
         RETURN
      CASE (1)
         IRNGT (1) = 1
         RETURN
      CASE DEFAULT
         CONTINUE
      END SELECT
!
!  Fill-in the index array, creating ordered couples
!
      DO IIND = 2, NVAL, 2
         IF (XVALT(IIND-1) < XVALT(IIND)) THEN
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         ELSE
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         END IF
      END DO
      IF (MODULO(NVAL, 2) /= 0) THEN
         IRNGT (NVAL) = NVAL
      END IF
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      DO
         IF (NVAL <= 4) EXIT
!
!   Loop on merges of A and B into C
!
         DO IWRKD = 0, NVAL - 1, 4
            IF ((IWRKD+4) > NVAL) THEN
               IF ((IWRKD+2) >= NVAL) EXIT
!
!   1 2 3
!
               IF (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) EXIT
!
!   1 3 2
!
               IF (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) THEN
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               ELSE
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               END IF
               EXIT
            END IF
!
!   1 2 3 4
!
            IF (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) CYCLE
!
!   1 3 x x
!
            IF (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) THEN
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               IF (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) THEN
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               ELSE
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               END IF
!
!   3 x x x
!
            ELSE
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               IF (XVALT(IRNG1) <= XVALT(IRNGT(IWRKD+4))) THEN
                  IRNGT (IWRKD+2) = IRNG1
                  IF (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) THEN
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  ELSE
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  END IF
               ELSE
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               END IF
            END IF
         END DO
!
!  The Cs become As and Bs
!
         LMTNA = 4
         EXIT
      END DO
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is DOubled.
!
      DO
         IF (2*LMTNA >= NVAL) EXIT
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         DO
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            IF (IWRKF >= NVAL) THEN
               IF (JINDA >= NVAL) EXIT
               IWRKF = NVAL
            END IF
            IINDA = 1
            IINDB = JINDA + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
            XVALA = XVALT (JWRKT(IINDA))
            XVALB = XVALT (IRNGT(IINDB))
!
            DO
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               IF (XVALA > XVALB) THEN
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  IF (IINDB > IWRKF) THEN
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     EXIT
                  END IF
                  XVALB = XVALT (IRNGT(IINDB))
               ELSE
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  IF (IINDA > LMTNA) EXIT! Only B still with unprocessed values
                  XVALA = XVALT (JWRKT(IINDA))
               END IF
!
            END DO
         END DO
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      END DO
!
!   Last merge of A and B into C, with removal of duplicates.
!
      IINDA = 1
      IINDB = LMTNA + 1
      NUNI = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
      IF (IINDB <= NVAL) THEN
        XTST = NEARLESS (Min(XVALT(JWRKT(1)), XVALT(IRNGT(IINDB))))
      ELSE
        XTST = NEARLESS (XVALT(JWRKT(1)))
      ENDIF
      DO IWRK = 1, NVAL
!
!  We still have unprocessed values in both A and B
!
         IF (IINDA <= LMTNA) THEN
            IF (IINDB <= NVAL) THEN
               IF (XVALT(JWRKT(IINDA)) > XVALT(IRNGT(IINDB))) THEN
                  IRNG = IRNGT (IINDB)
                  IINDB = IINDB + 1
               ELSE
                  IRNG = JWRKT (IINDA)
                  IINDA = IINDA + 1
               END IF
            ELSE
!
!  Only A still with unprocessed values
!
               IRNG = JWRKT (IINDA)
               IINDA = IINDA + 1
            END IF
         ELSE
!
!  Only B still with unprocessed values
!
            IRNG = IRNGT (IWRK)
         END IF
         IF (XVALT(IRNG) > XTST) THEN
            XTST = XVALT (IRNG)
            NUNI = NUNI + 1
            IRNGT (NUNI) = IRNG
         END IF
!
      END DO
!
      RETURN
!
END SUBROUTINE R_unirnk
SUBROUTINE I_unirnk (XVALT, IRNGT, NUNI)
! __________________________________________________________
!   UNIRNK = Merge-sort ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it discards indices that correspond to
!   duplicate entries.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      INTEGER, DIMENSION (:), INTENT (In) :: XVALT
      INTEGER, DIMENSION (:), INTENT (Out) :: IRNGT
      INTEGER, INTENT (Out) :: NUNI
! __________________________________________________________
      INTEGER, DIMENSION (SIZE(IRNGT)) :: JWRKT
      INTEGER :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
      INTEGER :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
      INTEGER :: XTST, XVALA, XVALB
!
!
      NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
      NUNI = NVAL

      SELECT CASE (NVAL)
      CASE (:0)
         RETURN
      CASE (1)
         IRNGT (1) = 1
         RETURN
      CASE DEFAULT
         CONTINUE
      END SELECT
!
!  Fill-in the index array, creating ordered couples
!
      DO IIND = 2, NVAL, 2
         IF (XVALT(IIND-1) < XVALT(IIND)) THEN
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         ELSE
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         END IF
      END DO
      IF (MODULO(NVAL, 2) /= 0) THEN
         IRNGT (NVAL) = NVAL
      END IF
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      DO
         IF (NVAL <= 4) EXIT
!
!   Loop on merges of A and B into C
!
         DO IWRKD = 0, NVAL - 1, 4
            IF ((IWRKD+4) > NVAL) THEN
               IF ((IWRKD+2) >= NVAL) EXIT
!
!   1 2 3
!
               IF (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) EXIT
!
!   1 3 2
!
               IF (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) THEN
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               ELSE
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               END IF
               EXIT
            END IF
!
!   1 2 3 4
!
            IF (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) CYCLE
!
!   1 3 x x
!
            IF (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) THEN
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               IF (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) THEN
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               ELSE
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               END IF
!
!   3 x x x
!
            ELSE
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               IF (XVALT(IRNG1) <= XVALT(IRNGT(IWRKD+4))) THEN
                  IRNGT (IWRKD+2) = IRNG1
                  IF (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) THEN
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  ELSE
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  END IF
               ELSE
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               END IF
            END IF
         END DO
!
!  The Cs become As and Bs
!
         LMTNA = 4
         EXIT
      END DO
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is DOubled.
!
      DO
         IF (2*LMTNA >= NVAL) EXIT
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         DO
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            IF (IWRKF >= NVAL) THEN
               IF (JINDA >= NVAL) EXIT
               IWRKF = NVAL
            END IF
            IINDA = 1
            IINDB = JINDA + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
            XVALA = XVALT (JWRKT(IINDA))
            XVALB = XVALT (IRNGT(IINDB))
!
            DO
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               IF (XVALA > XVALB) THEN
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  IF (IINDB > IWRKF) THEN
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     EXIT
                  END IF
                  XVALB = XVALT (IRNGT(IINDB))
               ELSE
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  IF (IINDA > LMTNA) EXIT! Only B still with unprocessed values
                  XVALA = XVALT (JWRKT(IINDA))
               END IF
!
            END DO
         END DO
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      END DO
!
!   Last merge of A and B into C, with removal of duplicates.
!
      IINDA = 1
      IINDB = LMTNA + 1
      NUNI = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
      IF (IINDB <= NVAL) THEN
        XTST = NEARLESS (Min(XVALT(JWRKT(1)), XVALT(IRNGT(IINDB))))
      ELSE
        XTST = NEARLESS (XVALT(JWRKT(1)))
      ENDIF
      DO IWRK = 1, NVAL
!
!  We still have unprocessed values in both A and B
!
         IF (IINDA <= LMTNA) THEN
            IF (IINDB <= NVAL) THEN
               IF (XVALT(JWRKT(IINDA)) > XVALT(IRNGT(IINDB))) THEN
                  IRNG = IRNGT (IINDB)
                  IINDB = IINDB + 1
               ELSE
                  IRNG = JWRKT (IINDA)
                  IINDA = IINDA + 1
               END IF
            ELSE
!
!  Only A still with unprocessed values
!
               IRNG = JWRKT (IINDA)
               IINDA = IINDA + 1
            END IF
         ELSE
!
!  Only B still with unprocessed values
!
            IRNG = IRNGT (IWRK)
         END IF
         IF (XVALT(IRNG) > XTST) THEN
            XTST = XVALT (IRNG)
            NUNI = NUNI + 1
            IRNGT (NUNI) = IRNG
         END IF
!
      END DO
!
      RETURN
!
END SUBROUTINE I_unirnk

FUNCTION D_nearless (XVAL) RESULT (D_nl)
!  Nearest value less than given value
! __________________________________________________________
      REAL (KIND=8), INTENT (In) :: XVAL
      REAL (KIND=8) :: D_nl
! __________________________________________________________
      D_nl = nearest (XVAL, -1.0_kdp)
      RETURN
!
END FUNCTION D_nearless
FUNCTION R_nearless (XVAL) RESULT (R_nl)
!  Nearest value less than given value
! __________________________________________________________
      REAL(KIND=4), INTENT (In) :: XVAL
      REAL(KIND=4) :: R_nl
! __________________________________________________________
      R_nl = nearest (XVAL, -1.0)
      RETURN
!
END FUNCTION R_nearless
FUNCTION I_nearless (XVAL) RESULT (I_nl)
!  Nearest value less than given value
! __________________________________________________________
      INTEGER, INTENT (In) :: XVAL
      INTEGER :: I_nl
! __________________________________________________________
      I_nl = XVAL - 1
      RETURN
!
END FUNCTION I_nearless

END MODULE mo_unirnk
