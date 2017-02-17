!**************************************************************************!
!***
!***  Computation of Sinks and Sources from the Gas Phase Chemistry
!***
!**************************************************************************!
!
  SUBROUTINE Henry(nc,iaF,ieF,yg,fg,ca,fa,t)
     USE mo_reac
     USE mo_met
     USE mo_mphys
     USE mo_control

     IMPLICIT NONE 

     INTEGER :: nc, iaF, ieF
     REAL(8) :: t
     REAL(8) :: yg(nc,ntgas), ca(nc,ntaqua,ieF) 
     REAL(8) :: fg(nc,ntgas), fa(nc,ntaqua,ieF) 
     REAL(8) :: Shade(nc), Rate(nc,ntFrac)

!--------------------------------------------------------------------------!
!---  internal varibles
!
     TYPE (reaction), POINTER :: current
     INTEGER :: iFrac, ihen

     REAL(8) :: ConvFac, term_diff, nue_otemp
     REAL(8) :: term_accom(nC)
     REAL(8) :: kmt(nc), w(nc), conc_gas(nc)

!--------------------------------------------------------------------------!
     INTERFACE
       SUBROUTINE rates(nc,nsp,iaF,ieF,yc,current,Rate)
          USE mo_reac
          USE mo_met
          USE mo_control
     
          INTEGER :: nc, nsp, iaF, ieF
          REAL(8) :: yc(nc,nsp,ieF)
          REAL(8) :: Rate(nc,ieF), Shade(nc)

          TYPE (reaction), POINTER :: current
       END SUBROUTINE rates
     END INTERFACE

!--------------------------------------------------------------------------!
!---  set pointer CURRENT
!
     IF ((ieF <= 0) .OR. (ntaqua <= 0)) RETURN

     IF(ASSOCIATED(henry_first)) THEN
        current => henry_first
     ELSE
        NULLIFY(current)
     END IF

!==========================================================================!
!===  Begin Main Loop
!==========================================================================!
!
     IF (ASSOCIATED(current))  THEN

!--  conversion of GasUnit into  [mol/l]
        ConvFac = ConvGas * 1.e-03

!--------------------------------------------------------------------------!
!---  Mass transfer between reactive species
!--------------------------------------------------------------------------!
!          
           
        DO ihen=1,nHenry 

!--- computation of reaction rates
           CALL Rates(nc,ntAqua,iaF,ieF,ca,current,Rate)

           current%w(:,:) = 0.E0
           current%v(:,:) = 0.E0

           conc_gas(:) = yg(:,current%so(1))    ! gas phase concentration

           IF (current%anz_p > 2) THEN
              term_diff   = 1/(3.E0*current%dparam(current%anz_p-1))
              nue_oTemp    = SQRT(8.E0 * 8.313E0 / Pi1 / current%dparam(current%anz_p-2)) ! RS: new version
              term_accom(:) = 4.E0 / (3.E0*current%dparam(current%anz_p)*nue_oTemp) / SQRT(temp(:))
           ELSE
              term_diff   = y_c1(current%so(1))    ! gas diffusion term
              term_accom(:) = y_c2(current%so(1)) / SQRT(temp(:))
           ENDIF

! 
!--------------------------------------------------------------------------!
           DO iFrac=iaF,ieF

!--  mass transfer coefficient
              IF (term_diff .ne. 0.0e0)  THEN   
                 kmt(:) = 1.0e0/(term_diff*cradius(:,iFrac)*cradius(:,iFrac) &
&                       + term_accom(:)*cradius(:,iFrac))
              ELSE
                 kmt(:) = dkmt
              END IF

!--  equilibrium term  ( Gas * ConvFac - Aqua/(RT*HenryConst))
              w(:) =  Rate(:,iFrac) * GasConst_R * temp(:)
              w(:) =  conc_gas(:) * ConvFac - ca(:,current%si(1),iFrac)/w(:)

!--  function values
              w(:) =  kmt(:) * w(:)

              fg(:,current%so(1)) = fg(:,current%so(1)) - 1.E-03 * lwc(:,iFrac)*w(:) / ConvFac
              current%w(:,iFrac)    = -  lwc(:,iFrac) * w(:) 

              fa(:,current%si(1),iFrac) = fa(:,current%si(1),iFrac) + w(:)
              current%v(:,iFrac)        = lwc(:,iFrac) * w(:)

           END DO

           IF (.NOT.ASSOCIATED(current%next)) EXIT
           current=>current%next
        END DO
        STOP 'initchem'
!--------------------------------------------------------------------------!
!---  Mass transfer for non-reactive species
!--------------------------------------------------------------------------!
! 
        DO ihen=nHenry+1,nreakhenry
           current%w(:,:) = 0.E0
           current%v(:,:) = 0.E0

           

           IF (ind_henry(ihen,1) == 0)  THEN
              IF (ASSOCIATED(current%factor)) THEN
                SELECT CASE (current%factor)

                 CASE ('$O2')                 
                        conc_gas(1:nc) = ((mO2*mair(1:nc))**current%fac_exp)*current%fac_A
                      
                 CASE ('$N2')                 
                        conc_gas(1:nc) = ((mN2*mair(1:nc))**current%fac_exp)*current%fac_A

                 CASE ('$H2')                 
                        conc_gas(1:nc) = ((mH2*mair(1:nc))**current%fac_exp)*current%fac_A

                 CASE DEFAULT
                    WRITE(*,*) 'HENRY: WRONG FACTOR FOR PASSIVE SPECIES'
                    STOP       'HENRY: WRONG FACTOR FOR PASSIVE SPECIES'
             
                END SELECT
              ELSE
                conc_gas(:)   = ykat(:,current%so(1))             ! see below
              ENDIF

              IF (current%anz_p <= 2) THEN
                 term_diff     = y_c1(ntGas+ntAqua+current%so(1))
                 term_accom(:) = y_c2(ntGas+ntAqua+current%so(1)) / SQRT(temp(:))
              ENDIF
   
           ELSE IF (ind_henry(ihen,2) == 0)  THEN
              conc_gas(:) = yg(:,current%so(1))    ! gas phase concentration

              IF (current%anz_p <= 2) THEN
                 term_diff   = y_c1(current%so(1))    ! gas diffusion term
                 term_accom(:) = y_c2(current%so(1)) / SQRT(temp(:))
              ENDIF

           ELSE
              STOP 'HENRY: Wrong index IND_HENRY for passive species!!'
           END IF

           IF (current%anz_p > 2) THEN
              term_diff   = 1/(3.E0*current%dparam(current%anz_p-1))    ! gas diffusion term

              nue_oTemp    = SQRT(8.E0 * 8.313E0 / Pi1 / current%dparam(current%anz_p-2))
              term_accom(:) = 4.E0 / (3.E0*current%dparam(current%anz_p)*nue_oTemp) / SQRT(temp(:))
           ENDIF

! 
!--- computation of reaction rates
           CALL Rates(nc,ntAqua,iaF,ieF,ca,current,Rate)

           DO iFrac=iaF,ieF

!--  mass transfer coefficient
              IF (term_diff .ne. 0.0e0)  THEN   
                 kmt(:) = 1.0e0/(term_diff*cradius(:,iFrac)*cradius(:,iFrac) &
&                       + term_accom(:)*cradius(:,iFrac))
              ELSE
                 kmt(:) = dkmt
              END IF

!--  equilibrium term  ( Gas * ConvFac - Aqua/(RT*HenryConst))
              w(:) =  Rate(:,iFrac) * GasConst_R * temp(:)
              w(:) =  conc_gas(:) * ConvFac - ca(:,current%si(1),iFrac)/w(:)

!--  function values
              w(:) =  kmt(:) * w(:)

              current%w(:,iFrac) = - lwc(:,iFrac) * w(:) 
              current%v(:,iFrac) =   lwc(:,iFrac) * w(:)

              IF (ind_henry(ihen,1) == 0)  THEN         
                 fa(:,current%si(1),iFrac) = fa(:,current%si(1),iFrac) + w(:)
              ELSE 
                 fg(:,current%so(1)) = fg(:,current%so(1)) - 1.E-03 * lwc(:,iFrac)*w(:) / ConvFac
              END IF
           END DO

           IF (.NOT.ASSOCIATED(current%next)) EXIT
           current=>current%next
        END DO

     END IF

!--------------------------------------------------------------------------!
!
  END SUBROUTINE Henry

