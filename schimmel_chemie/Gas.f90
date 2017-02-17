!**************************************************************************!
!***   
!***  Computation of Sinks and Sources from the Gas Phase Chemistry
!***   
!**************************************************************************!
!
SUBROUTINE Gas(nc,yg,fg,t)
  USE mo_reac
  USE mo_met
  USE mo_control

  IMPLICIT NONE 

  INTEGER ::  nc
  REAL(8) ::  yg(nc,ntgas)
  REAL(8) ::  fg(nc,ntgas)
  REAL(8) ::  Rate(nc,ntFrac) 
  REAL(8) ::  t

  TYPE (reaction), POINTER :: current

  REAL(8) :: w(nc)
  INTEGER :: l
  
  REAL(8) :: SumRate=0.0d0

!--------------------------------------------------------------------------!
  INTERFACE
    SUBROUTINE Rates(nc,nsp,iaF,ieF,yc,current,Rate)
       USE mo_reac
       USE mo_met
       USE mo_control
!   
       INTEGER :: nc, nsp, iaF,ieF
       REAL(8) :: yc(nc,nsp,ieF)
       REAL(8) :: Rate(nc,ieF)

       TYPE (reaction), POINTER :: current
    END SUBROUTINE rates
  END INTERFACE
!
!--------------------------------------------------------------------------!
!---  set pointer CURRENT
!
  IF(chi <= PiHalf .AND.ASSOCIATED(gphoto_first)) THEN
     current=>gphoto_first
  ELSE IF(ASSOCIATED(gconst_first)) THEN
     current=>gconst_first
  ELSE IF(ASSOCIATED(gas_first)) THEN
     current=>gas_first
  ELSE
     NULLIFY(current)
  END IF


!--------------------------------------------------------------------------!
!---  Begin Main Loop
!--------------------------------------------------------------------------!
!
  IF (ASSOCIATED(current)) THEN   
     DO
       !---  computation of reaction rates

        CALL Rates(nc,ntGas,1,1,yg,current,Rate)

!--------------------------------------------------------------------------!
!--- Checking for additional factors
!
        IF (ASSOCIATED(current%factor)) THEN  
           SELECT CASE (current%factor)

              CASE ('$O2')               
                    Rate(:,1) = Rate(:,1)*((mO2*mair(:))**current%fac_exp)*current%fac_A
                    
              CASE ('$H2O')                 
                    Rate(:,1) = Rate(:,1)*(mH2O(:)**current%fac_exp)*current%fac_A
     
              CASE ('$N2')                 
                    Rate(:,1) = Rate(:,1)*((mN2*mair(:))**current%fac_exp)*current%fac_A

              CASE ('$H2')                 
                    Rate(:,1) = Rate(:,1)*((mH2*mair(:))**current%fac_exp)*current%fac_A

              CASE ('$RO2')
                    Rate(:,1) = Rate(:,1) * SumRO2(:)

              CASE ('$M')
                    Rate(:,1) = Rate(:,1)*(mair(:)**current%fac_exp)*current%fac_A

              CASE ('$O2O2')
                    Rate(:,1) = Rate(:,1)*(((mO2*mair(:))**2)**current%fac_exp)*current%fac_A

              CASE ('$O2N2')
                    Rate(:,1) = Rate(:,1)*(((mO2*mair(:))*(mN2*mair(:)))**current%fac_exp)*current%fac_A

              CASE DEFAULT
                   WRITE(*,*) 'GasPhase: WRONG FACTOR FOR PASSIVE SPECIES ',current%factor
                   STOP       'GasPhase: WRONG FACTOR FOR PASSIVE SPECIES'

           END SELECT
        ENDIF

        w(:) = Rate(:,1)                 ! rate constant

        DO l=current%n_so_a+1,current%n_so            ! passive part               
           w(:) = w(:) * ykat(:,current%so(l))
        END DO


        DO l=1,current%n_so_a              
!          WRITE(*,*) 'WW',w(:),yg(:,current%so(l))
           w(:) = w(:) * yg(:,current%so(l))			!original
           
        ! ==============TEST====================
        !   w(:) = w(:) * 1.0d0
        ! ======================================
        END DO 
        
        ! ==============TEST====================
        !SumRate=SumRate+w(1)
        ! ======================================
		
        DO l=1,current%n_siso
           fg(:,current%reactant(l)%i_spc) = fg(:,current%reactant(l)%i_spc) &
&               + current%reactant(l)%d_koef*w(:)
        END DO

        current%w(:,1) = w(:)                 ! save for diagnose 
        
        ! ==============TEST====================
		!	WRITE(*,*) current%str_class, current%str_type, 'w',current%w(:,1)
		! ======================================

        IF (.NOT.ASSOCIATED(current%next)) EXIT
        current => current%next

     END DO
  END IF
     ! ======================================
	 ! 			     WRITE(*,*) 'STOP in Gas.f90',SumRate
	 !						STOP
	 ! ======================================
  		

!--------------------------------------------------------------------------!
END SUBROUTINE Gas

