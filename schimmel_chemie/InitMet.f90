!=================================================================
!  subroutine initmet() 
!=================================================================
!
  SUBROUTINE InitMet(Time) 
     USE mo_met
     USE mo_control

     IMPLICIT NONE

     REAL(8) :: Time

     INTEGER ::  i,il,it, iFrac
     REAL(8) ::  expo

    INTEGER, PARAMETER :: ntap11 = 11   &
&                        ,klines = 8

!==================================================
!===  Initialize Meteorology Data
!==================================================
!
!---  array allocation
     ALLOCATE(meteo(3))
     DO it=1,3
        ALLOCATE(meteo(it)%temp(nCell))
        ALLOCATE(meteo(it)%pres(nCell))
        ALLOCATE(meteo(it)%dens(nCell))
        ALLOCATE(meteo(it)%rhum(nCell))
        ALLOCATE(meteo(it)%ch2o(nCell))
     !
        ALLOCATE(meteo(it)%u(nCell))
        ALLOCATE(meteo(it)%w(nCell))
        ALLOCATE(meteo(it)%kz(nCell))
        ALLOCATE(meteo(it)%centrain(nCell))
     !
     !---  set zero
        meteo(it)%temp(:) = 0.E0
        meteo(it)%pres(:) = 0.E0
        meteo(it)%dens(:) = 0.E0
        meteo(it)%rhum(:) = 0.E0
        meteo(it)%ch2o(:) = 0.E0
     !
        meteo(it)%u(:)  = 0.E0
        meteo(it)%w(:)  = 0.E0
        meteo(it)%kz(:) = 0.E0
        meteo(it)%centrain(:) = 0.E0
     END DO

!---  set pointer
     u  => meteo(1)%u(:)
     ww => meteo(1)%w(:)
     kz => meteo(1)%kz(:)
     centrain => meteo(1)%centrain(:)

     temp   => meteo(2)%temp(:)
     pres   => meteo(2)%pres(:)
     densi  => meteo(2)%dens(:)
     rhum   => meteo(2)%rhum(:)
     ch2o   => meteo(2)%ch2o(:)

     ALLOCATE(mair(nCell))
     ALLOCATE(mH2O(nCell))

!======================================================
!===  Set Meteorology Pointer
!======================================================
!
     met_cur => meteo(1)
     met_new => meteo(2)
     met_old => meteo(3)
!
     met_index_new = 2
     met_index_old = 3
     met_time_new  = tmet
     met_time_old  = -1.d0
!
     FlagAqua = 0           ! liquid water switching flag 
     met_new%ModAqua = 1    ! aqueous phase on
!
!--------------------------------------------------
  END SUBROUTINE InitMet
