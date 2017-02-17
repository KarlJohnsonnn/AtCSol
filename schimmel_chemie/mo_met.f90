 MODULE mo_met

!--  photolysis
   REAL(8) :: cosz,chi
   REAL(8) :: tmet

!--  meteorology
   TYPE meteorology
      REAL(8) :: time
      REAL(8) :: ModAqua
      REAL(8), POINTER :: u(:),w(:),kz(:),temp(:),pres(:)
      REAL(8), POINTER :: rhum(:),ch2o(:),dens(:),centrain(:)
   END TYPE meteorology

   TYPE (meteorology), POINTER:: met_cur, met_old, met_new
   TYPE (meteorology), ALLOCATABLE, TARGET:: meteo(:)

   REAL(8), POINTER :: u(:),ww(:),kz(:),temp(:),pres(:),densi(:)
   REAL(8), POINTER :: rhum(:),ch2o(:),centrain(:)
   REAL(8), ALLOCATABLE ::  mH2O(:)
   !REAL(8), ALLOCATABLE :: mair(:), mH2O(:)
   REAL(8), PARAMETER :: mO2  = 0.2095e0 ,& ! proportion of O2
&                        mN2  = 0.7809e0 ,& ! proportion of N2
&                        mH2  = 0.e0        ! proportion of H2 (not used as passive specie, yet)

!--  supersaturation
   REAL(8), ALLOCATABLE :: SuperSat(:)

   INTEGER :: MetInput = 0       &
&            ,MetHour  = 0       &
&            ,ModAqua  = 0       &
&            ,FlagAqua = 0
   INTEGER :: met_index_old = 2,  met_index_new = 3

   REAL(8) :: MetPeriod = 0.E0
   REAL(8) :: met_time_old, met_time_new, met_time_next

 END MODULE mo_met
