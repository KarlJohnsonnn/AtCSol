 MODULE mo_mphys

!--  micophysics structures
   TYPE microphys
      INTEGER,POINTER :: mpAnf(:,:),mpEnd(:,:)
      INTEGER,POINTER :: mpStat(:,:)

      REAL(8) :: time
      REAL(8), POINTER :: Number(:,:)
      REAL(8), POINTER :: Mass(:,:)
      REAL(8), POINTER :: LwcSum(:)
      REAL(8), POINTER :: Lwc(:,:)
      REAL(8), POINTER :: cradius(:,:)
      REAL(8), POINTER :: DryMass(:,:)
      REAL(8), POINTER :: LwcFlux(:,:,:)
      REAL(8), POINTER :: ActFlux(:,:,:)
   END TYPE microphys

!--  time-dependent microphysics
   TYPE (microphys), POINTER:: mp_cur, mp_old, mp_new
   TYPE (microphys), ALLOCATABLE, TARGET:: mphys(:)

   INTEGER :: mp_index_old = 2,  mp_index_new = 3
   REAL(8) :: mp_time_old, mp_time_new, mp_time_next

!--  temporary pointers and variables
   INTEGER, POINTER :: mpAnf(:,:),mpEnd(:,:)
   INTEGER, POINTER :: mpStat(:,:)

   REAL(8), POINTER :: Number(:,:)
   REAL(8), POINTER :: Mass(:,:)
   REAL(8), POINTER :: LwcSum(:)
 !  REAL(8), POINTER :: Lwc(:,:)
   REAL(8), POINTER :: cradius(:,:)
   REAL(8), POINTER :: DryMass(:,:)
   REAL(8), POINTER :: LwcFlux(:,:,:)
   REAL(8), POINTER :: ActFlux(:,:,:)

!--  set boundaries of fractions
   REAL(8), POINTER :: BdRadius(:)

!--  number of soluted moles per particle/droplet
   REAL(8), POINTER :: Solute(:,:)

!--  number of ions per kg-Water
   REAL(8), POINTER :: NbIons(:,:)

!--  arrays for saving "fixed aerosol" concentrations for species
!--  which are not included in microphysics
   REAL(8), POINTER :: AeroConc(:)
 
!--  flags und control parameters
   INTEGER :: mpInput = 0       &
&            ,mpUnit  = 0       &
&           , mpFlag  = 0

!--  control parameters
   INTEGER :: FracScav = 1       &   ! minimum fraction for gas uptake
&            ,FracChem = 1           ! minimum fraction for chemistry
   REAL(8) :: Mass2Vol               ! water unit transfer  [kg/kg] ==> [l/m3]

 END MODULE mo_mphys
