MODULE NetCDF_Mod
!--- Netcdf Output
  USE Kind_Mod
  USE Meteo_Mod
  USE netcdf
  USE mo_reac
  USE mo_control
  USE mo_MPI
  !
  IMPLICIT NONE
  !
  !
  INTEGER, ALLOCATABLE :: Diag_varID(:), DiagERR_varID(:)
  INTEGER, ALLOCATABLE :: OutNetcdfSpc(:)
  CHARACTER(1),   ALLOCATABLE :: OutNetcdfPhase(:) !g=gas , a=aqua
  INTEGER              :: OutNetcdfANZ, OutNetcdfDIM
  REAL(dp)       :: altit               ! altitude
  
  INTEGER :: iStpNetCDF 
  INTEGER :: itime_NetCDF
  REAL(dp), ALLOCATABLE :: yNcdf(:)     ! current output vector
  INTEGER :: ncid

  INTEGER :: pH_ind
  !
  INTEGER :: x_varid,y_varid,z_varid,rec_varid, LWC_varid, traj_varid
  INTEGER :: StepSize_varid, Gassum_varid, Aquasum_varid, wetRadius_varid
  INTEGER :: zenith_varid, schwefel_varid, error_varid, Temperature_varid
  INTEGER :: MaxErrorSpc_varid, pH_varid

  CONTAINS

  SUBROUTINE InitNetCDF(StartNetcdf, EndNetcdf)
  !
  !
  !==================================================================
  !===  Initialization of Netcdf Output File
  !==================================================================
  !
  !
  ! -- Input parameter --
  !
  !
  REAL(dp) :: StartNetcdf      &  ! First Time for Netcdf Output
  &         , EndNetcdf           ! End Time of Netcdf Output
  !
  ! -- Output parameter --
  !
  ! -- Temporary variables --
  INTEGER :: jt
  !
  !    ============================================================
  !    Variable Attribute Names
  !    ============================================================
  CHARACTER (LEN = *), PARAMETER :: REC_NAME     = "time"
  CHARACTER (LEN = *), PARAMETER :: REC_LONGNAME = "time"
  CHARACTER (LEN = *), PARAMETER :: LON_NAME     = "lon"
  CHARACTER (LEN = *), PARAMETER :: LAT_NAME     = "lat"
  CHARACTER (LEN = *), PARAMETER :: Z_NAME       = "altitude"
  CHARACTER (LEN = *), PARAMETER :: TRAJ_NAME    = "trajectory"
  CHARACTER (LEN = *), PARAMETER :: LON_LONGNAME = "longitude"
  CHARACTER (LEN = *), PARAMETER :: LAT_LONGNAME = "latitude"
  CHARACTER (LEN = *), PARAMETER :: Z_LONGNAME   = "height above mean sea level"
!    ============================================================
!    Variable Attribute Units
!    ============================================================
  CHARACTER(60), ALLOCATABLE :: DIAG_UNITS(:)
  CHARACTER (LEN = *), PARAMETER :: NC_UNITS     = "units"
  CHARACTER (LEN = *), PARAMETER :: LON_UNITS    = "degrees_east"
  CHARACTER (LEN = *), PARAMETER :: LAT_UNITS    = "degrees_north"
  CHARACTER (LEN = *), PARAMETER :: Z_UNITS      = "m"
  CHARACTER (LEN = *), PARAMETER :: REC_UNITS    = "hours"! since 2000-01-01 00:00:00"
  CHARACTER (LEN = *), PARAMETER :: DRYMASS_UNITS= "g/m3"
  CHARACTER (LEN = *), PARAMETER :: MOLAL_UNITS  = "g/l"
  !    ============================================================
  !    Dimension ID
  !    ============================================================
  INTEGER, PARAMETER :: NDIMS  = 1
  INTEGER            :: dimIDs(NDIMS)
  INTEGER            :: rec_dimid
  CHARACTER(20)      :: traj_id
  !    ============================================================
  INTEGER, PARAMETER :: name_strlen=10
  LOGICAL            :: FileExist
  CHARACTER(10)      :: answer
  CHARACTER(60)      :: tmpName
  INTEGER            :: iDiagSpc
  INTEGER            :: strich
  !

  ! -- Allocate Netcdf Names etc.
  !
  ALLOCATE(Diag_LongName(OutNetcdfANZ))
  ALLOCATE(Diag_Name_Netcdf(OutNetcdfANZ),Diag_UNITS(OutNetcdfANZ))
  ALLOCATE(DiagERR_Name_Netcdf(SIZE(OutNetcdfspc)))
  ALLOCATE(Diag_varID(OutNetcdfANZ))
  ALLOCATE(DiagERR_varID(OutNetcdfANZ))
  Diag_varID=0
  DiagERR_varID=0

  jt=0
  iDiagSpc=0

  DO iDiagSpc=1,SIZE(OutNetcdfSpc)
    jt=jt+1
    tmpName = ADJUSTL(y_name(OutNetcdfSpc(iDiagSpc)))
    !
    IF (OutNetcdfPhase(iDiagSpc)=='g') THEN
      Diag_Name_Netcdf(jt) = tmpName
      Diag_LongName(jt)    = tmpName
      DIAG_UNITS(jt)       = "mol/m3"
      IF (Teq) DIAG_UNITS(jt) = 'mass fraction'
    ELSEIF (OutNetcdfPhase(iDiagSpc)=='a') THEN
      Diag_Name_Netcdf(jt) = TRIM(tmpName)//'_l'
      Diag_LongName(jt)    = TRIM(tmpName)//'AQUA'
      DIAG_UNITS(jt)       = "mol/l"
      jt=jt+1
      Diag_Name_Netcdf(jt) = TRIM(tmpName)//'_m3'
      Diag_LongName(jt)    = TRIM(tmpName)//'AIR'
      DIAG_UNITS(jt)       = "mol/m3"
      IF (tmpName == 'Hp' ) pH_ind = jt
    END IF
    DO
      IF (INDEX(Diag_Name_Netcdf(jt),'/')>0) THEN
        CALL check_name2(Diag_Name_Netcdf(jt))
      ELSE
        EXIT
      END IF
    END DO
    IF (tmpName(1:1)=='[') THEN
      CALL check_name3(Diag_Name_Netcdf(jt))
    END IF
    IF (OutNetcdfPhase(iDiagSpc)=='a') THEN
      strich=INDEX(Diag_Name_Netcdf(jt),'_')
      tmpName=tmpName(:strich-1)
    END IF
    DiagERR_Name_Netcdf(iDiagSpc)=TRIM(tmpName)//'_ERR'
  END DO
  !
  ! ============================================================
  ! --  create Netcdf file (for each size bin / fraction)
  ! ============================================================
  ! Create new NetCDF File
  !
  CALL check(  NF90_CREATE(  TRIM(NetCDFFile) &   ! NetCDF output file name
  &                        , NF90_CLOBBER     &   ! Overwrite existing file with the same name
  &                        , ncid     )    )      ! Returned netCDF ID
  !
  ! ============================================================
  ! --  define dimension variables
  ! ============================================================
  ! Define the dimensions. The record dimension is defined to have
  ! unlimited length - it can grow as needed. In this example it is
  ! the time dimension.
  CALL check(  NF90_DEF_DIM(  ncid           & ! NetCDF ID, from prev call NF90_CREATE 
  &                         , REC_NAME          & ! Time is the record dimension
  &                         , NF90_UNLIMITED    & ! space for "unlimited" time steps
  &                         , rec_dimID     )   ) ! returned dimension ID
  !
  ! ============================================================
  ! --  define coordinate variables
  ! ============================================================
  ! Define the coordinate variables. We will only define coordinate
  ! variables for lat and lon.  Ordinarily we would need to provide
  ! an array of dimension IDs for each variable's dimensions, but
  ! since coordinate variables only have one dimension, we can
  ! simply provide the address of that dimension ID (lat_dimid) and
  ! similarly for (lon_dimid).
  !dimIDs = (/ rec_dimid, x_dimid, y_dimid, z_dimid /)
  !
  dimIDs = (/ rec_dimID /)
  traj_id = 'name_strlen'
  CALL check(  NF90_DEF_VAR(  ncid, REC_NAME,       NF90_DOUBLE, dimIDs, rec_varid  ) )
  CALL check(  NF90_DEF_VAR(  ncid, TRIM(LON_NAME), NF90_DOUBLE, dimIDs, x_varid    ) )
  CALL check(  NF90_DEF_VAR(  ncid, TRIM(LAT_NAME), NF90_DOUBLE, dimIDs, y_varid    ) )
  CALL check(  NF90_DEF_VAR(  ncid, Z_NAME,         NF90_DOUBLE, dimIDs, z_varid    ) )
  CALL check(  NF90_DEF_VAR(  ncid, "trajectory",   NF90_CHAR,           traj_varid ) )
  !
  ! ============================================================
  ! --  define data variables
  ! ============================================================
  DO jt=1,OutNetcdfANZ
    CALL check(  NF90_DEF_VAR(  ncid, TRIM(Diag_Name_Netcdf(jt)), NF90_DOUBLE  &
    &                         , dimIDs,  Diag_varID(jt)                         ) )
  END DO
  DO jt=1,SIZE(DiagERR_Name_Netcdf)
    CALL check(  NF90_DEF_VAR(  ncid, TRIM(DiagERR_Name_Netcdf(jt)), NF90_DOUBLE  &
    &                         , dimIDs,  DiagERR_varID(jt)                         ) )
  END DO
  !
  ! lwc
  CALL check ( NF90_DEF_VAR( ncid, 'Step_Size' , NF90_DOUBLE, dimIDS, StepSize_varid ) )
  CALL check ( NF90_DEF_VAR( ncid, 'Zenith' , NF90_DOUBLE, dimIDS, zenith_varid ) )
  CALL check ( NF90_DEF_VAR( ncid, 'SchwefelSumme' , NF90_DOUBLE, dimIDS, schwefel_varid ) )
  CALL check ( NF90_DEF_VAR( ncid, 'loc_error' , NF90_DOUBLE, dimIDS, error_varid ) )
  CALL check ( NF90_DEF_VAR( ncid, 'loc_error_spc' , NF90_DOUBLE, dimIDS, MaxErrorSpc_varid ) )
  IF ( ns_GAS >0 ) THEN
    CALL check ( NF90_DEF_VAR( ncid, 'GasSum'    , NF90_DOUBLE, dimIDS, gassum_varid ) )
  END IF
  IF ( ns_AQUA>0 ) THEN
    CALL check ( NF90_DEF_VAR( ncid, 'LWC_Level' , NF90_DOUBLE, dimIDS, LWC_varid ) )
    CALL check ( NF90_DEF_VAR( ncid, 'wetRadius' , NF90_DOUBLE, dimIDS, wetRadius_varid ) )
    CALL check ( NF90_DEF_VAR( ncid, 'pH_Value' , NF90_DOUBLE, dimIDS, pH_varid ) )
    CALL check ( NF90_DEF_VAR( ncid, 'AquaSum'   , NF90_DOUBLE, dimIDS, aquasum_varid ) )
  END IF
  IF ( Teq ) THEN
    CALL check ( NF90_DEF_VAR( ncid, 'Temperature' , NF90_DOUBLE, dimIDS, Temperature_varid ) )
  END IF
  !
  ! ============================================================
  ! --  define attributes (name, unit...) of the variables
  ! ============================================================
  CALL check( NF90_PUT_ATT( ncid, traj_varid, "cf_role", "trajectory_id") )
  CALL check( NF90_PUT_ATT( ncid, MaxErrorSpc_varid, "Species", "ErrorVal") )
  ! x = longitude
  CALL check( NF90_PUT_ATT( ncid, x_varid, "standard_name", TRIM(LON_LONGNAME) ) )
  CALL check( NF90_PUT_ATT( ncid, x_varid, "long_name",     TRIM(LON_LONGNAME) ) )
  CALL check( NF90_PUT_ATT( ncid, x_varid, NC_UNITS,        TRIM(LON_UNITS)    ) )
  CALL check( NF90_PUT_ATT( ncid, x_varid, "axis",          "x"                ) )
  CALL check( NF90_PUT_ATT( ncid, x_varid, "_CoordinateAxisType", "lon"        ) )
  ! y = latitude
  CALL check( NF90_PUT_ATT( ncid, y_varid, "standard_name", TRIM(LAT_LONGNAME) ) )
  CALL check( NF90_PUT_ATT( ncid, y_varid, "long_name",     TRIM(LAT_LONGNAME) ) )
  CALL check( NF90_PUT_ATT( ncid, y_varid, NC_UNITS,        TRIM(LAT_UNITS)    ) )
  CALL check( NF90_PUT_ATT( ncid, y_varid, "axis",          "y"                ) )
  CALL check( NF90_PUT_ATT( ncid, y_varid, "_CoordinateAxisType", "lat"        ) )
  ! z = altitude
  CALL check( NF90_PUT_ATT( ncid, z_varid, "standard_name", "altitude" ) )
  CALL check( NF90_PUT_ATT( ncid, z_varid, "long_name",     Z_LONGNAME ) )
  CALL check( NF90_PUT_ATT( ncid, z_varid, NC_UNITS,        Z_UNITS    ) )
  CALL check( NF90_PUT_ATT( ncid, z_varid, "positive",      "up"       ) )
  CALL check( NF90_PUT_ATT( ncid, z_varid, "axis",          "z"        ) )
  CALL check( NF90_PUT_ATT( ncid, z_varid, "_CoordinateAxisType", "z"  ) )
  ! rec = time
  CALL check( NF90_PUT_ATT( ncid, rec_varid, "standard_name",       "time"    ) )
  CALL check( NF90_PUT_ATT( ncid, rec_varid, "long_name",           REC_NAME  ) )
  CALL check( NF90_PUT_ATT( ncid, rec_varid, NC_UNITS,              REC_UNITS ) )
  CALL check( NF90_PUT_ATT( ncid, rec_varid, "_CoordinateAxisType", "time"    ) )
  ! Diagnose species
  DO jt=1,OutNetcdfANZ
    CALL check( NF90_PUT_ATT(ncid, Diag_varID(jt), NC_UNITS, Diag_UNITS(jt)) )
    CALL check( NF90_PUT_ATT(ncid, Diag_varID(jt), "long_name", TRIM(Diag_LongName(jt))) )
    CALL check( NF90_PUT_ATT(ncid, Diag_varID(jt), "_CoordinateAxes", "time") )
  END DO
  DO jt=1,SIZE(DiagERR_Name_Netcdf)
    CALL check( NF90_PUT_ATT(ncid, DiagERR_varID(jt), NC_UNITS, ' ' ) ) 
    CALL check( NF90_PUT_ATT(ncid, DiagERR_varID(jt), "long_name",    &
    &               'Error - '//TRIM(DiagERR_Name_Netcdf(jt))//' each step') )
  END DO
  !
  ! stepsize
  CALL check( NF90_PUT_ATT(ncid, StepSize_varid, NC_UNITS, '[sec]' ) )  
  CALL check( NF90_PUT_ATT(ncid, StepSize_varid, "long_name", 'step size in [sec]') ) 
  CALL check( NF90_PUT_ATT(ncid, StepSize_varid, "_CoordinateAxes", "time") )
  IF ( ns_GAS>0 ) THEN
    ! Gas summe
    CALL check( NF90_PUT_ATT(ncid, Gassum_varid, NC_UNITS, '[molec/m3]' ) )  
    CALL check( NF90_PUT_ATT(ncid, Gassum_varid, "Long_name", 'sum gaseus conc [molec/cm3]') ) 
    CALL check( NF90_PUT_ATT(ncid, Gassum_varid, "_CoordinateAxes", "time") )
  END IF
  ! Zenith angle
  CALL check( NF90_PUT_ATT(ncid, wetRadius_varid, NC_UNITS, '[-]' ) )  
  CALL check( NF90_PUT_ATT(ncid, wetRadius_varid, "long_name", 'Zenith angle [-]') ) 
  CALL check( NF90_PUT_ATT(ncid, wetRadius_varid, "_CoordinateAxes", "time") )
  ! schwefel summe
  CALL check( NF90_PUT_ATT(ncid, schwefel_varid, NC_UNITS, '[molec/m3]' ) )  
  CALL check( NF90_PUT_ATT(ncid, schwefel_varid, "long_name", 'Schwefelkonzentration') ) 
  CALL check( NF90_PUT_ATT(ncid, schwefel_varid, "_CoordinateAxes", "time") )
  ! local error estimation
  CALL check( NF90_PUT_ATT(ncid, error_varid, NC_UNITS, '[-]' ) )  
  CALL check( NF90_PUT_ATT(ncid, error_varid, "long_name", 'local error estimation') ) 
  CALL check( NF90_PUT_ATT(ncid, error_varid, "_CoordinateAxes", "time") )
  ! species of max local error estimation
  CALL check( NF90_PUT_ATT(ncid, MaxErrorSpc_varid, NC_UNITS, '[-]' ) )  
  CALL check( NF90_PUT_ATT(ncid, MaxErrorSpc_varid, "long_name", 'species of max error val') ) 
  CALL check( NF90_PUT_ATT(ncid, MaxErrorSpc_varid, "_CoordinateAxes", "time") )
  ! pH value
  IF ( ns_AQUA>0 ) THEN
    ! lwc
    CALL check( NF90_PUT_ATT(ncid, LWC_varid, NC_UNITS, '[l/m3]' ) )  
    CALL check( NF90_PUT_ATT(ncid, LWC_varid, "long_name", '[liter/m3]') ) 
    CALL check( NF90_PUT_ATT(ncid, LWC_varid, "_CoordinateAxes", "time") )
    ! aqua summe
    CALL check( NF90_PUT_ATT(ncid, Aquasum_varid, NC_UNITS, '[molec/m3]' ) )  
    CALL check( NF90_PUT_ATT(ncid, Aquasum_varid, "Long_name", 'sum aqueus conc [molec/cm3]') ) 
    CALL check( NF90_PUT_ATT(ncid, Aquasum_varid, "_CoordinateAxes", "time") )
    ! wet dropplet radius
    CALL check( NF90_PUT_ATT(ncid, wetRadius_varid, NC_UNITS, '[m]' ) )  
    CALL check( NF90_PUT_ATT(ncid, wetRadius_varid, "long_name", 'wet droplett radius [m]') ) 
    CALL check( NF90_PUT_ATT(ncid, wetRadius_varid, "_CoordinateAxes", "time") )
    CALL check( NF90_PUT_ATT(ncid, pH_varid, NC_UNITS, '[-]' ) )  
    CALL check( NF90_PUT_ATT(ncid, pH_varid, "long_name", 'pH-Value ( = -log10[Hp] )') ) 
    CALL check( NF90_PUT_ATT(ncid, pH_varid, "_CoordinateAxes", "time") )
  END IF
  ! save temperature for Teq simulation
  IF ( Teq ) THEN
    CALL check( NF90_PUT_ATT(ncid, Temperature_varid, NC_UNITS, '[K]' ) )  
    CALL check( NF90_PUT_ATT(ncid, Temperature_varid, "long_name", 'Temperature in Kelvin') ) 
    CALL check( NF90_PUT_ATT(ncid, Temperature_varid, "_CoordinateAxes", "temp") )
  END IF
  ! ============================================================
  ! --  End define mode
  ! ============================================================
  !
  CALL check( NF90_ENDDEF( ncid ) )
  !
!  ! ============================================================
  ! -- Write the record variable into the netCDF file.
  ! ============================================================
  !
  !IF (MPI_ID==0) print *,"  *****************************************************  "
  !IF (MPI_ID==0) print *,"  *** SUCCESS writing netcdf-file ", TRIM(NetcdfFile),  "!"
  !IF (MPI_ID==0) print *,"  *****************************************************  "
  !
  ! ============================================================
  ! --  Close netcdf file
  ! ============================================================
  !
  CALL check( NF90_CLOSE( ncid ) )
  !      
  !
  ! ============================================================
  CONTAINS
  !  
  !SUBROUTINE check_name(NetcdfName)
  !  CHARACTER(60) :: NetcdfName
  !  CHARACTER(60) :: NetcdfName_new
  !  !
  !  IF ( (INDEX(ADJUSTL(NetcdfName),'[')>0) .AND.&
  !  &    (INDEX(ADJUSTL(NetcdfName),']')>0) ) THEN
  !    NetcdfName_new = TRIM(NetcdfName(INDEX(NetcdfName,'[')+1:&
  !    &                                INDEX(NetcdfName,']')-1))
  !    NetcdfName = TRIM(NetcdfName_new)
  !  END IF
  !END SUBROUTINE check_name
  !
  SUBROUTINE check_name2(NetcdfName)
    CHARACTER(60) :: NetcdfName
    CHARACTER(60) :: NetcdfName_new
    !
    IF ( (INDEX(ADJUSTL(NetcdfName),'/')>0) ) THEN
      NetcdfName_new = TRIM(NetcdfName(1:INDEX( NetcdfName,'/')-1) // '_' //  &
      &                     NetcdfName(INDEX( NetcdfName,'/')+1:)             ) 
      NetcdfName = TRIM(NetcdfName_new)
    END IF
  END SUBROUTINE check_name2
  !
  SUBROUTINE check_name3(NetcdfName)
    CHARACTER(60) :: NetcdfName
    CHARACTER(60) :: NetcdfName_new
    !
    IF ( (INDEX(TRIM(ADJUSTL(NetcdfName)),'['))==1) THEN
      NetcdfName_new = '_' // TRIM(NetcdfName)
      NetcdfName = TRIM(NetcdfName_new)
    END IF
  END SUBROUTINE check_name3
  !
  !
  SUBROUTINE check(STATUS)
    INTEGER, INTENT ( in) :: STATUS
    ! 
    IF(STATUS /= nf90_noerr) THEN
      WRITE(*,*) '  Error with NetCDF data ', trim(nf90_strerror(STATUS)), STATUS
      STOP 
    END IF
  END SUBROUTINE check
  ! 
END SUBROUTINE InitNetCDF
!
!
    SUBROUTINE SetOutputNcdf( y, yout, Time, actLWC)

!==================================================================
!===  Conversion of Values into Output Array 
!==================================================================
!


     REAL(dp), INTENT(IN)    ::  y(:)         ! Concentrations of Species
     REAL(dp), INTENT(INOUT) ::  yout(:)      ! Output Array
     REAL(dp) ::  Time               ! Model Time
     REAL(dp) :: actLWC

     !-- internal variable
     INTEGER :: jt,  idx, iDiagSpc
     LOGICAL :: NaN=.FALSE.
!==================================================================
!===  Saving Output
!==================================================================
!
     yout(:) = ZERO
     IF ( Teq ) yOut(OutNetcdfDIM)=y(nDIM)
     jt=0

! internal calculations all in molec/cm3
!
! mol2part converts from molec/cm3 to mol/m3 and vice versa
     DO iDiagSpc=1,SIZE(OutNetcdfspc)
       jt  = jt + 1
       idx = OutNetcdfSpc(iDiagSpc)
       NaN = ISNAN(y(idx))
       IF      (OutNetcdfPhase(iDiagSpc)=='a'.AND..NOT.NaN) THEN
         yout(jt)   = y(idx) / (actLWC * mol2part)  ! convert to mol/l water
         yout(jt+1) = y(idx) / (mol2part) ! convert to mol/m3 air
         jt = jt + 1
         !
       ELSE IF (OutNetcdfPhase(iDiagSpc)=='g'.AND..NOT.NaN) THEN
         IF (Teq) THEN
           !yout(jt) = y(idx) * 1.0d6       ! [mol/cm3] to [mol/m3]
           yout(jt) = y(idx)          ! in mass fractions
         ELSE
           yout(jt) = y(idx) / (mol2part) 
         END IF
         !
       ELSE
         WRITE(*,*) ' '
         WRITE(*,*) ' '
         WRITE(*,*) ' '
         WRITE(*,*) '  Error: Some value is NaN ! '
         WRITE(*,*) '  -------------------------- '
         WRITE(*,*) '  spc idx   =  ', idx
         WRITE(*,*) '  spc name  =  ', y_name(idx)
         WRITE(*,*) '  spc val   =  ', y(idx)
         STOP 
       END IF
     END DO
!-------------------------------------------------------------------------
  END SUBROUTINE SetOutputNcdf
  !
  !
  SUBROUTINE StepNetcdf( t, yout, time_ind, otherStuff, ErrInd, Schwefel,Error)
  !
  !==================================================================
  !===  Initialization of ASCII Output File
  !==================================================================
  !  
  INTEGER :: time_ind         ! step number

  ! other stuff
  REAL(dp) :: otherStuff(6)! (/actLWC, StepSize, Gassum, Aquasum, wetRadius, Zenith/)
  REAL(dp) :: Schwefel
  REAL(dp) :: Error
  REAL(dp) :: pH
  INTEGER :: ErrInd(1,1)
  !
  REAL(dp) :: t                 ! Current time [in seconds]
  !
  REAL(dp) :: yOut(OutNetcdfDIM)               ! Output Array
  REAL(dp) :: timeLoc
  !
  !-- internal variable
  INTEGER :: jt
  REAL(dp), PARAMETER :: y_Min = 1.E-40    ! minimum for logarithmic plot


  time_ind = time_ind + 1
  !
  ! ====================================================================================
  ! == Open Netcdf-File in write mode to modify data ===================================
  ! ====================================================================================
  
  CALL check( NF90_OPEN( TRIM(NetcdfFile), NF90_WRITE, ncid ) ) 

  ! ====================================================================================
  ! == Write new data to netcdf-file ===================================================
  ! ====================================================================================
  timeLoc = t/HOUR   ! time output in hours
  IF ( Teq ) timeLoc=t ! time in seconds for Teq mechanism
  CALL check( NF90_PUT_VAR( ncid, rec_varid, timeLoc, start=(/time_ind/) ) )
  CALL check( NF90_PUT_VAR( ncid, x_varid,   rlon,    start=(/time_ind/) ) )
  CALL check( NF90_PUT_VAR( ncid, y_varid,   rlat,    start=(/time_ind/) ) )
  CALL check( NF90_PUT_VAR( ncid, z_varid,   altit,   start=(/time_ind/) ) ) 
  !
  !----------------------------------------------------------------
  DO jt=1,OutNetcdfANZ
    !OutNetcdf(jt) = MAX(yout(jt),1.d-40)        ! for logarithmic plot
    CALL check( NF90_PUT_VAR( ncid, Diag_varID(jt), yOut(jt), start = (/time_ind/) ) )
  END DO
  CALL check( NF90_PUT_VAR( ncid, StepSize_varid,  otherStuff(2), start = (/time_ind/) ) )
  IF ( ns_GAS>0 ) THEN
    CALL check( NF90_PUT_VAR( ncid, Gassum_varid,    otherStuff(3), start = (/time_ind/) ) )
  END IF
  CALL check( NF90_PUT_VAR( ncid, zenith_varid,    otherStuff(6), start = (/time_ind/) ) )
  CALL check( NF90_PUT_VAR( ncid, schwefel_varid,  Schwefel, start = (/time_ind/) ) )
  CALL check( NF90_PUT_VAR( ncid, error_varid,     error,    start = (/time_ind/) ) )
  CALL check( NF90_PUT_VAR( ncid, MaxErrorSpc_varid, ErrInd, start = (/time_ind/) ) )
  IF ( ns_AQUA>0 ) THEN
    CALL check( NF90_PUT_VAR( ncid, Aquasum_varid,   otherStuff(4), start = (/time_ind/) ) )
    CALL check( NF90_PUT_VAR( ncid, wetRadius_varid, otherStuff(5), start = (/time_ind/) ) )
    CALL check( NF90_PUT_VAR( ncid, LWC_varid,       otherStuff(1), start = (/time_ind/) ) )
    pH = -LOG10( yOut(pH_ind) )
    CALL check( NF90_PUT_VAR( ncid, pH_varid, pH, start = (/time_ind/) ) )
  END IF
  IF ( Teq ) THEN
    CALL check( NF90_PUT_VAR( ncid, Temperature_varid, yOut(OutNetcdfDIM), start = (/time_ind/) ) )
  END IF
  ! ====================================================================================
  ! == Close netcdf-file ===============================================================
  ! ====================================================================================
  !
  CALL check( nf90_sync(ncid) )
  CALL check( nf90_close(ncid) )
  !
  CONTAINS
    SUBROUTINE check(STATUS)
      INTEGER, INTENT ( in) :: STATUS
      !
      IF(STATUS /= nf90_noerr) THEN
        PRINT *, '  Error StepNetCDF :: ' , trim(nf90_strerror(STATUS))
        STOP "Stopped at StepNetCDF"
      END IF
    END SUBROUTINE check
    !
  END SUBROUTINE StepNetcdf

END MODULE NetCDF_Mod

