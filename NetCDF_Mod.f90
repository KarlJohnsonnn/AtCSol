MODULE NetCDF_Mod
!--- Netcdf Output
  USE Kind_Mod
  USE Meteo_Mod
  USE netcdf
  USE Reac_Mod
  USE Control_Mod
  USE MPI_Mod
  !
  IMPLICIT NONE

  TYPE NetCDF_T
    INTEGER                   :: n_Out         ! number of diagnosis species
    INTEGER                   :: iTime         ! Time step counter
    REAL(dp)                  :: Time          ! Time value
    INTEGER,      ALLOCATABLE :: Spc_Pos(:)    ! Species Positions (indices)
    CHARACTER(1), ALLOCATABLE :: Spc_Phase(:)  ! Phase of species g=gas , a=aqua, s=solid, p=parti
    REAL(dp),     ALLOCATABLE :: Spc_Conc(:)   ! Species concentrations at iTime
    REAL(dp)                  :: Temperature   ! Temperature value at iTime
    REAL(dp)                  :: LWC           ! Liquid water content value at iTime
    REAL(dp)                  :: StepSize      ! Step size value at iTime
    REAL(dp)                  :: GasConc       ! Sum of all gaseous species at iTime
    REAL(dp)                  :: AquaConc      ! Sum of all aqueous species at iTime
    REAL(dp)                  :: SulfConc      ! Sum of all sulphuric species at iTime
    REAL(dp)                  :: Zenith        ! Zenith value at iTime
    REAL(dp),     ALLOCATABLE :: WetRadius(:)  ! Size of cloud dropletts
    INTEGER                   :: MaxErrorSpc   ! Index of species which has max error
    REAL(dp)                  :: ROWerror      ! Error value ROS step at iTime
  END TYPE NetCDF_T

  TYPE(NetCDF_T) :: NetCDF

  INTEGER  :: iStpNetCDF 
  INTEGER  :: ncid
  REAL(dp) :: altit               ! altitude
  !
  !

  INTEGER, ALLOCATABLE :: Diag_varid(:), WetRadius_varid(:), pH_varid(:)
  INTEGER :: x_varid,y_varid,z_varid,rec_varid, LWC_varid, traj_varid
  INTEGER :: StepSize_varid, Gassum_varid, Aquasum_varid
  INTEGER :: zenith_varid, schwefel_varid, error_varid, Temperature_varid
  INTEGER :: MaxErrorSpc_varid

  INTEGER, ALLOCATABLE :: pH_ind(:)


  CONTAINS

  !==================================================================
  !===  Initialization of Netcdf Output File
  !==================================================================
  SUBROUTINE InitNetCDF

    ! -- Temporary variables --
    INTEGER :: j, pos, iFr, i
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
    CHARACTER(60),  ALLOCATABLE    :: DIAG_UNITS(:)
    !CHARACTER(100), ALLOCATABLE    :: Diag_Name(:)
    CHARACTER(200), ALLOCATABLE    :: Diag_LongName(:)
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
    iStpNetCDF   = 1
    NetCDF%iTime = 0
    NetCDF%n_Out = nNcdfGas + 2*nFrac*nNcdfAqua + nNcdfSolid + nNcdfParti 
    
    ! output array containing diagnosis species (and temperature if ChemKin mechanism)
    ALLOCATE(NetCDF%Spc_Conc(NetCDF%n_Out), NetCDF%WetRadius(nFrac)) 

    ! -- Allocate Netcdf Names etc.
    !
    ALLOCATE( Diag_LongName(NetCDF%n_Out) , Diag_Name(NetCDF%n_Out)    &
    &       , Diag_UNITS(NetCDF%n_Out)    , Diag_varid(NetCDF%n_Out)   &
    &       , WetRadius_varid(nFrac) , pH_varid(nFrac) , pH_ind(nFrac) )
    
    Diag_varID      = 0
    WetRadius_varid = 0

    j = 0
    ! gaseous species
    ALLOCATE(iNCout_G(0))
    DO iDiagSpc = 1 , nNcdfGas
      j = j + 1
      tmpName = ADJUSTL(y_name(iNcdfGas(iDiagSpc)))
  
      Diag_Name(j)     = TRIM(tmpName)
      Diag_LongName(j) = TRIM(tmpName)
      IF ( UnitGas == 1 ) THEN
        DIAG_UNITS(j)    = "mol/m3"
      ELSE
        DIAG_UNITS(j)    = "molec/cm3"
      END IF
  
      CALL check_name_slash(Diag_Name(j))
      CALL check_name_bracket(Diag_Name(j))
      iNCout_G = [iNCout_G, j]
    END DO
  
    ! aqueous species
    ALLOCATE(iNCout_A_l(0),iNCout_A_m3(0))
    DO iFr = 1 , nFrac
      DO iDiagSpc = 1 , nNcdfAqua  
        j = j + 1
        WRITE(tmpName,'(A,I0)') TRIM(y_name(iNcdfAqua(iDiagSpc)))//'_',iFr
        Diag_Name(j)     = TRIM(tmpName)//'_l'
        Diag_LongName(j) = TRIM(tmpName)//'AQUA'
        DIAG_UNITS(j)    = "mol/l"
        iNCout_A_l  = [iNCout_A_l, j]
  
        j = j + 1
        Diag_Name(j)     = TRIM(tmpName)//'_m3'
        Diag_LongName(j) = TRIM(tmpName)//'AIR'
        DIAG_UNITS(j)    = "mol/m3"
        iNCout_A_m3 = [iNCout_A_m3, j]
  
        IF (INDEX(tmpName,'Hp_')>0) THEN
          pH_ind(iFr) = j
        END IF
  
        CALL check_name_slash(Diag_Name(j))
        CALL check_name_bracket(Diag_Name(j))
      END DO
    END DO


    ! solid species
    ALLOCATE(iNCout_S(0))
    DO iDiagSpc = 1 , nNcdfSolid
      j = j + 1
      tmpName = ADJUSTL(y_name(iNcdfSolid(iDiagSpc)))
      Diag_Name(j)     = TRIM(tmpName)
      Diag_LongName(j) = TRIM(tmpName)
      DIAG_UNITS(j)    = "mol/m3"
      iNCout_S = [iNCout_S, j]
    END DO
  
    ! partic species
    ALLOCATE(iNCout_P(0))
    DO iDiagSpc = 1 , nNcdfParti
      j = j + 1
      tmpName = ADJUSTL(y_name(iNcdfParti(iDiagSpc)))
      Diag_Name(j)     = TRIM(tmpName)
      Diag_LongName(j) = TRIM(tmpName)
      DIAG_UNITS(j)    = "mol/m3"
      iNCout_P = [iNCout_P, j]
    END DO


    CALL TikZ_init('OUTPUT/'//TRIM(BSP)//'_tikz.dat',Diag_Name)

    
    ! ============================================================
    ! --  create Netcdf file (for each size bin / fraction)
    ! ============================================================
    ! Create new NetCDF File
    IF (MPI_master) THEN
      
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
      DO j = 1 , NetCDF%n_Out
        CALL check( NF90_DEF_VAR( ncid, TRIM(Diag_Name(j)), NF90_DOUBLE, dimIDs, Diag_varid(j)) )
      END DO

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
        CALL check ( NF90_DEF_VAR( ncid, 'AquaSum'   , NF90_DOUBLE, dimIDS, aquasum_varid ) )

        DO iFr=1,nFrac
          WRITE(tmpName,'(A,I0)') 'wetRadius_',iFr
          CALL check ( NF90_DEF_VAR( ncid, TRIM(tmpName), NF90_DOUBLE, dimIDS, wetRadius_varid(iFr) ) )
          WRITE(tmpName,'(A,I0)') 'pH_Value_',iFr
          CALL check ( NF90_DEF_VAR( ncid, TRIM(tmpName), NF90_DOUBLE, dimIDS, pH_varid(iFr) ) )
        END DO
      END IF

      CALL check ( NF90_DEF_VAR( ncid, 'Temperature' , NF90_DOUBLE, dimIDS, Temperature_varid ) )
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
      DO j=1,NetCDF%n_Out
        CALL check( NF90_PUT_ATT(ncid, Diag_varid(j), NC_UNITS, Diag_UNITS(j)) )
        CALL check( NF90_PUT_ATT(ncid, Diag_varid(j), "long_name", TRIM(Diag_LongName(j))) )
        CALL check( NF90_PUT_ATT(ncid, Diag_varid(j), "_CoordinateAxes", "time") )
      END DO

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
        ! wet dropplet radius and pH-value
        DO iFr=1,nFrac
          CALL check( NF90_PUT_ATT(ncid, wetRadius_varid(iFr), NC_UNITS, '[m]' ) )  
          CALL check( NF90_PUT_ATT(ncid, wetRadius_varid(iFr), "long_name", 'wet droplett radius [m]') ) 
          CALL check( NF90_PUT_ATT(ncid, wetRadius_varid(iFr), "_CoordinateAxes", "time") )
          ! pH value
          CALL check( NF90_PUT_ATT(ncid, pH_varid(iFr), NC_UNITS, '[-]' ) )  
          CALL check( NF90_PUT_ATT(ncid, pH_varid(iFr), "long_name", 'pH-Value ( = -log10[Hp] )') ) 
          CALL check( NF90_PUT_ATT(ncid, pH_varid(iFr), "_CoordinateAxes", "time") )
        END DO
      END IF
      ! save temperature for Teq simulation
      CALL check( NF90_PUT_ATT(ncid, Temperature_varid, NC_UNITS, '[K]' ) )  
      CALL check( NF90_PUT_ATT(ncid, Temperature_varid, "long_name", 'Temperature in Kelvin') ) 
      CALL check( NF90_PUT_ATT(ncid, Temperature_varid, "_CoordinateAxes", "temp") )

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

      WRITE(*,'(10X,A)') 'Init NetCDF................... done'
    END IF
    !      
    !
    ! ============================================================
    CONTAINS
    !  
    ! Check if species name contains slashes, if -> replace all '/' with '_'
    SUBROUTINE check_name_slash(Name)
      CHARACTER(100) :: Name
      INTEGER        :: pos

      DO
        pos = INDEX(TRIM(ADJUSTL(Name)),'/')
        IF ( pos > 0 ) THEN
          Name = TRIM(Name(1:pos-1)//'_'//Name(pos+1:))
        ELSE
          EXIT
        END IF
      END DO
    END SUBROUTINE check_name_slash

    ! Check if the first character is a bracket '[', if -> replace '[' with '_['
    SUBROUTINE check_name_bracket(Name)
      CHARACTER(100) :: Name
      IF ( (INDEX(TRIM(ADJUSTL(Name)),'['))==1) Name = '_'//TRIM(ADJUSTL(Name))
    END SUBROUTINE check_name_bracket
    
    
    SUBROUTINE check(STATUS)
      INTEGER, INTENT ( in) :: STATUS
      ! 
      IF(STATUS /= nf90_noerr) THEN
        WRITE(*,*) '  Error with NetCDF data ', trim(nf90_strerror(STATUS)), STATUS
        STOP 
      END IF
    END SUBROUTINE check
  
END SUBROUTINE InitNetCDF
!
!
!SUBROUTINE SetOutputNcdf( y, yout, Time, actLWC)
SUBROUTINE SetOutputNcdf(NCDF,Time,StpSize,iERR,ERR,Conc,Temp)
  ! OUT:
  TYPE(NetCDF_T) :: NCDF

  REAL(dp), INTENT(IN) :: Conc(:)
	REAL(dp), INTENT(IN) :: Time, StpSize, ERR, Temp
	INTEGER,  INTENT(IN) :: iERR(1,1)

  !-- internal variable
	INTEGER :: j,  idx, iDiagSpc, iFr
	LOGICAL :: NaN=.FALSE.
  REAL(dp), PARAMETER    :: y_Min = 1.e-40_dp    
  REAL(dp), ALLOCATABLE :: tConc(:), tG(:), tA(:), tS(:), tP(:), tAF(:), dbg(:)

  !==================================================================
  !===  Saving Output
  !==================================================================

  DO j=1,NetCDF%n_Out
    IF ( ISNAN(Conc(j)) ) THEN
      WRITE(*,*); WRITE(*,*); WRITE(*,*)
      WRITE(*,'(10X,A)')        '  ERROR:  Species concentration is NaN ! '
      WRITE(*,'(10X,A)')        '  -------------------------------------- '
      WRITE(*,'(10X,A,I0)')     '  NCDF idx      =  ', j
      WRITE(*,'(10X,A,Es12.4)') '  Species val   =  ', Conc(j)
      WRITE(*,*); WRITE(*,*); WRITE(*,*)
      CALL DropOut()
    END IF
  END DO

 
  tConc = [ MAX(Conc,y_Min) ]             ! minimum for logarithmic plot
  IF ( ChemKin ) THEN
    tConc = MoleConc_to_MoleFr(tConc(1:nspc))
    tG = [ tConc(iNcdfGas) ]     ! in mole fractions [-]
    NCDF%Temperature = tConc(nDIM)
  ELSE
    IF ( UnitGas == 0 ) THEN
      tG = [ tConc(iNcdfGas)   ] !  molec/cm3
    ELSE
      tG = [ tConc(iNcdfGas)   / mol2part ] ! convert to mol/m3
    END IF
    tA = [ tConc(iNcdfAqua)  / mol2part ] 
    tS = [ tConc(iNcdfSolid) / mol2part ] 
    tS = [ tConc(iNcdfParti) / mol2part ] 
    NCDF%Temperature = Temp
  END IF

  IF ( hasAquaSpc ) tAF = GatherAquaFractions( tA ) 

  IF ( MPI_master ) THEN
    NCDF%Time        = Time
    NCDF%StepSize    = StpSize
    NCDF%MaxErrorSpc = iERR(1,1)
    NCDF%ROWerror    = ERR
    

    ! NCDF%Spc_Conc = 
    ! ------------------------------------------------------------------------------ !
    ! |   GAS   |   AQ_1  |   AQ_2   |   ....   |   AQ_n   |   SOLID   |   PARTI   | !   
    ! ------------------------------------------------------------------------------ !
   
    IF (hasGasSpc) NCDF%Spc_Conc(iNCout_G) = tG

    IF (hasAquaSpc) THEN
      NCDF%LWC       = pseudoLWC(Time)
      NCDF%WetRadius = (Pi34 * NCDF%LWC / Mode%Number)**(rTHREE) * 0.1_dp
      NCDF%Spc_Conc(iNCout_A_l)  = tAF/NCDF%LWC
      NCDF%Spc_Conc(iNCout_A_m3) = tAF
    END IF

    IF (hasSolidSpc) NCDF%Spc_Conc(iNCout_S) = tS
    IF (hasPartiSpc) NCDF%Spc_Conc(iNCout_P) = tP

    IF ( hasPhotoReac ) NCDF%Zenith = Zenith(Time)
  END IF

  CONTAINS
  
  FUNCTION MoleConc_to_MoleFr(MoleConc) RESULT(MoleFr)
    REAL(dp), ALLOCATABLE :: MoleFr(:)   ! Mole fraction [mol/mol]
    REAL(dp), INTENT(IN)  :: MoleConc(:) ! Mole concentration  [mol/cm3]

    MoleFr = MoleConc / SUM( MoleConc)

  END FUNCTION MoleConc_to_MoleFr

END SUBROUTINE SetOutputNcdf
  
  
  
  SUBROUTINE StepNetcdf( NCDF )
  !
  !==================================================================
  !===  Initialization of ASCII Output File
  !==================================================================
  !  
  ! IN:
  TYPE(NetCDF_T) :: NCDF

  ! other stuff
  REAL(dp) :: pH(nFrac)
  INTEGER :: ErrInd(1,1)
  
  REAL(dp) :: timeLoc
  !
  !-- internal variable
  INTEGER :: j, iFr


  IF ( MPI_master ) THEN
    NetCDF%iTime = NetCDF%iTime + 1
    !
    ! ====================================================================================
    ! == Open Netcdf-File in write mode to modify data ===================================
    ! ====================================================================================
    CALL check( NF90_OPEN( TRIM(NetcdfFile), NF90_WRITE, ncid ) ) 

    ! ====================================================================================
    ! == Write new data to netcdf-file ===================================================
    ! ====================================================================================
    timeLoc = NCDF%Time/HOUR          ! time output in hours
    IF ( Teq ) timeLoc=NCDF%Time      ! time in seconds for Teq mechanism
    CALL check( NF90_PUT_VAR( ncid, rec_varid, timeLoc, start=(/NCDF%iTime/) ) )
    CALL check( NF90_PUT_VAR( ncid, x_varid,   rlon,    start=(/NCDF%iTime/) ) )
    CALL check( NF90_PUT_VAR( ncid, y_varid,   rlat,    start=(/NCDF%iTime/) ) )
    CALL check( NF90_PUT_VAR( ncid, z_varid,   altit,   start=(/NCDF%iTime/) ) ) 
    
    !----------------------------------------------------------------
    ! Write concentration
    DO j=1,NCDF%n_Out
      CALL check( NF90_PUT_VAR( ncid, Diag_varid(j), NCDF%Spc_Conc(j), start = (/NCDF%iTime/) ) )
    END DO

    ! Write stepsize
    CALL check( NF90_PUT_VAR( ncid, StepSize_varid, NCDF%StepSize, start = (/NCDF%iTime/) ) )
    
    !IF ( ns_GAS>0 ) THEN
    !  CALL check( NF90_PUT_VAR( ncid, Gassum_varid,    otherStuff(3), start = (/time_ind/) ) )
    !END IF
    !CALL check( NF90_PUT_VAR( ncid, schwefel_varid,  Schwefel, start = (/time_ind/) ) )
    
    ! write zenith, error value stepssize control, index of max. error species
    CALL check( NF90_PUT_VAR( ncid, zenith_varid,      NCDF%Zenith,      start = (/NCDF%iTime/) ) )
    CALL check( NF90_PUT_VAR( ncid, error_varid,       NCDF%ROWerror,    start = (/NCDF%iTime/) ) )
    CALL check( NF90_PUT_VAR( ncid, MaxErrorSpc_varid, NCDF%MaxErrorSpc, start = (/NCDF%iTime/) ) )

    ! write pH-values, droplet radii, liquid water content
    IF ( ns_AQUA>0 ) THEN
      DO iFr=1,nFrac
        pH(iFr) = -LOG10( NCDF%Spc_Conc(pH_ind(iFr)) )
        CALL check( NF90_PUT_VAR( ncid, wetRadius_varid(iFr), NCDF%WetRadius(iFr), start = (/NCDF%iTime/) ) )
        CALL check( NF90_PUT_VAR( ncid, pH_varid(iFr), pH(iFr), start = (/NCDF%iTime/) ) )
      END DO
      CALL check( NF90_PUT_VAR( ncid, LWC_varid, NCDF%LWC, start = (/NCDF%iTime/) ) )
      !CALL check( NF90_PUT_VAR( ncid, Aquasum_varid,   NCDF%AquaConc, start = (/NCDF%iTime/) ) )
    END IF

    ! Write temperatue value
    CALL check( NF90_PUT_VAR( ncid, Temperature_varid, NCDF%Temperature, start = (/NCDF%iTime/) ) )

    ! ====================================================================================
    ! == Close netcdf-file ===============================================================
    ! ====================================================================================
    !
    IF (Teq) THEN 
      CALL TikZ_write( timeLoc , NCDF%Spc_Conc , temp=NCDF%Temperature )
    ELSE
      CALL TikZ_write( timeLoc , NCDF%Spc_Conc , LWC=NCDF%LWC, zen=NCDF%Zenith )
    END IF
    CALL check( nf90_sync(ncid) )
    CALL check( nf90_close(ncid) )

  END IF
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


  ! Writing Output in format for LaTex TikZ Package
  !
  SUBROUTINE TikZ_init(Filename,y_names)

    CHARACTER(*) :: Filename
    CHARACTER(100) :: y_names(:)
    INTEGER :: j
    OPEN(UNIT=TikZUnit,FILE=Filename,STATUS='UNKNOWN')
    
    IF (Teq) THEN
      WRITE(TikZUnit,'(*(A,2X))') 'time','temp',(TRIM(y_names(j)),j=1,SIZE(y_names))
    ELSE
      WRITE(TikZUnit,'(*(A,2X))') 'time','LWC','solar',(TRIM(y_names(j)),j=1,SIZE(y_names))
    END IF

  END SUBROUTINE TikZ_init

  SUBROUTINE TikZ_write(time,conc,temp,LWC,zen)
    REAL(dp) :: time, conc(:)
    REAL(dp), OPTIONAL :: temp, LWC, zen 
    INTEGER :: j
    
    IF (Teq) THEN
      WRITE(TikZUnit,'(*(Es18.12,2X))') time,temp,(conc(j), j=1,NetCDF%n_Out)
    ELSE
      IF ( zen > pi/TWO ) THEN
        zen = pi/TWO
      END IF
      zen = COS(zen)
      WRITE(TikZUnit,'(*(Es18.12,2X))') time , LWC , zen , (conc(j), j=1,NetCDF%n_Out)
      IF (time==tEnd .OR. time==tBegin) THEN
        BACKSPACE(TikZUnit)
        WRITE(TikZUnit,'(*(Es18.12,2X))') time , 0.0_dp , 4.0_dp*zen , (conc(j), j=1,NetCDF%n_Out) 
      END IF
    END IF

  END SUBROUTINE TikZ_write

  SUBROUTINE TikZ_finished()
    CLOSE(TikZUnit)
  END SUBROUTINE TikZ_finished
!

END MODULE NetCDF_Mod

