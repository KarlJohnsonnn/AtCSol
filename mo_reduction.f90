MODULE mo_reduction

  USE Kind_Mod,       ONLY: dp
  USE ChemSys_Mod,    ONLY: PositionSpeciesAll
  USE InputTool_Mod
  !USE mo_reac

  CONTAINS

    SUBROUTINE Read_Target_Spc(Idx,Names,FileName)
      ! OUT:
      INTEGER,       ALLOCATABLE :: Idx(:)
      CHARACTER(80), ALLOCATABLE :: Names(:)
      ! IN:
      CHARACTER(*)  :: FileName
      ! TEMP:
      INTEGER       :: i, iPos
      REAL(dp)      :: c1
      CHARACTER(20) :: SpeciesName
      LOGICAL       :: Back = .FALSE.

      CALL OpenIniFile(FileName)
      i=0
      DO
        CALL LineFile( Back,                       &
        &              Start1 = 'TARGET_SPC',      &
        &              End    = 'END_TARGET_SPC',  &
        &              Name1  = SpeciesName        )
        IF (Back) EXIT
        !
        IF ( PositionSpeciesAll(SpeciesName) > 0 ) i = i + 1
      END DO

      CALL RewindFile
      CALL ClearIniFile
     
      ALLOCATE( Names(i) , Idx(i) )
      ! read initial values
      i=0
      DO
        CALL LineFile( Back,                       &
        &              Start1 = 'TARGET_SPC',      &
        &              End    = 'END_TARGET_SPC',  &
        &              Name1  = SpeciesName        )
        IF (Back)   EXIT
       
        iPos = PositionSpeciesAll(SpeciesName)
        
        IF (iPos>0) THEN
          i = i + 1
          READ(SpeciesName,*) Names(i)
          Idx(i) = iPos
        END IF
      END DO
      CALL CloseIniFile
    END SUBROUTINE Read_Target_Spc


END MODULE mo_reduction
