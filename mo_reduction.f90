MODULE mo_reduction

  USE Kind_Mod
  USE ChemSys_Mod
  USE InputTool_Mod
  USE mo_reac

  CONTAINS

    SUBROUTINE Read_Target_Spc(FileName)

      CHARACTER(*)  :: FileName
      
      INTEGER       :: i, iPos
      REAL(dp)      :: c1
      CHARACTER(20) :: SpeciesName
      LOGICAL       :: Back = .FALSE.
      

      CALL OpenIniFile(FileName)
      i=0
      DO
        CALL LineFile( Back,                       &
        &              Start1 = 'BEGIN_TARGET_SPC',&
        &              End    = 'END_TARGET_SPC',  &
        &              Name1  = SpeciesName        )
        IF (Back) EXIT
        !
        IF ( PositionSpecies(SpeciesName) > 0 ) i = i + 1
      END DO

      CALL RewindFile
      CALL ClearIniFile
     
      ALLOCATE( Red_Names(i) , Red_Index(i) )
      ! read initial values
      i=0
      DO
        CALL LineFile( Back,                       &
        &              Start1 = 'BEGIN_TARGET_SPC',&
        &              End    = 'END_TARGET_SPC',  &
        &              Name1  = SpeciesName        )
        IF (Back)   EXIT
       
        iPos = PositionSpeciesAll(SpeciesName)
        
        IF (iPos>0) THEN
          i = i + 1
          READ(SpeciesName,*) Red_Names(i)
          Red_Index(i) = iPos
        END IF
      END DO
      CALL CloseIniFile
    END SUBROUTINE Read_Target_Spc


END MODULE mo_reduction
