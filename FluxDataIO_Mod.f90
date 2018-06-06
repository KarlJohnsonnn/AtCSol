MODULE FluxDataIO_Mod

  IMPLICIT NONE

  INTEGER, PARAMETER  :: dp = 8

  CHARACTER(80) :: FluxFile, FluxMetaFile, RedFile, TargetFile
  INTEGER       :: FluxUnit, FluxMetaUnit, RedUnit
  INTEGER       :: iStpFlux

  REAL(dp) :: TimeRedEnd, TimeRedStart
  REAL(dp) :: eps_red

  INTEGER        :: io_stat
  CHARACTER(100) :: io_msg


  CONTAINS

    ! module for writing flux data to unformatted file  


    SUBROUTINE InitFluxWriting


      !-----------------------------------------------------------------
      !--- NAMELISTS
      NAMELIST /SCENARIO/  TargetFile ,     &
      &                    TimeRedStart,    &
      &                    TimeRedEnd,      &
      &                    eps_red
		

      !--- Open run control file
      OPEN(UNIT=RedUnit,FILE=RedFile,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'opening run-file')

      !-----------------------------------------------------------------
      !---  SCENARIO
      !-----------------------------------------------------------------

      TargetFile   = ''
      TimeRedStart = 0.0d0
      TimeRedEnd   = 0.0d0
      eps_red      = 0.1_dp

      READ(RedUnit,SCENARIO,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading SCENARIO list')
      
      !--- Adjust Filenames
      CALL FileNameCheck(TargetFile,'TargetFile')
      TargetFile = ADJUSTL(TargetFile)

      IF (TimeRedStart >= TimeRedEnd) THEN
        WRITE(*,*) '  TStart >= TimeRedEnd ---> change timespan'
        STOP
      END IF

      iStpFlux = 0

    END SUBROUTINE InitFluxWriting

    SUBROUTINE OpenFile_wStream(UnitNr,FileName)
			INTEGER,      INTENT(IN) :: UnitNr
      CHARACTER(*), INTENT(IN) :: FileName
      OPEN(unit=UnitNr, file=FileName, status='replace', action='write', access='stream', iostat=io_stat)
      CALL file_err(FileName,io_stat)
      END SUBROUTINE OpenFile_wStream


      SUBROUTINE OpenFile_rStream(UnitNr,FileName)
      INTEGER,      INTENT(IN) :: UnitNr
      CHARACTER(*), INTENT(IN) :: FileName
      OPEN(unit=UnitNr, file=FileName, status='old', action='read', access='stream', iostat=io_stat)
      CALL file_err(FileName,io_stat)
    END SUBROUTINE OpenFile_rStream


    ! writing reaction rates , time and stepsize h to file via stream access
    SUBROUTINE StreamWriteFluxes(Rate,t,h)
      REAL(dp) :: Rate(:)
      REAL(dp) :: t , h
      
      INTEGER :: io_pos
      
      OPEN(unit=FluxUnit,      file=FluxFile,  status='old',   action='write', &
      &    position='append', access='stream', iostat=io_stat, iomsg=io_msg    )
      CALL file_err(FluxFile,io_stat,io_msg)
      INQUIRE(FluxUnit, POS=io_pos)
      WRITE(FluxUnit) Rate
      CLOSE(FluxUnit)
      
      iStpFlux   = iStpFlux + 1
      OPEN(unit=FluxMetaUnit, file=FluxMetaFile, status='old', action='write', position='append')
      WRITE(FluxMetaUnit,*) iStpFlux, io_pos ,t , h
      CLOSE(FluxMetaUnit)

    END SUBROUTINE StreamWriteFluxes


    SUBROUTINE file_err(FileName,io_status,io_message)
      CHARACTER(Len=*), INTENT(in) :: FileName
      INTEGER         , INTENT(in) :: io_status
      CHARACTER(Len=*), INTENT(in), OPTIONAL :: io_message
      IF (io_status /= 0) THEN
          WRITE(*,"(79('!'))")
          WRITE(*,'(A,I0)')    'ERROR operating on file:  '//TRIM(FileName)//'  with io status:  ',io_status 
          IF (PRESENT(io_message)) WRITE(*,'(A)')  'Message:  '//TRIM(io_message)
          WRITE(*,"(79('!'))")
          WRITE(*,*)'Exit ...'
          STOP
      END IF
    END SUBROUTINE file_err

    SUBROUTINE FileNameCheck(Name,miss)
      CHARACTER(*) :: Name
      CHARACTER(*) :: miss
      LOGICAL      :: ex

      INQUIRE(FILE=TRIM(Name), EXIST=ex)
      
      IF ( TRIM(Name) == '' .OR. .NOT.ex ) THEN
        WRITE(*,*); WRITE(*,*)
        WRITE(*,'(10X,A)') 'ERROR    Missing:  '//TRIM(miss)
        WRITE(*,'(10X,A)') '         FileName: '//TRIM(Name)
        WRITE(*,*); WRITE(*,*)
        STOP
      END IF
    END SUBROUTINE FileNameCheck

END MODULE FluxDataIO_Mod