MODULE Kind_Mod

  USE MPI

  IMPLICIT NONE
  
  INCLUDE 'dmumps_struc.h'
  
  INTEGER, PARAMETER :: RealKind=8

  INTEGER :: idxHp
   
  REAL(RealKind) :: Timer_Start=0.0d0, Timer_Finish=0.d0
  REAL(RealKind) :: Timer_MumpsA, Timer_MumpsE
  REAL(RealKind) :: Tspan(2)
  
  INTEGER :: NumFac=0

  REAL(RealKind) :: TimeFac=0.0d0
  REAL(RealKind) :: TimeSolve=0.0d0
  REAL(RealKind) :: TimeRates=0.0d0
  REAL(RealKind) :: TimeRateSend=0.0d0
  REAL(RealKind) :: TimeJac=0.0d0
  REAL(RealKind) :: Time_Read=0.0d0
  REAL(RealKind) :: TimeSymbolic=0.0d0
  REAL(RealKind) :: TimeNetCDF=0.0d0
  
  REAL(RealKind) :: TimeIntegrationA=0.0d0
  REAL(RealKind) :: TimeIntegrationE=0.0d0
  REAL(RealKind) :: TimeRateA=0.0d0
  REAL(RealKind) :: TimeRateE=0.0d0
  REAL(RealKind) :: TimeRateSendA=0.0d0
  REAL(RealKind) :: TimeRateSendE=0.0d0
  REAL(RealKind) :: TimeJacobianA=0.0d0
  REAL(RealKind) :: TimeJacobianE=0.0d0
  REAL(RealKind) :: TimeNetCDFA=0.0d0

  LOGICAL :: ok=.FALSE.
  LOGICAL :: ckTEMP=.FALSE.
  !LOGICAL :: MatrixPrint=.false.
  !LOGICAL :: ERRORLOG=.false.

END MODULE Kind_Mod
