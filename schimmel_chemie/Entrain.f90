!========================================================================
!===   Subroutine  Entrain  computes the entrainment rates 
!========================================================================
!
SUBROUTINE Entrain(nc,yGas,fGas, yAqua, fAqua)
   USE mo_reac
   USE mo_control
   USE mo_entrain
   USE mo_micfys

   IMPLICIT NONE
   INTEGER :: nc                            ! = 1 ... Number of grid cells
   REAL(8) :: yGas(ntGas), yAqua(ntAqua,ntFrac)
   REAL(8) :: fGas(ntGas), fAqua(ntAqua,ntFrac)

!---  internal variables 
   INTEGER :: na, ne
   INTEGER :: iFrac 

   REAL(8) :: SurSpec(ntAqua)

!-----------------------------------------------------------------------------------
!
!---  Entrainment of gas phase
   MuGas = MuGas0 * MuFac
   IF (ExGas >= 1)  THEN
      fGas(:) = fGas(:) - MuGas * (yGas(:)-SurGas(:))
   END IF

!---  Entrainment of aqueous phase
   IF (ExComp >= 1)  THEN
      DO iFrac=1,ntFrac
         SurSpec(:)     = 1.e3 * QS(iFrac) * SurComp(:,iFrac) / MolMass(:)
         fAqua(:,iFrac) = fAqua(:,iFrac) - MuPart * (yAqua(:,iFrac)-SurSpec(:))
      END DO
   END IF

!-----------------------------------------------------------------------------------
END SUBROUTINE Entrain

