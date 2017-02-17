!========================================================================
! subroutine jac_dep(n,df) coputes the Jacobian for deposition
!========================================================================
!
SUBROUTINE Depos(nc,y,f)
   USE mo_reac
   USE mo_met
   USE mo_mphys
   USE mo_control

   IMPLICIT NONE
   INTEGER :: nc, nonz           ! Anzahl der gew. Dgl (= INTLAY(neq))
   REAL(8) :: y(nc,nt)
   REAL(8) :: f(nc,nt)

!---  internal variables 
   INTEGER :: na, ne
   INTEGER :: iFrac, ic, ncDep

   REAL(8) :: DropRad, vd_aqua
   REAL(8) :: vd1(nc,nt)

!-----------------------------------------------------------------------------------
!---  Deposition
   IF (nDep == 1)  THEN
      ncDep = 1
   ELSE IF (nDep == 2)  THEN
      ncDep = nc
   ELSE
      RETURN
   END IF

   vd1(:,:) = 0.e0

   IF (nDepAqua == 1)  THEN
      na = ntFrac*ntAqua + 1
      ne = ntFrac*ntAqua + ntGas
      DO ic=1,ncDep
         vd1(ic,na:ne) = vd(1:ntGas)
      END DO

!--- fixed predefine droplet deposition velocity (from microphysics)
      na = 1
      ne = ntFrac*ntAqua 
      vd1(:,na:ne) = VdepDrop
   ELSE IF (nDepAqua >= 2)  THEN
      DO ic=1,ncDep
         na = ntFrac*ntAqua + 1
         ne = ntFrac*ntAqua + ntGas
         vd1(ic,na:ne) = vd(1:ntGas)

!--- sedimentation velocity
         na = 1
         ne = ntAqua
         DO iFrac=1,ntFrac
            DropRad = cradius(ic,iFrac)
            vd_aqua = 0.3e-4 * (DropRad/5)**(1.93e0+0.05e0*log(DropRad))
            vd1(ic,na:ne) = vd_aqua
            na = na + ntAqua
            ne = ne + ntAqua
         END DO
      END DO
   END IF
!
!--- Deposition
   f(1:ncDep,:) = f(1:ncDep,:) - vd1(1:ncDep,:) * y(1:ncDep,:)
!
END SUBROUTINE Depos

