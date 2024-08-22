! This file is part of the Kestrel software for simulations
! of sediment-laden Earth surface flows.
!
! Version 1.0
!
! Copyright 2023 Mark J. Woodhouse, Jake Langham, (University of Bristol).
!
! This program is free software: you can redistribute it and/or modify it 
! under the terms of the GNU General Public License as published by the Free 
! Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT 
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
! more details.
!
! You should have received a copy of the GNU General Public License along with 
! this program. If not, see <https://www.gnu.org/licenses/>. 


! This module defines closures for various physical processes that may be
! determined at runtime by the user, including the drag, morphodynamic processes
! and the inclusion of geometric corrections. It also specfies some useful
! functions (used elsewhere in the code) that return commonly computed flow
! observable, including the flow depth, speed, Shields number and so on.
module closures_module

   use set_precision_module, only: wp, pi
   use runsettings_module, only: RunSet
 
   implicit none

   private
   
   ! Helper functions to compute derived observables.
   public :: ComputeHn
   public :: FlowSquaredSpeedSlopeAligned
   public :: Density
 
   ! Main model closure declarations, followed by definitions for their interfaces.
   
   public :: GeometricCorrectionFactor, GeometricCorrectionFactor_gradin
   public :: IversonOuyangGeometricCorrectionFactor, NoGeometricCorrectionFactor
   public :: IversonOuyangGeometricCorrectionFactor_gradin, NoGeometricCorrectionFactor_gradin

   public :: DepositionClosure
   public :: NoDeposition, SimpleHinderedSettling, SpearmanManningHinderedSettling

   public :: ErosionClosure
   public :: SimpleErosion, FluidErosion, GranularErosion, MixedErosion, NoErosion

   public :: ErosionTransition
   public :: SmoothErosionTransition, StepErosionTransition, NoErosionTransition

   public :: DragClosure
   public :: ChezyDrag, CoulombDrag, VoellmyDrag, PouliquenDrag, Edwards2019Drag
   public :: ManningDrag, VariableDrag

   public :: fswitch
   public :: tanhSwitch, rat3Switch, cosSwitch
   public :: equalSwitch, zeroSwitch, oneSwitch, linearSwitch, stepSwitch

   public :: MorphoDamping
   public :: NoMorphoDamping, tanhMorphoDamping, rat3MorphoDamping

   public :: DesingularizeFunc
   public :: Desingularize_L1, Desingularize_L2, Desingularize_Linfty, Desingularize_step

   pointer :: GeometricCorrectionFactor
   interface
       pure function GeometricCorrectionFactor(RunParams, u) result(gam)
           import :: wp
           import :: RunSet
           type(RunSet), intent(in) :: RunParams
           real(kind=wp), dimension(:), intent(in) :: u
           real(kind=wp) :: gam
       end function GeometricCorrectionFactor
   end interface

   pointer :: GeometricCorrectionFactor_gradin
   interface
       pure function GeometricCorrectionFactor_gradin(dbdx, dbdy) result(gam)
           import :: wp
           real(kind=wp), intent(in) :: dbdx, dbdy
           real(kind=wp) :: gam
       end function GeometricCorrectionFactor_gradin
   end interface

   pointer :: DepositionClosure
   interface
       pure function DepositionClosure(RunParams, psi) result(deposition)
           import :: wp
           import :: RunSet
           type(RunSet), intent(in) :: RunParams
           real(kind=wp), intent(in) :: psi
           real(kind=wp) :: deposition
       end function DepositionClosure
   end interface

   pointer :: ErosionClosure
   interface
       pure function ErosionClosure(RunParams, uvect) result(erosion)
           import :: wp
           import :: RunSet
           type(RunSet), intent(in) :: RunParams
           real(kind=wp), dimension(:), intent(in) :: uvect
           real(kind=wp) :: erosion
       end function ErosionClosure
   end interface

   pointer :: ErosionTransition
   interface
       pure function ErosionTransition(RunParams, uvect) result(erosionRate)
           import :: wp
           import :: RunSet
           type(RunSet), intent(in) :: RunParams
           real(kind=wp), dimension(:), intent(in) :: uvect
           real(kind=wp) :: erosionRate
       end function ErosionTransition
   end interface

   pointer :: fswitch
   interface
       pure function fswitch(RunParams, psi) result(f)
           import :: wp
           import :: RunSet
           type(RunSet), intent(in) :: RunParams
           real(kind=wp), intent(in) :: psi
           real(kind=wp) :: f
       end function fswitch
   end interface

   pointer :: DragClosure
   interface
       pure function DragClosure(RunParams, uvect) result(friction)
           import :: wp
           import :: RunSet
           type(RunSet), intent(in) :: RunParams
           real(kind=wp), dimension(:), intent(in) :: uvect
           real(kind=wp) :: friction
       end function DragClosure
   end interface

   pointer :: MorphoDamping
   interface
       pure function MorphoDamping(RunParams, Hn) result(damping)
           import :: wp
           import :: RunSet
           type(RunSet), intent(in) :: RunParams
           real(kind=wp), intent(in) :: Hn
           real(kind=wp) :: damping
       end function MorphoDamping
   end interface

   pointer :: DesingularizeFunc
   interface
       pure function DesingularizeFunc(x,eps) result(x_recip)
            import :: wp
            real(kind=wp), intent(in) :: x
            real(kind=wp), intent(in) :: eps
            real(kind=wp) :: x_recip
       end function DesingularizeFunc
   end interface

contains

   ! Compute flow depth Hn from Hn = (w - b0 - bt)*gam. The calculation is
   ! broken into steps as typically bt << b0 ~ w, so loss of floating point
   ! precision can occur if done in one line.
   pure function ComputeHn(w, b0, bt, gam) result(Hn)
      implicit none

      real(kind=wp), intent(in) :: w
      real(kind=wp), intent(in) :: b0
      real(kind=wp), intent(in) :: bt
      real(kind=wp), intent(in) :: gam
      real(kind=wp) :: Hn

      Hn = -bt
      Hn = Hn + (w - b0)
      Hn = Hn * gam

   end function ComputeHn

   ! The magnitude of the flow speed, squared. If using geometric corrections,
   ! then this equals the horizontal component, plus a vertical component that
   ! we do not explicitly solve for, but may be obtained due to the assumption
   ! the velocity is always slope parallel. Specifically, the 3d velocity is
   ! (u, v, u*dbdx + v*dbdy) and we take the magnitude of this.
   pure function FlowSquaredSpeedSlopeAligned(RunParams, uvect) result(modu2)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: modu2

      real(kind=wp) :: u, v, dbdx, dbdy

      dbdx = uvect(RunParams%Vars%dbdx)
      dbdy = uvect(RunParams%Vars%dbdy)

      u = uvect(RunParams%Vars%u)
      if (RunParams%geometric_factors) then
         modu2 = u*u * (1.0_wp + dbdx*dbdx)
      else
         modu2 = u*u
      end if

      if (RunParams%isOneD) return

      v = uvect(RunParams%Vars%v)
      if (RunParams%geometric_factors) then
         modu2 = modu2 + v*v * (1.0_wp + dbdy*dbdy) + 2.0_wp*dbdx*dbdy*u*v
      else
         modu2 = modu2 + v*v
      end if
   end function FlowSquaredSpeedSlopeAligned

   ! Shields number for dilute flow.
   pure function FlowShieldsNumber(RunParams, uvect) result(shields)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: shields

      real(kind=wp) :: gred, modu2, chezy_friction

      ! g' -> g'_\perp = g' cos(theta)
      gred = RunParams%gred / GeometricCorrectionFactor(RunParams, uvect)

      modu2 = FlowSquaredSpeedSlopeAligned(RunParams, uvect)
      chezy_friction = RunParams%ChezyCo * modu2

      shields = chezy_friction / (gred * RunParams%SolidDiameter)
   end function FlowShieldsNumber

   ! Flow particle speed, u_p = sqrt(g' d) where g' = reduced gravity and d is
   ! the particle diamter.
   pure function FlowParticleSpeed(RunParams, uvect) result(u_p)
       implicit none

       type(RunSet), intent(in) :: RunParams
       real(kind=wp), dimension(:), intent(in) :: uvect
       real(kind=wp) :: u_p

       real(kind=wp) :: gred

       ! g' -> g'_\perp = g' cos(theta)
       gred = RunParams%gred / GeometricCorrectionFactor(RunParams, uvect)

       u_p = sqrt(gred * RunParams%SolidDiameter)
   end function FlowParticleSpeed

   ! Bulk density. This varies linearly with solid concentration psi on the
   ! interval 0<=psi<=1, but note psi is bounded to psi <= maximum packing.
   pure function Density(RunParams, psi) result(rho)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: rho

      real(kind=wp) :: rhow, rhos

      rhow = RunParams%rhow
      rhos = RunParams%rhos

      rho = rhow + (rhos - rhow) * psi
   end function Density

! -- Geometric correction factors --

   ! Return correction 'gamma' = sqrt(1 + dbdx^2 + dbdy^2) = 1 / cos(theta),
   ! originally suggested in Iverson & Ouyang, Rev. Geophys. (2015) to correct
   ! the bed evolution equation. If the shallow layer equations are derived by
   ! depth-integrating with respect to the bed-normal of an arbitrary surface,
   ! rather than a mean-slope-aligned frame, this factor appears elsewhere and
   ! is required for mass conservation. This is the perspective we take (see
   ! also https://arxiv.org/abs/2306.16185).
   pure function IversonOuyangGeometricCorrectionFactor(RunParams, u) result(gam)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp) :: gam

      real(kind=wp) :: dbdx, dbdy

      dbdx = u(RunParams%Vars%dbdx)
      dbdy = u(RunParams%Vars%dbdy)

      gam = sqrt(1.0_wp + dbdx*dbdx + dbdy*dbdy)

   end function IversonOuyangGeometricCorrectionFactor

   pure function IversonOuyangGeometricCorrectionFactor_gradin(dbdx, dbdy) result(gam)
      implicit none

      real(kind=wp), intent(in) :: dbdx, dbdy
      real(kind=wp) :: gam

      gam = sqrt(1.0_wp + dbdx*dbdx + dbdy*dbdy)

   end function IversonOuyangGeometricCorrectionFactor_gradin

   ! If geometric factors are switched off, set gamma = 1.
   pure function NoGeometricCorrectionFactor(RunParams, u) result(gam)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp) :: gam

      gam = 1.0_wp

   end function NoGeometricCorrectionFactor

   pure function NoGeometricCorrectionFactor_gradin(dbdx, dbdy) result(gam)
      implicit none

      real(kind=wp), intent(in) :: dbdx, dbdy
      real(kind=wp) :: gam

      gam = 1.0_wp

   end function NoGeometricCorrectionFactor_gradin

! -- Closures for the deposition law. --

   ! No deposition, return 0
   pure function NoDeposition(RunParams, psi) result(deposition)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: deposition

      deposition = 0.0_wp
   end function NoDeposition

   ! Simple quadratic hindered settling function, with zeros at 
   ! psi = 0 and psi = maximum packing
   pure function SimpleHinderedSettling(RunParams, psi) result(deposition)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: deposition

      deposition = psi * (1.0_wp - psi / RunParams%maxPack)
   end function SimpleHinderedSettling

   ! Hindered settling function from Spearman & Manning (2017)
   ! doi:10.1007/s10236-017-1034-7.
   pure function SpearmanManningHinderedSettling(RunParams, psi) result(deposition)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: deposition

      real(kind=wp) :: a, b

      a = 2.7_wp - 0.15_wp * RunParams%nsettling
      b = 0.62_wp * RunParams%nsettling - 1.46_wp
      deposition = psi * (1.0_wp - psi)**a * (1.0_wp - psi / RunParams%maxPack)**b
   end function SpearmanManningHinderedSettling

! -- Drag closures --
!
!    N.B. Because of the way the time stepper works (treating the drag
!    implicitly), these are required to define a 'friction' function F such that
!    the basal stress is \tau_x = -F \rho u/|(u,v)|.

   ! Chezy drag has friction function F = Cd * |u|^2.
   pure function ChezyDrag(RunParams, uvect) result(friction)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: modu2

      modu2 = FlowSquaredSpeedSlopeAligned(RunParams, uvect)
      
      friction = RunParams%ChezyCo * modu2
   end function ChezyDrag

   ! Coulomb drag has friction function F = \mu * g * Hn.
   pure function CoulombDrag(RunParams, uvect) result(friction)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: Hn
      real(kind=wp) :: gam
      real(kind=wp) :: g

      Hn = uvect(RunParams%Vars%Hn)
      gam = GeometricCorrectionFactor(RunParams, uvect)
      g = RunParams%g / gam
      
      friction = RunParams%CoulombCo * g * Hn
   end function CoulombDrag

   ! Voellmy drag is the sum of Chezy and Coulomb drag.
   pure function VoellmyDrag(RunParams, uvect) result(friction)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      friction = ChezyDrag(RunParams, uvect) + CoulombDrag(RunParams, uvect)
   end function VoellmyDrag

   ! Pouliquen drag has friction function F = \mu(I) * g * Hn, where mu(I) is a
   ! flow-dependent friction coefficient, I is the inertial number.
   ! See e.g. Pouliquen & Forterre, J. Fluid Mech. (2002).
   pure function PouliquenDrag(RunParams, uvect) result(friction)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: Hn
      real(kind=wp) :: modu2
      real(kind=wp) :: gam
      real(kind=wp) :: g
      real(kind=wp) :: mu

      Hn = uvect(RunParams%Vars%Hn)
      gam = GeometricCorrectionFactor(RunParams, uvect)
      g = RunParams%g / gam

      modu2 = FlowSquaredSpeedSlopeAligned(RunParams, uvect)

      mu = PouliquenFrictionCoefficient(RunParams, g, Hn, sqrt(modu2))

      if (modu2 > 0) then
         friction = mu * g * Hn
      else
         friction = 0.0_wp
      end if
   end function PouliquenDrag

   ! Pouliquen's intertial number dependent friction coefficient.
   pure function PouliquenFrictionCoefficient(RunParams, gcostheta, Hn, modu) result(mu)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: gcostheta
      real(kind=wp), intent(in) :: Hn
      real(kind=wp), intent(in) :: modu

      real(kind=wp) :: mu, mu1, mu2, beta, Fr, I

      mu1 = RunParams%PouliquenMinSlope
      mu2 = RunParams%PouliquenMaxSlope
      beta = RunParams%PouliquenBeta

      if (Hn > RunParams%heightThreshold) then
         ! Froude number
         Fr = modu / sqrt(gcostheta * Hn)
         ! Inertial number
         I = Fr * RunParams%SolidDiameter / Hn
         mu = mu1 + (mu2 - mu1) * I / (beta + I)
      else
         mu = mu1
      end if
   end function PouliquenFrictionCoefficient

   ! Extended 'Pouliquen'-like law featuring velocity-weakening behaviour at 
   ! low Froude number.
   ! See Edwards, Russell, Johnson, Gray, J. Fluid Mech. (2019).
   ! N.B. The Fr = 0 is omitted, since stopped regions require special 
   ! attention and this is difficult to implement in Kestrel currently.
   pure function Edwards2019Drag(RunParams, uvect) result(friction)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: Hn, modu, gam, gperp, mu, Fr
      real(kind=wp) :: betastar, beta, mu1, mu2, mu3, capgam, L, kappa

      Hn = uvect(RunParams%Vars%Hn)
      gam = GeometricCorrectionFactor(RunParams, uvect)
      gperp = RunParams%g / gam
      modu = sqrt(FlowSquaredSpeedSlopeAligned(RunParams, uvect))
      Fr = modu / sqrt(gperp * Hn)

      mu1 = RunParams%PouliquenMinSlope
      mu2 = RunParams%PouliquenMaxSlope
      mu3 = RunParams%PouliquenIntermediateSlope
      beta = RunParams%Pouliquenbeta
      betastar = RunParams%Edwards2019betastar
      kappa = RunParams%Edwards2019kappa
      capgam = RunParams%Edwards2019Gamma
      L = RunParams%SolidDiameter

      if (Fr > betastar) then
         friction = mu1 + (mu2 - mu1) / (1.0_wp + Hn * beta / (L * (Fr + capgam)))
      else
         friction = ((Fr / betastar)**kappa) * &
             (mu1 + (mu2 - mu1) / (1.0_wp + Hn * beta / (L * (betastar + capgam))) - &
             mu3 - (mu2 - mu1) / (1.0_wp + Hn / L)) + &
             mu3 + (mu2 - mu1) / (1.0_wp + Hn / L)
      end if

      friction = friction * gperp * Hn
   end function Edwards2019Drag

   ! Manning drag has friction function F = g * n^2 / (Hn^{1/3}) where n is the
   ! Manning coefficient.
   pure function ManningDrag(RunParams, uvect) result(friction)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: Hn, Hn_recip
      real(kind=wp) :: gam
      real(kind=wp) :: g
      real(kind=wp) :: ManningCo

      Hn = uvect(RunParams%Vars%Hn)
      gam = GeometricCorrectionFactor(RunParams, uvect)
      g = RunParams%g / gam

      ManningCo = RunParams%ManningCo

      Hn_recip = DesingularizeFunc(Hn, RunParams%heightThreshold*gam)

      friction = g * ManningCo*ManningCo * (Hn_recip**(1.0_wp/3.0_wp))

   end function ManningDrag

   ! Variable drag parameterization. This is a concentration-dependent
   ! variation from Chezy drag for dilute flows to Pouliquen drag for
   ! concentrated flow. The transition is specified by a switching function
   ! fSwitch.
   pure function VariableDrag(RunParams, uvect) result(friction)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: psi
      real(kind=wp) :: chezy_friction
      real(kind=wp) :: pouliquen_friction
      real(kind=wp) :: fc

      psi = uvect(RunParams%Vars%psi)
      
      chezy_friction = ChezyDrag(RunParams, uvect)
      pouliquen_friction = PouliquenDrag(RunParams, uvect)
      
      fc = fswitch(RunParams, psi)

      friction = chezy_friction * (1.0_wp - fc) + pouliquen_friction * fc
   end function VariableDrag

! -- Erosion closures --

   ! A fluid (Shields number) based erosion *without* a critical Shields number.
   ! The erosion rate is given by E = up * eps * S where up is the particle
   ! speed scale (= sqrt(g' d)) eps is the user defined erosion rate S is the
   ! Shields number.
   pure function SimpleErosion(RunParams, uvect) result(ero)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: ero

      real(kind=wp) :: shields

      shields = FlowShieldsNumber(RunParams, uvect)

      ero = RunParams%EroRate * shields
      ero = ero * FlowParticleSpeed(RunParams, uvect)
   end function SimpleErosion
   
   ! A fluid (Shields number) based erosion *with* a critical Shields number. The
   ! erosion rate is given by E = up * eps * (S - S_c) for S>S_c where up is the
   ! particle speed scale (= sqrt(g' d)) eps is the user defined erosion rate S
   ! is the Shields number, S_c is the critical shields.
   pure function FluidErosion(RunParams, uvect) result(ero)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: ero

      real(kind=wp) :: shields

      shields = FlowShieldsNumber(RunParams, uvect)
      if (shields > RunParams%CriticalShields) then
         ero = RunParams%EroRate * (shields - RunParams%CriticalShields)
   
         ero = ero * FlowParticleSpeed(RunParams, uvect)
      else
         ero = 0.0_wp
      end if
   end function FluidErosion

   ! A granular erosion based on a critical friction coefficient (from a
   ! 'neutral angle' so known as the neutral coefficient). The erosion rate is
   ! given by E = up * eps_g * (mu - mu_N) for mu>mu_N, where up is the particle
   ! speed scale (= sqrt(g' d)), eps is the user defined granular erosion rate,
   ! mu is the Pouliquen granular friction coefficient and mu_N is the neutral
   ! friction coefficient.
   pure function GranularErosion(RunParams, uvect) result(ero)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: ero

      real(kind=wp) :: Hn, modu2
      real(kind=wp) :: PouliquenMaxSlope, PouliquenMinSlope, PouliquenStaticSlope
      real(kind=wp) :: t1, gcostheta, mu, muNeutral

      Hn = uvect(RunParams%Vars%Hn)
      modu2 = FlowSquaredSpeedSlopeAligned(RunParams, uvect)

      PouliquenMinSlope = RunParams%PouliquenMinSlope
      PouliquenMaxSlope = RunParams%PouliquenMaxSlope
      t1 = tan(pi / 180.0_wp) ! t1 = tan(1 deg)
      PouliquenStaticSlope = (PouliquenMinSlope + t1) / (1.0_wp - PouliquenMinSlope * t1)
      ! granular erosion, eps_g*(mu - mu_n)
      gcostheta = RunParams%g / GeometricCorrectionFactor(RunParams, uvect)
      mu = PouliquenFrictionCoefficient(RunParams, gcostheta, Hn, sqrt(modu2))
      muNeutral = PouliquenMinSlope + (PouliquenStaticSlope - PouliquenMinSlope) /  &
         (1.0_wp + (Hn / 25.0_wp / RunParams%SolidDiameter)**2.0_wp)
      if (mu > muNeutral) then
         ero = RunParams%EroRateGranular * (mu - muNeutral)
         ero = ero*FlowParticleSpeed(RunParams, uvect)
      else
         ero = 0.0_wp
      end if
   end function GranularErosion

   ! Solid concentration weighted combination of fluid and granular erosion.
   ! The weighting is determined using the switching function.
   pure function MixedErosion(RunParams, uvect) result(ero)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: ero

      real(kind=wp) :: psi
      real(kind=wp) :: fc
      real(kind=wp) :: fluid_ero
      real(kind=wp) :: granular_ero

      psi = uvect(RunParams%Vars%psi)

      ! Get the weighting.
      fc = fswitch(RunParams, psi)

      fluid_ero = FluidErosion(RunParams, uvect)
      granular_ero = GranularErosion(RunParams, uvect)

      ero = (1.0_wp - fc) * fluid_ero + fc * granular_ero
   end function MixedErosion

   ! Set erosion rate to zero.
   pure function NoErosion(RunParams, uvect) result(ero)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: ero

      ero = 0.0_wp
   end function NoErosion

! -- Transition functions from erodible to inerodible at the maximum erosion
!    depth (RunParams%EroDepth). --
!
!    These functions specify how to transition when the elevation decrease
!    approaches the erosion depth.

   ! Rapid tanh transition.
   pure function SmoothErosionTransition(RunParams, uvect) result(erosionRate)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: erosionRate

      real(kind=wp) :: Bt
      real(kind=wp) :: eroDepth

      Bt = uvect(RunParams%Vars%Bt)

      eroDepth = RunParams%EroDepth

      erosionRate = 0.5_wp * (1.0_wp + tanh(1e5_wp * (Bt + eroDepth)))
   end function SmoothErosionTransition

   ! Step transition.
   pure function StepErosionTransition(RunParams, uvect) result(erosionRate)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: erosionRate

      real(kind=wp) :: Bt
      real(kind=wp) :: eroDepth

      Bt = uvect(RunParams%Vars%Bt)

      eroDepth = RunParams%EroDepth

      if (Bt < -eroDepth) then
         erosionRate = 0.0_wp
      else
         erosionRate = 1.0_wp
      end if
   end function StepErosionTransition

   ! No transition.
   pure function NoErosionTransition(RunParams, uvect) result(erosionRate)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: erosionRate

      erosionRate = 1.0_wp
   end function NoErosionTransition

! -- Switching functions for damping morphodynamics near critical erosion height
!    (RunParam%EroCriticalHeight). --
!
!    These functions define damping terms for Hn -> EroCriticalHeight
!    with damping = 0 => no morphodynamics;
!         damping = 1 => full morphodynamic (erosion/deposition) rates
!         0 < damping < 1 => reduced morphodynamic rates.


    ! No damping, return 1
    pure function NoMorphoDamping(RunParams, Hn) result(damping)
       implicit none

       type(RunSet), intent(in) :: RunParams
       real(kind=wp), intent(in) :: Hn
       real(kind=wp) :: damping

       damping = 1.0_wp
    end function NoMorphoDamping

    ! tanh form, this form ensures:
    !   damping -> as Hn -> 0,
    !   damping -> 1 as Hn -> infty (but damping ~ 1 for Hn>EroCriticalHeight)
    pure function tanhMorphoDamping(RunParams, Hn) result(damping)
       implicit none

       type(RunSet), intent(in) :: RunParams
       real(kind=wp), intent(in) :: Hn
       real(kind=wp) :: damping

       real(kind=wp) :: Hncrit

       Hncrit = RunParams%EroCriticalHeight
       damping = 0.5_wp * (1.0_wp + tanh(10.0_wp * log(Hn / Hncrit)))
    end function tanhMorphoDamping

    ! Cubic rational function with
    !   damping = 0 for Hn <= EroCriticalHeight,
    !   damping = 1 for Hn => 2 * EroCriticalHeight,
    pure function rat3MorphoDamping(RunParams, Hn) result(damping)
       implicit none

       type(RunSet), intent(in) :: RunParams
       real(kind=wp), intent(in) :: Hn
       real(kind=wp) :: damping

       real(kind=wp) :: Hncrit, t

       Hncrit = RunParams%EroCriticalHeight
       if (Hn < Hncrit) then
           damping = 0.0_wp
       else if (Hn > 2.0_wp * Hncrit) then
           damping = 1.0_wp
       else
           t = Hn / Hncrit - 1.0_wp
           damping = (t**3) / ((1.0_wp - t)**3 + t**3)
       end if
    end function rat3MorphoDamping

! -- Switching functions -- 
!    Chezy-Pouliquen switching term as a function of the solid fraction psi.

   ! tanh form.
   pure function tanhSwitch(RunParams, psi) result(f)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: f

      f = 0.5_wp * (1.0_wp + &
          tanh(RunParams%VoellmySwitchRate*(psi - RunParams%VoellmySwitchValue)))
   end function tanhSwitch

   ! Cubic rational function form.
   pure function rat3Switch(RunParams, psi) result(f)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: f

      real(kind=wp) :: a, b
      real(kind=wp) :: x

      a = RunParams%VoellmySwitchValue - 1.5_wp/RunParams%VoellmySwitchRate
      b = RunParams%VoellmySwitchValue + 1.5_wp/RunParams%VoellmySwitchRate

      if (psi<=a) then
          f = 0.0_wp
      elseif (psi>=b) then
          f = 1.0_wp
      else
          x = (psi-a)/(b-a)
          f = (x**3)/((1-x)**3 + x**3)
      end if
   end function rat3Switch

   ! Cosine form.
   pure function cosSwitch(RunParams, psi) result(f)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: f

      real(kind=wp) :: a, b
      real(kind=wp) :: x

      a = RunParams%VoellmySwitchValue - 0.25_wp*pi/RunParams%VoellmySwitchRate
      b = RunParams%VoellmySwitchValue + 0.25_wp*pi/RunParams%VoellmySwitchRate

      if (psi<=a) then
          f = 0.0_wp
      elseif (psi>=b) then
          f = 1.0_wp
      else
          x = (psi-a)/(b-a)
          f = 0.5_wp*(1.0_wp - cos(pi*x))
      end if
   end function cosSwitch

   ! 50--50 weighting.
   pure function equalSwitch(RunParams, psi) result(f)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: f
      
      f = 0.5_wp
   end function equalSwitch

   ! No switch, f(psi) = 0.
   pure function zeroSwitch(RunParams, psi) result(f)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: f
      
      f = 0.0_wp
   end function zeroSwitch

   ! Linear on 0<= psi <= maximum packing.
   pure function linearSwitch(RunParams, psi) result(f)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: f
      
      f = psi / RunParams%maxPack
   end function linearSwitch

   ! No switch, f(psi) = 1.
   pure function oneSwitch(RunParams, psi) result(f)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: f
      
      f = 1.0_wp
   end function oneSwitch

   ! Abrupt switch as RunParams%VoellmySwitchValue,
   ! f = 0 if psi < VoellmySwitchValue
   ! f = 1 if psi >= VoellmySwitchValue.
   pure function stepSwitch(RunParams, psi) result(f)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: f

      if (psi < RunParams%VoellmySwitchValue) then
         f = 0.0_wp
      else
         f = 1.0_wp
      end if
   end function stepSwitch

! -- Desingularization formulae -- 
!    We want to compute 1/x but avoid blow-up for small x.
!    These functions take x and a threshold, eps, as input and return
!    a 'desingularized' value for x_recip \approx 1/x.  We should have
!    x_recip = 1/x for x>eps, but a thresholded value for x < eps.
!    Commonly used desingularization formulae have the form
!    x_recip = alpha*x/(norm(x^2,max(x^2,eps^2)) for x<eps, with norm representing
!    the p-norm, and alpha = norm(1,1)
!    For examples:
!       Chertock et al (2015) uses the 1-norm,
!       Kurganov & Petrova (2007) uses the 2-norm,
!    and this could reasonably be continued, and culminate in the infinity-norm
!       x_recip = x/eps^2 for x<eps
!    The naming of the desingularization formulae below is based on this,
!    with the exception of Bollermann et al (2013) who use a step-function
!    with x_recip = 0 for x<eps -- we call this the Desingularize_step
    pure function Desingularize_L1(x,eps) result(x_recip)
        ! Desingularize using the L1 norm
        ! norm(x^2, max(x^2,eps^2)) = x^2 + max(x^2,eps^2)
        ! Chertock et al (2015) https://doi.org/10.1002/fld.4023
        implicit none
        real(kind=wp), intent(in) :: x
        real(kind=wp), intent(in) :: eps
        real(kind=wp) :: x_recip
        
        real(kind=wp) :: x2
        x2 = x * x

        x_recip = 2.0_wp*x/(x2 + max(x2,eps*eps))
        return
    end function Desingularize_L1

    pure function Desingularize_L2(x,eps) result(x_recip)
        ! Desingularize using the L2 norm for x<eps:
        ! norm(x^4, max(x^4,eps^4)) = sqrt(x^4 + max(x^4,eps^4))
        ! Kurganov & Petrova (2007) http://dx.doi.org/10.4310/CMS.2007.v5.n1.a6
        implicit none
        real(kind=wp), intent(in) :: x
        real(kind=wp), intent(in) :: eps
        real(kind=wp) :: x_recip

        real(kind=wp) :: x4
        x4 = x * x * x * x

        x_recip = sqrt(2.0_wp)*x/sqrt(x4 + max(x4,eps*eps*eps*eps))
        return
    end function Desingularize_L2

    pure function Desingularize_Linfty(x,eps) result(x_recip)
        ! Desingularize using the L-infty norm for x<eps:
        ! norm(x^2, max(x^2,eps^2)) = max(x^2,eps^2)
        implicit none
        real(kind=wp), intent(in) :: x
        real(kind=wp), intent(in) :: eps
        real(kind=wp) :: x_recip
        
        x_recip = x/(max(x*x,eps*eps))
        return
    end function Desingularize_Linfty

    pure function Desingularize_step(x,eps) result(x_recip)
        ! Desingularize using the step function
        ! x_recip = 1/x for x>eps, 0 otherwise
        ! Bollermann et al. (2013) https://doi.org/10.1007/s10915-012-9677-5
        ! It is likely that this formula will have VERY bad consequences -- when used to
        ! calculate drag, thin layers will become drag-free...
        implicit none
        real(kind=wp), intent(in) :: x
        real(kind=wp), intent(in) :: eps
        real(kind=wp) :: x_recip

        if (x>eps) then
            x_recip = 1.0_wp/x
        else
            x_recip = 0.0_wp
        end if
        return
    end function Desingularize_step

end module Closures_module
