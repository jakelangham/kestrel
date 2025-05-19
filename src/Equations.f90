! This file is part of the Kestrel software for simulations
! of sediment-laden Earth surface flows.
!
! Version v1.1.1
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


! This module contains routines for computing individual terms in the governing
! equations, as well as the characteristic wave speeds for the problem, which
! are needed for the numerical method.
module equations_module

   use set_precision_module, only: wp
   use grid_module, only: GridType
   use runsettings_module, only: RunSet
   use closures_module

   implicit none

   private
   public :: XConvectionFlux
   private :: XConvectionFlux_with_gam, XConvectionFlux_without_gam
   public :: YConvectionFlux
   private :: YConvectionFlux_with_gam, YConvectionFlux_without_gam
   public :: XHydrostaticFlux
   public :: YHydrostaticFlux
   public :: XDiffusionFlux
   public :: YDiffusionFlux
   public :: XMaxWaveSpeeds
   private :: XMaxWaveSpeeds_with_gam, XMaxWaveSpeeds_without_gam
   public :: XMinWaveSpeeds
   private :: XMinWaveSpeeds_with_gam, XMinWaveSpeeds_without_gam
   public :: YMaxWaveSpeeds
   private :: YMaxWaveSpeeds_with_gam, YMaxWaveSpeeds_without_gam
   public :: YMinWaveSpeeds
   private :: YMinWaveSpeeds_with_gam, YMinWaveSpeeds_without_gam
   public :: ErosionDepositionTerms
   public :: ExplicitSourceTerms
   public :: ImplicitSourceTerms

   interface XConvectionFlux
      module procedure :: XConvectionFlux_with_gam, XConvectionFlux_without_gam
   end interface

   interface YConvectionFlux
      module procedure :: YConvectionFlux_with_gam, YConvectionFlux_without_gam
   end interface

   interface XMaxWaveSpeeds
      module procedure :: XMaxWaveSpeeds_with_gam, XMaxWaveSpeeds_without_gam
   end interface

   interface XMinWaveSpeeds
      module procedure :: XMinWaveSpeeds_with_gam, XMinWaveSpeeds_without_gam
   end interface

   interface YMaxWaveSpeeds
      module procedure :: YMaxWaveSpeeds_with_gam, YMaxWaveSpeeds_without_gam
   end interface

   interface YMinWaveSpeeds
      module procedure :: YMinWaveSpeeds_with_gam, YMinWaveSpeeds_without_gam
   end interface

contains

   ! Compute convection fluxes in the x direction for each of the governing
   ! equations and store in the vector xFlux.
   pure subroutine XConvectionFlux_without_gam(RunParams, u, xFlux)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(:), intent(out) :: xFlux(size(RunParams%iFlux))

      integer :: iw, iu, irhoHnu, irhoHnv, iHnpsi, iHn, irho
      real(kind=wp) :: gam, ugam

      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u
      irho = RunParams%Vars%rho

      gam = GeometricCorrectionFactor(RunParams, u)
      ugam = u(iu) * gam
      xFlux(iw) = u(iHn) * ugam
      xFlux(iHnpsi) = u(iHnpsi) * ugam
      xFlux(irhoHnu) = u(irhoHnu) * u(iu)
      xFlux(irhoHnv) = u(irhoHnv) * u(iu)

   end subroutine XConvectionFlux_without_gam

   pure subroutine XConvectionFlux_with_gam(RunParams, u, gam, xFlux)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), intent(in) :: gam
      real(kind=wp), dimension(:), intent(out) :: xFlux(size(RunParams%iFlux))

      integer :: iw, iu, irhoHnu, irhoHnv, iHnpsi, iHn, irho
      real(kind=wp) :: ugam

      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u
      irho = RunParams%Vars%rho

      ugam = u(iu) * gam

      xFlux(iw) = u(iHn) * ugam
      xFlux(iHnpsi) = u(iHnpsi) * ugam
      xFlux(irhoHnu) = u(irhoHnu) * u(iu)
      xFlux(irhoHnv) = u(irhoHnv) * u(iu)

   end subroutine XConvectionFlux_with_gam

   ! Compute convection fluxes in the y direction for each of the governing
   ! equations and store in the vector yFlux.
   pure subroutine YConvectionFlux_without_gam(RunParams, u, yFlux)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(:), intent(out) :: yFlux(size(RunParams%iFlux))

      integer :: iw, iv, irhoHnu, irhoHnv, iHnpsi, iHn, irho
      real(kind=wp) :: gam, vgam

      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      iHn = RunParams%Vars%Hn
      iv = RunParams%Vars%v
      irho = RunParams%Vars%rho

      gam = GeometricCorrectionFactor(RunParams, u)
      vgam = u(iv) * gam
      yFlux(iw) = u(iHn) * vgam
      yFlux(iHnpsi) = u(iHnpsi) * vgam
      yFlux(irhoHnu) = u(irhoHnu) * u(iv)
      yFlux(irhoHnv) = u(irhoHnv) * u(iv)

   end subroutine YConvectionFlux_without_gam

   pure subroutine YConvectionFlux_with_gam(RunParams, u, gam, yFlux)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), intent(in) :: gam
      real(kind=wp), dimension(:), intent(out) :: yFlux(size(RunParams%iFlux))

      integer :: iw, iv, irhoHnu, irhoHnv, iHnpsi, iHn, irho
      real(kind=wp) :: vgam
      
      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      iHn = RunParams%Vars%Hn
      iv = RunParams%Vars%v
      irho = RunParams%Vars%rho

      vgam = u(iv) * gam

      yFlux(iw) = u(iHn) * vgam
      yFlux(iHnpsi) = u(iHnpsi) * vgam
      yFlux(irhoHnu) = u(irhoHnu) * u(iv)
      yFlux(irhoHnv) = u(irhoHnv) * u(iv)

   end subroutine YConvectionFlux_with_gam

   ! Compute hydrostatic pressure fluxes in the x direction for each of the
   ! governing equations and store in the vector xFlux.
   pure subroutine XHydrostaticFlux(RunParams, u, xFlux)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(:), intent(out) :: xFlux(size(RunParams%iFlux))

      real(kind=wp) :: g, hp
      integer :: iw, irhoHnu, irhoHnv, iHnpsi, iHn, ib0, ibt, irho

      g = RunParams%g

      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      iHn = RunParams%Vars%Hn
      ib0 = RunParams%Vars%b0
      ibt = RunParams%Vars%bt
      irho = RunParams%Vars%rho

      xFlux(iw) = 0.0_wp
      xFlux(iHnpsi) = 0.0_wp

      hp = -u(ibt)
      hp = hp + (u(iw) - u(ib0))
      xFlux(irhoHnu) = 0.5_wp * g * u(irho) * hp * hp
      xFlux(irhoHnv) = xFlux(irhoHnu)

   end subroutine XHydrostaticFlux

   ! Compute hydrostatic pressure fluxes in the y direction for each of the
   ! governing equations and store in the vector yFlux.
   pure subroutine YHydrostaticFlux(RunParams, u, yFlux)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(:), intent(out) :: yFlux(size(RunParams%iFlux))

      real(kind=wp) :: g, hp
      integer :: iw, irhoHnu, irhoHnv, iHnpsi, iHn, ib0, ibt, irho

      g = RunParams%g

      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      iHn = RunParams%Vars%Hn
      ib0 = RunParams%Vars%b0
      ibt = RunParams%Vars%bt
      irho = RunParams%Vars%rho

      yFlux(iw) = 0.0_wp
      yFlux(iHnpsi) = 0.0_wp

      hp = -u(ibt)
      hp = hp + (u(iw) - u(ib0))
      yFlux(irhoHnu) = 0.5_wp * g * u(irho) * hp * hp
      yFlux(irhoHnv) = yFlux(irhoHnu)

   end subroutine YHydrostaticFlux

   ! Compute diffusion fluxes in the x direction for each of the governing
   ! equations and store in the vector xFlux. These implement an eddy viscosity
   ! parametrisation.
   pure subroutine XDiffusionFlux(RunParams, u, dudx, xFlux)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(:), intent(in) :: dudx
      real(kind=wp), dimension(:), intent(out) :: xFlux(size(RunParams%iFlux))

      real(kind=wp) :: nu
      real(kind=wp) :: nurhoHn

      integer :: iw, iHn, irho, irhoHnu, irhoHnv, iHnpsi, iu, iv

      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u
      iv = RunParams%Vars%v
      irho = RunParams%Vars%rho

      xFlux(iw) = 0.0_wp
      xFlux(iHnpsi) = 0.0_wp

      if (u(iHn) < 0.0_wp) then ! just in case
         xFlux(irhoHnu) = 0.0_wp
         xFlux(irhoHnv) = 0.0_wp
      else
         nu = RunParams%EddyViscosity
         nurhoHn = nu*u(irho)*u(iHn)
         xFlux(irhoHnu) = nurhoHn*dudx(iu)
         xFlux(irhoHnv) = nurhoHn*dudx(iv)
      end if
   end subroutine XDiffusionFlux

   ! Compute diffusion fluxes in the y direction for each of the governing
   ! equations and store in the vector yFlux. These implement an eddy viscosity
   ! parametrisation.
   pure subroutine YDiffusionFlux(RunParams, u, dudy, yFlux)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: u
      real(kind=wp), dimension(:), intent(in) :: dudy
      real(kind=wp), dimension(:), intent(out) :: yFlux(size(RunParams%iFlux))

      real(kind=wp) :: nu
      real(kind=wp) :: nurhoHn

      integer :: iw, iHn, irho, irhoHnu, irhoHnv, iHnpsi, iu, iv

      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u
      iv = RunParams%Vars%v
      irho = RunParams%Vars%rho

      yFlux(iw) = 0.0_wp
      yFlux(iHnpsi) = 0.0_wp

      if (u(iHn) < 0.0_wp) then ! just in case
         yFlux(irhoHnu) = 0.0_wp
         yFlux(irhoHnv) = 0.0_wp
      else
         nu = RunParams%EddyViscosity
         nurhoHn = nu*u(irho)*u(iHn)
         yFlux(irhoHnu) = nurhoHn*dudy(iu)
         yFlux(irhoHnv) = nurhoHn*dudy(iv)
      end if
   end subroutine YDiffusionFlux

   ! Given a solution vector uvect, compute its maximum characteristic wave
   ! speed in the x direction, return in xMaxWS.
   pure function XMaxWaveSpeeds_without_gam(RunParams, uvect) result(xMaxWS)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: xMaxWS

      real(kind=wp) :: Hn, u, gam, by

      integer :: iHn, iu

      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u

      Hn = uvect(iHn)

      ! just in case
      if (Hn <= 0.0_wp) then
         Hn = 0.0_wp
      end if

      u = uvect(iu)

      if (RunParams%geometric_factors) then
         gam = GeometricCorrectionFactor(RunParams, uvect)
         by = uvect(RunParams%Vars%dbdy)
         xMaxWS = u + sqrt(RunParams%g*Hn*(1.0_wp + by*by)/(gam*gam*gam))
      else
         xMaxWS = u + sqrt(RunParams%g*Hn)
      end if
   end function XMaxWaveSpeeds_without_gam

   pure function XMaxWaveSpeeds_with_gam(RunParams, uvect, gam) result(xMaxWS)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp), intent(in) :: gam
      real(kind=wp) :: xMaxWS

      real(kind=wp) :: Hn, u, by

      integer :: iHn, iu

      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u

      Hn = uvect(iHn)

      ! just in case
      if (Hn <= 0.0_wp) then
         Hn = 0.0_wp
      end if

      u = uvect(iu)

      if (RunParams%geometric_factors) then
         by = uvect(RunParams%Vars%dbdy)
         xMaxWS = u + sqrt(RunParams%g*Hn*(1.0_wp + by*by)/(gam*gam*gam))
      else
         xMaxWS = u + sqrt(RunParams%g*Hn)
      end if
   end function XMaxWaveSpeeds_with_gam

   ! Given a solution vector uvect, compute its minimum characteristic wave
   ! speed in the x direction, return in xMinWS.
   pure function XMinWaveSpeeds_without_gam(RunParams, uvect) result(xMinWS)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: xMinWS

      real(kind=wp) :: Hn, u, gam, by

      integer :: iHn, iu

      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u

      Hn = uvect(iHn)

      ! just in case
      if (Hn <= 0.0_wp) then
         Hn = 0.0_wp
      end if

      u = uvect(iu)
      
      if (RunParams%geometric_factors) then
         gam = GeometricCorrectionFactor(RunParams, uvect)
         by = uvect(RunParams%Vars%dbdy)
         xMinWS = u - sqrt(RunParams%g*Hn*(1.0_wp + by*by)/(gam*gam*gam))
      else
         xMinWS = u - sqrt(RunParams%g*Hn)
      end if
   end function XMinWaveSpeeds_without_gam

   pure function XMinWaveSpeeds_with_gam(RunParams, uvect, gam) result(xMinWS)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp), intent(in) :: gam
      real(kind=wp) :: xMinWS

      real(kind=wp) :: Hn, u, by

      integer :: iHn, iu

      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u

      Hn = uvect(iHn)

      ! just in case
      if (Hn <= 0.0_wp) then
         Hn = 0.0_wp
      end if

      u = uvect(iu)

      if (RunParams%geometric_factors) then
         by = uvect(RunParams%Vars%dbdy)
         xMinWS = u - sqrt(RunParams%g*Hn*(1.0_wp + by*by)/(gam*gam*gam))
      else
         xMinWS = u - sqrt(RunParams%g*Hn)
      end if
   end function XMinWaveSpeeds_with_gam

   ! Given a solution vector uvect, compute its maximum characteristic wave
   ! speed in the y direction, return in yMaxWS.
   pure function YMaxWaveSpeeds_without_gam(RunParams, uvect) result(yMaxWS)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: yMaxWS

      real(kind=wp) :: Hn, v, gam, bx

      integer :: iHn, iv

      iHn = RunParams%Vars%Hn
      iv = RunParams%Vars%v

      Hn = uvect(iHn)

      ! just in case
      if (Hn <= 0.0_wp) then
         Hn = 0.0_wp
      end if

      v = uvect(iv)

      if (RunParams%geometric_factors) then
         gam = GeometricCorrectionFactor(RunParams, uvect)
         bx = uvect(RunParams%Vars%dbdx)
         yMaxWS = v + sqrt(RunParams%g*Hn*(1.0_wp + bx*bx)/(gam*gam*gam))
      else
         yMaxWS = v + sqrt(RunParams%g*Hn)
      end if
   end function YMaxWaveSpeeds_without_gam

   pure function YMaxWaveSpeeds_with_gam(RunParams, uvect, gam) result(yMaxWS)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp), intent(in) :: gam
      real(kind=wp) :: yMaxWS

      real(kind=wp) :: Hn, v, bx

      integer :: iHn, iv

      iHn = RunParams%Vars%Hn
      iv = RunParams%Vars%v

      Hn = uvect(iHn)

      ! just in case
      if (Hn <= 0.0_wp) then
         Hn = 0.0_wp
      end if

      v = uvect(iv)

      if (RunParams%geometric_factors) then
         bx = uvect(RunParams%Vars%dbdx)
         yMaxWS = v + sqrt(RunParams%g*Hn*(1.0_wp + bx*bx)/(gam*gam*gam))
      else
         yMaxWS = v + sqrt(RunParams%g*Hn)
      end if
   end function YMaxWaveSpeeds_with_gam

   ! Given a solution vector uvect, compute its minimum characteristic wave
   ! speed in the y direction, return in yMinWS.
   pure function YMinWaveSpeeds_without_gam(RunParams, uvect) result(yMinWS)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: yMinWS

      real(kind=wp) :: Hn, v, gam, bx

      integer :: iHn, iv

      iHn = RunParams%Vars%Hn
      iv = RunParams%Vars%v

      Hn = uvect(iHn)

      ! just in case
      if (Hn <= 0.0_wp) then
         Hn = 0.0_wp
      end if

      v = uvect(iv)

      if (RunParams%geometric_factors) then
         gam = GeometricCorrectionFactor(RunParams, uvect)
         bx = uvect(RunParams%Vars%dbdx)
         yMinWS = v - sqrt(RunParams%g*Hn*(1.0_wp + bx*bx)/(gam*gam*gam))
      else
         yMinWS = v - sqrt(RunParams%g*Hn)
      end if
   end function YMinWaveSpeeds_without_gam

   pure function YMinWaveSpeeds_with_gam(RunParams, uvect, gam) result(yMinWS)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp), intent(in) :: gam
      real(kind=wp) :: yMinWS

      real(kind=wp) :: Hn, v, bx

      integer :: iHn, iv

      iHn = RunParams%Vars%Hn
      iv = RunParams%Vars%v

      Hn = uvect(iHn)

      ! just in case
      if (Hn <= 0.0_wp) then
         Hn = 0.0_wp
      end if

      v = uvect(iv)

      if (RunParams%geometric_factors) then
         bx = uvect(RunParams%Vars%dbdx)
         yMinWS = v - sqrt(RunParams%g*Hn*(1.0_wp + bx*bx)/(gam*gam*gam))
      else
         yMinWS = v - sqrt(RunParams%g*Hn)
      end if
   end function YMinWaveSpeeds_with_gam

   ! Given a solution vector uvect, compute corresponding erosion and deposition
   ! rates, returning in E and D respectively.
   pure subroutine ErosionDepositionTerms(RunParams, uvect, gam, E, D)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp), intent(in) :: gam
      real(kind=wp), intent(out) :: E
      real(kind=wp), intent(out) :: D

      real(kind=wp) :: Hn, Hncrit, damping

      Hn = uvect(RunParams%Vars%Hn)

      E = Erosion(RunParams, uvect, gam)
      D = Deposition(RunParams, uvect)
        
      ! Damp morphodynamics at small length scales Hn ~ Hncrit (various
      ! different damping functions are given in Closures.f90).
      Hncrit = RunParams%EroCriticalHeight
      damping = MorphoDamping(RunParams, Hn)
      E = E * damping
      D = D * damping
   end subroutine ErosionDepositionTerms

   ! Given a solution vector uvect, compute the corresponding erosion rate,
   ! returning in ero.
   pure function Erosion(RunParams, uvect, gam) result(ero)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp), intent(in) :: gam
      real(kind=wp) :: ero

      ! Get erosion rate from chosen closure.
      ero = ErosionClosure(RunParams, uvect, gam)

      ! Alter erosion rate with flow depth and smoothing choice.
      ero = ero * ErosionTransition(RunParams, uvect)
      
   end function Erosion

   ! Given a solution vector uvect, compute the corresponding deposition rate,
   ! returning in dep.
   pure function Deposition(RunParams, uvect) result(dep)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: dep

      real(kind=wp) :: psi, alpha

      psi = uvect(RunParams%Vars%psi)

      if (psi >= RunParams%maxPack) then
         alpha = 0.0_wp
      else if (psi > 0.0_wp) then
         alpha = DepositionClosure(RunParams, psi)
      else
         alpha = 0.0_wp
      end if

      dep = RunParams%ws0 * alpha
   end function

   ! Compute (hydraulic) explicit source terms, given flow at (x,y), time t and
   ! with fields in uvect, for each of the governing equations. Store the result
   ! in the vector stvect.
   !
   ! The source terms that are time stepped explicitly are: flux sources and
   ! downslope gravitational forcing from the weight of the flow..
   pure subroutine ExplicitSourceTerms(RunParams, grid, t, x, y, uvect, gam, tile_has_source, stvect)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(in) :: grid
      real(kind=wp), intent(in) :: t
      real(kind=wp), intent(in) :: x
      real(kind=wp), intent(in) :: y
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp), intent(in) :: gam
      logical, intent(in) :: tile_has_source
      real(kind=wp), dimension(:), intent(out) :: stvect(size(uvect))

      integer :: nSrc
      real(kind=wp) :: srcX, srcY, srcR
      integer :: nFluxSeries

      real(kind=wp), dimension(:), allocatable :: Qf, rhoQf, psiQf
      real(kind=wp) :: Qt, rhoQt, psiQt, Mt
      real(kind=wp) :: Qa, Qb
      real(kind=wp) :: psia, psib, psif
      real(kind=wp) :: ta, tb, tfrac
      real(kind=wp) :: rho
      real(kind=wp) :: sourceProp
      integer :: J, K

      integer :: iw, irhoHnu, irhoHnv, iHnpsi, irho
      integer :: iHn, iu, iv, ib0, ibt

      real(kind=wp) :: heightThreshold
      real(kind=wp) :: g, hp_o_gam, dbdx, dbdy
      real(kind=wp) :: d2bdxx, d2bdyy, d2bdxy
      real(kind=wp) :: d2bdtx, d2bdty
      real(kind=wp) :: u, v
      real(kind=wp) :: curvatureTerm

      stvect(:) = 0.0_wp

      heightThreshold = RunParams%heightThreshold

      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      irho = RunParams%Vars%rho
      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u
      iv = RunParams%Vars%v
      ib0 = RunParams%Vars%b0
      ibt = RunParams%Vars%bt

      nSrc = RunParams%nSources

      Qt = 0.0_wp
      rhoQt = 0.0_wp
      psiQt = 0.0_wp
      Mt = 0.0_wp

      ! If there is a flux source present, compute the total flux of water and
      ! solids through the current cell, storing the results in Qt, psiQt
      ! respectively. The logic below determines the sum of all sources present,
      ! finding the flux value by interpolating it from the user-specified time
      ! series data sttored in RunParams%FluxSources(:)%{time,flux,psi,etc}.
      if (tile_has_source) then
         if (nSrc >= 1) then
            allocate (Qf(nSrc), rhoQf(nSrc), psiQf(nSrc), source=0.0_wp)
            ! Qf(:) = 0.0_wp
            ! rhoQf(:) = 0.0_wp
            ! psiQf(:) = 0.0_wp
            do J = 1, nSrc
               srcX = RunParams%FluxSources(J)%x
               srcY = RunParams%FluxSources(J)%y
               srcR = RunParams%FluxSources(J)%Radius
               if (((x - srcX)*(x - srcX) + (y - srcY)*(y - srcY)) < srcR*srcR) then
                  nFluxSeries = RunParams%FluxSources(J)%nFluxSeries
                  ! If just one element in the series, set constant flux that
                  ! applies after the given time.
                  if (nFluxSeries == 1) then
                     if (t < RunParams%FluxSources(J)%time(1) .or. &
                         (t == RunParams%FluxSources(J)%time(1) .and. &
                         grid%t < t)) then
                        Qf(J) = 0.0_wp
                        rhoQf(J) = 0.0_wp
                        psiQf(J) = 0.0_wp
                     else
                        Qf(J) = RunParams%FluxSources(J)%flux(1)
                        rho = Density(RunParams, RunParams%FluxSources(J)%psi(1))
                        rhoQf(J) = rho * Qf(J)
                        psiQf(J) = RunParams%FluxSources(J)%psi(1) * Qf(J)
                     end if
                  else
                     ! If we are outside the time series, flux is zero. The
                     ! final two conditions are needed to carefully treat any
                     ! discontinuities at the series' end points. For example,
                     ! at the end of the series (t = tn, say), we need to have
                     ! nonzero flux when grid%t < tn and zero flux afterwards.
                     ! If t = tn, we must distinguish between the multistep
                     ! method computing intermed3(t+dt) (which needs Q != 0) and
                     ! intermed0(t) at the next step (at which point Q = 0). In
                     ! the latter case, grid%t = t, which is what we're checking
                     ! for in the last condition. Similar reasoning applies for
                     ! the penultimate condition and the test above when
                     ! nFluxSeries == 1.
                     if (t < RunParams%FluxSources(J)%time(1) .or. &
                         t > RunParams%FluxSources(J)%time(nFluxSeries) .or. &
                         (t == RunParams%FluxSources(J)%time(1) .and. &
                         grid%t < t) .or. &
                         (t == RunParams%FluxSources(J)%time(nFluxSeries) .and. &
                         grid%t == t)) then
                        Qf(J) = 0.0_wp
                        rhoQf(J) = 0.0_wp
                        psiQf(J) = 0.0_wp
                     else
                        do K = 2, nFluxSeries
                           if ((t >= RunParams%FluxSources(J)%time(K - 1)) .and. &
                               (t <= RunParams%FluxSources(J)%time(K))) then
                              Qa = RunParams%FluxSources(J)%flux(K - 1)
                              psia = RunParams%FluxSources(J)%psi(K - 1)
                              ta = RunParams%FluxSources(J)%time(K - 1)
                              Qb = RunParams%FluxSources(J)%flux(K)
                              psib = RunParams%FluxSources(J)%psi(K)
                              tb = RunParams%FluxSources(J)%time(K)
                              tfrac = (t - ta) / (tb - ta)
                              Qf(J) = Qa + (Qb - Qa) * tfrac
                              psif = psia + (psib - psia) * tfrac
                              rho = Density(RunParams, psif)
                              rhoQf(J) = rho * Qf(J)
                              psiQf(J) = psif * Qf(J)
                           end if
                        end do
                     end if
                  end if
               end if
               if (RunParams%isOneD) then
                  sourceProp = RunParams%FluxSources(J)%NumCellsInSrc * grid%deltaX
                  Qf(J) = Qf(J)/sourceProp
                  rhoQf(J) = rhoQf(J)/sourceProp
                  psiQf(J) = psiQf(J)/sourceProp
               else
                  sourceProp = RunParams%FluxSources(J)%NumCellsInSrc * grid%deltaX * grid%deltaY
                  Qf(J) = Qf(J)/sourceProp
                  rhoQf(J) = rhoQf(J)/sourceProp
                  psiQf(J) = psiQf(J)/sourceProp
               end if
            end do
            Qt = Qt + sum(Qf)
            rhoQt = rhoQt + sum(rhoQf)
            psiQt = psiQt + sum(psiQf)
         end if
         if (allocated(Qf)) deallocate (Qf)
         if (allocated(rhoQf)) deallocate (rhoQf)
         if (allocated(psiQf)) deallocate (psiQf)
      end if

      dbdx = uvect(RunParams%Vars%dbdx)
      dbdy = uvect(RunParams%Vars%dbdy)

      ! Add on total fluxes of total volume and solid volume. These are subject
      ! to a geometric correction if gam != 1.
      stvect(iw) = stvect(iw) + Qt / (gam * gam)
      stvect(iHnpsi) = stvect(iHnpsi) + psiQt / gam

      if (uvect(iHn) > RunParams%heightThreshold) then

         ! Determine the gravitational forcing and add it on to the momentum
         ! equations.
         g = RunParams%g
         hp_o_gam = -uvect(ibt)
         hp_o_gam = hp_o_gam + (uvect(iw) - uvect(ib0))
         hp_o_gam = hp_o_gam / gam
         stvect(irhoHnu) = stvect(irhoHnu) - g*uvect(irho)*hp_o_gam*dbdx
         stvect(irhoHnv) = stvect(irhoHnv) - g*uvect(irho)*hp_o_gam*dbdy

         if (RunParams%curvature) then
            ! Determine the curvature terms and add it on to the momentum
            ! equations.
            d2bdxx = uvect(RunParams%Vars%d2bdxx)
            d2bdyy = uvect(RunParams%Vars%d2bdyy)
            d2bdxy = uvect(RunParams%Vars%d2bdxy)

            u = uvect(iu)
            if (RunParams%isOneD) then
               v = 0.0_wp
            else
               v = uvect(iv)
            end if

            if (RunParams%MorphodynamicsOn) then
               d2bdtx = uvect(RunParams%Vars%d2bdtx)
               if (RunParams%isOneD) then
                  d2bdty = 0.0_wp
               else
                  d2bdty = uvect(RunParams%Vars%d2bdty)
               end if
            else
               d2bdtx = 0.0_wp
               d2bdty = 0.0_wp
            end if

            curvatureTerm = uvect(irho)*hp_o_gam*(u*u*d2bdxx + 2.0_wp*u*v*d2bdxy + v*v*d2bdyy + u*d2bdtx + v*d2bdty)
            stvect(irhoHnu) = stvect(irhoHnu) - curvatureTerm * dbdx
            stvect(irhoHnv) = stvect(irhoHnv) - curvatureTerm * dbdy
        end if
    end if

   end subroutine ExplicitSourceTerms

   ! Compute (hydraulic) implicit source terms for each of the governing
   ! equations and store in stvect. In our current formulation, this is just the
   ! basal drag. Due to the way the implicit time stepping works [cf Chertock,
   ! Cui, Kurganov & Wu (2015)] this must be supplied as F / (Hn |u|), where u
   ! is the velocity vector and F is friction = drag / density.
   pure subroutine ImplicitSourceTerms(RunParams, uvect, friction, stvect)

      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp), intent(in) :: friction
      real(kind=wp), dimension(:), intent(out) :: stvect(size(uvect))

      integer :: irhoHnu, irhoHnv, iHn

      real(kind=wp) :: Hn, modu, hr
      real(kind=wp) :: rhs

      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHn = RunParams%Vars%Hn

      Hn = uvect(iHn)

      stvect(:) = 0.0_wp

      if (Hn > RunParams%heightThreshold) then
         
         modu = sqrt(FlowSquaredSpeedSlopeAligned(RunParams, uvect))

         if (modu > 1.0e-8_wp) then
            hr = 1.0_wp / Hn ! hr = 1/Hn calculated using desingularization
            rhs = -friction * hr / modu
            stvect(irhoHnu) = rhs
            stvect(irhoHnv) = rhs
         end if
      end if

   end subroutine ImplicitSourceTerms

end module equations_module
