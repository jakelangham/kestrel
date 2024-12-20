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


! This module contains routines involved in approximating the hydraulic operator
! of the shallow-layer governing equations using finite volumes. We use a
! central upwind scheme, adapted from [Chertock, Cui, Kurganov & Wu, Int. J.
! Numer. Meth. Fluids (2015)] to suit our system.
module hydraulic_rhs_module

   use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
   use set_precision_module, only: wp
   use grid_module, only: GridType, TileType, TileList
   use utilities_module, only: KahanSum
   use runsettings_module, only: RunSet
   use equations_module
   use limiters_module, only: limiter
   use closures_module, only: ComputeHn, GeometricCorrectionFactor, Density, DragClosure, DesingularizeFunc
   use messages_module, only: FatalErrorMessage
   use varstring_module

   implicit none

   private
   public :: CalculateHydraulicRHS
   public :: CalculateLimitedDerivsBoundary
   public :: CalculateLimitedDerivs
   public :: Reconstruct
   public :: CorrectSlopes
   public :: ComputeDesingularisedVariables
   public :: CalculateFluxes

contains

   ! Compute the hydraulic part of the right-hand side of the governing
   ! equations. (i.e. all non-morphodynamic terms that are not time
   ! derivatives.) The solution data over which the hydraulic operator H is to be
   ! applied is passed in the tileWorkspace array of tiles. All the calculations
   ! are carried out within tileWorkspace, which contains auxiliary storage for
   ! all the intermediate components involved in the calculation.
   !
   ! The procedure is as follows: we loop over the active tiles in tileWorkspace
   ! and compute each stage of the central-upwind discretisation of H in turn.
   ! This furnishes us with the flow speed, which we use to calculate a stable
   ! time step (in a way that depends on the particular substep of the time
   ! stepping scheme), output as advisedTimeStep. Finally, the various
   ! components of the RHS are summed in ConstructHydraulicRHS.
   subroutine CalculateHydraulicRHS(RunParams, grid, tileWorkspace, t, &
                                    substep, advisedTimeStep)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      type(TileType), dimension(:), intent(inout) :: tileWorkspace
      real(kind=wp), intent(in) :: t
      integer, intent(in) :: substep
      real(kind=wp), intent(out) :: advisedTimeStep

      type(TileList), pointer :: ActiveTiles
      real(kind=wp), dimension(:) :: unitCFLTimeStep(grid%nTiles)

      integer :: tt, ttk

      activeTiles => grid%activeTiles

      unitCFLTimeStep(:) = RunParams%maxdt

      ! Calculate limited spatial derivatives of w, rhoHnu, rhoHnv, Hnpsi.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call CalculateLimitedDerivs(RunParams, grid, RunParams%iFlux, tileWorkspace, ttk)
      end do

      ! Reconstruct cell boundary values of w, rhoHnu, rhoHnv, Hnpsi.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call Reconstruct(RunParams, grid, RunParams%iFlux, tileWorkspace, ttk)
      end do

      ! Correct slopes for w, Hnpsi to ensure positivity.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call CorrectSlopes(RunParams, grid, tileWorkspace, ttk)
      end do

      ! Compute derived variables Hn, u, v, psi and rho at cell centres.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call ComputeDesingularisedVariables(RunParams, tileWorkspace, ttk, .true.)
      end do

      ! Calculate limited spatial derivatives of Hn, u, v, psi, rho.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call CalculateLimitedDerivs(RunParams, grid, RunParams%iDesing, tileWorkspace, ttk)
      end do

      ! Reconstruct cell boundary values of Hn, u, v, psi, rho.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call Reconstruct(RunParams, grid, RunParams%iDesing, tileWorkspace, ttk)
      end do

      ! Calculate numerical flux terms, and also local propagation speeds which
      ! give the CFL time step.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call CalculateFluxes(RunParams, grid, RunParams%iFlux, tileWorkspace, &
                              ttk, unitCFLTimeStep(ttk))  ! Look at eqn input
      end do

      advisedTimeStep = ComputeAdvisedTimeStep(RunParams, grid, &
                                               substep, unitCFLTimeStep)

      ! Sum flux and source terms to complete the calculation.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call ConstructHydraulicRHS(RunParams, grid, tileWorkspace, ttk, t)
      end do
   end subroutine CalculateHydraulicRHS

   ! This contains all the logic for computing a suitable time step, i.e. one
   ! that is stable with respect to the CFL condition and the diffusive time
   ! scale set by eddy viscosity.
   function ComputeAdvisedTimeStep(RunParams, grid, substep, &
                                   unitCFLTimeStep) result(advisedTimeStep)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: substep ! i.e. stage of multistep timestepper
      real(kind=wp), dimension(:) :: unitCFLTimeStep(grid%nTiles)

      real(kind=wp) maxTimeStep, advisedTimeStep

      ! N.B. If eddy viscosity is present we need to limit the time step by the 
      ! diffusive time scale ~ dx^2 / nu
      maxTimeStep = min(RunParams%cfl * minval(unitCFLTimeStep), &
                        RunParams%diffusiveTimeScale)

      maxTimeStep = min(maxTimeStep, RunParams%maxdt)

      if (substep == 1) then
         ! Reasoning for this: sometimes as flow is speeding up, the
         ! requisite time step gets progressively smaller leading to a situation
         ! where the second call to H always wants to refine the time step,
         ! leading to a costly reset. By making the initial step more conservative
         ! the hope is to avoid this problem.
         advisedTimeStep = 0.9_wp * maxTimeStep

         ! Since this is the first substep, the advised step is the step that will
         ! actually be taken, unless overridden in later substeps.
         grid%dt = advisedTimeStep
      else
         advisedTimeStep = maxTimeStep
      end if

   end function ComputeAdvisedTimeStep

   ! Compute limited derivatives in x and y for the given variables, indexed by
   ! the variables array, over the tile indexed by tID. Derivatives in x and y
   ! are stored in the uLimX and uLimY arrays respectively. The limiter function
   ! (which is user-settable) tries to prevent oscillations whose lengthscales
   ! are comparable to deltaX, deltaY.
   subroutine CalculateLimitedDerivs(RunParams, grid, variables, tiles, tID)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(in) :: grid
      integer, dimension(:), intent(in) :: variables
      type(TileType), target, dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID

      real(kind=wp), dimension(:,:,:), pointer :: u

      integer :: i, j, k, d, nVars, nYPoints, nXPoints, neighbour

      nXPoints = RunParams%nXpertile
      nYPoints = RunParams%nYpertile
      nVars = size(variables)

      u => tiles(tID)%u

      ! Over all y-points and x-points, excluding East and West boundaries,
      ! calculate limited X-derivatives.
      do j = 1, nYPoints
         do i = 2, nXPoints - 1
            do k = 1, nVars
               d = variables(k)
               tiles(tID)%uLimX(d, i, j) = grid%deltaXRecip *  &
                  limiter(u(d, i + 1, j) - u(d, i, j),  &
                  u(d, i, j) - u(d, i - 1, j))
            end do
         end do
      end do

      ! Over all x-points and y-points, excluding North and South boundaries,
      ! calculate limited Y-derivatives.
      do j = 2, nYPoints - 1
         do i = 1, nXPoints
            do k = 1, nVars
               d = variables(k)
               tiles(tID)%uLimY(d, i, j) = grid%deltaYRecip *  &
                  limiter(u(d, i, j + 1) - u(d, i, j),  &
                  u(d, i, j) - u(d, i, j - 1))
            end do
         end do
      end do

      ! Now fill in the boundaries.
      if (.not. RunParams%isOneD) then
         ! Y-derivs on South boundary
         j = 1
         do i = 1, nXPoints
            do k = 1, nVars
               d = variables(k)
               call CalculateLimitedDerivsBoundary(RunParams, grid, d, i, j,  &
                  tiles, tID, 'S')
            end do
         end do

         ! Do south ghost tile too if it exists
         neighbour = grid%tileContainer(tID)%South
         if (grid%tileContainer(neighbour)%TileOn) then
            j = nYpoints
            do i = 1, nXPoints
               do k = 1, nVars
                  d = variables(k)
                  call CalculateLimitedDerivsBoundary(RunParams, grid, d, i, j,  &
                     tiles, neighbour, 'N')
               end do
            end do
         end if

         ! Y-derivs on North boundary
         j = nYPoints
         do i = 1, nXPoints
            do k = 1, nVars
               d = variables(k)
               call CalculateLimitedDerivsBoundary(RunParams, grid, d, i, j,  &
                  tiles, tID, 'N')
            end do
         end do

         ! Do north ghost tile too if it exists
         neighbour = grid%tileContainer(tID)%North
         if (grid%tileContainer(neighbour)%TileOn) then
            j = 1
            do i = 1, nXpoints
               do k = 1, nVars
                  d = variables(k)
                  call CalculateLimitedDerivsBoundary(RunParams, grid, d, i, j,  &
                     tiles, neighbour, 'S')
               end do
            end do
         end if
      end if

      ! X-derivs on West boundary
      i = 1
      do j = 1, nYPoints
         do k = 1, nVars
            d = variables(k)
            call CalculateLimitedDerivsBoundary(RunParams, grid, d, i, j,  &
               tiles, tID, 'W')
         end do
      end do

      ! Do West ghost tile too if it exists
      neighbour = grid%tileContainer(tID)%West
      if (grid%tileContainer(neighbour)%TileOn) then
         i = nXpoints
         do j = 1, nYpoints
            do k = 1, nVars
               d = variables(k)
               call CalculateLimitedDerivsBoundary(RunParams, grid, d, i, j,  &
                  tiles, neighbour, 'E')
            end do
         end do
      end if

      ! X-derivs on East boundary
      i = nXPoints
      do j = 1, nYPoints
         do k = 1, nVars
            d = variables(k)
            call CalculateLimitedDerivsBoundary(RunParams, grid, d, i, j,  &
               tiles, tID, 'E')
         end do
      end do

      ! Do East ghost tile too if it exists
      neighbour = grid%tileContainer(tID)%East
      if (grid%tileContainer(neighbour)%TileOn) then
         i = 1
         do j = 1, nYpoints
            do k = 1, nVars
               d = variables(k)
               call CalculateLimitedDerivsBoundary(RunParams, grid, d, i, j,  &
                  tiles, neighbour, 'W')
            end do
         end do
      end if

   end subroutine CalculateLimitedDerivs

   ! Compute the limited derivatives at the boundary of the tile (indexed by
   ! tID) for the variable d and at the (i, j)-th point. The variable bid passes 
   ! the direction of the boundary (NSEW).
   subroutine CalculateLimitedDerivsBoundary(RunParams, grid, d, i, j, tiles, &
                                             tID, bid)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(in) :: grid
      integer, intent(in) :: d, i, j
      type(TileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID
      character(len=1), intent(in) :: bid

      integer :: nextTile = 0, prevTile = 0

      real(kind=wp) :: uN, uS, uE, uW
      real(kind=wp) :: deltaXRecip, deltaYRecip

      deltaXRecip = grid%deltaXRecip
      deltaYRecip = grid%deltaYRecip

      select case (bid)
       case ('N')
         if (.not. RunParams%isOneD) then
            nextTile = tiles(tID)%North
            uN = tiles(nextTile)%u(d, i, 1)
            tiles(tID)%uLimY(d, i, j) = deltaYRecip *  &
               limiter(uN - tiles(tID)%u(d, i, j),   &
               tiles(tID)%u(d, i, j) -  &
               tiles(tID)%u(d, i, j - 1))
         else
            tiles(tID)%uLimY(d, i, j) = 0.0_wp
         end if
       case ('S')
         if (.not. RunParams%isOneD) then
            prevTile = tiles(tID)%South
            uS = tiles(prevTile)%u(d, i, RunParams%nYpertile)
            tiles(tID)%uLimY(d, i, j) = deltaYRecip *  &
               limiter(tiles(tID)%u(d, i, j + 1) -  &
               tiles(tID)%u(d, i, j),  &
               tiles(tID)%u(d, i, j) - uS)
         else
            tiles(tID)%uLimY(d, i, j) = 0.0_wp
         end if
       case ('E')
         nextTile = tiles(tID)%East
         uE = tiles(nextTile)%u(d, 1, j)
         tiles(tID)%uLimX(d, i, j) = deltaXRecip *  &
            limiter(uE - tiles(tID)%u(d, i, j),  &
            tiles(tID)%u(d, i, j) -  &
            tiles(tID)%u(d, i - 1, j))
       case ('W')
         prevTile = tiles(tID)%West
         uW = tiles(prevTile)%u(d, RunParams%nXpertile, j)
         tiles(tID)%uLimX(d, i, j) = deltaXRecip *  &
            limiter(tiles(tID)%u(d, i + 1, j) -  &
            tiles(tID)%u(d, i, j),  &
            tiles(tID)%u(d, i, j) - uW)
       case default
         call FatalErrorMessage('bid not recognized')
      end select

   end subroutine CalculateLimitedDerivsBoundary

   ! Reconstruct solution fields at the edges of cells, given their values at
   ! cell centres. The solution data to be reconstructed is assumed to be in
   ! tiles(tID)%u(dims,:,:), where dims is an indexing array. Let (x,y) denote
   ! the centre of cell (i,j) and let (dx,dy) be the numerical grid spacings. 
   ! Then the reconstructions are stored at:
   ! tiles(tID)%uMinusX(d, i, j), for (x - dx/2, y)
   ! tiles(tID)%uPlusX(d, i, j),  for (x + dx/2, y)
   ! tiles(tID)%uMinusY(d, i, j), for (x, y - dy/2)
   ! tiles(tID)%uPlusY(d, i, j),  for (x, y + dy/2), for each d in dims.
   subroutine Reconstruct(RunParams, grid, dims, tiles, tID)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout) :: grid
      integer, dimension(:), intent(in) :: dims
      type(TileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID 

      integer :: i, j, k, d, nd
      integer :: iw, irhoHnu, irhoHnv, iHnpsi
      integer :: iHn, iu, iv, ipsi, irho, ib0, ibt

      integer :: nextTile, prevTile

      integer :: nXPoints, nYPoints
      real(kind=wp) :: deltaX, deltaY
      real(kind=wp) :: gam

      nXPoints = RunParams%nXpertile
      nYPoints = RunParams%nYpertile
      nd = size(dims)

      deltaX = grid%deltaX
      deltaY = grid%deltaY

      iw = RunParams%Vars%w
      ib0 = RunParams%Vars%b0
      ibt = RunParams%Vars%bt
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      irho = RunParams%Vars%rho
      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u
      iv = RunParams%Vars%v
      ipsi = RunParams%Vars%psi

      ! Calculate reconstructed variables at each cell edge
      prevTile = tiles(tID)%West
      do k = 1, nd
         d = dims(k)
         tiles(tID)%uPlusX(d, 1, :) = tiles(tID)%u(d, 1, :) - &
            tiles(tID)%uLimX(d, 1, :) * 0.5_wp * deltaX
         tiles(tID)%uMinusX(d, 1, :) = tiles(prevTile)%u(d, nXPoints, :) + &
            tiles(prevTile)%uLimX(d, nXPoints, :) * 0.5_wp * deltaX
      end do
      do k = 1, nd
         d = dims(k)
         tiles(tID)%uPlusX(d, 2:nXPoints, :) = tiles(tID)%u(d, 2:nXPoints, :) - &
            tiles(tID)%uLimX(d, 2:nXPoints, :) * 0.5_wp * deltaX
         tiles(tID)%uMinusX(d, 2:nXPoints, :) = tiles(tID)%u(d, 1:nXPoints-1, :) + &
            tiles(tID)%uLimX(d, 1:nXPoints-1, :) * 0.5_wp * deltaX
      end do

      nextTile = tiles(tID)%East
      do k = 1, nd
         d = dims(k)
         tiles(tID)%uPlusX(d, nXPoints+1, :) = tiles(nextTile)%u(d, 1, :) - &
            tiles(nextTile)%uLimX(d, 1, :) * 0.5_wp * deltaX
         tiles(tID)%uMinusX(d, nXPoints+1, :) = tiles(tID)%u(d, nXPoints, :) + &
            tiles(tID)%uLimX(d, nXPoints, :) * 0.5_wp * deltaX
      end do

      if (.not. RunParams%isOneD) then
         prevTile = tiles(tID)%South
         do k = 1, nd
            d = dims(k)
            tiles(tID)%uPlusY(d, :, 1) = tiles(tID)%u(d, :, 1) - &
               tiles(tID)%uLimY(d, :, 1) * 0.5_wp * deltaY
            tiles(tID)%uMinusY(d, :, 1) = tiles(prevTile)%u(d, :, nYPoints) + &
               tiles(prevTile)%uLimY(d, :, nYPoints) * 0.5_wp * deltaY
         end do
         do k = 1, nd
            d = dims(k)
            tiles(tID)%uPlusY(d, :, 2:nYPoints) = tiles(tID)%u(d, :, 2:nYPoints) - &
               tiles(tID)%uLimY(d, :, 2:nYPoints) * 0.5_wp * deltaY
            tiles(tID)%uMinusY(d, :, 2:nYPoints) = tiles(tID)%u(d, :, 1:nYPoints-1) + &
               tiles(tID)%uLimY(d, :, 1:nYPoints-1) * 0.5_wp * deltaY
         end do

         nextTile = tiles(tID)%North
         do k = 1, nd
            d = dims(k)
            tiles(tID)%uPlusY(d, :, nYPoints+1) = tiles(nextTile)%u(d, :, 1) - &
               tiles(nextTile)%uLimY(d, :, 1) * 0.5_wp * deltaY
            tiles(tID)%uMinusY(d, :, nYPoints+1) = tiles(tID)%u(d, :, nYPoints) + &
               tiles(tID)%uLimY(d, :, nYPoints) * 0.5_wp * deltaY
         end do
      end if

      if (all(dims(1:4) .eq. RunParams%iFlux)) then
         return
      end if

      ! Don't 'reconstruct' Hn - we require that hp = w - b exactly at cell
      ! vertices (as well as centres) for numerical lake at rest to work.
      do i = 1, RunParams%nXpertile + 1
         do j = 1, RunParams%nYpertile
            ! n.b. uPlusX(b,i,j) = uMinusX(b,i,j), so it doesn't matter which
            ! side gam is computed from
            gam = GeometricCorrectionFactor(RunParams, tiles(tID)%uPlusX(:,i,j))
            tiles(tID)%uPlusX(iHn,i,j) = ComputeHn(tiles(tID)%uPlusX(iw,i,j), &
               tiles(tID)%uPlusX(ib0,i,j), tiles(tID)%uPlusX(ibt,i,j), gam)
            gam = GeometricCorrectionFactor(RunParams, tiles(tID)%uMinusX(:,i,j))
            tiles(tID)%uMinusX(iHn,i,j) = ComputeHn(tiles(tID)%uMinusX(iw,i,j), &
               tiles(tID)%uMinusX(ib0,i,j), tiles(tID)%uMinusX(ibt,i,j), gam)
         end do
      end do
      if (.not. RunParams%isOneD) then
         do i = 1, RunParams%nXpertile
            do j = 1, RunParams%nYpertile + 1
               gam = GeometricCorrectionFactor(RunParams, tiles(tID)%uPlusY(:,i,j))
               tiles(tID)%uPlusY(iHn,i,j) = ComputeHn(tiles(tID)%uPlusY(iw,i,j), &
                  tiles(tID)%uPlusY(ib0,i,j), tiles(tID)%uPlusY(ibt,i,j), gam)
               gam = GeometricCorrectionFactor(RunParams, tiles(tID)%uMinusY(:,i,j))
               tiles(tID)%uMinusY(iHn,i,j) = ComputeHn(tiles(tID)%uMinusY(iw,i,j), &
                  tiles(tID)%uMinusY(ib0,i,j), tiles(tID)%uMinusY(ibt,i,j), gam)
            end do
         end do
      end if

      ! Recompute reconstructions of flux variables so that they are
      ! commensurate with the corresponding reconstructions of H, u, v, etc.
      tiles(tID)%uPlusX(irhoHnu,:,:) = tiles(tID)%uPlusX(irho,:,:) *  &
         tiles(tID)%uPlusX(iHn,:,:) * tiles(tID)%uPlusX(iu,:,:)
      tiles(tID)%uMinusX(irhoHnu,:,:) = tiles(tID)%uMinusX(irho,:,:) *  &
         tiles(tID)%uMinusX(iHn,:,:) * tiles(tID)%uMinusX(iu,:,:)
      ! Could recalculate Hnpsi in the same way here, but to maintain Hnpsi > 0
      ! requires (cf Kurganov & Petrova, 2007)
      ! u(iHnpsi,i,j) = 0.5 * (uPlusX(iHnpsi,i,j) + uMinusX(iHnpsi,i+1,j)),
      ! which is true by construction but can be violated if we recompute.
!      tiles(tID)%uPlusX(iHnpsi,:,:) = tiles(tID)%uPlusX(ipsi,:,:) *  &
!         tiles(tID)%uPlusX(iHn,:,:)
!      tiles(tID)%uMinusX(iHnpsi,:,:) = tiles(tID)%uMinusX(ipsi,:,:) *  &
!         tiles(tID)%uMinusX(iHn,:,:)
      if (.not. RunParams%isOneD) then
        tiles(tID)%uPlusY(irhoHnu,:,:) = tiles(tID)%uPlusY(irho,:,:) *  &
           tiles(tID)%uPlusY(iHn,:,:) * tiles(tID)%uPlusY(iu,:,:)
        tiles(tID)%uMinusY(irhoHnu,:,:) = tiles(tID)%uMinusY(irho,:,:) *  &
           tiles(tID)%uMinusY(iHn,:,:) * tiles(tID)%uMinusY(iu,:,:)
        tiles(tID)%uPlusX(irhoHnv,:,:) = tiles(tID)%uPlusX(irho,:,:) *  &
           tiles(tID)%uPlusX(iHn,:,:) * tiles(tID)%uPlusX(iv,:,:)
        tiles(tID)%uMinusX(irhoHnv,:,:) = tiles(tID)%uMinusX(irho,:,:) *  &
           tiles(tID)%uMinusX(iHn,:,:) * tiles(tID)%uMinusX(iv,:,:)
        tiles(tID)%uPlusY(irhoHnv,:,:) = tiles(tID)%uPlusY(irho,:,:) *  &
           tiles(tID)%uPlusY(iHn,:,:) * tiles(tID)%uPlusY(iv,:,:)
        tiles(tID)%uMinusY(irhoHnv,:,:) = tiles(tID)%uMinusY(irho,:,:) *  &
           tiles(tID)%uMinusY(iHn,:,:) * tiles(tID)%uMinusY(iv,:,:)
        ! (cf above comment)
!        tiles(tID)%uPlusY(iHnpsi,:,:) = tiles(tID)%uPlusY(ipsi,:,:) *  &
!           tiles(tID)%uPlusY(iHn,:,:)
!        tiles(tID)%uMinusY(iHnpsi,:,:) = tiles(tID)%uMinusY(ipsi,:,:) *  &
!           tiles(tID)%uMinusY(iHn,:,:)
      end if
   end subroutine Reconstruct

   ! This routine applies a correction to the default reconstruction of w and
   ! Hnpsi at cell interfaces (see above) in order to ensure that no
   ! reconstructed points end up with w or Hnpsi < 0. We follow the well
   ! balanced approach of Chertock et al. (2015), though the logic ends up being
   ! rather complicated in order to account for tile boundaries and the need to
   ! make sure that calculations are independent of the tiling layout.
   subroutine CorrectSlopes(RunParams, grid, tiles, tID)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout) :: grid
      type(TileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID ! tileID

      integer :: i, j, prevTile, nextTile

      integer :: nXPoints, nYPoints
      integer :: iw, iHnpsi, ib0, ibt

      real(kind=wp) :: b_left, b_right, b_bottom, b_top

      nXPoints = RunParams%nXpertile
      nYPoints = RunParams%nYpertile
      iw = RunParams%Vars%w
      iHnpsi = RunParams%Vars%Hnpsi
      ib0 = RunParams%Vars%b0
      ibt = RunParams%Vars%bt

      ! First check the 'horizontal' reconstructions, ie uMinusX/uPlusX.
      do j = 1, nYPoints

         ! Left tile boundary
         prevTile = tiles(tID)%West
         if (.not. RunParams%isOneD) then
            ! If 2D, have to interpolate b_left/right in middle of cell interface
            b_left = InterpolateB(tiles(prevTile)%b0(nXpoints,j  ), tiles(prevTile)%b0(nXpoints,j+1), &
                                  tiles(prevTile)%bt(nXpoints,j  ), tiles(prevTile)%bt(nXpoints,j+1))
            b_right = InterpolateB(tiles(prevTile)%b0(nXpoints+1,j  ), tiles(prevTile)%b0(nXpoints+1,j+1), &
                                   tiles(prevTile)%bt(nXpoints+1,j  ), tiles(prevTile)%bt(nXpoints+1,j+1))
         else
            b_left = tiles(prevTile)%b0(nXpoints,1) + tiles(prevTile)%bt(nXpoints,1)
            b_right = tiles(prevTile)%b0(nXpoints+1,1) + tiles(prevTile)%bt(nXpoints+1,1)
         end if
         if ((tiles(prevTile)%uMinusX(iw,nXpoints+1,j) < b_right) .or.  &
            (tiles(prevTile)%uPlusX(iw,nXpoints,j) < b_left)) then
            tiles(tID)%uMinusX(iw,1,j) = tiles(prevTile)%u(iw,nXpoints,j) +  &
               0.5_wp * (b_right - b_left)
         else if (tiles(tID)%uMinusX(iw,1,j) /= tiles(prevTile)%uMinusX(iw,nXpoints+1,j) &
            .and. grid%tileContainer(prevTile)%isGhostTile .eqv. .false.) then
            ! This case covers the fact that if the left tile already corrected
            ! its uMinusX value, we need to make sure both copies at the tile
            ! interface are the same (since the previous test will fail to flag it)
            tiles(tID)%uMinusX(iw,1,j) = tiles(prevTile)%uMinusX(iw,nXpoints+1,j)
         end if

         if (tiles(tID)%uMinusX(iHnpsi,1,j) < 0.0_wp) then
            tiles(tID)%uMinusX(iHnpsi,1,j) = tiles(prevTile)%u(iHnpsi,nXPoints,j)
         end if

         do i = 1, nXPoints
            if (.not. RunParams%isOneD) then
               ! If 2D, have to interpolate b_left/right in middle of cell interface
               b_left = InterpolateB(tiles(tID)%b0(i,j  ), tiles(tID)%b0(i,j+1), &
                                     tiles(tID)%bt(i,j  ), tiles(tID)%bt(i,j+1))
               b_right = InterpolateB(tiles(tID)%b0(i+1,j  ), tiles(tID)%b0(i+1,j+1), &
                                      tiles(tID)%bt(i+1,j  ), tiles(tID)%bt(i+1,j+1))
            else
               b_left = tiles(tID)%b0(i,1) + tiles(tID)%bt(i,1)
               b_right = tiles(tID)%b0(i+1,1) + tiles(tID)%bt(i+1,1)
            end if

            if ((tiles(tID)%uMinusX(iw,i+1,j) < b_right) .or.  &
               (tiles(tID)%uPlusX(iw,i,j) < b_left)) then
               ! Replace reconstruction by setting w gradient equal to bed gradient
               tiles(tID)%uMinusX(iw,i+1,j) = tiles(tID)%u(iw,i,j) +  &
                  0.5_wp * (b_right - b_left)
               tiles(tID)%uPlusX(iw,i,j) = tiles(tID)%u(iw,i,j) +  &
                  0.5_wp * (b_left - b_right)
            end if

            ! Force reconstruction to have non-negative height at cell edges
            if ((tiles(tID)%uMinusX(iHnpsi,i+1,j) < 0.0_wp) .or.  &
               (tiles(tID)%uPlusX(iHnpsi,i,j) < 0.0_wp)) then
               tiles(tID)%uMinusX(iHnpsi,i+1,j) = tiles(tID)%u(iHnpsi,i,j)
               tiles(tID)%uPlusX(iHnpsi,i,j) = tiles(tID)%u(iHnpsi,i,j)
            end if
         end do

         ! Right tile boundary
         nextTile = tiles(tID)%East
         if (.not. RunParams%isOneD) then
            ! if 2D, have to interpolate b_left/right in middle of cell interface
            b_left = InterpolateB(tiles(nextTile)%b0(1,j  ), tiles(nextTile)%b0(1,j+1), &
                                  tiles(nextTile)%bt(1,j  ), tiles(nextTile)%bt(1,j+1))
            b_right = InterpolateB(tiles(nextTile)%b0(2,j  ), tiles(nextTile)%b0(2,j+1), &
                                   tiles(nextTile)%bt(2,j  ), tiles(nextTile)%bt(2,j+1))
         else
            b_left = tiles(nextTile)%b0(1,1) + tiles(nextTile)%bt(1,1)
            b_right = tiles(nextTile)%b0(2,1) + tiles(nextTile)%bt(2,1)
         end if

         if ((tiles(nextTile)%uMinusX(iw,2,j) < b_right) .or.  &
            (tiles(nextTile)%uPlusX(iw,1,j) < b_left)) then
            tiles(tID)%uPlusX(iw,nXpoints+1,j) = tiles(nextTile)%u(iw,1,j) +  &
               0.5_wp * (b_left - b_right)
         else if (tiles(tID)%uPlusX(iw,nXpoints+1,j) /= tiles(nextTile)%uPlusX(iw,1,j) &
            .and. grid%tileContainer(nextTile)%isGhostTile .eqv. .false.) then
            ! Likewise, this case covers the fact that if the right tile already
            ! corrected its uPlusX value, we need to make sure both copies at the tile
            ! interface are the same (since the previous test will fail to flag it).
            tiles(tID)%uPlusX(iw,nXpoints+1,j) = tiles(nextTile)%uPlusX(iw,1,j)
         end if
         if (tiles(tID)%uPlusX(iHnpsi,nXPoints+1,j) < 0.0_wp) then
            tiles(tID)%uPlusX(iHnpsi,nXPoints+1,j) = tiles(nextTile)%u(iHnpsi,1,j)
         end if
      end do

      if (.not. RunParams%isOneD) then
         do i = 1, nXPoints

            ! Bottom tile boundary
            prevTile = tiles(tID)%South
            b_bottom = InterpolateB(tiles(prevTile)%b0(i  ,nYpoints), tiles(prevTile)%b0(i+1,nYpoints), &
                                    tiles(prevTile)%bt(i  ,nYpoints), tiles(prevTile)%bt(i+1,nYpoints))
            b_top = InterpolateB(tiles(prevTile)%b0(i  ,nYpoints+1), tiles(prevTile)%b0(i+1,nYpoints+1), &
                                 tiles(prevTile)%bt(i  ,nYpoints+1), tiles(prevTile)%bt(i+1,nYpoints+1))
            if ((tiles(prevTile)%uMinusY(iw,i,nYpoints+1) < b_top) .or.  &
               (tiles(prevTile)%uPlusY(iw,i,nYpoints) < b_bottom)) then
               tiles(tID)%uMinusY(iw,i,1) = tiles(prevTile)%u(iw,i,nYpoints) +  &
                  0.5_wp * (b_top - b_bottom)
            else if (tiles(tID)%uMinusY(iw,i,1) /= tiles(prevTile)%uMinusY(iw,i,nYpoints+1) &
               .and. grid%tileContainer(prevTile)%isGhostTile .eqv. .false.) then
               tiles(tID)%uMinusY(iw,i,1) = tiles(prevTile)%uMinusY(iw,i,nYpoints+1)
            end if
            if (tiles(tID)%uMinusY(iHnpsi,i,1) < 0.0_wp) then
               tiles(tID)%uMinusY(iHnpsi,i,1) = tiles(prevTile)%u(iHnpsi,i,nYPoints)
            end if

            ! Interior points
            do j = 1, nYPoints
               b_bottom = InterpolateB(tiles(tID)%b0(i  ,j), tiles(tID)%b0(i+1,j), &
                                       tiles(tID)%bt(i  ,j), tiles(tID)%bt(i+1,j))
               b_top = InterpolateB(tiles(tID)%b0(i  ,j+1), tiles(tID)%b0(i+1,j+1), &
                                    tiles(tID)%bt(i  ,j+1), tiles(tID)%bt(i+1,j+1))
               if ((tiles(tID)%uMinusY(iw,i,j+1) < b_top) .or.  &
                  (tiles(tID)%uPlusY(iw,i,j) < b_bottom)) then
                  ! Replace reconstruction by setting w gradient equal to bed gradient
                  tiles(tID)%uMinusY(iw,i,j+1) = tiles(tID)%u(iw,i,j) +  &
                     0.5_wp * (b_top - b_bottom)
                  tiles(tID)%uPlusY(iw,i,j) = tiles(tID)%u(iw,i,j) +  &
                     0.5_wp * (b_bottom - b_top)
               end if

               if ((tiles(tID)%uMinusY(iHnpsi,i,j+1) < 0.0_wp) .or.  &
                  (tiles(tID)%uPlusY(iHnpsi,i,j) < 0.0_wp)) then
                  tiles(tID)%uMinusY(iHnpsi,i,j+1) = tiles(tID)%u(iHnpsi,i,j)
                  tiles(tID)%uPlusY(iHnpsi,i,j) = tiles(tID)%u(iHnpsi,i,j)
               end if
            end do

            ! Top tile boundary
            nextTile = tiles(tID)%North
            b_bottom = InterpolateB(tiles(nextTile)%b0(i  ,1), tiles(nextTile)%b0(i+1,1), &
                                    tiles(nextTile)%bt(i  ,1), tiles(nextTile)%bt(i+1,1))
            b_top = InterpolateB(tiles(nextTile)%b0(i  ,2), tiles(nextTile)%b0(i+1,2), &
                                 tiles(nextTile)%bt(i  ,2), tiles(nextTile)%bt(i+1,2))
            if ((tiles(nextTile)%uMinusY(iw,i,2) < b_top) .or.  &
               (tiles(nextTile)%uPlusY(iw,i,1) < b_bottom)) then
               tiles(tID)%uPlusY(iw,i,nYpoints+1) = tiles(nextTile)%u(iw,i,1) +  &
                  0.5_wp * (b_bottom - b_top)
            else if (tiles(tID)%uPlusY(iw,i,nYpoints+1) /= tiles(nextTile)%uPlusY(iw,i,1) &
               .and. grid%tileContainer(nextTile)%isGhostTile .eqv. .false.) then
               tiles(tID)%uPlusY(iw,i,nYpoints+1) = tiles(nextTile)%uPlusY(iw,i,1)
            end if

            if (tiles(tID)%uPlusY(iHnpsi,i,nYPoints+1) < 0.0_wp) then
               tiles(tID)%uPlusY(iHnpsi,i,nYPoints+1) = tiles(nextTile)%u(iHnpsi,i,1)
            end if

         end do
      end if

   end subroutine CorrectSlopes

   ! This is a helper function for CorrectSlopes that linearly interpolates b
   ! between two points carefully, making sure to do the addition of the small
   ! magnitude bt parts first.
   pure function InterpolateB(b0_0, b0_1, bt_0, bt_1) result(b)
      implicit none

      real(kind=wp), intent(in) :: b0_0, b0_1
      real(kind=wp), intent(in) :: bt_0, bt_1
      real(kind=wp) :: b

      b = bt_0 + bt_1
      b = b + b0_0
      b = b + b0_1
      b = b * 0.5_wp

   end function InterpolateB

   ! This routine computes derived variables Hn, u, v, psi, rho from the primary
   ! flow variables. It employs 'desingularisation' formulae to prevent
   ! unphysical blow-up when flow depths are lower than
   ! RunParams%heightThreshold.  It is called within both the hydraulic and
   ! morphodynamic parts of the time stepper. The computeVelocities flag exists
   ! so we can exclude u & v computation within the morphodynamic part (within
   ! which u & v are considered constant).
   recursive subroutine ComputeDesingularisedVariables(RunParams, tiles, tID, &
      computeVelocities)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(TileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID
      logical, intent(in) :: computeVelocities

      integer :: i, j, nXPoints, nYPoints

      real(kind=wp) :: rhow, rhos

      real(kind=wp) :: rhoHnu, rhoHnv, Hnpsi
      real(kind=wp) :: rho, Hn, u, v, psi
      real(kind=wp) :: Hneps, gam, Hneps_gam
      
      real(kind=wp) :: Hn_recip ! 1/Hn computed using Desingularization

      integer :: iw, irhoHnu, irhoHnv, iHnpsi
      integer :: iHn, iu, iv, ipsi, ib0, ibt, irho

      nXPoints = RunParams%nXpertile
      nYPoints = RunParams%nYpertile

      Hneps = RunParams%heightThreshold

      rhow = RunParams%rhow
      rhos = RunParams%rhos

      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      iHn = RunParams%Vars%Hn
      iu = RunParams%Vars%u
      iv = RunParams%Vars%v
      ipsi = RunParams%Vars%psi
      irho = RunParams%Vars%rho
      ib0 = RunParams%Vars%b0
      ibt = RunParams%Vars%bt

      do i = 1, nXPoints
         do j = 1, nYPoints
            rhoHnu = tiles(tID)%u(irhoHnu,i,j)
            Hnpsi = tiles(tID)%u(iHnpsi,i,j)

            gam = GeometricCorrectionFactor(RunParams, tiles(tID)%u(:,i,j))
            Hn = ComputeHn(tiles(tID)%u(iw,i,j), tiles(tID)%u(ib0,i,j), tiles(tID)%u(ibt,i,j), gam)

            Hneps_gam = Hneps * gam

            Hn_recip = DesingularizeFunc(Hn,Hneps_gam)

#if DEBUG_NEGATIVE_DEPTH==1 || DEBUG_NEGATIVE_DEPTH==2
            ! In principle this should never happen within the hydraulic
            ! stepper, so the check is only needed for debugging.
            if (Hn < 0.0_wp .and. abs(Hn) > 10.0_wp * epsilon(Hn) .and. computeVelocities) then
               write(stderr, *) "Hn < 0", tID, i, j, Hn, tiles(tID)%u(iw,i,j), &
                  tiles(tID)%u(ib0,i,j), tiles(tID)%u(ibt,i,j), gam
#if DEBUG_NEGATIVE_DEPTH==2
               call exit(1)
#endif
            end if
#endif

#if DEBUG_NEGATIVE_CONC==1 || DEBUG_NEGATIVE_CONC==2
            ! Likewise, this should never happen within the hydraulic 
            ! stepper, so the check is only needed for debugging.
            if (Hnpsi < 0.0_wp .and. abs(Hnpsi) > 10.0_wp * epsilon(Hnpsi) .and. computeVelocities) then
               write(stderr, *) "Hnpsi < 0", tID, i, j, Hnpsi, tiles(tID)%u(RunParams%Vars%psi,i,j), &
                     Hn, tiles(tID)%u(iw,i,j), tiles(tID)%u(ib0,i,j) + tiles(tID)%u(ibt,i,j), gam
#if DEBUG_NEGATIVE_CONC==2
               call exit(1)
#endif
            end if
#endif

            ! Within substeps of the morphodynamic operator, w & Hnpsi are 
            ! permitted to deviate s.t. Hn & psi would be computed as negative,
            ! but then later, either 1) w is 'topped up' by subsequent
            ! substep, 2) corrected by RedistributDeposition or 3) the time
            ! step is cancelled and refined. Either way, we need Hn & Hnpsi to
            ! stay non-negative for purposes of computing E & D.
            if (Hn < 0.0_wp) then
               Hn = 0.0_wp
            end if
            if (Hnpsi < 0.0_wp) then
               Hnpsi = 0.0_wp
            end if

            psi = min(Hnpsi * Hn_recip, RunParams%maxPack)
            rho = Density(RunParams, psi)

# if DEBUG_EXCESS_CONC==1 || DEBUG_EXCESS_CONC==2
            if ((Hn > 0.0_wp) .and. (Hn < Hneps)) then
               if (Hnpsi > 0.5_wp * RunParams%maxPack * Hn * (1.0_wp + Hneps * Hneps / Hn / Hn)) then
                  write(stderr, *) "Hnpsi too big ", tID, i, j, Hnpsi, &
                     0.5_wp * RunParams%maxPack * Hn * (1.0_wp + Hneps*Hneps / (Hn*Hn)), Hn, psi
#if DEBUG_EXCESS_CONC==2
                  call exit
#endif
               end if
            end if
#endif

            tiles(tID)%u(iHn,i,j) = Hn
            tiles(tID)%u(ipsi,i,j) = psi
            tiles(tID)%u(irho,i,j) = rho

            if (computeVelocities) then
               ! u = rhoHnu * (1/Hn) / rho with 1/Hn calculated using desingularization
               u = rhoHnu * Hn_recip / rho

               tiles(tID)%u(iu,i,j) = u

               if (.not. RunParams%isOneD) then
                  rhoHnv = tiles(tID)%u(irhoHnv, i, j)
                  ! v = rhoHnv * (1/Hn) / rho with 1/Hn calculated using desingularization
                  v = rhoHnv * Hn_recip / rho
                  tiles(tID)%u(iv,i,j) = v
               end if
            end if
         end do
      end do
   end subroutine ComputeDesingularisedVariables

   ! This routine determines the numerical approximations to the flux terms in
   ! the governing equations. The values are calculated via central-upwinding
   ! formulae that combine flux reconstructions at each interface
   ! (u{Plus,Minus}*{X,Y}), [pertaining to the distinct (Plus/Minus)
   ! reconstructions from either side], with their corresponding characteristic
   ! wave speeds (ws{Plus,Minus}).  For the full construction, see Chertock et
   ! al (2015) or e.g. Kurganov & Tadmor (2000).
   !
   ! For each equation describing the evolution of one the primary variables, we
   ! ultimately calculate the following numerical fluxes at the edges of each
   ! cell in tile tID: 
   ! h{X,Y}flux - advection terms, g{X,Y}flux - hydrostatic pressure and 
   ! p{X,Y}flux - diffusion fluxes (for eddy viscosity regularisation). 
   ! N.b. 1 - some of these fluxes will be zero in some cases, i.e. when the
   ! equation possesses no corresponding term.
   ! N.b. 2 - code particular to each equation is to be found in Equations.f90.
   subroutine CalculateFluxes(RunParams, grid, dims, &
                              tiles, tID, unitCFLTimeStep)
      implicit none

      type(RunSet), intent(in) :: RunParams
      integer, dimension(:), intent(in) :: dims
      type(GridType), intent(inout) :: grid
      type(TileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID
      real(kind=wp), intent(inout) :: unitCFLTimeStep

      integer :: i, j, k, d, nd

      integer :: nXPoints, nYPoints
      integer :: nextTile, prevTile
      real(kind=wp) :: deltaX, deltaY
      real(kind=wp) :: deltaXRecip, deltaYRecip
      real(kind=wp) :: gam_ratio

      real(kind=wp) :: wsPlus, wsMinus, aPos, aNeg, gamplus, gamneg
      real(kind=wp), dimension(:), allocatable :: uPlusConvectionFlux, uMinusConvectionFlux
      real(kind=wp), dimension(:), allocatable :: uPlusHydrostaticFlux, uMinusHydrostaticFlux
      real(kind=wp), dimension(:), allocatable :: uPlusDiffusionFlux, uMinusDiffusionFlux

      real(kind=wp) :: dif
      real(kind=wp) :: hX, hY
#if DEBUG_SPD==1 || DEBUG_SPD==2
      integer :: iu, iv

      iu = RunParams%Vars%u
      iv = RunParams%Vars%v
#endif

      nXPoints = RunParams%nXpertile
      nYPoints = RunParams%nYpertile

      nd = size(dims)

      deltaX = grid%deltaX
      deltaY = grid%deltaY
      deltaXRecip = grid%deltaXRecip
      deltaYRecip = grid%deltaYRecip

      wsPlus = 0.0_wp
      wsMinus = 0.0_wp

      allocate(uPlusConvectionFlux(nd), uMinusConvectionFlux(nd))
      allocate(uPlusHydrostaticFlux(nd), uMinusHydrostaticFlux(nd))
      allocate(uPlusDiffusionFlux(nd), uMinusDiffusionFlux(nd))

      unitCFLTimeStep = huge(unitCFLTimeStep)
      gam_ratio = 1.0_wp

      do j = 1,nYPoints
         do i = 1,nXPoints+1
            wsPlus =  xMaxWaveSpeeds(RunParams, tiles(tID)%uPlusX(:,i, j))
            wsMinus = xMaxWaveSpeeds(RunParams, tiles(tID)%uMinusX(:,i, j))
#if DEBUG_SPD == 1 || DEBUG_SPD == 2
            if (tiles(tID)%uPlusX(iU, i, j) > 50.0) then
               write(stderr, *) 'High x-velocity: uPlusX = ', tiles(tID)%uPlusX(iU,i, j), ' for tile ', tID, ' at cell index ', i, ', ', j
#if DEBUG_SPD==2
               call exit
#endif
            end if
            if (tiles(tID)%uMinusX(iU, i, j) > 50.0) then
               write(stderr, *) 'High x-velocity: uMinusX = ', tiles(tID)%uMinusX(iU,i, j), ' for tile ', tID, ' at cell index ', i, ', ', j
#if DEBUG_SPD==2
               call exit
#endif
            end if
#endif
            if (wsPlus > wsMinus) then
               aPos = wsPlus
            else
               aPos = wsMinus
            end if
            if (aPos < 0.0_wp) aPos = 0.0_wp 

            wsPlus =  xMinWaveSpeeds(RunParams, tiles(tID)%uPlusX(:,i, j))
            wsMinus = xMinWaveSpeeds(RunParams, tiles(tID)%uMinusX(:,i, j))
            if (wsPlus<wsMinus) then
               aNeg = wsPlus
            else
               aNeg = wsMinus
            end if
            if (aNeg > 0.0_wp) aNeg = 0.0_wp 

            if (aPos > epsilon(aPos)) then
               ! If gamma at cell interface > gamma at cell centre, this places
               ! a restriction on the time step in order to maintain positivity.
               if (i == 1) then
                  prevTile = tiles(tID)%West
                  gam_ratio = min(GeometricCorrectionFactor(RunParams, tiles(prevTile)%u(:, RunParams%nXpertile, j)) / &
                              GeometricCorrectionFactor(RunParams, tiles(tID)%uPlusX(:, i, j)), 1.0_wp)
               else
                  gam_ratio = min(GeometricCorrectionFactor(RunParams, tiles(tID)%u(:, i - 1, j)) / &
                              GeometricCorrectionFactor(RunParams, tiles(tID)%uPlusX(:, i, j)), 1.0_wp)
               end if
               unitCFLTimeStep = min(gam_ratio * gam_ratio * deltaX / aPos, unitCFLTimeStep)
            end if
            if (abs(aNeg) > epsilon(aNeg)) then ! changed from: if (aNeg > 1.0e-16_wp) then
               if (i == RunParams%nXpertile + 1) then
                  nextTile = tiles(tID)%East
                  gam_ratio = min(GeometricCorrectionFactor(RunParams, tiles(nextTile)%u(:, 1, j)) / &
                              GeometricCorrectionFactor(RunParams, tiles(tID)%uMinusX(:, i, j)), 1.0_wp)
               else
                  gam_ratio = min(GeometricCorrectionFactor(RunParams, tiles(tID)%u(:, i, j)) / &
                              GeometricCorrectionFactor(RunParams, tiles(tID)%uMinusX(:, i, j)), 1.0_wp)
               end if
               unitCFLTimeStep = min(gam_ratio * gam_ratio * deltaX / abs(aNeg), unitCFLTimeStep)
            end if

            call XConvectionFlux(RunParams,tiles(tID)%uPlusX(:,i,j),uPlusConvectionFlux)
            call XConvectionFlux(RunParams,tiles(tID)%uMinusX(:,i,j),uMinusConvectionFlux)
            call XHydrostaticFlux(RunParams,tiles(tID)%uPlusX(:,i,j),uPlusHydrostaticFlux)
            call XHydrostaticFlux(RunParams,tiles(tID)%uMinusX(:,i,j),uMinusHydrostaticFlux)
            if (i == RunParams%nXpertile + 1) then
               nextTile = tiles(tID)%East
               call XDiffusionFlux(RunParams,tiles(tID)%uPlusX(:,i,j),tiles(nextTile)%uLimX(:,1,j),uPlusDiffusionFlux)
               call XDiffusionFlux(RunParams,tiles(tID)%uMinusX(:,i,j),tiles(tID)%uLimX(:,i-1,j),uMinusDiffusionFlux)
            else if (i == 1) then
               prevTile = tiles(tID)%West
               call XDiffusionFlux(RunParams,tiles(tID)%uPlusX(:,i,j),tiles(tID)%uLimX(:,i,j),uPlusDiffusionFlux)
               call XDiffusionFlux(RunParams,tiles(tID)%uMinusX(:,i,j),tiles(prevTile)%uLimX(:,RunParams%nXpertile,j),uMinusDiffusionFlux)
            else
               call XDiffusionFlux(RunParams,tiles(tID)%uPlusX(:,i,j),tiles(tID)%uLimX(:,i,j),uPlusDiffusionFlux)
               call XDiffusionFlux(RunParams,tiles(tID)%uMinusX(:,i,j),tiles(tID)%uLimX(:,i-1,j),uMinusDiffusionFlux)
            end if

            dif = aPos - aNeg
            do k = 1, nd
               d = dims(k)
               if (dif < 1e-10_wp) then
                  tiles(tID)%hXFlux(d, i, j) = 0.0_wp
                  tiles(tID)%gXFlux(d, i, j) = 0.0_wp
                  tiles(tID)%pXFlux(d, i, j) = 0.0_wp
               else if (d == RunParams%Vars%w) then
                  ! This variable is time-stepped by proxy (through Hn*gam) and 
                  ! thus we need to use (Hn*gam)'s interfacial values for this flux
                  ! rather than w's. (Nb. gamplus=gamneg in our scheme!)
                  gamplus = GeometricCorrectionFactor(RunParams, tiles(tID)%uPlusX(:, i, j))
                  gamneg = GeometricCorrectionFactor(RunParams, tiles(tID)%uMinusX(:, i, j))
                  hX = tiles(tID)%uPlusX(RunParams%Vars%Hn, i, j) * gamplus - &
                       tiles(tID)%uMinusX(RunParams%Vars%Hn, i, j) * gamneg
                  hX = hX * aPos * aNeg
                  hX = hX  + (aPos * uMinusConvectionFlux(k) - aNeg * uPlusConvectionFlux(k))
                  hX = hX / dif
                  tiles(tID)%hXFlux(d, i, j) = hX
                  tiles(tID)%gXFlux(d, i, j) = 0.0_wp
                  tiles(tID)%pXFlux(d, i, j) = 0.0_wp
               else if (d == RunParams%Vars%Hnpsi) then
                  gamplus = GeometricCorrectionFactor(RunParams, tiles(tID)%uPlusX(:, i, j))
                  gamneg = GeometricCorrectionFactor(RunParams, tiles(tID)%uMinusX(:, i, j))
                  hX = tiles(tID)%uPlusX(d, i, j) * gamplus - tiles(tID)%uMinusX(d, i, j) * gamneg
                  hX = hX * aPos * aNeg
                  hX = hX  + (aPos * uMinusConvectionFlux(k) - aNeg * uPlusConvectionFlux(k))
                  hX = hX / dif
                  tiles(tID)%hXFlux(d, i, j) = hX
                  tiles(tID)%gXFlux(d, i, j) = 0.0_wp
                  tiles(tID)%pXFlux(d, i, j) = 0.0_wp
               else
                  hX = tiles(tID)%uPlusX(d, i, j) - tiles(tID)%uMinusX(d, i, j)
                  hX = hX * aPos * aNeg
                  hX = hX  + (aPos * uMinusConvectionFlux(k) - aNeg * uPlusConvectionFlux(k))
                  hX = hX / dif
                  tiles(tID)%hXFlux(d, i, j) = hX
                  tiles(tID)%gXFlux(d, i, j) =  (aPos*uMinusHydrostaticFlux(k) - aNeg*uPlusHydrostaticFlux(k))/dif
                  tiles(tID)%pXFlux(d, i, j) = 0.5_wp * (uPlusDiffusionFlux(k) + uMinusDiffusionFlux(k))
               end if
            end do
         end do
      end do

      wsPlus = 0.0_wp
      wsMinus = 0.0_wp
      aPos = 0.0_wp
      aNeg = 0.0_wp

      ! Do Y fluxes
      if (.not. RunParams%isOneD) then
         do j = 1, nYPoints + 1
            do i = 1, nXPoints
               wsPlus = yMaxWaveSpeeds(RunParams, tiles(tID)%uPlusY(:,i, j))
               wsMinus = yMaxWaveSpeeds(RunParams, tiles(tID)%uMinusY(:,i, j))
               if (wsPlus>wsMinus) then
                  aPos = wsPlus
               else
                  aPos = wsMinus
               end if
               if (aPos < 0.0_wp) aPos = 0.0_wp 

               wsPlus = yMinWaveSpeeds(RunParams, tiles(tID)%uPlusY(:,i, j))
               wsMinus = yMinWaveSpeeds(RunParams, tiles(tID)%uMinusY(:,i, j))
               if (wsPlus<wsMinus) then
                  aNeg = wsPlus
               else
                  aNeg = wsMinus
               end if
               if (aNeg > 0.0_wp) aNeg = 0.0_wp 

               if (aPos > epsilon(aPos)) then 
                  if (j == 1) then
                     prevTile = tiles(tID)%South
                     gam_ratio = min(GeometricCorrectionFactor(RunParams, tiles(prevTile)%u(:, i, RunParams%nYpertile)) / &
                                 GeometricCorrectionFactor(RunParams, tiles(tID)%uPlusY(:, i, j)), 1.0_wp)
                  else
                     gam_ratio = min(GeometricCorrectionFactor(RunParams, tiles(tID)%u(:, i, j - 1)) / &
                                 GeometricCorrectionFactor(RunParams, tiles(tID)%uPlusY(:, i, j)), 1.0_wp)
                  end if
                  unitCFLTimeStep = min(gam_ratio * gam_ratio * deltaY / aPos, unitCFLTimeStep)
               end if
               if (abs(aNeg) > epsilon(aNeg)) then 
                  if (j == RunParams%nYpertile + 1) then
                     ! look elsewhere for gam
                     nextTile = tiles(tID)%North
                     gam_ratio = min(GeometricCorrectionFactor(RunParams, tiles(nextTile)%u(:, i, 1)) / &
                                 GeometricCorrectionFactor(RunParams, tiles(tID)%uMinusY(:, i, j)), 1.0_wp)
                  else
                     gam_ratio = min(GeometricCorrectionFactor(RunParams, tiles(tID)%u(:, i, j)) / &
                                 GeometricCorrectionFactor(RunParams, tiles(tID)%uMinusY(:, i, j)), 1.0_wp)
                  end if
                  unitCFLTimeStep = min(gam_ratio * gam_ratio * deltaY / abs(aNeg), unitCFLTimeStep)
               end if

               call YConvectionFlux(RunParams,tiles(tID)%uPlusY(:,i,j),uPlusConvectionFlux)
               call YConvectionFlux(RunParams,tiles(tID)%uMinusY(:,i,j),uMinusConvectionFlux)
               call YHydrostaticFlux(RunParams,tiles(tID)%uPlusY(:,i,j),uPlusHydrostaticFlux)
               call YHydrostaticFlux(RunParams,tiles(tID)%uMinusY(:,i,j),uMinusHydrostaticFlux)
               if (j == RunParams%nYpertile + 1) then
                  nextTile = tiles(tID)%North
                  call YDiffusionFlux(RunParams,tiles(tID)%uPlusY(:,i,j),tiles(nextTile)%uLimY(:,i,1),uPlusDiffusionFlux)
                  call YDiffusionFlux(RunParams,tiles(tID)%uMinusY(:,i,j),tiles(tID)%uLimY(:,i,j-1),uMinusDiffusionFlux)
               else if (j == 1) then
                  prevTile = tiles(tID)%South
                  call YDiffusionFlux(RunParams,tiles(tID)%uPlusY(:,i,j),tiles(tID)%uLimY(:,i,j),uPlusDiffusionFlux)
                  call YDiffusionFlux(RunParams,tiles(tID)%uMinusY(:,i,j),tiles(prevTile)%uLimY(:,i,RunParams%nYpertile),uMinusDiffusionFlux)
               else
                  call YDiffusionFlux(RunParams,tiles(tID)%uPlusY(:,i,j),tiles(tID)%uLimY(:,i,j),uPlusDiffusionFlux)
                  call YDiffusionFlux(RunParams,tiles(tID)%uMinusY(:,i,j),tiles(tID)%uLimY(:,i,j-1),uMinusDiffusionFlux)
               end if

               dif = aPos - aNeg
               do k = 1, nd
                  d = dims(k)
                  if (dif < 1e-10_wp) then
                     tiles(tID)%hYFlux(d, i, j) = 0.0_wp
                     tiles(tID)%gYFlux(d, i, j) = 0.0_wp
                     tiles(tID)%pYFlux(d, i, j) = 0.0_wp
                  else if (d == RunParams%Vars%w) then
                     ! This variable is time-stepped by proxy (through Hn*gam) and 
                     ! thus we need to use (Hn*gam)'s interfacial values for this flux
                     ! rather than w's. (Nb. gamplus=gamneg in our scheme!)
                     gamplus = GeometricCorrectionFactor(RunParams, tiles(tID)%uPlusY(:, i, j))
                     gamneg = GeometricCorrectionFactor(RunParams, tiles(tID)%uMinusY(:, i, j))
                     hY = tiles(tID)%uPlusY(RunParams%Vars%Hn, i, j) * gamplus - &
                          tiles(tID)%uMinusY(RunParams%Vars%Hn, i, j) * gamneg
                     hY = hY * aPos * aNeg
                     hY = hY  + (aPos * uMinusConvectionFlux(k) - aNeg * uPlusConvectionFlux(k))
                     hY = hY / dif
                     tiles(tID)%hYFlux(d, i, j) = hY
                     tiles(tID)%gYFlux(d, i, j) = 0.0_wp
                     tiles(tID)%pYFlux(d, i, j) = 0.0_wp
                  else if (d == RunParams%Vars%Hnpsi) then
                     gamplus = GeometricCorrectionFactor(RunParams, tiles(tID)%uPlusY(:, i, j))
                     gamneg = GeometricCorrectionFactor(RunParams, tiles(tID)%uMinusY(:, i, j))
                     hY = tiles(tID)%uPlusY(d, i, j) * gamplus - tiles(tID)%uMinusY(d, i, j) * gamneg
                     hY = hY * aPos * aNeg
                     hY = hY  + (aPos * uMinusConvectionFlux(k) - aNeg * uPlusConvectionFlux(k))
                     hY = hY / dif
                     tiles(tID)%hYFlux(d, i, j) = hY
                     tiles(tID)%gYFlux(d, i, j) = 0.0_wp
                     tiles(tID)%pYFlux(d, i, j) = 0.0_wp
                  else
                     hY = tiles(tID)%uPlusY(d, i, j) - tiles(tID)%uMinusY(d, i, j)
                     hY = hY * aPos * aNeg
                     hY = hY  + (aPos * uMinusConvectionFlux(k) - aNeg * uPlusConvectionFlux(k))
                     hY = hY / dif
                     tiles(tID)%hYFlux(d, i, j) = hY
                     tiles(tID)%gYFlux(d, i, j) = (aPos*uMinusHydrostaticFlux(k) - aNeg*uPlusHydrostaticFlux(k))/dif
                     tiles(tID)%pYFlux(d, i, j) = 0.5_wp * (uPlusDiffusionFlux(k) + uMinusDiffusionFlux(k))
                  end if
               end do
            end do
         end do
      end if

   end subroutine CalculateFluxes

   ! Compute any source terms for the equations and combine these with the
   ! numerical fluxes to determine the right-hand side terms for time
   ! integration. The results are stored in 
   ! tiles(tID)%ddtExplicit(1:nFlux,:,:) and tiles(tID)%ddtImplicit(1:nFlux,:,:)
   ! which respectively record the parts on the RHS that are to be time stepped
   ! explicitly and implicity.
   subroutine ConstructHydraulicRHS(RunParams, grid, tiles, tID, t)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout) :: grid
      type(TileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID 
      real(kind=wp), intent(in) :: t

      integer :: i, j, k, d, nFlux, bt
      integer :: nXPoints, nYPoints

      real(kind=wp) :: deltaX, deltaY
      real(kind=wp) :: deltaXRecip, deltaYRecip
      real(kind=wp) :: Friction, gam, dbdx, dbdy

      real(kind=wp), dimension(:) :: STF(RunParams%nDimensions), STE(RunParams%nDimensions), STI(RunParams%nDimensions)
      real(kind=wp), dimension(:) :: gX_prefactors(RunParams%nDimensions), gY_prefactors(RunParams%nDimensions)

      nXPoints = RunParams%nXpertile
      nYPoints = RunParams%nYpertile

      nFlux = size(RunParams%iFlux)

      deltaX = grid%deltaX
      deltaY = grid%deltaY
      deltaXRecip = grid%deltaXRecip
      deltaYRecip = grid%deltaYRecip

      bt = RunParams%Vars%bt

      ! These premultiply the derivatives of the hydrostatic pressure fluxes.
      ! In particular, they are needed for the 2D surface gradient operator and
      ! are computed inside the loop below.
      gX_prefactors(:) = 0.0_wp
      gY_prefactors(:) = 0.0_wp

      if (.not.RunParams%isOneD) then
         do i = 1, nXPoints
            do j = 1, nYPoints
               gam = GeometricCorrectionFactor(RunParams, tiles(tID)%u(:, i, j))
               dbdx = tiles(tID)%u(RunParams%Vars%dbdx, i, j)
               dbdy = tiles(tID)%u(RunParams%Vars%dbdy, i, j)
               if (RunParams%geometric_factors) then
                  gX_prefactors(RunParams%Vars%rhoHnu) = (1.0_wp + dbdy * dbdy) / gam
                  gX_prefactors(RunParams%Vars%rhoHnv) = -dbdx * dbdy / gam
                  gY_prefactors(RunParams%Vars%rhoHnu) = -dbdx * dbdy / gam
                  gY_prefactors(RunParams%Vars%rhoHnv) = (1.0_wp + dbdx * dbdx) / gam
               else 
                  gX_prefactors(RunParams%Vars%rhoHnu) = 1.0_wp
                  gX_prefactors(RunParams%Vars%rhoHnv) = 0.0_wp
                  gY_prefactors(RunParams%Vars%rhoHnu) = 0.0_wp
                  gY_prefactors(RunParams%Vars%rhoHnv) = 1.0_wp
               end if
               do k = 1, nFlux
                  d = RunParams%iFlux(k)
                  if (d == RunParams%Vars%w) then
                     ! Divide by gamma to convert Hn back to hp again! (see also comment
                     ! in Equations.f90)
                     STF(d) = (tiles(tID)%hXFlux(d, i, j) - tiles(tID)%hXFlux(d, i + 1, j)) * deltaXRecip / (gam * gam) +  &
                        (tiles(tID)%hYFlux(d, i, j) - tiles(tID)%hYFlux(d, i, j + 1)) * deltaYRecip / (gam * gam)
                  else if (d == RunParams%Vars%Hnpsi) then
                     STF(d) = (tiles(tID)%hXFlux(d, i, j) - tiles(tID)%hXFlux(d, i + 1, j)) * deltaXRecip / gam +  &
                        (tiles(tID)%hYFlux(d, i, j) - tiles(tID)%hYFlux(d, i, j + 1)) * deltaYRecip / gam
                  else
                     STF(d) = KahanSum([tiles(tID)%hXFlux(d, i, j) - tiles(tID)%hXFlux(d, i + 1, j), &
                        (tiles(tID)%gXFlux(d, i, j) - tiles(tID)%gXFlux(d, i + 1, j)) * gX_prefactors(d), &
                        tiles(tID)%pXFlux(d, i + 1, j) - tiles(tID)%pXFlux(d, i, j)]) * deltaXRecip
                     STF(d) = STF(d) + KahanSum([tiles(tID)%hYFlux(d, i, j) - tiles(tID)%hYFlux(d, i, j + 1), &
                        (tiles(tID)%gYFlux(d, i, j) - tiles(tID)%gYFlux(d, i, j + 1)) * gY_prefactors(d),  &
                        tiles(tID)%pYFlux(d, i, j + 1) - tiles(tID)%pYFlux(d, i, j)]) * deltaYRecip
                  end if
               end do

               call ExplicitSourceTerms(RunParams, grid, t, tiles(tID)%x(i), tiles(tID)%y(j), &
                                        tiles(tID)%u(:,i,j), tiles(tID)%containsSource, STE)
               tiles(tID)%ddtExplicit(1:nFlux,i,j) = STF(1:nFlux) + STE(1:nFlux)

               Friction = DragClosure(RunParams, tiles(tID)%u(:,i,j))
               call ImplicitSourceTerms(RunParams, tiles(tID)%u(:,i,j), Friction, STI)
               tiles(tID)%ddtImplicit(1:nFlux,i,j) = STI(1:nFlux)
            end do
         end do
      else
         do i = 1, nXPoints
            gam = GeometricCorrectionFactor(RunParams, tiles(tID)%u(:, i, 1))
            do k = 1, nFlux
               d = RunParams%iFlux(k)
               if (d == RunParams%Vars%w) then
                  ! divide by gamma to convert H back to h again! (see also comment
                  ! in Equations.f90)
                  STF(d) = (tiles(tID)%hXFlux(d, i, 1) - tiles(tID)%hXFlux(d, i + 1, 1)) * deltaXRecip / (gam * gam)
               else if (d == RunParams%Vars%Hnpsi) then
                  STF(d) = (tiles(tID)%hXFlux(d, i, 1) - tiles(tID)%hXFlux(d, i + 1, 1)) * deltaXRecip / gam
               else
                  STF(d) = (tiles(tID)%hXFlux(d, i, 1) - tiles(tID)%hXFlux(d, i + 1, 1) +  &
                     (tiles(tID)%gXFlux(d, i, 1) - tiles(tID)%gXFlux(d, i + 1, 1)) / gam +  &
                     tiles(tID)%pXFlux(d, i + 1, 1) - tiles(tID)%pXFlux(d, i, 1)) * deltaXRecip
                  STF(d) = KahanSum([tiles(tID)%hXFlux(d, i, 1) - tiles(tID)%hXFlux(d, i + 1, 1), &
                     (tiles(tID)%gXFlux(d, i, 1) - tiles(tID)%gXFlux(d, i + 1, 1)) / gam,  &
                     tiles(tID)%pXFlux(d, i + 1, 1) - tiles(tID)%pXFlux(d, i, 1)]) * deltaXRecip
               end if
            end do

            call ExplicitSourceTerms(RunParams, grid, t, tiles(tID)%x(i), tiles(tID)%y(1), &
                                     tiles(tID)%u(:,i,1), tiles(tID)%containsSource, STE)
            tiles(tID)%ddtExplicit(1:nFlux,i,1) = STF(1:nFlux) + STE(1:nFlux)

            Friction = DragClosure(RunParams, tiles(tID)%u(:,i,1))
            call ImplicitSourceTerms(RunParams, tiles(tID)%u(:,i,1), Friction, STI)
            tiles(tID)%ddtImplicit(1:nFlux,i,1) = STI(1:nFlux)
         end do
      end if

   end subroutine ConstructHydraulicRHS

end module hydraulic_rhs_module
