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


! This module contains routines to calculate the morphodynamic part of the RHS
! of the governing equations. In particular, it computes the source term for the
! bed evolution equation only.  In our framework, morphodynamic updates for
! other variables depend linearly on the bed elevation update, so this is all
! that's needed for the time stepper (see also
! https://arxiv.org/abs/2306.16185).
!
! That is, we compute the RHS of
!
! db/dt = -gamma * M(u),
!
! where gamma is a geometric correction (if used), u represents the full vector
! of primary flow variables and 
!
! M = (E - D) / (1 - RunParams%BedPorosity), 
!
! with E and D being user-specified closures for flow-dependent erosion and
! deposition rates. Additionally, we damp morphodynamics to zero in the limit of
! low flow depths in order to avoid unphysical bulking from fast thin layers.
!
! The calculation is slightly complicated by the fact that b is stored on cell
! vertices, while the other flow variables are defined and time stepped at cell
! centres. So compute M(u) and interpolate it for the benefit of b.
! Consequently, this module contains other topographic interpolation routines
! that are needed elsewhere in the code and lack a natural home.
module morphodynamic_rhs_module

   use grid_module, only: GridType, TileType, TileList, IsActiveGhostTile
   use set_precision_module, only: wp
   use utilities_module, only: KahanSum, InVector
   use runsettings_module, only: RunSet
   use equations_module, only: ErosionDepositionTerms
   use closures_module, only: GeometricCorrectionFactor, GeometricCorrectionFactor_gradin_scalar, GeometricCorrectionFactor_gradin_array, DragClosure
   use limiters_module, only: limiter
   use blur_module, only: MultiBinom3Blur

   implicit none

   private
   public :: CalculateMorphodynamicRHS
   public :: ComputeCellCentredTopographicData
   public :: ComputeTopographicCurvatures
   public :: ComputeMorphodynamicCurvatures
   public :: ComputeInterfacialTopographicData
   public :: EqualiseTopographicBoundaryData

   interface ComputeCellCentredTopographicData
      module procedure :: ComputeCellCentredTopographicData_tileID, ComputeCellCentredTopographicData_ij
   end interface

contains

   ! Calculate -gamma * M(u) for every active tile in tiles. The result is
   ! stored in ddtExplicitBt (by BtSourceTerm), for use by the time stepper.
   subroutine CalculateMorphodynamicRHS(RunParams, grid, tiles)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      type(tileType), dimension(:), intent(inout) :: tiles

      type(TileList), pointer :: activeTiles

      real(kind=wp) :: Friction, Ero, Depo 
      real(kind=wp) :: HnW, HnE, HnS, HnN, Hneps
      real(kind=wp) :: gam

      integer :: tt, tID, ttW, ttE, ttS, ttN
      integer :: i, j, nd, bt, iHn

      activeTiles => grid%activeTiles
      nd = size(RunParams%iFlux)
      iHn = RunParams%Vars%Hn
      bt = RunParams%Vars%bt

      Hneps = RunParams%heightThreshold
      HnS = 0.0_wp
      HnN = 0.0_wp

      ! do tt = 1, ActiveTiles%Size
!$omp parallel do schedule(auto), default(none), &
!$omp private(tt, tID, i, j, ttW, HnW, ttE, HnE, ttS, HnS, ttN, HnN, gam, Friction, Ero, Depo), &
!$omp shared(ActiveTiles, RunParams, tiles, grid, iHn, Hneps, nd, GeometricCorrectionFactor, DragClosure)
      do tt = 1, ActiveTiles%Size
         tID = ActiveTiles%List(tt)
         do i = 1, RunParams%nXPertile
            do j = 1, RunParams%nYPertile
               if (i == 1) then
                  ttW = tiles(tID)%West
                  HnW = tiles(ttW)%u(iHn, RunParams%nXpertile, j)
               else
                  HnW = tiles(tID)%u(iHn, i - 1, j)
               end if
               if (i == RunParams%nXpertile) then
                  ttE = tiles(tID)%East
                  HnE = tiles(ttE)%u(iHn, 1, j)
               else
                  HnE = tiles(tID)%u(iHn, i + 1, j)
               end if

               if (.not. RunParams%isOneD) then
                  if (j == 1) then
                     ttS = tiles(tID)%South
                     HnS = tiles(ttS)%u(iHn, i, RunParams%nYpertile)
                  else
                     HnS = tiles(tID)%u(iHn, i, j - 1)
                  end if
                  if (j == RunParams%nYpertile) then
                     ttN = tiles(tID)%North
                     HnN = tiles(ttN)%u(iHn, i, 1)
                  else
                     HnN = tiles(tID)%u(iHn, i, j + 1)
                  end if
               end if

               ! If point is next to a dry cell (or is itself dry), no
               ! morphodynamics. This is meant to defend against a numerical
               ! ratcheting phenomenon that can occur at flow fronts in some
               ! scenarios, where erosion at the front leads to bulking and in
               ! turn, erosion at a neighbouring (previously dry) cell.
               if (tiles(tID)%u(iHn,i,j) < Hneps .or. HnW < Hneps .or. HnE < Hneps &
                  .or. ((.not. RunParams%isOneD) .and. (HnS < Hneps .or. HnN < Hneps))) then
                  tiles(tID)%ddtExplicit(1:nd,i,j) = 0.0_wp
                  tiles(tID)%EminusD(i,j) = 0.0_wp
               else
                  ! Compute E - D at each cell centre and save it.
                  gam = GeometricCorrectionFactor(RunParams, tiles(tID)%u(:,i,j))
                  Friction = DragClosure(RunParams, tiles(tID)%u(:,i,j), gam)
                  call ErosionDepositionTerms(RunParams, tiles(tID)%u(:,i,j), gam, Ero, Depo)
                  tiles(tID)%EminusD(i,j) = Ero - Depo
               end if
            end do
         end do
      end do
!$omp end parallel do

      call BtSourceTerm(RunParams, grid, tiles)

   end subroutine

   ! Take E - D on the cell centres (which we stored above) and interpolate to
   ! obtain the source term for bt which needs to be time stepped using the cell
   ! vertices as grid points.
   subroutine BtSourceTerm(RunParams, grid, tiles)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      type(tileType), dimension(:), intent(inout) :: tiles

      type(TileList), pointer :: activeTiles
      integer :: tt, tID, ttE, ttW, ttN, ttS, ttSW, ttSE, ttNW, ttNE
      integer :: i, j, nXpertile, nYpertile, idbdx, idbdy
      real(kind=wp) :: psib, dbdx, dbdy, gam

      activeTiles => grid%activeTiles

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      idbdx = RunParams%Vars%dbdx
      idbdy = RunParams%Vars%dbdy

      psib = 1.0_wp - RunParams%BedPorosity

   if (.not. RunParams%isOneD) then
!$omp parallel do schedule(auto), default(none), &
!$omp private(tt, tID, i, j, dbdx, dbdy, gam, ttW, ttE, ttS, ttN, ttSW, ttSE, ttNW, ttNE), &
!$omp shared(ActiveTiles, tiles, idbdx, idbdy, nXpertile, nYpertile, psib, GeometricCorrectionFactor_gradin_scalar)
      do tt = 1, ActiveTiles%Size
         tID = ActiveTiles%List(tt)
         
            ! interior
            do i = 2, nXpertile
               do j = 2, nYpertile
                  dbdx = 0.25_wp * KahanSum([tiles(tID)%u(idbdx, i - 1, j - 1), tiles(tID)%u(idbdx, i - 1, j), &
                     tiles(tID)%u(idbdx, i    , j - 1), tiles(tID)%u(idbdx, i    , j)])
                  dbdy = 0.25_wp * KahanSum([tiles(tID)%u(idbdy, i - 1, j - 1), tiles(tID)%u(idbdy, i - 1, j), &
                     tiles(tID)%u(idbdy, i    , j - 1), tiles(tID)%u(idbdy, i    , j)])
                  gam = GeometricCorrectionFactor_gradin_scalar(dbdx, dbdy)
                  tiles(tID)%ddtExplicitBt(i,j) = -0.25_wp * gam / psib *  &
                     KahanSum([tiles(tID)%EminusD(i - 1, j - 1), tiles(tID)%EminusD(i - 1, j), tiles(tID)%EminusD(i, j - 1), tiles(tID)%EminusD(i, j)])
               end do
            end do

            ! computations on edges/corners
            ttW = tiles(tID)%West
            ttE = tiles(tID)%East
            ttS = tiles(tID)%South
            ttN = tiles(tID)%North
            ttSW = tiles(tID)%SouthWest
            ttSE = tiles(tID)%SouthEast
            ttNW = tiles(tID)%NorthWest
            ttNE = tiles(tID)%NorthEast

            do j = 2, nYpertile
               ! west bdry
               dbdx = 0.25_wp * KahanSum([tiles(ttW)%u(idbdx, nXpertile, j - 1), tiles(ttW)%u(idbdx, nXpertile, j), &
                  tiles(tID)%u(idbdx, 1, j - 1), tiles(tID)%u(idbdx, 1, j)])
               dbdy = 0.25_wp * KahanSum([tiles(ttW)%u(idbdy, nXpertile, j - 1), tiles(ttW)%u(idbdy, nXpertile, j), &
                  tiles(tID)%u(idbdy, 1, j - 1), tiles(tID)%u(idbdy, 1, j)])
               gam = GeometricCorrectionFactor_gradin_scalar(dbdx, dbdy)
               tiles(tID)%ddtExplicitBt(1, j) = -0.25_wp * gam / psib *  &
                  KahanSum([tiles(ttW)%EminusD(nXpertile, j - 1), tiles(ttW)%EminusD(nXpertile, j), tiles(tID)%EminusD(1, j - 1) + tiles(tID)%EminusD(1, j)]) 
               ! east bdry
               dbdx = 0.25_wp * KahanSum([tiles(tID)%u(idbdx, nXpertile, j - 1), tiles(tID)%u(idbdx, nXpertile, j), &
                  tiles(ttE)%u(idbdx, 1, j - 1), tiles(ttE)%u(idbdx, 1, j)])
               dbdy = 0.25_wp * KahanSum([tiles(tID)%u(idbdy, nXpertile, j - 1), tiles(tID)%u(idbdy, nXpertile, j), &
                  tiles(ttE)%u(idbdy, 1, j - 1), tiles(ttE)%u(idbdy, 1, j)])
               gam = GeometricCorrectionFactor_gradin_scalar(dbdx, dbdy)
               tiles(tID)%ddtExplicitBt(nXpertile + 1, j) = -0.25_wp * gam / psib *  &
                  KahanSum([tiles(tID)%EminusD(nXpertile, j - 1), tiles(tID)%EminusD(nXpertile, j), tiles(ttE)%EminusD(1, j - 1), tiles(ttE)%EminusD(1, j)])
            end do
            do i = 2, nXpertile
               ! north bdry
               dbdx = 0.25_wp * KahanSum([tiles(tID)%u(idbdx, i - 1, nYpertile), tiles(ttN)%u(idbdx, i - 1, 1), &
                  tiles(tID)%u(idbdx, i, nYpertile), tiles(ttN)%u(idbdx, i, 1)])
               dbdy = 0.25_wp * KahanSum([tiles(tID)%u(idbdy, i - 1, nYpertile), tiles(ttN)%u(idbdy, i - 1, 1), &
                  tiles(tID)%u(idbdy, i, nYpertile), tiles(ttN)%u(idbdy, i, 1)])
               gam = GeometricCorrectionFactor_gradin_scalar(dbdx, dbdy)
               tiles(tID)%ddtExplicitBt(i, nYpertile + 1) = -0.25_wp * gam / psib *  &
                  KahanSum([tiles(tID)%EminusD(i - 1, nYpertile), tiles(ttN)%EminusD(i - 1, 1), tiles(tID)%EminusD(i, nYpertile), tiles(ttN)%EminusD(i, 1)])
               ! south bdry
               dbdx = 0.25_wp * KahanSum([tiles(ttS)%u(idbdx, i - 1, nYpertile), tiles(tID)%u(idbdx, i - 1, 1), &
                  tiles(ttS)%u(idbdx, i, nYpertile), tiles(tID)%u(idbdx, i, 1)])
               dbdy = 0.25_wp * KahanSum([tiles(ttS)%u(idbdy, i - 1, nYpertile), tiles(tID)%u(idbdy, i - 1, 1), &
                  tiles(ttS)%u(idbdy, i, nYpertile), tiles(tID)%u(idbdy, i, 1)])
               gam = GeometricCorrectionFactor_gradin_scalar(dbdx, dbdy)
               tiles(tID)%ddtExplicitBt(i, 1) = -0.25_wp * gam / psib * &
                  KahanSum([tiles(ttS)%EminusD(i - 1, nYpertile), tiles(tID)%EminusD(i - 1, 1), tiles(ttS)%EminusD(i, nYpertile), tiles(tID)%EminusD(i, 1)])
            end do
            ! corners
            ! SW
            dbdx = 0.25_wp * KahanSum([tiles(ttSW)%u(idbdx, nXpertile, nYpertile), tiles(ttW)%u(idbdx, nXpertile, 1), &
               tiles(ttS )%u(idbdx, 1, nYpertile), tiles(tID)%u(idbdx, 1, 1)])
            dbdy = 0.25_wp * KahanSum([tiles(ttSW)%u(idbdy, nXpertile, nYpertile), tiles(ttW)%u(idbdy, nXpertile, 1), &
               tiles(ttS )%u(idbdy, 1, nYpertile), tiles(tID)%u(idbdy, 1, 1)])
            gam = GeometricCorrectionFactor_gradin_scalar(dbdx, dbdy)
            tiles(tID)%ddtExplicitBt(1, 1) = -0.25_wp * gam / psib *  &
               KahanSum([tiles(ttSW)%EminusD(nXpertile, nYpertile), tiles(ttW)%EminusD(nXpertile, 1), tiles(ttS )%EminusD(1, nYpertile), tiles(tID)%EminusD(1, 1)])
            ! SE
            dbdx = 0.25_wp * KahanSum([tiles(ttS )%u(idbdx, nXpertile, nYpertile), tiles(tID)%u(idbdx, nXpertile, 1), &
               tiles(ttSE)%u(idbdx, 1, nYpertile), tiles(ttE)%u(idbdx, 1, 1)])
            dbdy = 0.25_wp * KahanSum([tiles(ttS )%u(idbdy, nXpertile, nYpertile), tiles(tID)%u(idbdy, nXpertile, 1), &
               tiles(ttSE)%u(idbdy, 1, nYpertile), tiles(ttE)%u(idbdy, 1, 1)])
            gam = GeometricCorrectionFactor_gradin_scalar(dbdx, dbdy)
            tiles(tID)%ddtExplicitBt(nXpertile + 1, 1) = -0.25_wp * gam / psib *  &
               KahanSum([tiles(ttS )%EminusD(nXpertile, nYpertile), tiles(tID)%EminusD(nXpertile, 1), tiles(ttSE)%EminusD(1, nYpertile), tiles(ttE)%EminusD(1, 1)])
            ! NW
            dbdx = 0.25_wp * KahanSum([tiles(ttW)%u(idbdx, nXpertile, nYpertile), tiles(ttNW)%u(idbdx, nXpertile, 1), &
               tiles(tID)%u(idbdx, 1, nYpertile), tiles(ttN )%u(idbdx, 1, 1)])
            dbdy = 0.25_wp * KahanSum([tiles(ttW)%u(idbdy, nXpertile, nYpertile), tiles(ttNW)%u(idbdy, nXpertile, 1), &
               tiles(tID)%u(idbdy, 1, nYpertile), tiles(ttN )%u(idbdy, 1, 1)])
            gam = GeometricCorrectionFactor_gradin_scalar(dbdx, dbdy)
            tiles(tID)%ddtExplicitBt(1, nYpertile + 1) = -0.25_wp * gam / psib *  &
               KahanSum([tiles(ttW)%EminusD(nXpertile, nYpertile), tiles(ttNW)%EminusD(nXpertile, 1), &
               tiles(tID)%EminusD(1, nYpertile), tiles(ttN )%EminusD(1, 1)])
            ! NE
            dbdx = 0.25_wp * KahanSum([tiles(tID)%u(idbdx, nXpertile, nYpertile), tiles(ttN )%u(idbdx, nXpertile, 1), &
               tiles(ttE)%u(idbdx, 1, nYpertile) + tiles(ttNE)%u(idbdx, 1, 1)])
            dbdy = 0.25_wp * KahanSum([tiles(tID)%u(idbdy, nXpertile, nYpertile), tiles(ttN )%u(idbdy, nXpertile, 1), &
               tiles(ttE)%u(idbdy, 1, nYpertile), tiles(ttNE)%u(idbdy, 1, 1)])
            gam = GeometricCorrectionFactor_gradin_scalar(dbdx, dbdy)
            tiles(tID)%ddtExplicitBt(nXpertile + 1, nYpertile + 1) = -0.25_wp * gam / psib *  &
               KahanSum([tiles(tID)%EminusD(nXpertile, nYpertile), tiles(ttN )%EminusD(nXpertile, 1), &
               tiles(ttE)%EminusD(1, nYpertile), tiles(ttNE)%EminusD(1, 1)])
         end do
!$omp end parallel do
      else
!$omp parallel do schedule(auto), default(none), &
!$omp private(tt, tID, dbdx, gam, ttW, ttE), &
!$omp shared(ActiveTiles, tiles, grid, idbdx, idbdy, nXpertile, nYpertile, psib, GeometricCorrectionFactor_gradin_scalar)
         do tt = 1, ActiveTiles%Size
            tID = ActiveTiles%List(tt)
            ! 1D is much simpler...
            do i = 2, nXpertile
               dbdx = 0.5_wp * (tiles(tID)%u(idbdx, i - 1, 1) + tiles(tID)%u(idbdx, i, 1))
               gam = GeometricCorrectionFactor_gradin_scalar(dbdx, 0.0_wp)
               tiles(tID)%ddtExplicitBt(i,1) = -0.5_wp * gam * (tiles(tID)%EminusD(i - 1, 1) + tiles(tID)%EminusD(i, 1)) / psib
            end do

            ! handle tile boundaries
            ttW = tiles(tID)%West
            if (ttW > 0 .and. grid%tileContainer(ttW)%TileOn) then
               dbdx = 0.5_wp * (tiles(ttW)%u(iDBDX, nXpertile, 1) + tiles(tID)%u(iDBDX, 1, 1))
               gam = GeometricCorrectionFactor_gradin_scalar(dbdx, 0.0_wp)
               tiles(tID)%ddtExplicitBt(1,1) = -0.5_wp * gam * (tiles(ttW)%EminusD(nXpertile, 1) + tiles(tID)%EminusD(1, 1)) / psib
            else
               dbdx = 0.5_wp * tiles(tID)%u(idbdx, 1, 1)
               gam = GeometricCorrectionFactor_gradin_scalar(dbdx, 0.0_wp)
               tiles(tID)%ddtExplicitBt(1,1) = -0.5_wp * gam * tiles(tID)%EminusD(1, 1) / psib
            end if

            ttE = tiles(tID)%East
            if (ttE > 0 .and. grid%tileContainer(ttE)%TileOn) then
               dbdx = 0.5_wp * (tiles(tID)%u(idbdx, nXpertile, 1) + tiles(ttE)%u(idbdx, 1, 1))
               gam = GeometricCorrectionFactor_gradin_scalar(dbdx, 0.0_wp)
               tiles(tID)%ddtExplicitBt(nXpertile+1,1) = -0.5_wp * gam * (tiles(tID)%EminusD(nXpertile, 1) + tiles(ttE)%EminusD(1, 1)) / psib
            else
               dbdx = 0.5_wp * tiles(tID)%u(idbdx, nXpertile, 1)
               gam = GeometricCorrectionFactor_gradin_scalar(dbdx, 0.0_wp)
               tiles(tID)%ddtExplicitBt(nXpertile+1,1) = -0.5_wp * gam * tiles(tID)%EminusD(nXpertile, 1) / psib
            end if
         end do
!$omp end parallel do
      end if

   end subroutine

   ! Interpolate b at cell centres, for tile tID.
   pure subroutine ComputeCellCentredTopographicData_tileID(RunParams, grid, tiles, tID)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout), target :: grid
      type(tileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID

      integer :: i, j, ib0, ibt, idbdx, idbdy
      integer :: nXpertile, nYpertile

      real(kind=wp) :: deltaXRecip, deltaYRecip
      real(kind=wp) :: dbdx, dbdy, b0_centre, bt_centre

      ib0 = RunParams%Vars%b0
      ibt = RunParams%Vars%bt
      idbdx = RunParams%Vars%dbdx
      idbdy = RunParams%Vars%dbdy
      
      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      deltaXRecip = grid%deltaXRecip
      deltaYRecip = grid%deltaYRecip

      if (.not. RunParams%isOneD) then
         do i = 1, nXpertile
            do j = 1, nYpertile

               b0_centre = 0.25_wp * KahanSum([tiles(tID)%b0(i, j), tiles(tID)%b0(i+1,j), &
                  tiles(tID)%b0(i, j+1), tiles(tID)%b0(i+1, j+1)])
               bt_centre = 0.25_wp * KahanSum([tiles(tID)%bt(i, j), tiles(tID)%bt(i+1,j), &
                  tiles(tID)%bt(i, j+1), tiles(tID)%bt(i+1, j+1)])
               tiles(tID)%u(ib0,i,j) = b0_centre
               tiles(tID)%u(ibt,i,j) = bt_centre

               dbdx = 0.5_wp * deltaXRecip *  &
                  KahanSum([tiles(tID)%b0(i+1, j  ), tiles(tID)%bt(i+1, j  ), &
                  -tiles(tID)%b0(i  , j  ), -tiles(tID)%bt(i  , j  ), &
                  tiles(tID)%b0(i+1, j+1), tiles(tID)%bt(i+1, j+1), &
                  -tiles(tID)%b0(i  , j+1), -tiles(tID)%bt(i  , j+1)])
               tiles(tID)%u(idbdx,i,j) = dbdx

               dbdy = 0.5_wp * deltaYRecip *  &
                  KahanSum([tiles(tID)%b0(i  , j+1), tiles(tID)%bt(i  , j+1), &
                  -tiles(tID)%b0(i  , j  ), -tiles(tID)%bt(i  , j  ), &
                  tiles(tID)%b0(i+1, j+1), tiles(tID)%bt(i+1, j+1), &
                  -tiles(tID)%b0(i+1, j  ), -tiles(tID)%bt(i+1, j  )])
               tiles(tID)%u(idbdy,i,j) = dbdy
            end do
         end do

      else
         do i = 1, nXpertile
            b0_centre = 0.5_wp * (tiles(tID)%b0(i, 1) + tiles(tID)%b0(i+1, 1))
            bt_centre = 0.5_wp * (tiles(tID)%bt(i, 1) + tiles(tID)%bt(i+1, 1))
            tiles(tID)%u(ib0,i,1) = b0_centre
            tiles(tID)%u(ibt,i,1) = bt_centre
            dbdx = deltaXRecip * KahanSum([tiles(tID)%b0(i+1, 1), tiles(tID)%bt(i+1, 1), &
               -tiles(tID)%b0(i  , 1), -tiles(tID)%bt(i  , 1)])
            tiles(tID)%u(idbdx,i,1) = dbdx
         end do
      end if

   end subroutine ComputeCellCentredTopographicData_tileID

   ! Interpolate b at cell centres, for cell i,j in tile tID.
   pure subroutine ComputeCellCentredTopographicData_ij(RunParams, grid, tiles, tID, i, j)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout), target :: grid
      type(tileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID, i, j

      integer :: ibt, idbdx, idbdy

      real(kind=wp) :: deltaXRecip, deltaYRecip
      real(kind=wp) :: dbdx, dbdy, bt_centre

      ibt = RunParams%Vars%bt
      idbdx = RunParams%Vars%dbdx
      idbdy = RunParams%Vars%dbdy

      deltaXRecip = grid%deltaXRecip
      deltaYRecip = grid%deltaYRecip

      if (.not. RunParams%isOneD) then
         bt_centre = 0.25_wp * KahanSum([tiles(tID)%bt(i, j  ), tiles(tID)%bt(i+1,   j), &
            tiles(tID)%bt(i, j+1), tiles(tID)%bt(i+1, j+1)])
         tiles(tID)%u(ibt,i,j) = bt_centre

         dbdx = 0.5_wp * deltaXRecip *  &
            KahanSum([tiles(tID)%b0(i+1, j  ), tiles(tID)%bt(i+1, j  ), &
            -tiles(tID)%b0(i  , j  ), -tiles(tID)%bt(i  , j  ), &
            tiles(tID)%b0(i+1, j+1), tiles(tID)%bt(i+1, j+1), &
            -tiles(tID)%b0(i  , j+1), -tiles(tID)%bt(i  , j+1)])
         tiles(tID)%u(idbdx,i,j) = dbdx

         dbdy = 0.5_wp * deltaYRecip *  &
            KahanSum([tiles(tID)%b0(i  , j+1), tiles(tID)%bt(i  , j+1), &
            -tiles(tID)%b0(i  , j  ), -tiles(tID)%bt(i  , j  ), &
            tiles(tID)%b0(i+1, j+1), tiles(tID)%bt(i+1, j+1), &
            -tiles(tID)%b0(i+1, j  ), -tiles(tID)%bt(i+1, j  )])
         tiles(tID)%u(idbdy,i,j) = dbdy
      else
         bt_centre = 0.5_wp * (tiles(tID)%bt(i, 1) + tiles(tID)%bt(i+1, 1))
         tiles(tID)%u(ibt,i,1) = bt_centre
         dbdx = deltaXRecip * KahanSum([tiles(tID)%b0(i+1, 1), tiles(tID)%bt(i+1, 1), &
            -tiles(tID)%b0(i  , 1), -tiles(tID)%bt(i  , 1)])
         tiles(tID)%u(idbdx,i,1) = dbdx
      end if
   end subroutine ComputeCellCentredTopographicData_ij

   ! Interpolate b at cell interfaces, for tile tID.
   pure subroutine ComputeInterfacialTopographicData(RunParams, grid, tiles, tID)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout), target :: grid
      type(tileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID

      integer :: i, j, prevTile, nextTile, southTile, northTile
      integer :: ibt, idbdx, idbdy, nXpertile, nYpertile
      real(kind=wp) :: bt, dbdx, dbdy, deltaXRecip, deltaYRecip

      ibt = RunParams%Vars%bt
      idbdx = RunParams%Vars%dbdx
      idbdy = RunParams%Vars%dbdy

      deltaXRecip = grid%deltaXRecip
      deltaYRecip = grid%deltaYRecip
      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      if (.not. RunParams%isOneD) then
         ! u{Plus,Minus}X
         ! x derivs
         do i = 2, nXpertile
            do j = 1, nYpertile
               dbdx = 0.25_wp * deltaXRecip * &
                  KahanSum([tiles(tID)%b0(i+1,j  ), tiles(tID)%bt(i+1,j  ), &
                     tiles(tID)%b0(i+1,j+1), tiles(tID)%bt(i+1,j+1), &
                     -tiles(tID)%b0(i-1,j  ), -tiles(tID)%bt(i-1,j  ), &
                     -tiles(tID)%b0(i-1,j+1), - tiles(tID)%bt(i-1,j+1)])
               tiles(tID)%uMinusX(idbdx,i,j) = dbdx
               tiles(tID)%uPlusX(idbdx,i,j) = dbdx
            end do
         end do

         ! lh boundary
         prevTile = tiles(tID)%West
         do j = 1, nYpertile
            dbdx = 0.25_wp * deltaXRecip * &
               KahanSum([tiles(tID)%b0(2,j  ), tiles(tID)%bt(2,j  ), &
                  tiles(tID)%b0(2,j+1), tiles(tID)%bt(2,j+1), &
                  -tiles(prevTile)%b0(nXpertile,j  ), -tiles(prevTile)%bt(nXpertile,j  ), &
                  -tiles(prevTile)%b0(nXpertile,j+1), -tiles(prevTile)%bt(nXpertile,j+1)])
            tiles(tID)%uMinusX(idbdx,1,j) = dbdx
            tiles(tID)%uPlusX(idbdx,1,j) = dbdx
         end do

         ! rh boundary
         nextTile = tiles(tID)%East
         do j = 1, nYpertile
            dbdx = 0.25_wp * deltaXRecip * &
               KahanSum([tiles(nextTile)%b0(2,j  ), tiles(nextTile)%bt(2,j  ), &
                  tiles(nextTile)%b0(2,j+1), tiles(nextTile)%bt(2,j+1), &
                  -tiles(tID)%b0(nXpertile,j  ), -tiles(tID)%bt(nXpertile,j  ), &
                  -tiles(tID)%b0(nXpertile,j+1), -tiles(tID)%bt(nXpertile,j+1)])
            tiles(tID)%uMinusX(idbdx,nXpertile+1,j) = dbdx
            tiles(tID)%uPlusX(idbdx,nXpertile+1,j) = dbdx
         end do

         ! y derivs & bt
         do i = 1, nXpertile + 1
            do j = 1, nYpertile
               bt = 0.5_wp * (tiles(tID)%bt(i,j) + tiles(tID)%bt(i,j+1))
               dbdy = deltaYRecip * &
                  KahanSum([tiles(tID)%b0(i,j+1), tiles(tID)%bt(i,j+1), &
                     -tiles(tID)%b0(i,j  ), -tiles(tID)%bt(i,j  )])
               tiles(tID)%uPlusX(ibt,i,j) = bt
               tiles(tID)%uMinusX(ibt,i,j) = bt
               tiles(tID)%uPlusX(idbdy,i,j) = dbdy
               tiles(tID)%uMinusX(idbdy,i,j) = dbdy
            end do
         end do

         ! u{Plus,Minus}Y
         ! x derivs & bt
         do i = 1, nXpertile
            do j = 1, nYpertile + 1
               bt = 0.5_wp * (tiles(tID)%bt(i,j) + tiles(tID)%bt(i+1,j))
               dbdx = deltaXRecip * &
                  KahanSum([tiles(tID)%b0(i+1,j), tiles(tID)%bt(i+1,j), &
                     -tiles(tID)%b0(i  ,j), -tiles(tID)%bt(i  ,j)])
               tiles(tID)%uPlusY(ibt,i,j) = bt
               tiles(tID)%uMinusY(ibt,i,j) = bt
               tiles(tID)%uPlusY(idbdx,i,j) = dbdx
               tiles(tID)%uMinusY(idbdx,i,j) = dbdx
            end do
         end do

         ! y derivs
         do i = 1, nXpertile
            do j = 2, nYpertile
               dbdy = 0.25_wp * deltaYRecip * &
                  KahanSum([tiles(tID)%b0(i+1,j+1), tiles(tID)%bt(i+1,j+1), &
                     tiles(tID)%b0(i  ,j+1), tiles(tID)%bt(i  ,j+1), &
                     -tiles(tID)%b0(i+1,j-1), -tiles(tID)%bt(i+1,j-1), &
                     -tiles(tID)%b0(i  ,j-1), -tiles(tID)%bt(i  ,j-1)])
               tiles(tID)%uPlusY(idbdy,i,j) = dbdy
               tiles(tID)%uMinusY(idbdy,i,j) = dbdy
            end do
         end do

         ! bottom bdry
         southTile = tiles(tID)%South
         do i = 1, nXpertile
            dbdy = 0.25_wp * deltaYRecip * &
               KahanSum([tiles(tID)%b0(i+1,2), tiles(tID)%bt(i+1,2), &
                  tiles(tID)%b0(i  ,2), tiles(tID)%bt(i  ,2), &
                  -tiles(southTile)%b0(i+1,nYpertile), -tiles(southTile)%bt(i+1,nYpertile), &
                  -tiles(southTile)%b0(i  ,nYpertile), - tiles(southTile)%bt(i  ,nYpertile)])
            tiles(tID)%uPlusY(idbdy,i,1) = dbdy
            tiles(tID)%uMinusY(idbdy,i,1) = dbdy
         end do

         ! top bdry
         northTile = tiles(tID)%North
         do i = 1, nXpertile
            dbdy = 0.25_wp * deltaYRecip * &
               KahanSum([tiles(northTile)%b0(i+1,2), tiles(northTile)%bt(i+1,2), &
                  tiles(northTile)%b0(i  ,2), tiles(northTile)%bt(i  ,2), &
                  -tiles(tID)%b0(i+1,nYpertile), -tiles(tID)%bt(i+1,nYpertile), &
                  -tiles(tID)%b0(i  ,nYpertile), -tiles(tID)%bt(i  ,nYpertile)])
            tiles(tID)%uPlusY(idbdy,i,nYpertile+1) = dbdy
            tiles(tID)%uMinusY(idbdy,i,nYpertile+1) = dbdy
         end do
      else
         do i = 2, nXpertile
            bt = tiles(tID)%bt(i,1)
            dbdx = 0.5_wp * deltaXRecip * &
               KahanSum([tiles(tID)%b0(i+1,1), tiles(tID)%bt(i+1,1), &
                  -tiles(tID)%b0(i-1,1), -tiles(tID)%bt(i-1,1)])
            tiles(tID)%uPlusX(ibt,i,1) = bt
            tiles(tID)%uMinusX(ibt,i,1) = bt
            tiles(tID)%uPlusX(idbdx,i,1) = dbdx
            tiles(tID)%uMinusX(idbdx,i,1) = dbdx
         end do

         ! lh boundary
         prevTile = tiles(tID)%West
         bt = tiles(tID)%bt(1,1)
         dbdx = 0.5_wp * deltaXRecip * &
            KahanSum([tiles(tID)%b0(2,1), tiles(tID)%bt(2,1), &
               -tiles(prevTile)%b0(nXpertile,1), -tiles(prevTile)%bt(nXpertile,1)])
         tiles(tID)%uPlusX(ibt,1,1) = bt
         tiles(tID)%uMinusX(ibt,1,1) = bt
         tiles(tID)%uPlusX(idbdx,1,1) = dbdx
         tiles(tID)%uMinusX(idbdx,1,1) = dbdx

         ! rh boundary
         nextTile = tiles(tID)%East
         bt = tiles(tID)%bt(nXpertile+1,1)
         dbdx = 0.5_wp * deltaXRecip * &
            KahanSum([tiles(nextTile)%b0(2,1), tiles(nextTile)%bt(2,1), &
               -tiles(tID)%b0(nXpertile,1), -tiles(tID)%bt(nXpertile,1)])
         tiles(tID)%uPlusX(ibt,RunParams%nXpertile+1,1) = bt
         tiles(tID)%uMinusX(ibt,RunParams%nXpertile+1,1) = bt
         tiles(tID)%uPlusX(idbdx,RunParams%nXpertile+1,1) = dbdx
         tiles(tID)%uMinusX(idbdx,RunParams%nXpertile+1,1) = dbdx
      end if
   end subroutine ComputeInterfacialTopographicData

   pure subroutine ComputeTopographicCurvatures(RunParams, grid, tiles, tID)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout), target :: grid
      type(tileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID

      integer :: idbdx, idbdy
      integer :: id2bdxx, id2bdyy, id2bdxy
      integer :: nX, nY
      integer :: ttW, ttE, ttS, ttN
      integer :: ttNW, ttNE, ttSW, ttSE

      real(kind=wp) :: deltaXRecip, deltaYRecip
      real(kind=wp), dimension(:,:), allocatable :: d2bdxx, d2bdyy, d2bdxy, d2bdyx

      real(kind=wp), dimension(:,:), allocatable :: dbdx, dbdy

      logical :: isActiveTile

      if (.not.RunParams%curvature) return

      isActiveTile = .not. tiles(tID)%isGhostTile

      if (.not. isActiveTile) return

      idbdx = RunParams%Vars%dbdx
      idbdy = RunParams%Vars%dbdy
      id2bdxx = RunParams%Vars%d2bdxx
      id2bdyy = RunParams%Vars%d2bdyy
      id2bdxy = RunParams%Vars%d2bdxy

      nX = RunParams%nXpertile
      nY = RunParams%nYpertile

      deltaXRecip = grid%deltaXRecip
      deltaYRecip = grid%deltaYRecip

      ttW = tiles(tID)%West
      ttE = tiles(tID)%East

      if (.not. RunParams%isOneD) then

         ! Allocate dbdx, dbdy -- take neighours 3 pixels wide on each edge to enable smoothing, if needed
         ! Note, the indices i=1...nXpertile correspond to this tile, others to neighbours
         allocate(dbdx(-2:RunParams%nXpertile+3,-2:RunParams%nYpertile+3), dbdy(-2:RunParams%nXpertile+3,-2:RunParams%nYpertile+3))

         ! computations on edges/corners
         ttS = tiles(tID)%South
         ttN = tiles(tID)%North
         ttSW = tiles(tID)%SouthWest
         ttSE = tiles(tID)%SouthEast
         ttNW = tiles(tID)%NorthWest
         ttNE = tiles(tID)%NorthEast

         !interior
         dbdx(1:nX,1:nY) = tiles(tID)%u(idbdx,:,:)
         dbdy(1:nX,1:nY) = tiles(tID)%u(idbdy,:,:)
         ! West
         dbdx(-2:0,1:nY) = tiles(ttW)%u(idbdx,nX-2:nX,:)
         dbdy(-2:0,1:nY) = tiles(ttW)%u(idbdy,nX-2:nX,:)
         ! East
         dbdx(nX+1:nX+3,1:nY) = tiles(ttE)%u(idbdx,1:3,:)
         dbdy(nX+1:nX+3,1:nY) = tiles(ttE)%u(idbdy,1:3,:)
         ! South
         dbdx(1:nX,-2:0) = tiles(ttS)%u(idbdx,:,nY-2:nY)
         dbdy(1:nX,-2:0) = tiles(ttS)%u(idbdy,:,nY-2:nY)
         ! North
         dbdx(1:nX,nY+1:nY+3) = tiles(ttN)%u(idbdx,:,1:3)
         dbdy(1:nX,nY+1:nY+3) = tiles(ttN)%u(idbdy,:,1:3)
         ! South-West
         dbdx(-2:0,-2:0) = tiles(ttSW)%u(idbdx,nX-2:nX,nY-2:nY)
         dbdy(-2:0,-2:0) = tiles(ttSW)%u(idbdy,nX-2:nX,nY-2:nY)
         ! North-West
         dbdx(-2:0,nY+1:nY+3) = tiles(ttNW)%u(idbdx,nX-2:nX,1:3)
         dbdy(-2:0,nY+1:nY+3) = tiles(ttNW)%u(idbdy,nX-2:nX,1:3)
         ! North-East
         dbdx(nX+1:nX+3,nY+1:nY+3) = tiles(ttNE)%u(idbdx,1:3,1:3)
         dbdy(nX+1:nX+3,nY+1:nY+3) = tiles(ttNE)%u(idbdy,1:3,1:3)
         ! South-East
         dbdx(nX+1:nX+3,-2:0) = tiles(ttSE)%u(idbdx,1:3,nY-2:nY)
         dbdy(nX+1:nX+3,-2:0) = tiles(ttSE)%u(idbdy,1:3,nY-2:nY)

         if (RunParams%nBlur>0) then
            dbdx = MultiBinom3Blur(dbdx, RunParams%nBlur)
            dbdy = MultiBinom3Blur(dbdy, RunParams%nBlur)
         end if

         allocate(d2bdxx(nX,nY), d2bdxy(nX,nY), d2bdyx(nX,nY), d2bdyy(nX,nY))

         d2bdxx(:,:) = (0.25_wp*(dbdx(2:nX+1,1:nY) - dbdx(0:nX-1,1:nY)) + 0.125_wp*(dbdx(3:nX+2,1:nY) - dbdx(-1:nX-2,1:nY))) * deltaXRecip
         d2bdyx(:,:) = (0.25_wp*(dbdy(2:nX+1,1:nY) - dbdy(0:nX-1,1:nY)) + 0.125_wp*(dbdy(3:nX+2,1:nY) - dbdy(-1:nX-2,1:nY))) * deltaXRecip
         d2bdyy(:,:) = (0.25_wp*(dbdy(1:nX,2:nY+1) - dbdy(1:nX,0:nY-1)) + 0.125_wp*(dbdy(1:nX,3:nY+2) - dbdy(1:nX,-1:nY-2))) * deltaYRecip
         d2bdxy(:,:) = (0.25_wp*(dbdx(1:nX,2:nY+1) - dbdx(1:nX,0:nY-1)) + 0.125_wp*(dbdx(1:nX,3:nY+2) - dbdx(1:nX,-1:nY-2))) * deltaYRecip

         tiles(tID)%u(id2bdxx,:,:) = d2bdxx(:,:)
         tiles(tID)%u(id2bdxy,:,:) = 0.5_wp*(d2bdyx(:,:) + d2bdxy(:,:))
         tiles(tID)%u(id2bdyy,:,:) = d2bdyy(:,:)

         deallocate(dbdx, dbdy)
         deallocate(d2bdxx, d2bdxy, d2bdyx, d2bdyy)
         
      else
         ! Allocate dbdx -- take neighours 3 pixels wide on each edge to enable smoothing, if needed
         ! Note, the indices i=1...nXpertile correspond to this tile, others to neighbours
         allocate(dbdx(-2:RunParams%nXpertile+3,-2:RunParams%nYpertile+3))

         !interior
         dbdx(1:nX,1) = tiles(tID)%u(idbdx,:,1)
         ! West
         dbdx(-2:0,1) = tiles(ttW)%u(idbdx,nX-2:nX,1)
         ! East
         dbdx(nX+1:nX+3,1) = tiles(ttE)%u(idbdx,1:3,1)

         if (RunParams%nBlur>0) then
            dbdx = MultiBinom3Blur(dbdx, RunParams%nBlur)
         end if

         tiles(tID)%u(id2bdxx,:,1) = (0.25_wp*(dbdx(2:nX+1,1) - dbdx(0:nX-1,1)) + 0.125_wp*(dbdx(3:nX+2,1) - dbdx(-1:nX-2,1))) * deltaXRecip
         deallocate(dbdx)

      end if

   end subroutine ComputeTopographicCurvatures

   pure subroutine ComputeMorphodynamicCurvatures(RunParams, grid, tiles, tID)
      implicit none
      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout), target :: grid
      type(tileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID

      integer :: idbdx, idbdy
      integer :: id2bdtx, id2bdty
      integer :: id2bdxx, id2bdyy, id2bdxy
      integer :: nX, nY
      integer :: ttW, ttE, ttS, ttN
      integer :: ttNW, ttNE, ttSW, ttSE

      real(kind=wp) :: deltaXRecip, deltaYRecip
      real(kind=wp), dimension(:,:), allocatable :: dbdx, dbdy
      real(kind=wp), dimension(:,:), allocatable :: dbdt
      real(kind=wp), dimension(:,:), allocatable :: d2bdtx, d2bdty
      real(kind=wp), dimension(:,:), allocatable :: d2bdxx, d2bdyy, d2bdxy, d2bdyx

      real(kind=wp), dimension(:,:), allocatable :: gam
      real(kind=wp) :: psib

      logical :: isActiveTile

      if (.not.RunParams%curvature) return

      isActiveTile = .not. tiles(tID)%isGhostTile

      if (.not. isActiveTile) return

      idbdx = RunParams%Vars%dbdx
      idbdy = RunParams%Vars%dbdy
      id2bdtx = RunParams%Vars%d2bdtx
      id2bdty = RunParams%Vars%d2bdty
      id2bdxx = RunParams%Vars%d2bdxx
      id2bdxy = RunParams%Vars%d2bdxy
      id2bdyy = RunParams%Vars%d2bdyy
            
      nX = RunParams%nXpertile
      nY = RunParams%nYpertile

      deltaXRecip = grid%deltaXRecip
      deltaYRecip = grid%deltaYRecip

      psib = 1.0_wp - RunParams%BedPorosity

      ttW = tiles(tID)%West
      ttE = tiles(tID)%East

      if (.not. RunParams%isOneD) then

         ! Allocate dbdx, dbdy -- take neighours 3 pixels wide on each edge to enable smoothing, if needed
         ! Note, the indices i=1...nXpertile correspond to this tile, others to neighbours
         allocate(dbdx(-2:RunParams%nXpertile+3,-2:RunParams%nYpertile+3))
         allocate(dbdy(-2:RunParams%nXpertile+3,-2:RunParams%nYpertile+3))
         allocate(dbdt(-2:RunParams%nXpertile+3,-2:RunParams%nYpertile+3))
         allocate(gam(-2:RunParams%nXpertile+3,-2:RunParams%nYpertile+3))

         ! computations on edges/corners
         ttS = tiles(tID)%South
         ttN = tiles(tID)%North
         ttSW = tiles(tID)%SouthWest
         ttSE = tiles(tID)%SouthEast
         ttNW = tiles(tID)%NorthWest
         ttNE = tiles(tID)%NorthEast

         !interior
         dbdx(1:nX,1:nY) = tiles(tID)%u(idbdx,:,:)
         dbdy(1:nX,1:nY) = tiles(tID)%u(idbdy,:,:)
         if (allocated(tiles(tID)%EminusD)) then
            dbdt(1:nX,1:nY) = -tiles(tID)%EminusD(:,:)
         else
            dbdt(1:nX,1:nY) = 0.0_wp
         end if
         
         ! West
         dbdx(-2:0,1:nY) = tiles(ttW)%u(idbdx,nX-2:nX,:)
         dbdy(-2:0,1:nY) = tiles(ttW)%u(idbdy,nX-2:nX,:)
         if (allocated(tiles(ttW)%EminusD)) then
            dbdt(-2:0,1:nY) = -tiles(ttW)%EminusD(nX-2:nX,:)
         else
            dbdt(-2:0,1:nY) = 0.0_wp
         end if
         
         ! East
         dbdx(nX+1:nX+3,1:nY) = tiles(ttE)%u(idbdx,1:3,:)
         dbdy(nX+1:nX+3,1:nY) = tiles(ttE)%u(idbdy,1:3,:)
         if (allocated(tiles(ttE)%EminusD)) then
            dbdt(nX+1:nX+3,1:nY) = -tiles(ttE)%EminusD(1:3,:)
         else
            dbdt(nX+1:nX+3,1:nY) = 0.0_wp
         end if

         ! South
         dbdx(1:nX,-2:0) = tiles(ttS)%u(idbdx,:,nY-2:nY)
         dbdy(1:nX,-2:0) = tiles(ttS)%u(idbdy,:,nY-2:nY)
         if (allocated(tiles(ttS)%EminusD)) then
            dbdt(1:nX,-2:0) = -tiles(ttS)%EminusD(:,nY-2:nY)
         else
            dbdt(1:nX,-2:0) = 0.0_wp
         end if
         
         ! North
         dbdx(1:nX,nY+1:nY+3) = tiles(ttN)%u(idbdx,:,1:3)
         dbdy(1:nX,nY+1:nY+3) = tiles(ttN)%u(idbdy,:,1:3)
         if (allocated(tiles(ttN)%EminusD)) then
            dbdt(1:nX,nY+1:nY+3) = -tiles(ttN)%EminusD(:,1:3)
         else
            dbdt(1:nX,nY+1:nY+3) = 0.0_wp
         end if
         
         ! South-West
         dbdx(-2:0,-2:0) = tiles(ttSW)%u(idbdx,nX-2:nX,nY-2:nY)
         dbdy(-2:0,-2:0) = tiles(ttSW)%u(idbdy,nX-2:nX,nY-2:nY)
         if (allocated(tiles(ttSW)%EminusD)) then
            dbdt(-2:0,-2:0) = -tiles(ttSW)%EminusD(nX-2:nX,nY-2:nY)
         else
            dbdt(-2:0,-2:0) = 0.0_wp
         end if
         
         ! North-West
         dbdx(-2:0,nY+1:nY+3) = tiles(ttNW)%u(idbdx,nX-2:nX,1:3)
         dbdy(-2:0,nY+1:nY+3) = tiles(ttNW)%u(idbdy,nX-2:nX,1:3)
         if (allocated(tiles(ttNW)%EminusD)) then
            dbdt(-2:0,nY+1:nY+3) = -tiles(ttNW)%EminusD(nX-2:nX,1:3)
         else
            dbdt(-2:0,nY+1:nY+3) = 0.0_wp
         end if
         
         ! North-East
         dbdx(nX+1:nX+3,nY+1:nY+3) = tiles(ttNE)%u(idbdx,1:3,1:3)
         dbdy(nX+1:nX+3,nY+1:nY+3) = tiles(ttNE)%u(idbdy,1:3,1:3)
         if (allocated(tiles(ttNE)%EminusD)) then
            dbdt(nX+1:nX+3,nY+1:nY+3) = -tiles(ttNE)%EminusD(1:3,1:3)
         else
            dbdt(nX+1:nX+3,nY+1:nY+3) = 0.0_wp
         end if
         
         ! South-East
         dbdx(nX+1:nX+3,-2:0) = tiles(ttSE)%u(idbdx,1:3,nY-2:nY)
         dbdy(nX+1:nX+3,-2:0) = tiles(ttSE)%u(idbdy,1:3,nY-2:nY)
         if (allocated(tiles(ttSE)%EminusD)) then
            dbdt(nX+1:nX+3,-2:0) = -tiles(ttSE)%EminusD(1:3,nY-2:nY)
         else
            dbdt(nX+1:nX+3,-2:0) = 0.0_wp
         end if

         ! Update dbdt to include factor gamma/psib
         gam = GeometricCorrectionFactor_gradin_array(dbdx, dbdy)
         dbdt(:,:) = dbdt(:,:) * gam(:,:)/psib

         if (RunParams%nBlur>0) then
            dbdx = MultiBinom3Blur(dbdx, RunParams%nBlur)
            dbdy = MultiBinom3Blur(dbdy, RunParams%nBlur)
            dbdt = MultiBinom3Blur(dbdt, RunParams%nBlur)
         end if

         allocate(d2bdtx(nX,nY), d2bdty(nX,nY))
         allocate(d2bdxx(nX,nY), d2bdxy(nX,nY), d2bdyx(nX,nY), d2bdyy(nX,nY))

         d2bdtx(:,:) = (0.25_wp*(dbdt(2:nX+1,1:nY) - dbdt(0:nX-1,1:nY)) + 0.125_wp*(dbdt(3:nX+2,1:nY) - dbdt(-1:nX-2,1:nY))) * deltaXRecip
         d2bdty(:,:) = (0.25_wp*(dbdt(1:nX,2:nY+1) - dbdt(1:nX,0:nY-1)) + 0.125_wp*(dbdt(1:nX,3:nY+2) - dbdt(1:nX,-1:nY-2))) * deltaYRecip

         d2bdxx(:,:) = (0.25_wp*(dbdx(2:nX+1,1:nY) - dbdx(0:nX-1,1:nY)) + 0.125_wp*(dbdx(3:nX+2,1:nY) - dbdx(-1:nX-2,1:nY))) * deltaXRecip
         d2bdyx(:,:) = (0.25_wp*(dbdy(2:nX+1,1:nY) - dbdy(0:nX-1,1:nY)) + 0.125_wp*(dbdy(3:nX+2,1:nY) - dbdy(-1:nX-2,1:nY))) * deltaXRecip
         d2bdyy(:,:) = (0.25_wp*(dbdy(1:nX,2:nY+1) - dbdy(1:nX,0:nY-1)) + 0.125_wp*(dbdy(1:nX,3:nY+2) - dbdy(1:nX,-1:nY-2))) * deltaYRecip
         d2bdxy(:,:) = (0.25_wp*(dbdx(1:nX,2:nY+1) - dbdx(1:nX,0:nY-1)) + 0.125_wp*(dbdx(1:nX,3:nY+2) - dbdx(1:nX,-1:nY-2))) * deltaYRecip

         tiles(tID)%u(id2bdtx,:,:) = d2bdtx(:,:)
         tiles(tID)%u(id2bdty,:,:) = d2bdty(:,:)

         tiles(tID)%u(id2bdxx,:,:) = d2bdxx(:,:)
         tiles(tID)%u(id2bdxy,:,:) = 0.5_wp*(d2bdyx(:,:) + d2bdxy(:,:))
         tiles(tID)%u(id2bdyy,:,:) = d2bdyy(:,:)

         deallocate(dbdx, dbdy, dbdt)
         deallocate(d2bdxx, d2bdxy, d2bdyx, d2bdyy, d2bdtx, d2bdty)
         
      else
         ! Allocate dbdx -- take neighours 3 pixels wide on each edge to enable smoothing, if needed
         ! Note, the indices i=1...nXpertile correspond to this tile, others to neighbours
         allocate(dbdx(-2:RunParams%nXpertile+3,1))
         allocate(dbdy(-2:RunParams%nXpertile+3,1)) ! Only needed to zero
         allocate(dbdt(-2:RunParams%nXpertile+3,1))

         dbdy(:,:) = 0.0_wp

         !interior
         dbdx(1:nX,1) = tiles(tID)%u(idbdx,:,1)
         if (allocated(tiles(tID)%EminusD)) then
            dbdt(1:nX,1) = -tiles(tID)%EminusD(:,1)
         else
            dbdt(1:nX,1) = 0.0_wp
         end if

         ! West
         dbdx(-2:0,1) = tiles(ttW)%u(idbdx,nX-2:nX,1)
         if (allocated(tiles(ttW)%EminusD)) then
            dbdt(-2:0,1) = -tiles(ttW)%EminusD(nX-2:nX,1)
         else
            dbdt(-2:0,1) = 0.0_wp
         end if

         ! East
         dbdx(nX+1:nX+3,1) = tiles(ttE)%u(idbdx,1:3,1)
         if (allocated(tiles(ttE)%EminusD)) then
            dbdt(nX+1:nX+3,1) = -tiles(ttE)%EminusD(1:3,1)
         else
            dbdt(nX+1:nX+3,1) = 0.0_wp
         end if

         ! Update dbdt to include factor gamma/psib
         gam = GeometricCorrectionFactor_gradin_array(dbdx, dbdy)
         dbdt(:,1) = dbdt(:,1) * gam(:,1)/psib

         if (RunParams%nBlur>0) then
            dbdx = MultiBinom3Blur(dbdx, RunParams%nBlur)
            dbdt = MultiBinom3Blur(dbdt, RunParams%nBlur)
         end if

         tiles(tID)%u(id2bdxx,:,1) = (0.25_wp*(dbdx(2:nX+1,1) - dbdx(0:nX-1,1)) + 0.125_wp*(dbdx(3:nX+2,1) - dbdx(-1:nX-2,1))) * deltaXRecip
         tiles(tID)%u(id2bdtx,:,1) = (0.25_wp*(dbdt(2:nX+1,1) - dbdt(0:nX-1,1)) + 0.125_wp*(dbdt(3:nX+2,1) - dbdt(-1:nX-2,1))) * deltaXRecip

         deallocate(dbdx, dbdy, dbdt)

      end if
   end subroutine ComputeMorphodynamicCurvatures

   ! For DEMS / SRTMS etc, GetHeights returns different values of b0 at the
   ! interface depending on which tile it's called from. These O(1e-10)
   ! discrepancies in interpolation can be enough to cause negative depths in
   ! the solver.
   ! Therefore, we adopt a convention of overwriting b0 in the right and top
   ! interfaces of each tile with its neighbour's values. Since every active
   ! tile is initialised with ghost tiles surrounding it, this routine should
   ! ensure that b0 and its derivatives are single-valued fields.
   pure subroutine EqualiseTopographicBoundaryData(RunParams, grid, tID)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout), target :: grid
      integer, intent(in) :: tID

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ghostTiles, activeTiles

      integer :: i, j, ib0, ttW, ttE, ttN, ttS, ttSW, ttNE
      logical :: SW_updated, NE_updated

      tileContainer => grid%tileContainer
      ghostTiles => grid%ghostTiles
      activeTiles => grid%activeTiles

      ib0 = RunParams%Vars%b0

      SW_updated = .false.
      NE_updated = .false.

      ! Check left tile
      ttW = tileContainer(tID)%West
      if (InVector(activeTiles%List, ttW) .or. IsActiveGhostTile(grid, ttW)) then
         do j = 1, RunParams%nYpertile
            tileContainer(ttW)%b0(RunParams%nXpertile + 1, j) = tileContainer(tID)%b0(1, j)
            tileContainer(ttW)%uPlusX(ib0, RunParams%nXpertile + 1, j) = tileContainer(tID)%uPlusX(ib0, 1, j)
            tileContainer(ttW)%uMinusX(ib0, RunParams%nXpertile + 1, j) = tileContainer(tID)%uMinusX(ib0, 1, j)
         end do
         call ComputeCellCentredTopographicData(RunParams, grid, tileContainer, ttW)
         if (.not. RunParams%isOneD) then
            ttSW = tileContainer(ttW)%South
            if (InVector(activeTiles%List, ttSW) .or. IsActiveGhostTile(grid, ttSW)) then
               tileContainer(ttSW)%b0(RunParams%nXpertile + 1, RunParams%nYpertile + 1) = &
                  tileContainer(tID)%b0(1, 1)
               SW_updated = .true.
               call ComputeCellCentredTopographicData(RunParams, grid, tileContainer, ttSW)
            end if
         end if
      end if
      ! Check right tile
      ttE = tileContainer(tID)%East
      if (InVector(activeTiles%List, ttE) .or. IsActiveGhostTile(grid, ttE)) then
         do j = 1, RunParams%nYpertile
            tileContainer(tID)%b0(RunParams%nXpertile + 1, j) = tileContainer(ttE)%b0(1, j)
            tileContainer(tID)%uPlusX(ib0, RunParams%nXpertile + 1, j) = tileContainer(ttE)%uPlusX(ib0, 1, j)
            tileContainer(tID)%uMinusX(ib0, RunParams%nXpertile + 1, j) = tileContainer(ttE)%uMinusX(ib0, 1, j)
         end do
         if (.not. RunParams%isOneD) then
            ttNE = tileContainer(ttE)%North
            if (InVector(activeTiles%List, ttNE) .or. IsActiveGhostTile(grid, ttNE)) then
               tileContainer(tID)%b0(RunParams%nXpertile + 1, RunParams%nYpertile + 1) = &
                  tileContainer(ttNE)%b0(1, 1)
               NE_updated = .true.
            end if
         end if
         call ComputeCellCentredTopographicData(RunParams, grid, tileContainer, tID)
      end if

      if (.not. RunParams%isOneD) then
         ! Check top tile
         ttN = tileContainer(tID)%North
         if (InVector(activeTiles%List, ttN) .or. IsActiveGhostTile(grid, ttN)) then
            do i = 1, RunParams%nXpertile
               tileContainer(tID)%b0(i, RunParams%nYpertile + 1) = tileContainer(ttN)%b0(i, 1)
               tileContainer(tID)%uPlusY(ib0, i, RunParams%nYpertile + 1) = tileContainer(ttN)%uPlusY(ib0, i, 1)
               tileContainer(tID)%uMinusY(ib0, i, RunParams%nYpertile + 1) = tileContainer(ttN)%uMinusY(ib0, i, 1)
            end do
            if (.not. NE_updated) then
               ttNE = tileContainer(ttN)%East
               if (InVector(activeTiles%List, ttNE) .or. IsActiveGhostTile(grid, ttNE)) then
                  tileContainer(tID)%b0(RunParams%nXpertile + 1, RunParams%nYpertile + 1) = &
                     tileContainer(ttNE)%b0(1, 1)
               end if
            end if
            call ComputeCellCentredTopographicData(RunParams, grid, tileContainer, tID)
         end if
         ! Check bottom tile
         ttS = tileContainer(tID)%South
         if (InVector(activeTiles%List, ttS) .or. IsActiveGhostTile(grid, ttS)) then
            do i = 1, RunParams%nXpertile
               tileContainer(ttS)%b0(i, RunParams%nYpertile + 1) = tileContainer(tID)%b0(i, 1)
               tileContainer(ttS)%uPlusY(ib0, i, RunParams%nYpertile + 1) = tileContainer(tID)%uPlusY(ib0, i, 1)
               tileContainer(ttS)%uMinusY(ib0, i, RunParams%nYpertile + 1) = tileContainer(tID)%uMinusY(ib0, i, 1)
            end do
            call ComputeCellCentredTopographicData(RunParams, grid, tileContainer, ttS)
            if (.not. SW_updated) then
               ttSW = tileContainer(ttS)%West
               if (InVector(activeTiles%List, ttSW) .or. IsActiveGhostTile(grid, ttSW)) then
                  tileContainer(ttSW)%b0(RunParams%nXpertile + 1, RunParams%nYpertile + 1) = &
                     tileContainer(tID)%b0(1, 1)
                  call ComputeCellCentredTopographicData(RunParams, grid, tileContainer, ttSW)
               end if
            end if
         end if

      end if ! isOneD

   end subroutine EqualiseTopographicBoundaryData

end module morphodynamic_rhs_module
