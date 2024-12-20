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
   use closures_module, only: GeometricCorrectionFactor_gradin, DragClosure
   use limiters_module, only: limiter

   implicit none

   private
   public :: CalculateMorphodynamicRHS
   public :: ComputeCellCentredTopographicData
   public :: ComputeTopographicCurvatures
   public :: ComputeInterfacialTopographicData
   public :: EqualiseTopographicBoundaryData

   interface ComputeCellCentredTopographicData
      module procedure :: ComputeCellCentredTopographicData_tileID, ComputeCellCentredTopographicData_ij
   end interface

   interface ComputeTopographicCurvatures
      module procedure :: ComputeTopographicCurvatures_tileID, ComputeTopographicCurvatures_ij
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

      integer :: tt, tID, ttW, ttE, ttS, ttN
      integer :: i, j, nd, bt, iHn

      activeTiles => grid%activeTiles
      nd = size(RunParams%iFlux)
      iHn = RunParams%Vars%Hn
      bt = RunParams%Vars%bt

      Hneps = RunParams%heightThreshold
      HnS = 0.0_wp
      HnN = 0.0_wp

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
                  Friction = DragClosure(RunParams, tiles(tID)%u(:,i,j))
                  call ErosionDepositionTerms(RunParams, tiles(tID)%u(:,i,j), Ero, Depo)
                  tiles(tID)%EminusD(i,j) = Ero - Depo
               end if
            end do
         end do
      end do

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

      do tt = 1, ActiveTiles%Size
         tID = ActiveTiles%List(tt)
         if (.not. RunParams%isOneD) then
            ! interior
            do i = 2, nXpertile
               do j = 2, nYpertile
                  dbdx = 0.25_wp * KahanSum([tiles(tID)%u(idbdx, i - 1, j - 1), tiles(tID)%u(idbdx, i - 1, j), &
                     tiles(tID)%u(idbdx, i    , j - 1), tiles(tID)%u(idbdx, i    , j)])
                  dbdy = 0.25_wp * KahanSum([tiles(tID)%u(idbdy, i - 1, j - 1), tiles(tID)%u(idbdy, i - 1, j), &
                     tiles(tID)%u(idbdy, i    , j - 1), tiles(tID)%u(idbdy, i    , j)])
                  gam = GeometricCorrectionFactor_gradin(dbdx, dbdy)
                  tiles(tID)%ddtExplicitBt(i,j) = -0.25_wp * gam / psib *  &
                     KahanSum([tiles(tID)%EminusD(i - 1, j - 1), tiles(tID)%EminusD(i - 1, j), tiles(tID)%EminusD(i, j - 1), tiles(tID)%EminusD(i, j)])
               end do
            end do

            ! computations on edges/corners
            ttW = tiles(tID)%West
            ttE = tiles(tID)%East
            ttS = tiles(tID)%South
            ttN = tiles(tID)%North
            ttSW = tiles(ttW)%South
            ttSE = tiles(ttE)%South
            ttNW = tiles(ttW)%North
            ttNE = tiles(ttE)%North

            do j = 2, nYpertile
               ! west bdry
               dbdx = 0.25_wp * KahanSum([tiles(ttW)%u(idbdx, nXpertile, j - 1), tiles(ttW)%u(idbdx, nXpertile, j), &
                  tiles(tID)%u(idbdx, 1, j - 1), tiles(tID)%u(idbdx, 1, j)])
               dbdy = 0.25_wp * KahanSum([tiles(ttW)%u(idbdy, nXpertile, j - 1), tiles(ttW)%u(idbdy, nXpertile, j), &
                  tiles(tID)%u(idbdy, 1, j - 1), tiles(tID)%u(idbdy, 1, j)])
               gam = GeometricCorrectionFactor_gradin(dbdx, dbdy)
               tiles(tID)%ddtExplicitBt(1, j) = -0.25_wp * gam / psib *  &
                  KahanSum([tiles(ttW)%EminusD(nXpertile, j - 1), tiles(ttW)%EminusD(nXpertile, j), tiles(tID)%EminusD(1, j - 1) + tiles(tID)%EminusD(1, j)]) 
               ! east bdry
               dbdx = 0.25_wp * KahanSum([tiles(tID)%u(idbdx, nXpertile, j - 1), tiles(tID)%u(idbdx, nXpertile, j), &
                  tiles(ttE)%u(idbdx, 1, j - 1), tiles(ttE)%u(idbdx, 1, j)])
               dbdy = 0.25_wp * KahanSum([tiles(tID)%u(idbdy, nXpertile, j - 1), tiles(tID)%u(idbdy, nXpertile, j), &
                  tiles(ttE)%u(idbdy, 1, j - 1), tiles(ttE)%u(idbdy, 1, j)])
               gam = GeometricCorrectionFactor_gradin(dbdx, dbdy)
               tiles(tID)%ddtExplicitBt(nXpertile + 1, j) = -0.25_wp * gam / psib *  &
                  KahanSum([tiles(tID)%EminusD(nXpertile, j - 1), tiles(tID)%EminusD(nXpertile, j), tiles(ttE)%EminusD(1, j - 1), tiles(ttE)%EminusD(1, j)])
            end do
            do i = 2, nXpertile
               ! north bdry
               dbdx = 0.25_wp * KahanSum([tiles(tID)%u(idbdx, i - 1, nYpertile), tiles(ttN)%u(idbdx, i - 1, 1), &
                  tiles(tID)%u(idbdx, i, nYpertile), tiles(ttN)%u(idbdx, i, 1)])
               dbdy = 0.25_wp * KahanSum([tiles(tID)%u(idbdy, i - 1, nYpertile), tiles(ttN)%u(idbdy, i - 1, 1), &
                  tiles(tID)%u(idbdy, i, nYpertile), tiles(ttN)%u(idbdy, i, 1)])
               gam = GeometricCorrectionFactor_gradin(dbdx, dbdy)
               tiles(tID)%ddtExplicitBt(i, nYpertile + 1) = -0.25_wp * gam / psib *  &
                  KahanSum([tiles(tID)%EminusD(i - 1, nYpertile), tiles(ttN)%EminusD(i - 1, 1), tiles(tID)%EminusD(i, nYpertile), tiles(ttN)%EminusD(i, 1)])
               ! south bdry
               dbdx = 0.25_wp * KahanSum([tiles(ttS)%u(idbdx, i - 1, nYpertile), tiles(tID)%u(idbdx, i - 1, 1), &
                  tiles(ttS)%u(idbdx, i, nYpertile), tiles(tID)%u(idbdx, i, 1)])
               dbdy = 0.25_wp * KahanSum([tiles(ttS)%u(idbdy, i - 1, nYpertile), tiles(tID)%u(idbdy, i - 1, 1), &
                  tiles(ttS)%u(idbdy, i, nYpertile), tiles(tID)%u(idbdy, i, 1)])
               gam = GeometricCorrectionFactor_gradin(dbdx, dbdy)
               tiles(tID)%ddtExplicitBt(i, 1) = -0.25_wp * gam / psib * &
                  KahanSum([tiles(ttS)%EminusD(i - 1, nYpertile), tiles(tID)%EminusD(i - 1, 1), tiles(ttS)%EminusD(i, nYpertile), tiles(tID)%EminusD(i, 1)])
            end do
            ! corners
            ! SW
            dbdx = 0.25_wp * KahanSum([tiles(ttSW)%u(idbdx, nXpertile, nYpertile), tiles(ttW)%u(idbdx, nXpertile, 1), &
               tiles(ttS )%u(idbdx, 1, nYpertile), tiles(tID)%u(idbdx, 1, 1)])
            dbdy = 0.25_wp * KahanSum([tiles(ttSW)%u(idbdy, nXpertile, nYpertile), tiles(ttW)%u(idbdy, nXpertile, 1), &
               tiles(ttS )%u(idbdy, 1, nYpertile), tiles(tID)%u(idbdy, 1, 1)])
            gam = GeometricCorrectionFactor_gradin(dbdx, dbdy)
            tiles(tID)%ddtExplicitBt(1, 1) = -0.25_wp * gam / psib *  &
               KahanSum([tiles(ttSW)%EminusD(nXpertile, nYpertile), tiles(ttW)%EminusD(nXpertile, 1), tiles(ttS )%EminusD(1, nYpertile), tiles(tID)%EminusD(1, 1)])
            ! SE
            dbdx = 0.25_wp * KahanSum([tiles(ttS )%u(idbdx, nXpertile, nYpertile), tiles(tID)%u(idbdx, nXpertile, 1), &
               tiles(ttSE)%u(idbdx, 1, nYpertile), tiles(ttE)%u(idbdx, 1, 1)])
            dbdy = 0.25_wp * KahanSum([tiles(ttS )%u(idbdy, nXpertile, nYpertile), tiles(tID)%u(idbdy, nXpertile, 1), &
               tiles(ttSE)%u(idbdy, 1, nYpertile), tiles(ttE)%u(idbdy, 1, 1)])
            gam = GeometricCorrectionFactor_gradin(dbdx, dbdy)
            tiles(tID)%ddtExplicitBt(nXpertile + 1, 1) = -0.25_wp * gam / psib *  &
               KahanSum([tiles(ttS )%EminusD(nXpertile, nYpertile), tiles(tID)%EminusD(nXpertile, 1), tiles(ttSE)%EminusD(1, nYpertile), tiles(ttE)%EminusD(1, 1)])
            ! NW
            dbdx = 0.25_wp * KahanSum([tiles(ttW)%u(idbdx, nXpertile, nYpertile), tiles(ttNW)%u(idbdx, nXpertile, 1), &
               tiles(tID)%u(idbdx, 1, nYpertile), tiles(ttN )%u(idbdx, 1, 1)])
            dbdy = 0.25_wp * KahanSum([tiles(ttW)%u(idbdy, nXpertile, nYpertile), tiles(ttNW)%u(idbdy, nXpertile, 1), &
               tiles(tID)%u(idbdy, 1, nYpertile), tiles(ttN )%u(idbdy, 1, 1)])
            gam = GeometricCorrectionFactor_gradin(dbdx, dbdy)
            tiles(tID)%ddtExplicitBt(1, nYpertile + 1) = -0.25_wp * gam / psib *  &
               KahanSum([tiles(ttW)%EminusD(nXpertile, nYpertile), tiles(ttNW)%EminusD(nXpertile, 1), &
               tiles(tID)%EminusD(1, nYpertile), tiles(ttN )%EminusD(1, 1)])
            ! NE
            dbdx = 0.25_wp * KahanSum([tiles(tID)%u(idbdx, nXpertile, nYpertile), tiles(ttN )%u(idbdx, nXpertile, 1), &
               tiles(ttE)%u(idbdx, 1, nYpertile) + tiles(ttNE)%u(idbdx, 1, 1)])
            dbdy = 0.25_wp * KahanSum([tiles(tID)%u(idbdy, nXpertile, nYpertile), tiles(ttN )%u(idbdy, nXpertile, 1), &
               tiles(ttE)%u(idbdy, 1, nYpertile), tiles(ttNE)%u(idbdy, 1, 1)])
            gam = GeometricCorrectionFactor_gradin(dbdx, dbdy)
            tiles(tID)%ddtExplicitBt(nXpertile + 1, nYpertile + 1) = -0.25_wp * gam / psib *  &
               KahanSum([tiles(tID)%EminusD(nXpertile, nYpertile), tiles(ttN )%EminusD(nXpertile, 1), &
               tiles(ttE)%EminusD(1, nYpertile), tiles(ttNE)%EminusD(1, 1)])
         else
            ! 1D is much simpler...
            do i = 2, nXpertile
               dbdx = 0.5_wp * (tiles(tID)%u(idbdx, i - 1, 1) + tiles(tID)%u(idbdx, i, 1))
               gam = GeometricCorrectionFactor_gradin(dbdx, 0.0_wp)
               tiles(tID)%ddtExplicitBt(i,1) = -0.5_wp * gam * (tiles(tID)%EminusD(i - 1, 1) + tiles(tID)%EminusD(i, 1)) / psib
            end do

            ! handle tile boundaries
            ttW = tiles(tID)%West
            if (ttW > 0 .and. grid%tileContainer(ttW)%TileOn) then
               dbdx = 0.5_wp * (tiles(ttW)%u(iDBDX, nXpertile, 1) + tiles(tID)%u(iDBDX, 1, 1))
               gam = GeometricCorrectionFactor_gradin(dbdx, 0.0_wp)
               tiles(tID)%ddtExplicitBt(1,1) = -0.5_wp * gam * (tiles(ttW)%EminusD(nXpertile, 1) + tiles(tID)%EminusD(1, 1)) / psib
            else
               dbdx = 0.5_wp * tiles(tID)%u(idbdx, 1, 1)
               gam = GeometricCorrectionFactor_gradin(dbdx, 0.0_wp)
               tiles(tID)%ddtExplicitBt(1,1) = -0.5_wp * gam * tiles(tID)%EminusD(1, 1) / psib
            end if

            ttE = tiles(tID)%East
            if (ttE > 0 .and. grid%tileContainer(ttE)%TileOn) then
               dbdx = 0.5_wp * (tiles(tID)%u(idbdx, nXpertile, 1) + tiles(ttE)%u(idbdx, 1, 1))
               gam = GeometricCorrectionFactor_gradin(dbdx, 0.0_wp)
               tiles(tID)%ddtExplicitBt(nXpertile+1,1) = -0.5_wp * gam * (tiles(tID)%EminusD(nXpertile, 1) + tiles(ttE)%EminusD(1, 1)) / psib
            else
               dbdx = 0.5_wp * tiles(tID)%u(idbdx, nXpertile, 1)
               gam = GeometricCorrectionFactor_gradin(dbdx, 0.0_wp)
               tiles(tID)%ddtExplicitBt(nXpertile+1,1) = -0.5_wp * gam * tiles(tID)%EminusD(nXpertile, 1) / psib
            end if
         end if
      end do
   end subroutine

   ! Interpolate b at cell centres, for tile tID.
   subroutine ComputeCellCentredTopographicData_tileID(RunParams, grid, tiles, tID)
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
   subroutine ComputeCellCentredTopographicData_ij(RunParams, grid, tiles, tID, i, j)
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
   subroutine ComputeInterfacialTopographicData(RunParams, grid, tiles, tID)
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

   subroutine ComputeTopographicCurvatures_tileID(RunParams, grid, tiles, tID)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout), target :: grid
      type(tileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID

      integer :: i, j
      integer :: ib
      integer :: idbdx, idbdy
      integer :: id2bdxx, id2bdyy, id2bdxy
      integer :: nX, nY
      integer :: ttW, ttE, ttS, ttN
      integer :: ttNW, ttNE, ttSW, ttSE

      real(kind=wp) :: deltaXRecip, deltaYRecip
      real(kind=wp) :: d2bdxx, d2bdyy, d2bdxy, d2bdyx

      logical :: isActiveTile

      if (.not.RunParams%curvature) return

      isActiveTile = .not. tiles(tID)%isGhostTile

      if (.not. isActiveTile) return

      ib = RunParams%Vars%b0
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

         ! interior
         do i = 2, nX-1
            do j = 2, nY-1

               d2bdxx = deltaXRecip*deltaXRecip * (tiles(tID)%u(ib,i+1,j) - 2.0_wp*tiles(tID)%u(ib,i,j) + tiles(tID)%u(ib,i-1,j))
               d2bdxy = 0.25_wp*deltaXRecip*deltaYRecip * (tiles(tID)%u(ib,i+1,j+1) - tiles(tID)%u(ib,i+1,j-1) - tiles(tID)%u(ib,i-1,j+1) + tiles(tID)%u(ib,i-1,j-1))
               d2bdyy = deltaYRecip*deltaYRecip * (tiles(tID)%u(ib,i,j+1) - 2.0_wp*tiles(tID)%u(ib,i,j) + tiles(tID)%u(ib,i,j-1))

               ! d2bdxx = deltaXRecip * limiter(tiles(tID)%u(idbdx,i+1,j)-tiles(tID)%u(idbdx,i,j), tiles(tID)%u(idbdx,i,j)-tiles(tID)%u(idbdx,i-1,j))
               ! d2bdyx = deltaXRecip * limiter(tiles(tID)%u(idbdy,i+1,j)-tiles(tID)%u(idbdy,i,j), tiles(tID)%u(idbdy,i,j)-tiles(tID)%u(idbdy,i-1,j))
               ! d2bdxy = deltaYRecip * limiter(tiles(tID)%u(idbdx,i,j+1)-tiles(tID)%u(idbdx,i,j), tiles(tID)%u(idbdx,i,j)-tiles(tID)%u(idbdx,i,j-1))
               ! d2bdyy = deltaYRecip * limiter(tiles(tID)%u(idbdy,i,j+1)-tiles(tID)%u(idbdy,i,j), tiles(tID)%u(idbdy,i,j)-tiles(tID)%u(idbdy,i,j-1))

               ! d2bdxx = ( tiles(tID)%u(idbdx, i+1, j+1) - tiles(tID)%u(idbdx, i-1, j+1) &
               !          + tiles(tID)%u(idbdx, i+1, j-1) - tiles(tID)%u(idbdx, i-1, j-1) ) * 0.25_wp * deltaXRecip
               ! d2bdyx = ( tiles(tID)%u(idbdy, i+1, j+1) - tiles(tID)%u(idbdy, i-1, j+1) &
               !          + tiles(tID)%u(idbdy, i+1, j-1) - tiles(tID)%u(idbdy, i-1, j-1) ) * 0.25_wp * deltaXRecip
               ! d2bdxy = ( tiles(tID)%u(idbdx, i+1, j+1) - tiles(tID)%u(idbdx, i+1, j-1) &
               !          + tiles(tID)%u(idbdx, i-1, j+1) - tiles(tID)%u(idbdx, i-1, j-1) ) * 0.25_wp * deltaYRecip
               ! d2bdyy = ( tiles(tID)%u(idbdy, i+1, j+1) - tiles(tID)%u(idbdy, i+1, j-1) &
               !          + tiles(tID)%u(idbdy, i-1, j+1) - tiles(tID)%u(idbdy, i-1, j-1) ) * 0.25_wp * deltaYRecip
               tiles(tID)%u(id2bdxx, i, j) = d2bdxx
               tiles(tID)%u(id2bdxy, i, j) = d2bdxy !0.5_wp*(d2bdxy + d2bdyx)
               tiles(tID)%u(id2bdyy, i, j) = d2bdyy
            end do
         end do

         ! computations on edges/corners
         ttS = tiles(tID)%South
         ttN = tiles(tID)%North
         ttSW = tiles(tID)%SouthWest
         ttSE = tiles(tID)%SouthEast
         ttNW = tiles(tID)%NorthWest
         ttNE = tiles(tID)%NorthEast

         do j = 2, nY-1
            ! west bdry: i = 1
            d2bdxx = deltaXRecip*deltaXRecip * (tiles(tID)%u(ib,2,j) - 2.0_wp*tiles(tID)%u(ib,1,j) + tiles(ttW)%u(ib,nX,j))
            d2bdxy = 0.25_wp*deltaXRecip*deltaYRecip * (tiles(tID)%u(ib,2,j+1) - tiles(tID)%u(ib,2,j-1) - tiles(ttW)%u(ib,nX,j+1) + tiles(ttW)%u(ib,nX,j-1))
            d2bdyy = deltaYRecip*deltaYRecip * (tiles(tID)%u(ib,1,j+1) - 2.0_wp*tiles(tID)%u(ib,1,j) + tiles(tID)%u(ib,1,j-1))

               ! d2bdxx = deltaXRecip * limiter(tiles(tID)%u(idbdx,2,j)-tiles(tID)%u(idbdx,1,j), tiles(tID)%u(idbdx,1,j)-tiles(ttW)%u(idbdx,nX,j))
               ! d2bdyx = deltaXRecip * limiter(tiles(tID)%u(idbdy,2,j)-tiles(tID)%u(idbdy,1,j), tiles(tID)%u(idbdy,1,j)-tiles(ttW)%u(idbdy,nX,j))
               ! d2bdxy = deltaYRecip * limiter(tiles(tID)%u(idbdx,1,j+1)-tiles(tID)%u(idbdx,1,j), tiles(tID)%u(idbdx,1,j)-tiles(tID)%u(idbdx,1,j-1))
               ! d2bdyy = deltaYRecip * limiter(tiles(tID)%u(idbdy,1,j+1)-tiles(tID)%u(idbdy,1,j), tiles(tID)%u(idbdy,1,j)-tiles(tID)%u(idbdy,1,j-1))

               ! d2bdxx = ( tiles(tID)%u(idbdx, 2, j+1) - tiles(ttW)%u(idbdx, nX, j+1) &
               !          + tiles(tID)%u(idbdx, 2, j-1) - tiles(ttW)%u(idbdx, nX, j-1) ) * 0.25_wp * deltaXRecip
               ! d2bdyx = ( tiles(tID)%u(idbdy, 2, j+1) - tiles(ttW)%u(idbdy, nX, j+1) &
               !          + tiles(tID)%u(idbdy, 2, j-1) - tiles(ttW)%u(idbdy, nX, j-1) ) * 0.25_wp * deltaXRecip
               ! d2bdxy = ( tiles(tID)%u(idbdx,  2, j+1) - tiles(tID)%u(idbdx,  2, j-1) &
               !          + tiles(ttW)%u(idbdx, nX, j+1) - tiles(ttW)%u(idbdx, nX, j-1) ) * 0.25_wp * deltaYRecip
               ! d2bdyy = ( tiles(tID)%u(idbdy,  2, j+1) - tiles(tID)%u(idbdy,  2, j-1) &
               !          + tiles(ttW)%u(idbdy, nX, j+1) - tiles(ttW)%u(idbdy, nX, j-1) ) * 0.25_wp * deltaYRecip
            
            tiles(tID)%u(id2bdxx, 1, j) = d2bdxx
            tiles(tID)%u(id2bdxy, 1, j) = d2bdxy !0.5_wp*(d2bdxy+d2bdyx)
            tiles(tID)%u(id2bdyy, 1, j) = d2bdyy

            ! east bdry: i = nX
            d2bdxx = deltaXRecip*deltaXRecip * (tiles(ttE)%u(ib,1,j) - 2.0_wp*tiles(tID)%u(ib,nX,j) + tiles(tID)%u(ib,nX-1,j))
            d2bdxy = 0.25_wp*deltaXRecip*deltaYRecip * (tiles(ttE)%u(ib,1,j+1) - tiles(ttE)%u(ib,1,j-1) - tiles(tID)%u(ib,nX-1,j+1) + tiles(tID)%u(ib,nX-1,j-1))
            d2bdyy = deltaYRecip*deltaYRecip * (tiles(tID)%u(ib,nX,j+1) - 2.0_wp*tiles(tID)%u(ib,nX,j) + tiles(tID)%u(ib,nX,j-1))
            
               ! d2bdxx = deltaXRecip * limiter(tiles(ttE)%u(idbdx,1,j)-tiles(tID)%u(idbdx,nX,j), tiles(tID)%u(idbdx,nX,j)-tiles(tID)%u(idbdx,nX-1,j))
               ! d2bdyx = deltaXRecip * limiter(tiles(ttE)%u(idbdy,1,j)-tiles(tID)%u(idbdy,nX,j), tiles(tID)%u(idbdy,nX,j)-tiles(tID)%u(idbdy,nX-1,j))
               ! d2bdxy = deltaYRecip * limiter(tiles(tID)%u(idbdx,nX,j+1)-tiles(tID)%u(idbdx,nX,j), tiles(tID)%u(idbdx,nX,j)-tiles(tID)%u(idbdx,nX,j-1))
               ! d2bdyy = deltaYRecip * limiter(tiles(tID)%u(idbdy,nX,j+1)-tiles(tID)%u(idbdy,nX,j), tiles(tID)%u(idbdy,nX,j)-tiles(tID)%u(idbdy,nX,j-1))

               ! d2bdxx = ( tiles(ttE)%u(idbdx, 1, j+1) - tiles(tID)%u(idbdx, nX-1, j+1) &
               !          + tiles(ttE)%u(idbdx, 1, j-1) - tiles(tID)%u(idbdx, nX-1, j-1) ) * 0.25_wp * deltaXRecip
               ! d2bdyx = ( tiles(ttE)%u(idbdy, 1, j+1) - tiles(tID)%u(idbdy, nX-1, j+1) &
               !          + tiles(ttE)%u(idbdy, 1, j-1) - tiles(tID)%u(idbdy, nX-1, j-1) ) * 0.25_wp * deltaXRecip
               ! d2bdxy = ( tiles(ttE)%u(idbdx,    1, j+1) - tiles(ttE)%u(idbdx,    1, j-1) &
               !          + tiles(tID)%u(idbdx, nX-1, j+1) - tiles(tID)%u(idbdx, nX-1, j-1) ) * 0.25_wp * deltaYRecip
               ! d2bdyy = ( tiles(ttE)%u(idbdy,    1, j+1) - tiles(ttE)%u(idbdy,    1, j-1) &
               !          + tiles(tID)%u(idbdy, nX-1, j+1) - tiles(tID)%u(idbdy, nX-1, j-1) ) * 0.25_wp * deltaYRecip
            
            tiles(tID)%u(id2bdxx, Nx, j) = d2bdxx
            tiles(tID)%u(id2bdxy, Nx, j) = d2bdxy !0.5_wp*(d2bdxy+d2bdyx)
            tiles(tID)%u(id2bdyy, Nx, j) = d2bdyy
         end do
         
         do i = 2, nX-1
            ! north bdry: j = nY
            d2bdxx = deltaXRecip*deltaXRecip * (tiles(tID)%u(ib,i+1,nY) - 2.0_wp*tiles(tID)%u(ib,i,nY) + tiles(tID)%u(ib,i-1,nY))
            d2bdxy = 0.25_wp*deltaXRecip*deltaYRecip * (tiles(ttN)%u(ib,i+1,1) - tiles(tID)%u(ib,i+1,nY-1) - tiles(ttN)%u(ib,i-1,1) + tiles(tID)%u(ib,i-1,nY-1))
            d2bdyy = deltaYRecip*deltaYRecip * (tiles(ttN)%u(ib,i,1) - 2.0_wp*tiles(tID)%u(ib,i,nY) + tiles(tID)%u(ib,i,nY-1))
            
               ! d2bdxx = deltaXRecip * limiter(tiles(tID)%u(idbdx,i+1,nY)-tiles(tID)%u(idbdx,i,nY), tiles(tID)%u(idbdx,i,nY)-tiles(tID)%u(idbdx,i-1,nY))
               ! d2bdyx = deltaXRecip * limiter(tiles(tID)%u(idbdy,i+1,nY)-tiles(tID)%u(idbdy,i,nY), tiles(tID)%u(idbdy,i,nY)-tiles(tID)%u(idbdy,i-1,nY))
               ! d2bdxy = deltaYRecip * limiter(tiles(ttN)%u(idbdx,i,1)-tiles(tID)%u(idbdx,i,nY), tiles(tID)%u(idbdx,i,nY)-tiles(tID)%u(idbdx,i,nY-1))
               ! d2bdyy = deltaYRecip * limiter(tiles(ttN)%u(idbdy,i,1)-tiles(tID)%u(idbdy,i,nY), tiles(tID)%u(idbdy,i,nY)-tiles(tID)%u(idbdy,i,nY-1))

               ! d2bdxx = ( tiles(ttN)%u(idbdx, i+1,    1) - tiles(ttN)%u(idbdx, i-1,    1) &
               !          + tiles(tID)%u(idbdx, i+1, nY-1) - tiles(tID)%u(idbdx, i-1, nY-1) ) * 0.25_wp * deltaXRecip
               ! d2bdyx = ( tiles(ttN)%u(idbdy, i+1,    1) - tiles(ttN)%u(idbdy, i-1,    1) &
               !          + tiles(tID)%u(idbdy, i+1, nY-1) - tiles(tID)%u(idbdy, i-1, nY-1) ) * 0.25_wp * deltaXRecip
               ! d2bdxy = ( tiles(ttN)%u(idbdx, i+1, 1) - tiles(tID)%u(idbdx, i+1, nY-1) &
               !          + tiles(ttN)%u(idbdx, i-1, 1) - tiles(tID)%u(idbdx, i-1, nY-1) ) * 0.25_wp * deltaYRecip
               ! d2bdyy = ( tiles(ttN)%u(idbdy, i+1, 1) - tiles(tID)%u(idbdy, i+1, nY-1) &
               !          + tiles(ttN)%u(idbdy, i-1, 1) - tiles(tID)%u(idbdy, i-1, nY-1) ) * 0.25_wp * deltaYRecip
            
            tiles(tID)%u(id2bdxx, i, nY) = d2bdxx
            tiles(tID)%u(id2bdxy, i, nY) = d2bdxy !0.5_wp*(d2bdxy+d2bdyx)
            tiles(tID)%u(id2bdyy, i, nY) = d2bdyy
            
            ! south bdry: j = 1
            d2bdxx = deltaXRecip*deltaXRecip * (tiles(tID)%u(ib,i+1,1) - 2.0_wp*tiles(tID)%u(ib,i,1) + tiles(tID)%u(ib,i-1,1))
            d2bdxy = 0.25_wp*deltaXRecip*deltaYRecip * (tiles(tID)%u(ib,i+1,2) - tiles(ttS)%u(ib,i+1,nY) - tiles(tID)%u(ib,i-1,2) + tiles(ttS)%u(ib,i-1,nY))
            d2bdyy = deltaYRecip*deltaYRecip * (tiles(tID)%u(ib,i,2) - 2.0_wp*tiles(tID)%u(ib,i,1) + tiles(ttS)%u(ib,i,nY))
            
               ! d2bdxx = deltaXRecip * limiter(tiles(tID)%u(idbdx,i+1,1)-tiles(tID)%u(idbdx,i,1), tiles(tID)%u(idbdx,i,1)-tiles(tID)%u(idbdx,i-1,1))
               ! d2bdyx = deltaXRecip * limiter(tiles(tID)%u(idbdy,i+1,1)-tiles(tID)%u(idbdy,i,1), tiles(tID)%u(idbdy,i,1)-tiles(tID)%u(idbdy,i-1,1))
               ! d2bdxy = deltaYRecip * limiter(tiles(tID)%u(idbdx,i,2)-tiles(tID)%u(idbdx,i,1), tiles(tID)%u(idbdx,i,1)-tiles(ttS)%u(idbdx,i,nY))
               ! d2bdyy = deltaYRecip * limiter(tiles(tID)%u(idbdy,i,2)-tiles(tID)%u(idbdy,i,1), tiles(tID)%u(idbdy,i,1)-tiles(ttS)%u(idbdy,i,nY))

               ! d2bdxx = ( tiles(tID)%u(idbdx, i+1,  2) - tiles(tID)%u(idbdx, i-1,  2) &
               !          + tiles(ttS)%u(idbdx, i+1, nY) - tiles(ttS)%u(idbdx, i-1, nY) ) * 0.25_wp * deltaXRecip
               ! d2bdyx = ( tiles(tID)%u(idbdy, i+1,  2) - tiles(tID)%u(idbdy, i-1,  2) &
               !          + tiles(ttS)%u(idbdy, i+1, nY) - tiles(ttS)%u(idbdy, i-1, nY) ) * 0.25_wp * deltaXRecip
               ! d2bdxy = ( tiles(tID)%u(idbdx, i+1, 2) - tiles(ttS)%u(idbdx, i+1, nY) &
               !          + tiles(tID)%u(idbdx, i-1, 2) - tiles(ttS)%u(idbdx, i-1, nY) ) * 0.25_wp * deltaYRecip
               ! d2bdyy = ( tiles(tID)%u(idbdy, i+1, 2) - tiles(ttS)%u(idbdy, i+1, nY) &
               !          + tiles(tID)%u(idbdy, i-1, 2) - tiles(ttS)%u(idbdy, i-1, nY) ) * 0.25_wp * deltaYRecip
            
            tiles(tID)%u(id2bdxx, i, 1) = d2bdxx
            tiles(tID)%u(id2bdxy, i, 1) = d2bdxy !0.5_wp*(d2bdxy+d2bdyx)
            tiles(tID)%u(id2bdyy, i, 1) = d2bdyy
         end do
         ! corners
         ! SW: i=1, j=1
         d2bdxx = deltaXRecip*deltaXRecip * (tiles(tID)%u(ib,2,1) - 2.0_wp*tiles(tID)%u(ib,1,1) + tiles(ttW)%u(ib,nX,1))
         d2bdxy = 0.25_wp*deltaXRecip*deltaYRecip * (tiles(tID)%u(ib,2,2) - tiles(ttS)%u(ib,2,nY) - tiles(ttW)%u(ib,nX,2) + tiles(ttSW)%u(ib,nX,nY))
         d2bdyy = deltaYRecip*deltaYRecip * (tiles(tID)%u(ib,1,2) - 2.0_wp*tiles(tID)%u(ib,1,1) + tiles(ttS)%u(ib,1,nY))
         
            ! d2bdxx = deltaXRecip * limiter(tiles(tID)%u(idbdx,2,1)-tiles(tID)%u(idbdx,1,1), tiles(tID)%u(idbdx,1,1)-tiles(ttW)%u(idbdx,nX,1))
            ! d2bdyx = deltaXRecip * limiter(tiles(tID)%u(idbdy,2,1)-tiles(tID)%u(idbdy,1,1), tiles(tID)%u(idbdy,1,1)-tiles(ttW)%u(idbdy,nX,j))
            ! d2bdxy = deltaYRecip * limiter(tiles(tID)%u(idbdx,1,2)-tiles(tID)%u(idbdx,1,1), tiles(tID)%u(idbdx,1,1)-tiles(ttS)%u(idbdx,1,nY))
            ! d2bdyy = deltaYRecip * limiter(tiles(tID)%u(idbdy,1,2)-tiles(tID)%u(idbdy,1,1), tiles(tID)%u(idbdy,1,1)-tiles(ttS)%u(idbdy,1,nY))

            ! d2bdxx = ( tiles(tID)%u(idbdx, 2,  2) - tiles( ttW)%u(idbdx, nX,  2) &
            !          + tiles(ttS)%u(idbdx, 2, nY) - tiles(ttSW)%u(idbdx, nX, nY) ) * 0.25_wp * deltaXRecip
            ! d2bdyx = ( tiles(tID)%u(idbdy, 2,  2) - tiles( ttW)%u(idbdy, nX,  2) &
            !          + tiles(ttS)%u(idbdy, 2, nY) - tiles(ttSW)%u(idbdy, nX, nY) ) * 0.25_wp * deltaXRecip
            ! d2bdxy = ( tiles(tID)%u(idbdx,  2, 2) - tiles( ttS)%u(idbdx,  2, nY) &
            !          + tiles(ttW)%u(idbdx, nX, 2) - tiles(ttSW)%u(idbdx, nX, nY) ) * 0.25_wp * deltaYRecip
            ! d2bdyy = ( tiles(tID)%u(idbdy,  2, 2) - tiles( ttS)%u(idbdy,  2, nY) &
            !          + tiles(ttW)%u(idbdy, nX, 2) - tiles(ttSW)%u(idbdy, nX, nY) ) * 0.25_wp * deltaYRecip
         
         tiles(tID)%u(id2bdxx, 1, 1) = d2bdxx
         tiles(tID)%u(id2bdxy, 1, 1) = d2bdxy !0.5_wp*(d2bdxy+d2bdyx)
         tiles(tID)%u(id2bdyy, 1, 1) = d2bdyy

         ! SE: i=nX, j=1
         d2bdxx = deltaXRecip*deltaXRecip * (tiles(ttE)%u(ib,1,1) - 2.0_wp*tiles(tID)%u(ib,nX,1) + tiles(tID)%u(ib,nX-1,1))
         d2bdxy = 0.25_wp*deltaXRecip*deltaYRecip * (tiles(ttE)%u(ib,1,2) - tiles(ttSE)%u(ib,1,nY) - tiles(tID)%u(ib,nX-1,2) + tiles(ttS)%u(ib,nX-1,nY))
         d2bdyy = deltaYRecip*deltaYRecip * (tiles(tID)%u(ib,nX,2) - 2.0_wp*tiles(tID)%u(ib,nX,1) + tiles(ttS)%u(ib,nX,nY))
         
            ! d2bdxx = deltaXRecip * limiter(tiles(ttE)%u(idbdx,1,1)-tiles(tID)%u(idbdx,nX,1), tiles(tID)%u(idbdx,nX,1)-tiles(tID)%u(idbdx,nX-1,1))
            ! d2bdyx = deltaXRecip * limiter(tiles(ttE)%u(idbdy,1,1)-tiles(tID)%u(idbdy,nX,1), tiles(tID)%u(idbdy,nX,1)-tiles(tID)%u(idbdy,nX-1,1))
            ! d2bdxy = deltaYRecip * limiter(tiles(tID)%u(idbdx,nX,2)-tiles(tID)%u(idbdx,nX,1), tiles(tID)%u(idbdx,nX,1)-tiles(ttS)%u(idbdx,nX,nY))
            ! d2bdyy = deltaYRecip * limiter(tiles(tID)%u(idbdy,nX,2)-tiles(tID)%u(idbdy,nX,1), tiles(tID)%u(idbdy,nX,1)-tiles(ttS)%u(idbdy,nX,nY-1))

            ! d2bdxx = ( tiles( ttE)%u(idbdx, 1,  2) - tiles(tID)%u(idbdx, nX-1, 2) &
            !          + tiles(ttSE)%u(idbdx, 1, nY) - tiles(ttS)%u(idbdx, nX-1, nY) ) * 0.25_wp * deltaXRecip
            ! d2bdyx = ( tiles( ttE)%u(idbdy, 1,  2) - tiles(tID)%u(idbdy, nX-1, 2) &
            !          + tiles(ttSE)%u(idbdy, 1, nY) - tiles(ttS)%u(idbdy, nX-1, nY) ) * 0.25_wp * deltaXRecip
            ! d2bdxy = ( tiles(ttE)%u(idbdx,    1, 2) - tiles(ttSE)%u(idbdx,    1, nY) &
            !          + tiles(tID)%u(idbdx, nX-1, 2) - tiles( ttS)%u(idbdx, nX-1, nY) ) * 0.25_wp * deltaYRecip
            ! d2bdyy = ( tiles(ttE)%u(idbdy,    1, 2) - tiles(ttSE)%u(idbdy,    1, nY) &
            !          + tiles(tID)%u(idbdy, nX-1, 2) - tiles( ttS)%u(idbdy, nX-1, nY) ) * 0.25_wp * deltaYRecip
         
         tiles(tID)%u(id2bdxx, nX, 1) = d2bdxx
         tiles(tID)%u(id2bdxy, nX, 1) = d2bdxy !0.5_wp*(d2bdxy+d2bdyx)
         tiles(tID)%u(id2bdyy, nX, 1) = d2bdyy

         ! NW: i=1, j=nY
         d2bdxx = deltaXRecip*deltaXRecip * (tiles(tID)%u(ib,2,nY) - 2.0_wp*tiles(tID)%u(ib,1,nY) + tiles(ttW)%u(ib,nX,nY))
         d2bdxy = 0.25_wp*deltaXRecip*deltaYRecip * (tiles(ttN)%u(ib,2,1) - tiles(tID)%u(ib,2,nY-1) - tiles(ttNW)%u(ib,nX,1) + tiles(ttW)%u(ib,nX,nY-1))
         d2bdyy = deltaYRecip*deltaYRecip * (tiles(ttN)%u(ib,1,1) - 2.0_wp*tiles(tID)%u(ib,1,nY) + tiles(tID)%u(ib,1,nY-1))
         
            ! d2bdxx = deltaXRecip * limiter(tiles(tID)%u(idbdx,2,nY)-tiles(tID)%u(idbdx,1,nY), tiles(tID)%u(idbdx,1,nY)-tiles(ttW)%u(idbdx,nX,nY))
            ! d2bdyx = deltaXRecip * limiter(tiles(tID)%u(idbdy,2,nY)-tiles(tID)%u(idbdy,2,nY), tiles(tID)%u(idbdy,1,nY)-tiles(ttW)%u(idbdy,nX,nY))
            ! d2bdxy = deltaYRecip * limiter(tiles(ttN)%u(idbdx,1,1)-tiles(tID)%u(idbdx,1,nY), tiles(tID)%u(idbdx,1,nY)-tiles(tID)%u(idbdx,1,nY-1))
            ! d2bdyy = deltaYRecip * limiter(tiles(ttN)%u(idbdy,1,1)-tiles(tID)%u(idbdy,1,nY), tiles(tID)%u(idbdy,1,nY)-tiles(tID)%u(idbdy,1,nY-1))
            
            ! d2bdxx = ( tiles(ttN)%u(idbdx, 2,    1) - tiles(ttNW)%u(idbdx, nX,    1) &
            !          + tiles(tID)%u(idbdx, 2, nY-1) - tiles( ttW)%u(idbdx, nX, nY-1) ) * 0.25_wp * deltaXRecip
            ! d2bdyx = ( tiles(ttN)%u(idbdy, 2,    1) - tiles(ttNW)%u(idbdy, nX,    1) &
            !          + tiles(tID)%u(idbdy, 2, nY-1) - tiles( ttW)%u(idbdy, nX, nY-1) ) * 0.25_wp * deltaXRecip
            ! d2bdxy = ( tiles( ttN)%u(idbdx,  2, 1) - tiles(tID)%u(idbdx,  2, nY-1) &
            !          + tiles(ttNW)%u(idbdx, nX, 1) - tiles(ttW)%u(idbdx, nX, nY-1) ) * 0.25_wp * deltaYRecip
            ! d2bdyy = ( tiles( ttN)%u(idbdy,  2, 1) - tiles(tID)%u(idbdy,  2, nY-1) &
            !          + tiles(ttNW)%u(idbdy, nX, 1) - tiles(ttW)%u(idbdy, nX, nY-1) ) * 0.25_wp * deltaYRecip
         
         tiles(tID)%u(id2bdxx, 1, nY) = d2bdxx
         tiles(tID)%u(id2bdxy, 1, nY) = d2bdxy !0.5_wp*(d2bdxy+d2bdyx)
         tiles(tID)%u(id2bdyy, 1, nY) = d2bdyy

         ! NE: i=nX, j=nY
         d2bdxx = deltaXRecip*deltaXRecip * (tiles(ttE)%u(ib,1,nY) - 2.0_wp*tiles(tID)%u(ib,nX,nY) + tiles(tID)%u(ib,nX-1,nY))
         d2bdxy = 0.25_wp*deltaXRecip*deltaYRecip * (tiles(ttNE)%u(ib,1,1) - tiles(ttE)%u(ib,1,nY-1) - tiles(ttN)%u(ib,nX-1,1) + tiles(tID)%u(ib,nX-1,nY-1))
         d2bdyy = deltaYRecip*deltaYRecip * (tiles(ttN)%u(ib,nX,1) - 2.0_wp*tiles(tID)%u(ib,nX,nY) + tiles(tID)%u(ib,nX,nY-1))
         
            ! d2bdxx = deltaXRecip * limiter(tiles(ttE)%u(idbdx,1,nY)-tiles(tID)%u(idbdx,nX,nY), tiles(tID)%u(idbdx,nX,nY)-tiles(tID)%u(idbdx,nX-1,nY))
            ! d2bdyx = deltaXRecip * limiter(tiles(ttE)%u(idbdy,1,nY)-tiles(tID)%u(idbdy,nY,nY), tiles(tID)%u(idbdy,nX,nY)-tiles(tID)%u(idbdy,nX-1,nY))
            ! d2bdxy = deltaYRecip * limiter(tiles(ttN)%u(idbdx,nX,1)-tiles(tID)%u(idbdx,nX,nY), tiles(tID)%u(idbdx,nX,nY)-tiles(tID)%u(idbdx,nX,nY-1))
            ! d2bdyy = deltaYRecip * limiter(tiles(ttN)%u(idbdy,nX,1)-tiles(tID)%u(idbdy,nX,nY), tiles(tID)%u(idbdy,nX,nY)-tiles(tID)%u(idbdy,nX,nY-1))

            ! d2bdxx = ( tiles(ttNE)%u(idbdx, 1,    1) - tiles(ttN)%u(idbdx, nX-1,    1) &
            !          + tiles( ttE)%u(idbdx, 1, nY-1) - tiles(tID)%u(idbdx, nX-1, nY-1) ) * 0.25_wp * deltaXRecip
            ! d2bdyx = ( tiles(ttNE)%u(idbdy, 1,    1) - tiles(ttN)%u(idbdy, nX-1,    1) &
            !          + tiles( ttE)%u(idbdy, 1, nY-1) - tiles(tID)%u(idbdy, nX-1, nY-1) ) * 0.25_wp * deltaXRecip
            ! d2bdxy = ( tiles(ttNE)%u(idbdx,    1, 1) - tiles(ttE)%u(idbdx,    1, nY-1) &
            !          + tiles( ttN)%u(idbdx, nX-1, 1) - tiles(tID)%u(idbdx, nX-1, nY-1) ) * 0.25_wp * deltaYRecip
            ! d2bdyy = ( tiles(ttNE)%u(idbdy,    1, 1) - tiles(ttE)%u(idbdy,    1, nY-1) &
            !          + tiles( ttN)%u(idbdy, nX-1, 1) - tiles(tID)%u(idbdy, nX-1, nY-1) ) * 0.25_wp * deltaYRecip
         
         tiles(tID)%u(id2bdxx, nX, nY) = d2bdxx
         tiles(tID)%u(id2bdxy, nX, nY) = d2bdxy !0.5_wp*(d2bdxy+d2bdyx)
         tiles(tID)%u(id2bdyy, nX, nY) = d2bdyy
         
      else

         do i=2,nX-1

            tiles(tID)%u(id2bdxx,i,1) = deltaXRecip*deltaXRecip * (tiles(tID)%u(ib,i+1,1) - 2.0_wp*tiles(tID)%u(ib,i,1) + tiles(tID)%u(ib,i-1,1))

            ! tiles(tID)%u(id2bdxx,i,1) = deltaXRecip * limiter(tiles(tID)%u(idbdx,i+1,1)-tiles(tID)%u(idbdx,i,1), tiles(tID)%u(idbdx,i,1)-tiles(tID)%u(idbdx,i-1,1))

            ! tiles(tID)%u(id2bdxx,i,1) = 0.5_wp*deltaXRecip * (tiles(tID)%u(idbdx, i+1, 1) - tiles(tID)%u(idbdx, i-1, 1))
         end do

         ! West bndy: i=1
         d2bdxx = deltaXRecip*deltaXRecip * (tiles(tID)%u(ib,2,1) - 2.0_wp*tiles(tID)%u(ib,1,1) + tiles(ttW)%u(ib,nX,1))

            ! d2bdxx = deltaXRecip * limiter(tiles(tID)%u(idbdx,2,1)-tiles(tID)%u(idbdx,1,1), tiles(tID)%u(idbdx,1,1)-tiles(ttW)%u(idbdx,nX,1))
            ! d2bdxx = ( tiles(tID)%u(idbdx, 2, 1) - tiles(ttW)%u(idbdx, nX, 1) ) * 0.5_wp * deltaXRecip
         tiles(tID)%u(id2bdxx,1,1) = d2bdxx

         ! East bndy: i=nX
         d2bdxx = deltaXRecip*deltaXRecip * (tiles(ttE)%u(ib,1,1) - 2.0_wp*tiles(tID)%u(ib,nX,1) + tiles(tID)%u(ib,nX-1,1))
            ! d2bdxx = deltaXRecip * limiter(tiles(ttE)%u(idbdx,1,1)-tiles(tID)%u(idbdx,nX,1), tiles(tID)%u(idbdx,nX,1)-tiles(tID)%u(idbdx,nX-1,1))
            ! d2bdxx = ( tiles(ttE)%u(idbdx, 1, 1) - tiles(tID)%u(idbdx, nX-1, 1) ) * 0.5_wp * deltaXRecip
         tiles(tID)%u(id2bdxx,nX,1) = d2bdxx

      end if
   end subroutine ComputeTopographicCurvatures_tileID

   subroutine ComputeTopographicCurvatures_ij(RunParams, grid, tiles, tID, i, j)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout), target :: grid
      type(tileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID, i, j

      integer :: ibt, idbdx, idbdy
      integer :: id2bdxx, id2bdyy, id2bdxy

      real(kind=wp) :: deltaXRecip, deltaYRecip
      real(kind=wp) :: d2bdxx, d2bdyy, d2bdxy

      if (.not.RunParams%curvature) return

      idbdx = RunParams%Vars%dbdx
      idbdy = RunParams%Vars%dbdy
      id2bdxx = RunParams%Vars%d2bdxx
      id2bdyy = RunParams%Vars%d2bdyy
      id2bdxy = RunParams%Vars%d2bdxy

      deltaXRecip = grid%deltaXRecip
      deltaYRecip = grid%deltaYRecip

      if (.not. RunParams%isOneD) then
         d2bdxx = deltaXRecip * ( &
            tiles(tID)%u(idbdx,i+1, j) - tiles(tID)%u(idbdx,i, j) &
         )
         tiles(tID)%u(id2bdxx, i, j) = d2bdxx

         d2bdyy = deltaYRecip * ( &
            tiles(tID)%u(idbdy,i, j+1) - tiles(tID)%u(idbdy,i, j) &
         )
         tiles(tID)%u(id2bdyy, i, j) = d2bdyy

         d2bdxy = 0.5_wp * deltaYRecip * &
            KahanSum([tiles(tID)%u(idbdx,i, j+1), -tiles(tID)%u(idbdx,i, j), &
               tiles(tID)%u(idbdx,i+1, j+1), -tiles(tID)%u(idbdx,i+1, j)])
         tiles(tID)%u(id2bdxy, i, j) = d2bdxy
      else
         d2bdxx = deltaXRecip * (tiles(tID)%u(idbdx, i+1, 1) - tiles(tID)%u(idbdx, i, 1))
         tiles(tID)%u(id2bdxx,i,1) = d2bdxx
      end if

   end subroutine ComputeTopographicCurvatures_ij


   ! For DEMS / SRTMS etc, GetHeights returns different values of b0 at the
   ! interface depending on which tile it's called from. These O(1e-10)
   ! discrepancies in interpolation can be enough to cause negative depths in
   ! the solver.
   ! Therefore, we adopt a convention of overwriting b0 in the right and top
   ! interfaces of each tile with its neighbour's values. Since every active
   ! tile is initialised with ghost tiles surrounding it, this routine should
   ! ensure that b0 and its derivatives are single-valued fields.
   subroutine EqualiseTopographicBoundaryData(RunParams, grid, tID)
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
