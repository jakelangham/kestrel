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


! This module provides utilities to work with tiles in the numerical
! simulation grid -- held in a GridType object.
!
! The grid is subdivided into tiles which can be activated as simulations
! progress, so that calculations are limited to 'active' tiles in the grid.
! Tiles in the grid are labelled by integers (origin at bottom left of grid
! and increasing in column-major order).
module update_tiles_module

   use set_precision_module, only: wp
   use messages_module, only: FatalErrorMessage
   use runsettings_module, only: RunSet
   use grid_module, only: GridCoords, GridToPhysical, GridType, OnDomainEdge, &
      TileID, TileInBounds, TileList, TileType, EastTile, WestTile, NorthTile, SouthTile
   use utilities_module, only: AddToOrderedVector, InVector, RemoveFromVector, AddToVector, WrapIndex
   use hydraulic_rhs_module, only: CalculateFluxes, CalculateLimitedDerivs, CalculateLimitedDerivsBoundary, ComputeDesingularisedVariables, &
      CorrectSlopes, Reconstruct
   use morphodynamic_rhs_module, only: ComputeCellCentredTopographicData, ComputeInterfacialTopographicData, EqualiseTopographicBoundaryData, ComputeTopographicCurvatures
   use dem_module, only: GetHeights
   use closures_module, only : Density, GeometricCorrectionFactor

   implicit none

   private
   public :: AddTile
   public :: AddTiles
   public :: AddToActiveTiles
   public :: ActivateTile
   public :: AllocateTile
   public :: ReconstructwAtEdges

contains

   ! This is essentially a wrapper routine providing a grid-agnostic
   ! interface for adding tiles to grids. 
   subroutine AddTile(grid, tile, RunParams)
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: tile
      type(RunSet), intent(in) :: RunParams
      
      if (.not. TileInBounds(grid, tile) .or. &
         (OnDomainEdge(grid, tile) .and. RunParams%bcs%s /= "periodic")) then
         if (RunParams%bcs%s == "halt") then
            call FatalErrorMessage("Error: tried to add a tile outside the domain." // &
              new_line('A') // " Try increasing the domain size or repositioning its location.")
         else ! if boundary conditions are not 'halt' this routine fails silently
            return
         end if
         return
      end if

      call AddToActiveTiles(grid, tile)
      call AllocateTile(RunParams, grid, tile)
      call ActivateTile(RunParams, grid, tile)
      call SetTileBoundaries(RunParams, grid, tile)
      call UpdateNeighbourTiles(grid, tile)
      grid%xmin = min(grid%tileContainer(tile)%xmin, grid%xmin)
      grid%xmax = max(grid%tileContainer(tile)%xmax, grid%xmax)
      grid%ymin = min(grid%tileContainer(tile)%ymin, grid%ymin)
      grid%ymax = max(grid%tileContainer(tile)%ymax, grid%ymax)
   end subroutine AddTile

   subroutine AddTiles(grid, tiles, RunParams)
      type(GridType), target, intent(inout) :: grid
      integer, dimension(:), intent(in) :: tiles
      type(RunSet), intent(in) :: RunParams

      integer :: t

      do t=1,size(tiles)
         call AddTile(grid, tiles(t), RunParams)
      end do
   end subroutine AddTiles

   ! Add tile to ActiveTiles list
   pure subroutine AddToActiveTiles(grid, k)

      implicit none

      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: k

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles, ghostTiles

      logical isGhostTile

      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles
      ghostTiles => grid%ghostTiles

      ! Is this the first active tile?
      if (.not. allocated(grid%ActiveTiles%List)) then
         allocate(grid%ActiveTiles%List(1))
         grid%ActiveTiles%List(1) = k
         grid%ActiveTiles%Size = 1
      else 
         ! Else we add only if it hasn't been added before.
         if (.not. InVector(ActiveTiles%List, k)) then
            call AddToOrderedVector(ActiveTiles%List, k)
            ActiveTiles%Size = ActiveTiles%Size + 1
         end if
            
         ! If it's a ghost tile, remove it.
         if (.not. allocated(ghostTiles%List)) return
         call RemoveFromVector(ghostTiles%List, k, isGhostTile)
         if (isGhostTile) then
            ghostTiles%size = ghostTiles%size - 1
            tileContainer(k)%isGhostTile = .false.
         end if
      end if

   end subroutine AddToActiveTiles

   ! Allocate all the data structures and set their default values
   subroutine AllocateTile(RunParams, grid, k)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: k

      integer :: d, nXpertile, nYpertile, nFlux

      type(tileType), dimension(:), pointer :: tileContainer

      tileContainer => grid%tileContainer

      tileContainer(k)%TileOn = .true.

      d = RunParams%nDimensions

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile
      nFlux = RunParams%nFlux

      if (.not. allocated(tileContainer(k)%Hnmax)) then
         allocate (tileContainer(k)%Hnmax(nXpertile, nYpertile, 2), source=0.0_wp)
        !  tileContainer(k)%Hnmax(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%umax)) then
         allocate (tileContainer(k)%umax(nXpertile, nYpertile, 2), source=0.0_wp)
        !  tileContainer(k)%umax(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%emax)) then
         allocate (tileContainer(k)%emax(nXpertile, nYpertile, 2), source=0.0_wp)
        !  tileContainer(k)%emax(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%dmax)) then
        allocate (tileContainer(k)%dmax(nXpertile, nYpertile, 2), source=0.0_wp)
        ! tileContainer(k)%dmax(:, :, :) = 0.0_wp
     end if
      if (.not. allocated(tileContainer(k)%psimax)) then
         allocate (tileContainer(k)%psimax(nXpertile, nYpertile, 2), source=0.0_wp)
        !  tileContainer(k)%psimax(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%tfirst)) then
         allocate (tileContainer(k)%tfirst(nXpertile, nYpertile), source=-1.0_wp)
        !  tileContainer(k)%tfirst(:, :) = -1.0_wp
      end if

      call AllocateU(RunParams, grid, k)

      ! Load topography
      call SetTileBoundaries(RunParams, grid, k)
      call AllocateTopographicData(RunParams, grid, k)
      call GetHeights(RunParams, grid, k)
      call ComputeCellCentredTopographicData(RunParams, grid, tileContainer, k)
      call EqualiseTopographicBoundaryData(RunParams, grid, k)

      if (.not. allocated(tileContainer(k)%hXFlux)) then
         allocate (tileContainer(k)%hXFlux(d, nXpertile + 1, nYpertile), source=0.0_wp)
        !  tileContainer(k)%hXFlux(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%hYFlux)) then
         allocate (tileContainer(k)%hYFlux(d, nXpertile, nYpertile + 1), source=0.0_wp)
        !  tileContainer(k)%hYFlux(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%gXFlux)) then
         allocate (tileContainer(k)%gXFlux(d, nXpertile + 1, nYpertile), source=0.0_wp)
        !  tileContainer(k)%gXFlux(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%gYFlux)) then
         allocate (tileContainer(k)%gYFlux(d, nXpertile, nYpertile + 1), source=0.0_wp)
        !  tileContainer(k)%gYFlux(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%pXFlux)) then
         allocate (tileContainer(k)%pXFlux(d, nXpertile + 1, nYpertile), source=0.0_wp)
        !  tileContainer(k)%pXFlux(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%pYFlux)) then
         allocate (tileContainer(k)%pYFlux(d, nXpertile, nYpertile + 1), source=0.0_wp)
        !  tileContainer(k)%pYFlux(:, :, :) = 0.0_wp
      end if

      if (.not. allocated(tileContainer(k)%ddtExplicit)) then
         allocate (tileContainer(k)%ddtExplicit(nFlux, nXpertile, nYpertile), source=0.0_wp)
        !  tileContainer(k)%ddtExplicit(:, :, :) = 0.0_wp
      end if

      if (.not. allocated(tileContainer(k)%ddtImplicit)) then
         allocate (tileContainer(k)%ddtImplicit(nFlux, nXpertile, nYpertile), source=0.0_wp)
        !  tileContainer(k)%ddtImplicit(:, :, :) = 0.0_wp
      end if

      if (RunParams%MorphodynamicsOn) then
         if (.not. allocated(tileContainer(k)%ddtExplicitBt)) then
            allocate (tileContainer(k)%ddtExplicitBt(nXpertile+1, nYpertile+1), source=0.0_wp)
            ! tileContainer(k)%ddtExplicitBt(:, :) = 0.0_wp
         end if
      end if

      tileContainer(k)%containsSource = .false.

   end subroutine AllocateTile

   pure subroutine AllocateU(RunParams, grid, k)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: k

      integer :: d, nXpertile, nYpertile

      type(tileType), dimension(:), pointer :: tileContainer

      tileContainer => grid%tileContainer
      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      d = RunParams%nDimensions

      if (.not. allocated(tileContainer(k)%u)) then
         allocate (tileContainer(k)%u(d, nXpertile, nYpertile), source=0.0_wp)
        !  tileContainer(k)%u(:, :, :) = 0.0_wp
         ! Even dry regions should have rho = 'at least' rhow since sometimes
         ! we need to divide by rho near fronts
         tileContainer(k)%u(RunParams%Vars%rho, :, :) = RunParams%rhow
      end if
      if (.not. allocated(tileContainer(k)%uLimX)) then
         allocate (tileContainer(k)%uLimX(d, nXpertile, nYpertile), source=0.0_wp)
        !  tileContainer(k)%uLimX(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%uLimY)) then
         allocate (tileContainer(k)%uLimY(d, nXpertile, nYpertile), source=0.0_wp)
        !  tileContainer(k)%uLimY(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%uPlusX)) then
         allocate (tileContainer(k)%uPlusX(d, nXpertile + 1, nYpertile), source=0.0_wp)
        !  tileContainer(k)%uPlusX(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%uMinusX)) then
         allocate (tileContainer(k)%uMinusX(d, nXpertile + 1, nYpertile), source=0.0_wp)
        !  tileContainer(k)%uMinusX(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%uPlusY)) then
         allocate (tileContainer(k)%uPlusY(d, nXpertile, nYpertile + 1), source=0.0_wp)
        !  tileContainer(k)%uPlusY(:, :, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%uMinusY)) then
         allocate (tileContainer(k)%uMinusY(d, nXpertile, nYpertile + 1), source=0.0_wp)
        !  tileContainer(k)%uMinusY(:, :, :) = 0.0_wp
      end if
      
   end subroutine AllocateU

   pure subroutine AllocateTopographicData(RunParams, grid, k)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: k

      integer :: ii, jj, nXpertile, nYpertile, grid_i, grid_j
      real(kind=wp) :: x, y

      type(tileType), dimension(:), pointer :: tileContainer

      tileContainer => grid%tileContainer

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      call GridCoords(k, grid, grid_i, grid_j)

      if (.not. allocated(tileContainer(k)%x)) then
         allocate (tileContainer(k)%x(nXpertile))
         do ii = 1, nXpertile
            call GridToPhysical(RunParams, grid, grid_i, grid_j, ii, 1, x, y)
            tileContainer(k)%x(ii) = x
         end do

         tileContainer(k)%xmin = minval(tileContainer(k)%x)
         tileContainer(k)%xmax = maxval(tileContainer(k)%x)
      end if
      if (.not. allocated(tileContainer(k)%y)) then
         allocate (tileContainer(k)%y(nYpertile))
         do jj = 1, nYpertile
            call GridToPhysical(RunParams, grid, grid_i, grid_j, 1, jj, x, y)
            tileContainer(k)%y(jj) = y
         end do

         tileContainer(k)%ymin = minval(tileContainer(k)%y)
         tileContainer(k)%ymax = maxval(tileContainer(k)%y)
      end if
      if (.not. allocated(tileContainer(k)%x_vertex)) then
         allocate (tileContainer(k)%x_vertex(nXpertile + 1))
         do ii = 1, nXpertile
            tileContainer(k)%x_vertex(ii) = &
               tileContainer(k)%x(ii) - 0.5_wp*RunParams%deltaX
         end do

         tileContainer(k)%x_vertex(nXpertile + 1) = &
            tileContainer(k)%x(nXpertile) + 0.5_wp*RunParams%deltaX
      end if
      if (.not. allocated(tileContainer(k)%y_vertex)) then
         allocate (tileContainer(k)%y_vertex(nYpertile + 1))
         do jj = 1, nYpertile
            tileContainer(k)%y_vertex(jj) = &
               tileContainer(k)%y(jj) - 0.5_wp*RunParams%deltaY
         end do

         tileContainer(k)%y_vertex(nYpertile + 1) = &
            tileContainer(k)%y(nYpertile) + 0.5_wp*RunParams%deltaY
      end if

      if (.not. allocated(tileContainer(k)%b0)) then
         allocate (tileContainer(k)%b0(nXpertile + 1, nYpertile + 1))
      end if
      if (.not. allocated(tileContainer(k)%bt)) then
         allocate (tileContainer(k)%bt(nXpertile + 1, nYpertile + 1), source=0.0_wp)
        !  tileContainer(k)%bt(:, :) = 0.0_wp
      end if
      if (.not. allocated(tileContainer(k)%EminusD)) then
         allocate (tileContainer(k)%EminusD(nXpertile, nYpertile), source=0.0_wp)
        !  tileContainer(k)%EminusD(:, :) = 0.0_wp
      end if

   end subroutine

   ! Setup all the data needed for computations - boundary derivatives etc.
   subroutine ActivateTile(RunParams, grid, k, Restart)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout), target :: grid
      integer, intent(in) :: k
      logical, intent(in), optional :: Restart

      real(kind=wp) :: unitCFLTimeStep_dummy

      type(tileType), dimension(:), pointer :: tileContainer

      logical :: initialise

      tileContainer => grid%tileContainer

      tileContainer(k)%isGhostTile = .false.

      initialise = .True.
      if (present(Restart)) then
         if (Restart) initialise = .False.
      end if

      if (initialise) then
         ! Initialise w = b, bt = 0
         grid%tileContainer(k)%u(RunParams%Vars%w, :, :) = grid%tileContainer(k)%u(RunParams%Vars%b0, :, :)
         grid%tileContainer(k)%u(RunParams%Vars%bt, :, :) = 0.0_wp
      end if

      call AddGhostTiles(RunParams, grid, k)

      if (RunParams%curvature) then
        call ComputeTopographicCurvatures(RunParams, grid, tileContainer, k)
      end if

      call ComputeInterfacialTopographicData(RunParams, grid, tileContainer, k)
      call CalculateLimitedDerivs(RunParams, grid, RunParams%iFlux, tileContainer, k)
      call Reconstruct(RunParams, grid, RunParams%iFlux, tileContainer, k)
      call CorrectSlopes(RunParams, grid, tileContainer, k)
      call ComputeDesingularisedVariables(RunParams, tileContainer, k, .true.)
      call CalculateLimitedDerivs(RunParams, grid, RunParams%iDesing, tileContainer, k)
      call Reconstruct(RunParams, grid, RunParams%iDesing, tileContainer, k)
      unitCFLTimeStep_dummy = 100000.0_wp
      call CalculateFluxes(RunParams, grid, RunParams%iFlux, tileContainer, &
                           k, unitCFLTimeStep_dummy)
   end subroutine ActivateTile

   ! Check the tiles surrounding tile k. If they are inactive, add them as
   ! ghost tiles (if they aren't already ghost tiles). We are presuming k is a
   ! newly activated tile.
   subroutine AddGhostTiles(RunParams, grid, k)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: k

      type(TileList), pointer :: ghostTiles
      type(TileType), dimension(:), pointer :: tileContainer

      integer :: i, tile, nGhosts
      integer, dimension(:), allocatable :: newGhostTiles
      integer :: d, nXpertile, nYpertile
      integer, dimension(:), allocatable :: ghostTiles_temp
      logical :: alreadyGhostTile

      tileContainer => grid%tileContainer
      ghostTiles => grid%ghostTiles

      nGhosts = 0
    !   allocate (newGhostTiles(8))
    !   newGhostTiles(:) = 0

      ! Work out which of the surrounding tiles need to be ghost tiles.
      do i = 1, grid%dim*2 + 4*(grid%dim - 1)
         ! corner tiles - needed for computing geometric info at extreme vertices
         if (i > 4) then
            if (RunParams%bcs%s == 'periodic') then
               exit
            end if
            tile = tileContainer(k)%cornertiles(i - 4)
         else
            tile = tileContainer(k)%neighbours(i)
         end if
         if (.not. TileInBounds(grid, tile)) then
            call exit
         end if

         alreadyGhostTile = (allocated(ghostTiles%List)) .and. &
                            (InVector(ghostTiles%List, tile))
         ! If tile is a already a ghost tile then it needs to fill in its w
         ! reconstruction data at the neighbouring edge
         if (alreadyGhostTile) then
            call ReconstructwAtEdges(RunParams, grid, tileContainer, tile)
            ! otherwise, if it's not a regular active tile, mark as new ghost
         else if ((.not. tileContainer(tile)%TileOn) .and. &
                  (.not. InVector(newGhostTiles, tile))) then
            nGhosts = nGhosts + 1
            call AddToVector(newGhostTiles, tile)
            ! newGhostTiles(nGhosts) = tile
         end if
      end do

      if (nGhosts == 0) then
         return ! Quick exit if no tiles to add
      end if

      ! Add them
      allocate (ghostTiles_temp(ghostTiles%size + nGhosts))
      if (ghostTiles%size > 0) then
         ghostTiles_temp(1:ghostTiles%size) = ghostTiles%List(:)
      end if
      ghostTiles_temp((ghostTiles%size + 1):(ghostTiles%size + nGhosts)) = &
         newGhostTiles(1:nGhosts)
      call move_alloc(ghostTiles_temp, ghostTiles%List)
      ghostTiles%size = ghostTiles%size + nGhosts

      ! Activate the new ghost tiles so that they are ready to use.
      d = RunParams%nDimensions
      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile
      do i = 1, nGhosts
         tile = newGhostTiles(i)

         call AllocateU(RunParams, grid, tile)
         call AllocateTopographicData(RunParams, grid, tile)

         call SetTileBoundaries(RunParams, grid, tile)
         call GetHeights(RunParams, grid, tile)
         call ComputeCellCentredTopographicData(RunParams, grid, tileContainer, tile)

         call UpdateGhostTileBCs(RunParams, grid, tile)
         call ReconstructwAtEdges(RunParams, grid, tileContainer, tile)

         grid%tileContainer(tile)%isGhostTile = .true.
      end do

      do i = 1, nGhosts
         ! this needs to 'see' the other new ghost tiles before running
         tile = newGhostTiles(i)
         call EqualiseTopographicBoundaryData(RunParams, grid, tile)
      end do

   end subroutine AddGhostTiles

   ! Set the boundary data for the tile - i.e. coordinates needed for
   ! navigating the grid and boundary condition types.
   ! N.B. For now, this is a quick paste-in of code that used to be in
   ! ActivateTile so that it can be used for ghost tiles too.
   pure subroutine SetTileBoundaries(RunParams, grid, ttk)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: ttk

      type(tileType), dimension(:), pointer :: tileContainer

      tileContainer => grid%tileContainer

      ! Populate neighbours array, ordering is W, E, N, S
      ! (W, E are first because we have to account for the 1D case)
      tileContainer(ttk)%neighbours(1) = WestTile(ttk)
      tileContainer(ttk)%neighbours(2) = EastTile(ttk)
      tileContainer(ttk)%neighbours(3) = NorthTile(ttk, grid)
      tileContainer(ttk)%neighbours(4) = SouthTile(ttk, grid)
      tileContainer(ttk)%cornertiles(1) = WestTile(SouthTile(ttk, grid))
      tileContainer(ttk)%cornertiles(2) = EastTile(SouthTile(ttk, grid))
      tileContainer(ttk)%cornertiles(3) = WestTile(NorthTile(ttk, grid))
      tileContainer(ttk)%cornertiles(4) = EastTile(NorthTile(ttk, grid))
      if (RunParams%bcs%s == "periodic") then
         call SetPeriodicBoundaryConds(grid, ttk)
      end if

      ! for when directions have to be referenced individually
      tileContainer(ttk)%West => tileContainer(ttk)%neighbours(1)
      tileContainer(ttk)%East => tileContainer(ttk)%neighbours(2)
      if (.not. RunParams%isOneD) then
         tileContainer(ttk)%North => tileContainer(ttk)%neighbours(3)
         tileContainer(ttk)%South => tileContainer(ttk)%neighbours(4)
         tileContainer(ttk)%SouthWest => tileContainer(ttk)%cornertiles(1)
         tileContainer(ttk)%SouthEast => tileContainer(ttk)%cornertiles(2)
         tileContainer(ttk)%NorthWest => tileContainer(ttk)%cornertiles(3)
         tileContainer(ttk)%NorthEast => tileContainer(ttk)%cornertiles(4)
      end if

      call UpdateNeighbourTiles(grid, ttk)

   end subroutine SetTileBoundaries

   ! For a tile in the get, workout if the neigbour tiles are active
   pure subroutine UpdateNeighbourTiles(grid, ttk)
      implicit none
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: ttk

      type(tileType), dimension(:), pointer :: tileContainer

      tileContainer => grid%tileContainer

      if (associated(tileContainer(ttk)%West)) tileContainer(ttk)%WestOn => tileContainer(tileContainer(ttk)%West)%TileOn
      if (associated(tileContainer(ttk)%East)) tileContainer(ttk)%EastOn => tileContainer(tileContainer(ttk)%East)%TileOn
      if (associated(tileContainer(ttk)%North)) tileContainer(ttk)%NorthOn => tileContainer(tileContainer(ttk)%North)%TileOn
      if (associated(tileContainer(ttk)%South)) tileContainer(ttk)%SouthOn => tileContainer(tileContainer(ttk)%South)%TileOn
   end subroutine UpdateNeighbourTiles

   ! Determine if tile is on the edge of the domain and if so, set its
   ! neighbours to wrap around.
   pure subroutine SetPeriodicBoundaryConds(grid, k)
      implicit none

      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: k

      integer :: grid_i, grid_j
      integer :: grid_i_plus, grid_i_minus
      integer :: grid_j_plus, grid_j_minus

      type(tileType), dimension(:), pointer :: tileContainer

      tileContainer => grid%tileContainer
      call GridCoords(k, grid, grid_i, grid_j)

      grid_i_plus  = WrapIndex(grid_i+1, grid%nXtiles)
      grid_i_minus = WrapIndex(grid_i-1, grid%nXtiles)
      grid_j_plus  = WrapIndex(grid_j+1, grid%nYtiles)
      grid_j_minus = WrapIndex(grid_j-1, grid%nYtiles)

      tileContainer(k)%neighbours(1) = TileID(grid, grid_i_minus, grid_j) ! W
      tileContainer(k)%neighbours(2) = TileID(grid, grid_i_plus, grid_j)  ! E
      tileContainer(k)%neighbours(3) = TileID(grid, grid_i, grid_j_plus)  ! N
      tileContainer(k)%neighbours(4) = TileID(grid, grid_i, grid_j_minus) ! S

      tileContainer(k)%cornertiles(1) = TileID(grid, grid_i_minus, grid_j_minus) ! SW
      tileContainer(k)%cornertiles(2) = TileID(grid, grid_i_plus, grid_j_minus) ! SE
      tileContainer(k)%cornertiles(3) = TileID(grid, grid_i_minus, grid_j_plus) ! NW
      tileContainer(k)%cornertiles(4) = TileID(grid, grid_i_plus, grid_j_plus) ! NE

   end subroutine SetPeriodicBoundaryConds

   ! Update the u-values for a ghost tile.
   ! N.B. This routine does not check or care whether k actually is a ghost
   ! tile.
   subroutine UpdateGhostTileBCs(RunParams, grid, k)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: k

      type(TileType), dimension(:), pointer :: tileContainer
      real(kind=wp), dimension(:, :, :), pointer :: uGhost

      tileContainer => grid%tileContainer
      uGhost => tileContainer(k)%u

      if (OnDomainEdge(grid, k)) then
         call SetDomainBoundaryData(RunParams, grid, k)
      else
         call SetDefaultTileData(RunParams, grid, k)
      end if
   end subroutine UpdateGhostTileBCs

   ! Take a tile k, assumed to be at the domain edge and write boundary data
   ! to its fields.
   ! Typical case is that this is a ghost tile needed to enforce conditions
   ! at domain edge.
   subroutine SetDomainBoundaryData(RunParams, grid, k)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: k

      integer :: i, j
      real(kind=wp), dimension(:, :, :), pointer :: u
      real(kind=wp) :: Hnval, hpval, uval, vval, psival, rho

      real(kind=wp) :: w

      u => grid%tileContainer(k)%u

      select case (RunParams%bcs%s)
       case ('periodic')
         return ! nothing to do
       case ('halt', 'sponge')
         call SetDefaultTileData(RunParams, grid, k)
       case ('dirichlet')
         Hnval = RunParams%bcsHnval
         uval = RunParams%bcsuval
         vval = RunParams%bcsvval
         psival = RunParams%bcspsival
         rho = Density(RunParams, psival)

         ! these only really need to update the outside grid pts
         do i = 1, RunParams%nXpertile
            do j = 1, RunParams%nYpertile
               hpval = Hnval/GeometricCorrectionFactor(RunParams, u(:, i, j))
            !    u(RunParams%Vars%w, i, j) = hpval + u(RunParams%Vars%b0, i, j) + u(RunParams%Vars%bt, i, j)
               w = u(RunParams%Vars%bt, i, j) + hpval
               w = w + u(RunParams%Vars%b0, i, j)
               u(RunParams%Vars%w, i, j) = w
            end do
         end do
         u(RunParams%Vars%rhoHnu, :, :) = rho*Hnval*uval
         u(RunParams%Vars%rhoHnv, :, :) = rho*Hnval*vval
         u(RunParams%Vars%Hnpsi, :, :) = Hnval*psival
         u(RunParams%Vars%Hn, :, :) = Hnval
         u(RunParams%Vars%u, :, :) = uval
         u(RunParams%Vars%v, :, :) = vval
         u(RunParams%Vars%psi, :, :) = psival
         u(RunParams%Vars%rho, :, :) = rho
       case default
         call FatalErrorMessage('Unrecognised boundary condition')
      end select
   end subroutine SetDomainBoundaryData

   ! Initialise a tile with 'default data', i.e. all vars 0
   ! apart from rho = rhow, w = b0
   subroutine SetDefaultTileData(RunParams, grid, k)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: k

      real(kind=wp), dimension(:, :, :), pointer :: u
      integer :: nVars

      u => grid%tileContainer(k)%u

      nVars = size(RunParams%iFlux) + size(RunParams%iDesing)
      u(1:nVars, :, :) = 0.0_wp

      u(RunParams%Vars%rho, :, :) = RunParams%rhow
      u(RunParams%Vars%w, :, :) = u(RunParams%Vars%b0, :, :)
   end subroutine SetDefaultTileData

   ! For the benefit of CorrectSlopes (see CCKWSolver.f90), the ghost tiles
   ! need certain w reconstructions that depend on the cells adjacent to the
   ! boundary.
   subroutine ReconstructwAtEdges(RunParams, grid, tiles, tID)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout) :: grid
      type(tileType), dimension(:), intent(inout) :: tiles
      integer, intent(in) :: tID ! tileID

      integer :: i, j, iw, neighbour, nXpertile, nYpertile
      real(kind=wp) :: deltaX, deltaY

      deltaX = grid%deltaX
      deltaY = grid%deltaY
      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      iw = RunParams%Vars%w

      neighbour = grid%tileContainer(tID)%West

      if (TileInBounds(grid, neighbour)) then
         if (tiles(neighbour)%TileOn) then
            do j = 1, nYpertile
               call CalculateLimitedDerivsBoundary(RunParams, grid, iw, &
                                                   1, j, tiles, tID, 'W')
            end do

            tiles(tID)%uPlusX(iw, 1, :) = tiles(tID)%u(iw, 1, :) - &
               tiles(tID)%uLimX(iw, 1, :) * 0.5_wp * deltaX
            tiles(tID)%uMinusX(iw, 2, :) = tiles(tID)%u(iw, 1, :) + &
               tiles(tID)%uLimX(iw, 1, :) * 0.5_wp * deltaX
         end if
      end if

      neighbour = grid%tileContainer(tID)%East
      if (TileInBounds(grid, neighbour)) then
         if (tiles(neighbour)%TileOn) then
            do j = 1, nYpertile
               call CalculateLimitedDerivsBoundary(RunParams, grid, iw, &
                                                   nXpertile, j, tiles, tID, 'E')
            end do

            tiles(tID)%uPlusX(iw, nXpertile, :) = tiles(tID)%u(iw, nXpertile, :) - &
               tiles(tID)%uLimX(iw, nXpertile, :) * 0.5_wp * deltaX
            tiles(tID)%uMinusX(iw, nXpertile+1, :) = tiles(tID)%u(iw, nXpertile, :) + &
               tiles(tID)%uLimX(iw, nXpertile, :) * 0.5_wp * deltaX
         end if
      end if

      if (.not. RunParams%isOneD) then
         neighbour = grid%tileContainer(tID)%South

         if (TileInBounds(grid, neighbour)) then
            if (tiles(neighbour)%TileOn) then
               do i = 1, nXpertile
                  call CalculateLimitedDerivsBoundary(RunParams, grid, iw, &
                                                      i, 1, tiles, tID, 'S')
               end do

               tiles(tID)%uPlusY(iw, :, 1) = tiles(tID)%u(iw, :, 1) - &
                  tiles(tID)%uLimY(iw, :, 1) * 0.5_wp * deltaY
               tiles(tID)%uMinusY(iw, :, 2) = tiles(tID)%u(iw, :, 1) + &
                  tiles(tID)%uLimY(iw, :, 1) * 0.5_wp * deltaY
            end if
         end if

         neighbour = grid%tileContainer(tID)%North
         if (TileInBounds(grid, neighbour)) then
            if (tiles(neighbour)%TileOn) then
               do i = 1, nXpertile
                  call CalculateLimitedDerivsBoundary(RunParams, grid, iw, &
                                                      i, nYpertile, tiles, tID, 'N')
               end do

               tiles(tID)%uPlusY(iw, :, nYpertile) = tiles(tID)%u(iw, :, nYpertile) - &
                  tiles(tID)%uLimY(iw, :, nYpertile) * 0.5_wp * deltaY
               tiles(tID)%uMinusY(iw, :, nYpertile+1) = tiles(tID)%u(iw, :, nYpertile) + &
                  tiles(tID)%uLimY(iw, :, nYpertile) * 0.5_wp * deltaY
            end if
         end if
      end if
   end subroutine ReconstructwAtEdges

end module update_tiles_module
