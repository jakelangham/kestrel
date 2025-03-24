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


! This module defines the primary structure for storing data related to the
! numerical simulation grid -- GridType. It also defines associated helper
! routines for working with the grid.
!
! The grid is a rectangular domain divided into regular tiles (in 1 or 2D),
! whose number and dimension are user-defined. Each tile is further subdivided
! into a user-defined number of regular finite volumes (or cells), over which
! the governing equations are discretised. We structure the data in this way to
! give us convenient way to time step only over cell patches that actually
! contain flow, by marking tiles that contain (or are about to contain) flow as
! 'active'. Tiles have their own data structure -- TileType -- which holds
! solution data and auxiliary data for computing derivatives and fluxes. We
! record any tile adjacent to an active tile as a 'ghost' tile.  These tiles
! contain data for enforcing boundary conditions (if needed) and computing
! spatial derivatives at the edges of active tiles.
module grid_module

   use runsettings_module, only: RunSet
   use set_precision_module, only: wp
   use utilities_module, only: pair, InVector, KahanSum

   implicit none

   private
   public :: GridType
   public :: TileType, TileList
   public :: InitialiseGrid
   public :: GridCoords
   public :: TileID
   public :: TileInBounds
   public :: NorthTile, SouthTile, EastTile, WestTile
   public :: OnDomainEdge
   public :: GridToPhysical
   public :: IsActiveGhostTile

   type TileType
      logical :: TileOn ! Is the tile 'active', i.e. part of the simulation?
      logical :: isGhostTile
      logical :: containsSource

      ! Navigation data.
      ! Locations of adjacent tiles, ordered: W, E, N, S
      integer :: neighbours(4)
      logical :: activeneighbours(4)
      ! Locations of diagonal neighbours, ordered: SW, SE, NW, NE
      integer :: cornertiles(4)
      ! Pointers to the neighbours array for easy reading when they need to be
      ! referenced individually
      integer, pointer :: North, South, East, West
      integer, pointer :: NorthWest, NorthEast, SouthWest, SouthEast
      logical, pointer :: NorthOn, SouthOn, EastOn, WestOn

      ! x, y values at cell centres.
      real(kind=wp), dimension(:), allocatable :: x, y
      ! x, y values at cell vertices.
      real(kind=wp), dimension(:), allocatable :: x_vertex, y_vertex
      ! Tile extent.
      real(kind=wp) :: xmin, xmax, ymin, ymax
      ! Pair of latitude/longitude coordinates, stored if tile is georeferenced.
      type(pair), dimension(:,:), allocatable :: latlon

      ! The main solution vector, containing the values of observables at cell
      ! centres. Indexed as u(d, i, j) where d = field, i = horizontal tile
      ! coordinate, j = vertical tile coordinate.
      real(kind=wp), dimension(:,:,:), allocatable :: u 
      ! Topographic data stored on cell vertices. b0 is the original bed
      ! elevation at each vertex, measured vertically, while bt is the total
      ! change in bed elevation due to morphodynamics. (so b = b0 + bt gives the
      ! total bed height)
      real(kind=wp), dimension(:,:), allocatable :: b0
      real(kind=wp), dimension(:,:), allocatable :: bt
      ! Morphodynamic term (erosion - deposition).
      real(kind=wp), dimension(:,:), allocatable :: EminusD

      ! Storage for the maximum value of various fields over the simulation, at
      ! each cell, plus the time at which the maximum ocurred.
      real(kind=wp), dimension(:,:,:), allocatable :: Hnmax ! Max depth
      real(kind=wp), dimension(:,:,:), allocatable :: umax ! Max speed
      real(kind=wp), dimension(:,:,:), allocatable :: emax ! Max erosion
      real(kind=wp), dimension(:,:,:), allocatable :: dmax ! Max deposit
      real(kind=wp), dimension(:,:,:), allocatable :: psimax ! Max solid frac

      ! Time of first inundation (for each cell in tile).
      real(kind=wp), dimension(:,:), allocatable :: tfirst

      ! Storage for cell centred data computed during evalutaion of the
      ! hydraulic operator.
      !
      ! Limited derivatives
      real(kind=wp), dimension(:,:,:), allocatable :: uLimX, uLimY
      ! Reconstructed variables at cell edges.
      real(kind=wp), dimension(:,:,:), allocatable :: uPlusX, uMinusX
      real(kind=wp), dimension(:,:,:), allocatable :: uPlusY, uMinusY
      ! Numerical fluxes.
      real(kind=wp), dimension(:,:,:), allocatable :: hXFlux, hYFlux
      real(kind=wp), dimension(:,:,:), allocatable :: gXFlux, gYFlux
      real(kind=wp), dimension(:,:,:), allocatable :: pXFlux, pYFlux

      ! Arrays for time stepping updates (RHS data) at each substep. 
      real(kind=wp), dimension(:,:,:), allocatable :: ddtExplicit
      real(kind=wp), dimension(:,:,:), allocatable :: ddtImplicit
      real(kind=wp), dimension(:,:), allocatable :: ddtExplicitBt
   end type TileType

   type TileList
      integer, dimension(:), allocatable :: List
      integer :: size
   end type TileList

   type GridType
      ! Is this a 1D or 2D grid?
      integer :: dim
      ! Number of tiles in each direction and total, for a full grid that covers
      ! the whole domain.
      integer :: nXtiles, nYtiles, nTiles
      ! Distance between points in x and y directions (and their reciprocals).
      real(kind=wp) :: deltaX, deltaY, deltaXRecip, deltaYRecip
      ! Tile dimensions.
      real(kind=wp) :: Xtilesize, Ytilesize

      ! Associated tile data for the grid (containing the numerical solution
      ! fields).
      type(TileType), dimension(:), allocatable :: tileContainer
      type(TileList) :: activeTiles
      type(TileList) :: ghostTiles

      ! Extent of active tiles
      real(kind=wp) :: xmin, xmax, ymin, ymax

      ! Data for time stepping the grid.
      real(kind=wp) t ! Current simulation time
      real(kind=wp) t0 ! Time at beginning of substepping
      real(kind=wp) dt ! Time step

      ! These are auxiliary arrays of tiles used to store solution data at the
      ! intermediate substeps of the time stepping scheme.
      type(tileType), dimension(:), allocatable :: intermed0, intermed1
      type(tileType), dimension(:), allocatable :: intermed2, intermed3

   end type GridType

contains

   ! This should be called after ReadInputFile, in order to populate RunParams
   ! with the correct settings, but *before* we load initial/source conditions
   ! and run the simulation. Some of the default values set here may be
   ! overwritten by later routines.
   subroutine InitialiseGrid(RunParams, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(out) :: grid

      ! Grid may already exist, in which case, clear memory.
      if (allocated(grid%tileContainer)) deallocate(grid%tileContainer)
      if (allocated(grid%activeTiles%List)) deallocate(grid%activeTiles%List)

      if (RunParams%isOneD) then
         grid%dim = 1
      else
         grid%dim = 2
      end if

      grid%xmin =  huge(1.0_wp)
      grid%xmax = -huge(1.0_wp)
      grid%ymin =  huge(1.0_wp)
      grid%ymax = -huge(1.0_wp)

      grid%nTiles = RunParams%nTiles
      grid%nXtiles = RunParams%nXtiles
      grid%nYtiles = RunParams%nYtiles
      grid%deltaX = RunParams%deltaX
      grid%deltaY = RunParams%deltaY
      grid%deltaXRecip = RunParams%deltaXRecip
      grid%deltaYRecip = RunParams%deltaYRecip
      grid%Xtilesize = RunParams%Xtilesize
      grid%Ytilesize = RunParams%Ytilesize

      ! Empty grid has nTiles which are all off by default 
      allocate(grid%tileContainer(grid%nTiles))
      grid%tileContainer(:)%TileOn = .false.
      grid%activeTiles%Size = 0
      grid%ghostTiles%Size = 0

      ! Simulation time defaults to starting at t = 0.
      grid%t = 0.0_wp
      ! Default initial time step (n.b. almost always overridden by adaptive
      ! step calculated in ComputeHydraulicRHS)
      grid%dt = 1.0e-5_wp

      ! Allocate time stepping arrays.
      if (allocated(grid%intermed0)) deallocate(grid%intermed0)
      if (allocated(grid%intermed1)) deallocate(grid%intermed1)
      if (allocated(grid%intermed2)) deallocate(grid%intermed2)
      if (allocated(grid%intermed3)) deallocate(grid%intermed3)
      allocate(grid%intermed0(grid%nTiles))
      allocate(grid%intermed1(grid%nTiles))
      allocate(grid%intermed2(grid%nTiles))
      allocate(grid%intermed3(grid%nTiles))

   end subroutine InitialiseGrid

   ! Row co-ordinates of tile in grid.
   pure function TileRow(tile, grid) result(row)
      implicit none

      integer, intent(in) :: tile
      type(GridType), intent(in) :: grid
      integer :: row

      row = (tile - 1) / grid%nXtiles + 1

   end function

   ! Zero-indexed column co-ordinates of tile in grid.
   pure function TileCol(tile, grid) result(col)
      implicit none

      integer, intent(in) :: tile
      type(GridType), intent(in) :: grid
      integer :: col

      col = mod(tile - 1, grid%nXtiles) + 1

   end function

   ! Return grid coords of a given tile with respect to the grid.
   pure subroutine GridCoords(tileID, grid, i, j)
      implicit none

      integer, intent(in) :: tileID
      type(GridType), intent(in) :: grid
      integer, intent(out) :: i, j

      i = TileCol(tileID, grid)
      j = TileRow(tileID, grid)

   end subroutine GridCoords

   ! Return a unique integer for each tile in the grid.
   pure function TileID(grid, grid_i, grid_j) result(id)
      implicit none

      type(GridType), intent(in) :: grid
      integer, intent(in) :: grid_i, grid_j
      integer :: id

      id = grid_i + (grid_j - 1) * grid%nXtiles
   end function TileID

   pure function TileInBounds(grid, tileID) result(res)
      implicit none
      type(GridType), target, intent(in) :: grid
      integer, intent(in) :: tileID
      logical :: res

      res = .true.
      if (tileID < 1 .or. tileID > grid%nTiles) res = .false.
   end function TileInBounds

   pure function NorthTile(tileID, grid) result(N)
      implicit none

      integer, intent(in) :: tileID
      type(GridType), intent(in) :: grid

      integer :: N

      N = tileID + grid%nXtiles
   end function

   pure function SouthTile(tileID, grid) result(S)
      implicit none

      integer, intent(in) :: tileID
      type(GridType), intent(in) :: grid

      integer :: S

      S = tileID - grid%nXtiles
   end function

   pure function EastTile(tileID) result(E)
      implicit none

      integer, intent(in) :: tileID

      integer :: E

      E = tileID + 1
   end function

   pure function WestTile(tileID) result(W)
      implicit none

      integer, intent(in) :: tileID

      integer :: W

      W = tileID - 1
   end function

   pure function OnDomainEdge(grid, tileID) result(on)
      implicit none

      type(GridType), intent(in) :: grid
      integer, intent(in) :: tileID

      logical :: on
      integer :: i, j

      call GridCoords(tileID, grid, i, j)

      on = (i == 1 .or. i == grid%nXtiles)
      on = on .or. (grid%nYtiles > 1 .and. (j == 1 .or. j == grid%nYtiles))
   end function

   ! Return location (x, y) in physical space of the point at
   ! (grid_i, grid_j) in grid co-ords and (tile_i, tile_j) in tile co-ords.
   pure subroutine GridToPhysical(RunParams, grid, grid_i, grid_j,  &
                             tile_i, tile_j, x, y)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(in) :: grid
      integer, intent(in) :: grid_i, grid_j, tile_i, tile_j
      real(kind=wp), intent(out) :: x, y

      x = -0.5_wp * RunParams%xSize + grid%deltaX *  &
         ((grid_i - 1.0_wp) * RunParams%nXpertile + (tile_i - 0.5_wp))
      y = -0.5_wp * RunParams%ySize + grid%deltaY *  &
         ((grid_j - 1.0_wp) * RunParams%nYpertile + (tile_j - 0.5_wp))

   end subroutine GridToPhysical

   ! This routine exists to check if a tile is being used as a ghost tile. The
   ! rather convoluted procedure is neeeded because when we add the ghost tiles
   ! in a group, the List structure gets updated before all the tiles are fully
   ! ready to be used. In other situations we may get away with a simpler check.
   function IsActiveGhostTile(grid, tileID) result(active)
      implicit none

      type(GridType), target, intent(inout) :: grid
      integer, intent(in) :: tileID

      logical :: active

      active = .false.
      if (allocated(grid%ghostTiles%List)) then
         if (InVector(grid%ghostTiles%List, tileID)) then
            if (grid%tileContainer(tileID)%isGhostTile) then
               active = .true.
            end if
         end if
      end if
   end function IsActiveGhostTile

end module grid_module
