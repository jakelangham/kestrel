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


! This module contains routines for correcting overzealous deposition, which can
! lead to negative flow depths if Hn drops too far. This is essentially an
! outcome of finite time steps being too large. Therefore, in many cases we
! prefer to refine the time step rather than make a correction, but this is not
! always possible, or feasible. In these cases, we redistribute some of the
! solids that a cell deposits to its neighbours. Cells marked for deposition are
! stored in an ordered linked list, so that the redistribution procedure occurs
! in an order that is independent of the tiling layout.
module redistribute_module

   use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
   use set_precision_module, only: wp
   use grid_module, only: GridType, TileType
   use runsettings_module, only: RunSet
   use closures_module, only: GeometricCorrectionFactor, ComputeHn
   use utilities_module, only: KahanSum
   use messages_module, only: InfoMessage, WarningMessage
   use morphodynamic_rhs_module, only: ComputeCellCentredTopographicData

   implicit none

   private
   public :: RedistList
   public :: AddToRedistList
   public :: ExcessDeposition
   public :: RedistributeGrid

   ! List data structure for storing nodes that require redistribution of
   ! the deposition. Ordered from least to most in terms of 'how much' excess
   ! deposition the cell experienced. Idea is to redistribute in this order, 
   ! which is independent of tile size.
   type RedistNode
      integer :: ttk
      integer :: i, j
      real(kind=wp) :: excess_dep
      type(RedistNode), pointer :: next => null()
   end type RedistNode

   type RedistList
      type(RedistNode), pointer :: head => null()
   end type

contains

   ! Mark the cell at tile ttk, index (i, j) for redistribution by adding it to
   ! the list structure. The input excess_dep is the depth of extra deposition
   ! that needs to be corrected for.
   subroutine AddToRedistList(list, ttk, i, j, excess_dep)
      implicit none

      type(RedistList), intent(inout) :: list
      integer, intent(in) :: ttk, i, j
      real(kind=wp), intent(in) :: excess_dep
      
      type(RedistNode), pointer :: p => null()
      type(RedistNode), pointer :: newNode => null()

      call CreateNode(newNode, ttk, i, j, excess_dep)

      if (.not. associated(list%head)) then
         list%head => newNode
      else if (list%head%excess_dep > excess_dep) then
         ! add at head node
         newNode%next => list%head
         list%head => newNode
      else
         ! else find position to add node
         p => list%head
         do while (.true.)
            if (.not. associated(p%next)) exit
            if (p%next%excess_dep > excess_dep) exit
            p => p%next
         end do

         ! add newNode after p in list
         newNode%next => p%next
         p%next => newNode
      end if

   end subroutine AddToRedistList

   pure subroutine CreateNode(node, ttk, i, j, excess_dep)
      implicit none

      type(RedistNode), pointer, intent(inout) :: node
      integer, intent(in) :: ttk, i, j
      real(kind=wp), intent(in) :: excess_dep

      allocate(node)
      node%ttk = ttk
      node%i = i
      node%j = j
      node%excess_dep = excess_dep
      node%next => null()
   end subroutine CreateNode

   pure subroutine FreeList(list)
      implicit none

      type(RedistList), intent(inout) :: list

      type(RedistNode), pointer :: p, pnext

      if (.not. associated(list%head)) return

      p => list%head
      do while (.true.)
         pnext => p%next
         deallocate(p)
         if (.not. associated(pnext)) exit
         p => pnext
      end do
   end subroutine FreeList

   ! For debugging only.
   subroutine PrintRedistList(list)
      implicit none

      type(RedistList), intent(inout) :: list
      type(RedistNode), pointer :: p

      if (.not. associated(list%head)) then
         return
      else 
         p => list%head
         do while (associated(p))
            write(stdout,*) p%ttk, p%i, p%j, p%excess_dep
            p => p%next
         end do
      end if
   end subroutine PrintRedistList

   ! Determine if a cell has deposited more solids content than in possesses.
   ! Outputs: excess_dep - amount of excess deposition
   !          deltaBt - change in bed height
   !          psiold - solids fraction prior to time step
   pure subroutine ExcessDeposition(RunParams, tilesBefore, tilesAfter, &
                               ttk, i, j, excess_dep, deltaBt, psiold)
      implicit none

      type(RunSet), intent(inout) :: RunParams
      type(TileType), dimension(:), pointer, intent(in) :: tilesBefore, tilesAfter
      integer, intent(in) :: ttk, i, j
      real(kind=wp), intent(out) :: excess_dep, deltaBt, psiold

      integer :: iw, ib0, ibt, ipsi, iHnpsi
      real(kind=wp) :: gamold, gamnew, Hn_old, Hn_new

      iw = RunParams%Vars%w
      ib0 = RunParams%Vars%b0
      ibt = RunParams%Vars%bt
      ipsi = RunParams%Vars%psi
      iHnpsi = RunParams%Vars%Hnpsi

      gamold = GeometricCorrectionFactor(RunParams, tilesBefore(ttk)%u(:,i,j))
      gamnew = GeometricCorrectionFactor(RunParams, tilesAfter(ttk)%u(:,i,j))
      Hn_old = ComputeHn(tilesBefore(ttk)%u(iw, i, j), &
         tilesBefore(ttk)%u(ib0, i, j), tilesBefore(ttk)%u(ibt, i, j), gamold)
      Hn_new = ComputeHn(tilesAfter(ttk)%u(iw, i, j), &
         tilesAfter(ttk)%u(ib0, i, j), tilesAfter(ttk)%u(ibt, i, j), gamnew)
      deltaBt = tilesAfter(ttk)%u(ibt, i, j) - tilesBefore(ttk)%u(ibt, i, j)
      psiold = tilesBefore(ttk)%u(ipsi, i, j)

      ! This is the amount we deposited that exceeds the maximum solid content.
      ! Usually we want to correct the deposition by this amount so that exactly
      ! the max is deposited.
      excess_dep = -(tilesBefore(ttk)%u(iHnpsi, i, j) * gamold /  &
                     (1.0_wp - RunParams%BedPorosity) - deltaBt)

      ! If excess_dep < -Hn_new then the cell will still have negative h if all its
      ! solids are deposited. (This is ultimately because other cells contribute
      ! to its computed height... and potentially because the cell might have
      ! greater than the max solids content.) Therefore we must make a larger
      ! correction.
      excess_dep = max(excess_dep, -(Hn_old * gamold - deltaBt))

   end subroutine ExcessDeposition

   ! Loop over all the cells in the RedistList and call RedistributeCell, which
   ! does the actual correction. We recalculate the excess deposition for each
   ! cell, since it can change if neighbouring cells are redistributed.
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
   subroutine RedistributeGrid(RunParams, grid, redistcells, refineTimeStep)
#else
   pure subroutine RedistributeGrid(RunParams, grid, redistcells, refineTimeStep)
#endif
      implicit none

      type(RunSet), intent(inout) :: RunParams
      type(GridType), target, intent(inout) :: grid
      type(RedistList), intent(inout) :: redistcells
      logical, intent(inout) :: refineTimeStep

      type(RedistNode), pointer :: cell, cellnext
      logical :: redist_success

      type(tileType), dimension(:), pointer :: intermed0, intermed3
      real(kind=wp) :: excess_dep, deltaBt, psiold

      intermed0 => grid%intermed0
      intermed3 => grid%intermed3

      cell => redistcells%head
      do while (.true.)
         call ExcessDeposition(RunParams, intermed0, intermed3, &
                               cell%ttk, cell%i, cell%j, excess_dep, deltaBt, psiold)
         if (excess_dep > epsilon(excess_dep)) then
            call RedistributeCell(RunParams, grid, &
               cell%ttk, cell%i, cell%j, excess_dep, redist_success)
            if (.not. redist_success) then
               refineTimeStep = .true.
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
               call InfoMessage('Refining at failed Redistribute')
#if DEBUG_TIMESTEP==2
               call exit
#endif
#endif
               exit
            end if
         end if
         cellnext => cell%next
         if (.not. associated(cellnext)) exit
         cell => cellnext
      end do

      call FreeList(redistcells)
   end subroutine RedistributeGrid

   ! Correct for 'excess deposition' by modifying adjacent vertices such that
   ! the morphodynamic update to the cell centre changes from bt to bt - corr
   ! (and update w and rhoHc in a commensurate way).
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
   subroutine RedistributeCell(RunParams, grid, ttk, i, j, corr, success)
#else
   pure subroutine RedistributeCell(RunParams, grid, ttk, i, j, corr, success)
#endif
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout) :: grid
      integer, intent(in) :: ttk, i, j
      real(kind=wp), intent(in) :: corr !correction
      logical, intent(out) :: success

      real(kind=wp) :: b_diff(4), sum_b_diff, delta
      real(kind=wp) :: Hnpsiold, Hnold, Hngampsi_o_psib_old, gamold, gamnew, db, w
      real(kind=wp) :: tol, adjustment, adjustment_psi, discrepancy
      integer :: ci, cj, k, l
      integer :: depi(4), depj(4)
      integer :: cttk, N, nextTile, tileSouth, tileNorth
      integer :: ib0, ibt, iw, iHn, iHnpsi

      ib0 = RunParams%Vars%b0
      ibt = RunParams%Vars%bt
      iw = RunParams%Vars%w
      iHn = RunParams%Vars%Hn
      iHnpsi = RunParams%Vars%Hnpsi

      N = 0
      sum_b_diff = 0.0_wp
      b_diff(:) = 0.0_wp

      ! Scan through vertices of cell (i, j) checking whether each vertex
      ! was net depositional in the most recent time step and marking it
      ! for correction if so.
      b_diff(N + 1) = grid%intermed3(ttk)%bt(i, j) - grid%intermed0(ttk)%bt(i, j)
      if (b_diff(N + 1) > 0.0_wp) then
         N = N + 1
         depi(N) = i; depj(N) = j
         sum_b_diff = sum_b_diff + b_diff(N)
      end if
      b_diff(N + 1) = grid%intermed3(ttk)%bt(i+1, j) - grid%intermed0(ttk)%bt(i+1, j)
      if (b_diff(N + 1) > 0.0_wp) then
         N = N + 1
         depi(N) = i + 1; depj(N) = j
         sum_b_diff = sum_b_diff + b_diff(N)
      end if
      b_diff(N + 1) = grid%intermed3(ttk)%bt(i, j+1) - grid%intermed0(ttk)%bt(i, j+1)
      if (b_diff(N + 1) > 0.0_wp) then
         N = N + 1
         depi(N) = i; depj(N) = j + 1
         sum_b_diff = sum_b_diff + b_diff(N)
      end if
      b_diff(N + 1) = grid%intermed3(ttk)%bt(i+1, j+1) - grid%intermed0(ttk)%bt(i+1, j+1)
      if (b_diff(N + 1) > 0.0_wp) then
         N = N + 1
         depi(N) = i + 1; depj(N) = j + 1
         sum_b_diff = sum_b_diff + b_diff(N)
      end if

      if (N == 0) then
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
         call WarningMessage("Trying to correct excess deposition within a cell " &
            // "that borders no depositional vertices. " &
            // "(This should never happen.)")
         call WarningMessage("Excess deposition = ", corr, fmt='(a,2P,G0)')
#endif

        !  call exit
         success = .false.
         return
      end if

      ! Correction needed for each depositional vertex is
      ! db -> db * (1 - delta)
      delta = 4.0_wp * corr / sum_b_diff

      ! However, if correction is going to bring Hn very close to zero, the
      ! limitations of finite precision arithmetic mean we might nevertheless
      ! end up with a negative depth. We err on the side of caution and over
      ! correct in this case.
      Hnold = grid%intermed0(ttk)%u(iHn, i, j)
      gamold = GeometricCorrectionFactor(RunParams, grid%intermed0(ttk)%u(:, i, j))
      db = grid%intermed3(ttk)%u(ibt,i,j) - grid%intermed0(ttk)%u(ibt,i,j)

      ! Aim is to make sure corrections bring depths up to something commensurate
      ! with this tol, which tries to (over)estimate error in computing h = w - b
      tol = epsilon(grid%intermed3(ttk)%u(iw, i, j)) * grid%intermed3(ttk)%u(iw, i, j) * 10.0_wp

      adjustment = 0.0_wp
      discrepancy = KahanSum([Hnold * gamold, -db, corr])
      if (abs(discrepancy) < tol) then 
         ! This adjustment keeps delta < 1.
         adjustment = KahanSum([tol, -Hnold * gamold, db])
         adjustment = adjustment * (4.0_wp / sum_b_diff)
         adjustment = adjustment - delta
         adjustment = max(adjustment, 0.0_wp)
      end if

      ! Also have to do an analogous check for Hnpsi positivity.
      Hngampsi_o_psib_old = grid%intermed0(ttk)%u(iHnpsi, i, j) * gamold / (1.0_wp - RunParams%BedPorosity)
      tol = epsilon(Hngampsi_o_psib_old) * 10.0_wp
      discrepancy = KahanSum([Hngampsi_o_psib_old, -db, corr])
      if (abs(discrepancy) < tol) then 
         adjustment_psi = KahanSum([tol, -Hngampsi_o_psib_old, db])
         adjustment_psi = adjustment_psi * (4.0_wp / sum_b_diff)
         adjustment_psi = adjustment_psi - delta
         adjustment = max(adjustment, adjustment_psi)
      end if

      ! Correction adjusted accordingly (if need be).
      delta = delta + adjustment

      ! Apply corrections to bt at vertices.
      do k = 1, N
         grid%intermed3(ttk)%bt(depi(k), depj(k)) = &
            grid%intermed3(ttk)%bt(depi(k), depj(k)) - delta * b_diff(k)

         ! If at a tile boundary, the bt update must match the vertex on the
         ! adjacent tile.
         if (depi(k) == 1) then
            nextTile = grid%intermed3(ttk)%West
            grid%intermed3(nextTile)%bt(RunParams%nXpertile + 1, depj(k)) = grid%intermed3(ttk)%bt(depi(k), depj(k))
         else if (depi(k) == RunParams%nXpertile + 1) then
            nextTile = grid%intermed3(ttk)%East
            grid%intermed3(nextTile)%bt(1, depj(k)) = grid%intermed3(ttk)%bt(depi(k), depj(k))
         end if
         if (.not. RunParams%isOneD) then ! upper/lower bdries
            if (depj(k) == 1) then
               nextTile = grid%intermed3(ttk)%South
               grid%intermed3(nextTile)%bt(depi(k), RunParams%nYpertile + 1) = grid%intermed3(ttk)%bt(depi(k), depj(k))
               if (depi(k) == 1) then ! bottom left corner
                  tileSouth = grid%intermed3(ttk)%South
                  nextTile = grid%intermed3(tileSouth)%West
                  grid%intermed3(nextTile)%bt(RunParams%nXpertile + 1, RunParams%nYpertile + 1) = grid%intermed3(ttk)%bt(depi(k), depj(k))
               end if
               if (depi(k) == RunParams%nXpertile + 1) then ! bottom right corner
                  tileSouth = grid%intermed3(ttk)%South
                  nextTile = grid%intermed3(tileSouth)%East
                  grid%intermed3(nextTile)%bt(1, RunParams%nYpertile + 1) = grid%intermed3(ttk)%bt(depi(k), depj(k))
               end if
            end if
            if (depj(k) == RunParams%nYpertile + 1) then
               nextTile = grid%intermed3(ttk)%North
               grid%intermed3(nextTile)%bt(depi(k), 1) = grid%intermed3(ttk)%bt(depi(k), depj(k))
               if (depi(k) == 1) then ! top left corner
                  tileNorth = grid%intermed3(ttk)%North
                  nextTile = grid%intermed3(tileNorth)%West
                  grid%intermed3(nextTile)%bt(RunParams%nXpertile + 1, 1) = grid%intermed3(ttk)%bt(depi(k), depj(k))
               end if
               if (depi(k) == RunParams%nXpertile + 1) then ! top right corner
                  tileNorth = grid%intermed3(ttk)%North
                  nextTile = grid%intermed3(tileNorth)%East
                  grid%intermed3(nextTile)%bt(1, 1) = grid%intermed3(ttk)%bt(depi(k), depj(k))
               end if
            end if
         end if
      end do

      ! Scan through cell centres, updating the interpolated bt value.
      ! N.B. Up to 3^D cells may be affected, so we have to be careful about
      ! boundaries.
      if (i == 1) then
         cttk = grid%tileContainer(ttk)%West
         ci = RunParams%nXpertile
      else
         cttk = ttk
         ci = i - 1
      end if
      if (.not. RunParams%isOneD) then
         do k = -1, 1
            if (j == 1) then
               cttk = grid%tileContainer(cttk)%South
               cj = RunParams%nYpertile
            else
               cj = j - 1
            end if
            do l = -1, 1
               call ComputeCellCentredTopographicData(RunParams, grid, grid%intermed3, cttk, ci, cj)
               ! Update mass variables linearly dependent on bt.
               db = grid%intermed3(cttk)%u(ibt, ci, cj) - grid%intermed0(cttk)%u(ibt, ci, cj)
               gamold = GeometricCorrectionFactor(RunParams, grid%intermed0(cttk)%u(:, ci, cj))
               gamnew = GeometricCorrectionFactor(RunParams, grid%intermed3(cttk)%u(:, ci, cj))
               Hnold = grid%intermed0(cttk)%u(iHn, ci, cj)

               w = grid%intermed3(cttk)%u(ibt, ci, cj)
               w = w + (Hnold * gamold / gamnew - db / gamnew) / gamnew
               w = w + grid%intermed3(cttk)%u(ib0, ci, cj)
               grid%intermed3(cttk)%u(iw, ci, cj) = w
               Hnpsiold = grid%intermed0(cttk)%u(iHnpsi, ci, cj)
               grid%intermed3(cttk)%u(iHnpsi, ci, cj) = Hnpsiold * gamold / gamnew - (1.0_wp - RunParams%BedPorosity) * db / gamnew

               cj = cj + 1
               if (cj == RunParams%nYpertile + 1 .and. l < 1) then
                  cttk = grid%tileContainer(cttk)%North
                  cj = 1
               end if
            end do
            if (j + 1 > RunParams%nYpertile) then
               cttk = grid%tileContainer(cttk)%South ! move back down
            end if
            ci = ci + 1
            if (ci == RunParams%nXpertile + 1) then
               cttk = grid%tileContainer(cttk)%East
               ci = 1
            end if
         end do
      else ! 1D part
         do k = -1, 1
            call ComputeCellCentredTopographicData(RunParams, grid, grid%intermed3, cttk, ci, 1)

            ! Update mass variables linearly dependent on bt.
            db = grid%intermed3(cttk)%u(ibt, ci, 1) - grid%intermed0(cttk)%u(ibt, ci, 1)
            gamold = GeometricCorrectionFactor(RunParams, grid%intermed0(cttk)%u(:, ci, 1))
            gamnew = GeometricCorrectionFactor(RunParams, grid%intermed3(cttk)%u(:, ci, 1))
            Hnold = grid%intermed0(cttk)%u(iHn, ci, 1)
            w = - db / gamnew / gamnew
            w = w + Hnold * gamold / gamnew / gamnew
            w = w + grid%intermed3(cttk)%u(ibt, ci, 1)
            w = w + grid%intermed3(cttk)%u(ib0, ci, 1)
            grid%intermed3(cttk)%u(iw, ci, 1) = w
            Hnpsiold = grid%intermed0(cttk)%u(iHnpsi, ci, 1)
            grid%intermed3(cttk)%u(iHnpsi, ci, 1) = Hnpsiold * gamold / gamnew - (1.0_wp - RunParams%BedPorosity) * db / gamnew

            ci = ci + 1
            if (ci == RunParams%nXpertile + 1) then
               cttk = grid%tileContainer(cttk)%East
               ci = 1
            end if
         end do
      end if
      success = .true.

   end subroutine RedistributeCell

end module redistribute_module
