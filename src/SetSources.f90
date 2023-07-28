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


! This module creates the interface between user specified input sources
! and the numerical solver.  Flux sources, caps and cube sources are read from
! the input file in InitConds.f90 and stored in a type(RunSet)::RunParams object.
! In this module, the public LoadSourceConditions routine using RunParams and the
! numerical grid to create representation of the sources in the numerical solver.
module set_sources_module

   use set_precision_module, only: wp
   use messages_module, only: FatalErrorMessage
   use grid_module, only: GridToPhysical, GridType, TileID, TileList, TileType
   use runsettings_module, only: RunSet
   use closures_module, only : GeometricCorrectionFactor
   use update_tiles_module, only: AddTile, ReconstructwAtEdges

   implicit none

   private
   public :: LoadSourceConditions

contains

   ! Load source conditions held in RunParams into variable values
   ! on the grid.
   ! InOut: RunParams - a RunSet object containing source conditions, updated
   !        grid - a GridType object containing the numerical grid, updated
   subroutine LoadSourceConditions(RunParams, grid)
      implicit none

      type(RunSet), intent(inout) :: RunParams
      type(GridType), target, intent(inout) :: grid

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      integer :: iw, iHn, irhoHnu, irhoHnv, iHnpsi, irho, ib0, iu, ipsi
      integer :: i, j, k, n
      integer :: ii, jj, kk, tt, ttk

      integer :: nCaps, nCubes, nSrcs

      real(kind=wp) :: x, y, Hn
      real(kind=wp) :: capX, capY, capR, capHn, cappsi, capu, capv
      real(kind=wp) :: cubeX, cubeY, cubeL, cubeW, cubeHn, cubepsi
      real(kind=wp) :: srcX, srcY, srcR
      real(kind=wp) :: R2

      real(kind=wp) :: rho, rhow, rhos, rho_orig, Hn_orig, gam, hp

      ! These pointers and variables simplify notation.
      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      iw = RunParams%Vars%w
      iHn = RunParams%Vars%Hn
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      irho = RunParams%Vars%rho
      ib0 = RunParams%Vars%b0
      iu = RunParams%Vars%u
      ipsi = RunParams%Vars%psi

      rhow = RunParams%rhow
      rhos = RunParams%rhos

      nCaps = RunParams%nCaps
      nCubes = RunParams%nCubes
      nSrcs = RunParams%nSources

      if (nSrcs >= 1) then
         RunParams%FluxSources(:)%NumCellsInSrc = 0
      end if

      if (RunParams%isOneD) then
         do i = 1, RunParams%nXtiles
            k = TileID(grid, i, 1)

            do ii = 1, RunParams%nXpertile
               call GridToPhysical(RunParams, grid, i, 1, ii, 1, x, y)

               if (nCaps >= 1 .and. (RunParams%Restart .eqv. .false.)) then
                  do kk = 1, nCaps
                     capX = RunParams%CapSources(kk)%x
                     capR = RunParams%CapSources(kk)%Radius
                     capHn = RunParams%CapSources(kk)%Height
                     cappsi = RunParams%CapSources(kk)%psi
                     capu = RunParams%CapSources(kk)%u
                     call CheckSolidFractionInBounds(cappsi, RunParams)

                     rho = rhow + (rhos - rhow) * cappsi
                     Hn = capHn

                     R2 = (x-capX)*(x-capX)
                     if (R2 <= capR*capR) then
                        if (.not. tileContainer(k)%TileOn) then
                           call AddTile(grid, k, RunParams)
                        end if
                        gam = GeometricCorrectionFactor(RunParams, tileContainer(k)%u(:,ii,1))

                        ! 2nd check accounts for if tile failed to add
                        ! (e.g. if it lies outside the domain)
                        if (tileContainer(k)%TileOn) then
                           Hn_orig = tileContainer(k)%u(iHn,ii,1)
                           rho_orig = tileContainer(k)%u(irho,ii,1)
                           select case (RunParams%CapSources(kk)%Shape)
                            case ('flat')
                              tileContainer(k)%u(iw,ii,1) = tileContainer(k)%u(iw,ii,1) + capHn / gam
                              tileContainer(k)%u(iHn,ii,1) = tileContainer(k)%u(iHn,ii,1) + capHn
                              tileContainer(k)%Hnmax(ii,:,1) = tileContainer(k)%Hnmax(ii,:,1) + capHn
                              tileContainer(k)%psimax(ii,:,1) = tileContainer(k)%psimax(ii,:,1) + cappsi
                              tileContainer(k)%u(irhoHnu,ii,:) = tileContainer(k)%u(irhoHnu,ii,:) + rho*capHn*capu
                              tileContainer(k)%u(iu,ii,:) = tileContainer(k)%u(iu,ii,:) + capu
                              tileContainer(k)%u(ipsi,ii,:) = tileContainer(k)%u(ipsi,ii,:) + cappsi
                              tileContainer(k)%u(iHnpsi,ii,:) = tileContainer(k)%u(iHnpsi,ii,:) + cappsi*capHn
                            case ('para')
                              tileContainer(k)%u(iw,ii,1) = tileContainer(k)%u(iw,ii,1) + &
                                 capHn * (1.0_wp - R2 / capR / capR) / gam
                              tileContainer(k)%u(iHn,ii,1) = tileContainer(k)%u(iHn,ii,1) + &
                                 capHn * (1.0_wp - R2 / capR / capR)
                              tileContainer(k)%Hnmax(ii,:,1) = tileContainer(k)%Hnmax(ii,:,1) + &
                                 capHn * (1.0_wp - R2 / capR / capR)
                              tileContainer(k)%psimax(ii,:,1) = tileContainer(k)%psimax(ii,:,1) + cappsi
                              tileContainer(k)%u(irhoHnu,ii,:) = tileContainer(k)%u(irhoHnu,ii,:) + rho*capHn*capU
                              tileContainer(k)%u(iHnpsi,ii,:) = tileContainer(k)%u(iHnpsi,ii,:) + &
                                 cappsi*capHn*(1.0_wp-R2/capR/capR)
                            case ('level')
                              ! here, capHn actually refers to w
                              hp = capHn - tileContainer(k)%u(ib0,ii,1)
                              if (hp > 0.0_wp) then
                                 tileContainer(k)%u(iw,ii,:) = tileContainer(k)%u(iw,ii,:) + hp
                                 tileContainer(k)%Hnmax(ii,:,1) = tileContainer(k)%Hnmax(ii,:,1) + hp * gam
                                 tileContainer(k)%u(iHn,ii,:) = tileContainer(k)%u(iHn,ii,:) + hp * gam
                                 tileContainer(k)%u(irhoHnu,ii,:) = tileContainer(k)%u(irhoHnu,ii,:) + rho*hp*gam*capu
                                 tileContainer(k)%u(iHnpsi,ii,:) = tileContainer(k)%u(iHnpsi,ii,:) + cappsi*hp*gam
                                 tileContainer(k)%psimax(ii,:,1) = tileContainer(k)%psimax(ii,:,1) + cappsi
                              end if
                           end select
                           tileContainer(k)%u(irho,ii,1) = (rho_orig * Hn_orig + rho * Hn) / (Hn_orig + Hn)
                        end if

                     end if
                  end do
               end if

               if (nCubes >= 1 .and. (RunParams%Restart .eqv. .false.)) then
                  do kk = 1, nCubes
                     cubeX = RunParams%CubeSources(kk)%x
                     cubeL = RunParams%CubeSources(kk)%Length
                     cubeW = RunParams%CubeSources(kk)%Width
                     cubeHn = RunParams%CubeSources(kk)%Height
                     cubepsi = RunParams%CubeSources(kk)%psi
                     call CheckSolidFractionInBounds(cubepsi, RunParams)

                     rho = rhow + (rhos - rhow) * cubepsi
                     Hn = cubeHn

                     if (abs(x-cubeX) <= 0.5_wp*cubeL) then
                        if (.not. tileContainer(k)%TileOn) then
                           call AddTile(grid, k, RunParams)
                        end if
                        gam = GeometricCorrectionFactor(RunParams, tileContainer(k)%u(:,ii,1))

                        if (tileContainer(k)%TileOn) then
                           Hn_orig = tileContainer(k)%u(iHn,ii,1)
                           rho_orig = tileContainer(k)%u(irho,ii,1)
                           select case (RunParams%CubeSources(kk)%Shape)

                            case ('level')
                              ! in this case cubeHn actually refers to w
                              hp = cubeHn - tileContainer(k)%u(ib0,ii,1)
                              Hn = hp * gam
                              if (Hn > 0.0_wp) then
                                 tileContainer(k)%u(iw,ii,:) = tileContainer(k)%u(iw,ii,:) + &
                                    (cubeHn-tileContainer(k)%u(ib0,ii,:))
                                 tileContainer(k)%u(iHn,ii,:) = tileContainer(k)%u(iHn,ii,:) + Hn
                                 tileContainer(k)%Hnmax(ii,:,1) = tileContainer(k)%Hnmax(ii,:,1) + Hn
                                 tileContainer(k)%u(iHnpsi,ii,:) = tileContainer(k)%u(iHnpsi,ii,:) + cubepsi*Hn
                              end if
                            case ('flat')
                              tileContainer(k)%u(iw,ii,:) = tileContainer(k)%u(iw,ii,:) + cubeHn / gam
                              tileContainer(k)%u(iHn,ii,:) = tileContainer(k)%u(iHn,ii,:) + cubeHn
                              tileContainer(k)%Hnmax(ii,:,1) = tileContainer(k)%Hnmax(ii,:,1) + cubeHn
                              tileContainer(k)%u(iHnpsi,ii,:) = tileContainer(k)%u(iHnpsi,ii,:) + cubepsi*cubeHn
                           end select
                           tileContainer(k)%u(irho,ii,1) = (rho_orig * Hn_orig + rho * Hn) / (Hn_orig + Hn)
                        end if

                     end if
                  end do
               end if

               if (nSrcs >= 1) then ! flux sources
                  do kk = 1, nSrcs
                     srcX = RunParams%FluxSources(kk)%x
                     srcR = RunParams%FluxSources(kk)%Radius
                     do n = 1, RunParams%FluxSources(kk)%nFluxSeries
                        call CheckSolidFractionInBounds(RunParams%FluxSources(kk)%psi(n), RunParams)
                     end do

                     R2 = (x-srcX)*(x-srcX)
                     if (R2 <= srcR*srcR) then
                        if (.not. tileContainer(k)%TileOn) then
                           call AddTile(grid, k, RunParams)
                        end if
                        tileContainer(k)%containsSource = .true.
                        RunParams%FluxSources(kk)%NumCellsInSrc = RunParams%FluxSources(kk)%NumCellsInSrc + 1
                     end if
                  end do
               end if

            end do
         end do
      else
         do i = 1, RunParams%nXtiles
            do j = 1, RunParams%nYtiles
               k = TileID(grid, i, j)
               do ii = 1, RunParams%nXpertile
                  do jj = 1, RunParams%nYpertile

                     call GridToPhysical(RunParams, grid, i, j, ii, jj, x, y)

                     if (nCaps >= 1 .and. (RunParams%Restart .eqv. .false.)) then
                        do kk = 1, nCaps
                           capX = RunParams%CapSources(kk)%x
                           capY = RunParams%CapSources(kk)%y
                           capR = RunParams%CapSources(kk)%Radius
                           capHn = RunParams%CapSources(kk)%Height
                           cappsi = RunParams%CapSources(kk)%psi
                           capu = RunParams%CapSources(kk)%u
                           capv = RunParams%CapSources(kk)%v
                           call CheckSolidFractionInBounds(cappsi, RunParams)

                           rho = rhow + (rhos - rhow) * cappsi
                           Hn = capHn

                           R2 = (x-capX)*(x-capX)+(y-capY)*(y-capY)
                           if (R2 <= capR*capR) then
                              if (.not. tileContainer(k)%TileOn) then
                                 call AddTile(grid, k, RunParams)
                              end if
                              gam = GeometricCorrectionFactor(RunParams, tileContainer(k)%u(:,ii,jj))

                              if (tileContainer(k)%TileOn) then
                                 Hn_orig = tileContainer(k)%u(iHn,ii,jj)
                                 rho_orig = tileContainer(k)%u(irho,ii,jj)
                                 select case (RunParams%CapSources(kk)%Shape)
                                  case ('flat')
                                    tileContainer(k)%u(iw,ii,jj) = tileContainer(k)%u(iw,ii,jj) + capHn / gam
                                    tileContainer(k)%u(iHn,ii,jj) = tileContainer(k)%u(iHn,ii,jj) + capHn
                                    tileContainer(k)%Hnmax(ii,jj,1) = tileContainer(k)%Hnmax(ii,jj,1) + capHn
                                    tileContainer(k)%u(irhoHnu,ii,jj) = tileContainer(k)%u(irhoHnu,ii,jj) + rho*capHn*capu
                                    tileContainer(k)%u(irhoHnv,ii,jj) = tileContainer(k)%u(irhoHnv,ii,jj) + rho*capHn*capv
                                    tileContainer(k)%u(iHnpsi,ii,jj) = tileContainer(k)%u(iHnpsi,ii,jj) + cappsi*capHn
                                  case ('para')
                                    tileContainer(k)%u(iw,ii,jj) = tileContainer(k)%u(iw,ii,jj) + &
                                       capHn * (1.0_wp - R2 / capR / capR) / gam
                                    tileContainer(k)%u(iHn,ii,jj) = tileContainer(k)%u(iHn,ii,jj) + &
                                       capHn * (1.0_wp - R2 / capR / capR)
                                    tileContainer(k)%Hnmax(ii,jj,1) = tileContainer(k)%Hnmax(ii,jj,1) + &
                                       capHn * (1.0_wp - R2 / capR / capR)
                                    tileContainer(k)%u(irhoHnu,ii,jj) = tileContainer(k)%u(irhoHnu,ii,jj) + rho*capHn*capu
                                    tileContainer(k)%u(irhoHnv,ii,jj) = tileContainer(k)%u(irhoHnv,ii,jj) + rho*capHn*capv
                                    tileContainer(k)%u(iHnpsi,ii,jj) = tileContainer(k)%u(iHnpsi,ii,jj) + &
                                       cappsi * capHn * (1.0_wp - R2 / capR / capR)
                                  case ('level')
                                    ! in this case cubeHn actually refers to w
                                    hp = capHn - tileContainer(k)%u(ib0,ii,jj)
                                    Hn = hp * gam
                                    if (Hn > 0.0_wp) then
                                       tileContainer(k)%u(iw,ii,jj) = tileContainer(k)%u(iw,ii,jj) + hp
                                       tileContainer(k)%u(iHn,ii,jj) = tileContainer(k)%u(iHn,ii,jj) + Hn
                                       tileContainer(k)%Hnmax(ii,jj,1) = tileContainer(k)%Hnmax(ii,jj,1) + Hn
                                       tileContainer(k)%u(irhoHnu,ii,jj) = tileContainer(k)%u(irhoHnu,ii,jj) + rho*Hn*capu
                                       tileContainer(k)%u(irhoHnv,ii,jj) = tileContainer(k)%u(irhoHnv,ii,jj) + rho*Hn*capv
                                       tileContainer(k)%u(iHnpsi,ii,jj) = tileContainer(k)%u(iHnpsi,ii,jj) + cappsi*Hn
                                    end if
                                 end select
                                 tileContainer(k)%u(irho,ii,jj) = (rho_orig * Hn_orig + rho * Hn) / (Hn_orig + Hn)
                              end if
                           end if
                        end do
                     end if

                     if (nCubes >= 1 .and. (RunParams%Restart .eqv. .false.)) then ! cuboid block I.C
                        do kk = 1, nCubes
                           cubeX = RunParams%CubeSources(kk)%x
                           cubeY = RunParams%CubeSources(kk)%y
                           cubeL = RunParams%CubeSources(kk)%Length
                           cubeW = RunParams%CubeSources(kk)%Width
                           cubeHn = RunParams%CubeSources(kk)%Height
                           cubepsi = RunParams%CubeSources(kk)%psi
                           call CheckSolidFractionInBounds(cubepsi, RunParams)

                           ! n.b. formula is for user-specified volume fraction,
                           ! not concentration
                           rho = rhow + (rhos - rhow) * cubepsi
                           Hn = cubeHn

                           if ((abs(x-CubeX) <= 0.5_wp*cubeL).and.(abs(y-CubeY) <= 0.5_wp*cubeW)) then
                              if (.not. tileContainer(k)%TileOn) then
                                 call AddTile(grid, k, RunParams)
                              end if
                              gam = GeometricCorrectionFactor(RunParams, tileContainer(k)%u(:,ii,jj))

                              if (tileContainer(k)%TileOn) then
                                 Hn_orig = tileContainer(k)%u(iHn,ii,jj)
                                 rho_orig = tileContainer(k)%u(irho,ii,jj)
                                 select case (RunParams%CubeSources(kk)%Shape)

                                  case ('level')
                                    ! in this case cubeHn actually refers to w
                                    hp = cubeHn - tileContainer(k)%u(ib0,ii,jj)
                                    Hn = hp * gam
                                    if (Hn > 0.0_wp) then
                                       tileContainer(k)%u(iw,ii,jj) = tileContainer(k)%u(iw,ii,jj) + hp
                                       tileContainer(k)%u(iHn,ii,jj) = tileContainer(k)%u(iHn,ii,jj) + Hn
                                       tileContainer(k)%Hnmax(ii,jj,1) = tileContainer(k)%Hnmax(ii,jj,1) + Hn
                                       tileContainer(k)%u(iHnpsi,ii,jj) = tileContainer(k)%u(iHnpsi,ii,jj) + Hn * cubepsi
                                       tileContainer(k)%u(irhoHnu,ii,jj) = tileContainer(k)%u(irhoHnu,ii,jj) + 0.0_wp
                                       tileContainer(k)%u(irhoHnv,ii,jj) = tileContainer(k)%u(irhoHnv,ii,jj) + 0.0_wp
                                    end if
                                  case ('flat')
                                    tileContainer(k)%u(iw,ii,jj) = tileContainer(k)%u(iw,ii,jj) + cubeHn / gam
                                    tileContainer(k)%u(iHn,ii,jj) = tileContainer(k)%u(iHn,ii,jj) + cubeHn
                                    tileContainer(k)%Hnmax(ii,jj,1) = tileContainer(k)%Hnmax(ii,jj,1) + cubeHn
                                    tileContainer(k)%u(irhoHnu,ii,jj) = tileContainer(k)%u(irhoHnu,ii,jj) + 0.0_wp
                                    tileContainer(k)%u(irhoHnv,ii,jj) = tileContainer(k)%u(irhoHnv,ii,jj) + 0.0_wp
                                    tileContainer(k)%u(iHnpsi,ii,jj) = tileContainer(k)%u(iHnpsi,ii,jj) + cubepsi*cubeHn
                                 end select
                                 tileContainer(k)%u(irho,ii,jj) = (rho_orig * Hn_orig + rho * Hn) / (Hn_orig + Hn)
                              end if
                           end if
                        end do
                     end if

                     if (nSrcs >= 1) then ! flux sources
                        do kk = 1, nSrcs
                           srcX = RunParams%FluxSources(kk)%x
                           srcY = RunParams%FluxSources(kk)%y
                           srcR = RunParams%FluxSources(kk)%Radius
                           do n = 1, RunParams%FluxSources(kk)%nFluxSeries
                              call CheckSolidFractionInBounds(RunParams%FluxSources(kk)%psi(n), RunParams)
                           end do

                           R2 = (x-srcX)*(x-srcX)+(y-srcY)*(y-srcY)
                           if (R2 <= srcR*srcR) then
                              if (.not. tileContainer(k)%TileOn) then
                                 call AddTile(grid, k, RunParams)
                              end if
                              tileContainer(k)%containsSource = .true.
                              RunParams%FluxSources(kk)%NumCellsInSrc = RunParams%FluxSources(kk)%NumCellsInSrc + 1
                           end if
                        end do
                     end if

                  end do
               end do

            end do
         end do
      end if

      ! Need to update ghost tile w construction in case we just loaded w
      ! data that's directly next to a ghost tile (E.g. cap source bordering 
      ! Dirichlet boundary)
      do tt = 1, grid%ghostTiles%size
         ttk = grid%ghostTiles%List(tt)
         call ReconstructwAtEdges(RunParams, grid, grid%tileContainer, ttk)
      end do

   end subroutine LoadSourceConditions

   ! Check that user-input solid fractions are appropriate, i.e.
   ! 0 <= psi <= psi_b
   subroutine CheckSolidFractionInBounds(psi, RunParams)
      implicit none

      real(kind=wp), intent(in) :: psi
      type(RunSet), intent(in) :: RunParams

      if (psi < 0.0_wp .or. psi > (1.0_wp - RunParams%BedPorosity)) then
         call FatalErrorMessage('User-input solid fraction out of bounds')
      end if
   end subroutine CheckSolidFractionInBounds

end module set_sources_module
