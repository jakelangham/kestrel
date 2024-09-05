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


! This module assigns values to variables in RunSet associated
! with the model domain, as read from the "Domain" block of an input file.
! Some input validation is performed, producing Warning or FatalError messages.
module domain_settings_module

   use set_precision_module, only: wp
   use messages_module, only: WarningMessage, FatalErrorMessage, InputLabelUnrecognized
   use varstring_module, only: varString
   use runsettings_module, only: RunSet

   implicit none

   private
   public :: DomainSettings_Set

! default values
   character(len=4), parameter :: bcs_d = 'halt'

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set RunParams related to Domain block in input file.
   ! Inputs: DomainLabels - Labels for Domain block variables;
   !         DomainValues - The values associated with the labels.
   ! Output: RunParams [type RunSet] - The run settings information.
   subroutine DomainSettings_Set(DomainLabels,DomainValues,RunParams)
      type(varString), dimension(:), intent(in) :: DomainLabels
      type(varString), dimension(:), intent(in) :: DomainValues
      type(RunSet), intent(inout) :: RunParams

      type(varString) :: label
      type(varString) :: bcs

      logical :: set_Lat, set_Lon
      logical :: set_nXPoints, set_nYPoints, set_regionSizeX, set_regionSizeY
      logical :: set_nXTiles, set_nYTiles
      logical :: set_nXpertile, set_nYpertile
      logical :: set_Xtilesize, set_Ytilesize
      logical :: set_boundaryConds
      logical :: set_bcs_Hnval, set_bcs_uval, set_bcs_vval, set_bcs_psival

      integer :: J, N

      N = size(DomainLabels)

      set_Lat=.FALSE.
      set_Lon=.FALSE.
      set_nXPoints=.FALSE.
      set_nYPoints=.FALSE.
      set_regionSizeX=.FALSE.
      set_regionSizeY=.FALSE.
      set_nXTiles=.FALSE.
      set_nYTiles=.FALSE.
      set_nXpertile=.FALSE.
      set_nYpertile=.FALSE.
      set_Xtilesize=.FALSE.
      set_Ytilesize=.FALSE.
      set_boundaryConds=.FALSE.
      set_bcs_Hnval=.FALSE.
      set_bcs_uval=.FALSE.
      set_bcs_vval=.FALSE.
      set_bcs_psival=.FALSE.

      do J=1,N
         label = DomainLabels(J)%to_lower()
         select case (label%s)
          case ('lat','latitude')
            set_Lat=.TRUE.
            RunParams%Lat = DomainValues(J)%to_real()
          case ('lon','longitude')
            set_Lon=.TRUE.
            RunParams%Lon = DomainValues(J)%to_real()
          case ('nxtiles')
            set_nXTiles=.TRUE.
            RunParams%nXTiles = DomainValues(J)%to_int()
          case ('nytiles')
            set_nYTiles=.TRUE.
            RunParams%nYTiles = DomainValues(J)%to_int()
          case ('nxpertile')
            set_nXpertile=.TRUE.
            RunParams%nXpertile = DomainValues(J)%to_int()
          case ('nypertile')
            set_nYpertile=.TRUE.
            RunParams%nYpertile = DomainValues(J)%to_int()
          case ('xtilesize')
            if (set_Ytilesize) then
               call WarningMessage("The 'Domain' block contains both Xtilesize and Ytilesize.  Only one should be set.  " &
                  // "Ytilesize is given first, so ignoring Xtilesize.")
            else
               set_Xtilesize=.TRUE.
               RunParams%Xtilesize = DomainValues(J)%to_real()
            end if
          case ('ytilesize')
            if (set_Xtilesize) then
               call WarningMessage("The 'Domain' block contains both Xtilesize and Ytilesize.  Only one should be set.  " &
                  // "Xtilesize is given first, so ignoring Ytilesize.")
            else
               set_Ytilesize=.TRUE.
               RunParams%Ytilesize = DomainValues(J)%to_real()
            end if
          case ('boundary conditions')
            bcs = DomainValues(J)%to_lower()
            select case (bcs%s)
             case ('halt','periodic','dirichlet')
               set_boundaryConds=.TRUE.
               RunParams%bcs = bcs
             case ('sponge')
               ! call FatalErrorMessage("Boundary condition 'sponge' is not supported yet")
               set_boundaryConds=.TRUE.
               RunParams%SpongeLayer=.TRUE.
               RunParams%bcs = bcs
             case ('reflect')
               call FatalErrorMessage("This boundary condition 'reflect' is not supported yet")
             case default
               call InputLabelUnrecognized(DomainLabels(J)%s)
            end select
          case ('boundary hn')
            set_bcs_Hnval=.TRUE.
            RunParams%bcsHnval = DomainValues(J)%to_real()
          case ('boundary u')
            set_bcs_uval=.TRUE.
            RunParams%bcsuval = DomainValues(J)%to_real()
          case ('boundary v')
            set_bcs_vval=.TRUE.
            RunParams%bcsvval = DomainValues(J)%to_real()
          case ('boundary psi')
            set_bcs_psival=.TRUE.
            RunParams%bcspsival = DomainValues(J)%to_real()
          case default
            call InputLabelUnrecognized(DomainLabels(J)%s)
         end select
      end do

      if (.not.set_Lat) call FatalErrorMessage("The 'Domain' block in input file "// trim(RunParams%InputFile%s) // new_line('A') &
         // " does not contain the required 'Lat' variable.")
      if (.not.set_Lon) call FatalErrorMessage("The 'Domain' block in input file "// trim(RunParams%InputFile%s) // new_line('A') &
         // " does not contain the required 'Lon' variable.")

      RunParams%LatLon%first = RunParams%Lat
      RunParams%LatLon%second = RunParams%Lon

      if (.not.set_nXTiles) call FatalErrorMessage("The 'Domain' block in input file "// trim(RunParams%InputFile%s) // new_line('A') &
         // " does not contain the required 'nXTiles' variable.")
      if (.not.set_nYTiles) call FatalErrorMessage("The 'Domain' block in input file "// trim(RunParams%InputFile%s) // new_line('A') &
         // " does not contain the required 'nYTiles' variable.")
      if (.not.set_nXpertile) call FatalErrorMessage("The 'Domain' block in input file "// trim(RunParams%InputFile%s) // new_line('A') &
         // " does not contain the required 'nXpertile' variable.")
      if (.not.set_nYpertile) call FatalErrorMessage("The 'Domain' block in input file "// trim(RunParams%InputFile%s) // new_line('A') &
         // " does not contain the required 'nYpertile' variable.")

      RunParams%nTiles = RunParams%nXTiles*RunParams%nYTiles

      if ((.not.set_Xtilesize).AND.(.not.set_Ytilesize)) then
         call FatalErrorMessage("The 'Domain' block in input file "// trim(RunParams%InputFile%s) // new_line('A') &
            // " does not contain either Xtilesize or Ytilesize.  One of these mustbe set.")
      else
         if (set_Xtilesize) then
            RunParams%Ytilesize = RunParams%Xtilesize*real(RunParams%nYpertile,kind=wp)/real(RunParams%nXpertile,kind=wp)
         elseif (set_Ytilesize) then
            RunParams%Xtilesize = RunParams%Ytilesize*real(RunParams%nXpertile,kind=wp)/real(RunParams%nYpertile,kind=wp)
         end if
      end if

      if (.not. set_boundaryConds) then
         RunParams%bcs = varString(bcs_d)
         call WarningMessage("'Boundary Conditions' are not given. Using boundary conditions " // RunParams%bcs%s)
      end if

      if (RunParams%bcs%s == "dirichlet") then
         if (.not. set_bcs_Hnval) then
            RunParams%bcsHnval = 0.0_wp
         end if
         if (.not. set_bcs_uval) then
            RunParams%bcsuval = 0.0_wp
         end if
         if (.not. set_bcs_vval) then
            RunParams%bcsvval = 0.0_wp
         end if
         if (.not. set_bcs_psival) then
            RunParams%bcspsival = 0.0_wp
         end if
      else if (set_bcs_Hnval .or. set_bcs_uval .or.  &
         set_bcs_vval .or. set_bcs_psival) then
         call WarningMessage("Boundary values will be ignored since not relevant" // &
            " for chosen boundary type.")
      end if

      RunParams%xSize = RunParams%nXtiles*RunParams%Xtilesize
      RunParams%ySize = RunParams%nYtiles*RunParams%Ytilesize

      RunParams%nXPoints = RunParams%nXpertile*RunParams%nXtiles
      RunParams%nYPoints = RunParams%nYpertile*RunParams%nYtiles

      if (RunParams%nYtiles*RunParams%nYpertile==1) then
         RunParams%isOneD = .TRUE.
      else
         RunParams%isOneD = .FALSE.
      end if

      RunParams%deltaX = RunParams%Xtilesize / real(RunParams%nXpertile,kind=wp)
      RunParams%deltaY = RunParams%Ytilesize / real(RunParams%nYpertile,kind=wp)
      RunParams%deltaXRecip = 1.0_wp / RunParams%deltaX
      RunParams%deltaYRecip = 1.0_wp / RunParams%deltaY
      
      return

   end subroutine DomainSettings_Set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end Module domain_settings_module
