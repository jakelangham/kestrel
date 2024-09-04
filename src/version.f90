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

module version_module

   use, intrinsic :: iso_c_binding, only: c_int, c_char

   interface
      subroutine get_version_tag(str, length) bind(C, name="get_version_tag_c")
         import :: c_int, c_char
         implicit none
         character(kind=c_char), dimension(*), intent(out) :: str
         integer(c_int), intent(out) :: length
      end subroutine get_version_tag
   end interface

contains

   function GetVersion() result(tag)
      implicit none
      character(len=:), allocatable :: tag
      character(kind=c_char), dimension(:), allocatable :: str_buffer
      integer(c_int) :: length

      ! Allocate a buffer large enough to hold the C string
      allocate(str_buffer(100))
      
      call get_version_tag(str_buffer, length)
      if (length > 0) then
         allocate(character(len=length) :: tag)
         tag = transfer(str_buffer(1:length), tag)
      else
         tag = "Unknown"
      end if

      deallocate(str_buffer)
      
   end function GetVersion

end module version_module
