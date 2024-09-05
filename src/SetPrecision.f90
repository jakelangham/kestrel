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


! This main purpose of this module is to define the precision of real numbers.
! The simulation requires double precision, and to ensure interoperability with
! C/C++ libraries (GDAL, Proj) we use c_double from iso_c_bindings.
!
! Additionally, the frequently used double precision value of pi is defined.
module set_precision_module

   use, intrinsic :: iso_c_binding

   implicit none

   private
   public :: wp, pi

   integer, parameter :: wp = kind(1.0_c_double)

   real(kind=wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

end module set_precision_module
