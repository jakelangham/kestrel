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


! This module contains routines to interpolate two-dimensional data.
!
! Currently only Bicubic interpolation is implemented, with a generic interface to
! a function Bicubic that implements interpolation of integer and real valued input data.
! The output is real valued.
! A required matrix m is saved as a private module variable.

! Bicubic interpolation:
! Assume function takes the form:
!    f(x,y) = sum(i=0,...,3) sum(j=0,...,3) a_ij x^i y^j        for (x,y) in (0,1)x(0,1)
! a_ij constants determined by known f[], the f values, derivatives, and cross derivative at (0,0), (0,1), (1,0) and (1,1)
!
! Given assumed form for f(x,y), can get formulae for the components of f[] in terms of the a_ij
!      e.g. f_x(0,0) = sum(i=0,...,3) sum(j=0,...,3) a_ij i x^(i-1) y^j     at (x=0, y=0)
!                 = sum(i=0,...,3) sum(j=0,...,3) a_ij i 0^(i-1) 0^j
!                 = a_10                                         since 0^k = 0 if k>0
!
! 15 other similar equations, for f, f_x, f_y, f_xy at each of (0,0), (0,1), (1,0), (1,1)
! Thus we know f[] in terms of a matrix applied to the length-16 vector a_ij
! Invert this 16x16 matrix to get m[][], which determines the vector a_ij from the known f[]
! We don't actually know the derivative parts of f[], but can approximate these by finite differencing
module interpolate_2d_module

   use iso_c_binding
   use set_precision_module, only: wp

   implicit none

   private
   public :: Bicubic

   interface Bicubic
      module procedure Bicubic_int, Bicubic_int16, Bicubic_r
   end interface Bicubic

   real(kind=wp), dimension(:,:), parameter :: m(16,16) = reshape([ &
      1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! m(1,:)
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! m(2,:)
      -3.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! m(3,:)
      2.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! m(4,:)
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! m(5,:)
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! m(6,:)
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, & ! m(7,:)
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, & ! m(8,:)
      -3.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! m(9,:)
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, & ! m(10,:)
      9.0_wp, -9.0_wp, -9.0_wp, 9.0_wp, 6.0_wp, 3.0_wp, -6.0_wp, -3.0_wp, 6.0_wp, -6.0_wp, 3.0_wp, -3.0_wp, 4.0_wp, 2.0_wp, 2.0_wp, 1.0_wp, & ! m(11,:)
      -6.0_wp, 6.0_wp, 6.0_wp, -6.0_wp, -4.0_wp, -2.0_wp, 4.0_wp, 2.0_wp, -3.0_wp, 3.0_wp, -3.0_wp, 3.0_wp, -2.0_wp, -1.0_wp, -2.0_wp, -1.0_wp, & ! m(12,:)
      2.0_wp, -2.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! m(13,:)
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, -2.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, & ! m(14,:)
      -6.0_wp, 6.0_wp, 6.0_wp, -6.0_wp, -3.0_wp, -3.0_wp, 3.0_wp, 3.0_wp, -4.0_wp, 4.0_wp, -2.0_wp, 2.0_wp, -2.0_wp, -2.0_wp, -1.0_wp, -1.0_wp, & ! m(15,:)
      4.0_wp, -4.0_wp, -4.0_wp, 4.0_wp, 2.0_wp, 2.0_wp, -2.0_wp, -2.0_wp, 2.0_wp, -2.0_wp, 2.0_wp, -2.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp & ! m(16,:)
      ], [16,16], order=[2,1])

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Bicubic interpolation for a surface z = f(x,y) where
   ! z is real valued but passed as integer valued data.
   ! This is necessary as digital elevation maps are often
   ! provided as integer-truncated data.
   ! Data is passed as an integer two-dimensional array.
   ! The coordinates x and y are assumed to be defined on the
   ! interval 1 <= x <= Nx, 1 <= y <= Ny where Nx and Ny are the
   ! sizes of the passed data array
   ! Inputs: real :: x - x-coordinate
   !         real :: y - y-coordinate
   !         integer :: z(Nx, Ny) - integer truncated data
   ! Returns: real :: ans - interpolated result ans = z(x,y)
   pure function Bicubic_int(x,y,z) result(ans)

      implicit none

      real(kind=wp), intent(in) :: x, y
      integer, dimension(:,:), intent(in) :: z
      real(kind=wp) :: ans

      integer :: i
      integer :: x0, x1, y0, y1

      real(kind=wp) :: t, u
      real(kind=wp), dimension(16) :: fVec
      real(kind=wp), dimension(16) :: aVec
      real(kind=wp), dimension(4,4) :: a

      ! Get integer (pixel) bounds x0, x1 with x0 <= x <= x1
      x0 = floor(x)
      x1 = ceiling(x)
      ! Get integer (pixel) bounds y0, y1 with y0 <= y <= y1
      y0 = floor(y)
      y1 = ceiling(y)

      ! Retrieve values at corners
      fVec(1) = real(z(x0,y0),kind=wp) ! f(0,0)
      fVec(2) = real(z(x1,y0),kind=wp) ! f(1,0)
      fVec(3) = real(z(x0,y1),kind=wp) ! f(0,1)
      fVec(4) = real(z(x1,y1),kind=wp) ! f(1,1)
      ! x-derivative at corners
      fVec(5) = (real(z(x1,y0),kind=wp)-real(z(x0-1,y0),kind=wp))/real(x1-x0+1,kind=wp) ! fx(0,0)
      fVec(6) = (real(z(x1+1,y0),kind=wp)-real(z(x0,y0),kind=wp))/real(x1+1-x0,kind=wp) ! fx(1,0)
      fVec(7) = (real(z(x1,y1),kind=wp)-real(z(x0-1,y1),kind=wp))/real(x1-x0+1,kind=wp) ! fx(0,1)
      fVec(8) = (real(z(x1+1,y1),kind=wp)-real(z(x0,y1),kind=wp))/real(x1+1-x0,kind=wp) ! fx(1,1)
      ! y-derivative at corners
      fVec(9)  = (real(z(x0,y1),kind=wp)-real(z(x0,y0-1),kind=wp))/real(y1-y0+1,kind=wp) ! fy(0,0)
      fVec(10) = (real(z(x1,y1),kind=wp)-real(z(x1,y0-1),kind=wp))/real(y1-y0+1,kind=wp) ! fy(1,0)
      fVec(11) = (real(z(x0,y1+1),kind=wp)-real(z(x0,y0),kind=wp))/real(y1+1-y0,kind=wp) ! fy(0,1)
      fVec(12) = (real(z(x1,y1+1),kind=wp)-real(z(x1,y0),kind=wp))/real(y1+1-y0,kind=wp) ! fy(1,1)
      ! xy derivatives at corners
      fVec(13) = ((real(z(x1,y1),kind=wp)-real(z(x0-1,y1),kind=wp))-(real(z(x1,y0-1),kind=wp)-real(z(x0-1,y0-1),kind=wp)))/real(x1-x0+1,kind=wp)/real(y1-y0+1,kind=wp) ! fxy(0,0)
      fVec(14) = ((real(z(x1+1,y1),kind=wp)-real(z(x0,y1),kind=wp))-(real(z(x1+1,y0-1),kind=wp)-real(z(x0,y0-1),kind=wp)))/real(x1+1-x0,kind=wp)/real(y1-y0+1,kind=wp) ! fxy(1,0)
      fVec(15) = ((real(z(x1,y1+1),kind=wp)-real(z(x0-1,y1+1),kind=wp))-(real(z(x1,y0),kind=wp)-real(z(x0-1,y0),kind=wp)))/real(x1-x0+1,kind=wp)/real(y1+1-y0,kind=wp) ! fxy(0,1)
      fVec(16) = ((real(z(x1+1,y1+1),kind=wp)-real(z(x0,y1+1),kind=wp))-(real(z(x1+1,y0),kind=wp)-real(z(x0,y0),kind=wp)))/real(x1+1-x0,kind=wp)/real(y1+1-y0,kind=wp) ! fxy(1,1)

      aVec = matmul(m,fVec)

      a=reshape(aVec,(/4,4/),order=(/2,1/))

      t=(x-x0)/real(x1-x0,kind=wp)
      u=(y-y0)/real(y1-y0,kind=wp)

      ans=0.0_wp
      do i=4,1,-1
         ans=t*ans+((a(i,4)*u+a(i,3))*u+a(i,2))*u+a(i,1)
      end do

      return

   end function Bicubic_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   pure function Bicubic_int16(x,y,z) result(ans)

      implicit none

      real(kind=wp), intent(in) :: x, y
      integer(kind=c_int16_t), dimension(:,:), intent(in) :: z
      real(kind=wp) :: ans

      integer :: i
      integer :: x0, x1, y0, y1

      real(kind=wp) :: t, u
      real(kind=wp), dimension(16) :: fVec
      real(kind=wp), dimension(16) :: aVec
      real(kind=wp), dimension(4,4) :: a

      x0 = floor(x)
      x1 = ceiling(x)
      y0 = floor(y)
      y1 = ceiling(y)

      ! values at corners
      fVec(1) = real(z(x0,y0),kind=wp) ! f(0,0)
      fVec(2) = real(z(x1,y0),kind=wp) ! f(1,0)
      fVec(3) = real(z(x0,y1),kind=wp) ! f(0,1)
      fVec(4) = real(z(x1,y1),kind=wp) ! f(1,1)
      ! x-derivative at corners
      fVec(5) = (real(z(x1,y0),kind=wp)-real(z(x0-1,y0),kind=wp))/real(x1-x0+1,kind=wp) ! fx(0,0)
      fVec(6) = (real(z(x1+1,y0),kind=wp)-real(z(x0,y0),kind=wp))/real(x1+1-x0,kind=wp) ! fx(1,0)
      fVec(7) = (real(z(x1,y1),kind=wp)-real(z(x0-1,y1),kind=wp))/real(x1-x0+1,kind=wp) ! fx(0,1)
      fVec(8) = (real(z(x1+1,y1),kind=wp)-real(z(x0,y1),kind=wp))/real(x1+1-x0,kind=wp) ! fx(1,1)
      ! y-derivative at corners
      fVec(9)  = (real(z(x0,y1),kind=wp)-real(z(x0,y0-1),kind=wp))/real(y1-y0+1,kind=wp) ! fy(0,0)
      fVec(10) = (real(z(x1,y1),kind=wp)-real(z(x1,y0-1),kind=wp))/real(y1-y0+1,kind=wp) ! fy(1,0)
      fVec(11) = (real(z(x0,y1+1),kind=wp)-real(z(x0,y0),kind=wp))/real(y1+1-y0,kind=wp) ! fy(0,1)
      fVec(12) = (real(z(x1,y1+1),kind=wp)-real(z(x1,y0),kind=wp))/real(y1+1-y0,kind=wp) ! fy(1,1)
      ! xy derivatives at corners
      fVec(13) = ((real(z(x1,y1),kind=wp)-real(z(x0-1,y1),kind=wp))-(real(z(x1,y0-1),kind=wp)-real(z(x0-1,y0-1),kind=wp)))/real(x1-x0+1,kind=wp)/real(y1-y0+1,kind=wp) ! fxy(0,0)
      fVec(14) = ((real(z(x1+1,y1),kind=wp)-real(z(x0,y1),kind=wp))-(real(z(x1+1,y0-1),kind=wp)-real(z(x0,y0-1),kind=wp)))/real(x1+1-x0,kind=wp)/real(y1-y0+1,kind=wp) ! fxy(1,0)
      fVec(15) = ((real(z(x1,y1+1),kind=wp)-real(z(x0-1,y1+1),kind=wp))-(real(z(x1,y0),kind=wp)-real(z(x0-1,y0),kind=wp)))/real(x1-x0+1,kind=wp)/real(y1+1-y0,kind=wp) ! fxy(0,1)
      fVec(16) = ((real(z(x1+1,y1+1),kind=wp)-real(z(x0,y1+1),kind=wp))-(real(z(x1+1,y0),kind=wp)-real(z(x0,y0),kind=wp)))/real(x1+1-x0,kind=wp)/real(y1+1-y0,kind=wp) ! fxy(1,1)

      aVec = matmul(m,fVec)

      a=reshape(aVec,(/4,4/),order=(/2,1/))

      t=(x-x0)/real(x1-x0,kind=wp)
      u=(y-y0)/real(y1-y0,kind=wp)

      ans=0.0_wp
      do i=4,1,-1
         ans=t*ans+((a(i,4)*u+a(i,3))*u+a(i,2))*u+a(i,1)
      end do

      return

   end function Bicubic_int16

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function Bicubic_r(x,y,z) result(ans)

      implicit none

      real(kind=wp), intent(in) :: x, y
      real(kind=wp), dimension(:,:), intent(in) :: z
      real(kind=wp) :: ans

      integer :: i
      integer :: x0, x1, y0, y1

      real(kind=wp) :: t, u
      real(kind=wp), dimension(16) :: fVec
      real(kind=wp), dimension(16) :: aVec
      real(kind=wp), dimension(4,4) :: a

      x0 = floor(x)
      x1 = ceiling(x)
      y0 = floor(y)
      y1 = ceiling(y)

      ! values at corners
      fVec(1) = z(x0,y0) ! f(0,0)
      fVec(2) = z(x1,y0) ! f(1,0)
      fVec(3) = z(x0,y1) ! f(0,1)
      fVec(4) = z(x1,y1) ! f(1,1)
      ! x-derivative at corners
      fVec(5) = (z(x1,y0)-z(x0-1,y0))/real(x1-x0+1,kind=wp) ! fx(0,0)
      fVec(6) = (z(x1+1,y0)-z(x0,y0))/real(x1+1-x0,kind=wp) ! fx(1,0)
      fVec(7) = (z(x1,y1)-z(x0-1,y1))/real(x1-x0+1,kind=wp) ! fx(0,1)
      fVec(8) = (z(x1+1,y1)-z(x0,y1))/real(x1+1-x0,kind=wp) ! fx(1,1)
      ! y-derivative at corners
      fVec(9)  = (z(x0,y1)-z(x0,y0-1))/real(y1-y0+1,kind=wp) ! fy(0,0)
      fVec(10) = (z(x1,y1)-z(x1,y0-1))/real(y1-y0+1,kind=wp) ! fy(1,0)
      fVec(11) = (z(x0,y1+1)-z(x0,y0))/real(y1+1-y0,kind=wp) ! fy(0,1)
      fVec(12) = (z(x1,y1+1)-z(x1,y0))/real(y1+1-y0,kind=wp) ! fy(1,1)
      ! xy derivatives at corners
      fVec(13) = ((z(x1,y1)-z(x0-1,y1))-(z(x1,y0-1)-z(x0-1,y0-1)))/real(x1-x0+1,kind=wp)/real(y1-y0+1,kind=wp) ! fxy(0,0)
      fVec(14) = ((z(x1+1,y1)-z(x0,y1))-(z(x1+1,y0-1)-z(x0,y0-1)))/real(x1+1-x0,kind=wp)/real(y1-y0+1,kind=wp) ! fxy(1,0)
      fVec(15) = ((z(x1,y1+1)-z(x0-1,y1+1))-(z(x1,y0)-z(x0-1,y0)))/real(x1-x0+1,kind=wp)/real(y1+1-y0,kind=wp) ! fxy(0,1)
      fVec(16) = ((z(x1+1,y1+1)-z(x0,y1+1))-(z(x1+1,y0)-z(x0,y0)))/real(x1+1-x0,kind=wp)/real(y1+1-y0,kind=wp) ! fxy(1,1)

      aVec = 0.0_wp

      aVec = matmul(m,fVec)

      a=reshape(aVec,(/4,4/),order=(/2,1/))

      if (x==x0) then
        t = 0.0_wp
      else
         t=(x-x0)/real(x1-x0,kind=wp)
      end if
      if (y==y0) then
         u = 0.0_wp
      else
         u=(y-y0)/real(y1-y0,kind=wp)
      end if

      ans=0.0_wp
      do i=4,1,-1
         ans= t*ans &
              +( (a(i,4)*u+a(i,3))*u + a(i,2) )*u &
              + a(i,1)
      end do

      return

   end function Bicubic_r

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module interpolate_2d_module
