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


! This module contains a collection of useful subroutines.
!
! We use the term 'vector' to refer to a 1d-array which can change size
! as elements are added or removed.
! AddToVector allows an element to be added.
! InVector checks whether a value is contained in a vector
! RemoveFromVector allows an element to be removed.
!
! An 'ordered vector' is a vector in which elements are sorted in accending order.
! AddToOrderedVector allows an element to be added to the ordered vector in the 
! correct position
!
! CheckPath determines if a varString contains a valid path
! CheckFileExists determines a filename pass as a varString exists
!
! KahanAdd and KahanSum implement Kahan compensated summation.
! KahanAdd updates a running sum using Kahan summation.
! KahanSum uses Kahan summation to add elements of a real 1d array.
!
! pair defines a type with two real numbers, stored as pair%first, pair%second.
module utilities_module
   
   use set_precision_module, only: wp
   use messages_module, only: FatalErrorMessage
   use varstring_module, only: varString

   implicit none

   private
   public :: Int2String
   public :: AddToVector, InVector, RemoveFromVector
   public :: AddToOrderedVector
   public :: PathTrail, CheckFileExists
   public :: KahanAdd, KahanSum
   public :: pair

   interface AddToVector
      module procedure AddToVector_r, AddToVector_i, AddToVector_rvec, AddToVector_ivec, AddToVector_varString, AddToVector_varString_string
   end interface AddToVector

   interface AddToOrderedVector
      module procedure AddToOrderedVector_i
   end interface AddToOrderedVector

   interface RemoveFromVector
      module procedure RemoveFromVector_i
   end interface

   interface InVector
      module procedure InVector_i, InVector_r
   end interface InVector

   interface CheckFileExists
      module procedure check_file_exists_str, check_file_exists_vstr
   end interface CheckFileExists

   type pair
      ! a two component vector
      real(kind=wp) :: first
      real(kind=wp) :: second
   end type pair

contains

   ! Convert an integer to a character string
   ! Input: n - integer to convert
   ! Returns: intString - integer as a string
   function Int2String(n) result(intString)
      integer, intent(in) :: n
      character(len=:), allocatable :: intString

      logical :: neg
      real(kind=4) :: logn

      character(len=20) :: iwdth
      character(len=20) :: nstr

      if (n==0) then
         WRITE(nstr,"(I1)") n
      elseif (n<0) then
         neg=.TRUE.
         logn=log10(-1.0*n)
         allocate(character(len=floor(logn)+1) :: intString)
         WRITE(iwdth,*) floor(logn)+1
         WRITE(nstr,"(I" // ADJUSTL(iwdth) // "." // ADJUSTL(iwdth) // ")") n
      else
         neg=.FALSE.
         logn=log10(1.0*n)
         allocate(character(len=floor(logn)) :: intString)
         WRITE(iwdth,*) floor(logn)+1
         WRITE(nstr,"(I" // ADJUSTL(iwdth) // ")") n
      end if

      intString = trim(nstr)

      return
   end function Int2String

   pure subroutine AddToVector_r(vector,val)
      real(kind=wp), dimension(:), intent(inout), allocatable :: vector
      real(kind=wp), intent(in) :: val

      real(kind=wp), dimension(:), allocatable :: temp

      integer :: N

      if (allocated(vector)) then
         N = size(vector)
         allocate(temp(N+1))
         temp(1:N) = vector
         temp(N+1) = val
         call move_alloc(temp,vector)
      else
         allocate(vector(1))
         vector(1) = val
      end if

      return

   end subroutine AddToVector_r

   pure subroutine AddToVector_i(vector,val)

      implicit none

      integer, dimension(:), intent(inout), allocatable :: vector
      integer, intent(in) :: val

      integer, dimension(:), allocatable :: temp

      integer :: N

      if (allocated(vector)) then
         N = size(vector)
         allocate(temp(N+1))
         temp(1:N) = vector
         temp(N+1) = val
         call move_alloc(temp,vector)
      else
         allocate(vector(1))
         vector(1) = val
      end if

      return

   end subroutine AddToVector_i

   pure subroutine AddToVector_varString(vector,val)

      implicit none

      type(varString), dimension(:), intent(inout), allocatable :: vector
      type(varString), intent(in) :: val

      type(varString), dimension(:), allocatable :: temp

      integer :: N

      if (allocated(vector)) then
         N = size(vector)
         allocate(temp(N+1))
         temp(1:N) = vector
         temp(N+1) = val
         call move_alloc(temp,vector)
      else
         allocate(vector(1))
         vector(1) = val
      end if

   end subroutine AddToVector_varString

   pure subroutine AddToVector_varString_string(vector,val)

      implicit none

      type(varString), dimension(:), intent(inout), allocatable :: vector
      character(*), intent(in) :: val

      type(varString), dimension(:), allocatable :: temp

      integer :: N

      if (allocated(vector)) then
         N = size(vector)
         allocate(temp(N+1))
         temp(1:N) = vector
         temp(N+1) = varString(val)
         call move_alloc(temp,vector)
      else
         allocate(vector(1))
         vector(1) = varString(val)
      end if

   end subroutine AddToVector_varString_string

   pure subroutine AddToVector_rvec(vector,val)

      implicit none

      real(kind=wp), dimension(:), intent(inout), allocatable :: vector
      real(kind=wp), dimension(:), intent(in) :: val

      real(kind=wp), dimension(:), allocatable :: temp

      integer :: Nold, Nval, Nnew

      Nval = size(val)

      if (allocated(vector)) then
         Nold = size(vector)
         Nnew = Nold+Nval

         allocate(temp(Nnew))
         temp(1:Nold) = vector
         temp(Nold+1:Nnew) = val
         call move_alloc(temp,vector)
      else
         allocate(vector(Nval))
         vector(:) = val(:)
      end if

   end subroutine AddToVector_rvec

   pure subroutine AddToVector_ivec(vector,val)

      implicit none

      integer, dimension(:), intent(inout), allocatable :: vector
      integer, dimension(:), intent(in) :: val

      integer, dimension(:), allocatable :: temp

      integer :: Nold, Nval, Nnew

      Nval = size(val)

      if (allocated(vector)) then
         Nold = size(vector)
         Nnew = Nold+Nval

         allocate(temp(Nnew))
         temp(1:Nold) = vector
         temp(Nold+1:Nnew) = val
         call move_alloc(temp,vector)
      else
         allocate(vector(Nval))
         vector(:) = val(:)
      end if

   end subroutine AddToVector_ivec

   ! Assuming numerically ordered vector, add val in correct position
   recursive subroutine AddToOrderedVector_i(vector, val)

      implicit none

      integer, dimension(:), intent(inout), allocatable :: vector
      integer, intent(in) :: val

      integer, dimension(:), allocatable :: temp

      integer :: N, i, pos

      if (.not. allocated(vector)) then
        allocate(vector(1))
        vector(1) = val
        return
      end if

      N = size(vector)

      pos = 1
      do i = 1, N
         if (val < vector(i)) then
            exit
         end if
         pos = pos + 1
      end do

      allocate(temp(N+1))
      temp(1:pos-1) = vector(1:pos-1)
      temp(pos) = val
      temp(pos+1:N+1) = vector(pos:N)

      call move_alloc(temp,vector)

   end subroutine AddToOrderedVector_i

   ! N.B.! This only removes the first occurence of val from the vector.
   ! Return value is true if something was actually removed.
   function RemoveFromVector_i(vector, val) result(removed)
      implicit none

      integer, dimension(:), allocatable, intent(inout) :: vector
      integer, intent(in) :: val

      integer, dimension(:), allocatable :: temp

      integer :: N, i, j
      logical :: removed

      N = size(vector)

      j = 0
      removed = .false.
      do i = 1, N
         if (vector(i) /= val) then
            j = j + 1
            vector(j) = vector(i)
         else
            removed = .true.
         end if
      end do

      if (removed) then
         allocate(temp(j))
         temp = vector(1:j)
         call move_alloc(temp, vector)
      end if

   end function RemoveFromVector_i

   recursive function InVector_r(vector,val) result(invec)

      implicit none

      real(kind=wp), dimension(:), intent(in) :: vector
      real(kind=wp), intent(in) :: val
      logical :: invec

      integer :: N, I

      N = size(vector)

      invec = .FALSE.
      do I = 1,N
         if (vector(I) .eq. val) then
            invec = .TRUE.
            exit
         end if
      end do

   end function InVector_r

   recursive function InVector_i(vector,val) result(invec)

      implicit none

      integer, dimension(:), allocatable, intent(in) :: vector
      integer, intent(in) :: val
      logical :: invec

      integer :: N, I

      if (.not. allocated(vector)) then
         invec = .FALSE.
         return
      end if

      N = size(vector)

      invec = .FALSE.
      do I = 1,N
         if (vector(I) .eq. val) then
            invec = .TRUE.
            exit
         end if
      end do

   end function InVector_i

   pure function PathTrail(path) result(valid_path)
      type(varString), intent(in) :: path
      type(varString) :: valid_path
      ! character(len=:), allocatable :: tmp_path

      integer :: N
      character(len=1) :: last_char

      N = len(path%s)
      last_char = path%s(N:N)

      if (last_char /= '/') then
         valid_path = path + '/'
      else
         valid_path = path
      end if

      return

   end function PathTrail

   function check_file_exists_str(fname) result(exist)
      character(len=*), intent(in) :: fname
      logical :: exist

      inquire(file=trim(fname), exist=exist)

   end function check_file_exists_str

   function check_file_exists_vstr(fname) result(exist)
      type(varString), intent(in) :: fname
      logical :: exist

      inquire(file=fname%s, exist=exist)

   end function check_file_exists_vstr

   ! Compute the sum of input array num using Kahan compensated summation
   ! @cite https://en.wikipedia.org/wiki/Kahan_summation_algorithm
   pure function KahanSum(nums) result(s)
      real(kind=wp), dimension(:), intent(in) :: nums
      real(kind=wp) :: s
      real(kind=wp) :: c

      integer :: i, N

      N = size(nums)
      s = 0.0_wp
      c = 0.0_wp

      do i=1,N
         call KahanAdd(nums(i), s, c)
      end do

   end function KahanSum

   ! Incremental Kahan summation
   ! Add value x to current sum s with compensated increment c
   ! Inputs s and c are updated and returned.  On initiation, they should both be zero.
   ! @cite https://en.wikipedia.org/wiki/Kahan_summation_algorithm
   pure subroutine KahanAdd(x, s, c)
      real(kind=wp), intent(in) :: x
      real(kind=wp), intent(inout) :: s
      real(kind=wp), intent(inout) :: c
      real(kind=wp) :: y, t ! temporary values
      y = x - c
      t = s + y
      c = (t - s) - y ! N.B. parenthesis is crucial
      s = t
   end subroutine KahanAdd

end module utilities_module
