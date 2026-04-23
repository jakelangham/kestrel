! This file is part of the Kestrel software for simulations
! of sediment-laden Earth surface flows.
!
! Version v1.1.2
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


! This module defines a class for a poor-mans dictionary (Dict) of key=value pairs.
! The Dict has allocatable 1D-arrays for the keys and values.
! There is no efficient searching, only looping through the array, so the Dict is not
! intended as a data-structure for repeated use.  Instead, it is a convenient for temporary
! storage for initialization.
! The keys and values are held as varStrings (see varStringClass.f90), but values can be set and
! retrieved as integer, real, logical, character or varString.
!
! The Dict is empty on instantiation, and data is added using the overloaded append subroutine, e.g.
!     type(Dict) :: d
!     call d%append('a', 1)
!     call d%append('b', 2.0_wp)
!     call d%append('c', .TRUE.)
!     call d%append('d', 'four')
!
! The length of the Dict is the number of keys (which must equal the number of values), and is found
! using the function len, e.g.
!     type(Dict) :: d
!     call d%append('a', 1)
!     call d%append('b', 2.0_wp)
!     print *, d%len()  ! return 2
! If no keys/values are set, len returns 0.
!
! Values are retrieved from the Dict using the overloaded get subroutine, e.g.
!     type(Dict) :: d
!     integer :: a_value
!     call d%append('a', 1)
!     call d%get('a', a_value)
!     print *, a  ! return 1
! Note, there is an attempt to cast the varString value to the data type of the return value, e.g.
!     type(Dict) :: d
!     real(kind=wp) :: a_value
!     call d%append('a', 1)
!     call d%get('a', a_value)
!     print *, a_value  ! return 1.0
!
! The append method skips if the key already exists, e.g.
!     type(Dict) :: d
!     integer :: a_value
!     call d%append('a', 1)
!     call d%append('a', 2)
!     call d%get('a', a_value)
!     print *, a_value  ! return 1
!
! The existence of a key in the Dict can be queried using the has_key function, e.g.
!     type(Dict) :: d
!     logical :: key_in_d
!     call d%append('a', 1)
!     key_in_d = d%has_key('a')  ! return T
!     key_in_d = d%has_key('b')  ! return F
!
! Dicts can be combined using addition, e.g.
!     type(Dict) :: d1, d2, d
!     call d1%append('a', 1)
!     call d2%append('b', 2)
!     d = d1 + d2    ! d has keys ('a', 'b') with values (1, 2)
! Note, if keys of d1 have precedence over the keys of d2.

module dict_module

   use set_precision_module, only : wp
   use messages_module, only: FatalErrorMessage, WarningMessage
   use varstring_module, only: varString
   use utilities_module, only: AddToVector, InVector, IndexInVector

   implicit none
   private
   public :: Dict

   type :: Dict
      type(varString), dimension(:), allocatable :: keys
      type(varString), dimension(:), allocatable :: values
   contains
      procedure, pass(this) :: append_str_str, append_vstr_vstr, append_str_vstr, append_vstr_str
      procedure, pass(this) :: append_str_int, append_vstr_int
      procedure, pass(this) :: append_str_real, append_vstr_real
      generic, public :: append => append_str_str, append_vstr_str, &
                                append_str_vstr, append_vstr_vstr, &
                                append_str_int, append_vstr_int, &
                                append_str_real, append_vstr_real

      procedure, pass(this) :: get_real_str, get_real_vstr
      procedure, pass(this) :: get_int_str, get_int_vstr
      procedure, pass(this) :: get_logical_str, get_logical_vstr
      procedure, pass(this) :: get_vstring_str, get_vstring_vstr
      generic, public :: get => get_vstring_str, get_vstring_vstr, &
         get_real_str, get_real_vstr, &
         get_int_str, get_int_vstr, &
         get_logical_str, get_logical_vstr

      procedure, pass(this) :: get_real_str_default, get_real_vstr_default
      procedure, pass(this) :: get_int_str_default, get_int_vstr_default
      procedure, pass(this) :: get_logical_str_default, get_logical_vstr_default
      procedure, pass(this) :: get_vstring_str_default, get_vstring_vstr_default
      generic, public :: get_or_default => get_vstring_str_default, get_vstring_vstr_default, &
         get_real_str_default, get_real_vstr_default, &
         get_int_str_default, get_int_vstr_default, &
         get_logical_str_default, get_logical_vstr_default

      procedure, pass(this) :: has_key_str, has_key_vstr
      generic, public :: has_key => has_key_str, has_key_vstr

      procedure, pass(this) :: dict_length
      generic, public :: len => dict_length

      procedure, pass(this) :: combine
      generic, public :: operator(+) => combine

      final :: destructor
   end type Dict

   interface write(unformatted)
      module procedure dict_write_ufmt
   end interface

contains

   subroutine append_str_str(this, key, value, warn)
      class(Dict), intent(inout) :: this
      character(len=*), intent(in) :: key, value
      logical, intent(in), optional :: warn

      logical :: w

      if (present(warn)) then
         w = warn
      else
         w = .FALSE.
      end if
      
      if (this%has_key(key)) then
         if (w) call WarningMessage('key ' // key // ' already in dictionary.  Skipping')
         return
      end if

      call AddToVector(this%keys, varString(key))
      call AddToVector(this%values, varString(value))

   end subroutine append_str_str

   subroutine append_vstr_vstr(this, key, value, warn)
      class(Dict), intent(inout) :: this
      type(varString), intent(in) :: key, value
      logical, intent(in), optional :: warn

      logical :: w

      if (present(warn)) then
         w = warn
      else
         w = .FALSE.
      end if
      
      if (this%has_key(key)) then
         if (w) call WarningMessage('key ' // key%s // ' already in dictionary.  Skipping')
         return
      end if

      call AddToVector(this%keys, key)
      call AddToVector(this%values, value)

   end subroutine append_vstr_vstr

   subroutine append_vstr_str(this, key, value, warn)
      class(Dict), intent(inout) :: this
      type(varString), intent(in) :: key
      character(len=*), intent(in) :: value
      logical, intent(in), optional :: warn

      logical :: w

      if (present(warn)) then
         w = warn
      else
         w = .FALSE.
      end if
      
      if (this%has_key(key)) then
         if (w) call WarningMessage('key ' // key%s // ' already in dictionary.  Skipping')
         return
      end if

      call AddToVector(this%keys, key)
      call AddToVector(this%values, varString(value))

   end subroutine append_vstr_str

   subroutine append_str_vstr(this, key, value, warn)
      class(Dict), intent(inout) :: this
      character(len=*), intent(in) :: key
      type(varString), intent(in) :: value
      logical, intent(in), optional :: warn

      logical :: w

      if (present(warn)) then
         w = warn
      else
         w = .FALSE.
      end if
      
      if (this%has_key(key)) then
         if (w) call WarningMessage('key ' // key // ' already in dictionary.  Skipping')
         return
      end if

      call AddToVector(this%keys, varString(key))
      call AddToVector(this%values, value)

   end subroutine append_str_vstr

   subroutine append_str_real(this, key, value, warn)
      class(Dict), intent(inout) :: this
      character(len=*), intent(in) :: key
      real(kind=wp), intent(in) :: value
      logical, intent(in), optional :: warn

      logical :: w

      if (present(warn)) then
         w = warn
      else
         w = .FALSE.
      end if
      
      if (this%has_key(key)) then
         if (w) call WarningMessage('key ' // key // ' already in dictionary.  Skipping')
         return
      end if

      call AddToVector(this%keys, varString(key))
      call AddToVector(this%values, varString(value))

   end subroutine append_str_real

   subroutine append_vstr_real(this, key, value, warn)
      class(Dict), intent(inout) :: this
      type(varString), intent(in) :: key
      real(kind=wp), intent(in) :: value
      logical, intent(in), optional :: warn

      logical :: w

      if (present(warn)) then
         w = warn
      else
         w = .FALSE.
      end if
      
      if (this%has_key(key)) then
         if (w) call WarningMessage('key ' // key%s // ' already in dictionary.  Skipping')
         return
      end if

      call AddToVector(this%keys, key)
      call AddToVector(this%values, varString(value))

   end subroutine append_vstr_real

   subroutine append_str_int(this, key, value, warn)
      class(Dict), intent(inout) :: this
      character(len=*), intent(in) :: key
      integer, intent(in) :: value
      logical, intent(in), optional :: warn

      logical :: w

      if (present(warn)) then
         w = warn
      else
         w = .FALSE.
      end if
      
      if (this%has_key(key)) then
         if (w) call WarningMessage('key ' // key // ' already in dictionary.  Skipping')
         return
      end if

      call AddToVector(this%keys, varString(key))
      call AddToVector(this%values, varString(value))

   end subroutine append_str_int

   subroutine append_vstr_int(this, key, value, warn)
      class(Dict), intent(inout) :: this
      type(varString), intent(in) :: key
      integer, intent(in) :: value
      logical, intent(in), optional :: warn

      logical :: w

      if (present(warn)) then
         w = warn
      else
         w = .FALSE.
      end if
      
      if (this%has_key(key)) then
         if (w) call WarningMessage('key ' // key%s // ' already in dictionary.  Skipping')
         return
      end if

      call AddToVector(this%keys, key)
      call AddToVector(this%values, varString(value))

   end subroutine append_vstr_int

   function has_key_str(this, key) result(has_key)
      class(Dict), intent(in) :: this
      character(len=*), intent(in) :: key
      logical :: has_key

      has_key = InVector(this%keys, varString(key))
      
   end function has_key_str

   function has_key_vstr(this, key) result(has_key)
      class(Dict), intent(in) :: this
      type(varString), intent(in) :: key
      logical :: has_key

      has_key = InVector(this%keys, key)
      
   end function has_key_vstr

   subroutine get_real_str(this, key, val, report)
      class(Dict), intent(in) :: this
      character(len=*), intent(in) :: key
      real(kind=wp), intent(out) :: val
      logical, intent(in), optional :: report
      integer :: i

      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage('Key ' // key // ' not found in dictionary. Skipping.')
         return
      end if
         
      i = IndexInVector(this%keys, varString(key))

      val = this%values(i)%to_real()
      
   end subroutine get_real_str

   subroutine get_real_str_default(this, key, val, default, report)
      class(Dict), intent(in) :: this
      character(len=*), intent(in) :: key
      real(kind=wp), intent(out) :: val
      real(kind=wp), intent(in) :: default
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage("Key " // key // " not found. Using default value ", default)
         val = default
         return
      end if

      i = IndexInVector(this%keys, varString(key))

      val = this%values(i)%to_real()
      
   end subroutine get_real_str_default

   subroutine get_real_vstr(this, key, val, report)
      class(Dict), intent(in) :: this
      type(varString), intent(in) :: key
      real(kind=wp), intent(out) :: val
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage('Key ' // key%s // ' not found in dictionary. Skipping.')
         return
      end if

      i = IndexInVector(this%keys, key)

      val = this%values(i)%to_real()
      
   end subroutine get_real_vstr

   subroutine get_real_vstr_default(this, key, val, default, report)
      class(Dict), intent(in) :: this
      type(varString), intent(in) :: key
      real(kind=wp), intent(out) :: val
      real(kind=wp), intent(in) :: default
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage("Key " // key%s // " not found. Using default value ", default)
         val = default
         return
      end if

      i = IndexInVector(this%keys, key)

      val = this%values(i)%to_real()
      
   end subroutine get_real_vstr_default

   subroutine get_int_str(this, key, val, report)
      class(Dict), intent(in) :: this
      character(len=*), intent(in) :: key
      integer, intent(out) :: val
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage('Key ' // key // ' not found in dictionary. Skipping.')
         return
      end if

      i = IndexInVector(this%keys, varString(key))

      val = this%values(i)%to_int()
      
   end subroutine get_int_str

   subroutine get_int_str_default(this, key, val, default, report)
      class(Dict), intent(in) :: this
      character(len=*), intent(in) :: key
      integer, intent(out) :: val
      integer, intent(in) :: default
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage("Key " // key // " not found. Using default value ", default)
         val = default
         return
      end if

      i = IndexInVector(this%keys, varString(key))

      val = this%values(i)%to_int()
      
   end subroutine get_int_str_default

   subroutine get_int_vstr(this, key, val, report)
      class(Dict), intent(in) :: this
      type(varString), intent(in) :: key
      integer, intent(out) :: val
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage('Key ' // key%s // ' not found in dictionary. Skipping.')
         return
      end if

      i = IndexInVector(this%keys, key)

      val = this%values(i)%to_int()
      
   end subroutine get_int_vstr

   subroutine get_int_vstr_default(this, key, val, default, report)
      class(Dict), intent(in) :: this
      type(varString), intent(in) :: key
      integer, intent(out) :: val
      integer, intent(in) :: default
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage("Key " // key%s // " not found. Using default value ", default)
         val = default
         return
      end if

      i = IndexInVector(this%keys, key)

      val = this%values(i)%to_int()
      
   end subroutine get_int_vstr_default

   subroutine get_logical_str(this, key, val, report)
      class(Dict), intent(in) :: this
      character(len=*), intent(in) :: key
      logical, intent(out) :: val
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage('Key ' // key // ' not found in dictionary. Skipping.')
         return
      end if

      i = IndexInVector(this%keys, varString(key))

      val = this%values(i)%to_logical()
      
   end subroutine get_logical_str

   subroutine get_logical_str_default(this, key, val, default, report)
      class(Dict), intent(in) :: this
      character(len=*), intent(in) :: key
      logical, intent(out) :: val
      logical, intent(in) :: default
      logical, intent(in), optional :: report

      type(varString) :: default_vstr
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (default) then
            default_vstr = varString("True")
         else
            default_vstr = varString("False")
         end if
         if (r) call WarningMessage("Key " // key // " not found. Using default value " // default_vstr%s)
         val = default
         return
      end if

      i = IndexInVector(this%keys, varString(key))

      val = this%values(i)%to_logical()
      
   end subroutine get_logical_str_default

   subroutine get_logical_vstr(this, key, val, report)
      class(Dict), intent(in) :: this
      type(varString), intent(in) :: key
      logical, intent(out) :: val
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage('Key ' // key%s // ' not found in dictionary. Skipping.')
         return
      end if

      i = IndexInVector(this%keys, key)

      val = this%values(i)%to_logical()
      
   end subroutine get_logical_vstr

   subroutine get_logical_vstr_default(this, key, val, default, report)
      class(Dict), intent(in) :: this
      type(varString), intent(in) :: key
      logical, intent(out) :: val
      logical, intent(in) :: default
      logical, intent(in), optional :: report

      type(varString) :: default_vstr
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (default) then
            default_vstr = varString("True")
         else
            default_vstr = varString("False")
         end if
         if (r) call WarningMessage("Key " // key%s // " not found. Using default value " // default_vstr%s)
         val = default
         return
      end if

      i = IndexInVector(this%keys, key)

      val = this%values(i)%to_logical()
      
   end subroutine get_logical_vstr_default

   subroutine get_vstring_str(this, key, val, report)
      class(Dict), intent(in) :: this
      character(len=*), intent(in) :: key
      type(varString), intent(out) :: val
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage('Key ' // key // ' not found in dictionary. Skipping.')
         return
      end if

      i = IndexInVector(this%keys, varString(key))

      val = this%values(i)
      
   end subroutine get_vstring_str

   subroutine get_vstring_str_default(this, key, val, default, report)
      class(Dict), intent(in) :: this
      character(len=*), intent(in) :: key
      type(varString), intent(out) :: val
      type(varString), intent(in) :: default
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage("Key " // key // " not found. Using default value " // default%s)
         val = default
         return
      end if

      i = IndexInVector(this%keys, varString(key))

      val = this%values(i)
      
   end subroutine get_vstring_str_default

   subroutine get_vstring_vstr(this, key, val, report)
      class(Dict), intent(in) :: this
      type(varString), intent(in) :: key
      type(varString), intent(out) :: val
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage('Key ' // key%s // ' not found in dictionary. Skipping.')
         return
      end if

      i = IndexInVector(this%keys, key)

      val = this%values(i)
      
   end subroutine get_vstring_vstr

   subroutine get_vstring_vstr_default(this, key, val, default, report)
      class(Dict), intent(in) :: this
      type(varString), intent(in) :: key
      type(varString), intent(out) :: val
      type(varString), intent(in) :: default
      logical, intent(in), optional :: report
      integer :: i
      logical :: r

      if (present(report)) then
         r = report
      else
         r = .false.
      end if

      if (.not. allocated(this%keys)) then
         call FatalErrorMessage('Dictionary has no entries')
      end if

      if (.not. this%has_key(key)) then
         if (r) call WarningMessage("Key " // key%s // " not found. Using default value " // default%s)
         val = default
         return
      end if

      i = IndexInVector(this%keys, key)

      val = this%values(i)
      
   end subroutine get_vstring_vstr_default

   pure function dict_length(this) result(len)
      class(Dict), intent(in) :: this
      integer :: len

      if (.not. allocated(this%keys)) then
         len = 0
         return
      end if

      len = size(this%keys)
      return
   end function dict_length

   subroutine destructor(this)
      implicit none
      type(Dict), intent(inout) :: this

      if (allocated(this%keys)) then
         deallocate(this%keys)
      end if
      if (allocated(this%values)) then
         deallocate(this%values)
      end if
   end subroutine destructor

   subroutine dict_write_ufmt(this, unit, iostat, iomsg)
      class(Dict), intent(in) :: this
      integer, intent(in) :: unit
      integer, intent(out) :: iostat
      character(len=*), intent(inout) :: iomsg

      integer j, N

      N = this%len()

      do j=1,N
         write(unit, fmt=*, iostat=iostat, iomsg=iomsg) this%keys(j)%s, ' = ', this%values(j)%s
      end do

   end subroutine dict_write_ufmt

   function combine(this, that) result(out)
      class(Dict), intent(in) :: this
      class(Dict), intent(in) :: that
      type(Dict) :: out

      integer :: J

      do J=1,this%len()
         call out%append(this%keys(J), this%values(J), warn=.FALSE.)
      end do
      do J=1,that%len()
         call out%append(that%keys(J), that%values(J), warn=.FALSE.)
      end do
   end function combine

end module dict_module
   