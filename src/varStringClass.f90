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


! This module defines a class for variable length strings (varString)
!
! A varString is a convenient character string for Fortran.
! A varString has a single attribute (%s) containing the character string.
! The varString constructor assigns an arbitrary length character string to the %s attribute,
! so a varString v can be constructed from a character string c="test" using
!     type(varString) :: v
!     character(len=4) :: c
!     c = "test"
!     v = varString(c)
!
! A varString can be copied to a new varString using the = assignment, e.g.
!     type(varString) :: old, new
!     old = varString('Hello world')
!     new = old
!     print *, new%s ! prints "Hello world"
!
! Two varStrings a and b can be concatenated using the + operator,
! so that c = a + b creates a varstring c, e.g.
!     type(varString) :: a, b, c
!     a = varString("Hello,")
!     b = varString("World")
!     c = a + b ! gives c%s = "Hello, World"
!
! Two varStrings can be checked for equality using the == comparison, e.g.
!     type(varString) :: a, b, c
!     a = varString("Hello,")
!     b = varString("World")
!     c = a
!     print *, a==b ! prints "F"
!     print *, a==c ! prints "T"
!
! Methods defined for a varString v are:
!  N = v%len()
!     the length of the character string v%s.
!     Return: integer :: N
!
!  cstr = v%to_cstring()
!     convert a varString to a c-compatible string (NULL terminated array or length-1 c_char characters)
!     Return: character(len=1,kind=c_char) :: cstr(v%len()+1)
!
!  ext = v%get_extension()
!     get extension (substring following last '.' character) from varString
!     Return: varString :: ext
!
!  in_v = v%contains(x)
!     check if v contains a given substring
!     Inputs: character(len=*) OR type(varString :: x;
!     Return: logical :: in_v;  .TRUE. if x is found in v, otherwise in_v=.FALSE.
!
!  is_num = v%is_numeric()
!     check if v contains only numerical characters (0--9, ., +, - and e)
!     Return: logical :: is_num; .TRUE. if so, otherwise .FALSE.
!
!  nums = v%get_numeric()
!     get numerical characters (0--9, ., +, - and e) from the varString
!     Return: type(varString) :: nums
!
!  call v%split(at, first, remain)
!     split a varString at specified substring
!     Inputs: character(len=*) or type(varString) :: at -- substring at which to split
!     Returns: type(varString) :: first -- characters preceeding x; empty if x is not in v%s
!              type(varString), optional :: remain -- if present, characters following x
!
!  call v%get_between(from, to, first, remain)
!     get the part of a varstring between specified substrings passed as a varStrings.
!     Inputs:  character(len=*) or type(varString) :: from - string marking start of portion to extract
!              character(len=*) or type(varString) :: to - string marking end of portion to extract
!     Returns: type(varString) :: first - the first substring of 'this' between the strings 'from' and 'to'
!              optional varString remain - the substring of 'this' following the first occurence of string 'to'
!
!  v_new = v%adjustl()
!     Implement adjustl (adjust-left, i.e. removing initial whitespace)
!     Note, as with adjustl, leading whitespace is removed
!     but v_new%s will have trailing whitespace to be the same length as v%s
!     Returns: type(varString) :: v_new
!
!  v_new = v%trim()
!     trim (i.e. remove training whitespace)
!     Returns: type(varString) :: v_new
!
!  a = v%first_char()
!     get first character from varString
!     Returns: character(len=1) :: a
!
!  a = v%drop_comments()
!  a = v%drop_comments("%")
!     get varString ignoring characters following optional
!     comment character (default "#")
!     Returns: type(varString) :: a
!
!  N = v%count_substring(",")
!     find number of times substring "," appears in varString
!     Returns: integer :: N
!
!  a = v%to_lower()
!     convert varString to lower case
!     Returns: type(varString) :: a
!
!  a = v%to_upper()
!     convert varString to upper case
!     Returns: type(varString) :: a
!
!  N = v%to_int()
!     convert varString to integer
!     Returns: integer :: N
!
!  x = v%to_real()
!  x = v%to_real('(f8.3)')
!     convert varString to real
!     If the string is input as a fraction "a/b", the real is calculated.
!     Optional format string can be passed.
!     Returns: real(kind=wp) :: x
!     
!  x = v%in_list(list)
!     check whether varString is in a list (1d array) of varStrings
!     Input: type(varString), dimension(:) :: list
!     Returns: logical :: x
!
!  call v%read_list(list, delim)
!     convert varString containing list items separated by delim character
!     into a list (1d-array) of type of list
!     Inputs: character(len=1), optional :: delim [default delim=","]
!     Outputs: type(varString), dimension(:), allocatable :: list
!              integer, dimension(:), allocatable :: list
!              real(kind=wp), dimension(:), allocatable :: list
!
!  call v%read_Set(list)
!     convert varString containing a set -- string of the form "(a,b,c,...z)"
!     into a list (1d-array) of type list
!     Outputs: type(varString), dimension(:), allocatable :: list
!              integer, dimension(:), allocatable :: list
!              real(kind=wp), dimension(:), allocatable :: list
!
!  The module also defines the public subroutine ReadFileLine.
!  This reads a line of an open file into a varString.
module varstring_module

   use, intrinsic :: iso_c_binding
   use, intrinsic :: iso_fortran_env, only: CSZ=>character_storage_size
   use set_precision_module, only : wp

   implicit none

   private
   public :: varString
   public :: ReadFileLine

   type varString
      character(len=:), allocatable :: s
   contains
      
      procedure :: get_extension => varString_get_extension ! extracts an extension (substring following '.')

      ! make a copy (=)
      procedure, pass(this) :: varString_copy
      generic, public :: assignment(=) => varString_copy

      ! concatenate (+)
      procedure, pass(this) :: add_varstring ! concatenate two varStrings
      procedure, pass(this) :: add_string_left ! concatenate a varString and a character string
      procedure, pass(this) :: add_string_right ! concatenate a character string and a varString
      generic, public :: operator(+) => add_varstring, add_string_left, add_string_right ! Overload + operator for varStrings

      ! check equality (==)
      procedure, pass(this) :: equal_varstring ! test equality of  two varStrings
      procedure, pass(this) :: equal_varstring_left ! test equality of varString and character string
      procedure, pass(this) :: equal_varstring_right ! test equality of character string and varString
      generic, public :: operator(==) => equal_varstring, equal_varstring_left, equal_varstring_right ! Overload == comparison operator for varStrings

      ! check non-equality (/=)
      procedure, pass(this) :: not_equal_varstring ! test non-equality of  two varStrings
      procedure, pass(this) :: not_equal_varstring_left ! test non-equality of varString and character string
      procedure, pass(this) :: not_equal_varstring_right ! test non-equality of character string and varString
      generic, public :: operator(/=) => not_equal_varstring, not_equal_varstring_left, not_equal_varstring_right ! Overload /= comparison operator for varStrings

      ! Methods:

      ! len [integer function]
      ! get string length
      procedure :: len => varString_length

      ! to_cstring [character(kind=c_char) function]
      ! converts a varString into a c-compatible string 
      ! (i.e. an array of lenth-1 character strings of type c_char)
      procedure, pass(this) :: varString_to_cstring
      generic, public :: to_cstring => varString_to_cstring

      ! contains [logical function]
      ! check for substring
      procedure, pass(this) :: varstring_contains_string ! determine whether a varString contains a substring passed as a character string
      procedure, pass(this) :: varstring_contains_varstring ! determine whether a varString contains a substring passed as a varString
      generic, public :: contains => varstring_contains_string, varstring_contains_varstring ! overload "contains" to determine whether a character string contains a substring

      ! is_numeric [logical function]
      ! check if numeric
      procedure, public, pass(this) :: is_numeric => varstring_isnumeric ! determine whether a varString contains only numerical characters

      ! get_numeric [varString function]
      ! get numeric characters
      procedure, public, pass(this) :: get_numeric => varstring_get_numeric ! get numerical characters from a varString

      ! split [subroutine]
      ! split at character
      procedure, pass(this) :: varstring_split_string ! split a varString at a specified character passed as a character string
      procedure, pass(this) :: varstring_split_varstring ! split a varString at a specified character passed as a varString
      generic, public :: split => varstring_split_string, varstring_split_varstring ! define method "split" to split a varString using at character passed either as character string or varString

      ! get_between [subroutine]
      ! get string between two substrings
      procedure, pass(this) :: varstring_between_string ! find varString between two substrings passed as character strings
      procedure, pass(this) :: varstring_between_varstring ! find varString between two substrings passed as varStrings
      generic, public :: get_between => varstring_between_string ! find varString between two substrings passed as either character strings or varStrings

      ! adjustl [varString function]
      ! remove leading whitespace
      procedure, pass(this) :: varstring_adjustl ! remove leading whitespace
      generic, public :: adjustl => varstring_adjustl ! define adjustl for varStrings

      ! trim [varString function]
      ! remove trailing whitespace
      procedure, pass(this) :: varstring_trim ! remove trailing whitespace
      generic, public :: trim => varstring_trim ! define trim for varStrings

      ! first_char [character(len=1) function]
      ! get first character
      procedure, pass(this) :: varstring_first_char ! get first character of a varString
      generic, public :: first_char => varstring_first_char ! get first character of a varString

      ! drop_comments [varString function]
      ! remove comments following substring ()
      procedure, pass(this) :: varstring_drop_comments ! strip comments out from a varString
      generic, public :: drop_comments => varstring_drop_comments

      ! count_substring [integer function]
      ! count number of occurrences of substring
      procedure, pass(this) :: varstring_count_varstring
      procedure, pass(this) :: varstring_count_substring
      generic, public :: count_substring => varstring_count_varstring, varstring_count_substring

      ! to_lower [varString function]
      ! convert to lower case
      procedure, public, pass(this) :: to_lower => varstring_to_lower

      ! to_upper [varString function]
      ! convert to upper case
      procedure, public, pass(this) :: to_upper => varstring_to_upper

      ! to_int
      ! convert to int
      procedure, public, pass(this) :: to_int => varstring_to_int

      ! convert to real (to_real)
      procedure, pass(this) :: varstring_to_real
      procedure, pass(this) :: varstring_to_real_fmt
      generic, public :: to_real => varstring_to_real, varstring_to_real_fmt

      ! check if this varstring is in a list (1d array) of varstrings
      procedure, public, pass(this) :: in_list => varstring_in_list

      ! read a list (delimiter separated varString) to array
      procedure, pass(this) :: varstring_list_to_varstrings
      procedure, pass(this) :: varstring_list_to_ints
      procedure, pass(this) :: varstring_list_to_reals
      generic, public :: read_list => varstring_list_to_varstrings, varstring_list_to_ints, varstring_list_to_reals

      ! read a set -- varString of the form (a,b,c,d,...) -- to array
      procedure, pass(this) :: varstring_set_to_varstrings
      procedure, pass(this) :: varstring_set_to_ints
      procedure, pass(this) :: varstring_set_to_reals
      generic, public :: read_set => varstring_set_to_varstrings, varstring_set_to_ints, varstring_set_to_reals


   end type varString

   ! Interface to varString constructor
   interface varString
      procedure constructor_char, constructor_real, constructor_int
   end interface varString

contains

   ! Create a varString from a inferred length character string
   ! Inputs: character string str
   ! Returns: varString this - varString containing str in attribute this%s
   pure function constructor_char(str, trim_str) result(this)
      character(len=*), intent(in) :: str
      logical, intent(in), optional :: trim_str
      type(varString) :: this

      logical :: do_trim

      if (present(trim_str)) then
         do_trim = trim_str
      else
         do_trim = .FALSE.
      end if
      
      if (do_trim) then
         this%s = trim(str)
      else
         this%s = str
      end if
   end function constructor_char

   pure function constructor_int(val) result(this)
      integer, intent(in) :: val
      type(varString) :: this

      integer :: isize
      character(len=:), allocatable :: valStr

      isize = storage_size(val)/CSZ

      if (allocated(valStr)) deallocate(valStr)
      allocate(character(len=isize) :: valStr)
      valStr(:) = transfer(val, valStr)

      this%s = valStr

   end function constructor_int

   pure function constructor_real(val) result(this)
      real(kind=wp), intent(in) :: val
      type(varString) :: this

      integer :: isize
      character(len=:), allocatable :: valStr

      isize = storage_size(val)/CSZ

      if (allocated(valStr)) deallocate(valStr)
      allocate(character(len=isize) :: valStr)
      valStr(:) = transfer(val, valStr)

      this%s = valStr

   end function constructor_real
   
   ! Private method to copy a varString
   !  new_varString = old_varString
   pure subroutine varString_copy(this, from)
      class(varString), intent(inout) :: this
      class(varString), intent(in) :: from
      this%s = from%s
   end subroutine varString_copy

   ! Private method to concatenate this varString with another
   !  out = this%add_varstring(txt)
   ! Inputs: varString txt
   ! Returns: varString out - out = this + txt
   pure function add_varstring(this, txt) result(out)
      class(varString), intent(in) :: this
      type(varString), intent(in) :: txt
      type(varString) :: out

      character(len=len(this%s)) :: tmp0
      character(len=len(txt%s)) :: tmp1

      tmp0 = this%s
      tmp1 = txt%s

      out = varString(tmp0 // tmp1)

   end function add_varstring

   ! Private method to concatenate this varString with a character string
   !  out = this%add_varstring(txt)
   ! Inputs: character(len=:) txt
   ! Returns: varString out - out = this + txt
   pure function add_string_left(this, txt) result(out)
      class(varString), intent(in) :: this
      character(len=*), intent(in) :: txt
      type(varString) :: out

      character(len=len(this%s)) :: tmp0

      tmp0 = this%s

      out = varString(tmp0 // txt)

   end function add_string_left

   ! Private method to concatenate this varString with a character string
   !  out = this%add_varstring(txt)
   ! Inputs: character(len=:) txt - string to concatenate with this
   ! Returns: varString out - out = txt + this
   pure function add_string_right(txt, this) result(out)
      character(len=*), intent(in) :: txt
      class(varString), intent(in) :: this
      type(varString) :: out

      character(len=len(this%s)) :: tmp0

      tmp0 = this%s

      out = varString(txt // tmp0)

   end function add_string_right

   ! Private method to compare this varString with varString txt
   !  isequal = this == txt
   pure function equal_varstring(this, txt) result(isequal)
      class(varString), intent(in) :: this
      type(varString), intent(in) :: txt
      logical isequal

      isequal = (this%s == txt%s)
   end function equal_varstring

   ! Private method to compare this varString with character string txt
   !  isequal = this == txt
   pure function equal_varstring_left(this, txt) result(isequal)
      class(varString), intent(in) :: this
      character(len=*), intent(in) :: txt
      logical isequal

      isequal = (this%s == txt)
   end function equal_varstring_left

   ! Private method to compare character string txt with this varString
   !  isequal = this == txt
   pure function equal_varstring_right(txt, this) result(isequal)
      character(len=*), intent(in) :: txt
      class(varString), intent(in) :: this
      logical isequal

      isequal = (this%s == txt)
   end function equal_varstring_right

   ! Private method to compare this varString with varString txt
   !  notisequal = this /= txt
   pure function not_equal_varstring(this, txt) result(isnotequal)
      class(varString), intent(in) :: this
      type(varString), intent(in) :: txt
      logical isnotequal

      isnotequal = (this%s /= txt%s)
   end function not_equal_varstring

   ! Private method to compare this varString with character string txt
   !  notisequal = this /= txt
   pure function not_equal_varstring_left(this, txt) result(notisequal)
      class(varString), intent(in) :: this
      character(len=*), intent(in) :: txt
      logical notisequal

      notisequal = (this%s /= txt)
   end function not_equal_varstring_left

   ! Private method to compare character string txt with this varString
   !  isequal = this == txt
   pure function not_equal_varstring_right(txt, this) result(isnotequal)
      character(len=*), intent(in) :: txt
      class(varString), intent(in) :: this
      logical isnotequal

      isnotequal = (this%s /= txt)
   end function not_equal_varstring_right

   ! Determine the length of the string in this%s
   ! Called as N = this%len()
   ! Inputs: 
   ! Returns: int N - the length of the character string in this%s
   pure function varString_length(this) result(N)
      class(varString), intent(in) :: this
      integer :: N

      N = len(this%s)

   end function varString_length

   ! Determine if varString contains given substring
   ! Called as out = this%contains(str)
   ! Inputs: character string str - substring
   ! Returns: logical out - out = .TRUE. if str is a substring of this%s, otherwise .FALSE.
   pure function varString_contains_string(this, str) result(out)
      ! Check if varstring contains substring
      class(varString), intent(in) :: this
      character(len=*), intent(in) :: str
      logical :: out

      if (index(this%s, str)==0) then
         out = .false.
      else
         out = .true.
      end if
   end function varString_contains_string

   ! Determine if varString contains given substring
   ! Called as out = this%contains(str)
   ! Inputs: varString str - substring
   ! Returns: logical out - out = .TRUE. if str is a substring of this%s, otherwise .FALSE.
   pure function varString_contains_varstring(this, str) result(out)
      ! Check if varstring contains substring
      class(varString), intent(in) :: this
      type(varString), intent(in) :: str
      logical :: out

      if (index(this%s, str%s)==0) then
         out = .false.
      else
         out = .true.
      end if
   end function varString_contains_varstring

   ! Determine if varstring contains only numeric characters.
   ! Note numeric characters include .,+,- and e
   ! Called as out = this%is_numeric()
   ! Inputs:
   ! Returns: logical out - out=.TRUE. if this%s contains only numeric characters, otherwise .FALSE.
   pure function varstring_isnumeric(this) result(out)
      class(varString), intent(in) :: this
      logical :: out

      integer :: J

      do J=1,len_trim(this%s)
         if (.not.((lle('0',this%s(J:J)).and.lle(this%s(J:J),'9')) &
            .or.(this%s(J:J)=='.') &
            .or.(this%s(J:J)=='-') &
            .or.(this%s(J:J)=='+') &
            .or.(this%s(J:J)=='e'))) then
            out = .false.
            return
         end if
      end do
      out = .true.
   end function varstring_isnumeric

   ! Find the numerical characters that appear together in varString
   ! (with no gaps and no alphabetic characters, other than e, within).
   ! Inputs: 
   ! Returns: varString out - numerical characters in this%s
   pure function varstring_get_numeric(this) result(out)
      class(varString), intent(in) :: this
      type(varString) :: out

      integer :: J
      integer :: K

      out = varString('')

      K=0
      do J=1,this%len()
         if ((lle('0',this%s(J:J)).and.lle(this%s(J:J),'9')) &
            .or.(this%s(J:J)=='.') &
            .or.(this%s(J:J)=='-') &
            .or.(this%s(J:J)=='+') &
            .or.(this%s(J:J)=='e') &
            ) then
            K = K+1
            out = out + this%s(J:J)
         else
            if (K>0) exit
         end if
      end do
   end function varstring_get_numeric

   ! Split a varstring at a specified substring substring passed as a string
   ! Called as "call this%split(at, first, remain)"
   ! Inputs: character string at - substring at which to split
   ! Returns: varString first - the substring of 'this' preceeding the string 'at'
   !          optional varString remain - the substring of 'this' following the string 'at' 
   pure subroutine varString_split_string(this, at, first, remain)
      class(varString), intent(in) :: this
      character(len=*), intent(in) :: at
      type(varString), intent(out) :: first
      type(varString), intent(out), optional :: remain

      integer j

      if (.not. this%contains(at)) then
         first = this
         if (present(remain)) remain=varString('')
      else
         j = scan(this%s, at)
         first = varString(this%s(1:j-1))
         first = first%adjustl()
         first = first%trim()
         if (present(remain)) then
            remain = varString(adjustl(this%s(j+1:)))
            remain = remain%adjustl()
            remain = remain%trim()
         end if
      end if

   end subroutine varString_split_string

   ! Split a varstring at a specified substring substring passed as a varString
   ! Called as "call this%split(at, first, remain)"
   ! Inputs:  varString at - substring at which to split
   ! Returns: varString first - the substring of 'this' preceeding the string 'at'
   !          optional varString remain - the substring of 'this' following the string 'at' 
   pure subroutine varString_split_varstring(this, at, first, remain)
      ! Check if varstring contains substring
      class(varString), intent(in) :: this
      type(varString), intent(in) :: at
      type(varString), intent(out) :: first
      type(varString), intent(out), optional :: remain

      integer j

      if (.not. this%contains(at%s)) then
         first = this
         if (present(remain)) remain=varString('')
      else
         j = scan(this%s, at%s)
         first = varString(this%s(1:j))
         first = first%adjustl()
         first = first%trim()
         if (present(remain)) then
            remain = varString(this%s(j+1:))
            remain = remain%adjustl()
            remain = remain%trim()
         end if
      end if

   end subroutine varString_split_varstring

   ! Get portion of varstring between specified substrings passed as a character string
   ! Called as "call this%get_between(from, to, first, remain)"
   ! Inputs:  character string from - mark start of portion to extract
   !          character string to - mark end of portion to extract
   ! Returns: varString first - the substring of 'this' between the strings 'from' and 'to'
   !          optional varString remain - the substring of 'this' following the string 'to' 
   pure subroutine varString_between_string(this, from, to, first, remain)
      class(varString), intent(in) :: this
      character(len=*), intent(in) :: from
      character(len=*), intent(in) :: to
      type(varString), intent(out) :: first
      type(varString), intent(out), optional :: remain

      integer i0, i1

      if (.not. this%contains(from)) then
         first = varString('')
         if (present(remain)) remain = this
      elseif (.not. this%contains(to)) then
         first = varString('')
         if (present(remain)) remain = this
      else
         i0 = scan(this%s, from)
         i1 = scan(this%s, to)

         first = varString(this%s(i0+1:i1-1))
         if (present(remain)) remain = varString(this%s(i1:))
      end if

   end subroutine varString_between_string

   ! Get portion of varstring between specified substrings passed as a varStrings
   ! Called as "call this%get_between(from, to, first, remain)"
   ! Inputs:  varString from - mark start of portion to extract
   !          varString to - mark end of portion to extract
   ! Returns: varString first - the substring of 'this' between the strings 'from' and 'to'
   !          optional varString remain - the substring of 'this' following the string 'to' 
   pure subroutine varString_between_varString(this, from, to, first, remain)
      class(varString), intent(in) :: this
      type(varString), intent(in) :: from
      type(varString), intent(in) :: to
      type(varString), intent(out) :: first
      type(varString), intent(out), optional :: remain

      integer i0, i1

      if (.not. this%contains(from)) then
         first = varString('')
         if (present(remain)) remain = this
      elseif (.not. this%contains(to)) then
         first = varString('')
         if (present(remain)) remain = this
      else
         i0 = scan(this%s, from%s)
         i1 = scan(this%s, to%s)

         first = varString(this%s(i0+1:i1-1))
         if (present(remain)) remain = varString(this%s(i1+1:))
      end if

   end subroutine varString_between_varString

   ! Convert a varString to a c-compatible string
   ! Called as cstr = this%to_cstring()
   ! Returns: character(len=1, kind=c_char) :: cstr(v%len()+1)
   pure function varString_to_cstring(this) result(cString)

      implicit none

      class(varString), intent(in) :: this
      character(len=1,kind=c_char) :: cString(len(this%s)+1)

      integer :: N, i

      N = this%len()
      do i = 1, N
         cString(i) = this%s(i:i)
      end do
      cString(N + 1) = C_NULL_CHAR

   end function varString_to_cstring

   ! Get extension (substring following last '.' character) from varString
   ! Called as ext = this%get_extension()
   ! Returns: varString :: ext
   pure function varString_get_extension(this) result(ext)
      class(varString), intent(in) :: this
      type(varString) :: ext
      type(varString) :: str
      type(varString) :: first
      type(varString) :: remain

      logical :: split

      if (.not. this%contains('.')) then
         ext = varString('')
         return
      end if

      str = this
      split = .True.
      do
         call str%split('.', first, remain)
         if (.not. remain%contains('.')) exit
         str = remain
      end do

      ext = remain

   end function varString_get_extension

   ! Implement adjustl intrinsic for varString
   ! Called as adjusted = this%adjustl()
   ! Returns: varString
   pure function varstring_adjustl(this) result(adjusted)
      class(varString), intent(in) :: this
      type(varString) :: adjusted

      adjusted = varString(adjustl(this%s))
   end function varstring_adjustl

   ! Implement trim intrinsic for varString
   ! Called as trimmed = this%trim()
   ! Returns: varString
   pure function varstring_trim(this) result(trimmed)
      class(varString), intent(in) :: this
      type(varString) :: trimmed

      trimmed = varString(trim(this%s))
   end function varstring_trim

   ! Get first character for a varString
   ! Called as c = this%first_char()
   ! Returns: character(len=1) :: c
   pure function varstring_first_char(this) result(c)
      class(varString), intent(in) :: this
      character(len=1) :: c

      c = this%s(1:1)
   end function varstring_first_char

   ! Drop comments from a varString specified by character
   ! Called as out = this%drop_comments(comment_char)
   ! Inputs: optional character string comment_char - character marking comments (default = #)
   ! Returns: varString out - the varString preceeding the comment character
   pure function varString_drop_comments(this, comment_char_in) result(out)
      class(varString), intent(in) :: this
      character(len=1), intent(in), optional :: comment_char_in
      type(varString) :: out

      character(len=1) :: comment_char

      if (present(comment_char_in)) then
         comment_char = comment_char_in
      else
         comment_char = "#"
      end if

      call this%split(comment_char, out)

      out = out%trim()
      
   end function varString_drop_comments

   ! Count number of occurrences of substring s in a varString specified by character
   ! Called as N = this%count_substring(str)
   ! Inputs: character string comment_char - character marking comments (default = #)
   ! Returns: integer N - number of times substring occurs
   pure function varstring_count_varstring(this, s) result(N)
      class(varString), intent(in) :: this
      type(varString), intent(in) :: s
      integer :: N
      integer :: p, posn

      N = 0
      if (s%len() == 0) return
      p = 1
      do
         posn = index(this%s(p:), s%s)
         if(posn == 0) return
         N = N + 1
         p = p + posn + s%len()-1
      end do

   end function varstring_count_varstring

   ! Count number of occurrences of substring s in a varString specified by character
   ! Called as N = this%count_substring(str)
   ! Inputs: character string comment_char - character marking comments (default = #)
   ! Returns: integer N - number of times substring occurs
   pure function varstring_count_substring(this, s) result(N)
      class(varString), intent(in) :: this
      character(*), intent(in) :: s
      integer :: N
      integer :: p, posn

      N = 0
      if (len(s) == 0) return
      p = 1
      do
         posn = index(this%s(p:), s)
         if(posn == 0) return
         N = N + 1
         p = p + posn + len(s)-1
      end do

   end function varstring_count_substring

   ! Convert this to lower case
   ! called as str = this%to_lower()
   ! Input: this - varString
   ! Returns: str - varString
   pure function varstring_to_lower(this) result(str)
      class(varString), intent(in) :: this
      type(varString) :: str
      integer :: i
      character(len=len(this%s)) :: str_tmp

      do i = 1, this%len()
         select case(this%s(i:i))
          case("A":"Z")
            str_tmp(i:i) = achar(iachar(this%s(i:i))+32)
          case default
            str_tmp(i:i) = this%s(i:i)
         end select
      end do
      str = varString(str_tmp)
   end function varstring_to_lower

   ! Convert this to upper case
   ! called as str = this%to_upper()
   ! Input: this - varString
   ! Returns: str - varString
   pure function varstring_to_upper(this) result(str)
      class(varString), intent(in) :: this
      type(varString) :: str
      integer :: i
      character(len=len(this%s)) :: str_tmp

      do i = 1, this%len()
         select case(this%s(i:i))
          case("a":"z")
            str_tmp(i:i) = achar(iachar(this%s(i:i))-32)
          case default
            str_tmp(i:i) = this%s(i:i)
         end select
      end do
      str = varString(str_tmp)
   end function varstring_to_upper

   ! Convert this containing an integer to an integer.
   ! Called as val = this%to_int()
   ! If the string is blank, a fatal error occurs.
   ! Input: this - varstring
   ! Output: val - integer
   function varstring_to_int(this) result(val)
      class(varString), intent(in) :: this
      integer :: val

      integer :: stat
      type(varString) :: msg

      read(this%s, *, iostat=stat) val
      if (stat.ne.0) then
         msg = varString('In varstring to_int, could not read varString ')
         msg = msg + this%s
         msg = msg + ' as integer. Read status =' + varString(stat)
         error stop msg%s
      end if

   end function varstring_to_int

   ! Convert this containing an real number to a real.
   ! Called as val = this%to_real()
   ! If the string is blank, a fatal error occurs.
   ! If the string is input as a fraction "a/b", the real is calculated.
   ! Input: this - varstring
   ! Output: val - integer
   function varstring_to_real(this) result(val)
      class(varString), intent(in) :: this
      real(kind=wp) :: val

      type(varString) :: numStr
      type(varString) :: denStr
      real(kind=wp) :: num, den

      type(varString) :: msg
      integer :: stat

      if (this%contains('/')) then
         
         call this%split('/', numStr, denStr)
         read(numStr%s, *,iostat=stat) num
         read(denStr%s, *,iostat=stat) den
         val = num/den
      else
         read(this%s, *, iostat=stat) val
      end if

      if (stat.ne.0) then
         msg = varString('In varstring to_real, could not read varString ')
         msg = msg + this%s
         msg = msg + ' as real. Read status =' + varString(stat)
         error stop msg%s
      end if

   end function varstring_to_real

   ! Convert this containing an real number to a real using a format specifier
   ! Called as val = this%to_real()
   ! If the string is blank, a fatal error occurs.
   ! Inputs: this - varstring
   !         fmt - format specifier in standard notation (e.g. fmt='(f5.3)')
   ! Output: val - real
   function varstring_to_real_fmt(this, fmt) result(val)
      class(varString), intent(in) :: this
      character(len=*), intent(in) :: fmt
      real(kind=wp) :: val

      integer :: stat
      type(varString) :: msg

      read(this%s, fmt=fmt, iostat=stat) val

      if (stat.ne.0) then
         msg = varString('In varstring to_real, could not read varString ')
         msg = msg + this%s
         msg = msg + ' as real. Read status =' + varString(stat)
         error stop msg%s
      end if

   end function varstring_to_real_fmt

   ! Check if this is a member of a list of varStrings
   ! Called as inlist = this%in_list(vlist)
   ! Inputs: this - varstring
   !         vlist - a 1d array of varStrings
   ! Returns: inlist - logical: inlist = .TRUE. if this is in vlist, otherwise .FALSE.
   function varstring_in_list(this, vlist) result(inlist)
      class(varString), intent(in) :: this
      type(varString), dimension(:), intent(in) :: vlist
      logical :: inlist

      integer J

      do J=1,size(vlist)
         if (this == vlist(J)) then
            inlist = .TRUE.
            return
         end if
      end do

      inlist = .FALSE.
   end function varstring_in_list

   ! Convert a varString containing a list (delimiter separated strings)
   ! into a 1d array of varStrings.
   ! Called as call this%to_list(vals, delimiter=",")
   ! Inputs: this - varString
   !         delimiter - optional character - delimiter separating items [default ","]
   ! In/Out: vals - 1d array of integers, reallocated if allocated on input
   subroutine varstring_list_to_varstrings(this, vals, delimiter)

      class(varString), intent(in) :: this
      character(len=1), intent(in), optional :: delimiter
      type(varString), dimension(:), allocatable, intent(out) :: vals

      character(len=1) :: delim

      type(varString) :: list, remain

      integer :: N, J

      if (present(delimiter)) then
         delim = delimiter
      else
         delim = ","
      end if

      ! Get number of items
      N = this%count_substring(delim)+1
      if (allocated(vals)) deallocate(vals)
      allocate(vals(N))
      if (N==1) then
         vals(1) = this
         return
      end if

      list = this

      do J=1,N-1
         call list%split(delim, vals(J), remain=remain)
         list = remain
      end do
      vals(N) = list

   end subroutine varstring_list_to_varstrings

   ! Convert a varString containing a list of integers (delimiter separated strings)
   ! into a 1d array of integers.
   ! Called as call this%to_list(vals, delimiter=",")
   ! Inputs: this - varString
   !         delimiter - optional character - delimiter separating items [default ","]
   ! In/Out: vals - 1d array of integers, reallocated if allocated on input
   subroutine varstring_list_to_ints(this, vals, delimiter)

      class(varString), intent(in) :: this
      character(len=1), intent(in), optional :: delimiter
      integer, dimension(:), allocatable, intent(out) :: vals

      character(len=1) :: delim

      type(varString) :: list, item, remain

      integer :: N, J

      if (present(delimiter)) then
         delim = delimiter
      else
         delim = ","
      end if

      ! Get number of items
      N = this%count_substring(delim)+1

      if (allocated(vals)) deallocate(vals)
      allocate(vals(N))

      if (N==1) then
         vals(1) = this%to_int()
         return
      end if

      list = this

      do J=1,N-1
         call list%split(delim, item, remain=remain)
         vals(J) = item%to_int()
         list = remain
      end do
      vals(N) = list%to_int()

   end subroutine varstring_list_to_ints

   ! Convert a varString containing a list of real numbers (delimiter separated strings)
   ! into a 1d array of reals.
   ! Called as call this%to_list(vals, delimiter=",")
   ! Inputs: this - varString
   !         delimiter - optional character - delimiter separating items [default ","]
   ! In/Out: vals - 1d array of real(kind=wp), reallocated if allocated on input
   subroutine varstring_list_to_reals(this, vals, delimiter)

      class(varString), intent(in) :: this
      character(len=1), intent(in), optional :: delimiter
      real(kind=wp), dimension(:), allocatable, intent(out) :: vals

      character(len=1) :: delim

      type(varString) :: list, item, remain

      integer :: N, J

      if (present(delimiter)) then
         delim = delimiter
      else
         delim = ","
      end if

      ! Get number of items
      N = this%count_substring(delim)+1

      if (allocated(vals)) deallocate(vals)
      allocate(vals(N))

      if (N==1) then
         vals(1) = this%to_real()
         return
      end if

      list = this

      do J=1,N-1
         call list%split(delim, item, remain=remain)
         vals(J) = item%to_real()
         list = remain
      end do
      vals(N) = list%to_real()

   end subroutine varstring_list_to_reals

   ! Convert a varString containing a set -- comma-separated string of the form (a,b,c,...,z) --
   ! into a 1d array of varStrings
   ! Called as call this%read_set(vals)
   ! Inputs: this - varString
   ! In/Out: vals - 1d array of varStrings, reallocated if allocated on input
   subroutine varstring_set_to_varstrings(this, vals)

      class(varString), intent(in) :: this
      type(varString), dimension(:), allocatable, intent(out) :: vals

      type(varString) :: list
      type(varString) :: msg

      if (.not. (this%contains('(') .and. this%contains(')'))) then
         msg = varString('varString ')
         msg = msg + this%s
         msg = msg + 'is not recognized as a set'
         error stop msg%s
      end if

      call this%get_between('(', ')', list)

      call list%read_list(vals, delimiter=',')
      
   end subroutine varstring_set_to_varstrings

   ! Convert a varString containing a set -- comma-separated string of the form (a,b,c,...,z) --
   ! of integers into a 1d array of integers
   ! Called as call this%read_set(vals)
   ! Inputs: this - varString
   ! In/Out: vals - 1d array of integers, reallocated if allocated on input
   subroutine varstring_set_to_ints(this, vals)

      class(varString), intent(in) :: this
      integer, dimension(:), allocatable, intent(out) :: vals

      type(varString) :: list
      type(varString) :: msg

      if (.not. (this%contains('(') .and. this%contains(')'))) then
         msg = varString('varString ')
         msg = msg + this%s
         msg = msg + 'is not recognized as a set'
         error stop msg%s
      end if

      call this%get_between('(', ')', list)

      call list%read_list(vals, delimiter=',')
      
   end subroutine varstring_set_to_ints

   ! Convert a varString containing a set -- comma-separated string of the form (a,b,c,...,z) --
   ! of reals into a 1d array of reals
   ! Called as call this%read_set(vals)
   ! Inputs: this - varString
   ! In/Out: vals - 1d array of reals, reallocated if allocated on input
   subroutine varstring_set_to_reals(this, vals)

      class(varString), intent(in) :: this
      real(kind=wp), dimension(:), allocatable, intent(out) :: vals

      type(varString) :: list
      type(varString) :: msg

      if (.not. (this%contains('(') .and. this%contains(')'))) then
         msg = varString('varString ')
         msg = msg + this%s
         msg = msg + 'is not recognized as a set'
         error stop msg%s
      end if

      call this%get_between('(', ')', list)

      call list%read_list(vals, delimiter=',')
      
   end subroutine varstring_set_to_reals

   ! Read a line of an opened file into a varString
   ! The file is passed as an integer unit and must have been opened
   ! Note, returns on io error
   ! Inputs: integer :: unit - the file unit number for an open file
   !         integer, optional :: maxlen - maximum number of characters to read.
   !                                       If not present, defaults to huge(1)
   ! Returns: type(varString) :: vstring - a varString containing the line from file
   !          integer :: iostat - the iostat from the file read, =0 for no problems
   subroutine ReadFileLine(unit, vstring, iostat, maxlen)
      integer, intent(in) :: unit
      type(varString), intent(out) :: vstring
      integer, intent(out) :: iostat
      integer, intent(in), optional :: maxlen

      integer, parameter :: buffer_size = 100
      character(len=buffer_size) :: buffer
      integer :: next_read_length
      integer :: num_read
      integer :: num_to_read

      if (present(maxlen)) then
         num_to_read = maxlen
      else
         num_to_read = huge(1)
      end if

      vstring = varString('')

      do
         if (num_to_read <= 0) exit
         next_read_length = min(buffer_size, num_to_read)
         read(unit, fmt='(A)', advance='NO', eor=9999, size=num_read, iostat=iostat) buffer(1:next_read_length)
         if (iostat /= 0) return
         vstring = vstring + buffer(1:next_read_length)
         num_to_read = num_to_read - next_read_length
      end do
      return
      9999 vstring = vstring + buffer(1:num_read)
      iostat = 0
   end subroutine ReadFileLine
    
end module varstring_module
