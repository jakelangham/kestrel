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