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


! A variety of routines for printing messages, comprising:
!
! fatal error messages (which terminate the program), written to stderr warnings
! (might be a problem, so we let the user know), written to stdout info messages
! (information the user might be interested in), written to stdout time step
! messages (completed time steps, and time step revision), written to stdout
module messages_module

   use, intrinsic :: iso_fortran_env, only : stdout=>output_unit, stderr=>error_unit
   use set_precision_module, only: wp
   
   implicit none

   private
   public :: FatalErrorMessage, WarningMessage, InfoMessage
   public :: TimestepMessage, RevisedTimestepMessage
   public :: InputLabelUnrecognized

   interface FatalErrorMessage
      module procedure FatalErrorMessage_b, FatalErrorMessage_m, FatalErrorMessage_r, FatalErrorMessage_rv, FatalErrorMessage_i, FatalErrorMessage_is
   end interface

   interface WarningMessage
      module procedure WarningMessage_m, WarningMessage_r, WarningMessage_rv, WarningMessage_i
   end interface

   interface InfoMessage
      module procedure InfoMessage_m, InfoMessage_r
   end interface InfoMessage

contains

   ! Fatal Error message with no information.
   ! Terminates with status 1.
   subroutine FatalErrorMessage_b
      write(stderr, *) '                                    '
      write(stderr, *) '************************************'
      write(stderr, *) '************************************'
      write(stderr, *) 'FATAL ERROR :                       '
      write(stderr, *) 'Quitting.'
      write(stderr, *) '************************************'
      write(stderr, *) '************************************'
      write(stderr, *) '                                    '

      call exit(1) ! Fatal error return with status 1

   end subroutine FatalErrorMessage_b

   ! Fatal Error message with information.
   ! Input: Message - reported message
   ! Terminates with status 1.
   subroutine FatalErrorMessage_m(Message)
      character(*), intent(in) :: Message

      write(stderr, *) '                                    '
      write(stderr, *) '************************************'
      write(stderr, *) '************************************'
      write(stderr, *) 'FATAL ERROR :                       '
      write(stderr, *) trim(Message)
      write(stderr, *) 'Quitting.'
      write(stderr, *) '************************************'
      write(stderr, *) '************************************'
      write(stderr, *) '                                    '

      call exit(1) ! Fatal error return with status 1

   end subroutine FatalErrorMessage_m

   ! Fatal Error message with information and real value.
   ! Input: Message - reported message
   ! Terminates with status 1.
   subroutine FatalErrorMessage_r(Message, val)

      implicit none

      character(*), intent(in) :: Message
      real(kind=wp), intent(in) :: val

      write(stderr,*) '                                    '
      write(stderr,*) '************************************'
      write(stderr,*) '************************************'
      write(stderr,*) 'FATAL ERROR :                       '
      write(stderr,*) trim(Message), val
      write(stderr,*) 'Quitting.'
      write(stderr,*) '************************************'
      write(stderr,*) '************************************'
      write(stderr,*) '                                    '

      call exit(1) ! Fatal error return with status 1

   end subroutine FatalErrorMessage_r

   ! Fatal Error message with information and integer value.
   ! Inputs: Message - reported message
   !         val - integer value appended to message
   ! Terminates with status 1.
   subroutine FatalErrorMessage_i(Message, val)

      character(*), intent(in) :: Message
      integer, intent(in) :: val

      write(stderr,*) '                                    '
      write(stderr,*) '************************************'
      write(stderr,*) '************************************'
      write(stderr,*) 'FATAL ERROR :                       '
      write(stderr,*) trim(Message), val
      write(stderr,*) 'Quitting.'
      write(stderr,*) '************************************'
      write(stderr,*) '************************************'
      write(stderr,*) '                                    '

      call exit(1) ! Fatal error return with status 1

   end subroutine FatalErrorMessage_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Fatal Error message with information. integer value and more 
   ! information.
   ! Input: Message1 - first part of reported message
   !        val - integer value in message
   !        Message2 - second part of reported message
   ! Terminates with status 1.
   subroutine FatalErrorMessage_is(Message1, val, Message2)
      character(*), intent(in) :: Message1
      integer, intent(in) :: val
      character(*), intent(in) :: Message2

      write(stderr,*) '                                    '
      write(stderr,*) '************************************'
      write(stderr,*) '************************************'
      write(stderr,*) 'FATAL ERROR :                       '
      write(stderr,*) trim(Message1), val, trim(Message2)
      write(stderr,*) 'Quitting.'
      write(stderr,*) '************************************'
      write(stderr,*) '************************************'
      write(stderr,*) '                                    '

      call exit(1) ! Fatal error return with status 1

   end subroutine FatalErrorMessage_is

   ! Fatal Error message with information and real valued array.
   ! Inputs: Message - reported message
   !         val - 1d array of reals printed after Message
   ! Terminates with status 1.
   subroutine FatalErrorMessage_rv(Message, val)
      character(*), intent(in) :: Message
      real(kind=wp), dimension(:), intent(in) :: val

      write(stderr,*) '                                    '
      write(stderr,*) '************************************'
      write(stderr,*) '************************************'
      write(stderr,*) 'FATAL ERROR :                       '
      write(stderr,*) trim(Message)
      write(stderr,*) val
      write(stderr,*) 'Quitting.'
      write(stderr,*) '************************************'
      write(stderr,*) '************************************'
      write(stderr,*) '                                    '

      call exit(1) ! Fatal error return with status 1

   end subroutine FatalErrorMessage_rv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Warning message with information.
   ! Input: Message - reported message.
   subroutine WarningMessage_m(Message)
      character(*), intent(in) :: Message

      write(stdout,*) 'Warning: ' // Message
      return
   end subroutine WarningMessage_m

   ! Warning message with information and real value.
   ! Optional format string can be passed
   ! Inputs: Message - reported message
   !         val - real number to append
   !         fmt [optional] - format string in usual syntax
   subroutine WarningMessage_r(Message, val, fmt)
      character(*), intent(in) :: Message
      real(kind=wp), intent(in) :: val
      character(*), intent(in), optional :: fmt

      if (present(fmt)) then
         write(stdout,fmt=fmt) 'Warning: ' // trim(Message), val
      else
         write(stdout,*) 'Warning: ' // trim(Message), val
      end if
      return
   end subroutine WarningMessage_r

   ! Warning message with information and integer value.
   ! Inputs: Message - reported message
   !         val - integer to append
   subroutine WarningMessage_i(Message, val)
      character(*), intent(in) :: Message
      integer, intent(in) :: val

      write(stdout,*) 'Warning: ' // trim(Message), val      
      return
   end subroutine WarningMessage_i

   ! Warning message with information and real valued 1d-array.
   ! Optional format string can be passed
   ! Inputs: Message - reported message
   !         val - 1d-array of real numbers to append
   !         fmt [optional] - format string in usual syntax
   subroutine WarningMessage_rv(Message, val, fmt)
      character(*), intent(in) :: Message
      real(kind=wp), dimension(:), intent(in) :: val
      character(*), intent(in), optional :: fmt

      if (present(fmt)) then
         write(stdout,fmt=fmt) 'Warning: ' // trim(Message), val
      else
         write(stdout,*) 'Warning: ' // trim(Message), val
      end if
      return
   end subroutine WarningMessage_rv

   ! Info message with information.
   ! Input: Message - reported message
   subroutine InfoMessage_m(Message)
      character(*), intent(in) :: Message

      write(stdout,*) Message
      return
   end subroutine InfoMessage_m

   ! Info message with information and real value.
   ! Optional format string can be passed
   ! Inputs: Message - reported message
   !         val - real number to append
   !         fmt [optional] - format string in usual syntax
   subroutine InfoMessage_r(Message, val, fmt)
      character(*), intent(in) :: Message
      real(kind=wp), intent(in) :: val
      character(len=*), intent(in), optional:: fmt

      if (present(fmt)) then
         write(stdout,fmt=fmt) Message, val
      else
         write(stdout,*) Message, val
      end if
      return
   end subroutine InfoMessage_r

   ! Write timestep info to stdout
   ! Inputs: t - current time
   !         dt - time step
   subroutine TimestepMessage(t, dt)
      real(kind=wp), intent(in) :: t
      real(kind=wp), intent(in) :: dt
   
      write(stdout, fmt='(a4, f20.12, a7, f20.12)') &
         't = ', t, ', dt = ', dt

      return
   end subroutine TimestepMessage

   ! Write revised timestep info to stdout
   ! Inputs: t - revised next time
   !         dt - new time step
   subroutine RevisedTimestepMessage(t, dt)
      real(kind=wp), intent(in) :: t
      real(kind=wp), intent(in) :: dt
   
      write (stdout, fmt='(a16,f13.7,a10,f10.7)') &
         'revised nextT = ', t, ' new dt = ', dt
      return
   end subroutine RevisedTimestepMessage

   ! Display warning that an input label is not recognized.
   ! Input: Label - label that is not recognized.
   subroutine InputLabelUnrecognized(Label)
      character(*), intent(in) :: Label

      write(stderr,*) 'Warning: Input ' // trim(Label) // ' not recognized. Skipping.'
      return
   end subroutine InputLabelUnrecognized

end module messages_module
