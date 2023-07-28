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


! This module defines limiters for computing derivatives of potentially
! discontinous fields.  The chosen limiter is determined at runtime by the user.
!
! The core hydrodynamic solver algorithm (in HydraulicsRHS.f90) uses slope limiters
! to approximate derivatives while accounting for discontinuities that occur in 
! hyperbolic systems of conservation laws.  The slope limiters in general take
! as inputs forward and backward finite difference approximations of a derivative 
! at a point and using a limiter function determine the derivative from these.
!
! The generic 'limiter' function takes as inputs a and b, the forward and
! backward differences, respectively, and returns a derivative approximate computed
! from these.
!
! The limiter function is called frequently, and therefore efficiencies in the
! specific limiter implementation can result in cumulative benefits in overall
! run time.
!
! Selection of limiters occurs in SolverSettings.f90
!
! Currently implemented limiters are:
! - MinMod1: minmod limiter
! - MinMod2: [recommended] generalized minmod limiter with parameter theta=1.3
! - vanAlbada: van Albada smooth limiter with parameter epsilon=0
! - WENO: weighted essentially non-oscillatory limiter
! - LimiterNone: no limiter, gives central difference
module limiters_module

   use set_precision_module, only: wp

   implicit none

   private

   public :: limiter
   public :: MinMod1
   public :: MinMod2
   public :: vanAlbada
   public :: WENO
   public :: LimiterNone
   
   ! Generic limiter function
   ! Various limiters are selected via this interface by Solver_Set.
   ! This must be a pure function.
   ! Inputs: a - forward slope approximation
   !         b - backward slope approximation
   ! Returns: c - limited slope approximation
   pointer :: limiter
   interface
      pure function limiter(a, b) result(c)
         import :: wp
         real(kind=wp), intent(in) :: a, b
         real(kind=wp) :: c
      end function limiter
   end interface

contains

   ! minmod limiter defined as
   ! minmod(a, b) = (sgn(a) + sgn(b))*min(abs(a), abs(b))/2
   ! where sgn(a) is the signum function.
   ! This formula can be more simply implemented.
   pure function MinMod1(a, b) result(c)
      real(kind=wp), intent(in) :: a, b
      real(kind=wp) :: c
        
      if (a*b<=0.0_wp) then ! either (i) a>0, b<0; (ii) a<0, b>0; (iii) a and/or b = 0
         c = 0.0_wp
      else
         if (a > 0.0_wp) then ! a>0, b>0
            c = min(a, b)
         else  ! a<0 so b<0
            c = max(a, b)
         end if
      end if
      return
   end function MinMod1

   ! minmod2 limiter is a generalization of minmod and is defined as
   ! minmod2(a, b) = sgn(a) * min(theta*abs(a), abs(a+b)/2, theta*abs(b)), if sgn(a)=sgn(b)
   !               = 0, if sgn(a)+sgn(b) = 0
   ! The parameter 1<=theta<=2.  Here we impose theta = 1.3 which we have found to work well.
   ! This formula can be more simply implemented.
   ! @cite https://doi.org/10.1002/fld.4023
   pure function MinMod2(a, b) result(c)
      real(kind=wp), intent(in) :: a, b
      real(kind=wp) :: c
      real(kind=wp), parameter :: theta = 1.3_wp

      if (a*b<=0.0_wp) then ! either (i) a>0, b<0; (ii) a<0, b>0; (iii) a and/or b = 0
         c = 0.0_wp
      else
         if (a > 0.0_wp) then ! a>0, b>0
            c = min(theta * a, min(theta * b, 0.5_wp * (a + b)))
         else  ! a<0 so b<0
            c = max(theta * a, max(theta * b, 0.5_wp * (a + b)))
         end if
      end if

   end function MinMod2

   ! no limiter
   ! Use central approximation c = (a + b)/2
   pure function LimiterNone(a, b) result(c)
      real(kind=wp), intent(in) :: a,b
      real(kind=wp) :: c

      c = 0.5_wp*(a+b)
      return

   end function LimiterNone

   ! van Albada limiter
   ! Smooth limiter.
   ! Tends to central difference in smooth regions
   ! Across discontinuities the averaged slope is biased to smallest
   ! value of the two one-sided slopes.
   ! In general van Albada limiter has a parameter.
   ! Here we use a version with parameter set to zero.
   ! @cite https://doi.org/10.48550/arXiv.2304.00437
   pure function vanAlbada(a, b) result(c)
      real(kind=wp), intent(in) :: a,b
      real(kind=wp) :: c

      real(kind=wp) :: den

      den = a*a + b*b

      if (den==0.0_wp) then
         c = 0.0_wp
      else
         c = (a*a*b + a*b*b)/den
      end if

   end function

   pure function weno_func(a) result(c)

      implicit none

      real(kind=wp), intent(in) :: a
      real(kind=wp) :: c
      real(kind=wp), parameter :: eps = 1.0e-6_wp
      real(kind=wp) :: epsaa

      epsaa = a*a+eps
      c = 1.0_wp/(epsaa*epsaa)

   end function weno_func

   pure function WENO(a, b) result(c)

      implicit none

      real(kind=wp), intent(in) :: a,b
      real(kind=wp) :: c
      real(kind=wp) :: wa, wb

      wa = weno_func(a)
      wb = weno_func(b)

      c = (wa*a + wb*b)/(wa + wb)

   end function WENO

end module limiters_module
