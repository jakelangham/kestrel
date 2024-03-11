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


! This module defines analytical formulae for topographic elevation.
! The chosen formula is determined at runtime by the user.
! 
! An interface to a generic TopogFunc is defined the specific choices
! is set in TopogSettings.f90.
!
! TopogFunc is a pure subroutine of the form
!   subroutine TopogFunc(RunParams, x, y, b0)
! with inputs: type(RunSet) :: RunParams - parameters of the model
!              real(kind=wp), dimension(:) :: x - the Easting coordinate
!              real(kind=wp), dimension(:) :: y - the Northing coordinate
! and output: real(kind=wp), dimension(:,:) :: b0(size(x),size(y)) - the topographic elevation at x, y
! Note, for 1d simulations, size(y)=1.
! Each specific subroutine has the same form.
!
! Currently implemented formulae are:
!   x2slopes - Two slopes connected by a circular arc
!   usgs - Parameterization of the USGS flume
!   flume - General flume based on the USGS flume
!   xbislope - Two slopes with smooth transition
!   flat - Flat topography
!   xslope - Uniform slope along x
!   yslope - Uniform slope along y
!   xyslope - Uniform slope along x-y
!   xsinslope - one-dimensional sinusoidal variation
!   xysinslope - two-dimensional sinusoidal variation
!   xhump - one-dimensional cosine hump on a flat topography
!   xtanh - one-dimensional tanh on a flat topography
!   xparab - one-dimensional parabola
!   xyparab - two-dimensional parabola
module topog_funcs_module

   use set_precision_module, only: wp, pi
   use runsettings_module, only: RunSet

   implicit none

   private
   public :: TopogFunc ! Generic template for topography formulae
   public :: x2slopes ! Two slopes connected by a circular arc
   public :: usgs ! Parameterization of the USGS flume
   public :: flume ! General flume based on the USGS flume
   public :: quadratic_flume ! Flume with runout and quadratic cross-section
   public :: xbislope ! Two slopes with smooth transition
   public :: flat ! Flat topography
   public :: xslope ! Uniform slope along x
   public :: yslope ! Uniform slope along y
   public :: xyslope ! Uniform slope along x-y
   public :: xsinslope ! one-dimensional sinusoidal variation
   public :: xysinslope ! two-dimensional sinusoidal variation
   public :: xhump ! one-dimensional cosine hump on a flat topography
   public :: xtanh ! one-dimensional tanh on a flat topography
   public :: xparab ! one-dimensional parabola
   public :: xyparab ! two-dimensional parabola

   pointer :: TopogFunc
   interface
       pure subroutine TopogFunc(RunParams, x, y, b0)
           import :: wp
           import :: RunSet
           type(RunSet), intent(in) :: RunParams
           real(kind=wp), dimension(:), intent(in) :: x
           real(kind=wp), dimension(:), intent(in) :: y
           real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))
       end subroutine TopogFunc
   end interface

contains

   ! Two slopes connected by a circular arc
   ! Three parameters required:
   !  alpha -- slope on left-hand-side ; passed in RunParams%TopogFuncParams(1)
   !  beta -- slope on right-hand-side; passed in RunParams%TopogFuncParams(2)
   !  R -- radius of connecting cicular arc; passed in RunParams%TopogFuncParams(3)
   pure subroutine x2slopes(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))

      real(kind=wp) :: alpha, beta, R
      real(kind=wp) :: xc0, zc0, x1, z1, x2, z2

      integer :: ii

      alpha = RunParams%TopogFuncParams(1)
      beta = RunParams%TopogFuncParams(2)
      R = RunParams%TopogFuncParams(3)
      
      ! Determine centre of connecting circular arc
      xc0 = (sqrt(1.0_wp+alpha*alpha) - sqrt(1.0_wp+beta*beta))*R/(alpha-beta)
      zc0 = (alpha*sqrt(1.0_wp+beta*beta) - beta*sqrt(1.0_wp+alpha*alpha))*R/(alpha-beta)
      ! Determine points where arc joins slopes on left (x1,z1) and right (x2,z2)
      x1 = xc0 - alpha*R/sqrt(1.0_wp+alpha*alpha)
      z1 = zc0 - R/sqrt(1.0_wp+alpha*alpha)
      x2 = xc0 - beta*R/sqrt(1.0_wp+beta*beta)
      z2 = zc0 - R/sqrt(1.0_wp+beta*beta)

      do ii=1,size(x)
         if (x(ii)<x1) then
            b0(ii,:) = -alpha*x(ii)
         elseif (x(ii)>x2) then
            b0(ii,:) = -beta*x(ii)
         else
            b0(ii,:) = zc0 - sqrt(R*R - (x(ii)-xc0)*(x(ii)-xc0))
         end if
      end do
      return
   end subroutine x2slopes

   ! Parameterization of the USGS flume
   ! This has slope of 31 deg for x<0, and slope 2.4 deg for x>x1>0
   ! that are connected by a smooth cosh curve section.
   ! Note x1 is determined to ensure smooth connection.
   ! The channel is confined by walls for x<8.5 m, that are represented as tanh profile humps
   ! Two parameters required:
   !  wallH -- height of confining walls; passed in RunParams%TopogFuncParams(1)
   !  sigma -- 'width' of tanh walls; passed in RunParams%TopogFuncParams(2)
   pure subroutine usgs(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))

      real(kind=wp) :: alpha, xc0, zc0, x1, sigma, wallH
      integer :: ii, jj
   
      ! Set parameters
      alpha = 8.5_wp/(asinh(-tan(4.0_wp*pi/180_wp))-asinh(-tan(31.0_wp*pi/180_wp)))
      xc0 = - alpha*asinh(-tan(31.0_wp*pi/180.0_wp))
      zc0 = -alpha*cosh((-xc0)/alpha)
      x1 = xc0 + alpha*asinh(-tan(2.4_wp*pi/180.0_wp))
      wallH = RunParams%TopogFuncParams(1)
      sigma = RunParams%TopogFuncParams(2)

      ! Build topography along x
      do ii=1,size(x)
         if (x(ii)<0.0_wp) then
            b0(ii,:) = -tan(31.0_wp*pi/180.0_wp)*x(ii)
         elseif (x(ii)>x1) then
            b0(ii,:) = zc0 + alpha*cosh((x1-xc0)/alpha) - tan(2.4_wp*pi/180.0_wp)*(x(ii)-x1)
         else
            b0(ii,:) = zc0 + alpha*cosh((x(ii)-xc0)/alpha)
         end if
      end do

      ! Add confining walls
      do ii=1,size(x)
         if (x(ii)<8.5_wp) then
            do jj=1,size(y)
               b0(ii,jj) = b0(ii,jj) + 0.5_wp*wallH*(tanh(sigma*(y(jj)-1.0_wp)) &
                  - tanh(sigma*(y(jj)-2.0_wp)) + tanh(sigma*(y(jj)+2.0_wp)) - tanh(sigma*(y(jj)+1.0_wp)))
            end do
         end if
      end do
      return
   end subroutine usgs

   ! Parameterization of the general flume
   ! This has slope theta0 for x<0, and slope theta1 for x>x1>0
   ! that are connected by a smooth cosh curve section.
   ! Note x1 is determined to ensure smooth connection.
   ! The channel is confined by walls for x<xwall, that are represented as tanh profile humps
   !  parameters required:
   !  theta0 -- slope (in degrees) for x<0; passed in RunParams%TopogFuncParams(1)
   !  theta1 -- slope (in degrees) for x>x1; passed in RunParams%TopogFuncParams(2)
   !  xwall -- x coordinate for end of wall; passed in RunParams%TopogFuncParams(3)
   !  wallW -- width of confining walls; passed in RunParams%TopogFuncParams(4)
   !  wallH -- height of confining walls; passed in RunParams%TopogFuncParams(5)
   !  sigma -- 'width' of tanh walls; passed in RunParams%TopogFuncParams(6)
   pure subroutine flume(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))

      real(kind=wp) :: theta0, theta1, xwall, wallW, wallH, sigma
      real(kind=wp) :: alpha, xc0, zc0, x1
      integer :: ii, jj

      ! Get parameters
      theta0 = RunParams%TopogFuncParams(1)
      theta1 = RunParams%TopogFuncParams(2)
      xwall = RunParams%TopogFuncParams(3)
      wallW = RunParams%TopogFuncParams(4)
      wallH = RunParams%TopogFuncParams(5)
      sigma = RunParams%TopogFuncParams(6)
   
      ! Set parameters
      alpha = 8.5_wp/(asinh(-tan(4.0_wp*pi/180_wp))-asinh(-tan(theta0*pi/180_wp)))
      xc0 = - alpha*asinh(-tan(theta0*pi/180.0_wp))
      zc0 = -alpha*cosh((-xc0)/alpha)
      x1 = xc0 + alpha*asinh(-tan(theta1*pi/180.0_wp))
      
      ! Build topography along x
      do ii=1,size(x)
         if (x(ii)<0.0_wp) then
            b0(ii,:) = -tan(theta0*pi/180.0_wp)*x(ii)
         elseif (x(ii)>x1) then
            b0(ii,:) = zc0 + alpha*cosh((x1-xc0)/alpha) - tan(theta1*pi/180.0_wp)*(x(ii)-x1)
         else
            b0(ii,:) = zc0 + alpha*cosh((x(ii)-xc0)/alpha)
         end if
      end do

      ! Add confining walls
      do ii=1,size(x)
         if (x(ii)<xwall) then
            do jj=1,size(y)
               b0(ii,jj) = b0(ii,jj) + 0.5_wp*wallH*(tanh(sigma*(y(jj)-0.5_wp*wallW)) &
                  - tanh(sigma*(y(jj)-1.5_wp*wallW)) + tanh(sigma*(y(jj)+1.5_wp*wallW)) - tanh(sigma*(y(jj)+0.5_wp*wallW)))
            end do
         end if
      end do
      return
   end subroutine flume

   ! Parameterization of a quadratic flume
   ! This has slope theta0 for x<0, and slope theta1 for x>x1>0
   ! that are connected by a smooth cosh curve section.
   ! Note x1 is determined to ensure smooth connection.
   ! The channel cross-section is quadratic b(x,y) = -theta0*x + a*y^2
   !  parameters required:
   !  theta0 -- slope (in degrees) for x<0; passed in RunParams%TopogFuncParams(1)
   !  theta1 -- slope (in degrees) for x>x1; passed in RunParams%TopogFuncParams(2)
   !  xwall -- x coordinate for end of wall; passed in RunParams%TopogFuncParams(3)
   !  a -- quadratic coefficient; passed in RunParams%TopogFuncParams(4)
   pure subroutine quadratic_flume(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))

      real(kind=wp) :: theta0, theta1, xwall, a
      real(kind=wp) :: alpha, xc0, zc0, x1
      integer :: ii, jj

      ! Get parameters
      theta0 = RunParams%TopogFuncParams(1)
      theta1 = RunParams%TopogFuncParams(2)
      xwall = RunParams%TopogFuncParams(3)
      a = RunParams%TopogFuncParams(4)
      
      ! Set parameters
      alpha = 8.5_wp/(asinh(-tan(4.0_wp*pi/180_wp))-asinh(-tan(theta0*pi/180_wp)))
      xc0 = - alpha*asinh(-tan(theta0*pi/180.0_wp))
      zc0 = -alpha*cosh((-xc0)/alpha)
      x1 = xc0 + alpha*asinh(-tan(theta1*pi/180.0_wp))
      
      ! Build topography along x
      do ii=1,size(x)
         if (x(ii)<0.0_wp) then
            b0(ii,:) = -tan(theta0*pi/180.0_wp)*x(ii)
         elseif (x(ii)>x1) then
            b0(ii,:) = zc0 + alpha*cosh((x1-xc0)/alpha) - tan(theta1*pi/180.0_wp)*(x(ii)-x1)
         else
            b0(ii,:) = zc0 + alpha*cosh((x(ii)-xc0)/alpha)
         end if
      end do

      ! Add confining walls
      do ii=1,size(x)
         if (x(ii)<xwall) then
            do jj=1,size(y)
               b0(ii,jj) = b0(ii,jj) + a*y(jj)*y(jj)
            end do
         end if
      end do
      return
   end subroutine quadratic_flume

   ! Two slopes with smooth transition
   ! Three parameters required:
   !  phi1 -- angle of slope in degrees on left-hand-side; passed in RunParams%TopogFuncParams(1)
   !  phi2 -- angle of slope in degrees on right-hand-side; passed in RunParams%TopogFuncParams(2)
   !  lam -- length scale of connecting smooth curve; passed in RunParams%TopogFuncParams(3)
   pure subroutine xbislope(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))
   
      real(kind=wp) :: phi1, phi2, lam, a1, a2
      integer :: ii

      ! Set parameters
      phi1 = RunParams%TopogFuncParams(1) * pi / 180.0_wp
      phi2 = RunParams%TopogFuncParams(2) * pi / 180.0_wp
      lam = RunParams%TopogFuncParams(3)
      a1 = tan(phi1)
      a2 = tan(phi2)

      ! Build topography
      do ii = 1, size(x)
         b0(ii,:) = -0.5_wp * (a1 + a2) * x(ii) + &
            0.5_wp * (a1 - a2) * lam * log(cosh(x(ii) / lam))
      end do
   return
   end subroutine xbislope

   ! Flat topography
   pure subroutine flat(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))
      
      b0(:,:) = 0.0_wp
      return
   end subroutine flat

   ! Uniform slope along x
   ! One parameter required:
   !  slope -- slope; passed in RunParams%TopogFuncParams(1)
   pure subroutine xslope(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))

      real(kind=wp) :: slope
      integer :: ii
      slope = RunParams%TopogFuncParams(1)
   
      do ii=1,size(x)
         b0(ii,:) = slope*x(ii)
      end do
      return
   end subroutine xslope

   ! Uniform slope along y
   ! One parameter required:
   !  slope -- slope; passed in RunParams%TopogFuncParams(1)
   pure subroutine yslope(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))

      real(kind=wp) :: slope
      integer :: jj
      slope = RunParams%TopogFuncParams(1)
      
      do jj=1,size(y)
         b0(:,jj) = slope*y(jj)
      end do
      return
   end subroutine yslope

   ! Uniform slope along x-y
   ! Two parameter required:
   !  xgrad -- slope along x; passed in RunParams%TopogFuncParams(1)
   !  ygrad -- slope along y; passed in RunParams%TopogFuncParams(2)
   pure subroutine xyslope(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))

      real(kind=wp) :: xgrad, ygrad
      integer :: ii, jj
      xgrad = RunParams%TopogFuncParams(1)
      ygrad = RunParams%TopogFuncParams(2)
   
      do ii=1,size(x)
         do jj=1,size(y)
            b0(ii,jj) = xgrad*x(ii) + ygrad*y(jj)
         end do
      end do
      return
   end subroutine xyslope

   ! one-dimensional sinusoidal variation
   ! One parameter required:
   !  eps -- amplitude of the sinusoid; passed in RunParams%TopogFuncParams(1)
   pure subroutine xsinslope(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))
   
      real(kind=wp) :: eps, Lx, twopi_o_Lx
      integer :: ii

      eps = RunParams%TopogFuncParams(1) ! previously slope*eps
      Lx = RunParams%Xtilesize * RunParams%nXtiles
      twopi_o_Lx = 2.0_wp * pi / Lx
      do ii=1,size(x)
         b0(ii,:) = eps * sin(x(ii) * twopi_o_Lx)
      end do
      return
   end subroutine xsinslope

   ! two-dimensional sinusoidal variation
   ! One parameter required:
   !  eps -- amplitude of the sinusoid; passed in RunParams%TopogFuncParams(1)
   pure subroutine xysinslope(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))
   
      real(kind=wp) :: eps, Lx, Ly, twopi_o_Lx, twopi_o_Ly
      integer :: ii, jj

      eps = RunParams%TopogFuncParams(1) ! previously slope*eps
   
      Lx = RunParams%Xtilesize * RunParams%nXtiles
      Ly = RunParams%Ytilesize * RunParams%nYtiles
      twopi_o_Lx = 2.0_wp * pi / Lx
      twopi_o_Ly = 2.0_wp * pi / Ly
      do ii=1,size(x)
         do jj=1,size(y)
            b0(ii,jj) = eps * sin(x(ii) * twopi_o_Lx) * sin(y(jj) * twopi_o_Ly)
         end do
      end do
      return
   end subroutine xysinslope

   ! one-dimensional cosine hump on a flat topography
   ! Two parameters required:
   !  A -- amplitude of hump; passed in RunParams%TopogFuncParams(1)
   !  L -- length scale of hump; passed in RunParams%TopogFuncParams(2)
   pure subroutine xhump(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))

      real(kind=wp) :: A, L
      integer :: ii

      A = RunParams%TopogFuncParams(1)
      L = RunParams%TopogFuncParams(2)

      do ii=1,size(x)
         if ((x(ii)>-L) .and. (x(ii)<L)) then
            b0(ii,:) = 0.5_wp*A*(1.0_wp+cos(pi*x(ii)/L))
         else
            b0(ii,:) = 0.0_wp
         end if
      end do
      return
   end subroutine xhump

   ! one-dimensional tanh on a flat topography b = A*(1 + tanh((x-x0)/L))
   ! Three parameters required:
   !  x0 -- centre of tanh transition; passed in RunParams%TopogFuncParams(1)
   !  A -- amplitude of change in topography; passed in RunParams%TopogFuncParams(2)
   !  L -- length scale of transition; passed in RunParams%TopogFuncParams(3)
   pure subroutine xtanh(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))

      real(kind=wp) :: x0, A, L
      integer :: ii

      x0 = RunParams%TopogFuncParams(1) 
      A = RunParams%TopogFuncParams(2)
      L = RunParams%TopogFuncParams(3)
   
      do ii=1,size(x)
         b0(ii,:) = A*(1.0_wp+tanh((x(ii)-x0)/L))
      end do
      return
   end subroutine xtanh

   ! one-dimensional parabola b = A*x^2
   ! One parameter required:
   !  A -- coefficient; passed in RunParams%TopogFuncParams(1)
   pure subroutine xparab(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))
      
      real(kind=wp) :: A
      integer :: ii

      A = RunParams%TopogFuncParams(1)
      
      do ii=1,size(x)
         b0(ii,:) = A*x(ii)*x(ii)
      end do
      return
   end subroutine xparab

   ! two-dimensional parabola b = A*x^2 + B*y^2
   ! Two parameters required:
   !  A -- coefficient of x^2; passed in RunParams%TopogFuncParams(1)
   !  B -- coefficient of y^2; passed in RunParams%TopogFuncParams(2)
   pure subroutine xyparab(RunParams, x, y, b0)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: x
      real(kind=wp), dimension(:), intent(in) :: y
      real(kind=wp), dimension(:,:), intent(out) :: b0(size(x),size(y))
      
      real(kind=wp) :: A, B
      integer :: ii, jj

      A = RunParams%TopogFuncParams(1)
      B = RunParams%TopogFuncParams(2)

      do ii=1,size(x)
         do jj=1,size(y)
            b0(ii,jj) = A*x(ii)*x(ii) + B*y(jj)*y(jj)
         end do
      end do
      return
   end subroutine xyparab

end module topog_funcs_module
