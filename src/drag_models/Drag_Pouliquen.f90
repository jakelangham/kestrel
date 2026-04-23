submodule (drag_types) drag_pouliquen

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: FlowSquaredSpeedSlopeAligned, GeometricCorrectionFactor

   implicit none

   ! -- Drag closure: Pouliquen drag --
   !
   !    Pouliquen drag has friction function F = \mu(I) * g * Hn.
   !    with Pouliquen's intertial number (I) dependent friction coefficient, mu(I).
   !
   !    Parameters:
   !        pouliquen_min > 0
   !        pouliquen_max > 0
   !        pouliquen_beta > 0
   !        solid_diameter > 0
   !
contains

   ! Constructor
   module function pouliquen_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(PouliquenDrag) :: this
      this%name = varString("Pouliquen")
      if (params%has_key('pouliquen min')) then
         call params%get('pouliquen min', this%pouliquen_min)
      else
         call WarningMessage('Pouliquen Min not given.  Using default value ', this%pouliquen_min)
      end if
      if (params%has_key('pouliquen max')) then
         call params%get('pouliquen max', this%pouliquen_max)
      else
         call WarningMessage('Pouliquen Max not given.  Using default value ', this%pouliquen_max)
      end if
      if (params%has_key('pouliquen beta')) then
         call params%get('pouliquen beta', this%pouliquen_beta)
      else
         call WarningMessage('Pouliquen Beta not given.  Using default value ', this%pouliquen_beta)
      end if
      if (params%has_key('solid diameter')) then
         call params%get('solid diameter', this%solid_diameter)
      else
         call FatalErrorMessage('Required parameters Solid Diameter not given.')
      end if

      ! Parameters
      call this%parameters%append("pouliquen min", this%pouliquen_min)
      call this%parameters%append("pouliquen max", this%pouliquen_max)
      call this%parameters%append("pouliquen beta", this%pouliquen_beta)
   end function pouliquen_constructor

   ! Validator
   module subroutine pouliquen_validate(this)
      class(PouliquenDrag), intent(in) :: this
      
      if (this%pouliquen_min<=0) then
         call FatalErrorMessage('Pouliquen Min must be positive.  Received', this%pouliquen_min)
      end if

      if (this%pouliquen_max<=0) then
         call FatalErrorMessage('Pouliquen max must be positive.  Received', this%pouliquen_max)
      end if

      if (this%pouliquen_beta<=0) then
         call FatalErrorMessage('Pouliquen beta must be positive.  Received', this%pouliquen_beta)
      end if

      if (this%solid_diameter<=0) then
         call FatalErrorMessage('Solid Diameter must be positive.  Received', this%solid_diameter)
      end if
   end subroutine pouliquen_validate

   ! Friction
   module pure function pouliquen_friction(this, RunParams, uvect) result(friction)
      class(PouliquenDrag), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: Hn
      real(kind=wp) :: modu2
      real(kind=wp) :: gam
      real(kind=wp) :: g
      real(kind=wp) :: mu

      Hn = uvect(RunParams%Vars%Hn)
      gam = GeometricCorrectionFactor(RunParams, uvect)
      g = RunParams%g / gam

      modu2 = FlowSquaredSpeedSlopeAligned(RunParams, uvect)

      mu = this%friction_coefficient(RunParams%heightThreshold, g, Hn, sqrt(modu2))

      friction = mu * g * Hn
      
   end function pouliquen_friction

   ! Friction coefficient
   module pure function pouliquen_friction_coefficient(this, heightThreshold, gcostheta, Hn, modu) result(mu)
      class(PouliquenDrag), intent(in) :: this
      real(kind=wp), intent(in) :: heightThreshold
      real(kind=wp), intent(in) :: gcostheta
      real(kind=wp), intent(in) :: Hn
      real(kind=wp), intent(in) :: modu
      real(kind=wp) :: mu

      real(kind=wp) :: Fr, I

      if (Hn > heightThreshold) then
         ! Froude number
         Fr = modu / sqrt(gcostheta * Hn)
         ! Inertial number
         I = Fr * this%solid_diameter / Hn
         mu = this%pouliquen_min + (this%pouliquen_max - this%pouliquen_min) * I / (this%pouliquen_beta + I)
      else
         mu = this%pouliquen_min
      end if
   end function pouliquen_friction_coefficient

end submodule drag_pouliquen