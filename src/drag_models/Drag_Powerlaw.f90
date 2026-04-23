submodule (drag_types) drag_powerlaw

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: FlowSquaredSpeedSlopeAligned, Reciprocal

   implicit none

   ! -- Drag closure: Powerlaw drag --
   !
   !    Powerlaw drag has friction function F = C * (|u|/H)^n.
   !
   !    Parameters:
   !        power_law_co > 0
   !        power_law_power > 0
   !
contains

   ! Constructor
   module function powerlaw_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(PowerlawDrag) :: this
      this%name = varString("Powerlaw")
      if (params%has_key('powerlaw co')) then
         call params%get('powerlaw co', this%power_law_co)
      else
         call WarningMessage('Powerlaw Co not given.  Using default value ', this%power_law_co)
      end if
      if (params%has_key('powerlaw_power')) then
         call params%get('powerlaw_power', this%power_law_power)
      else
         call WarningMessage('Powerlaw Power not given.  Using default value ', this%power_law_power)
      end if

      ! Parameters
      call this%parameters%append("power law co", this%power_law_co)
      call this%parameters%append("power law power", this%power_law_power)
   end function powerlaw_constructor

   ! Validator
   module subroutine powerlaw_validate(this)
      class(PowerlawDrag), intent(in) :: this

      if (this%power_law_co<0) then
         call FatalErrorMessage('Power law Co must be positive.  Received', this%power_law_co)
      end if
      
   end subroutine powerlaw_validate

   ! Friction
   module pure function powerlaw_friction(this, RunParams, uvect) result(friction)
      class(PowerlawDrag), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: modu, Hn, Hneps, hr

      
      modu = sqrt(FlowSquaredSpeedSlopeAligned(RunParams, uvect))
      Hn = uvect(RunParams%Vars%Hn)
      Hneps = RunParams%heightThreshold

      hr = Reciprocal(Hn, RunParams%heightThreshold)
      friction = this%power_law_co * (modu * hr)**this%power_law_power

   end function powerlaw_friction

end submodule drag_powerlaw